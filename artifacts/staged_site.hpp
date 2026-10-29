#ifndef VG_STAGED_SITE_HPP_INCLUDED
#define VG_STAGED_SITE_HPP_INCLUDED

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "descent_policy.hpp"
#include "linkage_model.hpp"
#include "panel_lookup.hpp"
#include "read_likelihood_caller.hpp"
#include "site_values.hpp"
#include "../src/snarl_caller.hpp"

namespace vg::multipass {

using namespace std;

/**
 * Where a site sits on the reference: on a named path at an offset, or, off the reference, at a
 * stand-in position taken from its parent.
 */
struct ReferencePlacement {
    /// The reference path the site is on.
    ///
    /// A site with `no_reference` takes the path of its nearest ancestor that a reference path
    /// visits.
    string path_name;
    /// The offset added to positions on `path_name`.
    int64_t offset = 0;
    /// Whether no reference path visits the bounds of the site.
    ///
    /// Such a site has no REF or POS, so it gets no line. It is still genotyped and filed in the
    /// linkage model.
    bool no_reference = false;
    /// For a site with `no_reference`, the position that stands in for its own.
    ///
    /// It is the reference start of the parent plus `NestingPlacement::parent_offset`.
    int64_t position_from_parent = 0;

    /**
     * Get the locus the linkage model orders the site by.
     *
     * A site the reference path visits is placed where its first boundary node starts on the
     * path, or at its own base for a site inside a node. A site with `no_reference` is placed at
     * `position_from_parent`.
     *
     * @param net The site in `sites`, read forward along `path_name`.
     */
    SiteLocus locus(const PathPositionHandleGraph& graph,
                    const EditSnarlDecompositionOverlay& sites, const net_handle_t& net) const;
};

/// The alleles a site's record spells, by traversal, with the sites under it spliced in (see
/// `AlleleSplicer`).
///
/// A traversal the splice did not change is left as the site has it.
struct SplicedAlleles {
    vector<Traversal> travs;
    AlleleSequences sequences;
};

/// The genotyping result of one site from the direct pass, kept until the render builds the
/// site's records.
///
/// The ploidy of a nested chain depends on the genotype chosen for its parent, which is known only
/// after the direct pass. So the result holds the answer at both ploidies (see
/// `alt_ploidy_info`), and the render can write the record at whichever ploidy the linkage pass
/// gives the site.
///
/// The passes find a site by its record key.
struct StagedSite {
    // The name of the site, and its place in the nesting tree.

    /// The VCF `ID` of the site's records, without the suffix of a block record.
    ///
    /// A site the decomposition holds is named by its bounds (`vg::site_name`). Other sites have
    /// IDs of their own.
    string id;
    /// The hash of `id` (`vg::record_key_of`), which the tables that follow a site from pass to
    /// pass are keyed by.
    size_t record_key = 0;
    /**
     * The place of the site in the nesting tree, or the default for a site staged as top level.
     *
     * Its `level` is capped at 255, the deepest level the linkage model has.
     *
     * The linkage pass revises some fields of a nested chain from the genotype chosen for its
     * parent. It moves `parent_offset` to where the chain starts along that genotype, so that an
     * off-reference chain is placed along an allele of it. It computes `parent_crossing` and
     * `crossing_known` again for the children of a chain it revises or first files. Under the
     * linkage model it decides `reported_inline` again, which holds back the line of the site as
     * `reference.no_reference` does.
     *
     * Its `no_reference` says that an ancestor chain is off the reference, which can hold even
     * when a reference path visits this site (see `reference`).
     */
    NestingPlacement placement;
    // The site in the graph.

    /// The site's net handle in the overlay of sites, read forward along the reference path.
    ///
    /// This may be the decomposition's orientation turned round.
    net_handle_t net;
    /// The candidate alleles.
    CandidateAlleles alleles;
    /// The sites nested in this one, for its symbolic alleles.
    SiteChildren children;

    /// Return true if the site has no child sites.
    ///
    /// A site holds a child site wherever one of its child chains has two nodes, which makes an
    /// entry in `children`.
    bool is_leaf() const { return children.entries.empty(); }
    /// Where the site is on the reference.
    ReferencePlacement reference;

    // The call of the direct pass, and what the linkage pass decided about it.

    /// The genotype from the direct pass, before the linkage model.
    vector<int> genotype;
    /// The ploidy of `genotype`.
    ///
    /// For a nested chain, it is the number of the parent's called alleles that cross the chain,
    /// or the parent's ploidy when none does. `score` also holds the answer at the other ploidy.
    int ploidy = 2;
    /// The genotyper's call, which the record is written from.
    unique_ptr<SnarlCaller::CallInfo> call_info;
    /// The object `call_info` holds, typed as a read-likelihood score, or null if another genotyper
    /// called the site.
    ///
    /// `call_info` owns it, so `set_call` and `set_score` set the two together.
    SiteScore* score = nullptr;

    /// Give the site the genotyper's call `new_call`.
    ///
    /// `typed` is null or `new_call` as a read-likelihood score.
    void set_call(unique_ptr<SnarlCaller::CallInfo> new_call, SiteScore* typed);
    /// Give the site a read-likelihood score as its call.
    void set_score(unique_ptr<SiteScore> typed);
    /// Take the site's read-likelihood score out, leaving it no call.
    ///
    /// @pre `score` is not null.
    unique_ptr<SiteScore> take_score();
    /// Whether the chosen genotype of the parent, or of an ancestor, leaves this chain out.
    ///
    /// The chain and its descendants are then absent from the sample, and are not revised or
    /// written. Each linkage pass decides this again, so a chain dropped in one pass can come back
    /// in the next.
    bool dropped = false;

    /// Get `lookup.alleles(alleles)`, computing it on the first call.
    ///
    /// The traversals do not change after the direct pass, so the answer is kept, which saves each
    /// re-genotyping round the GBWT lookups.
    const vector<int>& panel_alleles(const PanelLookup& lookup);
    /// The kept answer of `panel_alleles`, once it is computed.
    optional<vector<int>> panel_cache;

    /// The alleles the site's record spells, if they are spliced from the chosen alleles of the
    /// sites under it. The render sets them while it writes the site's records, and then
    /// frees them.
    optional<SplicedAlleles> spliced;
};

/// The sources the passes read for what a staged site does not hold.
///
/// The caller that stages the sites sets them once.
struct SiteReader {
    /// The graph the sites are in.
    const PathPositionHandleGraph* graph = nullptr;
    /// The overlay of sites that the staged sites' net handles are in.
    const EditSnarlDecompositionOverlay* sites = nullptr;
    /// The read-likelihood genotyper that called the sites, or null where another genotyper did.
    const ReadLikelihoodSnarlCaller* genotyper = nullptr;
};

/**
 * The table of every staged site of a run, from the direct pass to the render.
 *
 * The direct pass files each site into a list of the thread that genotyped it, so filing takes no
 * lock. A top-level site, which the linkage pass will not revise, goes to a render queue. A nested
 * site goes to a nested list.
 *
 * The first linkage pass gathers the nested lists into one, and indexes it by parent. Every
 * linkage pass revises that list in place. Once the passes are done, `hand_off` moves each nested
 * site that gets a line of its own into the render queues. The render builds the records of each
 * queue on one thread.
 */
class StagedSiteTable {
public:
    /// Start staging, with one render queue and one nested list per thread.
    ///
    /// Call it before the direct pass, whose threads file sites without a lock.
    void start(size_t threads);

    /// Return true if staging is on.
    bool active() const { return !queues.empty(); }

    /// File a top-level site into the calling thread's render queue.
    void add_top_level(StagedSite&& site);

    /// File a nested site into the calling thread's nested list.
    void add_nested(StagedSite&& site);

    /// Gather the nested lists of the threads into one, and index it by parent if this added any
    /// site.
    ///
    /// Each linkage pass calls it. After the first, there is nothing left to gather.
    ///
    /// @return The gathered list, whose indices stay fixed until `hand_off`.
    vector<StagedSite>& gather_nested();

    /// Get the gathered nested list.
    vector<StagedSite>& nested() { return nested_sites; }

    /// Get the indices in `nested()` of the sites whose parent is `parent_key`, or null when there
    /// are none.
    const vector<size_t>* children_of(size_t parent_key) const;

    /// Return true if any gathered nested site has been indexed.
    bool has_children() const { return !children.empty(); }

    /// Index every staged site by record key.
    ///
    /// Where a key is staged both in a render queue and in the nested list, take the site in the
    /// render queue.
    unordered_map<size_t, StagedSite*> by_key();

    /// Visit each staged parent of a gathered nested site, with the indices of its children.
    ///
    /// Visit parents before their children, by the parent's level, then by its record key.
    ///
    /// @param by_key This table's `by_key()`.
    void for_each_parent_top_down(
        const unordered_map<size_t, StagedSite*>& by_key,
        const function<void(const StagedSite& parent, const vector<size_t>& children)>& visit)
        const;

    /**
     * List the sites a pass reads, in the order every pass visits them.
     *
     * The order is the render queues in turn, then the gathered nested list. Leave out the nested
     * sites that are `dropped` or `placement.reported_inline`, which the hand-off will hold back,
     * so that the list matches what is rendered. Both flags can change between linkage passes, so
     * each call tests them again.
     *
     * @param with_off_reference Include sites with no reference path. They get no line, but they
     * are genotyped, get anchors, and have a meaningful strand.
     */
    vector<StagedSite*> in_order(bool with_off_reference = false);

    /// Visit every staged site, dropped or not, one at a time, in the order of `in_order`.
    void for_each(const function<void(StagedSite&)>& visit);

    /// The counts of sites that `hand_off` held back.
    struct HandOff {
        /// Nested sites an enclosing block's ALT spells.
        size_t inline_unrendered = 0;
        /// Nested sites with no reference path.
        size_t no_ref_unrendered = 0;
    };

    /// Move every nested site that gets a line of its own into the render queues, spread over them
    /// in turn.
    ///
    /// Discard the `dropped` sites. Keep the `placement.reported_inline` sites and those with no
    /// reference path in `unrendered()`. Call it once, after the last linkage pass.
    HandOff hand_off();

    /// Count the render queues.
    size_t queue_count() const { return queues.size(); }

    /// Get render queue `q`.
    vector<StagedSite>& queue(size_t q) { return queues[q]; }

    /// Count the sites in the render queues.
    size_t queued_count() const;
    /// Get the nested sites that `hand_off` kept out of the render queues.
    ///
    /// The records of their parents spell their chosen alleles. Each keeps its walks and genotype,
    /// but not its scores.
    vector<StagedSite>& unrendered() { return held_back; }
    /// Free the sites the hand-off kept, once the render is done.
    void release_unrendered() { vector<StagedSite>().swap(held_back); }

private:
    /// The nested sites that `hand_off` kept out of the render queues.
    vector<StagedSite> held_back;
    /// The render queues, one per thread: the top-level sites, and after the hand-off also the
    /// nested ones that get a line of their own.
    vector<vector<StagedSite>> queues;
    /// The nested sites, one list per thread, as the direct pass files them.
    vector<vector<StagedSite>> nested_lists;
    /// The nested sites, gathered.
    vector<StagedSite> nested_sites;
    /// The indices in `nested_sites` of the children of each parent, by the record key of the
    /// parent.
    unordered_map<size_t, vector<size_t>> children;
};

}

#endif
