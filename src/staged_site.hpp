#ifndef VG_STAGED_SITE_HPP_INCLUDED
#define VG_STAGED_SITE_HPP_INCLUDED

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "panel_lookup.hpp"
#include "site_genotyper.hpp"
#include "site_values.hpp"
#include "snarl_caller.hpp"

namespace vg {

using namespace std;

/// A staged site: one site's genotyping result from the direct pass, kept until the render
/// builds the site's records.
///
/// A nested chain's ploidy depends on its parent's chosen genotype, which is known only after
/// the direct pass, so the result is kept rather than computed again. The `CallInfo` has the
/// answer at both ploidies (see `alt_ploidy_info`), so the record can be rendered at whichever
/// ploidy the linkage pass gives the site.
///
/// The site is named by its record key alone. The fields that describe it in the graph, its
/// bounds, its candidate walks and what it holds, are grouped apart from the rest.
struct StagedSite {
    // Identity, and the site's place in the nesting tree.

    /// See `VCFOutputCaller::record_key_of`.
    size_t record_key = 0;
    /// The record key of the site this one is nested in, or 0 for a top-level site.
    size_t parent_record_key = 0;
    /// See `NestingPlacement::chain_key`.
    size_t chain_key = 0;
    /// The site's level.
    uint8_t level = 0;
    /// See `NestingPlacement::parent_crossing`.
    uint64_t parent_crossing = 0;
    /// False when the parent has more than 64 candidate traversals, too many for
    /// `parent_crossing`. A 0 mask then means unknown, and the linkage pass leaves the chain at the
    /// ploidy the direct pass gave it. The linkage pass computes the mask again when it revises
    /// or first records the parent.
    bool crossing_known = true;

    // The site in the graph: its bounds, its candidate walks, its child chains, and the chain it
    // is in.

    /// Oriented forward along the reference path, which may be the decomposition's orientation
    /// turned round.
    SiteBounds bounds;
    vector<Traversal> travs;
    int ref_trav_idx = -1;
    /// The sites nested in this one, for its symbolic alleles.
    SiteChildren children;
    /// Whether the site has no child sites.
    bool leaf = true;
    /// The chain the site is in, when the decomposition knows the site (`in_chain`). A nested
    /// site's line is held back when its parent's blocks report its chain.
    ChildChain chain;
    bool in_chain = false;

    // Where the site is on the reference.

    string ref_path_name;
    int ref_offset = 0;
    /// This snarl has no reference path, so no line can be written for it, since REF and POS are
    /// undefined. It is still genotyped and recorded in the linkage model.
    bool no_reference = false;
    /// For a snarl with no reference path: its parent's reference start plus `chain_offset`,
    /// standing in for the position it lacks.
    int64_t position_from_parent = 0;
    /// `NestingPlacement::parent_offset` for this chain. The direct pass takes it from the parent's
    /// direct call, and the linkage pass computes it again from the parent's chosen genotype, so
    /// that an off-reference chain is placed along an allele of its parent's chosen genotype.
    size_t chain_offset = 0;

    // The direct pass's call, and what the linkage pass decided about it.

    /// The direct pass's genotype, before the linkage model, and its ploidy. A nested chain's
    /// ploidy here is the number of the parent's called alleles that cross it, or the parent's
    /// ploidy when none does; `score` also holds the answer at the other ploidy.
    vector<int> genotype;
    int ploidy = 2;
    /// The genotyper's call, which the record is written from.
    unique_ptr<SnarlCaller::CallInfo> call_info;
    /// `call_info` as the read-likelihood genotyper's score, or null where another genotyper
    /// called the site. `call_info` owns it, so the two are set together, by `set_call` or
    /// `set_score`.
    SiteScore* score = nullptr;

    /// Give the site the genotyper's call `info`, and `typed`, which is null or `info` as the
    /// read-likelihood genotyper's score.
    void set_call(unique_ptr<SnarlCaller::CallInfo> info, SiteScore* typed);
    /// Give the site a read-likelihood score as its call.
    void set_score(unique_ptr<SiteScore> typed);
    /// Take the site's read-likelihood score out, leaving it no call. `score` must not be null.
    unique_ptr<SiteScore> take_score();
    /// See `NestingPlacement::reported_inline`. Its line is held back, as for
    /// `no_reference`. The direct pass tests the parent's direct call; under the linkage model the
    /// linkage pass tests again with the parent's chosen genotype, which the parent's blocks are
    /// built from.
    bool reported_inline = false;
    /// Set when the parent's chosen genotype, or an ancestor's, does not carry this chain, so
    /// the chain and its descendants do not exist in the sample and are not revised or
    /// written. Each linkage pass decides it again, so a chain dropped in one pass can come
    /// back in the next.
    bool dropped = false;

    /// `lookup.alleles(travs)`, computed on the first call and kept: the traversals do not change
    /// after the direct pass, and each re-genotyping round would otherwise repeat the GBWT
    /// lookups. The direct pass can fill `panel_cache` first.
    const vector<int>& panel_alleles(const PanelLookup& lookup);
    vector<int> panel_cache;
    bool panel_cached = false;
};

/// Where the passes read what a staged site does not hold, set once by the caller that stages
/// the sites.
struct SiteReader {
    /// The graph the sites are in.
    const PathPositionHandleGraph* graph = nullptr;
    /// The read-likelihood genotyper that called the sites, or null where another genotyper did.
    const SiteGenotyper* genotyper = nullptr;
    /// A walk's sequence.
    function<string(const Traversal&)> spell;
    /// A site's ID, as its records name it.
    function<string(const SiteBounds&)> name;
};

/**
 * Every staged site of a run, from the direct pass to the render.
 *
 * The direct pass files each site from the thread that genotyped it, into that thread's own
 * list, so filing takes no lock. A top-level site, which the linkage pass will not revise, goes
 * to a render queue; a nested site goes to a nested list. The first linkage pass gathers the
 * nested lists into one, which every linkage pass then revises in place, and indexes it by
 * parent. Once the passes are done, `hand_off` moves each nested site that gets a line of its
 * own into the render queues, and the render builds each queue's records on one thread.
 */
class StagedSiteTable {
public:
    /// Start staging, with one render queue and one nested list per thread. Called before the
    /// direct pass, since its threads file into these without a lock.
    void start(size_t threads);

    /// Whether staging is on.
    bool active() const { return !queues.empty(); }

    /// File a top-level site into the calling thread's render queue.
    void add_top_level(StagedSite&& site);

    /// File a nested site into the calling thread's nested list.
    void add_nested(StagedSite&& site);

    /// Gather the threads' nested lists into one, and index it by parent if this added any site.
    /// Returns the gathered list, whose indices stay fixed until `hand_off`. Each linkage pass
    /// calls it; after the first, there is nothing left to gather.
    vector<StagedSite>& gather_nested();

    /// The gathered nested list.
    vector<StagedSite>& nested() { return nested_sites; }

    /// The indices in `nested()` of the sites whose parent is `parent_key`, or null when there
    /// are none.
    const vector<size_t>* children_of(size_t parent_key) const;

    /// Whether any gathered nested site has been indexed.
    bool has_children() const { return !children.empty(); }

    /// Every staged site by record key. Where a key is staged both in a render queue and in the
    /// nested list, the render queue's site wins.
    unordered_map<size_t, StagedSite*> by_key();

    /// Visit each staged parent of a gathered nested site, with the indices of its children,
    /// parents before their children: by the parent's level, then by its record key. `by_key`
    /// is this table's `by_key()`.
    void for_each_parent_top_down(
        const unordered_map<size_t, StagedSite*>& by_key,
        const function<void(const StagedSite& parent, const vector<size_t>& children)>& visit)
        const;

    /// The sites a pass reads, in the order every pass visits them: the render queues in turn,
    /// then the gathered nested list. Between a linkage pass and the hand-off it leaves out the
    /// nested sites the hand-off will hold back, so that it matches what is rendered: a dropped
    /// one, which the parent's chosen genotype does not carry, and a `reported_inline` one,
    /// which an enclosing block's ALT spells. Both are tested on each call, since both can
    /// change between linkage passes. With `with_off_reference`, sites with no reference path
    /// are included: they cannot be rendered, having no REF or POS, but they are genotyped, get
    /// anchors, and have a meaningful strand.
    vector<StagedSite*> in_order(bool with_off_reference = false);

    /// Visit every staged site, dropped or not, one at a time, in the order of `in_order`.
    void for_each(const function<void(StagedSite&)>& visit);

    /// What `hand_off` held back.
    struct HandOff {
        /// Nested sites an enclosing block's ALT spells.
        size_t inline_unrendered = 0;
        /// Nested sites with no reference path.
        size_t no_ref_unrendered = 0;
    };

    /// Move every nested site that gets a line of its own into the render queues, spread over
    /// them in turn, and discard the rest: a dropped site, which the sample does not carry, a
    /// `reported_inline` one and one with no reference path. Called once, after the last
    /// linkage pass.
    HandOff hand_off();

    /// How many render queues there are.
    size_t queue_count() const { return queues.size(); }

    /// Render queue `q`.
    vector<StagedSite>& queue(size_t q) { return queues[q]; }

    /// How many sites the render queues hold.
    size_t queued_count() const;

private:
    /// The render queues, one per thread: the top-level sites, and after the hand-off also the
    /// nested ones that get a line of their own.
    vector<vector<StagedSite>> queues;
    /// The nested sites, one list per thread, as the direct pass files them.
    vector<vector<StagedSite>> nested_lists;
    /// The nested sites, gathered.
    vector<StagedSite> nested_sites;
    /// Parent record key -> the indices in `nested_sites` of its children.
    unordered_map<size_t, vector<size_t>> children;
};

}

#endif
