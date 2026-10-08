#ifndef VG_FLOW_CALLER_HPP_INCLUDED
#define VG_FLOW_CALLER_HPP_INCLUDED

#include <atomic>
#include <iostream>
#include <algorithm>
#include <array>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include <gbwt/cached_gbwt.h>
#include "handle.hpp"
#include "linkage_model.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "anchor.hpp"
#include "read_phasing.hpp"
#include "regenotype.hpp"
#include "snarl_caller.hpp"
#include "symbolic_allele.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"
#include "vcf_genotype_likelihoods.hpp"
#include "graph_caller.hpp"
#include "vcf_output_caller.hpp"
#include "gaf_output_caller.hpp"

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/// A set of traversals through a child snarl that are consistent with
/// a single parent allele. Multiple traversals can exist if the child
/// has internal variation within a shared region.
using TraversalSet = vector<SnarlTraversal>;

/// One TraversalSet per parent allele (index matches parent genotype).
/// For a diploid parent with genotype [0,1], element 0 contains traversals
/// consistent with parent allele 0, element 1 with parent allele 1.
using ChildTraversalSets = vector<TraversalSet>;

/// Counts of what nested descent did in one run: how deep it went, and how many child chains it
/// skipped or recorded, and why. A member of each FlowCaller, so that runs count separately;
/// `mutable` there because the counting paths are const.
struct DescentCounters {
    /// How deep the descent went, by depth.
    std::atomic<size_t> depth_hist[16] = {};
    /// Children skipped because no reference path passes through them and off-reference descent
    /// is off.
    std::atomic<size_t> skipped_no_ref{0};
    /// Children that no reference path passes through, descended into because off-reference
    /// descent is on.
    std::atomic<size_t> off_reference{0};
    /// Sites that no reference path passes through, given an entry in the linkage model but no
    /// line.
    std::atomic<size_t> no_ref_recorded{0};
    /// Off-reference chains by copy number: 0, 1, 2.
    std::atomic<size_t> no_ref_copies[3] = {};
    /// Children that no called parent allele crosses in the direct pass. They are genotyped and kept,
    /// and the linkage pass decides from the parent's chosen genotype whether the sample has them.
    std::atomic<size_t> skipped_no_copy{0};
    /// Children that a called traversal enters more than once. Only the first entry counts, both
    /// for ploidy and for distance: one traversal crossing a chain twice is one strand carrying two
    /// copies, not two strands, so it does not make the chain ploidy 2.
    std::atomic<size_t> child_multi_crossing{0};
};

/**
 * FlowCaller: takes each snarl's candidate traversals from a TraversalFinder and genotypes them
 * with its SnarlCaller (support-based, or ReadLikelihoodSnarlCaller under --read-likelihood). It
 * works on any graph. With the flow traversal finder it does not report cyclic traversals;
 * haplotype enumeration can. The `nested` constructor flag is --top-down, which genotypes each
 * child against traversal sets derived from its parent's called alleles. Nested calling, the
 * descent into child chains in call_snarl_internal, is turned on instead by
 * set_symbolic_collapsing.
 *
 * With the linkage model or nested calling, calling runs in passes. The direct pass
 * genotypes every site from its own reads and stages it (see PendingRecord). Rounds
 * follow. Each round's linkage pass chooses the genotypes one level at a time, parents
 * before their children (see run_linkage_pass), and read phasing then re-decides the
 * phases. From round 2 on, re-genotyping first corrects the likelihoods from the phase
 * (see phase_and_regenotype). Finally the render builds each staged site's records once,
 * from its settled genotype, the one the last round chose (see render_retained_records).
 */
class FlowCaller : public GraphCaller, public VCFOutputCaller, public GAFOutputCaller {
public:
    /// Original constructor for non-nested mode
    FlowCaller(const PathPositionHandleGraph& graph,
               SupportBasedSnarlCaller& snarl_caller,
               SnarlManager& snarl_manager,
               const string& sample_name,
               TraversalFinder& traversal_finder,
               const vector<string>& ref_paths,
               const vector<size_t>& ref_path_offsets,
               const vector<int>& ref_path_ploidies,
               AlignmentEmitter* aln_emitter,
               bool traversals_only,
               bool gaf_output,
               size_t trav_padding,
               bool genotype_snarls,
               const pair<size_t, size_t>& allele_length_range);

    /// Extended constructor for nested mode with star alleles
    FlowCaller(const PathPositionHandleGraph& graph,
               SupportBasedSnarlCaller& snarl_caller,
               SnarlManager& snarl_manager,
               const string& sample_name,
               TraversalFinder& traversal_finder,
               const vector<string>& ref_paths,
               const vector<size_t>& ref_path_offsets,
               const vector<int>& ref_path_ploidies,
               AlignmentEmitter* aln_emitter,
               bool traversals_only,
               bool gaf_output,
               size_t trav_padding,
               bool genotype_snarls,
               const pair<size_t, size_t>& allele_length_range,
               bool nested,
               bool star_allele);

    virtual ~FlowCaller();

    /// GraphCaller::call_top_level_snarls, followed, when progress messages are on, by a report
    /// of what nested descent did.
    virtual void call_top_level_snarls(const HandleGraph& graph,
                                       RecurseType recurse_type = RecurseOnFail);

    virtual bool call_snarl(const Snarl& snarl);

    /// A site's place on the reference: a contig, and an offset along it. Positions on one contig
    /// order its sites, and the difference of two positions is the distance in bases between the
    /// sites, except where a position is a stand-in (see `position`). The linkage model relies on
    /// both: it builds its top-level linkage chains from the sites of one contig, orders a chain's
    /// sites by position and measures the gaps between them from it, and names a phase set
    /// (FORMAT/PS) by the position of the chain's first site.
    struct SiteLocus {
        /// The contig as the VCF names it: the locus part of a PanSN path name, so `chr20` for
        /// `CHM13#0#chr20`, or the path name itself when it is not a PanSN name.
        string contig;
        /// A 0-based offset along the contig. For a site the reference path passes through, it is
        /// where the site's first boundary node starts on that path. It must not depend on which
        /// alleles the site's records carry, so it is not their POS, which trimming the alleles
        /// can move. For a site that no reference path passes through, it is a stand-in (see
        /// off_reference_site_locus), and its difference from another position is a distance only
        /// for a site of the same chain, or for the chain's parent.
        size_t position = 0;
    };

    /// The locus of a site that the reference path `ref_path_name` passes through. `ref_offset`
    /// is added to every position along that path, to place the path on its contig.
    SiteLocus site_locus(const Snarl& snarl, const string& ref_path_name, int ref_offset) const;

    /// The locus of a site that no reference path passes through, and that so has no position of
    /// its own. `ref_path_name` is the reference path through the site's nearest ancestor on a
    /// reference path, and `stand_in_position` is that ancestor's position plus how far along the
    /// ancestor's allele the site's chain starts (`PendingRecord::position_from_parent`).
    SiteLocus off_reference_site_locus(const string& ref_path_name,
                                       int64_t stand_in_position) const;

    /// Record the site in the linkage model when it is genotyped, rather than when its line is
    /// written, since the linkage pass reads the collector before any line is written. The emitted
    /// allele map and whether a line was written are supplied later, by `set_allele_map`. The direct pass
    /// does not call this for a retained chain with a reference path (see
    /// `NestedContext::retain_only`); the linkage pass records that chain if the sample carries it.
    ///
    /// `ref_path_name` and `ref_offset` give the site's locus as for `site_locus`. `no_reference`
    /// marks a site that no reference path passes through: `ref_path_name` is then the reference
    /// path through the site's nearest ancestor on a reference path, `position_from_parent` is the
    /// site's stand-in position (see `off_reference_site_locus`), and `ref_offset` is not used.
    /// Otherwise `position_from_parent` is not used.
    ///
    /// Returns whether the site was recorded. If it was and `panel_out` is given, the panel
    /// alleles looked up for it, `panel_lookup.alleles(travs)`, are moved to `panel_out`, so that
    /// the staged record can keep them rather than look them up again.
    bool record_site(const Snarl& snarl, const vector<SnarlTraversal>& travs,
                     const vector<int>& trav_genotype,
                     const unique_ptr<SnarlCaller::CallInfo>& call_info, int ref_trav_idx,
                     const string& ref_path_name, int ref_offset,
                     bool no_reference = false, int64_t position_from_parent = 0,
                     vector<int>* panel_out = nullptr);

    /// The frequency exponent a site should decode with: `--hp-prior` at a run-length site, or -1
    /// for the model's own. Reads the traversals' sequences only when `--hp-prior` is on.
    double site_freq_prior(const vector<SnarlTraversal>& travs, int ref_trav_idx) const;

    /// Decide every heterozygous site's phase from the reads, and change the chosen phase to
    /// match, so that the GT order, the anchor slot column and the mosaic all follow from it.
    /// Genotypes are not changed. On FlowCaller because it needs the staged sites, which hold
    /// the per-read evidence.
    void apply_read_phasing();

    /// Re-score every retained site's genotype likelihoods with the reads' phase.
    ///
    /// Runs after `apply_read_phasing`, which supplies `phase_sites` and `phase_flips`. With
    /// --regeno-passes above 1 the corrected likelihoods replace each site's own, and GQ is
    /// recomputed from them where the best genotype changed and is the direct pass's elsewhere. Returns
    /// true if any site's corrected best genotype differs from its called one.
    bool apply_regenotyping();

    /// Feed the corrected likelihoods back to the linkage layer and run the linkage pass again.
    ///
    /// The correction changes likelihoods, not genotypes: the linkage model still decides, as in
    /// round 1.
    void rerun_linkage_pass();

    /// Read phasing, then rounds of re-genotyping, each followed by the linkage pass and read phasing
    /// again, as far as they are turned on. Does nothing unless read phasing is on.
    void phase_and_regenotype();

    /// Build the records of every staged site once, from its chosen genotype, and collect the
    /// site's anchors.
    void render_retained_records();

    /// Whether this snarl has no children, resolved through the manager's own copy. See the
    /// implementation for why the obvious `children_of(&snarl)` is not safe here.
    bool snarl_is_leaf(const Snarl& snarl) const;

    /// Stage every site during the direct pass, and write its records only after the linkage pass
    /// has chosen its genotype. A nested chain's ploidy, its strand, and whether it has a record at
    /// all then come from its parent's chosen genotype, and a parent is chosen before its children,
    /// so a child's evidence cannot change its parent. Sizes the per-thread queues, so it must be
    /// called before calling starts.
    void set_stage_records(bool defer);

    /// How many staged sites `render_records` holds, reported under --progress.
    size_t render_record_count() const;


    /// The linkage pass: choose the genotypes one level at a time. The linkage model chooses
    /// each level's genotypes; then each child chain of the next level takes the ploidy its
    /// parent's chosen genotype gives it, from the answers the direct pass kept at both
    /// ploidies, and a chain the parent does not carry is dropped with everything inside it.
    /// Once every level is done, it decides which chains an enclosing block spells
    /// (`PendingRecord::reported_inline`). Does nothing unless staging is on (see
    /// `set_stage_records`).
    void run_linkage_pass();

    /// Move every nested chain the linkage pass kept into the render's queues, and collect anchors for
    /// those that get no line. Separate from `run_linkage_pass`, which re-genotyping runs again,
    /// because moving the records and collecting their anchors must happen once.
    void hand_off_deferred_records();



    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

    /// See max_snarl_edges. Zero removes the limit.
    void set_max_snarl_edges(size_t edges) {
        max_snarl_edges = edges ? edges : numeric_limits<size_t>::max();
    }

protected:

    /// Add the record steps this caller needs to `record_steps`: phasing from the linkage model,
    /// the GL layout of the read-likelihood genotyper, block records, and telling the linkage model
    /// each site's allele numbering. Each does nothing when its part is turned off.
    void install_record_steps();

    /// Report what nested descent did: the depth histogram, and how many children it skipped and
    /// why. Does nothing in a run without nested descent.
    void report_descent_instrumentation() const;

    /// See `DescentCounters`.
    mutable DescentCounters descent_counters;

    /// the graph
    const PathPositionHandleGraph& graph;

    /// the traversal finder
    TraversalFinder& traversal_finder;

    /// keep track of the reference paths
    vector<string> ref_paths;
    unordered_set<string> ref_path_set;

    /// keep track of offsets in the reference paths
    map<string, size_t> ref_offsets;
    
    /// keep traco of the ploidies (todo: just one map for all path stuff!!)
    map<string, int> ref_ploidies;

    /// Do not genotype a snarl with more edges than this, including those of nested snarls
    /// (--max-snarl-edges). `call_top_level_snarls` then genotypes the snarl's children as if they
    /// were top-level snarls. No limit until `vg call` sets one.
    size_t max_snarl_edges = numeric_limits<size_t>::max();

    /// alignment emitter. if not null, traversals will be output here and
    /// no genotyping will be done
    AlignmentEmitter* alignment_emitter;

    /// toggle whether to genotype or just output the traversals
    bool traversals_only;

    /// toggle whether to output vcf or gaf
    bool gaf_output;

    /// toggle whether to genotype every snarl
    /// (by default, uncalled snarls are skipped, and coordinates are flattened
    ///  out to minimize variant size -- this turns all that off)
    bool genotype_snarls;

    /// clamp calling to alleles of a given length range
    /// more specifically, a snarl is only called if
    /// 1) its largest allele is >= allele_length_range.first and
    /// 2) all alleles are < allele_length_range.second
    pair<size_t, size_t> allele_length_range;

    /// --- Nested mode members ---

    /// --top-down: after a snarl is called, genotype each child against the traversal sets its
    /// called alleles allow (see call_snarl_internal).
    bool nested = false;

    /// use * alleles for spanning haplotypes that don't traverse nested sites
    bool star_allele = false;

    /// A staged site: one site's genotyping result from the direct pass, kept until the render
    /// builds the site's records.
    ///
    /// A nested chain's ploidy depends on its parent's chosen genotype, which is known only after
    /// the direct pass, so the result is kept rather than computed again. The `CallInfo` has the
    /// answer at both ploidies (see `alt_ploidy_info`), so the record can be rendered at whichever
    /// ploidy the linkage pass gives the site. `snarl` is held by value because
    /// `call_snarl_internal` may work on a flipped copy.
    struct PendingRecord {
        Snarl snarl;
        string ref_path_name;
        int ref_offset = 0;
        vector<SnarlTraversal> travs;
        int ref_trav_idx = -1;
        /// The direct pass's genotype, before the linkage model, and its ploidy. A nested chain's
        /// ploidy here is the number of the parent's called alleles that cross it, or the parent's
        /// ploidy when none does; `call_info` also holds the answer at the other ploidy.
        vector<int> genotype;
        int ploidy = 2;
        unique_ptr<SnarlCaller::CallInfo> call_info;
        size_t record_key = 0;
        size_t parent_record_key = 0;
        /// See NestedContext::chain_key.
        size_t chain_key = 0;
        /// See NestedContext::reported_inline. Its line is held back, as for `no_reference`. The
        /// direct pass tests the parent's direct call; under the linkage model the linkage pass
        /// tests again with the parent's chosen genotype, which the parent's blocks are built from.
        bool reported_inline = false;
        /// This snarl has no reference path, so no line can be written for it, since REF and POS are
        /// undefined. It is still genotyped and recorded in the linkage model.
        bool no_reference = false;
        /// For a snarl with no reference path: its parent's reference start plus `chain_offset`,
        /// standing in for the position it lacks.
        int64_t position_from_parent = 0;
        /// `NestedContext::parent_offset` for this chain. The direct pass takes it from the parent's
        /// direct call, and the linkage pass computes it again from the parent's chosen genotype, so
        /// that an off-reference chain is placed along an allele of its parent's chosen genotype.
        size_t chain_offset = 0;
        /// See NestedContext::parent_crossing.
        uint64_t parent_crossing = 0;
        /// False when the parent has more than 64 candidate traversals, too many for
        /// `parent_crossing`. A 0 mask then means unknown, and the linkage pass leaves the chain at the
        /// ploidy the direct pass gave it. The linkage pass computes the mask again when it revises
        /// or first records the parent.
        bool crossing_known = true;
        /// The site's level.
        uint8_t level = 0;
        /// Set when the parent's chosen genotype, or an ancestor's, does not carry this chain, so
        /// the chain and its descendants do not exist in the sample and are not revised or
        /// written. Each linkage pass decides it again, so a chain dropped in one pass can come
        /// back in the next.
        bool dropped = false;
        /// `panel_lookup.alleles(travs)`, computed once: the traversals do not change after the
        /// direct pass, and each re-genotyping round would otherwise repeat the GBWT lookups.
        vector<int> panel_cache;
        bool panel_cached = false;

    };

    /// The staged sites a pass should look at, wherever they currently are: between a linkage
    /// pass and the hand-off, nested chains are in `deferred_pending` and the rest in
    /// `render_records`.
    ///
    /// With `for_phasing`, chains with no reference path are included: they cannot be rendered,
    /// having no REF or POS, but they are genotyped, get anchors, and have a meaningful strand.
    vector<PendingRecord*> records_for_render(bool for_phasing = false);

    /// `panel_lookup.alleles` for a record, computed once and kept. See
    /// `PendingRecord::panel_cache`.
    const vector<int>& cached_panel_alleles(PendingRecord& rec);

    /// The chosen pair and ploidy per record, `{trav_first, trav_second, ploidy}`, for measuring
    /// whether a re-genotyping round changed anything.
    unordered_map<size_t, std::array<int, 3>> chosen_snapshot();
    /// How many records have a different chosen pair or ploidy in snapshot `after` than in
    /// snapshot `before`, counting a chain that gained or lost a chosen answer as moved.
    static size_t chosen_changed(const unordered_map<size_t, std::array<int, 3>>& before,
                                 const unordered_map<size_t, std::array<int, 3>>& after);
    /// A digest of a snapshot that does not depend on order, for spotting a state the rounds have
    /// reached before, which means they are cycling.
    static size_t snapshot_digest(const unordered_map<size_t, std::array<int, 3>>& snap);

    /// The nested chains' staged sites, which every linkage pass reads, merged out of
    /// `pending_records` by the first and kept until `hand_off_deferred_records` moves them to the
    /// render. A member because the linkage pass runs once per round.
    vector<PendingRecord> deferred_pending;
    /// How many times the linkage pass has run.
    size_t linkage_passes_run = 0;


    /// See set_stage_records.
    bool stage_records = false;

    /// The nested chains' staged sites, filled per thread during the direct pass.
    vector<vector<PendingRecord>> pending_records;

    /// The top-level staged sites, and after the hand-off also the nested ones that get a line of
    /// their own.
    /// A top-level site's ploidy comes from the contig or the BED, so the linkage pass never revises
    /// it, though the linkage model still chooses its genotype. Separate from
    /// `pending_records`, which `run_linkage_pass` moves out and clears, and whose index groups
    /// records by parent. Read through `records_for_render`.
    vector<vector<PendingRecord>> render_records;


    /// Make a top-level site's `PendingRecord` from its genotype, moving `call_info` into it.
    /// `travs` is left empty, because descent still reads the traversals; the caller moves them in
    /// once descent is done. Returns null, and leaves `call_info` alone, when staging is off (see
    /// `set_stage_records`).
    unique_ptr<PendingRecord> stage_render_record(const Snarl& snarl,
                                                 const vector<int>& trav_genotype, int ref_trav_idx,
                                                 unique_ptr<SnarlCaller::CallInfo>& call_info,
                                                 const string& ref_path_name, int ref_offset,
                                                 int ploidy);


    /// The genotype the linkage model chose for a staged site, or the direct pass's genotype
    /// where the model chose none.
    vector<int> chosen_genotype_for(const PendingRecord& rec) const;

    /// The gqn column's value for this record: the direct pass's `gq_fraction`, unless the linkage model
    /// changed the call, in which case the signed value recomputed for the chosen genotype. NaN,
    /// written as `.`, where there is no value: no gap to normalise, or a moved call whose margin
    /// cannot be recomputed.
    double anchor_gqn_for(const PendingRecord& rec, const vector<int>& chosen) const;

    /// Turn a staged site into anchors, if anchors are being written, with the phase order, the
    /// haploid slot, the leaf test and the gqn derived from it. Called once per staged site as the
    /// sites are rendered: just before its line is written, or, for a site with no line, by the
    /// hand-off. A site with no reference position still gets anchors, since a pin is placed by
    /// node ID, which is why this is not part of `emit_variant`. The genotype is a parameter
    /// because the render passes the chosen pair.
    void collect_anchors_for_record(const PendingRecord& rec, const vector<int>& genotype);

    /// How many nested chains are staged, over all threads.
    size_t pending_record_count() const;


    /// Internal implementation of call_snarl that accepts parent context for nested mode
    /// When nested=true, this recursively calls children after processing the current snarl
    /// @param parent_ref_path_name Reference path from parent (for off-reference snarls)
    /// @param parent_ref_interval Reference interval from parent
    /// @param parent_child_trav_sets If non-null, contains one TraversalSet per parent allele.
    ///                               Each set contains all traversals through this child that are
    ///                               consistent with that parent allele, and the child's genotype
    ///                               takes one allele from each set. --top-down (FlowCaller's
    ///                               `nested` mode) passes them; nested calling (the descent that
    ///                               set_symbolic_collapsing turns on) and -A pass null.
    /// @param ploidy_override If >= 0, the ploidy to genotype this snarl at, instead of the
    ///                        contig's or the --ploidy-bed region's. Nested calling passes the
    ///                        number of the parent's called alleles that cross the child, or the
    ///                        parent's ploidy when none does.
    bool call_snarl_internal(const Snarl& snarl,
                             const string& parent_ref_path_name,
                             pair<size_t, size_t> parent_ref_interval,
                             const ChildTraversalSets* parent_child_trav_sets = nullptr,
                             int ploidy_override = -1);

    /// Where each node ID is visited in one traversal, in ascending order, excluding visits to
    /// snarls. Built once per traversal per snarl, so that testing each child does not scan the
    /// whole traversal again.
    using TraversalNodeIndex = unordered_map<nid_t, vector<int>>;
    static TraversalNodeIndex index_traversal_nodes(const SnarlTraversal& trav);

    /// How many of the called parent alleles cross this child snarl, capped at `cap`.
    ///
    /// A traversal crosses the child when the child's start and end both appear in it in order, so
    /// a traversal that touches both boundaries on unrelated excursions does not count. One
    /// allele crossing a chain more than once, as in a cycle or tandem duplication, counts once,
    /// since the caller assumes ploidy 1 or 2; this is logged.
    int child_ploidy(const vector<TraversalNodeIndex>& visits, const vector<int>& genotype,
                     const Snarl& child, int cap) const;

    /// How many times one traversal crosses `child`, by the same in-order rule child_ploidy uses.
    static int crossings_of_child(const TraversalNodeIndex& visits, const Snarl& child);

public:
    /// Where along `trav` the child chain is first entered, as a visit index, or -1 if `trav` does
    /// not cross it, by the rule `crossings_of_child` uses.
    static int offset_of_child(const SnarlTraversal& trav, const Snarl& child);

    /// A nested site's place in the nesting tree: its record key, its parent's, and its level.
    struct NestedLink {
        size_t key = 0;
        size_t parent = 0;
        uint8_t level = 0;
    };

    /// Keep each nested site's strand pointing at the parent strand that carries it, after the
    /// sites in `flips` had their chosen pair swapped. A ploidy-1 nested site names one of its
    /// parent's two strands in `nested_strand` and holds its haplotype in the slot of that number,
    /// so where the parent's strands swapped, both move to the other strand. Each site in `links`
    /// is looked up in `phased` through `phase_index`. Returns how many strands moved.
    static size_t cascade_nested_strands(vector<LinkageCollector::PhaseCall>& phased,
                                         const std::unordered_map<size_t, size_t>& phase_index,
                                         vector<NestedLink> links,
                                         const unordered_set<size_t>& flips);

    /// How far along `trav`, in bases, the child chain is entered: the total length of the nodes
    /// visited before it, or -1 if `trav` does not cross it. It gives an off-reference chain its
    /// place along its parent (see `NestedContext::parent_offset`).
    int64_t base_offset_of_child(const SnarlTraversal& trav, const Snarl& child) const;

    /// `base_offset_of_child` along the first traversal of `genotype` that crosses `child`, or 0
    /// when none does.
    size_t offset_along_genotype(const vector<SnarlTraversal>& travs, const vector<int>& genotype,
                                 const Snarl& child) const;

    /// `base_offset_of_child` for every child of one traversal, by lookup. Calling
    /// `base_offset_of_child` once per child scans the traversal once per child, which a parent
    /// with many children and a long traversal makes quadratic.
    struct ChildOffsets {
        ChildOffsets(const HandleGraph& graph, const SnarlTraversal& trav);
        /// The same answer as `base_offset_of_child(trav, child)`.
        int64_t base_offset(const Snarl& child) const;
        /// The visit indices of each node, ascending. Child-snarl visits are left out, as
        /// `offset_of_child` skips them.
        unordered_map<nid_t, vector<int>> visits_of;
        /// The bases of the node visits before each visit index; one longer than the traversal.
        vector<int64_t> bases_before;
    };

    /// `offset_along_genotype`, answered from `offsets`, which holds a `ChildOffsets` per
    /// traversal and is filled as traversals are first used.
    size_t offset_along_genotype(const vector<SnarlTraversal>& travs, const vector<int>& genotype,
                                 const Snarl& child,
                                 unordered_map<const SnarlTraversal*, ChildOffsets>& offsets) const;
protected:

    /// The crossing mask: bit i is set where `travs[i]` crosses `child`. Indexed by traversal, not
    /// by VCF allele, since it is tested against the parent's chosen traversals. Returns 0 and
    /// sets `*known` to false when there are more than 64 traversals, so that the caller can tell
    /// unknown from "no traversal crosses".
    static uint64_t child_crossing_mask(const vector<TraversalNodeIndex>& visits,
                                        const Snarl& child, bool* known = nullptr);

    /// Find all traversals through a child snarl that are consistent with a parent traversal.
    /// "Consistent" means the child's entry/exit points match what's in the parent traversal.
    /// Uses the traversal finder to enumerate all valid paths through the child.
    /// @param parent_trav The parent traversal defining entry/exit constraints
    /// @param child The child snarl to find traversals through
    /// @return Set of traversals through child, empty if parent doesn't traverse child
    TraversalSet find_child_traversal_set(const SnarlTraversal& parent_trav,
                                          const Snarl& child) const;

    /// Extract the portion of a parent traversal that spans a child snarl (single traversal).
    /// This is a simpler version used when we only need one traversal from the parent.
};

}

#endif
