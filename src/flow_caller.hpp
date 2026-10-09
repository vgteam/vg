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
#include "staged_site.hpp"
#include "read_strand_table.hpp"
#include "round_history.hpp"
#include "record_renderer.hpp"
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
 * genotypes every site from its own reads and stages it (see StagedSite). Rounds
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

    /// Feed the corrected likelihoods back to the linkage layer and run the linkage pass again.
    ///
    /// The correction changes likelihoods, not genotypes: the linkage model still decides, as in
    /// round 1.
    void rerun_linkage_pass();

    /// Read phasing (see `ReadPhaser`), then rounds of re-genotyping (see `GenotypeRescorer`), each
    /// followed by the linkage pass and read phasing again, as far as they are turned on. The
    /// rounds stop when the correction moves no site's direct call, or when the chosen genotypes
    /// stop changing, return to an earlier round's (see `RoundHistory`), or reach --regeno-passes
    /// rounds.
    void phase_and_regenotype();

    /// Build the records of every staged site once, from its chosen genotype, and collect the
    /// site's anchors (see `RecordRenderer`).
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


    /// The linkage pass over the staged sites (see `GenotypeLinker::link`). Once every level is
    /// done, it decides which chains an enclosing block spells (`StagedSite::reported_inline`)
    /// from the chosen genotypes. Does nothing unless staging is on (see `set_stage_records`).
    void run_linkage_pass();




    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

    /// Genotype through `genotyper`, the read-likelihood genotyper this caller's `SnarlCaller` is,
    /// with explicit ploidies, so that the passes after the direct pass read its typed scores.
    /// Not owned.
    void set_site_genotyper(ReadLikelihoodSnarlCaller& genotyper);

    /// See max_snarl_edges. Zero removes the limit.
    void set_max_snarl_edges(size_t edges) {
        max_snarl_edges = edges ? edges : numeric_limits<size_t>::max();
    }

protected:

    /// Add the record steps this caller needs to `record_steps`: phasing from the linkage model,
    /// the GL layout of the read-likelihood genotyper, block records, and telling the linkage model
    /// each site's allele numbering. Each does nothing when its part is turned off.
    void install_record_steps();

    /// Configure the widgets of the passes with what they read from this caller.
    void install_widgets();

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

    /// Every staged site, while staging is on (see `set_stage_records`). A top-level site's ploidy
    /// comes from the contig or the BED, so the linkage pass never revises it, though the linkage
    /// model still chooses its genotype.
    StagedSiteTable staged_sites;

    /// What the reads say about each site's strands, from read phasing. Re-genotyping and the
    /// anchors read it.
    ReadStrandTable read_strands;

    /// The temper re-genotyping fitted in its first round.
    TemperFit temper_fit;

    /// Builds the staged sites' records once the passes are done.
    RecordRenderer record_renderer;

    /// Lists the child chains nested descent genotypes under each site.
    ChildPlacer child_placer;

    /// The read-likelihood genotyper, or null where the `SnarlCaller` is another (see
    /// `set_site_genotyper`).
    unique_ptr<SiteGenotyper> site_genotyper;

    /// Genotype one site at `ploidies`: through `site_genotyper` where there is one, and otherwise
    /// through the `SnarlCaller`, which takes the ploidy alone. `score` is set to the call info as
    /// the read-likelihood genotyper's score, or to null.
    pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> genotype_site(
        const Snarl& site, const vector<SnarlTraversal>& travs, int ref_trav_idx,
        const Ploidies& ploidies, const string& ref_path_name, pair<size_t, size_t> ref_range,
        SiteScore*& score);

    /// Make a top-level site's `StagedSite` from its genotype, moving `call_info` into it.
    /// `travs` is left empty, because descent still reads the traversals; the caller moves them in
    /// once descent is done. Returns null, and leaves `call_info` alone, when staging is off (see
    /// `set_stage_records`).
    unique_ptr<StagedSite> stage_render_record(const Snarl& snarl,
                                                 const vector<int>& trav_genotype, int ref_trav_idx,
                                                 unique_ptr<SnarlCaller::CallInfo>& call_info,
                                                 SiteScore* score,
                                                 const string& ref_path_name, int ref_offset,
                                                 int ploidy);




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
    /// @param placement Where the snarl sits in the nesting tree: the default for a top-level
    ///                  snarl, and what `ChildPlacer::place` gave a child. A --top-down child
    ///                  takes its parent's.
    bool call_snarl_internal(const Snarl& snarl,
                             const string& parent_ref_path_name,
                             pair<size_t, size_t> parent_ref_interval,
                             const ChildTraversalSets* parent_child_trav_sets,
                             int ploidy_override, const NestingPlacement& placement);


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
