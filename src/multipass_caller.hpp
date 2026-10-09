#ifndef VG_MULTIPASS_CALLER_HPP_INCLUDED
#define VG_MULTIPASS_CALLER_HPP_INCLUDED

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "candidate_finder.hpp"
#include "child_placer.hpp"
#include "graph_caller.hpp"
#include "read_likelihood_caller.hpp"
#include "read_strand_table.hpp"
#include "record_renderer.hpp"
#include "round_history.hpp"
#include "site_genotyper.hpp"
#include "staged_site.hpp"
#include "tree_genotyper.hpp"
#include "vcf_output_caller.hpp"

namespace vg {

using namespace std;

/**
 * MultiPassCaller: the caller for every --read-likelihood run. Its one method, `call`, runs the
 * passes in order:
 *
 * 1. The direct pass genotypes every site from its own reads and stages it (see `TreeGenotyper`
 *    and `StagedSite`).
 * 2. Rounds follow. Each round's linkage pass chooses the genotypes one level at a time, parents
 *    before their children (see `GenotypeLinker`), and read phasing then re-decides the phases
 *    (see `ReadPhaser`). From round 2 on, re-genotyping first corrects the likelihoods from the
 *    phase (see `GenotypeRescorer`). The rounds stop when the correction moves no site's direct
 *    call, when the chosen genotypes stop changing or return to an earlier round's, or at the
 *    round cap.
 * 3. The render builds each staged site's records once, from the genotype the last round chose
 *    (see `RecordRenderer`), and the anchors are written.
 *
 * It is not a `GraphCaller`: it walks the snarls itself, and its genotyper is not a `SnarlCaller`
 * that answers for one site at a time. It writes its records through the `VCFOutputCaller` it is,
 * which also holds the linkage model, read phasing, block emission and the anchor and mosaic
 * files, configured from call_main.
 */
class MultiPassCaller : public VCFOutputCaller {
public:
    /// Call `graph`'s snarls, as `snarl_manager` decomposes it, from the candidate traversals
    /// `traversal_finder` gives, genotyping each with `genotyper`. The reference paths, their
    /// offsets and ploidies, `genotype_snarls` (-a) and `allele_length_range` (-c and -C) are as
    /// for `FlowCaller`. `top_down` and `star_allele` are --top-down and -Y. Nothing is owned.
    MultiPassCaller(const PathPositionHandleGraph& graph, ReadLikelihoodSnarlCaller& genotyper,
                    SnarlManager& snarl_manager, const string& sample_name,
                    TraversalFinder& traversal_finder, const vector<string>& ref_paths,
                    const vector<size_t>& ref_path_offsets, const vector<int>& ref_path_ploidies,
                    bool genotype_snarls, const pair<size_t, size_t>& allele_length_range,
                    bool top_down, bool star_allele);

    virtual ~MultiPassCaller() = default;

    /// Run every pass, and add the records to this caller's VCF buffer. The direct pass visits the
    /// snarls as `GraphCaller::call_top_level_snarls` does with `recurse_type`.
    /// `after_direct_pass` runs once the direct pass is done and before the linkage pass; the
    /// passes after it read only what the direct pass kept, not reads.
    void call(GraphCaller::RecurseType recurse_type, const function<void()>& after_direct_pass);

    /// Do not genotype a snarl with more edges than this, including those of nested snarls
    /// (--max-snarl-edges). The walk then genotypes the snarl's children as if they were top-level
    /// snarls. Zero removes the limit, which is also the default.
    void set_max_snarl_edges(size_t edges) { candidates.set_max_snarl_edges(edges); }

    /// Batch the direct pass's parallel jobs by windows of node IDs, as
    /// `GraphCaller::set_snarl_batching` does.
    void set_snarl_batching(size_t window_size) { snarl_batch_window = window_size; }

    /// Toggle progress messages.
    void set_show_progress(bool show_progress) { this->show_progress = show_progress; }

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

private:

    /// Add the record steps this caller needs to `record_steps`: phasing from the linkage model,
    /// the GL layout of the read-likelihood genotyper, block records, and telling the linkage model
    /// each site's allele numbering. Each does nothing when its part is turned off.
    void install_record_steps();

    /// Configure the widgets of the passes with what they read from this caller.
    void install_widgets();

    /// Report what nested descent did: the depth histogram, and how many children it skipped and
    /// why. Does nothing in a run without nested descent.
    void report_descent_instrumentation() const;

    /// The linkage pass over the staged sites (see `GenotypeLinker::link`). Once every level is
    /// done, it decides which chains an enclosing block spells (`StagedSite::reported_inline`)
    /// from the chosen genotypes.
    void run_linkage_pass();

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

    /// the graph
    const PathPositionHandleGraph& graph;

    /// The read-likelihood genotyper, which also writes its fields into each record.
    ReadLikelihoodSnarlCaller& snarl_caller;

    /// Our snarls
    SnarlManager& snarl_manager;

    /// keep track of the reference paths
    vector<string> ref_paths;
    unordered_set<string> ref_path_set;

    /// keep track of offsets in the reference paths
    map<string, size_t> ref_offsets;

    /// keep track of the ploidies
    map<string, int> ref_ploidies;

    /// toggle whether to genotype every snarl
    /// (by default, uncalled snarls are skipped, and coordinates are flattened
    ///  out to minimize variant size -- this turns all that off)
    bool genotype_snarls;

    /// --top-down: genotype each child against the traversal sets its parent's called alleles
    /// allow (see `TreeGenotyper`).
    bool top_down = false;

    /// use * alleles for spanning haplotypes that don't traverse nested sites
    bool star_allele = false;

    /// See set_snarl_batching.
    size_t snarl_batch_window = 0;

    /// Toggle progress messages
    bool show_progress = false;

    /// Finds each site's reference path and candidate traversals. Declared after the members it
    /// reads.
    CandidateFinder candidates;

    /// Genotypes each site from its reads.
    SiteGenotyper site_genotyper;

    /// See `DescentCounters`.
    mutable DescentCounters descent_counters;

    /// Every staged site. A top-level site's ploidy comes from the contig or the BED, so the
    /// linkage pass never revises it, though the linkage model still chooses its genotype.
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

    /// The direct pass for one top-level site.
    TreeGenotyper tree_genotyper;
};

}

#endif
