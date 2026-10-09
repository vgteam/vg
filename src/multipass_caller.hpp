#ifndef VG_MULTIPASS_CALLER_HPP_INCLUDED
#define VG_MULTIPASS_CALLER_HPP_INCLUDED

#include <atomic>
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
#include "anchor.hpp"
#include "block_records.hpp"
#include "candidate_finder.hpp"
#include "child_placer.hpp"
#include "genotype_linker.hpp"
#include "genotype_rescorer.hpp"
#include "graph_caller.hpp"
#include "linkage_model.hpp"
#include "mosaic_writer.hpp"
#include "panel_lookup.hpp"
#include "phase_table.hpp"
#include "ploidy_regions.hpp"
#include "read_likelihood_caller.hpp"
#include "read_phaser.hpp"
#include "read_strand_table.hpp"
#include "record_renderer.hpp"
#include "round_history.hpp"
#include "site_genotyper.hpp"
#include "site_walker.hpp"
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
 * that answers for one site at a time. It writes its records through a `VCFOutputCaller` it holds,
 * to which it adds its own header lines and record and writing steps. It holds the linkage model,
 * read phasing, block emission and the anchor and mosaic files, configured from call_main.
 */
class MultiPassCaller {
public:
    /// Call `graph`'s sites, as `decomposition` gives them, from the candidate traversals
    /// `traversal_finder` gives, genotyping each with `genotyper`, and write the records through
    /// `output`, whose steps it sets. The reference paths, their offsets and ploidies,
    /// `genotype_snarls` (-a) and `allele_length_range` (-c and -C) are as for `FlowCaller`.
    /// `top_down` and `star_allele` are --top-down and -Y. Nothing is owned.
    MultiPassCaller(const PathPositionHandleGraph& graph, ReadLikelihoodSnarlCaller& genotyper,
                    const SnarlDecomposition& decomposition, VCFOutputCaller& output,
                    const string& sample_name, TraversalFinder& traversal_finder,
                    const vector<string>& ref_paths, const vector<size_t>& ref_path_offsets,
                    const vector<int>& ref_path_ploidies, bool genotype_snarls,
                    const pair<size_t, size_t>& allele_length_range, bool top_down,
                    bool star_allele);

    /// Run every pass, and add the records to the output's buffer. The direct pass visits the
    /// sites as `SiteWalker::walk` does with `recurse_type`.
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

    /// The header: the output's, with the read-likelihood genotyper's lines and this caller's.
    string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                      const vector<size_t>& contig_length_overrides = {}) const;

    /// Per-region ploidy overrides (--ploidy-bed).
    void set_ploidy_regions(PloidyRegions regions) { ploidy_regions = std::move(regions); }

    /// Nested calling: after `TreeGenotyper` genotypes a snarl, it descends into the snarl's child
    /// chains and genotypes them. It also turns on symbolic collapsing: a called traversal whose
    /// symbolic allele equals the reference traversal's is written as the reference allele, since
    /// it differs from the reference only inside child chains, whose own records report those
    /// differences.
    void set_nested_calling(bool on);

    /// Genotype and record chains that no reference path passes through, rather than skipping them.
    ///
    /// Such a chain has no REF or POS, so no record can be written for it, but it still takes part
    /// in the linkage model and gets anchors.
    void set_off_reference_nesting(bool on) { off_reference_nesting = on; }

    /// Write one record per difference block between the reference and each called strand's
    /// symbolic allele, instead of one record per snarl (--atomize-blocks).
    void set_atomize_blocks(bool on) { block_records.set_enabled(on); }

    /// Record a compact entry per site while calling, so that the linkage model can re-decide the
    /// genotypes afterwards. Neither pointer is owned; a null collector turns the model off.
    ///
    /// The GBWT and the haplotype of each of its sequences give the panel; see PanelLookup.
    void set_linkage(LinkageCollector* collector, const gbwt::GBWT* gbwt,
                     const vector<size_t>* sequence_to_haplotype);

    /// Write phased genotypes (`0|1`) and FORMAT/PS, from the chosen phase: the linkage model's
    /// Viterbi path, as read phasing reordered it when read phasing is on. Has no effect where the
    /// linkage model does not run.
    void set_emit_phasing(bool on) { this->emit_phasing = on; }

    /// `--min-confidence`, so that a record whose GQN the linkage model recomputes is marked
    /// against the same threshold.
    void set_linkage_min_confidence(double threshold) {
        this->linkage_min_confidence = threshold;
    }

    /// Write assembly anchors to `path`; see AnchorCollector::configure.
    void set_anchors_out(const string& path, const AnchorParams& params,
                         const string& graph_name, const string& reads_source,
                         double mismap_min) {
        anchor_collector.configure(path, params, graph_name, reads_source, mismap_min);
    }

    /// Turn on read phasing (--read-phasing); see read_phasing.hpp. Needs the linkage model, whose
    /// phase it changes.
    void set_read_phasing(bool on, const ReadPhasingParams& params) {
        read_phaser.configure(on, params);
    }

    /// Turn on re-genotyping from the phase (--regenotype); see regenotype.hpp. Needs read phasing,
    /// which gives each read its strand log-odds.
    void set_regenotype(bool on, const RegenotypeParams& params, size_t passes,
                        const string& ledger) {
        rescorer.configure(on, params, passes, ledger);
    }

    /// Where and how to write the mosaic. A path turns phasing on.
    void set_mosaic_out(MosaicParams params) {
        if (!params.path.empty()) {
            // The mosaic is the phasing, so phasing is on.
            this->emit_phasing = true;
        }
        mosaic_writer.set_params(std::move(params));
    }

private:

    /// The steps the output adds to `site`'s record: symbolic collapsing and the block count with
    /// nested calling, phasing from the linkage model, the GL layout of the read-likelihood
    /// genotyper, block records, and telling the linkage model the site's allele numbering. Each
    /// does nothing when its part is turned off. The steps read `site`, which must outlive them.
    VCFOutputCaller::SiteRecordSteps record_steps(const StagedSite& site);

    /// Set the output's header lines and writing steps this caller needs: the phase set and block
    /// lines of the header; before the records are written, the mosaic and the phasing report;
    /// on each line, the quality the linkage model moved; and after them, the refusals and the
    /// block report.
    void install_writer_steps();

    /// Write the mosaic, once every record exists, and report the phasing; separate from
    /// resolution because it needs to know which sites have a VCF line.
    void finalise_linkage_outputs();

    /// True when walk `trav_idx` of `site` takes the same route through the site as the reference
    /// walk and differs only inside child chains. Always false when nested calling is off.
    bool is_symbolically_reference(const StagedSite& site, int trav_idx, int ref_trav_idx) const;

    /// Phase a record's genotype from the linkage model's phase call for the site, as
    /// `SiteHooks::phase` does. The phased genotype must be a permutation of the record's own, so
    /// that phasing cannot change a genotype; a call that is not is counted as declined.
    int64_t phase_record_genotype(size_t record_key, const vector<int>& site_genotype,
                                  const map<int, int>& trav_to_allele, string& gt) const;

    /// A phase set is named by a position on its contig, so two contigs can share a name. Read
    /// phasing, re-genotyping and the anchors tell phase sets apart by an id instead, which
    /// `phase_set_id` gives each (contig, phase set) pair on first use.
    map<pair<string, size_t>, size_t> phase_set_ids;
    /// The id of a (contig, phase set) pair, the same for the whole run.
    size_t phase_set_id(const string& contig, size_t phase_set);

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


    /// the graph
    const PathPositionHandleGraph& graph;

    /// The read-likelihood genotyper, which also writes its fields into each record.
    ReadLikelihoodSnarlCaller& snarl_caller;

    /// The sites, as a snarl decomposition gives them.
    const SnarlDecomposition& decomposition;

    /// Walks the sites for the direct pass.
    SiteWalker walker;

    /// Where the records go.
    VCFOutputCaller& output;

    /// The sample the records are for.
    string sample_name;

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

    /// See set_ploidy_regions.
    PloidyRegions ploidy_regions;

    /// See set_nested_calling.
    bool nested_calling = false;

    /// See set_off_reference_nesting.
    bool off_reference_nesting = false;

    /// Every phased site. The linkage model fills it, the mosaic reads it, the linkage pass looks
    /// up a parent's chosen pair in it, and each record reads its phase from it as it is
    /// rendered.
    PhaseTable phase_table;

    /// Which allele each panel haplotype carries at a site. Empty without the linkage model.
    PanelLookup panel_lookup;
    /// The linkage model, and the linkage pass that chooses genotypes with it.
    GenotypeLinker linker;

    /// Records whose genotype the linkage model changed but whose quality fields could not be
    /// found on the line (no sample column, or FORMAT and sample columns of different lengths).
    /// They keep the per-site GQ.
    std::atomic<size_t> quality_declined{0};

    /// Phases refused while rendering because the record's genotype was not a permutation of the
    /// phased pair.
    mutable std::atomic<size_t> phase_declined{0};

    /// See set_linkage_min_confidence.
    double linkage_min_confidence = 0.0;

    /// Whether to emit phased GT and FORMAT/PS.
    bool emit_phasing = false;

    /// Writes the mosaic file, if one was asked for.
    MosaicWriter mosaic_writer;
    /// Writes sites as their difference blocks (see `set_atomize_blocks`), and counts block
    /// emission for the report.
    BlockRecordWriter block_records;

    /// Read phasing, if it is on.
    ReadPhaser read_phaser;

    /// Re-genotyping from the phase, if it is on.
    GenotypeRescorer rescorer;

    /// Collects the anchors and writes the anchor file, if one was asked for.
    AnchorCollector anchor_collector;

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
