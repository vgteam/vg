#ifndef VG_MULTIPASS_CALLER_HPP_INCLUDED
#define VG_MULTIPASS_CALLER_HPP_INCLUDED

#include <atomic>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../src/handle.hpp"
#include "../src/snarls.hpp"
#include "../src/traversal_finder.hpp"
#include "anchor.hpp"
#include "block_records.hpp"
#include "candidate_finder.hpp"
#include "allele_splicer.hpp"
#include "descent_policy.hpp"
#include "edit_snarl_decomposition_overlay.hpp"
#include "fill_pins.hpp"
#include "genotype_linker.hpp"
#include "genotype_rescorer.hpp"
#include "../src/graph_caller.hpp"
#include "linkage_model.hpp"
#include "mosaic_writer.hpp"
#include "off_panel_ledger.hpp"
#include "panel_lookup.hpp"
#include "phase_table.hpp"
#include "ploidy_regions.hpp"
#include "read_likelihood_caller.hpp"
#include "read_consensus_step.hpp"
#include "read_phaser.hpp"
#include "read_strand_table.hpp"
#include "record_renderer.hpp"
#include "round_history.hpp"
#include "site_walker.hpp"
#include "staged_site.hpp"
#include "tree_genotyper.hpp"
#include "../src/site_emitter.hpp"
#include "../src/sorted_vcf_writer.hpp"

namespace vg::multipass {

using namespace std;



/**
 * The caller for every run that genotypes from read likelihoods.
 *
 * Its one method, `call`, runs the passes in order:
 *
 * 1. The direct pass genotypes every site from its own reads and stages it (see `TreeGenotyper`
 *    and `StagedSite`).
 * 2. Rounds follow. Each round's linkage pass chooses the genotypes one level at a time, parents
 *    before their children (see `GenotypeLinker`). Read phasing then decides the phases again
 *    (see `ReadPhaser`). From round 2 on, re-genotyping first corrects the likelihoods from the
 *    phase (see `GenotypeRescorer`). The rounds stop when the correction moves no site's direct
 *    call, when the chosen genotypes stop changing or return to an earlier round's, or at the
 *    round cap.
 * 3. The render builds each staged site's records once, from the genotype the last round chose
 *    (see `RecordRenderer`).
 * 4. The anchor file and the mosaic are written, and the passes are reported.
 *
 * With fill anchors on, the fill points are placed before the direct pass. The reads crossing them
 * are placed from each window the read source fetches whole during the pass, and from the other
 * windows, which are fetched after it while the read source is still open. The render marks which
 * called walks cross them.
 *
 * It walks the sites of a `SnarlDecomposition` itself, rather than being a `GraphCaller`. It
 * writes its records through a `SortedVCFWriter`, to which it adds its own header lines and
 * record steps.
 */
class MultiPassCaller {
public:
    /**
     * Make a caller of the sites of `graph`, as `decomposition` gives them.
     *
     * It genotypes each site with `genotyper`, from the candidate traversals `traversal_finder`
     * gives. It writes the records through `output`, whose steps it sets, on the paths of
     * `ref_paths`. Nothing else is owned.
     *
     * It skips a site unless its longest traversal is at least `allele_length_range.first` long
     * and none is longer than `allele_length_range.second`. With `genotype_snarls` it writes every
     * site, and trims alleles only as far as the boundary nodes, so that every sample gets the
     * same alleles.
     */
    MultiPassCaller(const PathPositionHandleGraph& graph, ReadLikelihoodSnarlCaller& genotyper,
                    const SnarlDecomposition& decomposition, SortedVCFWriter& output,
                    const string& sample_name, SiteTraversalFinder& traversal_finder,
                    ReferencePathTable ref_paths, bool genotype_snarls,
                    const pair<size_t, size_t>& allele_length_range);

    /**
     * Run every pass and add the records to the output's buffer, then write the anchors and the
     * mosaic and report on the passes.
     *
     * @param recurse_type Which sites the direct pass goes on to after the top level. It visits
     * the children of every site (`RecurseAlways`), of every site that `TreeGenotyper::genotype`
     * returns false for (`RecurseOnFail`), or of none (`RecurseNever`), and so on down.
     * @param release_reads Frees what the read source keeps for fetching reads. It runs once
     * the direct pass is done, before the linkage pass, and again after the read-consensus step
     * if that ran, which fetches reads. The other passes after the direct pass read only what it
     * kept, not reads.
     */
    void call(GraphCaller::RecurseType recurse_type, const function<void()>& release_reads);

    /// Get the emitter that writes the record of each staged site, so that its settings can be
    /// changed.
    SiteEmitter& emitter() { return site_emitter; }

    /// Get the sites this caller called: those of the snarl decomposition it was given, and the
    /// edit sites laid over it.
    const EditSnarlDecompositionOverlay& site_decomposition() const { return sites; }

    /// Return the ID of the records of a site of `site_decomposition()`, read the way `site`
    /// traverses it.
    ///
    /// A snarl is named from its boundary nodes. An edit site has an ID of its own, or none
    /// (an empty ID) if `call` wrote no record for it.
    string record_id_of_site(const net_handle_t& site) const;

    /// Do not genotype a site with more edges than this, including those of nested sites.
    ///
    /// The walk then genotypes the site's children as if they were top-level sites. Zero removes
    /// the limit, which is also the default.
    void set_max_snarl_edges(size_t edges) { candidates.set_max_snarl_edges(edges); }

    /// Batch the top-level sites of the direct pass into one parallel job per window of this many
    /// node IDs (see `SiteScheduler`). At 0, each top-level site is its own job.
    void set_snarl_batching(size_t window_size) { snarl_batch_window = window_size; }

    /// Toggle progress messages.
    void set_show_progress(bool show_progress) { this->show_progress = show_progress; }

    /// Add the header lines of this caller's own fields to the output. Call it once, after the
    /// settings that add fields (phasing, off-panel variants and block records).
    void add_header_lines();

    /// Set the ploidy of regions, in place of the ploidy of their contigs.
    void set_ploidy_regions(PloidyRegions regions) { ploidy_regions = std::move(regions); }

    /**
     * Turn nested calling on or off.
     *
     * With nested calling, after `TreeGenotyper` genotypes a site, it genotypes the site's child
     * chains. Nested calling also writes a called traversal whose symbolic allele equals the
     * reference traversal's as the reference allele. It differs from the reference only inside
     * child chains, whose own records report those differences.
     */
    void set_nested_calling(bool on);

    /// Genotype and file chains that the reference walk does not cross, rather than skipping
    /// them.
    ///
    /// Such a chain has no REF or POS, so it gets no record, but it still takes part in the
    /// linkage model and gets anchors.
    void set_off_reference_nesting(bool on) { off_reference_nesting = on; }

    /// Set whether groups of each site's walks, from `candidate_groups`, choose its genotype. The
    /// walks then choose within the chosen groups.
    void set_group_candidates(bool on) { group_candidates = on; }

    /// Write one record per difference block between the reference and the symbolic allele of
    /// each called strand, instead of one record per site.
    void set_atomize_blocks(bool on) { block_records.set_enabled(on); }

    /**
     * Set whether the top-level sites are written as flat records: the base-level differences
     * of each strand from the reference, over a run of neighbouring sites, with the called
     * alleles of the sites nested in them spliced in (see `plan_flat_records` and
     * `BlockRecordWriter::write_flat`). The sites nested in them then get no lines of their
     * own. Only the records change: the phase, the anchors and the mosaic are what they are
     * with this off. Off by default.
     *
     * A site whose strands cannot be aligned to its reference is written as it is with this
     * off, and so are the sites nested in it.
     *
     * With the environment variable `VG_CHECK_FLAT_RECORDS` set, the flat records of each run
     * are put back into its reference and compared with its strands, and the report counts the
     * runs checked and the runs that failed.
     */
    void set_flatten_records(bool on);

    /// Set the most bases that the nodes from one site's end bound to the next site's start
    /// bound, both included, can have for the two to be in one run of flat records. With 0,
    /// every site is a run of its own.
    void set_flatten_distance(size_t bases) { flatten_distance = bases; }

    /// File each site in the linkage model while calling, so that it can choose the genotypes
    /// afterwards.
    ///
    /// The GBWT and the haplotype of each of its sequences give the panel (see PanelLookup). No
    /// pointer is owned. A null collector turns the model off.
    void set_linkage(LinkageCollector* collector, const gbwt::GBWT* gbwt,
                     const vector<size_t>* sequence_to_haplotype);

    /// Write phased genotypes (`0|1`) and FORMAT/PS from the chosen phase.
    ///
    /// The chosen phase is the linkage model's Viterbi path, as read phasing reordered it when read
    /// phasing is on. This has no effect without the linkage model.
    void set_emit_phasing(bool on) { this->emit_phasing = on; }

    /// Set the GQN below which a record whose GQN the linkage model recomputes gets FILTER
    /// `lowconf`, or 0 for none.
    ///
    /// Pass the threshold given to `ReadLikelihoodSnarlCaller::set_min_confidence`, so that these
    /// records are marked against the same threshold as the rest.
    void set_linkage_min_confidence(double threshold) {
        this->linkage_min_confidence = threshold;
    }

    /// Write assembly anchors to `path`; see AnchorCollector::configure.
    ///
    /// Fill anchors are written too where `params.fill_spacing` is above 0, their reads scored by
    /// `fill_scores`.
    void set_anchors_out(const string& path, const AnchorParams& params,
                         const string& graph_name, const string& reads_source,
                         double mismap_min, const FillScoreParams& fill_scores) {
        anchor_collector.configure(path, params, graph_name, reads_source, mismap_min);
        fill_anchors.configure(path.empty() ? 0 : params.fill_spacing, fill_scores,
                               params.min_read_score);
    }

    /// Turn read phasing on or off (see `ReadPhaser`).
    ///
    /// It needs the linkage model, whose phase it changes.
    void set_read_phasing(bool on, const ReadPhasingParams& params) {
        read_phaser.configure(on, params);
    }

    /// Turn re-genotyping from the phase on or off (see `GenotypeRescorer`).
    ///
    /// It needs read phasing, which gives each read its strand log-odds.
    void set_regenotype(bool on, const RegenotypeParams& params, size_t passes) {
        rescorer.configure(on, params, passes);
    }

    /// Set the read-consensus step (see `ReadConsensusStep`), which runs where
    /// `options.enabled`.
    void set_read_consensus(const ReadConsensusOptions& options) {
        read_consensus.set_options(options);
    }

    /// Get what the read-consensus step did.
    const ReadConsensusCounts& read_consensus_counts() const { return read_consensus.counts(); }

    /// Set the source the reads of the sites come from.
    ///
    /// It can hold a wide site's reads for the sites below it (see `SiteReadSource::hold_reads`).
    void set_read_source(const SiteReadSource* reads) {
        read_source = reads;
    }

    /// Turn on off-panel detection.
    ///
    /// `reads` counts each read's substitutions as it fetches them (see
    /// `SiteReadSource::count_edits`).
    void set_off_panel(const SiteReadSource* reads, const EditCallBar& bar) {
        off_panel.configure(reads, bar);
    }

    /// Set where and how to write the mosaic.
    ///
    /// A path turns phasing on, since the mosaic is the phasing.
    void set_mosaic_out(MosaicParams params) {
        if (!params.path.empty()) {
            this->emit_phasing = true;
        }
        mosaic_writer.set_params(std::move(params));
    }

private:

    /**
     * The writer of the record of a staged site.
     *
     * It finds the staged site of each called site by its ID, among the sites the render is
     * writing. It does symbolic collapsing with nested calling, phasing from the linkage model,
     * block records, the quality fields of a genotype the linkage model moved, and telling the
     * linkage model the site's allele numbering. Each does nothing when its part is turned off.
     */
    class StagedSiteEmitter : public SiteEmitter {
    public:
        /// Make an emitter for the sites of `caller`, writing into `output`. Both must outlive it.
        StagedSiteEmitter(MultiPassCaller& caller, SortedVCFWriter& output);

    protected:
        bool same_as_reference(const CalledSite& site, int trav) override;
        int64_t phase(const CalledSite& site, const vector<int>& genotype,
                      const map<int, int>& trav_to_allele, string& gt) override;
        bool withholds_record(const CalledSite& site) override;
        int write_blocks(const PathPositionHandleGraph& graph, const CalledSite& site,
                         const SiteRecord& record) override;
        void finish_record(const CalledSite& site, vcflib::Variant& record) override;
        void on_site_done(const CalledSite& site, const map<int, int>& trav_to_allele,
                          bool has_line) override;

    private:
        MultiPassCaller& caller;
    };

    /// Get the staged site that the render is writing as `site`.
    const StagedSite& rendered_site(const CalledSite& site) const;

    /// Write the records of a staged site at `genotype`, during the render, and return how they
    /// numbered its alleles, or nothing where the site got no line.
    RecordRenderer::WrittenAlleles write_site(StagedSite& site, const vector<int>& genotype);

    /// Return whether the records `site` has without `--flatten-records` are kept out of the
    /// output: those of a site written as flat records, of a site nested in one, and of an edit
    /// site between two sites of a run.
    ///
    /// Such a site is still rendered, with its lines dropped, so that the linkage model, the
    /// anchors and the mosaic are told what they are told without `--flatten-records`.
    bool withholds_lines(const StagedSite& site) const {
        return site.flat.has_value() || flat_withheld.contains(site.record_key);
    }

    /// A run of sites that `--flatten-records` writes together, before its strands are spelled.
    struct FlatRunPlan {
        /// The sites that write the run's records: neighbouring top-level snarls of one chain,
        /// in the order of the chain, or one top-level edit site.
        vector<StagedSite*> sites;
        /// The nodes from each site's end bound to the next site's start bound, both included.
        vector<vector<handle_t>> junctions;
        /// The top-level edit sites in the nodes of the junctions.
        vector<StagedSite*> between;
    };

    /**
     * Return whether the strands of the top-level site `site` can be spelled for flat records,
     * counting why not.
     *
     * They can where the site has a genotype and a phase whose strands carry it, and its
     * resolved reference is its reference allele, which its records are placed by.
     */
    bool can_flatten(const StagedSite& site) const;

    /// Return whether the phase of the site with `phase` orders its two strands, so that its
    /// record claims a phase.
    bool orders_strands(const LinkageCollector::PhaseCall& phase) const;

    /**
     * Group the top-level sites `roots` into the runs that are written as flat records.
     *
     * A run is a longest series of snarls that follow each other in a chain, where at most
     * `flatten_distance` bases lie in the nodes from one's end bound to the next one's start
     * bound, one phase orders the strands of all of them, the reference path reads them one
     * after another, and the run's reference has no more than `flat_run_max_bases` bases. The
     * edit sites in the nodes between two sites of a run are in the run, and two snarls with an
     * edit site between them that cannot be are in no run together. Every other site that can
     * be written as flat records is a run of its own.
     *
     * @param can Whether each site of `roots` can be written as flat records.
     */
    vector<FlatRunPlan> flat_runs(const vector<StagedSite*>& roots, const vector<char>& can) const;

    /**
     * Spell each strand along the run of `plan`, align it to the run's reference, and give each
     * site of the run its part of the run's records, in `StagedSite::flat`.
     *
     * A strand is spelled by the splicer, from the first site's start bound to the last site's
     * end bound, with every site of the run in it.
     *
     * @param refused Set to why not, where false is returned.
     * @return False, giving no site anything, if a strand cannot be spelled or aligned.
     */
    bool align_flat_run(const FlatRunPlan& plan, FlatRefusal& refused) const;

    /// Return the objects the read-consensus step works with.
    ReadConsensusStep::Parts read_consensus_parts();

    /// Run the read-consensus step, if it is on, and choose the genotypes and phase the reads
    /// again if it gave any site an allele.
    void add_read_consensus_alleles();

    /**
     * Decide, before any record is written, which top-level sites are written as flat records,
     * in which runs, and so which other sites get no line.
     *
     * The line of a site nested in a site of a run, or of an edit site between two sites of a
     * run, is left out only once the run is known to align. A run of several sites that does
     * not is written as a run for each of its sites. A site alone that does not is written as
     * it is without flat records, with the sites nested in it.
     */
    void plan_flat_records();

    /// Give `record` the quality fields of its chosen genotype, if the linkage model moved its
    /// site's genotype (see ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype).
    void finish_moved_record(vcflib::Variant& record);

    /// Write the mosaic and report the phasing, once every record exists.
    ///
    /// Both need to know which sites have a VCF line.
    void finalise_linkage_outputs();

    /// Return true if walk `trav_idx` of `site` takes the same route through the site as the
    /// reference walk and differs only inside child chains.
    ///
    /// Return false when nested calling is off.
    bool is_symbolically_reference(const StagedSite& site, int trav_idx, int ref_trav_idx) const;

    /**
     * Phase a record's genotype from the linkage model's phase call for the site.
     *
     * Replace `gt` with the phased genotype, and return the phase set to write as PS, or -1 for
     * none. The phased genotype must be a permutation of the record's own, so that phasing cannot
     * change a genotype. A call that is not is counted as declined, and leaves `gt` unchanged.
     */
    int64_t phase_record_genotype(size_t record_key, const vector<int>& site_genotype,
                                  const map<int, int>& trav_to_allele, string& gt) const;

    /// The ID of each (contig, phase set) pair, by the pair.
    ///
    /// A phase set is named by a position on its contig, so two contigs can share a name. Read
    /// phasing, re-genotyping and the anchors tell phase sets apart by this ID instead.
    map<pair<string, size_t>, size_t> phase_set_ids;
    /// Get the ID of a (contig, phase set) pair, giving it one on first use.
    size_t phase_set_id(const string& contig, size_t phase_set);

    /// Configure the widgets of the passes with what they read from this caller.
    void install_widgets();

    /// The read source of the sites, or null when the caller has none of its own.
    const SiteReadSource* read_source = nullptr;
    /// Report what nested descent did: the depth histogram, and how many children it skipped and
    /// why.
    ///
    /// Do nothing in a run without nested descent.
    void report_descent_instrumentation() const;

    /// Run the linkage pass over the staged sites (see `GenotypeLinker::link`).
    void run_linkage_pass();

    /// Decide again which chains an enclosing block spells (`StagedSite::reported_inline`), from
    /// the chosen genotypes the render builds each parent's blocks from.
    ///
    /// A chain whose two strands carry different alleles is not spelled by a record that spells
    /// one traversal for both strands, so its line is not held back.
    ///
    /// @return The number of chains whose flag changed.
    size_t rederive_reported_inline();

    /// Report the memory kept for the staged sites, and what the linkage pass did.
    void report_linkage_pass(const GenotypeLinker::PassCounts& counts, size_t inline_rederived);

    /// Feed the corrected likelihoods back to the linkage layer and run the linkage pass again.
    ///
    /// The correction changes likelihoods, not genotypes. The linkage model still decides, as in
    /// round 1.
    void rerun_linkage_pass();

    /// Run read phasing (see `ReadPhaser`), if the linkage model is on.
    void phase_reads();

    /// Correct the likelihoods from the phase (see `GenotypeRescorer`), if the linkage model is
    /// on.
    ///
    /// @return True if the correction moved any site's direct call.
    bool rescore();

    /**
     * Run read phasing, then rounds of re-genotyping, each followed by the linkage pass and read
     * phasing again, as far as they are turned on.
     *
     * The rounds stop when the correction moves no site's direct call, or when the chosen
     * genotypes stop changing, return to an earlier round's (see `RoundHistory`), or reach the
     * round cap that `set_regenotype` sets.
     */
    void phase_and_regenotype();

    /// Run round `round` of re-genotyping, and report it.
    ///
    /// @return False if the rounds should stop.
    bool run_regenotyping_round(size_t round, RoundHistory& history);

    /// Build the records of every staged site once, from its chosen genotype, and collect the
    /// site's anchors (see `RecordRenderer`).
    void render_retained_records();

    const PathPositionHandleGraph& graph;

    /// The read-likelihood genotyper, which also writes its fields into each record.
    ReadLikelihoodSnarlCaller& snarl_caller;

    /// The sites: those of the snarl decomposition the caller was given, and the edit sites laid
    /// over it.
    EditSnarlDecompositionOverlay sites;

    /// Walks the sites for the direct pass.
    SiteWalker walker;

    /// Where the records go.
    SortedVCFWriter& output;

    /// Writes the record of each staged site into `output`.
    StagedSiteEmitter site_emitter;

    /// The staged sites the render writes, by record key, while it runs.
    unordered_map<size_t, StagedSite*> rendered_sites;

    /// The record ID of each edit site the render wrote, by its canonical net handle.
    unordered_map<net_handle_t, string> edit_record_ids;

    /// The sample the records are for.
    string sample_name;

    /// The reference paths, with their offsets and ploidies.
    ReferencePathTable ref_paths;

    /// Whether to write every site.
    ///
    /// By default, uncalled sites are skipped, and alleles are trimmed to make variants small.
    /// This turns both off.
    bool genotype_snarls;

    /// The width, in node IDs, of the windows the direct pass batches top-level sites by, or 0
    /// for no batching.
    size_t snarl_batch_window = 0;

    /// Whether to report progress to stderr.
    bool show_progress = false;

    /// The ploidy of regions, in place of the ploidy of their contigs.
    PloidyRegions ploidy_regions;

    /// Whether to genotype the child chains of each genotyped site, and write a called traversal
    /// whose symbolic allele equals the reference traversal's as the reference allele.
    bool nested_calling = false;

    /// Whether to genotype and file chains that the reference walk does not cross, rather than
    /// skip them.
    bool off_reference_nesting = false;

    /// Whether groups of each site's walks, from `candidate_groups`, choose its genotype.
    bool group_candidates = false;

    /// Every phased site.
    ///
    /// The linkage model fills it, the mosaic reads it, the linkage pass looks up a parent's
    /// chosen pair in it, and each record reads its phase from it as it is rendered.
    PhaseTable phase_table;

    /// Spells the alleles of each site with the chosen alleles of the sites under it, from
    /// `phase_table`.
    AlleleSplicer splicer;
    /// Whether the top-level sites are written as flat records (see `set_flatten_records`).
    bool flatten_records = false;
    /// See `set_flatten_distance`.
    size_t flatten_distance = 500;
    /// The most bases the reference of a run of several sites can have.
    static constexpr size_t flat_run_max_bases = 50000;
    /// The record keys of the sites in the render queues whose lines are kept out of the
    /// output, because the flat records of a run spell them: the sites nested in a site of the
    /// run, and the edit sites between two of its sites.
    unordered_set<size_t> flat_withheld;

    /// Which allele each panel haplotype carries at a site. Empty without the linkage model.
    PanelLookup panel_lookup;
    /// The linkage model, and the linkage pass that chooses genotypes with it.
    GenotypeLinker linker;

    /// The number of records whose genotype the linkage model changed but which have no FORMAT
    /// fields for the sample.
    ///
    /// They keep the per-site GQ.
    std::atomic<size_t> quality_declined{0};

    /// The number of phases refused while rendering because the record's genotype was not a
    /// permutation of the phased pair.
    mutable std::atomic<size_t> phase_declined{0};

    /// The GQN below which a record whose GQN the linkage model recomputes gets FILTER `lowconf`,
    /// or 0 for none.
    double linkage_min_confidence = 0.0;

    /// Whether to write phased GT and FORMAT/PS.
    bool emit_phasing = false;

    /// Writes the mosaic file, if one was asked for.
    MosaicWriter mosaic_writer;
    /// Writes sites as their difference blocks (see `set_atomize_blocks`), and counts the blocks
    /// for the report.
    BlockRecordWriter block_records;

    /// Read phasing, if it is on.
    ReadPhaser read_phaser;

    /// Re-genotyping from the phase, if it is on.
    GenotypeRescorer rescorer;
    /// The step that adds alleles built from the reads, if it is on.
    ReadConsensusStep read_consensus;
    /// Frees what the read source keeps for fetching reads (see `call`).
    function<void()> release_reads;

    /// Collects the anchors and writes the anchor file, if one was asked for.
    AnchorCollector anchor_collector;

    /// The fill points and the reads that cross them, if fill anchors are made.
    FillAnchors fill_anchors;

    /// Records what became of each off-panel candidate, and writes the summary.
    OffPanelLedger off_panel;

    /// Finds each site's reference path and candidate traversals. Declared after the members it
    /// reads.
    CandidateFinder candidates;

    /// Tallies what nested descent did.
    mutable DescentCounters descent_counters;

    /// Every staged site.
    ///
    /// A top-level site's ploidy comes from the contig or the BED, so the linkage pass never
    /// revises it, though the linkage model still chooses its genotype.
    StagedSiteTable staged_sites;

    /// What the reads say about the strands of each site, from read phasing.
    ///
    /// Re-genotyping and the anchors read it.
    ReadStrandTable read_strands;

    /// The temper re-genotyping fitted in its first round.
    TemperFit temper_fit;

    /// Builds the staged sites' records once the passes are done. Declared after the tables it
    /// reads.
    RecordRenderer record_renderer;

    /// Lists the child chains nested descent genotypes under each site.
    DescentPolicy descent_policy;

    /// Runs the direct pass for one top-level site.
    TreeGenotyper tree_genotyper;
};

}

#endif
