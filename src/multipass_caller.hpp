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

namespace vg::multipass {

using namespace std;

/**
 * A splicer that spells the alleles of a site with the chosen alleles of the sites staged under
 * it.
 *
 * A record spells each called ALT on the first strand of the site's phase that carries it. On
 * that strand, it spells every site staged under the site by its own chosen allele, recursively.
 * The splicer reads the chosen alleles from the phase table as the render has it.
 */
class AlleleSplicer {
public:
    /// The alleles a site's record spells, by traversal.
    ///
    /// A traversal the splice did not change is left as the site has it.
    struct SplicedAlleles {
        vector<Traversal> travs;
        AlleleSequences sequences;
    };

    /// One change to the sequence of one step of a spliced walk.
    struct StepEdit {
        /// The index of the step.
        size_t step;
        /// The offset in the step's sequence where the change starts.
        size_t at;
        /// The number of bases the change removes.
        size_t removed;
        /// The bases the change adds.
        string added;
        /// The record key of the child site that makes the change.
        size_t key;
    };

    /// Make a splicer over `graph` that reads the chosen alleles from `phase_table`.
    ///
    /// Neither is owned.
    AlleleSplicer(const HandleGraph& graph, const PhaseTable& phase_table);

    /// Index the staged sites under each site of `staged_sites`.
    ///
    /// Index nothing unless `phasing` is on and the phase table has its render phases, since
    /// without a phase there is no strand to splice on. The sites must not move until `clear()`.
    void index(StagedSiteTable& staged_sites, bool phasing);

    /// Forget the index.
    void clear();

    /// Return true if any staged site is indexed under the site with record key `record_key`.
    bool has_children(size_t record_key) const;

    /// Fill `out` with the alleles the record of `site` spells.
    ///
    /// @return False if no allele changed.
    bool splice(const StagedSite& site, SplicedAlleles& out) const;

private:
    /// The chosen allele of one child site on one strand, as the steps that replace the parent's
    /// visits [begin, end).
    struct Replacement {
        int begin;
        int end;
        Traversal walk;
        StepSequences steps;
        /// Whether the child lies inside nodes, so that its allele changes the bases of the
        /// parent's steps rather than replacing them.
        bool inside_node;
        /// The record key of the child.
        size_t key;
    };

    /**
     * Spell allele `trav` of `site` on strand `strand` of the site's phase, with each child staged
     * under it spelled by its chosen allele on that strand, recursively.
     *
     * @param haploid Whether the site has one strand, which all its one-copy children are on.
     * @param walk Set to the spliced walk.
     * @param steps Set to the sequence of each step of `walk`.
     * @return False, leaving `walk` and `steps` empty, if no child changes the allele.
     */
    bool spliced_allele(const StagedSite& site, int trav, int strand, bool haploid,
                        Traversal& walk, StepSequences& steps) const;

    /**
     * Make the replacement of a child site's chosen allele on one strand of its parent's walk.
     *
     * @param first_visit The index in `parent_walk` of the first visit to each node.
     * @return False if the child has no chosen allele on the strand, or the parent's walk does
     * not cross it.
     */
    bool chosen_replacement(const StagedSite& child, int strand, bool haploid,
                            const Traversal& parent_walk,
                            const unordered_map<nid_t, int>& first_visit, Replacement& rep) const;

    /// Spell `parent_walk` into `walk` and `steps`, with each of `routes` in place of the visits
    /// it covers.
    ///
    /// A route that overlaps one before it is skipped.
    ///
    /// @param own The step sequences of the parent's allele, or null if its walk spells it.
    /// @param routes Sorted by where they begin.
    /// @return True if any route changed the allele.
    bool splice_routes(const Traversal& parent_walk, const StepSequences* own,
                       const vector<Replacement>& routes, Traversal& walk,
                       StepSequences& steps) const;

    /// Find the changes that children inside nodes make to the steps of the spliced walk `walk`.
    ///
    /// Each child changes the bases its allele differs from its nodes in, where `walk` first
    /// visits its first node and then the rest of its walk, read the way `walk` reads them.
    vector<StepEdit> step_edits_of(const vector<Replacement>& edits, const Traversal& walk) const;

    /// Make the edits to the sequences of the steps, right to left within each step.
    ///
    /// Going right to left, an edit that changes the length does not move the offsets of those
    /// before it. An edit that overlaps one already made to its right is skipped.
    ///
    /// @return True if any sequence changed.
    static bool apply_step_edits(vector<StepEdit> step_edits, StepSequences& steps);

    const HandleGraph& graph;
    /// Gives the chosen allele of each site on each strand.
    const PhaseTable& phase_table;
    /// The staged sites under each site, by the record key of the parent.
    unordered_map<size_t, vector<const StagedSite*>> children;
};

/**
 * The record of what became of each off-panel edit candidate, which writes the dump and the
 * summary.
 *
 * The read source counts the candidates as it fetches reads. The direct pass files the outcome of
 * each candidate the call bar admits, from whichever thread genotyped it. Once every record is
 * written, the ledger reports the counts, and writes each candidate with its outcome to the dump,
 * placed on the reference.
 */
class OffPanelLedger {
public:
    /// Make a ledger that places candidates on `graph`, on the reference paths of `ref_paths`,
    /// shifted by their offsets.
    ///
    /// Neither is owned.
    OffPanelLedger(const PathPositionHandleGraph& graph, const ReferencePathTable& ref_paths);

    /// Turn off-panel detection on, for the candidates `reads` counts.
    ///
    /// `bar` decides which candidates become edit sites, and which edit sites are staged. With a
    /// `dump` path, every candidate is written there at the end.
    void configure(const SiteReadSource* reads, string dump, const EditCallBar& bar);

    /// Return whether off-panel detection is on.
    bool is_enabled() const { return candidate_source != nullptr; }

    /// Start filing outcomes, from thread numbers below `threads`.
    void open(size_t threads);

    /// Get the read source that counts the candidates, or null when detection is off.
    const SiteReadSource* source() const { return candidate_source; }

    /// Get the call bar.
    const EditCallBar& bar() const { return call_bar; }

    /// Get the outcomes, for the direct pass to file into.
    ///
    /// @return Null before `open`.
    EditOutcomes* outcomes() { return outcomes_filed.get(); }

    /**
     * Report what became of the candidates, and write every candidate to the dump, if there is
     * one.
     *
     * Do nothing when detection is off.
     *
     * @param nodes_after_pass The number of nodes whose top-level edit sites were genotyped after
     * the direct pass, for the report.
     */
    void write(size_t nodes_after_pass) const;

private:
    /// Report what became of the candidates.
    void report(size_t nodes_after_pass) const;

    /// Find where a node's base lies on the reference.
    ///
    /// @return The contig and the 1-based position of the base at forward offset `offset`, or an
    /// empty contig when no reference path visits the node.
    pair<string, int64_t> reference_position(nid_t node, uint32_t offset) const;

    const PathPositionHandleGraph& graph;
    /// The reference paths the candidates are placed on.
    const ReferencePathTable& ref_paths;

    /// The read source that counts the candidates, or null when detection is off.
    const SiteReadSource* candidate_source = nullptr;
    /// Where to write the dump, or empty for nowhere.
    string dump_path;
    /// Which candidates become edit sites, and which edit sites are staged.
    EditCallBar call_bar;
    /// What became of each candidate the call bar admitted.
    unique_ptr<EditOutcomes> outcomes_filed;
};

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
 * It walks the sites of a `SnarlDecomposition` itself, rather than being a `GraphCaller`. It
 * writes its records through a `VCFOutputCaller`, to which it adds its own header lines and
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
     * same alleles. With `top_down` it genotypes each child against the traversal sets its
     * parent's called alleles allow. With `star_allele` it gives a parent allele that skips a
     * child a star allele rather than a missing one.
     */
    MultiPassCaller(const PathPositionHandleGraph& graph, ReadLikelihoodSnarlCaller& genotyper,
                    const SnarlDecomposition& decomposition, VCFOutputCaller& output,
                    const string& sample_name, TraversalFinder& traversal_finder,
                    ReferencePathTable ref_paths, bool genotype_snarls,
                    const pair<size_t, size_t>& allele_length_range, bool top_down,
                    bool star_allele);

    /**
     * Run every pass and add the records to the output's buffer, then write the anchors and the
     * mosaic and report on the passes.
     *
     * @param recurse_type Which sites the direct pass visits, as for `SiteWalker::walk`.
     * @param after_direct_pass Runs once the direct pass is done, before the linkage pass. The
     * passes after it read only what the direct pass kept, not reads.
     */
    void call(GraphCaller::RecurseType recurse_type, const function<void()>& after_direct_pass);

    /// Get the tree of the sites `call` wrote records for, and of the sites enclosing them.
    ///
    /// The output's nesting INFO tags are computed from it. It is valid while this caller exists.
    unique_ptr<SiteTree> site_tree();

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

    /// Make the VCF header, with the lines of the read-likelihood genotyper and of this caller.
    string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                      const vector<size_t>& contig_length_overrides = {}) const;

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

    /// Write one record per difference block between the reference and the symbolic allele of
    /// each called strand, instead of one record per site.
    void set_atomize_blocks(bool on) { block_records.set_enabled(on); }

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
    void set_anchors_out(const string& path, const AnchorParams& params,
                         const string& graph_name, const string& reads_source,
                         double mismap_min) {
        anchor_collector.configure(path, params, graph_name, reads_source, mismap_min);
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
    void set_regenotype(bool on, const RegenotypeParams& params, size_t passes,
                        const string& ledger) {
        rescorer.configure(on, params, passes, ledger);
    }

    /// Set the source the reads of the sites come from.
    ///
    /// It can hold a wide site's reads for the sites below it (see `SiteReadSource::hold_reads`).
    void set_read_source(const SiteReadSource* reads) {
        read_source = reads;
    }

    /// Turn on off-panel detection.
    ///
    /// `reads` counts each read's substitutions as it fetches them (see
    /// `SiteReadSource::count_edits`). With a `dump` path, every candidate is written there at the
    /// end.
    void set_off_panel(const SiteReadSource* reads, const string& dump, const EditCallBar& bar) {
        off_panel.configure(reads, dump, bar);
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
     * Make the steps the output adds to the record of `site`.
     *
     * They do symbolic collapsing and the block count with nested calling, phasing from the
     * linkage model, the GL layout of the read-likelihood genotyper, block records, the quality
     * fields of a genotype the linkage model moved, and telling the linkage model the site's
     * allele numbering. Each does nothing when its part is turned off. The steps read `site`,
     * which must outlive them.
     *
     * @param spliced If not null, the alleles the site's record spells, which its blocks spell
     * too.
     */
    VCFOutputCaller::SiteRecordSteps record_steps(
        const StagedSite& site, const AlleleSplicer::SplicedAlleles* spliced = nullptr);

    /// Give `record` the quality fields of its chosen genotype, if the linkage model moved its
    /// site's genotype (see ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype).
    void finish_moved_record(vcflib::Variant& record);

    /// Set the header lines this caller adds to the output: the phase set and block lines.
    void install_header_steps();

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
     * `record_steps` makes this the `phase` step of each record. The phased genotype must be a
     * permutation of the record's own, so that phasing cannot change a genotype. A call that is
     * not is counted as declined.
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
    /// Without the linkage model every chosen genotype is the direct call, which the direct pass
    /// tested, so there is nothing to redo.
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

    /// The sites, as a snarl decomposition gives them.
    const SnarlDecomposition& decomposition;

    /// Walks the sites for the direct pass.
    SiteWalker walker;

    /// Where the records go.
    VCFOutputCaller& output;

    /// The sample the records are for.
    string sample_name;

    /// The reference paths, with their offsets and ploidies.
    ReferencePathTable ref_paths;

    /// Whether to write every site.
    ///
    /// By default, uncalled sites are skipped, and alleles are trimmed to make variants small.
    /// This turns both off.
    bool genotype_snarls;

    /// Genotype each child against the traversal sets its parent's called alleles allow (see
    /// `TreeGenotyper`).
    bool top_down = false;

    /// Whether to give a star allele, rather than a missing one, to a parent allele that skips a
    /// child, in top-down calling.
    bool star_allele = false;

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

    /// Every phased site.
    ///
    /// The linkage model fills it, the mosaic reads it, the linkage pass looks up a parent's
    /// chosen pair in it, and each record reads its phase from it as it is rendered.
    PhaseTable phase_table;

    /// Spells the alleles of each site with the chosen alleles of the sites under it, from
    /// `phase_table`.
    AlleleSplicer splicer;

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

    /// Collects the anchors and writes the anchor file, if one was asked for.
    AnchorCollector anchor_collector;

    /// Records what became of each off-panel candidate, and writes the dump and the summary.
    /// Declared after the reference offsets it reads.
    OffPanelLedger off_panel;

    /// Finds each site's reference path and candidate traversals. Declared after the members it
    /// reads.
    CandidateFinder candidates;

    /// Genotypes each site from its reads.
    SiteGenotyper site_genotyper;

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
    ChildPlacer child_placer;

    /// Runs the direct pass for one top-level site.
    TreeGenotyper tree_genotyper;
};

}

#endif
