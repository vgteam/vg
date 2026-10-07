#ifndef VG_VCF_OUTPUT_CALLER_HPP_INCLUDED
#define VG_VCF_OUTPUT_CALLER_HPP_INCLUDED

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

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/// Special marker value for star alleles in genotype vectors.
/// A star allele (*) represents a haplotype that spans a nested site in the
/// parent but doesn't have a defined traversal at the child level.
constexpr int STAR_ALLELE_MARKER = -2;

/// Special marker value for missing alleles in genotype vectors.
/// Used when a parent allele doesn't traverse a child snarl and star_allele
/// mode is disabled. Outputs as '.' in VCF to maintain consistent ploidy.
constexpr int MISSING_ALLELE_MARKER = -1;

/// Counters for the mosaic writer. A member of each VCFOutputCaller, so that runs count
/// separately; `mutable` there because the writing paths are const.
struct MosaicCounters {
    /// Runs with no position to walk from. Panel haplotypes are often fragments, so this is
    /// reported rather than expected to be zero.
    std::atomic<size_t> unwalkable{0};
    /// Of those, the ones that are only a head: the run could be walked from a later site, so the
    /// walkable rest is written separately.
    std::atomic<size_t> head_clipped{0};
    /// Run boundaries across which the first run's haplotype could be followed, and those across
    /// which it could not, which a reference fill or a new fragment has to cover.
    std::atomic<size_t> extended{0}, gap_left{0}, patched{0};
    /// Rows whose own haplotype does not span them, rewritten as a reference substitution.
    std::atomic<size_t> row_to_ref{0};
    /// Run boundaries between a parent and a child snarl, at the child's own boundary nodes.
    std::atomic<size_t> nested_enter{0}, nested_leave{0};
    /// Rows the current direction could not walk but the other direction could, at an inversion.
    std::atomic<size_t> direction_broken{0}, extended_left{0};
};

/// Counters for block emission. A member of each VCFOutputCaller, as MosaicCounters is.
struct AtomizeCounters {
    /// Sites that reached `tally_atomize`, so that the report can tell "nothing refused" from
    /// "never ran".
    std::atomic<size_t> sites{0};
    std::atomic<size_t> site_unresolvable{0};  // flip_snarl left projection with no symbols
    std::atomic<size_t> site_reversed{0};      // resolved only via the reversed pairing
    /// Sites written as blocks, and the lines they produced.
    std::atomic<size_t> split_sites{0}, split_lines{0};
    /// Chains whose own record is not written because a block's ALT already spells them.
    std::atomic<size_t> child_inlined{0};
    /// Why `emit_block_records` declined a site, by refusal point. Each means the site's single
    /// record is written instead.
    std::atomic<size_t> refuse[13] = {};
};

/**
 * Helper class that VCF writers can inherit from, for the common code to output sorted VCF.
 *
 * It also holds the state of the linkage model, read phasing, the anchor file, the mosaic file
 * and block emission, which only FlowCaller uses.
 */
class VCFOutputCaller {
public:
    /// Where a buffered VCF record sorts: by contig, POS, `id` (the ID column), then `block`.
    /// Several records can share a contig and POS, such as a nested site and its parent, and
    /// `std::sort` is not stable, so the output order is the same on every run only when a caller
    /// gives the records at one position distinct (`id`, `block`) pairs.
    struct BufferedRecordKey {
        string contig;
        size_t position = 0;
        string id;
        /// Which difference block of its snarl this record is; 0 for a record written for a
        /// whole snarl. Two blocks of one snarl can land on the same POS, such as a deletion on
        /// one strand next to an insertion on the other, and share an ID, so the block number
        /// keeps the order total.
        size_t block = 0;
    };

    /// Strict weak ordering on BufferedRecordKey. Public, so that a unit test can check the
    /// ordering directly.
    static bool buffered_record_key_less(const BufferedRecordKey& a, const BufferedRecordKey& b);

    VCFOutputCaller(const string& sample_name);

    virtual ~VCFOutputCaller();

    /// Write the vcf header (version and contigs and basic info)
    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides) const;

    /// Add a variant to our buffer
    /// Returns false if the variant line length exceeds VCFOutputCaller::max_vcf_line_length
    bool add_variant(vcflib::Variant& var, size_t block = 0) const;

    /**
     * Per-region ploidy overrides, from a BED of `CHROM START END PLOIDY`.
     *
     * `-d` and `--ploidy-regex` set ploidy per contig, which cannot express a contig whose copy
     * number changes along it, such as a male sample's chrX, which is haploid except in the
     * pseudoautosomal regions.
     *
     * The CHROM column matches the contig name as it appears in the output VCF -- the locus part
     * of a PanSN path name, so `chrX` rather than `CHM13#0#chrX`. Intervals are BED half-open and
     * 0-based, and a position no interval covers keeps the contig's ploidy from `-d` or
     * `--ploidy-regex`.
     *
     * Overlapping intervals are an error, since a BED that says two things about one base has no
     * correct reading.
     */
    void set_ploidy_regions(const string& bed_path);


    /// Ploidy at this reference position, or `fallback` where no interval covers it. `position` is
    /// a 0-based offset along the contig, as in the BED.
    int region_ploidy(const string& ref_path_name, size_t position, int fallback) const;

    /// region_ploidy for a snarl whose reference interval begins at `interval_start`: the first
    /// base of its first boundary node, which is the record's POS less 1 before the record's
    /// alleles are trimmed. Returns `fallback` when no BED is loaded.
    int ploidy_at(const string& ref_path_name, int64_t interval_start, int64_t ref_offset,
                  int fallback) const;

    /// Record a compact entry per site while calling, so that the linkage model can re-decide the
    /// genotypes afterwards. Neither pointer is owned; a null collector turns the model off.
    ///
    /// The GBWT gives the panel: which allele each haplotype carries at a site, found by asking
    /// which haplotypes take each traversal. It is the GBWT that haplotype enumeration draws the
    /// candidate alleles from.
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

    /// Write assembly anchors to `path`, with `params` deciding which sites and reads qualify.
    void set_anchors_out(const string& path, const AnchorParams& params,
                         const string& graph_name, const string& reads_source,
                         double mismap_min) {
        this->anchor_path = path;
        this->anchor_params = params;
        if (this->anchor_params.counters == nullptr) {
            // A caller that did not supply counters, such as a unit test, counts into this
            // instance's own, so that every use below can dereference without a check.
            this->anchor_params.counters = &this->owned_anchor_counters;
        }
        this->anchor_graph_name = graph_name;
        this->anchor_reads_source = reads_source;
        this->anchor_mismap_min = mismap_min;
        if (!path.empty()) {
            // One queue per OpenMP thread, since the passes that fill it are parallel.
            this->anchor_writer = make_unique<AnchorWriter>((size_t)max(1, get_thread_count()));
        }
    }

    /// Turn one chosen site into anchors, if anchors are being written. Called once per staged
    /// site as the sites are rendered: just before its line is written, or, for a site with no
    /// line, by the hand-off. A site with no reference position still gets anchors, since a pin
    /// is placed by node ID, which is why this is not part of `emit_variant`.
    ///
    /// `is_leaf` is supplied by the caller, since the snarl manager is on GraphCaller. `gqn` is
    /// the value for the anchor's gqn column, from `FlowCaller::anchor_gqn_for`; NaN is written as
    /// `.`.
    void collect_anchors_for(const Snarl& snarl, const vector<int>& genotype, int haploid_slot,
                             const unique_ptr<SnarlCaller::CallInfo>& call_info, bool is_leaf,
                             double gqn, size_t record_key);


    /// The chosen pair in phase order, for the anchors, which take each slot from the order of
    /// the pair they are given.
    ///
    /// `LinkageCollector::chosen_traversals` returns a sorted pair; the phase is in
    /// `render_phases`. The pair is swapped when the record's PhaseCall names the same two
    /// traversals in the other order, and returned unchanged otherwise: no phasing, no PhaseCall,
    /// or a PhaseCall naming other traversals.
    vector<int> phase_ordered_genotype(size_t record_key, const vector<int>& genotype) const;

    /// Which strand a one-allele genotype sits on, for the anchors: 0 or 1.
    ///
    /// A nested chain at ploidy 1 is one strand of its parent, the one `nested_strand` names;
    /// `emit_variant` writes it as `a|.` or `.|a`, and this keeps the anchor's slot the same.
    /// Returns 0 when there is no phasing, no entry, a ploidy other than 1, no nested strand, or a
    /// phase that names a different allele from the chosen one.
    int phase_haploid_slot(size_t record_key, const vector<int>& genotype) const;

    /// Turn on read phasing (--read-phasing); see read_phasing.hpp. Needs the linkage model, whose
    /// phase it changes.
    void set_read_phasing(bool on, const ReadPhasingParams& params) {
        read_phasing = on;
        read_phasing_params = params;
    }

    /// Turn on re-genotyping from the phase (--regenotype); see regenotype.hpp. Needs read phasing,
    /// which gives each read its strand log-odds.
    void set_regenotype(bool on, const RegenotypeParams& params, size_t passes,
                        const string& ledger) {
        regenotype = on;
        regenotype_params = params;
        regenotype_passes = passes;
        regenotype_ledger = ledger;
    }

    /// Write the anchor file and report the counters. Does nothing unless anchors are on.
    void write_anchors();

    /// Whether the anchors need each snarl's leaf status, which has a cost to find (see
    /// `FlowCaller::snarl_is_leaf`).
    bool anchors_want_leaf_test() const {
        return !anchor_path.empty() && anchor_params.leaf_only;
    }

    /// Where to write the mosaic, if anywhere. Turns phasing on.
    ///
    /// `reference_paths` are the full names of the reference paths called against, such as
    /// `CHM13#0#chr20`. The rows give only the contig as the VCF names it, and a graph can hold
    /// several references, so the header lists them.
    void set_mosaic_out(const string& path, const string& graph_name,
                        const vector<string>& haplotype_names = {},
                        const vector<string>& reference_paths = {},
                        bool patch_gaps = true, bool keep_nested = true,
                        bool connect_unexplained = true) {
        this->mosaic_path = path;
        this->mosaic_graph_name = graph_name;
        this->mosaic_haplotype_names = haplotype_names;
        this->mosaic_reference_paths = reference_paths;
        this->mosaic_patch_gaps = patch_gaps;
        this->mosaic_keep_nested = keep_nested;
        this->mosaic_connect_unexplained = connect_unexplained;
        if (!path.empty()) {
            // The mosaic is the phasing, so phasing is on.
            this->emit_phasing = true;
        }
    }

    /// Write the buffered records. It adds the nesting INFO tags, sorts the records, runs the
    /// linkage model if nothing has (`resolve_linkage`), and writes the mosaic
    /// (`finalise_linkage_outputs`). Then it writes each record, rewriting GQ, GQN and FILTER on
    /// those whose genotype the linkage model changed (`LinkageCollector::moved_quality`). Usable
    /// once. `snarl_manager` is needed if
    /// `include_nested` is true.
    void write_variants(ostream& out_stream, const SnarlManager* snarl_manager = nullptr);

    /// Run vcffixup from vcflib
    void vcf_fixup(vcflib::Variant& var) const;

    /// Add a translation map
    void set_translation(const unordered_map<nid_t, pair<string, size_t>>* translation);

    /// Assume writing nested snarls is enabled
    void set_nested(bool nested);

    /// Genotype and record chains that no reference path passes through, rather than skipping them.
    ///
    /// Such a chain has no REF or POS, so no record can be written for it, but it still takes part
    /// in the linkage model and gets anchors.
    void set_off_reference_nesting(bool on) { off_reference_nesting = on; }

    /// How deep in non-reference sequence each gRef contig sits, by contig name. INFO/CH is at
    /// least this for a record on that contig.
    void set_gref_levels(map<string, int> levels);

    /// Enable post-genotyping merging of near-identical called ALT alleles, so that a 1/2 call of
    /// two effectively-identical alleles collapses to 1/1 with a single ALT.  Uses the same
    /// similarity metric and the same core-length gate as "vg deconstruct -L/--cluster-min-len" (a
    /// length-weighted Jaccard, except that a pure deletion is scored against the site -- see
    /// weighted_traversal_similarity).  The gate is applied to the alleles each tool emits, and
    /// those sets differ, so the two can disagree at a given site:
    /// similarity is >= threshold to merge, and min_len > 0 restricts merging to sites whose
    /// core length reaches min_len bp (see allele_core_length).
    /// A threshold of 1.0 (the default) disables merging entirely.
    void set_allele_merge(double threshold, int64_t min_len);

    /// The set of reference contigs that actually have a record.  Reads the sort keys of the
    /// output buffer, so it costs nothing (no decompression) and does not need the snarl tree.
    /// Only meaningful once calling is finished and before write_variants() drains the buffer.
    unordered_set<string> get_output_contigs() const;

    /// Remove ##contig lines whose ID is not in keep, leaving every other line alone.
    /// A reference contig that produced no record is not worth declaring: with a gref cover
    /// most contigs are fragments, and on a human chromosome a third of them carry nothing.
    string prune_header_contigs(const string& header, const unordered_set<string>& keep) const;

    /// Turn on symbolic collapsing: a called traversal whose symbolic allele equals the reference
    /// traversal's is written as the reference allele, since it differs from the reference only
    /// inside child chains, whose own records report those differences. The manager is not owned
    /// and must outlive this caller.
    ///
    /// In FlowCaller a non-null manager also turns on nested calling: after call_snarl_internal
    /// calls a snarl, it descends into the snarl's child chains and genotypes them.
    void set_symbolic_collapsing(const SnarlManager* manager) { this->symbolic_manager = manager; }

    /// Write one record per difference block between the reference and each called strand's
    /// symbolic allele, instead of one record per snarl (--atomize-blocks). Does nothing on the
    /// calling paths that cannot support it.
    void set_atomize_blocks(bool on) { this->atomize_blocks = on; }

protected:

    /// Whether `child` is already reported by this snarl's own records, because every called strand
    /// crosses it only inside a difference block whose ALT spells the route through it, so that
    /// block emission reports each variant once. This can happen only when no called allele is the
    /// reference allele; a chain that no reference path passes through is handled separately by the
    /// caller.
    bool chain_reported_inline(const Snarl& snarl, const vector<SnarlTraversal>& travs,
                               const vector<int>& genotype, int ref_trav_idx,
                               const Snarl& child) const;

    /// The parts of that test that do not depend on the child: the site's symbolic projection,
    /// each called ALT's projection, and the difference blocks between the reference and each ALT.
    /// Built once per snarl rather than once per child, since the edit-distance alignment in
    /// `symbolic_diff` is the same for every child.
    struct ChainInlineContext {
        /// False when the answer is false for every child: indices out of range, an empty genotype,
        /// an unresolvable site, or the reference among the called alleles.
        bool usable = false;
        SymbolicAllele sref;
        struct Alt {
            SymbolicAllele salt;
            vector<DiffBlock> blocks;
        };
        /// One entry per called allele that is in range and not the reference, in genotype order.
        vector<Alt> alts;
    };

    /// Build the child-independent half of the rule. See ChainInlineContext.
    ChainInlineContext build_chain_inline_context(const Snarl& snarl,
                                                  const vector<SnarlTraversal>& travs,
                                                  const vector<int>& genotype,
                                                  int ref_trav_idx) const;

    /// The part of the test that depends on the child. Gives the same result as the five-argument
    /// form.
    bool chain_reported_inline(const ChainInlineContext& ctx, const Snarl& child) const;

    /// True when this called traversal takes the same route through the snarl as the reference and
    /// differs only inside child chains. Always false when symbolic collapsing is off.
    bool is_symbolically_reference(const vector<SnarlTraversal>& called_traversals,
                                   int trav_idx, int ref_trav_idx, const Snarl& snarl) const;

    /// The parent of the nested chain being genotyped, for the duration of that call.
    ///
    /// Thread-local rather than a parameter: descent runs on the calling thread, so the context is
    /// set just before the child call and cleared after it, and no other thread sees it.
    struct NestedContext {
        /// Exactly one of the parent's called alleles crosses this chain (see
        /// `LinkageCollector::SiteContext::nested`).
        bool one_copy = false;
        size_t parent_record_key = 0;
        /// The crossing mask: one bit per parent candidate traversal, set where that traversal
        /// crosses this chain. It is indexed by traversal rather than by VCF allele, since the
        /// linkage model chooses the parent's genotype as a pair of traversals, and the VCF alleles are
        /// chosen only when the parent's line is written. When the mask cannot be computed (more
        /// than 64 traversals), it is 0 and `crossing_known` is false. The linkage pass reads it to find
        /// how many copies of the chain the parent's chosen genotype carries, and on which
        /// strand.
        uint64_t parent_crossing = 0;
        /// Set where no called parent allele reaches the chain, and only when staging is on (see
        /// `set_stage_records`) and the linkage model runs (without it, such a chain is not
        /// genotyped). The chain is genotyped anyway, at the parent's ploidy, because the linkage
        /// model may still move the parent onto an allele that does reach it. Inherited by its
        /// children, which are genotyped at their own provisional ploidy.
        ///
        /// In the direct pass the chain is staged, not written, and not recorded in the linkage model.
        /// The exception is a snarl whose own boundaries are on no reference path: that is recorded
        /// whatever this flag says, and never gets a line.
        ///
        /// The linkage pass decides what happens to it from the parent's chosen pair. If the pair
        /// carries no copy, the chain and everything under it are dropped; so is a chain that no
        /// candidate traversal of the parent crosses (a crossing mask of 0). If the pair carries
        /// some, the chain is recorded at that many copies and rendered; where the direct pass scored no
        /// genotype at that ploidy, it is rendered at the parent's ploidy instead, unrecorded. A
        /// parent the linkage model gave no phase call is read at its own chosen genotype. Where
        /// the crossing mask is unknown (`crossing_known` false), the linkage pass cannot compare the
        /// chain with a chosen pair at all, so the chain is never dropped on its parent's account
        /// (a dropped ancestor still removes it), and is rendered at the ploidy it was genotyped at,
        /// unrecorded.
        bool retain_only = false;
        /// Where this chain starts along the first of the parent's called traversals that crosses
        /// it, in bases, plus the parent's own `parent_offset`. Added to the parent's reference
        /// start, it gives an off-reference chain a position of its own, so that its sites are
        /// ordered as that traversal visits them and the distance between two of them is known.
        /// The linkage pass computes it again from the parent's chosen genotype
        /// (`PendingRecord::chain_offset`).
        size_t parent_offset = 0;
        /// Permission to genotype a chain that no reference path passes through. Inherited, since
        /// everything under such a chain is also off the reference. Whether a given snarl has a
        /// reference path is still checked for each snarl, since a descendant's boundaries may lie
        /// on a reference path even when its parent's do not.
        bool no_reference = false;
        /// Under block emission, an enclosing block's ALT already spells this chain's variation, so
        /// its own line would repeat it. The chain is still genotyped and recorded, and its line is
        /// held back when records are rendered. Inherited.
        bool reported_inline = false;
        /// Identifies the chain being descended into, from its boundary nodes. The linkage model
        /// groups a chain's sites by it, and chains under one parent have no transitions between
        /// them.
        size_t chain_key = 0;
        /// False when the parent has more than 64 candidate traversals, too many for the crossing
        /// mask. A 0 mask then means unknown rather than "no allele crosses".
        bool crossing_known = true;
    };
    static thread_local NestedContext nested_context;

    /// Whether the last emit on this thread described a real record, so that a child descending from
    /// it can tell that its parent's traversals are available. Thread-local and read just after the
    /// emit, as nested_context is.
    static thread_local bool last_emit_valid;

    /// Snarl hierarchy for symbolic collapsing, or null to compare alleles by sequence alone.
    const SnarlManager* symbolic_manager = nullptr;

    /// Whether to decompose a snarl into one record per difference block.
    bool atomize_blocks = false;

    /// The level of the site being recorded now: its depth among the nested chains, which
    /// decides when in a linkage pass its genotype is chosen. Thread-local because the direct pass
    /// runs in parallel, and each descent saves, increments and restores it on its own thread.
    static thread_local size_t current_level;

    /// Resolve one level of the linkage model. `last` marks a linkage pass's final
    /// level.
    void resolve_linkage_level(size_t level, bool last);

    /// Write the mosaic, once every record exists; separate from resolution because it needs to
    /// know which sites have a VCF line.
    void finalise_linkage_outputs();

    /// Resolve the linkage model, if nothing has resolved it yet: one pass per level, from 0 to
    /// the deepest, which is looked up again after each pass since a pass can add a deeper chain.
    /// Safe to call more than once, so `write_variants` can call it unconditionally.
    void resolve_linkage();

    /// Every phased site, in the order the linkage model produced them. The mosaic reads this, and
    /// the linkage pass looks up a parent's chosen pair in it.
    vector<LinkageCollector::PhaseCall> linkage_phased;
    bool linkage_resolved = false;
    /// Time spent in the linkage model, over all levels, for the report.
    double linkage_seconds = 0.0;
    size_t linkage_changed = 0;

    /// The linkage model's collector. Not owned.
    LinkageCollector* linkage_collector = nullptr;
    const gbwt::GBWT* linkage_gbwt = nullptr;
    const vector<size_t>* linkage_sequence_to_haplotype = nullptr;
    /// Panel size, so `panel_alleles` sizes its row by the haplotypes rather than by the GBWT's
    /// sequence count, which is larger under a gRef cover.
    size_t linkage_panel_size = 0;

    /// One cache of decompressed GBWT records per thread, for the panel lookups, since adjacent
    /// snarls share most of their records. Per thread because `CachedGBWT` has no locking. Sized
    /// in `set_linkage`, so that nothing allocates in the parallel region.
    mutable vector<gbwt::CachedGBWT> linkage_gbwt_cache;

    /// The node ID at which each thread's cache was started. `CachedGBWT` only grows, so when a
    /// thread's site is more than a fixed span of node IDs from this, its cache is cleared and
    /// started again there.
    mutable vector<nid_t> linkage_gbwt_cache_origin;

    /// Records whose genotype the linkage model changed but whose quality fields could not be
    /// found on the line (no sample column, or FORMAT and sample columns of different lengths).
    /// They keep the per-site GQ.
    mutable std::atomic<size_t> quality_declined{0};

    /// See `set_off_reference_nesting`.
    bool off_reference_nesting = false;

    /// Counters for the mosaic writer; see `MosaicCounters`.
    mutable MosaicCounters mosaic_counters;
    /// Counters for block emission; see `AtomizeCounters`.
    mutable AtomizeCounters atomize_counters;

    /// Print the block-emission counters.
    void report_atomize_instrumentation() const;

    /// `linkage_phased` keyed by record, copied by `build_render_phases` when phasing is emitted.
    /// Each record reads its phase from it as it is rendered. Keyed by record rather than by
    /// (contig, POS), since POS depends on which alleles the line carries.
    std::unordered_map<size_t, LinkageCollector::PhaseCall> render_phases;

    /// Read phasing: whether it is on, its parameters, and its counters. See read_phasing.hpp.
    bool read_phasing = false;
    ReadPhasingParams read_phasing_params;
    ReadPhasingCounters read_phasing_counters;

    /// A phase set is named by a position on its contig, so two contigs can share a name. Read
    /// phasing, re-genotyping and the anchors tell phase sets apart by an id instead, which
    /// `phase_set_id` gives each (contig, phase set) pair on first use.
    map<pair<string, size_t>, size_t> phase_set_ids;
    /// The id of a (contig, phase set) pair, the same for the whole run.
    size_t phase_set_id(const string& contig, size_t phase_set);

    /// Each diploid heterozygous site's per-read evidence, reduced to its chosen pair, as read
    /// phasing last built it. Re-genotyping and the anchors read it.
    vector<PhaseSite> phase_sites;
    /// The record keys whose chosen pair read phasing last reversed. Their sites' contributions
    /// to a read's strand log-odds enter with the opposite sign.
    unordered_set<size_t> phase_flips;

    /// Each read's strand log-odds, built once from `phase_sites` and `phase_flips` after the last
    /// read-phasing pass, for the anchors. Positive means strand 0 of the read's phase set, which
    /// is GT field 0 and anchor slot 0. Keyed by read alone, as a homozygous site, which has no
    /// PhaseSite, needs.
    LambdaTable render_lambda;
    /// `phase_sites` indexed by record, so that a site's own contribution can be subtracted. Points
    /// into `phase_sites`, which must not be rebuilt afterwards.
    unordered_map<size_t, const PhaseSite*> render_lambda_site;
    /// Each site's phase set, by record, since a read's strand is usable only at sites of the phase
    /// set it was found in (see `read_strand_usable`).
    unordered_map<size_t, size_t> render_lambda_phase_set;
    /// The fitted temper. Zero, when no fit was possible, makes every read's tempered strand
    /// log-odds zero, so no read counts as placed.
    double render_lambda_temper = 0.0;
    double render_lambda_ceiling = 1.0;

    /// Re-genotyping: whether it is on, its parameters, and its counters. See regenotype.hpp.
    bool regenotype = false;
    RegenotypeParams regenotype_params;
    RegenotypeCounters regenotype_counters;
    /// The most linkage passes (--regeno-passes). With 1 the correction is computed and reported
    /// but not applied.
    size_t regenotype_passes = 2;
    /// Where to write the ledger (--regeno-ledger): one line per site whose best genotype the
    /// correction changes. Empty for none.
    string regenotype_ledger;

    /// Phases refused while rendering because the record's genotype was not a permutation of the
    /// phased pair.
    mutable std::atomic<size_t> phase_declined{0};

    /// Copy `linkage_phased` into `render_phases`, keyed by record, when phasing is emitted.
    /// `render_retained_records` calls it on every staged run, after `build_render_lambda` and
    /// before any record is rendered. If read phasing ran, `linkage_phased` already carries its
    /// swaps.
    void build_render_phases();

    /// Fill `render_lambda` from the chosen phasing. Called just before `build_render_phases`,
    /// when `phase_sites` and `phase_flips` are final.
    void build_render_lambda();

    /// This read's tempered strand log-odds, leaving out `record_key`, so that a site does not judge
    /// its own reads. Positive names slot 0. Zero means none: no table, no other contributing site,
    /// or no fitted temper. NaN means the read has a strand that is not usable in the site's phase
    /// set (see `read_strand_usable`). Re-genotyping does not use this, and gives such a read 0.
    double read_strand_log_odds(size_t record_key, std::string_view read_name) const;

    /// See set_linkage_min_confidence.
    double linkage_min_confidence = 0.0;

    /// Whether to emit phased GT and FORMAT/PS.
    bool emit_phasing = false;

    /// Destination for the anchor file, and what qualifies for it. The writer is created once the
    /// thread count is known.
    string anchor_path;
    AnchorParams anchor_params;
    /// Used only when the run supplied no counters; see set_anchors_out.
    AnchorCounters owned_anchor_counters;
    string anchor_graph_name;
    string anchor_reads_source;
    unique_ptr<AnchorWriter> anchor_writer;
    /// The mismap floor that bounds a read's anchor confidence, for the file's header.
    double anchor_mismap_min = 0.0;

    /// Destination for the mosaic file, and the graph it is to be read against.
    string mosaic_path;
    string mosaic_graph_name;
    /// Panel index -> "sample#phase", the unit the linkage model works in: a haplotype stored as
    /// several GBWT paths is one haplotype. With a row's contig, that is enough to find its paths.
    /// The index means nothing outside this run, so the header writes the whole mapping.
    vector<string> mosaic_haplotype_names;

    /// Full reference path names the run called against; see set_mosaic_out.
    vector<string> mosaic_reference_paths;
    /// Fill a gap across which no panel haplotype can be followed with the reference, so that a
    /// strand stays one walk. On by default; the fill is marked `ref` in the file.
    bool mosaic_patch_gaps = true;
    /// Include nested sites in the runs, so that a switch of haplotype at a nested site starts a
    /// new row. On by default. Off leaves nested sites out, so a strand follows its enclosing
    /// site's haplotype through them, and its walk need not spell the nested sites' called alleles.
    /// Recorded in the file's #nested header.
    bool mosaic_keep_nested = true;
    /// Carry the flanking haplotype through a stretch the panel cannot explain, rather than
    /// writing an unwalkable row and breaking the path. On by default; it gives up the called
    /// alleles across those sites for a contiguous path.
    bool mosaic_connect_unexplained = true;

    /// How far a mosaic walk may run before it is abandoned. Walks go only in a direction already
    /// established, so this limits a long run rather than a wrong-way search.
    static const size_t MOSAIC_WALK_LIMIT = 1u << 17;
    /// Follow `hap` from an oriented node to `to_node`, and report where it arrives.
    ///
    /// The caller gives the direction. A GBWT stores each path in both orientations, so where a
    /// haplotype visits a node has two answers, while where it gets to along a known walk has one.
    bool mosaic_follow(gbwt::edge_type start, int64_t to_node, gbwt::node_type* out_end) const;
    /// Where `hap` sits at one oriented node, or invalid.
    gbwt::edge_type mosaic_position_at(gbwt::node_type node, size_t hap) const;
    /// (oriented node, haplotype) -> GBWT position. The mosaic is written serially, so one map with
    /// no locking is enough; it lives only for that pass.
    mutable std::unordered_map<uint64_t, gbwt::edge_type> mosaic_position_cache;
    gbwt::edge_type mosaic_gbwt_position(int64_t node_id, size_t hap) const;

    /// Collapse the per-site phasing into runs, stretches of sites over which a strand copies one
    /// haplotype, and write them.
    void write_mosaic(const vector<LinkageCollector::PhaseCall>& phasing) const;

    /// Which allele of `travs` each panel haplotype carries, or -1 where it does not traverse the
    /// site. Asks the GBWT which haplotypes take each traversal.
    vector<int> panel_alleles(const HandleGraph& graph,
                              const vector<SnarlTraversal>& travs) const;

    /// add a traversal to the VCF info field in the format of a GFA W-line or GAF path
    void add_allele_path_to_info(const HandleGraph* graph, vcflib::Variant& v, int allele,
                                 const Traversal& trav, bool reversed, bool one_based) const;
    /// legacy version of above
    void add_allele_path_to_info(vcflib::Variant& v, int allele, const SnarlTraversal& trav, bool reversed, bool one_based) const;
    
    
    /// convert a traversal into an allele string
    string trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const;

    /// Convert a SnarlTraversal to the handle vector the clustering code works on.  Returns false
    /// (leaving out_trav unspecified) if the traversal cannot be represented: fewer than two visits
    /// (the "*" placeholder pushed for a star allele), or a visit carrying a child Snarl rather
    /// than a node, which NestedFlowCaller produces via SnarlGraph::embed_snarl.  (LegacyCaller
    /// expands its children into node visits in top_down_genotype, so it never reaches here.)
    static bool snarl_traversal_to_handles(const HandleGraph& graph, const SnarlTraversal& trav,
                                           Traversal& out_trav);

    /// The CORE LENGTH of a variant: the length of the longest allele after stripping the prefix
    /// and the suffix that every non-"*" allele shares.  This is the single definition of "how big
    /// is this variant" behind --cluster-min-len in BOTH vg call and vg deconstruct.  It is
    /// invariant to how much shared flanking context a caller keeps in its allele strings, which is
    /// the point: vg call flattens down to an anchor base while vg deconstruct emits the whole
    /// snarl interior, so a raw string length answers differently for the same variant.
    /// Consequences, all intended:
    ///   - the anchor base flatten_common_allele_ends must leave on every indel is a shared prefix,
    ///     so it is stripped: a 49bp indel measures 49, not 50.
    ///   - REF participates, so a pure deletion measures the deleted length.  A maximum over ALTs
    ///     alone measures 1 for a deletion of any size.
    ///   - "*" is a marker, not sequence, so it is excluded from both the affixes and the maximum.
    ///     That also neutralizes flatten_common_allele_ends being a no-op whenever a "*" is
    ///     present -- without -a because min_allele_len becomes 1 and max_flatten_len decrements to
    ///     0, and with -a because "*" matches no base at the first offset compared.  Either way the
    ///     un-flattened boundary sequence is common to every real allele, so it is stripped here.
    /// Note this measures the SPAN of the variant, not the size of any one event inside it: a
    /// haplotype differing from the reference at two bases 59bp apart has a core length of 60.
    static int64_t allele_core_length(const vector<string>& alleles);

    /// Merge near-identical called ALT alleles in an already populated variant. Must run after
    /// SnarlCaller::update_vcf_info and flatten_common_allele_ends, so that both see every allele:
    /// merging earlier would drop the absorbed allele's reads from AD, DP and the Poisson caller's
    /// total_other_support term. Rewrites the allele-indexed fields (alleles/alt, AT, AD, GL, GT,
    /// MAD) and records the merge in INFO/MAT. Returns true if anything merged.
    ///
    /// `gl_layout` is the order in which the caller that produced this record wrote its GL, which
    /// cannot be recovered from the record.
    bool merge_similar_alleles(const PathPositionHandleGraph& graph,
                               const vector<SnarlTraversal>& site_traversals,
                               vector<int>& site_genotype,
                               const string& sample_name,
                               vcflib::Variant& out_variant,
                               GLLayout gl_layout) const;

    /// print a vcf variant
    /// return value is taken from add_variant (see above)
    bool emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                      const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                      const string& ref_path_name, int ref_offset, bool genotype_snarls, int ploidy,
                      function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_string = nullptr);

    /// get the interval of a snarl from our reference path using the PathPositionHandleGraph interface
    /// the bool is true if the snarl's backward on the path
    /// first returned value -1 if no traversal found 
    tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> get_ref_interval(const PathPositionHandleGraph& graph, const Snarl& snarl,
                                                                                 const string& ref_path_name) const;

    /// used for making gaf traversal names
    pair<string, int64_t> get_ref_position(const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name,
                                           int64_t ref_path_offset) const;

    /// clean up the alleles to not share common prefixes / suffixes
    /// if len_override given, just do that many bases without thinking
    void flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override) const;

    /// Split a finished site record into one record per difference block and file them.
    ///
    /// Returns the number of lines written, or -1 when it declines, in which case the site record
    /// is written as it is. `site` must be the finished record, after update_vcf_info and
    /// flattening, since every field a block does not redefine is taken from it.
    /// `trav_to_allele` maps each called traversal to its allele in `site`. It declines a site
    /// whose alleles `merge_similar_alleles` merged (`alleles_merged`), since the site then numbers
    /// its alleles differently from that map, and its blocks would spell the merged alleles apart.
    int emit_block_records(const PathPositionHandleGraph& graph, const Snarl& snarl,
                           const vector<SnarlTraversal>& called_traversals,
                           const vector<int>& genotype, int ref_trav_idx,
                           const string& sample_name, const vcflib::Variant& site,
                           const map<int, int>& trav_to_allele, int64_t site_position,
                           GLLayout gl_layout, bool genotype_snarls, bool alleles_merged) const;

    /// print a snarl in a consistent form like >3435<12222
    /// if in_brackets set to true,  do (>3435<12222) instead (this is only used for nested caller)
    // The nesting INFO headers (LV/CH/PS/RC/RS/RD), for both vg call and vg deconstruct.
    //
    // One definition on purpose.  These used to be written out verbatim in two places, and
    // drifted: 54bfd0f2d corrected the CH description in graph_caller.cpp while deconstructor.cpp
    // -- the copy deconstruct actually emits -- kept the text that commit's own message called
    // false, so the released VCFs carried the wrong one.
    static string nesting_info_headers();

    string print_snarl(const HandleGraph* grpah, const handle_t& snarl_start, const handle_t& snarl_end, bool in_brackets = false) const;
    /// legacy version of above
    string print_snarl(const Snarl& snarl, bool in_brackets = false) const;
    /// The same as above, but print the snarl as if its orientation has been flipped
    string print_flipped_snarl(const Snarl& snarl, bool in_brackets = false) const;
    /// What the three above print, from the snarl's two boundary visits.
    string print_snarl(nid_t start_id, bool start_backward, nid_t end_id, bool end_backward,
                       bool in_brackets) const;

    /// A site's record key: the hash of the printed snarl, which is also the record's ID column.
    /// It identifies the site everywhere: in the linkage model, in the phasing and in the staged
    /// sites.
    ///
    /// `write_variants` finds a buffered line's key by hashing its ID column, so the key must be the
    /// hash of that string. It survives `--translation`, where both sides print the translated
    /// form. One function, so that every caller and the recovery in `write_variants` agree.
    size_t record_key_of(const Snarl& snarl) const;

    /// do the opposite of above
    /// So a string that looks like AACT(>12<17)TTT would invoke the callback three times with
    /// ("AACT", Snarl), ("", Snarl(12,-17)), ("TTT", Snarl(12,-17))
    /// The parameters are to be treated as unions:  A sequence fragment if non-empty, otherwise a snarl
    void scan_snarl(const string& allele_string, function<void(const string&, Snarl&)> callback) const;

    // update the PS and LV tags in the output buffer (called in write_variants if include_nested is true)
    void update_nesting_info_tags(const SnarlManager* snarl_manager);
    
    /// output vcf
    mutable vcflib::VariantCallFile output_vcf;

    /// Sample name
    string sample_name;

    /// output buffers (1/thread) (for sorting) variants stored as strings (and position key pairs)
    /// because vcflib::Variant in-memory struct so huge
    mutable vector<vector<pair<BufferedRecordKey, string>>> output_variants;

    /// Reference interval of a site that was visited but not emitted, because every traversal
    /// through it was the reference (or absent) and so it had no variant to report.  Such a site
    /// is invisible to the RC/RS/RD walk, which only sees sites that reached the VCF, and a record
    /// nested under one would otherwise have no reference coordinate to point at.  Common in gref
    /// graphs, where the parent of an island of non-reference sequence is often a large snarl that
    /// only the reference and its own gref copy span.
    ///
    /// Keyed by snarl name as print_snarl() spells it, which is how record IDs and chrom_of_name
    /// are keyed too.  One buffer per thread, like output_variants, merged in
    /// update_nesting_info_tags().
    struct SuppressedRef {
        string chrom;
        size_t pos;
        size_t ref_len;
    };
    mutable vector<unordered_map<string, SuppressedRef>> suppressed_ref_info;

    /// print up to this many uncalled alleles when doing ref-genotpes in -a mode
    size_t max_uncalled_alleles = 5;

    /// Contig name -> ploidy overrides, sorted by start and not overlapping. Empty unless
    /// --ploidy-bed was given. See set_ploidy_regions.
    struct PloidyRegion {
        size_t start;   ///< 0-based, inclusive
        size_t end;     ///< 0-based, exclusive
        int ploidy;
    };
    unordered_map<string, vector<PloidyRegion>> ploidy_regions;

    // optional node translation to apply to snarl names in variant IDs
    const unordered_map<nid_t, pair<string, size_t>>* translation;

    // need to write LV/PS info tags
    bool include_nested;
    /// Contig name -> gRef nesting level, the least INFO/CH for that contig. Empty without a
    /// cover.
    map<string, int> gref_levels;

    // post-genotyping ALT merging (vg call -L / --cluster-min-len).  Deliberately NOT named
    // cluster_threshold / cluster_min_allele_len: Deconstructor derives from this class and already
    // declares both for its own pre-allele-string clustering, and -Wshadow is silent when a derived
    // member shadows a base one.
    double allele_merge_threshold = 1.0;
    int64_t allele_merge_min_len = 0;

    // prevent giant variants
    static const int64_t max_vcf_line_length = 2000000000;
};

}

#endif
