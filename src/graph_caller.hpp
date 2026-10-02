#ifndef VG_GRAPH_CALLER_HPP_INCLUDED
#define VG_GRAPH_CALLER_HPP_INCLUDED

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

/// A set of traversals through a child snarl that are consistent with
/// a single parent allele. Multiple traversals can exist if the child
/// has internal variation within a shared region.
using TraversalSet = vector<SnarlTraversal>;

/// One TraversalSet per parent allele (index matches parent genotype).
/// For a diploid parent with genotype [0,1], element 0 contains traversals
/// consistent with parent allele 0, element 1 with parent allele 1.
using ChildTraversalSets = vector<TraversalSet>;

/// Counters for nested descent, reported under --progress. Members of the caller, so that each
/// run counts separately; `mutable` because the counting and reporting paths are const.
struct DescentCounters {
    /// How deep the descent went, by depth.
    std::atomic<size_t> depth_hist[16] = {};
    /// Children skipped because no reference path passes through them, and those recorded anyway
    /// because off-reference descent is on.
    std::atomic<size_t> skipped_no_ref{0};
    std::atomic<size_t> off_reference{0};
    std::atomic<size_t> no_ref_recorded{0};
    /// Off-reference chains by copy number: 0, 1, 2.
    std::atomic<size_t> no_ref_copies[3] = {};
    /// Children that no called parent allele crosses in the sweep. They are genotyped and kept,
    /// and the barrier decides from the parent's settled genotype whether the sample has them.
    std::atomic<size_t> skipped_no_copy{0};
    /// Children that a called traversal enters more than once. Only the first entry counts, both
    /// for ploidy and for distance: one traversal crossing a chain twice is one strand carrying two
    /// copies, not two strands, so it does not make the chain ploidy 2.
    std::atomic<size_t> child_multi_crossing{0};
};

/**
 * GraphCaller: Use the snarl decomposition to call snarls in a graph
 */
class GraphCaller {
public:

    enum RecurseType { RecurseOnFail, RecurseAlways, RecurseNever };
   
    GraphCaller(SnarlCaller& snarl_caller,
                SnarlManager& snarl_manager);

    virtual ~GraphCaller();

    /// Run call_snarl() on every top-level snarl in the manager.
    /// For any that return false, try the children, etc. (when recurse_on_fail true)
    /// Snarls are processed in parallel
    virtual void call_top_level_snarls(const HandleGraph& graph, RecurseType recurse_type = RecurseOnFail);

    /// Report what nested descent did: the depth histogram, and how many children it skipped and
    /// why. Does nothing in a run without nested descent.
    void report_descent_instrumentation() const;

    /// Counters for that report; see `DescentCounters`.
    mutable DescentCounters descent_counters;

    /// For every chain, cut it up into pieces using max_edges and max_trivial to cap the size of each piece
    /// then make a fake snarl for each chain piece and call it.  If a fake snarl fails to call,
    /// It's child chains will be recursed on (if selected)_
    virtual void call_top_level_chains(const HandleGraph& graph,
                                       size_t max_edges,
                                       size_t max_trivial,
                                       RecurseType recurise_type = RecurseOnFail);

    /// Call a given snarl, and print the output to out_stream
    virtual bool call_snarl(const Snarl& snarl) = 0;

    /// toggle progress messages
    void set_show_progress(bool show_progress);

    /// Visit top-level snarls in node-ID order, grouped into windows of window_size node IDs,
    /// instead of the default arbitrary order, so that a read source that fetches by node-ID
    /// window fetches each window once. Off by default; `vg call` turns it on for such a source.
    void set_node_id_ordering(bool ordered, size_t window_size);

protected:

    /// Break up a chain into bits that we want to call using size heuristics
    vector<Chain> break_chain(const HandleGraph& graph, const Chain& chain, size_t max_edges, size_t max_trivial);
    
protected:

    /// Our Genotyper
    SnarlCaller& snarl_caller;

    /// Our snarls
    SnarlManager& snarl_manager;

    /// See set_node_id_ordering.
    bool node_id_ordering = false;
    size_t node_id_window = 256;

    /// Toggle progress messages
    bool show_progress;
};

/// The order in which a caller wrote its `Number=G` GL vector. The two orders differ from three
/// alleles up: `PoissonSupportSnarlCaller` writes i-major (`for i; for j = i..n`), while
/// `ReadLikelihoodSnarlCaller` writes the VCF specification's colexicographic order. At n=3 they
/// differ at indices 2 and 3, `(1,1)` against `(0,2)`, so code that reindexes a GL vector must
/// know which it holds.
enum class GLLayout {
    IMajor,
    Colexicographic,
};

/// Index of genotype (i, j), i <= j, in a `Number=G` vector of the given layout.
size_t gl_genotype_index(size_t i, size_t j, size_t n_alleles, GLLayout layout);

/// Max-marginal fold of a diploid GL vector onto a smaller allele set.
///
/// `new_index[a]` is the allele `a` becomes. When several old alleles map to one new allele, their
/// genotypes merge, and the merged genotype takes the best likelihood among them. Separate from
/// `merge_similar_alleles` so that the unit tests can check it without a graph.
vector<double> fold_genotype_likelihoods(const vector<double>& old_gl,
                                         const vector<int>& new_index,
                                         size_t n_new, GLLayout layout);

/// Where a buffered VCF record sorts.
///
/// (contig, POS) is not a total order: several records can share a position, such as a nested
/// site under its parent, and `std::sort` is not stable, so ties would come out in an order that
/// varies between runs. `id` breaks the tie. On the FlowCaller path it is the snarl's own boundary
/// nodes (`print_snarl(snarl, false)`), so `vg call` on a graph is reproducible. On the
/// VCFGenotyper path it comes from the input VCF and is often ".", so `vg call -v` is not.
struct BufferedRecordKey {
    string contig;
    size_t position = 0;
    string id;
    /// Which difference block of its snarl this record is; 0 for a record written for a whole
    /// snarl. Two blocks of one snarl can land on the same POS, such as a deletion on one strand
    /// next to an insertion on the other, and share an ID, so the block number keeps the order
    /// total.
    size_t block = 0;
};

/// Strict weak ordering on BufferedRecordKey. A free function so that a unit test can check the
/// ordering directly.
bool buffered_record_key_less(const BufferedRecordKey& a, const BufferedRecordKey& b);



/// Counters for the mosaic writer; members of the caller, as for `DescentCounters`.
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

/// Counters for block emission; members of the caller, as for `DescentCounters`.
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
    std::atomic<size_t> refuse[10] = {};
};

/// Rewrite GQ, GQN and FILTER on one rendered VCF line whose genotype the linkage model changed,
/// so that they describe the settled genotype the line carries. GQ is -10 log10(1 - posterior)
/// times the direct call's `gq_factor`, capped at GQI. GQN is the settled genotype's margin in GL
/// over the best other genotype, divided by the direct call's achievable gap and multiplied by its
/// explained share, held within [-1, 1]. FILTER is decided again from that GQN against
/// `linkage_min_confidence`. Returns false if the line could not be parsed.
bool apply_linkage_quality(string& line, const LinkageCollector::MovedQuality& moved,
                           double linkage_min_confidence);

class VCFOutputCaller {
public:
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

    /// Write phased genotypes (`0|1`) and FORMAT/PS, from the settled phase: the linkage model's
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

    /// Turn one settled site into anchors, if anchors are being written. Called once per staged
    /// record as the records are rendered: just before its line is written, or, for a record with
    /// no line, by the hand-off. A record with no reference position still gets anchors, since a
    /// pin is placed by node ID, which is why this is not part of `emit_variant`.
    ///
    /// `is_leaf` is supplied by the caller, since the snarl manager is on GraphCaller. `gqn` is
    /// the value for the anchor's gqn column, from `FlowCaller::anchor_gqn_for`; NaN is written as
    /// `.`.
    void collect_anchors_for(const Snarl& snarl, const vector<int>& genotype, int haploid_slot,
                             const unique_ptr<SnarlCaller::CallInfo>& call_info, bool is_leaf,
                             double gqn, size_t record_key);


    /// The settled pair in phase order, for the anchors, which take each slot from the order of
    /// the pair they are given.
    ///
    /// `LinkageCollector::settled_traversals` returns a sorted pair; the phase is in
    /// `render_phases`. The pair is swapped when the record's PhaseCall names the same two
    /// traversals in the other order, and returned unchanged otherwise: no phasing, no PhaseCall,
    /// or a PhaseCall naming other traversals.
    vector<int> phase_ordered_genotype(size_t record_key, const vector<int>& genotype) const;

    /// Which strand a one-allele genotype sits on, for the anchors: 0 or 1.
    ///
    /// A nested chain at ploidy 1 is one strand of its parent, the one `nested_strand` names;
    /// `emit_variant` writes it as `a|.` or `.|a`, and this keeps the anchor's slot the same.
    /// Returns 0 when there is no phasing, no entry, a ploidy other than 1, no nested strand, or a
    /// phase that names a different allele from the settled one.
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
        /// linkage model settles the parent on a pair of traversals, and the VCF alleles are
        /// chosen only when the parent's line is written. When the mask cannot be computed (more
        /// than 64 traversals), it is 0 and `crossing_known` is false. The barrier reads it to find
        /// how many copies of the chain the parent's settled genotype carries, and on which
        /// strand.
        uint64_t parent_crossing = 0;
        /// Set where no called parent allele reaches the chain, and only when records are staged
        /// for the barrier and the linkage model runs (without it, such a chain is not genotyped).
        /// The chain is genotyped anyway, at the parent's ploidy, because the linkage model may
        /// still move the parent onto an allele that does reach it. Inherited by its children,
        /// which are genotyped at their own provisional ploidy.
        ///
        /// In the sweep the chain is staged, not written, and not recorded in the linkage model.
        /// The exception is a snarl whose own boundaries are on no reference path: that is recorded
        /// whatever this flag says, and never gets a line.
        ///
        /// The barrier decides what happens to it from the parent's settled pair. If the pair
        /// carries no copy, the chain and everything under it are dropped; so is a chain that no
        /// candidate traversal of the parent crosses (a crossing mask of 0). If the pair carries
        /// some, the chain is recorded at that many copies and rendered; where the sweep scored no
        /// genotype at that ploidy, it is rendered at the parent's ploidy instead, unrecorded. A
        /// parent the linkage model gave no phase call is read at its own settled genotype. Where
        /// the crossing mask is unknown (`crossing_known` false), the barrier cannot compare the
        /// chain with a settled pair at all, so the chain is never dropped on its parent's account
        /// (a dropped ancestor still removes it), and is rendered at the ploidy it was genotyped at,
        /// unrecorded.
        bool retain_only = false;
        /// Where this chain starts along the first of the parent's called traversals that crosses
        /// it, in bases, plus the parent's own `parent_offset`. Added to the parent's reference
        /// start, it gives an off-reference chain a position of its own, so that its sites are
        /// ordered as that traversal visits them and the distance between two of them is known.
        /// The barrier computes it again from the parent's settled genotype
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

    /// The generation of the site being recorded now: its depth among the nested chains, which
    /// decides the linkage pass that settles it. Thread-local because the sweep runs in parallel,
    /// and each descent saves, increments and restores it on its own thread.
    static thread_local size_t current_generation;

    /// Resolve one generation of the linkage model. `last` marks a barrier pass's final
    /// generation.
    void resolve_linkage_generation(size_t generation, bool last);

    /// Write the mosaic, once every record exists; separate from resolution because it needs to
    /// know which sites have a VCF line.
    void finalise_linkage_outputs();

    /// Resolve the linkage model, if nothing has resolved it yet: one pass per generation, from 0 to
    /// the deepest, which is looked up again after each pass since a pass can add a deeper chain.
    /// Safe to call more than once, so `write_variants` can call it unconditionally.
    void resolve_linkage();

    /// Every phased site, in the order the linkage model produced them. The mosaic reads this, and
    /// the barrier looks up a parent's settled pair in it.
    vector<LinkageCollector::PhaseCall> linkage_phased;
    bool linkage_resolved = false;
    /// Time spent in the linkage model, over all generations, for the report.
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

    /// Each diploid heterozygous site's per-read evidence, reduced to its settled pair, as read
    /// phasing last built it. Re-genotyping and the anchors read it.
    vector<PhaseSite> phase_sites;
    /// The record keys whose settled pair read phasing last reversed. Their sites' contributions
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
    /// The most barrier passes (--regeno-passes). With 1 the correction is computed and reported
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

    /// Fill `render_lambda` from the settled phasing. Called just before `build_render_phases`,
    /// when `phase_sites` and `phase_flips` are final.
    void build_render_lambda();

    /// This read's tempered strand log-odds, leaving out `record_key`, so that a site does not judge
    /// its own reads. Positive names slot 0. Zero means none: no table, no other contributing site,
    /// or no fitted temper. NaN means the read has a strand that is not usable in the site's phase
    /// set (see `read_strand_usable`). Re-genotyping does not use this, and gives such a read 0.
    double read_strand_log_odds(size_t record_key, const string& read_name) const;

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
    int emit_block_records(const PathPositionHandleGraph& graph, const Snarl& snarl,
                           const vector<SnarlTraversal>& called_traversals,
                           const vector<int>& genotype, int ref_trav_idx,
                           const string& sample_name, const vcflib::Variant& site,
                           const map<int, int>& trav_to_allele, int64_t site_position,
                           GLLayout gl_layout, bool genotype_snarls) const;

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

    /// A site's record key: the hash of the printed snarl, which is also the record's ID column.
    /// It identifies the site everywhere: in the linkage model, in the phasing and in the staged
    /// records.
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

    /// output buffers (1/thread) (for sorting)
    /// variants stored as strings (and position key pairs) because vcflib::Variant in-memory struct so huge
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

/**
 * Helper class for outputing snarl traversals as GAF
 */
class GAFOutputCaller {
public:
    /// The emitter object is created and owned by external forces
    GAFOutputCaller(AlignmentEmitter* emitter, const string& sample_name, const vector<string>& ref_paths,
                    size_t trav_padding);
    virtual ~GAFOutputCaller();

    /// print the GAF traversals
    void emit_gaf_traversals(const PathHandleGraph& graph, const string& snarl_name,
                             const vector<SnarlTraversal>& travs,
                             int64_t ref_trav_idx,
                             const string& ref_path_name, int64_t ref_path_position,
                             const TraversalSupportFinder* support_finder = nullptr);

    /// print the GAF genotype
    void emit_gaf_variant(const PathHandleGraph& graph, const string& snarl_name,
                          const vector<SnarlTraversal>& travs,
                          const vector<int>& genotype,
                          int64_t ref_trav_idx,
                          const string& ref_path_name, int64_t ref_path_position,
                          const TraversalSupportFinder* support_finder = nullptr);
    
    /// pad a traversal with (first found) reference path, adding up to trav_padding to each side
    SnarlTraversal pad_traversal(const PathHandleGraph& graph, const SnarlTraversal& trav) const;
    
protected:
    
    AlignmentEmitter* emitter;

    /// Sample name
    string gaf_sample_name;

    /// Add padding from reference paths to traversals to make them at least this long
    /// (only in emit_gaf_traversals(), not emit_gaf_variant)
    size_t trav_padding = 0;

    /// Reference paths are used to pad out traversals.  If there are none, then first path found is used
    unordered_set<string> ref_paths;

};

/**
 * VCFGenotyper : Genotype variants in a given VCF file
 */
class VCFGenotyper : public GraphCaller, public VCFOutputCaller, public GAFOutputCaller {
public:
    VCFGenotyper(const PathHandleGraph& graph,
                 SnarlCaller& snarl_caller,
                 SnarlManager& snarl_manager,
                 vcflib::VariantCallFile& variant_file,
                 const string& sample_name,
                 const vector<string>& ref_paths,
                 const vector<int>& ref_path_ploidies,
                 FastaReference* ref_fasta,
                 FastaReference* ins_fasta,
                 AlignmentEmitter* aln_emitter,
                 bool traversals_only,
                 bool gaf_output,
                 size_t trav_padding);

    virtual ~VCFGenotyper();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

protected:

    /// get path positions bounding a set of variants
    tuple<string, size_t, size_t>  get_ref_positions(const vector<vcflib::Variant*>& variants) const;

    /// munge out the contig lengths from the VCF header
    virtual unordered_map<string, size_t> scan_contig_lengths() const;

protected:

    /// the graph
    const PathHandleGraph& graph;

    /// input VCF to genotype, must have been loaded etc elsewhere
    vcflib::VariantCallFile& input_vcf;

    /// traversal finder uses alt paths to map VCF alleles from input_vcf
    /// back to traversals in the snarl
    VCFTraversalFinder traversal_finder;

    /// toggle whether to genotype or just output the traversals
    bool traversals_only;

    /// toggle whether to output vcf or gaf
    bool gaf_output;

    /// the ploidies
    unordered_map<string, int> path_to_ploidy;
};


/**
 * LegacyCaller : Preserves (most of) the old vg call logic by using 
 * the RepresentativeTraversalFinder to recursively find traversals
 * through arbitrary sites.   
 */
class LegacyCaller : public GraphCaller, public VCFOutputCaller {
public:
    LegacyCaller(const PathPositionHandleGraph& graph,
                 SupportBasedSnarlCaller& snarl_caller,
                 SnarlManager& snarl_manager,
                 const string& sample_name,
                 const vector<string>& ref_paths = {},
                 const vector<size_t>& ref_path_offsets = {},
                 const vector<int>& ref_path_ploidies = {});

    virtual ~LegacyCaller();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

protected:

    /// recursively genotype a snarl
    /// todo: can this be pushed to a more generic class? 
    pair<vector<SnarlTraversal>, vector<int>> top_down_genotype(const Snarl& snarl, TraversalFinder& trav_finder, int ploidy,
                                                                const string& ref_path_name, pair<size_t, size_t> ref_interval) const;
    
    /// we need the reference traversal for VCF, but if the ref is not called, the above method won't find it. 
    SnarlTraversal get_reference_traversal(const Snarl& snarl, TraversalFinder& trav_finder) const;

    /// re-genotype output of top_down_genotype.  it may give slightly different results as
    /// it's working with fully-defined traversals and can exactly determine lengths and supports
    /// it will also make sure the reference traversal is in the beginning of the output
    tuple<vector<SnarlTraversal>, vector<int>, unique_ptr<SnarlCaller::CallInfo>> re_genotype(const Snarl& snarl,
                                                                                              TraversalFinder& trav_finder,
                                                                                              const vector<SnarlTraversal>& in_traversals,
                                                                                              const vector<int>& in_genotype,
                                                                                              int ploidy,
                                                                                              const string& ref_path_name,
                                                                                              pair<size_t, size_t> ref_interval) const;

    /// check if a site can be handled by the RepresentativeTraversalFinder
    bool is_traversable(const Snarl& snarl);

    /// look up a path index for a site and return its name too
    pair<string, PathIndex*> find_index(const Snarl& snarl, const vector<PathIndex*> path_indexes) const;

protected:

    /// the graph
    const PathPositionHandleGraph& graph;
    /// non-vg inputs are converted into vg as-needed, at least until we get the
    /// traversal finding ported
    bool is_vg;

    /// The old vg call traversal finder.  It is fairly efficient but daunting to maintain.
    /// We keep it around until a better replacement is implemented.  It is *not* compatible
    /// with the Handle Graph API because it relise on PathIndex.  We convert to VG as
    /// needed in order to use it. 
    RepresentativeTraversalFinder* traversal_finder;
    /// Needed by above (only used when working on vg inputs -- generated on the fly otherwise)
    vector<PathIndex*> path_indexes;

    /// keep track of the reference paths
    vector<string> ref_paths;

    /// keep track of offsets in the reference paths
    map<string, size_t> ref_offsets;

    /// keep track of ploidies in the reference paths
    map<string, int> ref_ploidies;

    /// Tuning

    /// How many nodes should we be willing to look at on our path back to the
    /// primary path? Keep in mind we need to look at all valid paths (and all
    /// combinations thereof) until we find a valid pair.
    int max_search_depth = 1000;
    /// How many search states should we allow on the DFS stack when searching
    /// for traversals?
    int max_search_width = 1000;
    /// What's the maximum number of bubble path combinations we can explore
    /// while finding one with maximum support?
    size_t max_bubble_paths = 100;

};

/**
 * FlowCaller: takes each snarl's candidate traversals from a TraversalFinder and
 * genotypes them with its SnarlCaller (support-based, or ReadLikelihoodSnarlCaller
 * under --read-likelihood). It works on any graph. With the flow traversal finder it
 * does not report cyclic traversals; haplotype enumeration can. The `nested` constructor flag is --top-down, which genotypes each
 * child against traversal sets derived from its parent's called alleles. Nested
 * calling, the descent into child chains in call_snarl_internal, is turned on
 * instead by set_symbolic_collapsing.
 *
 * With the linkage model or nested calling, calling runs in stages. The sweep
 * genotypes every site from its own reads and stages a record for it (see
 * PendingRecord). The barrier settles the genotypes, parents before their children
 * (see run_barrier). Read phasing and re-genotyping may then change the phase and
 * the genotypes (see phase_and_regenotype). Finally every staged record is
 * rendered once, from its settled genotype (see render_retained_records).
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

    virtual bool call_snarl(const Snarl& snarl);

    /// Where a site sits for the linkage model, which orders sites and measures distances by it:
    /// the contig as the VCF names it, and the position where the site's first boundary node
    /// starts on the reference path. A record's POS can move once its alleles are trimmed, so it
    /// is not used. The model identifies a site by its record key instead.
    pair<string, size_t> site_ref_key(const Snarl& snarl, const string& ref_path_name,
                                      int ref_offset, bool no_reference = false,
                                      int64_t position_from_parent = 0) const;

    /// Record the site in the linkage model when it is genotyped, rather than when its line is
    /// written, since the barrier reads the collector before any line is written. The emitted
    /// allele map and whether a line was written are supplied later, by `set_allele_map`. The sweep
    /// does not call this for a retained chain with a reference path (see
    /// `NestedContext::retain_only`); the barrier records that chain if the sample carries it.
    void record_site(const Snarl& snarl, const vector<SnarlTraversal>& travs,
                     const vector<int>& trav_genotype,
                     const unique_ptr<SnarlCaller::CallInfo>& call_info, int ref_trav_idx,
                     const string& ref_path_name, int ref_offset,
                     bool no_reference = false, int64_t position_from_parent = 0);

    /// The frequency exponent a site should decode with: `--hp-prior` at a run-length site, or -1
    /// for the model's own. Reads the traversals' sequences only when `--hp-prior` is on.
    double site_freq_prior(const vector<SnarlTraversal>& travs, int ref_trav_idx) const;

    /// Decide every heterozygous site's phase from the reads, and change the settled phase to
    /// match, so that the GT order, the anchor slot column and the mosaic all follow from it.
    /// Genotypes are not changed. On FlowCaller because it needs the staged records, which hold
    /// the per-read evidence.
    void apply_read_phasing();

    /// Re-score every retained site's genotype likelihoods with the reads' phase.
    ///
    /// Runs after `apply_read_phasing`, which supplies `phase_sites` and `phase_flips`. With
    /// --regeno-passes above 1 the corrected likelihoods replace each site's own, and GQ is
    /// recomputed from them where the best genotype changed and is the sweep's elsewhere. Returns
    /// true if any site's corrected best genotype differs from its called one.
    bool apply_regenotyping();

    /// Feed the corrected likelihoods back to the linkage layer and settle again.
    ///
    /// The correction changes likelihoods, not genotypes: the linkage model still decides, as on
    /// the first pass.
    void regenotype_resettle();

    /// Read phasing, then rounds of re-genotyping, each followed by the barrier and read phasing
    /// again, as far as they are turned on. Does nothing unless read phasing is on.
    void phase_and_regenotype();

    /// Write every staged record once, from its settled genotype, and collect its anchors.
    void render_retained_records();

    /// Whether this snarl has no children, resolved through the manager's own copy. See the
    /// implementation for why the obvious `children_of(&snarl)` is not safe here.
    bool snarl_is_leaf(const Snarl& snarl) const;

    /// Stage every record during the sweep, and write it only after the barrier has settled its
    /// genotype. A nested chain's ploidy, its strand, and whether it has a record at all then come
    /// from its parent's settled genotype, and a parent is settled before its children, so a
    /// child's evidence cannot change its parent. Sizes the per-thread queues, so it must be
    /// called before calling starts.
    void set_settle_after_sweep(bool defer);

    /// How many records are staged for the render, reported under --progress.
    size_t render_record_count() const;


    /// The barrier: settle the genotypes one generation at a time. Each generation's linkage pass
    /// settles its sites; then each child chain of the next generation takes the ploidy its
    /// parent's settled genotype gives it, from the answers the sweep kept at both ploidies, and
    /// a chain the parent does not carry is dropped with everything inside it. Once every
    /// generation is settled, it decides which chains an enclosing block spells
    /// (`PendingRecord::reported_inline`). Does nothing unless records are staged (see
    /// `set_settle_after_sweep`).
    void run_barrier();

    /// Move every nested chain the barrier kept into the render's queues, and collect anchors for
    /// those that get no line. Separate from `run_barrier`, which re-genotyping runs again,
    /// because moving the records and collecting their anchors must happen once.
    void hand_off_deferred_records();



    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

    /// See max_snarl_edges. Zero removes the limit.
    void set_max_snarl_edges(size_t edges) {
        max_snarl_edges = edges ? edges : numeric_limits<size_t>::max();
    }

protected:

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

    /// A site's genotyping result, kept from the sweep until its record is rendered.
    ///
    /// A nested chain's ploidy depends on its parent's settled genotype, which is known only after
    /// the sweep, so the result is kept rather than computed again. The `CallInfo` has the answer at
    /// both ploidies (see `alt_ploidy_info`), so the record can be rendered at whichever ploidy the
    /// barrier settles on. `snarl` is held by value because `call_snarl_internal` may work on a
    /// flipped copy.
    struct PendingRecord {
        Snarl snarl;
        string ref_path_name;
        int ref_offset = 0;
        vector<SnarlTraversal> travs;
        int ref_trav_idx = -1;
        /// The sweep's genotype, before the linkage model, and its ploidy. A nested chain's
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
        /// sweep tests the parent's direct call; under the linkage model the barrier tests again
        /// with the parent's settled genotype, which the parent's blocks are built from.
        bool reported_inline = false;
        /// This snarl has no reference path, so no line can be written for it, since REF and POS are
        /// undefined. It is still genotyped and recorded in the linkage model.
        bool no_reference = false;
        /// For a snarl with no reference path: its parent's reference start plus `chain_offset`,
        /// standing in for the position it lacks.
        int64_t position_from_parent = 0;
        /// `NestedContext::parent_offset` for this chain. The sweep takes it from the parent's
        /// direct call, and the barrier computes it again from the parent's settled genotype, so
        /// that an off-reference chain is placed along the allele its parent settles on.
        size_t chain_offset = 0;
        /// See NestedContext::parent_crossing.
        uint64_t parent_crossing = 0;
        /// False when the parent has more than 64 candidate traversals, too many for
        /// `parent_crossing`. A 0 mask then means unknown, and the barrier leaves the chain at the
        /// ploidy the sweep gave it. The barrier computes the mask again when it revises or first
        /// records the parent.
        bool crossing_known = true;
        /// The site's generation.
        uint8_t generation = 0;
        /// Set when the settled parent, or an ancestor, does not carry this chain, so the chain and
        /// its descendants do not exist in the sample and are not revised or written. Each barrier
        /// pass decides it again, so a chain dropped in one pass can come back in the next.
        bool dropped = false;
        /// `panel_alleles(graph, travs)`, computed once: the traversals do not change after the sweep,
        /// and each re-genotyping round would otherwise repeat the GBWT lookups.
        vector<int> panel_cache;
        bool panel_cached = false;

    };

    /// The staged records a pass should look at, wherever they currently are: between a barrier
    /// pass and the hand-off, nested chains are in `deferred_pending` and the rest in
    /// `render_records`.
    ///
    /// With `for_phasing`, chains with no reference path are included: they cannot be rendered,
    /// having no REF or POS, but they are genotyped, get anchors, and have a meaningful strand.
    vector<PendingRecord*> records_for_render(bool for_phasing = false);

    /// `panel_alleles` for a record, computed once and kept. See `PendingRecord::panel_cache`.
    const vector<int>& cached_panel_alleles(PendingRecord& rec);

    /// The settled pair and ploidy per record, `{trav_first, trav_second, ploidy}`, for measuring
    /// whether a re-genotyping round changed anything.
    unordered_map<size_t, std::array<int, 3>> settled_snapshot();
    /// How many records settled differently from `before`, counting a chain that gained or lost a
    /// settled answer as moved.
    size_t settled_changed(const unordered_map<size_t, std::array<int, 3>>& before);
    /// A digest of a snapshot that does not depend on order, for spotting a state the rounds have
    /// reached before, which means they are cycling.
    static size_t snapshot_digest(const unordered_map<size_t, std::array<int, 3>>& snap);

    /// The nested chains the barrier settles, merged out of `pending_records` on the first pass and
    /// kept until `hand_off_deferred_records` moves them to the render. A member because the
    /// barrier runs once per re-genotyping round.
    vector<PendingRecord> deferred_pending;
    /// How many times the barrier has run.
    size_t barrier_passes_run = 0;


    /// See set_settle_after_sweep.
    bool settle_after_sweep = false;

    /// The nested chains' staged records, filled per thread during the sweep.
    vector<vector<PendingRecord>> pending_records;

    /// The top-level sites' staged records, and after the hand-off every record to be rendered.
    /// A top-level site's ploidy comes from the contig or the BED, so the barrier never revises
    /// it, though the linkage model still settles its genotype. Separate from
    /// `pending_records`, which `run_barrier` moves out and clears, and whose index groups
    /// records by parent. Read through `records_for_render`.
    vector<vector<PendingRecord>> render_records;


    /// Stage what a top-level site's record is rendered from. Takes the genotype and `call_info`
    /// now; the caller adds the traversals once descent, which still reads them, is done.
    unique_ptr<PendingRecord> stage_render_record(const Snarl& snarl,
                                                 const vector<int>& trav_genotype, int ref_trav_idx,
                                                 unique_ptr<SnarlCaller::CallInfo>& call_info,
                                                 const string& ref_path_name, int ref_offset,
                                                 int ploidy);


    /// The genotype the linkage model settled on for a staged record, or the sweep's where it
    /// settled none. Both anchor-collection paths use it.
    vector<int> settled_genotype_for(const PendingRecord& rec) const;

    /// The gqn column's value for this record: the sweep's `gq_fraction`, unless the linkage model
    /// changed the call, in which case the signed value recomputed for the settled genotype. NaN,
    /// written as `.`, where there is no value: no gap to normalise, or a moved call whose margin
    /// cannot be recomputed.
    double anchor_gqn_for(const PendingRecord& rec, const vector<int>& settled) const;

    /// `collect_anchors_for` for a staged record, with the phase order, the haploid slot and the
    /// leaf test derived from it. The genotype is a parameter because the render passes the
    /// settled pair.
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

    /// How far along `trav`, in bases, the child chain is entered: the total length of the nodes
    /// visited before it, or -1 if `trav` does not cross it. It gives an off-reference chain its
    /// place along its parent (see `NestedContext::parent_offset`).
    int64_t base_offset_of_child(const SnarlTraversal& trav, const Snarl& child) const;

    /// `base_offset_of_child` along the first traversal of `genotype` that crosses `child`, or 0
    /// when none does.
    size_t offset_along_genotype(const vector<SnarlTraversal>& travs, const vector<int>& genotype,
                                 const Snarl& child) const;
protected:

    /// The crossing mask: bit i is set where `travs[i]` crosses `child`. Indexed by traversal, not
    /// by VCF allele, since it is tested against the parent's settled traversals. Returns 0 and
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

class SnarlGraph;

/**
 * NestedFlowCaller : DEPRECATED - Use FlowCaller with nested=true instead.
 *
 * Uses any traversals finder (ex, FlowTraversalFinder) to find
 * traversals, and calls those based on how much support they have.
 * Should work on any graph but will not report cyclic traversals.
 * This class is being replaced by FlowCaller's nested mode.
 */
class NestedFlowCaller : public GraphCaller, public VCFOutputCaller, public GAFOutputCaller {
public:
    NestedFlowCaller(const PathPositionHandleGraph& graph,
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
                     bool genotype_snarls);

    virtual ~NestedFlowCaller();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

protected:

    /// stuff we remember for each snarl call, to be used when genotyping its parent
    struct CallRecord {
        vector<SnarlTraversal> travs;
        vector<pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>>> genotype_by_ploidy;
        string ref_path_name;
        pair<int64_t, int64_t> ref_path_interval;
        int ref_trav_idx; // index of ref paths in CallRecord::travs
    };
    typedef map<Snarl, CallRecord, NestedCachedPackedTraversalSupportFinder::snarl_less> CallTable;

    /// update the table of calls for each child snarl (and the input snarl)
    bool call_snarl_recursive(const Snarl& managed_snarl, int ploidy,
                              const string& parent_ref_path_name, pair<size_t, size_t> parent_ref_path_interval,
                              CallTable& call_table);

    /// emit the vcf of all reference-spanning snarls
    /// The call_table needs to be completely resolved
    bool emit_snarl_recursive(const Snarl& managed_snarl, int ploidy,
                              CallTable& call_table);

    /// transform the nested allele string from something like AAC<6_10>TTT to
    /// a proper string by recursively resolving the nested snarls into alleles
    string flatten_reference_allele(const string& nested_allele, const CallTable& call_table) const;
    string flatten_alt_allele(const string& nested_allele, int allele, int ploidy, const CallTable& call_table) const;

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

    /// until we support nested snarls, cap snarl size we attempt to process
    size_t max_snarl_shallow_size = 50000;

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

    /// a hook into the snarl_caller's nested support finder
    NestedCachedPackedTraversalSupportFinder& nested_support_finder;
};


/** Simplification of a NetGraph that ignores chains.  It is designed only for
    traversal finding.  Todo: generalize NestedFlowCaller to the point where we 
    can remove this and use NetGraph instead */
class SnarlGraph : virtual public HandleGraph {
public:
    // note: can only deal with one snarl "level" at a time
    SnarlGraph(const HandleGraph* backing_graph, SnarlManager& snarl_manager, vector<const Snarl*> snarls);

    // go from node to snarl (first val false if not a snarl)
    pair<bool, handle_t> node_to_snarl(handle_t handle) const;

    // go from edge to snarl (first val false if not a virtual edge)
    tuple<bool, handle_t, edge_t> edge_to_snarl_edge(edge_t edge) const;

    // replace a snarl node with an actual snarl in the traversal
    void embed_snarl(Visit& visit);
    void embed_snarls(SnarlTraversal& traversal);

    // replace a refpath through the snarl with the actual snarl in the traversal
    // todo: this is a bed of a hack
    void embed_ref_path_snarls(SnarlTraversal& traversal);

    ////////////////////////////////////////////////////////////////////////////
    // Handle-based interface (which is all identical to backing graph)
    ////////////////////////////////////////////////////////////////////////////
    bool has_node(nid_t node_id) const;
    handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
    nid_t get_id(const handle_t& handle) const;
    bool get_is_reverse(const handle_t& handle) const;
    handle_t flip(const handle_t& handle) const;
    size_t get_length(const handle_t& handle) const;
    std::string get_sequence(const handle_t& handle) const;    
    size_t get_node_count() const;
    nid_t min_node_id() const;
    nid_t max_node_id() const;
    
protected:

    bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// this is the only function that's changed to do anything different from the backing graph:
    /// it is changed to "pass through" snarls by pretending there are edges from into snarl starts out of ends and
    /// vice versa.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;    

    /// the backing graph
    const HandleGraph* backing_graph;

    /// the snarl manager
    SnarlManager& snarl_manager;

    /// the snarls (indexed both ways).  flag is true for original orientation
    unordered_map<handle_t, pair<handle_t, bool>> snarls;
};


}

#endif
