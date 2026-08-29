#ifndef VG_GRAPH_CALLER_HPP_INCLUDED
#define VG_GRAPH_CALLER_HPP_INCLUDED

#include <atomic>
#include <iostream>
#include <algorithm>
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
#include "snarl_caller.hpp"
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

    /// Report what symbolic descent did: the depth histogram, and how many children it skipped and
    /// why. Stage 0 instrumentation for post-linkage nested descent; see the counters in
    /// graph_caller.cpp. Inert in a run with no symbolic descent.
    void report_descent_instrumentation() const;

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

    /// Visit top-level snarls in node-ID order, grouped into windows of window_size
    /// node IDs, instead of the default arbitrary order.
    ///
    /// For read sources that fetch by node-ID range: ordered access lets them touch
    /// each window exactly once and release it, rather than re-querying per site. Off
    /// by default, because it changes the visit order of code the default caller
    /// shares and that path's output must stay byte-identical.
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

/// Which order a caller wrote its `Number=G` GL vector in.
///
/// Two orders are live in this tree and they differ from three alleles up.
/// `PoissonSupportSnarlCaller` writes i-major -- `for i; for j = i..n` -- while
/// `ReadLikelihoodSnarlCaller` writes the VCF spec's colexicographic order. At n=3 they disagree at
/// exactly indices 2 and 3, `(1,1)` against `(0,2)`. Anything that reindexes a GL vector after the
/// fact has to know which it is holding; the allele-merge fold did not, and assumed i-major in a
/// comment that named the Poisson caller as the only writer.
enum class GLLayout {
    IMajor,
    Colexicographic,
};

/// Index of genotype (i, j), i <= j, in a `Number=G` vector of the given layout.
size_t gl_genotype_index(size_t i, size_t j, size_t n_alleles, GLLayout layout);

/// Max-marginal fold of a diploid GL vector onto a smaller allele set.
///
/// `new_index[a]` is the allele `a` becomes. Several old alleles mapping to one new allele means
/// their genotype classes collapse, and the merged class takes the best of its members -- the merged
/// allele scores as whichever of the alleles it absorbed fit the reads best.
///
/// Extracted from `merge_similar_alleles` so the layout handling is testable without a graph, a
/// traversal set or a `vcflib::Variant`. That matters because the defect it fixes is invisible in
/// output: on chr20, 0 of 57 merged records violate the "no genotype is strictly more likely than
/// the called one" invariant, so the transposition changes the emitted likelihoods without changing
/// which genotype is the argmax.
vector<double> fold_genotype_likelihoods(const vector<double>& old_gl,
                                         const vector<int>& new_index,
                                         size_t n_new, GLLayout layout);

/// Where a buffered VCF record sorts.
///
/// (contig, POS) alone is not a total order: several records legitimately share a position -- a
/// nested site under its parent, two snarls flattened onto the same anchor base -- and `std::sort`
/// is not stable, so ties came out in whatever order the per-thread buffers happened to concatenate
/// in. Two runs of one binary on chr20 differed in 72 record pairs that way, which makes every
/// "output must not move" gate unevaluable.
///
/// `id` breaks the tie and is intrinsic to the site rather than to the run: on the FlowCaller path
/// it is `print_snarl(snarl, false)`, the snarl's own boundary nodes. It is NOT intrinsic on the
/// VCFGenotyper path, where it comes from the input VCF and is commonly ".", so this makes
/// `vg call` on a graph reproducible and does not make `vg call -v` reproducible.
struct BufferedRecordKey {
    string contig;
    size_t position = 0;
    string id;
    /// Which difference block of its snarl this record is. Zero for every record a snarl emits
    /// whole, so the ordinary case sorts exactly as it did.
    ///
    /// The sort needs this rather than a suffixed ID. Two blocks of one snarl CAN land on the same
    /// POS -- a deletion on one haplotype abutting an insertion on the other, both anchored on the
    /// same base -- and with a shared ID the comparator would then be non-antisymmetric, which
    /// makes std::sort input-order dependent. That is the defect this header already records at
    /// the tie-break comment above, measured once at 72 differing record pairs between two runs.
    size_t block = 0;
};

/// Strict weak ordering on BufferedRecordKey. A free function rather than a lambda so a unit test
/// can assert the ordering property directly, which is the only place the totality of the key is
/// checked -- the in-tree TAP fixtures produce no ties at all.
bool buffered_record_key_less(const BufferedRecordKey& a, const BufferedRecordKey& b);

/// Stage 2 of planning/symbolic-diff-decomposition.md: what a symbolic diff would decompose,
/// measured from inside the caller rather than from INFO/AT offline. Output-neutral.
///
/// Must be called at the END of write_variants, not from the descent report. Most records are
/// retained and rendered after the calling sweep, so anything printed with the descent report
/// describes only the sites that emitted inline -- which on chr20 is none of them. That mistake
/// has been made here once already and produced a plausible-looking meaningless number.
void report_atomize_instrumentation();

/**
 * Helper class that vcf writers can inherit from to for some common code to output sorted VCF
 */
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
     * `-d` and `--ploidy-regex` set ploidy per *contig*, which cannot express the cases that
     * actually arise. A male sample's chrX is haploid except in the pseudoautosomal regions,
     * where X and Y recombine and two copies are present; before this, calling it meant running
     * chrX twice at different ploidies and splicing the two VCFs on the PAR boundaries. The same
     * shape appears wherever copy number varies along a contig.
     *
     * The CHROM column matches the contig name as it appears in the output VCF -- the locus part
     * of a PanSN path name, so `chrX` rather than `CHM13#0#chrX`. Intervals are BED half-open and
     * 0-based, and a position no interval covers keeps the contig's ploidy from `-d` or
     * `--ploidy-regex`.
     *
     * Overlapping intervals are an error rather than a silent last-one-wins: a BED saying two
     * things about one base has no correct reading, and picking one would move the ambiguity into
     * the output instead of into the error message.
     */
    void set_ploidy_regions(const string& bed_path);

    /// True when any override was loaded. Lookups are skipped when false, so a run without the
    /// option pays nothing for this.
    bool has_ploidy_regions() const { return !ploidy_regions.empty(); }

    /// Ploidy at this reference position, or `fallback` where no interval covers it. `position` is
    /// a 0-based offset along the contig, the same coordinate system as the BED and as the VCF POS
    /// the site will be emitted at.
    int region_ploidy(const string& ref_path_name, size_t position, int fallback) const;

    /// region_ploidy for a snarl whose reference interval begins at `interval_start`, applying the
    /// same offset arithmetic emit_variant uses for POS. Callers pass the interval they already
    /// computed, so no path lookup is repeated; returns `fallback` untouched when no BED is
    /// loaded, which keeps a default run free of this entirely.
    int ploidy_at(const string& ref_path_name, int64_t interval_start, int64_t ref_offset,
                  int fallback) const;

    /// Collect a compact record per site during calling, so genotypes can be re-decided by linkage
    /// once calling is done. Neither pointer is owned; a null collector disables the pass.
    ///
    /// The GBWT is needed because the panel matrix -- which allele each haplotype carries at a
    /// site -- comes from asking which haplotypes traverse each allele. That only exists under -z,
    /// where the alleles are haplotype-derived in the first place; without it there is no panel to
    /// link against.
    void set_linkage(LinkageCollector* collector, const gbwt::GBWT* gbwt,
                     const vector<size_t>* sequence_to_haplotype);

    /// Emit phased genotypes (`0|1`) and a FORMAT/PS phase set, from the linkage layer's Viterbi
    /// path. Off by default: it changes the GT of every record, which a naive parser may treat
    /// differently, so it rides with the feature rather than appearing unasked.
    void set_emit_phasing(bool on) { this->emit_phasing = on; }

    /// Where to write the run-length-encoded mosaic, if anywhere. Implies phasing.
    ///
    /// `reference_paths` are the full path names the run called against, e.g. `CHM13#0#chr20`.
    /// They are recorded because the segment rows carry only the *locus* part of the name -- the
    /// contig as the VCF spells it -- and a graph may hold more than one reference. The HPRC
    /// graphs do: their GBZ tags name both `CHM13` and `GRCh38` as reference samples, and both
    /// appear in the panel, so a bare `chr20` does not say which assembly a position is measured
    /// against. Without this the coordinates are ambiguous and nothing in the file says so.
    void set_mosaic_out(const string& path, const string& graph_name,
                        const vector<string>& haplotype_names = {},
                        const vector<string>& reference_paths = {}) {
        this->mosaic_path = path;
        this->mosaic_graph_name = graph_name;
        this->mosaic_haplotype_names = haplotype_names;
        this->mosaic_reference_paths = reference_paths;
        if (!path.empty()) {
            this->emit_phasing = true;
        }
    }

    /// Sort then write variants in the buffer
    /// snarl_manager needed if include_nested is true
    void write_variants(ostream& out_stream, const SnarlManager* snarl_manager = nullptr);

    /// Run vcffixup from vcflib
    void vcf_fixup(vcflib::Variant& var) const;

    /// Add a translation map
    void set_translation(const unordered_map<nid_t, pair<string, size_t>>* translation);

    /// Assume writing nested snarls is enabled
    void set_nested(bool nested);

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
    /// traversal's is emitted as the reference allele, because it differs from the reference only
    /// inside child chains and those differences belong to those chains' own records.
    ///
    /// Without it such a traversal becomes a long ALT differing from REF at a handful of bases.
    /// 90.6% of vg's same-length structural false positives are that shape, and 55,222 of its
    /// autosomal SNV false negatives sit inside one. See the nested-calling design note in the
    /// companion evaluation repository (docs/nested-calling-design.md).
    ///
    /// The manager is not owned and must outlive this caller.
    void set_symbolic_collapsing(const SnarlManager* manager) { this->symbolic_manager = manager; }

    /// Emit one record per difference block between the reference and each called haplotype's
    /// symbolic allele, instead of one record per snarl. See
    /// planning/symbolic-diff-decomposition.md. Off by default, and refused outright rather than
    /// silently declined for the configurations it cannot serve.
    void set_atomize_blocks(bool on) { this->atomize_blocks = on; }

protected:

    /// True when this called traversal takes the same route through the snarl as the reference and
    /// differs only inside child chains. False whenever symbolic collapsing is off, so the default
    /// path is unchanged.
    /// Whether `child` is already reported by this snarl's own records, because every called
    /// haplotype crosses it only inside a difference block whose ALT spells the route through it.
    ///
    /// This is the exactly-once rule from planning/symbolic-diff-decomposition.md. Its population
    /// is narrow: a chain the reference does not cross is not descended into at all (the caller
    /// gates on that separately), and a genotype carrying the reference allele matches the chain by
    /// definition. What is left is a snarl called with no reference allele where the alignment
    /// matched the chain on neither haplotype -- loops and reorderings.
    bool chain_reported_inline(const Snarl& snarl, const vector<SnarlTraversal>& travs,
                               const vector<int>& genotype, int ref_trav_idx,
                               const Snarl& child) const;

    bool is_symbolically_reference(const vector<SnarlTraversal>& called_traversals,
                                   int trav_idx, int ref_trav_idx, const Snarl& snarl) const;

    /// Which parent a nested call hangs off, for the duration of that call.
    ///
    /// A nested site has ploidy 1 exactly when one called parent allele crosses the child chain and
    /// the other does not, so the strand it belongs to is determined by the parent's phase rather
    /// than estimated. The linkage layer needs to know the parent's record key and which of its two
    /// genotype slots did the crossing.
    ///
    /// Thread-local rather than a parameter threaded through emit_variant's already long signature:
    /// descent happens synchronously on the calling thread, so the context is set immediately before
    /// the child call and cleared after it, and no other thread can observe it.
    struct NestedContext {
        bool active = false;
        size_t parent_record_key = 0;
        /// The parent traversal this child hangs off, or -1 when the parent does not carry it on
        /// exactly one. The traversal itself, not the strand it sits at: the strand is whichever
        /// haplotype that traversal is phased onto, looked up when the parent is phased.
        /// One bit per parent VCF allele, set where that allele crosses this child chain. Zero
        /// means descent could not express it -- more than 64 alleles at the parent -- and must be
        /// read as unknown rather than as none.
        ///
        /// Linkage uses it to check that the ploidy this child was called at survives the parent's
        /// *final* genotype, which the slot above cannot answer: the slot names one strand, where the
        /// question is whether both called parent alleles cross the child, neither does, or one does.
        /// See LinkageCollector::resolve.
        uint64_t parent_crossing = 0;
        /// True when this chain, and everything under it, must be genotyped but not emitted.
        ///
        /// Set where no called parent allele reaches the chain. It is visited anyway -- its parent's
        /// genotype is not settled yet, and linkage may still move the parent onto an allele that does
        /// reach it -- but nothing about it may reach the VCF until the barrier says the sample
        /// carries it. Inherited by its own children: a chain whose existence depends on an ancestor
        /// the sample may not have cannot be emitted either.
        bool retain_only = false;
        /// Permission to genotype a chain the reference does not cross. Inherited, like
        /// `retain_only`: everything under such a chain is also off the reference. It is only
        /// PERMISSION -- whether a given snarl actually has a reference path is re-derived per
        /// invocation from the graph, because a descendant's own boundaries may sit on a declared
        /// reference path even when its parent's did not.
        bool no_reference = false;
        /// The index of the chain member being descended into, within its chain. Set beside the
        /// other per-child facts before the recursive call, and read by that child for its own
        /// entry -- the same way `parent_crossing` travels.
        /// Identity of the chain being descended into, from its boundary pair. The decode
        /// groups on this: sibling chains have no transition between them, so a chain must be
        /// distinguishable from its siblings and nothing more.
        size_t chain_key = 0;
        int chain_index = -1;
        /// False when the crossing mask could not be computed: the parent emitted nothing on this
        /// invocation, or carries too many alleles for a 64-bit mask. child_crossing_mask returns 0
        /// for "unknown", and the barrier must not read that 0 as "no allele crosses" -- so the
        /// distinction travels beside the mask instead of being conflated with it.
        bool crossing_known = true;
    };
    static thread_local NestedContext nested_context;

    /// What the last emit_variant call on this thread turned each called traversal into, and how
    /// many alleles the record ended up with.
    ///
    /// Symbolic descent needs a child's crossing pattern in VCF allele space, because that is the
    /// space GT, GL and the linkage layer are written in, and only emit_variant knows the mapping --
    /// alleles are deduplicated by string and symbolically-reference traversals all collapse onto
    /// allele 0. Thread-local and read immediately after the emit that filled it, for the same
    /// reason nested_context is.
    struct EmittedAlleles {
        bool valid = false;
        size_t num_alleles = 0;
        map<int, int> trav_to_allele;
        /// The contig and POS the record was actually written with -- after
        /// flatten_common_allele_ends has moved POS -- which is the key add_variant files the line
        /// under. A linkage entry recorded from anything else (get_ref_position, say) sits at the
        /// pre-flatten position and write_variants' lookup never finds it.
        string contig;
        size_t position = 0;
        /// Where add_variant buffered the line, so the barrier can retract or replace that exact
        /// line instead of re-identifying it by hashing the ID column and counting GT separators.
        /// buffer_thread is -1 when the emit buffered nothing: the line was rejected for length, or
        /// the record flattened to nothing and emit_variant returned true without writing.
        int buffer_thread = -1;
        size_t buffer_index = 0;
    };
    static thread_local EmittedAlleles last_emitted;

    /// Snarl hierarchy for symbolic collapsing, or null to compare alleles by sequence alone.
    const SnarlManager* symbolic_manager = nullptr;

    /// Whether to decompose a snarl into one record per difference block.
    bool atomize_blocks = false;

    /// Which resolve pass will settle the site being recorded right now: its depth in the nested
    /// tree, since a chain's ploidy depends on its parent's settled genotype.
    ///
    /// Thread-local, and it has to be. Descent is synchronous on the calling thread but the sweep is
    /// parallel over node-ID windows, so as a plain member the save-increment-restore around each
    /// recursive call races between threads and drifts upward without bound: chr20 reported 124
    /// generations for a tree six deep, ran 124 resolve passes instead of 7, and clamped the linkage
    /// layer against generations that do not exist.
    static thread_local size_t current_generation;

    /// Resolve one generation of the linkage pass, accumulating into the three maps below.
    ///
    /// `last` gates the reporting, the mosaic, and building the mosaic -- all of which want
    /// the whole accumulated call set rather than one generation of it.
    void resolve_linkage_generation(size_t generation, bool last);

    /// Build the patch index and write the mosaic, once every record exists.
    ///
    /// Separate from resolution because both read whether a site has a VCF line, and with the record
    /// built after the decision that is unknown while resolving -- every entry still says unemitted.
    void finalise_linkage_outputs();

    /// Resolve the linkage pass if nothing has resolved it yet, as a single generation-0 pass.
    /// Idempotent, so `write_variants` can call it unconditionally and a driver that already ran the
    /// generations itself is not undone.
    void resolve_linkage();

    /// Accumulated across generations, consumed by `write_variants` when it patches the buffer.
    /// Several records can legitimately share a (contig, position) -- a nested child whose first
    /// variant base coincides with its parent's anchor base, for one -- so each position holds every
    /// call made there and the patch is matched to its line by record_key (the hash of the ID
    /// column). A map keyed by position alone silently kept only the last writer, so one of the two
    /// records lost its phasing and, when their GT strings coincided, the survivor was applied to
    /// both lines.
    /// Keyed on (contig, record_key), not (contig, POS).
    ///
    /// POS is the wrong key and has now cost two bugs. `flatten_common_allele_ends` advances POS by
    /// the prefix every allele shares, so the position a record is filed under is a function of its
    /// emitted allele set -- and anything that changes that set, or that computes the position before
    /// flattening, files the patch where nothing looks for it. `respecify` had to be given a contig
    /// and position for exactly this reason, and moving `record()` earlier would have reintroduced it
    /// from the other side. The failure is also near-invisible: `phase_declined`
    /// count patches that were found and refused, so a patch nobody looks up reports zero.
    ///
    /// `record_key` is `hash(print_snarl(snarl, false))` and `id_key()` hashes the ID column, which
    /// is that same string, so the two already agree. The contig stays in the key because one snarl
    /// ID can appear on more than one contig in a multi-reference run.
    /// Every phased site, in the order the model produced them. The mosaic reads this, and deferred
    /// descent looks a parent's settled allele pair up in it.
    vector<LinkageCollector::PhaseCall> linkage_phased;
    bool linkage_resolved = false;
    /// Set while the barrier re-emits, so the second pass does not add a second linkage entry for a
    /// site the barrier has already respecified.
    bool suppress_linkage_record = false;
    /// Totals for the one-line report, summed over however many generations ran.
    double linkage_seconds = 0.0;
    size_t linkage_changed = 0;

    /// Linkage pass state. Not owned.
    LinkageCollector* linkage_collector = nullptr;
    const gbwt::GBWT* linkage_gbwt = nullptr;
    const vector<size_t>* linkage_sequence_to_haplotype = nullptr;

    /// One decompressed-record cache per thread, for the panel lookups.
    ///
    /// Profiling put `panel_alleles` at roughly three and a half times the self cost of
    /// `score_read_against_allele` -- the per-read likelihood inner loop it exists to annotate --
    /// and three of its five heaviest frames were GBWT record lookup and decompression rather
    /// than the DA sampling intrinsic to `locate`. Every site re-walks records from scratch, and
    /// adjacent snarls share most of them.
    ///
    /// Per thread, not shared: `CachedGBWT` is a mutable cache with no internal locking and the
    /// caller is parallel over snarls. Sized once in `set_linkage` so nothing allocates inside
    /// the parallel region.
    mutable vector<gbwt::CachedGBWT> linkage_gbwt_cache;

    /// The node-ID neighbourhood each thread's cache currently covers. CachedGBWT only ever grows
    /// -- the gbwt header recommends short-lived instances -- and with node-ID-ordered windows a
    /// thread never revisits a retired window, so without eviction every decompressed record the
    /// contig ever touched stays resident. When a thread's site moves past this anchor by more than
    /// the fetch-window width, its cache is cleared: adjacent snarls still share records, which is
    /// all the cache was measured to help with, and residency is bounded to about one window.
    mutable vector<nid_t> linkage_gbwt_cache_anchor;

    /// Patches declined specifically because they named an allele beyond the record's ALT list, as
    /// opposed to declined for describing a genotype the record does not carry. The two mean
    /// different things -- the first is the traversal/VCF numbering gap, the second is a patch that
    /// found the wrong record -- so they are counted apart.

    /// Quality rewrites the record refused -- a malformed FORMAT, essentially. Counted rather than
    /// silent, because a record that keeps its per-site GQ where the model moved its genotype is
    /// differently calibrated and nothing else would say so.
    mutable std::atomic<size_t> quality_declined{0};

    /// Phase by record key, built after the barrier and read while each record is rendered.
    ///
    /// Keyed on the record key rather than on (contig, POS) because POS is not an identity here: the
    /// flattened position depends on which alleles the line carries, so it is not known until the
    /// record exists.
    std::unordered_map<size_t, LinkageCollector::PhaseCall> render_phases;

    /// Records phased while being rendered, and phases refused because the record did not carry a
    /// permutation of the phased pair.
    mutable std::atomic<size_t> phased_records{0};
    mutable std::atomic<size_t> phase_declined{0};

    /// Fill `render_phases` from the resolved phasing. Called between the barrier and the render.
    void build_render_phases();


    ///
    /// Returns false when the line is left untouched: the genotype being replaced is not the one
    /// this change was derived from, or the replacement names an allele the record has no ALT for.
    /// Refusing is always safe -- the line keeps its per-site call, which is a real answer -- and
    /// the count is reported, because a patch that silently does nothing is how the linkage layer
    /// once dropped every haploid correction it made.
    /// Rewrite a rendered record's quality from the linkage posterior: GQ becomes the phred
    /// complement, discounted by the explained-read share and capped at GQI, GQN is blanked and a
    /// stale `lowconf` cleared. The genotype is untouched -- the line already carries the settled
    /// one, because it was built from it.
    bool apply_linkage_quality(string& line, double posterior, double explained_share) const;

    /// Rewrite one emitted line's GT into phased form and attach its phase set. Applied after
    /// `apply_linkage_change`, so it phases the genotype that is actually emitted. Returns false
    /// when the line is left untouched, for the same reasons.

    /// Whether to emit phased GT and FORMAT/PS.
    bool emit_phasing = false;

    /// Destination for the mosaic file, and the graph it is to be read against.
    string mosaic_path;
    string mosaic_graph_name;
    /// Panel index -> "sample#phase", which is the unit the linkage model works in: a haplotype
    /// present in several GBWT fragments is one haplotype, so the name deliberately identifies a
    /// (sample, phase) pair rather than any single GBWT path. Combined with a segment's contig
    /// that is enough to find the paths again.
    ///
    /// The index itself is internal -- it is assigned in GBWT metadata order by the run that
    /// produced the file and means nothing outside it -- so the header emits the whole mapping and
    /// the file is self-describing.
    vector<string> mosaic_haplotype_names;

    /// Full reference path names the run called against; see set_mosaic_out.
    vector<string> mosaic_reference_paths;

    /// GBWT position of `hap`'s fragment at node `node_id`, or gbwt::invalid_edge() if the
    /// haplotype does not traverse it.
    ///
    /// This is what makes a segment walkable without an r-index. A consumer given only a node and
    /// a haplotype *name* has to turn the name into a position, and that is `locate()`: without an
    /// r-index or dense DA sampling it walks back to a sample, or scans the path from its start.
    /// Given the position instead, reconstruction is `extract({node, offset})` and forward `LF()`,
    /// with no search at all. Resolving it once here costs the writer a few thousand lookups and
    /// saves every reader an index larger than the graph -- this project's chr20 r-index is 86 MB
    /// against a 77 MB GBZ.
    ///
    /// Both orientations are tried because the recorded node IDs are bare: the snarl boundary is
    /// stored without the orientation the haplotype traverses it in.
    gbwt::edge_type mosaic_gbwt_position(int64_t node_id, size_t hap) const;

    /// Collapse the per-site phasing into segments and write them.
    ///
    /// Run-length encoding is the whole point: measured on chr20 against a 34-haplotype panel,
    /// about 2% of sites are switch points, so 105k sites become roughly 2k segments and the file
    /// is ~140 KB rather than ~7 MB. Whole-genome that is a few megabytes against a few hundred.
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

    /// Merge near-identical called ALT alleles in an already-populated variant.  Must run AFTER
    /// SnarlCaller::update_vcf_info and after flatten_common_allele_ends, so that the genotyper and
    /// the allele-flattening both see the full pre-merge allele set: merging earlier drops the
    /// absorbed allele's reads from AD/DP and from the Poisson caller's total_other_support term.
    /// Rewrites the allele-indexed fields (alleles/alt, AT, AD, GL, GT, MAD) and records what was
    /// merged in the MAT info field.  Returns true if anything merged.
    /// `gl_layout` must be the order the caller that produced this record's GL actually wrote it
    /// in. There is no way to recover it from the record.
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

    /// Decompose a fully built site record into one record per difference block and file them.
    ///
    /// Returns the number of lines written, or -1 meaning "the site record stands" -- which is
    /// every case this declines, so declining reproduces today's output exactly rather than
    /// approximately. `site` must be the finished record, after update_vcf_info and flattening,
    /// because every field a block does not redefine is inherited from it.
    int emit_block_records(const PathPositionHandleGraph& graph, const Snarl& snarl,
                           const vector<SnarlTraversal>& called_traversals,
                           const vector<int>& genotype, int ref_trav_idx,
                           const string& sample_name, const vcflib::Variant& site,
                           const map<int, int>& trav_to_allele, int64_t site_position,
                           GLLayout gl_layout, bool genotype_snarls) const;

    /// print a snarl in a consistent form like >3435<12222
    /// if in_brackets set to true,  do (>3435<12222) instead (this is only used for nested caller)
    string print_snarl(const HandleGraph* grpah, const handle_t& snarl_start, const handle_t& snarl_end, bool in_brackets = false) const;
    /// legacy version of above
    string print_snarl(const Snarl& snarl, bool in_brackets = false) const;
    /// The same as above, but print the snarl as if its orientation has been flipped
    string print_flipped_snarl(const Snarl& snarl, bool in_brackets = false) const;

    /// The linkage layer's identity for a site.
    ///
    /// It has to be the hash of the printed snarl and not something cheaper, because there is a
    /// seventh producer that cannot follow a change: `write_variants` recovers a buffered line's
    /// key by re-hashing the line's own ID column, which `emit_variant` set from `print_snarl`.
    /// That is what lets a compressed record be matched to its linkage entry with nothing carried
    /// alongside it, and it is what makes the key survive `--translation`, where both sides print
    /// the translated form.
    ///
    /// Packing the two boundary node IDs directly was tried and is about 0.15 s faster on chr20.
    /// It also takes that run's 19,472 linkage-quality patches to zero, because the key in the
    /// file cannot be repacked. Kept as one function instead, so the six call sites and the
    /// recovery in `write_variants` are a single statement apart rather than six copies of an
    /// expression that has to agree with a string.
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

    /// print up to this many uncalled alleles when doing ref-genotpes in -a mode
    size_t max_uncalled_alleles = 5;

    /// Contig name -> ploidy overrides, sorted by start and guaranteed non-overlapping by
    /// set_ploidy_regions. Empty unless --ploidy-bed was given. See set_ploidy_regions.
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
 * FlowCaller : Uses any traversals finder (ex, FlowTraversalFinder) to find
 * traversals, and calls those based on how much support they have.
 * Should work on any graph but will not
 * report cyclic traversals.  Supports nested calling when enabled with the
 * nested flag, recursively processing child snarls.
 * Designed to replace LegacyCaller, as it should miss fewer obviously
 * good traversals, and is not dependent on old protobuf-based structures.
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

    /// Defer symbolic descent until the linkage pass has settled each parent's genotype.
    ///
    /// Descent decides a child's ploidy, its strand, and whether it is called at all from its
    /// parent's genotype. Inline, that genotype is the *pre-linkage* one, and linkage then rewrites
    /// parents -- so a child can end up at a ploidy its own parent contradicts, and a child the
    /// settled parent does carry can go uncalled because the pre-linkage parent did not. Deferring
    /// makes the decision from the settled genotype instead, which removes both by construction.
    ///
    /// Greedy: a parent is settled before any of its children exists, so a child's evidence cannot
    /// move its parent. That is the trade, and it is deliberate.
    ///
    /// Needs the linkage layer and phasing, since the settled allele pair comes from the phasing.
    /// Sizes the per-thread queues, so call it before calling starts.
    /// Emit the records staged during the sweep, in one pass.
    ///
    /// Stage 9 of planning/decide-then-render.md, and deliberately the smallest step that proves the
    /// mechanism: the pass runs immediately after the sweep and BEFORE the barrier, so `record()` is
    /// still called before the linkage layer resolves and needs no decoupling. Output must not move
    /// at all, which is only checkable because the buffer sort became total in `abba6f288`.
    ///
    /// Stage 10 moves this pass after `resolve_linkage` and renders from the settled genotype; that
    /// is when `record()` has to leave `emit_variant`, since the barrier would otherwise resolve an
    /// empty collector.
    /// Record the site into the linkage layer at the point it is genotyped, not the point it is
    /// emitted.
    ///
    /// The barrier resolves the collector, and the render pass runs after the barrier, so a
    /// `record()` inside `emit_variant` would leave the collector empty when the barrier reads it --
    /// every nested revision, retraction and gain would disappear. Two of the arguments it used to
    /// take cannot exist here (the emitted allele map, and whether a line was written) and are
    /// supplied afterwards by `set_allele_map`.
    /// The pre-flatten (contig, position) a site is filed under in the linkage layer.
    pair<string, size_t> site_ref_key(const Snarl& snarl, const string& ref_path_name,
                                      int ref_offset, bool no_reference = false,
                                      int64_t anchor_position = 0) const;

    void record_site(const Snarl& snarl, const vector<SnarlTraversal>& travs,
                     const vector<int>& trav_genotype,
                     const unique_ptr<SnarlCaller::CallInfo>& call_info,
                     const string& ref_path_name, int ref_offset,
                     bool no_reference = false, int64_t anchor_position = 0);

    void render_retained_records();

    void set_defer_nested_descent(bool defer);

    /// How many render records are staged. Reported under --progress; the assertion that this is
    /// non-zero on a graph with no nesting at all is what proves the top-level branch stages too.
    size_t render_record_count() const;


    /// Resolve, descend, repeat until nothing is queued. One linkage pass per level of the snarl
    /// tree that descent actually reaches -- six on chr20, with 99.6% of the descents in the first
    /// three. Does nothing unless deferral is on, and leaves the linkage pass resolved either way,
    /// so `write_variants` needs no knowledge of which mode ran.
    void run_deferred_descent();


    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

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

    /// until we support nested snarls, cap snarl size we attempt to process
    size_t max_snarl_edges = 10000;

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

    /// enable recursive calling of child snarls
    bool nested = false;

    /// use * alleles for spanning haplotypes that don't traverse nested sites
    bool star_allele = false;

    /// One nested chain's genotyping, kept until the barrier can say what ploidy it should have.
    ///
    /// This is what makes the single sweep possible. A chain's ploidy comes from its parent's settled
    /// genotype, which is not known while the reads are resident, and the previous design went back to
    /// the reads once per level to find out -- five sweeps of a contig where there had been one, and
    /// 48.8% more reads fetched. Keeping the genotyping result instead costs 3.18 kB a chain,
    /// measured: 87 MB for chr20's 27,404 chains against a 3.6 GB peak.
    ///
    /// The `CallInfo` carries both ploidies' answers (see `alt_ploidy_info`), so whichever ploidy the
    /// barrier settles on can be rendered from what is here, with no re-reading and no re-scoring.
    /// `snarl` is held by value because `call_snarl_internal` works on a local copy it may have
    /// flipped.
    struct PendingRecord {
        Snarl snarl;
        string ref_path_name;
        int ref_offset = 0;
        vector<SnarlTraversal> travs;
        int ref_trav_idx = -1;
        /// The genotype and the ploidy it was genotyped at: the parent's *pre-linkage* answer.
        vector<int> genotype;
        int ploidy = 2;
        unique_ptr<SnarlCaller::CallInfo> call_info;
        size_t record_key = 0;
        size_t parent_record_key = 0;
        /// The chain's column in the alignment of its parent's two settled traversals, computed at
        /// the barrier, or -1. Carried here because the barrier can compute it before the chain has
        /// a layer entry to write it to, and `record` then takes it as an argument.
        int align_rank = -1;
        /// This snarl's index within its chain, from the decomposition. Unlike `align_rank` this is
        /// known at descent, so it needs no replay -- it is carried only because the barrier's
        /// record() fallback builds the entry and must be given it.
        /// Identity of this snarl's chain, from its boundary pair; the decode's group key.
        size_t chain_key = 0;
        int chain_index = -1;
        /// This snarl has no reference path of its own, so no line may be written for it: REF and
        /// POS are undefined. It is genotyped and recorded into the linkage layer regardless, which
        /// is the whole point -- the reads reach it through node-ID ranges, not coordinates.
        bool no_reference = false;
        /// The parent's reference start, standing in for a position this snarl does not have. Used
        /// wherever a key or an ordering coordinate is needed; never as a distance.
        int64_t anchor_position = 0;
        /// Whether the parent's settled traversal crossed this snarl's chain backward, from
        /// `SymbolicStep::backward` on the matched chain symbol. A barrier-time fact like
        /// `align_rank`, so it travels the same way.
        bool chain_backward = false;
        /// One bit per parent *candidate traversal*, set where that traversal crosses this chain.
        uint64_t parent_crossing = 0;
        /// False when parent_crossing could not be computed (the parent emitted nothing during the
        /// sweep, or has too many alleles for the mask). The barrier recomputes the mask when it
        /// re-emits the parent; until then an unknown mask must not be read as "nothing crosses".
        bool crossing_known = true;
        /// The parent traversal this chain hangs off, or -1 when the settled parent does not carry
        /// it on exactly one. Set at descent and re-derived whenever the barrier looks at the chain,
        /// so it always names a traversal in the parent's current settled pair.
        uint8_t generation = 0;
        /// Set when the settled parent turns out not to carry this chain, directly or through an
        /// ancestor. Such a chain does not exist in the sample, so neither do its descendants, and
        /// nothing below it may be revised or emitted afterwards.
        bool dropped = false;
        /// Whether a record was written during the sweep. False where no called parent allele reached
        /// it, which is exactly the population the barrier may turn into a call.
        bool emitted = false;
        /// Where the sweep buffered this record's line (see EmittedAlleles::buffer_thread), so the
        /// barrier can retract or replace exactly that line. -1 when nothing was buffered.
        int buffer_thread = -1;
        size_t buffer_index = 0;
    };

    bool defer_nested_descent = false;

    /// Filled per thread during the sweep, which is parallel over node-ID windows.
    vector<vector<PendingRecord>> pending_records;

    /// The same retention, for snarls the barrier does not revise -- top-level ones, whose ploidy
    /// comes from the contig or the BED and is never re-derived from a parent.
    ///
    /// A second container rather than one, because `run_deferred_descent` MOVES `pending_records`
    /// out and clears it: anything left there is gone by the time a render pass could read it. And
    /// because its `children_of` index buckets by `parent_record_key`, so every top-level record
    /// would land under key 0 and be walked as a sibling set.
    ///
    /// Nothing reads this yet. It exists so the memory cost of building the record after the
    /// genotype is decided is paid, and measured, in a commit whose output is byte-identical --
    /// which is the one risk in that change a reviewer cannot check by reading.
    vector<vector<PendingRecord>> render_records;


    /// Stage the inputs a record could be rendered from, for a snarl the barrier will not revise.
    /// Consumes `call_info` but deliberately NOT the traversals: descent still reads those. The
    /// caller completes the record after descent, exactly as it does for a nested chain.
    unique_ptr<PendingRecord> stage_render_record(const Snarl& snarl,
                                                 const vector<int>& trav_genotype, int ref_trav_idx,
                                                 unique_ptr<SnarlCaller::CallInfo>& call_info,
                                                 const string& ref_path_name, int ref_offset,
                                                 int ploidy, bool emitted);


    /// How many nested chains were retained across all threads.
    size_t pending_record_count() const;


    /// Internal implementation of call_snarl that accepts parent context for nested mode
    /// When nested=true, this recursively calls children after processing the current snarl
    /// @param parent_ref_path_name Reference path from parent (for off-reference snarls)
    /// @param parent_ref_interval Reference interval from parent
    /// @param parent_child_trav_sets If non-null, contains one TraversalSet per parent allele.
    ///                               Each set contains all traversals through this child that are
    ///                               consistent with that parent allele. The child genotypes by
    ///                               picking the best pair (one from each set) based on read support.
    /// `ploidy_override` >= 0 forces the ploidy this snarl is genotyped at, instead of taking it
    /// from the contig or the --ploidy-bed region. Symbolic nested calling uses it to hand a child
    /// the number of called parent alleles that actually cross it: a chain crossed by one allele
    /// and deleted by the other is haploid there, whatever the contig's ploidy says.
    bool call_snarl_internal(const Snarl& snarl,
                             const string& parent_ref_path_name,
                             pair<size_t, size_t> parent_ref_interval,
                             const ChildTraversalSets* parent_child_trav_sets = nullptr,
                             int ploidy_override = -1);

    /// How many of the called parent alleles cross this child snarl, capped at `cap`.
    ///
    /// Crossing means the child's start and end both appear in the parent traversal *in order*.
    /// Testing for them independently, as find_child_traversal_set does, counts a traversal that
    /// happens to touch both boundaries on unrelated excursions.
    ///
    /// A single allele crossing a chain more than once -- a cycle or a tandem duplication -- would
    /// give a per-haplotype copy number above one. v1 caps rather than modelling that, and says so
    /// in the log, because the rest of the caller assumes ploidy in {1, 2}.
    int child_ploidy(const vector<SnarlTraversal>& travs, const vector<int>& genotype,
                     const Snarl& child, int cap) const;

    /// How many times one traversal crosses `child`, by the same in-order rule child_ploidy uses.
    static int crossings_of_child(const SnarlTraversal& trav, const Snarl& child);

    /// One bit per VCF allele, set where a traversal that became that allele crosses `child`.
    ///
    /// Returns 0 -- unknown -- if any allele index is beyond a 64-bit mask, and sets `*known` to
    /// false there, so a caller can tell that 0 apart from "no allele crosses" instead of silently
    /// conflating the two. Several traversals can share an allele index; symbolic collapsing puts
    /// every same-route traversal on allele 0, and those agree on their crossings by construction,
    /// so any disagreement can only come from two distinct routes that spell the same sequence.
    /// The bit is set if any of them crosses.
    static uint64_t child_crossing_mask(const vector<SnarlTraversal>& travs,
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
