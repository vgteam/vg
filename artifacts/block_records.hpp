#ifndef VG_BLOCK_RECORDS_HPP_INCLUDED
#define VG_BLOCK_RECORDS_HPP_INCLUDED

#include <array>
#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "../src/handle.hpp"
#include "../src/snarls.hpp"
#include "string_aligner.hpp"
#include "symbolic_allele.hpp"
#include "vcf_genotype_likelihoods.hpp"
#include "vcf_record.hpp"

namespace vg::multipass {

using namespace std;

/// Why `BlockRecordWriter::write` declined a site, so that the site's single record is written
/// instead. The value indexes `AtomizeCounters::refuse` and the report's list of reasons.
enum class AtomizeRefusal : size_t {
    NoGenotype,
    NoReferenceTraversal,
    Unresolvable,
    EmptyReferenceProjection,
    StepMapLength,
    NoDifferenceBlocks,
    NoAnchorBase,
    AnchorWithoutSequence,
    ReferenceBasesOnly,
    SameAsSiteRecord,
    AllelesMerged,
    StrandsDisagree,
    ChainCrossedTwice,
    BoundChanged,
    /// The number of reasons.
    Count
};

/// Counters for block emission, for the report.
struct AtomizeCounters {
    /// Sites that reached `BlockRecordWriter::count_site`, so that the report can tell "nothing
    /// refused" from "never ran".
    std::atomic<size_t> sites{0};
    /// Sites whose snarl does not resolve, so that their projection has no chain symbols.
    std::atomic<size_t> site_unresolvable{0};
    /// Sites that resolve only as the reversal of a snarl.
    std::atomic<size_t> site_reversed{0};
    /// Sites written as blocks.
    std::atomic<size_t> split_sites{0};
    /// Lines written for sites written as blocks.
    std::atomic<size_t> split_lines{0};
    /// Chains whose own record is not written because the ALT of a block already spells them.
    std::atomic<size_t> child_inlined{0};
    /// Sites `BlockRecordWriter::write` declined, by reason.
    std::array<std::atomic<size_t>, (size_t)AtomizeRefusal::Count> refuse{};
};

/// One place where a strand of a run of sites differs from the run's reference: a longest run of
/// columns of their alignment that are not matches.
struct FlatDifference {
    /// The offset in the reference of the first base the strand replaces.
    size_t ref_begin;
    /// The offset in the reference after the last base the strand replaces, which is `ref_begin`
    /// where the strand only inserts bases.
    size_t ref_end;
    /// The bases the strand has in their place, which are none where it only deletes bases.
    string alt;

    /// Return whether a record of the difference alone needs the reference base before it, as VCF
    /// requires where the reference or the strand has no bases.
    bool needs_anchor() const { return ref_begin == ref_end || alt.empty(); }
};

/// One flat record of a run, before it is written: a difference, and the strands that have it.
struct FlatCut {
    FlatDifference difference;
    /// Whether each strand has the difference.
    vector<bool> on_strand;

    /// Return the offset in the reference of the first base the record covers, which is the
    /// base before the difference where a record of it needs one.
    int begin() const {
        return (int)difference.ref_begin - (difference.needs_anchor() ? 1 : 0);
    }
};

/// One strand of a run of sites, aligned base by base to the run's reference.
struct FlatStrand {
    /// The strand's sequence over the run (see `AlleleSplicer::sequence_along`).
    string sequence;
    /// Where the strand differs from the reference, in reference order. No two touch.
    vector<FlatDifference> differences;
};

/**
 * The strands of a run of neighbouring top-level sites, in the order of the sites' phase, aligned
 * to the run's reference. A site on its own is a run of one. Flat records are cut from it (see
 * `BlockRecordWriter::write_flat`).
 */
struct FlatRun {
    /// The ID that the run's records share, before their indices: that of the run's one site,
    /// or the ID of a site from the start bound of its first site to the end bound of its last.
    string id;
    /// The reference: what the reference allele of each site and the nodes between the sites
    /// spell, from the first base of the first site's reference allele.
    string reference;
    /// The strands, strand 0 first.
    vector<FlatStrand> strands;
    /// The run's records, in reference order (see `cut_differences`).
    vector<FlatCut> cuts;
    /// Whether the called genotype of a site of the run, or of a site nested in one, has an
    /// allele that is a consensus of reads.
    bool consensus_called = false;

    /**
     * Set `cuts` from the differences of the strands.
     *
     * A difference that several strands have, with the same reference bases and the same bases
     * in their place, is one record of those strands. Every other difference is a record of its
     * own strand, whatever the other strands have at its reference bases.
     */
    void cut_differences();

    /// The records that the sites of the run have written, and the number of sites still to
    /// write theirs, for the check of the run's records. A site adds its records under the lock.
    mutable std::mutex check_lock;
    mutable vector<vcflib::Variant> check_records;
    mutable size_t check_waiting = 0;
};

/// The part of the flat records of a run that one site of the run writes: the records that
/// start in its reference allele, or after it and before the next site of the run.
struct FlatSite {
    /// The run. Its sites share it.
    shared_ptr<const FlatRun> run;
    /// The offset in the run's reference of the first base of the site's reference allele,
    /// which is where the site's record is placed.
    size_t origin = 0;
    /// The traversal of the site that each strand carries. A flat record that the site writes
    /// takes its AD and GL from the site allele of the traversal of each strand.
    vector<int> travs;
    /// The index in `run->cuts` of the first record the site writes, and their number.
    size_t first_cut = 0;
    size_t cut_count = 0;
};

/// Why a site is not written as flat records, so that it and the sites nested in it are written
/// as they are without `--flatten-records`, or why a run of several sites is not written as one,
/// so that each of its sites is written as a run of its own. The value indexes
/// `FlatCounters::refuse` and `FlatCounters::run_refuse`, and the report's list of reasons.
enum class FlatRefusal : size_t {
    NoGenotype,
    NoPhase,
    PhaseDiffers,
    ReferenceDiffers,
    AlignerLimit,
    NoAnchorBase,
    /// A strand cannot be spelled along the run.
    Unspelled,
    /// The number of reasons.
    Count
};

/// Counters for flat records, for the report.
struct FlatCounters {
    /// Runs written as flat records, and the sites in them.
    std::atomic<size_t> runs{0};
    std::atomic<size_t> sites{0};
    /// Of the runs, those of two or more sites.
    std::atomic<size_t> long_runs{0};
    /// Flat records written.
    std::atomic<size_t> lines{0};
    /// Runs whose strands all spell the reference, which write nothing.
    std::atomic<size_t> reference_runs{0};
    /// Nested sites with no line because the flat records of an enclosing site spell them.
    std::atomic<size_t> withheld{0};
    /// Edit sites in the nodes between two sites of a run, with no line because the run's flat
    /// records spell them.
    std::atomic<size_t> withheld_between{0};
    /// Runs whose flat records were checked against their strands.
    std::atomic<size_t> checked{0};
    /// Checked runs whose flat records do not spell their strands.
    std::atomic<size_t> check_failed{0};
    /// Sites not written as flat records, by reason.
    std::array<std::atomic<size_t>, (size_t)FlatRefusal::Count> refuse{};
    /// Runs of two or more sites written as a run for each site, by reason.
    std::array<std::atomic<size_t>, (size_t)FlatRefusal::Count> run_refuse{};
};

/**
 * A writer that writes a site as block records instead of one site record, and decides which
 * child chains those blocks already report.
 *
 * A block record covers a cluster of difference blocks between the reference and the symbolic
 * alleles of the called strands. Deciding which chains the blocks report lets each variant be
 * reported once.
 *
 * Block emission needs nested calling, without which no site has child chains to project. The
 * writer counts what it does, for the report. The counters are atomic, so the const methods can
 * be called from many threads.
 */
class BlockRecordWriter {
public:
    /// Set whether nested calling is on. Off turns block emission and the site counts off.
    void set_nested(bool nested_on) { this->nested = nested_on; }

    /// Turn block emission on or off. Off by default.
    void set_enabled(bool enabled_on) { this->enabled = enabled_on; }

    /// Return whether block emission was turned on, for the VCF header.
    bool is_enabled() const { return enabled; }

    /**
     * The parts of the test of `chain_reported_inline` that are the same for every child of a
     * site.
     *
     * These are the projection of the reference, the projection of each called ALT, and the
     * difference blocks between them. They are built once per site, since aligning them is the
     * costly part.
     */
    struct ChainInlineContext {
        /// Whether the test can say yes for any child. It cannot when block emission is off, the
        /// genotype is empty, the site does not resolve, or a strand carries the reference.
        bool usable = false;
        /// The projection of the reference traversal.
        SymbolicAllele sref;
        /// A called ALT, aligned to the reference.
        struct Alt {
            /// The projection of the ALT.
            SymbolicAllele salt;
            /// The difference blocks between `sref` and `salt`.
            vector<DiffBlock> blocks;
        };
        /// The called ALTs, in genotype order.
        vector<Alt> alts;
    };

    /// Build the parts of the inline test that are the same for every child of a site with
    /// children `children` and genotype `genotype` over `travs`.
    ChainInlineContext chain_inline_context(const HandleGraph& graph,
                                            const SiteChildren& children,
                                            const vector<Traversal>& travs,
                                            const vector<int>& genotype,
                                            int ref_trav_idx) const;

    /// Return whether the block records of the site already report the child chain `chain`.
    ///
    /// They do when every called strand crosses the chain only inside a difference block, whose
    /// ALT spells the route through it. A chain the reference does not cross is for the caller to
    /// handle.
    bool chain_reported_inline(const HandleGraph& graph, const ChainInlineContext& ctx,
                               const ChildChain& chain) const;

    /// Restart the count of chains reported inline, before a pass decides them all again.
    void restart_inline_count() { counters.child_inlined = 0; }

    /// Count a site with children `children` for the report, before its record is built, by
    /// whether the decomposition knows it.
    void count_site(const SiteChildren& children, const vector<Traversal>& travs,
                    int ref_trav_idx) const;

    /**
     * Write a site as block records, giving each to `add_line` with its index. Return the number
     * of lines `add_line` accepted, or -1 if the site record is to be written instead.
     *
     * The site record is written instead when block emission is off or declines the site. The
     * lines are written for `sample_name`, with nodes named by `translation` if it is not null.
     *
     * `record` must be the finished site record, after the fields of the snarl caller are written
     * and the alleles are trimmed, since each block record copies every field it does not
     * redefine. A site whose alleles were merged is declined, since it then numbers its alleles
     * differently from `record.trav_to_allele`, and its blocks would spell the merged alleles
     * apart.
     *
     * @throws std::out_of_range if an entry of `genotype` is neither an allele marker nor the
     * number of a traversal in `called_traversals`.
     */
    int write(const PathPositionHandleGraph& graph, const SiteChildren& children,
              const vector<Traversal>& called_traversals, const vector<int>& genotype,
              int ref_trav_idx, const string& sample_name, const NodeTranslation* translation,
              const SiteRecord& record, GLLayout gl_layout, bool genotype_snarls,
              const function<bool(vcflib::Variant&, size_t)>& add_line,
              const AlleleSequences* sequences = nullptr) const;

    /// Set the costs and limits of the aligner that `align_flat` uses.
    void set_flat_aligner(const StringAlignerParams& params) { flat_aligner = params; }

    /// Set whether `write_flat` checks each run's records with `flat_records_spell`, once every
    /// site of the run has written its own, and counts the result for the report.
    void set_flat_check(bool on) { flat_check = on; }

    /**
     * Align each strand of a run of sites to the run's reference `reference`.
     *
     * `strands` gives the sequence of each strand of the run, in the order of the phase. Each
     * strand is aligned by itself at the least cost (see `align_global`), and its gaps are
     * moved as far left as they go without leaving the first base of the run.
     *
     * @param refused Set to why, where nothing is returned.
     * @return Nothing if the aligner reaches its limit, or if a difference that needs the
     * reference base before it starts at the first base of the run.
     */
    shared_ptr<FlatRun> align_flat(string reference, vector<string> strands,
                                   FlatRefusal& refused) const;

    /// Count a site that is not written as flat records for `reason`.
    void count_flat_refusal(FlatRefusal reason) const { ++flat_counters.refuse[(size_t)reason]; }

    /// Count a run of two or more sites that is not written as one run for `reason`.
    void count_flat_run_refusal(FlatRefusal reason) const {
        ++flat_counters.run_refuse[(size_t)reason];
    }

    /// Count a run of `sites` sites that is to be written as flat records.
    void count_flat_run(const FlatRun& run, size_t sites) const;

    /// Count the sites left without a line because the flat records of a run spell them:
    /// `nested` sites nested in a site of a run, and `between` edit sites in the nodes between
    /// two sites of a run.
    void count_flat_withheld(size_t nested, size_t between) const {
        flat_counters.withheld += nested;
        flat_counters.withheld_between += between;
    }

    /**
     * Write the flat records of a run that the site of `flat` writes, giving each to `add_line`
     * with its index, and return the number of lines `add_line` accepted.
     *
     * A flat record is one difference of a strand from the run's reference. A difference is a
     * longest run of bases a strand does not match: a run of mismatches is one record, and so
     * is a gap, with the reference base before it, or a gap with the mismatches that touch it.
     * The GT lists the strands in order. A strand that has the difference has allele 1, so a
     * difference on both strands is homozygous. A strand that lacks it has allele 0, whatever
     * it has at those reference bases: another difference of that strand there is a record of
     * its own, so two records of a run can cover the same reference bases. A run whose strands
     * all spell the reference writes nothing.
     *
     * A record is placed by its offset in the reference from the site's first base, which is
     * at `record.unflattened_position`. It copies every field of `record`, the finished record
     * of the site, that it does not redefine, as a block record does: it takes the ID of the
     * run with its index among the run's records appended, INFO/SB with that index and the
     * number of the run's records, the site's PS where the site's GT is phased, and AD and GL
     * through the site allele of each strand's traversal. It has INFO/FL, and "." as the AT of
     * each allele, since no walk spells its alleles.
     */
    int write_flat(const PathPositionHandleGraph& graph, const FlatSite& flat,
                   const string& sample_name, const NodeTranslation* translation,
                   const SiteRecord& record, GLLayout gl_layout,
                   const function<bool(vcflib::Variant&, size_t)>& add_line) const;

    /**
     * Return whether `records`, in position order, spell every strand of `run`.
     *
     * Each strand is rebuilt from the reference of `run`, whose first base is at position
     * `first_position`, by putting in place of the REF of each record whose GT gives the strand
     * an ALT that ALT. The entries of the GT of `sample_name` are taken as the strands in
     * order, whichever separator it has.
     */
    static bool flat_records_spell(const FlatRun& run, const vector<vcflib::Variant>& records,
                                   int64_t first_position, const string& sample_name);

    /// Print the counters to stderr. Print nothing when no site was counted.
    void report() const;

private:
    /// Count a site declined for `reason`, and return -1, which `write` returns for it.
    int refuse(AtomizeRefusal reason) const;

    /// Print the declined sites to stderr, by reason, if there are any.
    void report_refusals() const;

    /// Whether nested calling is on.
    bool nested = false;
    /// Whether block emission is on.
    bool enabled = false;
    /// What the writer has done, for the report.
    mutable AtomizeCounters counters;
    /// The costs and limits of the aligner that `align_flat` uses. A mismatch costs 40, and a gap
    /// of L bases the lesser of 78 + 7L and 162 + 3L. Those are the scores minimap2 aligns two
    /// assemblies of one species with (asm5: match 1, mismatch -19, gaps -39 - 3L and -81 - L),
    /// as the costs of an alignment from end to end in which a match costs nothing, doubled so
    /// that they are whole numbers.
    StringAlignerParams flat_aligner{.mismatch = 40,
                                     .gap_open = 78,
                                     .gap_extend = 7,
                                     .long_gap_open = 162,
                                     .long_gap_extend = 3};
    /// Whether `write_flat` checks the records it writes.
    bool flat_check = false;
    /// What the writer has done with flat records, for the report.
    mutable FlatCounters flat_counters;

    /// Print the counters of flat records to stderr, if any site was counted.
    void report_flat() const;
};

}

#endif
