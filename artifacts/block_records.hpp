#ifndef VG_BLOCK_RECORDS_HPP_INCLUDED
#define VG_BLOCK_RECORDS_HPP_INCLUDED

#include <array>
#include <atomic>
#include <functional>
#include <string>
#include <vector>

#include "../src/handle.hpp"
#include "../src/snarls.hpp"
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
};

}

#endif
