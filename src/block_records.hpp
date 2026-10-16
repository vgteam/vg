#ifndef VG_BLOCK_RECORDS_HPP_INCLUDED
#define VG_BLOCK_RECORDS_HPP_INCLUDED

#include <atomic>
#include <functional>
#include <string>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"
#include "symbolic_allele.hpp"
#include "vcf_genotype_likelihoods.hpp"
#include "vcf_record.hpp"

namespace vg {
namespace multipass {

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
    std::atomic<size_t> site_unresolvable{0};  // flip_snarl left projection with no symbols
    std::atomic<size_t> site_reversed{0};      // resolved only via the reversed pairing
    /// Sites written as blocks, and the lines they produced.
    std::atomic<size_t> split_sites{0}, split_lines{0};
    /// Chains whose own record is not written because a block's ALT already spells them.
    std::atomic<size_t> child_inlined{0};
    /// Sites `BlockRecordWriter::write` declined, by reason.
    std::atomic<size_t> refuse[(size_t)AtomizeRefusal::Count] = {};
};

/**
 * Writes a site as one record per difference block between the reference and each called
 * strand's symbolic allele, instead of one record for the whole site, and decides which child
 * chains those blocks already report, so that each variant is reported once.
 *
 * Configured with whether nested calling is on, without which no site has child chains to
 * project, and whether block emission is on. Each site's children come with the site. Counts what
 * it does, for the report; the counters are atomic, so the const methods can be called from many
 * threads.
 */
class BlockRecordWriter {
public:
    /// Whether nested calling is on. Off turns block emission and the site counts off.
    void set_nested(bool nested) { this->nested = nested; }

    /// Turn block emission on or off. Off by default.
    void set_enabled(bool enabled) { this->enabled = enabled; }

    /// Whether block emission was turned on, for the VCF header.
    bool is_enabled() const { return enabled; }

    /// The parts of the inline test that do not depend on the child: the site's symbolic
    /// projection, each called ALT's projection, and the difference blocks between the reference
    /// and each ALT. Built once per site rather than once per child, since the edit-distance
    /// alignment in `symbolic_diff` is the same for every child.
    struct ChainInlineContext {
        /// False when the answer is false for every child: block emission off, indices out of
        /// range, an empty genotype, an unresolvable site, or the reference among the called
        /// alleles.
        bool usable = false;
        SymbolicAllele sref;
        struct Alt {
            SymbolicAllele salt;
            vector<DiffBlock> blocks;
        };
        /// One entry per called allele that is in range and not the reference, in genotype order.
        vector<Alt> alts;
    };

    /// Build the child-independent half of the inline test for a site in `graph` with children
    /// `children` and genotype `genotype` over `travs`. See ChainInlineContext.
    ChainInlineContext chain_inline_context(const HandleGraph& graph,
                                            const SiteChildren& children,
                                            const vector<Traversal>& travs,
                                            const vector<int>& genotype,
                                            int ref_trav_idx) const;

    /// Whether the child chain `chain` is already reported by the site's own block records,
    /// because every called strand crosses it only inside a difference block whose ALT spells the
    /// route through it. This can happen only when no called allele is the reference allele; a
    /// chain that no reference path passes through is handled separately by the caller.
    bool chain_reported_inline(const HandleGraph& graph, const ChainInlineContext& ctx,
                               const ChildChain& chain) const;

    /// Restart the count of chains reported inline, before a pass decides them all again.
    void restart_inline_count() { counters.child_inlined = 0; }

    /// Count a site with children `children` for the report, before its record is built, by
    /// whether the decomposition knows it. It changes no output.
    void count_site(const SiteChildren& children, const vector<Traversal>& travs,
                    int ref_trav_idx) const;

    /// Write a site as its difference blocks, giving each block line to `add_line` with its block
    /// index. The lines are written for `sample_name`, with nodes named by `translation` if it is
    /// not null. Returns the number of lines `add_line` accepted, or -1 when block emission is off or
    /// declines the site, in which case the site record is to be written as it is.
    ///
    /// `record` must be the finished site record, after the snarl caller's fields are written and
    /// the alleles are flattened, since every field a block does not redefine is taken from it.
    /// A site whose alleles were merged is declined, since it then numbers its alleles differently
    /// from `record.trav_to_allele`, and its blocks would spell the merged alleles apart.
    int write(const PathPositionHandleGraph& graph, const SiteChildren& children,
              const vector<Traversal>& called_traversals, const vector<int>& genotype,
              int ref_trav_idx, const string& sample_name, const NodeTranslation* translation,
              const SiteRecord& record, GLLayout gl_layout, bool genotype_snarls,
              const function<bool(vcflib::Variant&, size_t)>& add_line,
              const AlleleSequences* sequences = nullptr) const;

    /// Print the counters to stderr. Prints nothing when no site was counted.
    void report() const;

private:
    bool nested = false;
    bool enabled = false;
    mutable AtomizeCounters counters;
};

}
}

#endif
