#ifndef VG_CHILD_PLACER_HPP_INCLUDED
#define VG_CHILD_PLACER_HPP_INCLUDED

#include <atomic>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "handle.hpp"
#include <vg/vg.pb.h>

namespace vg {

using namespace std;

class BlockRecordWriter;
class SnarlManager;

/// Counts of what nested descent did in one run: how deep it went, and how many child chains it
/// skipped or recorded, and why. One per run, so that runs count separately.
struct DescentCounters {
    /// How deep the descent went, by depth.
    std::atomic<size_t> depth_hist[16] = {};
    /// Children skipped because no reference path passes through them and off-reference descent
    /// is off.
    std::atomic<size_t> skipped_no_ref{0};
    /// Children that no reference path passes through, descended into because off-reference
    /// descent is on.
    std::atomic<size_t> off_reference{0};
    /// Sites that no reference path passes through, given an entry in the linkage model but no
    /// line.
    std::atomic<size_t> no_ref_recorded{0};
    /// Off-reference chains by copy number: 0, 1, 2.
    std::atomic<size_t> no_ref_copies[3] = {};
    /// Children that no called parent allele crosses in the direct pass. They are genotyped and kept,
    /// and the linkage pass decides from the parent's chosen genotype whether the sample has them.
    std::atomic<size_t> skipped_no_copy{0};
    /// Children that a called traversal enters more than once. Only the first entry counts, both
    /// for ploidy and for distance: one traversal crossing a chain twice is one strand carrying two
    /// copies, not two strands, so it does not make the chain ploidy 2.
    std::atomic<size_t> child_multi_crossing{0};
};

/**
 * Where a site sits in the nesting tree: the place its parent's genotype gives it. A top-level
 * site has the default placement.
 */
struct NestingPlacement {
    /// Exactly one of the parent's called alleles crosses this chain (see
    /// `LinkageCollector::SiteContext::nested`).
    bool one_copy = false;
    size_t parent_record_key = 0;
    /// The site's level: its depth among the nested chains, which decides when in a linkage pass
    /// its genotype is chosen. 0 for a top-level site.
    size_t level = 0;
    /// The crossing mask: one bit per parent candidate traversal, set where that traversal
    /// crosses this chain. It is indexed by traversal rather than by VCF allele, since the
    /// linkage model chooses the parent's genotype as a pair of traversals, and the VCF alleles are
    /// chosen only when the parent's line is written. When the mask cannot be computed (more
    /// than 64 traversals), it is 0 and `crossing_known` is false. The linkage pass reads it to find
    /// how many copies of the chain the parent's chosen genotype carries, and on which
    /// strand.
    uint64_t parent_crossing = 0;
    /// Set where no called parent allele reaches the chain, and only when staging is on (see
    /// `FlowCaller::set_stage_records`) and the linkage model runs (without it, such a chain is
    /// not genotyped). The chain is genotyped anyway, at the parent's ploidy, because the linkage
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
    /// (`StagedSite::chain_offset`).
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

/**
 * Nested descent's placement: lists the child chains a genotyped site's called alleles cross,
 * each with its place in the nesting tree and the ploidy to genotype it at.
 *
 * Its static members say how a parent's traversals cross its child chains: which traversals
 * cross a child, and how far along a traversal the child starts. A traversal crosses a child when
 * it visits one of the child's boundary nodes and then the other, so a traversal that touches
 * both boundaries on unrelated excursions does not count. Visits to child snarls are left out.
 */
class ChildPlacer {
public:
    /// Where each node ID is visited in one traversal, in ascending order, excluding visits to
    /// snarls. Built once per traversal per snarl, so that testing each child does not scan the
    /// whole traversal again.
    using TraversalNodeIndex = unordered_map<nid_t, vector<int>>;
    static TraversalNodeIndex index_traversal_nodes(const SnarlTraversal& trav);

    /// How many times one traversal crosses `child`. One allele crossing a chain more than once,
    /// as in a cycle or tandem duplication, crosses it that many times.
    static int crossings_of_child(const TraversalNodeIndex& visits, const Snarl& child);

    /// Where along `trav` the child chain is first entered, as a visit index, or -1 if `trav` does
    /// not cross it, by the rule `crossings_of_child` uses.
    static int offset_of_child(const SnarlTraversal& trav, const Snarl& child);

    /// How far along one traversal, in bases, each child is entered: the total length of the
    /// nodes visited before it, or -1 if the traversal does not cross it. Answered by lookup, so
    /// that a parent with many children and a long traversal is not scanned once per child.
    struct ChildOffsets {
        ChildOffsets(const HandleGraph& graph, const SnarlTraversal& trav);
        /// The bases before the visit `offset_of_child(trav, child)` gives, or -1 where that is -1.
        int64_t base_offset(const Snarl& child) const;
        /// The visit indices of each node, ascending. Child-snarl visits are left out, as
        /// `offset_of_child` skips them.
        unordered_map<nid_t, vector<int>> visits_of;
        /// The bases of the node visits before each visit index; one longer than the traversal.
        vector<int64_t> bases_before;
    };

    /// How far along the first traversal of `genotype` that crosses `child` the child is entered,
    /// in bases, or 0 when none does. `offsets` holds a `ChildOffsets` per traversal of `travs`,
    /// filled as traversals are first used.
    static size_t offset_along_genotype(const HandleGraph& graph,
                                        const vector<SnarlTraversal>& travs,
                                        const vector<int>& genotype, const Snarl& child,
                                        unordered_map<const SnarlTraversal*, ChildOffsets>& offsets);

    /// The crossing mask: bit i is set where the traversal `visits[i]` indexes crosses `child`.
    /// Indexed by traversal, not by VCF allele, since it is tested against the parent's chosen
    /// traversals. Returns 0 and sets `*known` to false when there are more than 64 traversals,
    /// so that the caller can tell unknown from "no traversal crosses".
    static uint64_t child_crossing_mask(const vector<TraversalNodeIndex>& visits,
                                        const Snarl& child, bool* known = nullptr);

    /// A child chain to genotype, as `place` lists it.
    struct Placed {
        /// The child, as the snarl manager holds it.
        const Snarl* snarl = nullptr;
        /// The ploidy to genotype it at: how many of the parent's called alleles cross it, or the
        /// parent's ploidy when none does, the most copies a child can have.
        int ploidy = 0;
        NestingPlacement placement;
    };

    /// Place children in `graph`, as `manager` nests them, testing each against the difference
    /// blocks of `blocks`, and counting what descent does in `counters`. None is owned.
    void configure(const HandleGraph* graph, const SnarlManager* manager,
                   const BlockRecordWriter* blocks, DescentCounters* counters);

    /// The child chains to genotype under a genotyped site, in the manager's order, each placed
    /// under the site. The site is `site`, named `site_key`, with candidate traversals `travs`,
    /// the reference among them at `ref_trav_idx`, called at `genotype` and `ploidy`, and placed
    /// at `placement`.
    ///
    /// A chain the reference does not cross is left out, unless `off_reference` lets descent
    /// genotype it with no line. A chain no called allele crosses is left out, unless
    /// `keep_uncrossed`: then it is genotyped at the parent's ploidy and retained, since the
    /// linkage model may still move the parent onto an allele that crosses it.
    vector<Placed> place(const Snarl& site, size_t site_key, const vector<SnarlTraversal>& travs,
                         const vector<int>& genotype, int ref_trav_idx, int ploidy,
                         const NestingPlacement& placement, bool off_reference,
                         bool keep_uncrossed) const;

    /// How many of the called parent alleles cross this child snarl, capped at `cap`.
    ///
    /// One allele crossing a chain more than once, as in a cycle or tandem duplication, counts
    /// once, since the caller assumes ploidy 1 or 2; this is counted.
    int child_ploidy(const vector<TraversalNodeIndex>& visits, const vector<int>& genotype,
                     const Snarl& child, int cap) const;

    /// How far along `trav`, in bases, the child chain is entered: the total length of the nodes
    /// visited before it, or -1 if `trav` does not cross it, by the rule `offset_of_child` uses.
    /// It gives an off-reference chain its place along its parent (see
    /// `NestingPlacement::parent_offset`).
    int64_t base_offset_of_child(const SnarlTraversal& trav, const Snarl& child) const;

    /// `base_offset_of_child` along the first traversal of `genotype` that crosses `child`, or 0
    /// when none does.
    size_t offset_along_genotype(const vector<SnarlTraversal>& travs, const vector<int>& genotype,
                                 const Snarl& child) const;

private:
    const HandleGraph* graph = nullptr;
    const SnarlManager* manager = nullptr;
    const BlockRecordWriter* blocks = nullptr;
    DescentCounters* counters = nullptr;
};

}

#endif
