#ifndef VG_CHILD_PLACER_HPP_INCLUDED
#define VG_CHILD_PLACER_HPP_INCLUDED

#include <atomic>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

#include "block_records.hpp"
#include "handle.hpp"
#include "site_values.hpp"
#include <vg/vg.pb.h>

namespace vg {

using namespace std;

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
    /// Set where no called parent allele reaches the chain, and only when the linkage model runs
    /// (without it, such a chain is not genotyped). The chain is genotyped anyway, at the parent's ploidy, because the linkage
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
    /// Where each node ID is visited in one walk, in ascending order. Built once per walk per
    /// site, so that testing each child does not scan the whole walk again.
    using TraversalNodeIndex = unordered_map<nid_t, vector<int>>;
    static TraversalNodeIndex index_traversal_nodes(const HandleGraph& graph,
                                                    const Traversal& walk);

    /// How many times one walk crosses `child`. One allele crossing a chain more than once, as in
    /// a cycle or tandem duplication, crosses it that many times.
    static int crossings_of_child(const HandleGraph& graph, const TraversalNodeIndex& visits,
                                  const SiteBounds& child);

    /// Where along `walk` the child is first entered, as a handle index, or -1 if `walk` does not
    /// cross it, by the rule `crossings_of_child` uses.
    static int offset_of_child(const HandleGraph& graph, const Traversal& walk,
                               const SiteBounds& child);

    /// How far along one walk, in bases, each child is entered: the total length of the nodes
    /// visited before it, or -1 if the walk does not cross it. Answered by lookup, so that a
    /// parent with many children and a long walk is not scanned once per child.
    struct ChildOffsets {
        ChildOffsets(const HandleGraph& graph, const Traversal& walk);
        /// The bases before the handle `offset_of_child(graph, walk, child)` gives, or -1 where
        /// that is -1.
        int64_t base_offset(const HandleGraph& graph, const SiteBounds& child) const;
        /// The handle indices of each node, ascending.
        unordered_map<nid_t, vector<int>> visits_of;
        /// The bases of the handles before each handle index; one longer than the walk.
        vector<int64_t> bases_before;
    };

    /// How far along the first walk of `genotype` that crosses `child` the child is entered, in
    /// bases, or 0 when none does. `offsets` holds a `ChildOffsets` per walk of `travs`, filled as
    /// walks are first used.
    static size_t offset_along_genotype(
        const HandleGraph& graph, const vector<Traversal>& travs, const vector<int>& genotype,
        const SiteBounds& child, unordered_map<const Traversal*, ChildOffsets>& offsets);

    /// The crossing mask: bit i is set where the walk `visits[i]` indexes crosses `child`.
    /// Indexed by candidate walk, not by VCF allele, since it is tested against the parent's
    /// chosen walks. Returns 0 and sets `*known` to false when there are more than 64 walks, so
    /// that the caller can tell unknown from "no walk crosses".
    static uint64_t child_crossing_mask(const HandleGraph& graph,
                                        const vector<TraversalNodeIndex>& visits,
                                        const SiteBounds& child, bool* known = nullptr);

    /// A child chain to genotype, as `place` gives it.
    struct Placed {
        /// The child, as the snarl manager holds it.
        const Snarl* snarl = nullptr;
        /// The chain the child belongs to.
        ChildChain chain;
        /// The ploidy to genotype it at: how many of the parent's called alleles cross it, or the
        /// parent's ploidy when none does, the most copies a child can have.
        int ploidy = 0;
        NestingPlacement placement;
    };

    /// Place children in `graph`, as `manager` nests them, testing each against the difference
    /// blocks of `blocks`, and counting what descent does in `counters`. None is owned.
    void configure(const HandleGraph* graph, const SnarlManager* manager,
                   const BlockRecordWriter* blocks, DescentCounters* counters);

    /// Call `visit` with each child chain to genotype under a genotyped site, in the manager's
    /// order, each placed under the site. The site is `site`, with children `children`, named
    /// `site_key`, with candidate walks `travs`, the reference among them at `ref_trav_idx`, called
    /// at `genotype` and `ploidy`, and placed at `placement`. Each child is placed just before its visit, so the
    /// visits run in the order, and among the work, that descent has always genotyped children in.
    ///
    /// A chain the reference does not cross is left out, unless `off_reference` lets descent
    /// genotype it with no line. A chain no called allele crosses is left out, unless
    /// `keep_uncrossed`: then it is genotyped at the parent's ploidy and retained, since the
    /// linkage model may still move the parent onto an allele that crosses it.
    void place(const Snarl& site, const SiteChildren& children, size_t site_key,
               const vector<Traversal>& travs,
               const vector<int>& genotype, int ref_trav_idx, int ploidy,
               const NestingPlacement& placement, bool off_reference, bool keep_uncrossed,
               const function<void(const Placed& child)>& visit) const;

    /// How many of the called parent alleles cross this child snarl, capped at `cap`.
    ///
    /// One allele crossing a chain more than once, as in a cycle or tandem duplication, counts
    /// once, since the caller assumes ploidy 1 or 2; this is counted.
    int child_ploidy(const vector<TraversalNodeIndex>& visits, const vector<int>& genotype,
                     const SiteBounds& child, int cap) const;

    /// How far along `walk`, in bases, the child is entered: the total length of the nodes
    /// visited before it, or -1 if `walk` does not cross it, by the rule `offset_of_child` uses.
    /// It gives an off-reference chain its place along its parent (see
    /// `NestingPlacement::parent_offset`).
    int64_t base_offset_of_child(const Traversal& walk, const SiteBounds& child) const;

    /// `base_offset_of_child` along the first walk of `genotype` that crosses `child`, or 0 when
    /// none does.
    size_t offset_along_genotype(const vector<Traversal>& travs, const vector<int>& genotype,
                                 const SiteBounds& child) const;

private:
    /// Whether the site's own blocks already report `child`, a child snarl of the site, as
    /// `BlockRecordWriter::chain_reported_inline` decides for the child's chain.
    bool reported_inline(const BlockRecordWriter::ChainInlineContext& ctx,
                         const Snarl& child) const;

    const HandleGraph* graph = nullptr;
    const SnarlManager* manager = nullptr;
    const BlockRecordWriter* blocks = nullptr;
    DescentCounters* counters = nullptr;
};

}

#endif
