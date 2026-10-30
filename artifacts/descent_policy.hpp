#ifndef VG_DESCENT_POLICY_HPP_INCLUDED
#define VG_DESCENT_POLICY_HPP_INCLUDED

#include <array>
#include <atomic>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

#include "block_records.hpp"
#include "child_crossing_finder.hpp"
#include "../src/handle.hpp"
#include "site_values.hpp"
#include <vg/vg.pb.h>

namespace vg::multipass {

using namespace std;

/// The tallies of what nested descent did in one run: how deep it went, and how many child chains
/// it skipped or recorded, and why.
struct DescentCounters {
    /// The number of children placed at each level, for levels below 16.
    std::array<std::atomic<size_t>, 16> depth_hist{};
    /// Children skipped because the reference walk does not cross them and off-reference descent
    /// is off.
    std::atomic<size_t> skipped_no_ref{0};
    /// Children that the reference walk does not cross, descended into because off-reference
    /// descent is on.
    std::atomic<size_t> off_reference{0};
    /// Sites that no reference path passes through, given an entry in the linkage model but no
    /// line.
    std::atomic<size_t> no_ref_recorded{0};
    /// Sites with no reference path, by number of called copies: 0, 1, and 2 or more.
    std::array<std::atomic<size_t>, 3> no_ref_copies{};
    /// Children that no called parent allele crosses in the direct pass.
    ///
    /// With the linkage model they are genotyped and kept, and the linkage pass decides from the
    /// parent's chosen genotype whether the sample has them.
    std::atomic<size_t> skipped_no_copy{0};
    /// Times a called traversal crossed a child more than once.
    ///
    /// Only the first crossing counts, for ploidy and for distance. One traversal crossing a chain
    /// twice is one strand with two copies, so it does not make the chain ploidy 2.
    std::atomic<size_t> child_multi_crossing{0};
};

/**
 * The place of a site in the nesting tree, which its parent's genotype gives it.
 *
 * A top-level site has the default placement.
 */
struct NestingPlacement {
    /// Whether exactly one of the parent's called alleles crosses this chain (see
    /// `LinkageCollector::SiteContext::nested`).
    bool one_copy = false;
    /// The record key of the parent.
    size_t parent_record_key = 0;
    /// The level of the site among the nested chains, 0 for a top-level site.
    ///
    /// It decides when in a linkage pass the site's genotype is chosen.
    size_t level = 0;
    /**
     * The crossing mask, with one bit per candidate traversal of the parent, set where that
     * traversal crosses this chain.
     *
     * It is indexed by traversal, because the linkage model chooses the parent's genotype as a
     * pair of traversals. The linkage pass reads it to find how many copies of the chain the
     * parent's chosen genotype carries, and on which strand. It is 0 when `crossing_known` is
     * false.
     */
    uint64_t parent_crossing = 0;
    /**
     * Whether no called parent allele crosses the chain, under the linkage model.
     *
     * Without the linkage model, such a chain is not genotyped. With it, the chain is genotyped
     * anyway, at the parent's ploidy, because the linkage model may still move the parent onto an
     * allele that crosses it. Its children inherit the flag, and are genotyped at their own
     * ploidy.
     *
     * In the direct pass the chain is staged, but not written or filed in the linkage model. A
     * site whose own bounds are on no reference path is filed anyway, and never gets a line.
     *
     * The linkage pass decides the fate of the chain from the parent's chosen pair:
     * - If the pair carries no copy, or the crossing mask is 0, it drops the chain and everything
     *   under it.
     * - If the pair carries copies, it files the chain at that many copies and renders it. If the
     *   direct pass scored no genotype at that ploidy, it renders the chain at the parent's ploidy
     *   instead, unfiled.
     * - If `crossing_known` is false, it cannot compare the chain with the pair. It renders the
     *   chain at the ploidy it was genotyped at, unfiled, and drops it only with an ancestor.
     *
     * A parent that the linkage model gave no phase call is read at its own chosen genotype.
     */
    bool retain_only = false;
    /**
     * The offset of the chain along its parent, in bases.
     *
     * It is where the chain starts along the first of the parent's called traversals that crosses
     * it, plus the parent's own `parent_offset`. Added to the parent's reference start, it gives
     * an off-reference chain a position of its own. Then its sites are ordered as that traversal
     * visits them, and the distance between two of them is known. The linkage pass computes it
     * again from the parent's chosen genotype (`StagedSite::placement`).
     */
    size_t parent_offset = 0;
    /// Whether the chain is off the reference, and so may be genotyped with no line.
    ///
    /// Children inherit it, since everything under such a chain is off the reference too. The
    /// candidate finder still looks for a reference path through each site, since a descendant's
    /// bounds may lie on a reference path even when its parent's do not.
    bool no_reference = false;
    /// Whether, under block emission, an enclosing block's ALT already spells the variation in
    /// this chain.
    ///
    /// The chain is still genotyped and filed, and its line is held back when records are
    /// rendered. Children inherit it.
    bool reported_inline = false;
    /// The key of the chain (`chain_key_of`).
    ///
    /// The linkage model groups the sites of a chain by it, and puts no transitions between
    /// chains under one parent.
    size_t chain_key = 0;
    /// Whether `parent_crossing` is known.
    ///
    /// It is unknown when the parent has more than 64 candidate traversals, too many for the mask.
    bool crossing_known = true;
};

/**
 * The rules for nested descent: which child chains of a genotyped site to genotype next, and how.
 *
 * It lists the child chains that a genotyped site's called alleles cross, each with its place in
 * the nesting tree and the ploidy to genotype it at. It reads a site's children from the snarl
 * decomposition, and finds where the site's walks cross them with `ChildCrossingFinder`.
 */
class DescentPolicy {
public:
    /// A child chain to genotype, as `place` gives it.
    struct Placed {
        /// The child, as its net handle in the overlay of sites.
        net_handle_t site;
        /// The chain the child belongs to.
        ChildChain chain;
        /// The ploidy to genotype it at.
        ///
        /// It is the number of the parent's called alleles that cross it, or the parent's ploidy
        /// when none does, which is the most copies a child can have.
        int ploidy = 0;
        /// The place of the child in the nesting tree.
        NestingPlacement placement;
    };

    /// Place children in `graph`, as `decomposition` nests them.
    ///
    /// Test each against the difference blocks of `blocks`, and count what descent does in
    /// `counters`. None is owned.
    void configure(const HandleGraph* graph, const EditSnarlDecompositionOverlay* decomposition,
                   const BlockRecordWriter* blocks, DescentCounters* counters);

    /// The sites nested in a site.
    struct Nested {
        /// What the records of the site need to know of its children.
        SiteChildren children;
        /// The child sites, in the order of the decomposition.
        vector<net_handle_t> sites;
    };

    /// Get the sites nested in `site`, which the caller has with start bound `called_start` and
    /// end bound `called_end`: read the way the decomposition orients it, or turned round.
    Nested nested_in(const net_handle_t& site, const handle_t& called_start,
                     const handle_t& called_end) const;

    /// Find the place of node `id` in the decomposition (see `placement_of_node`).
    NodePlacement placement_of(nid_t id) const;

    /// A genotyped site, as the policy needs it to place the site's children.
    struct Parent {
        /// The site, as its net handle in the overlay of sites.
        net_handle_t site;
        /// The sites nested in it.
        const Nested& nested;
        /// Its record key.
        size_t record_key;
        /// Its candidate walks.
        const vector<Traversal>& travs;
        /// Its called genotype, as indices into `travs`.
        const vector<int>& genotype;
        /// The index in `travs` of the reference walk, or -1 for none.
        int ref_trav_idx;
        /// The ploidy it was called at.
        int ploidy;
        /// Its place in the nesting tree.
        const NestingPlacement& placement;
    };

    /**
     * Call `visit` with each child chain to genotype under `parent`, in the order of the
     * decomposition.
     *
     * Place each child just before its visit, so the placement work and the genotyping of
     * children interleave.
     *
     * @param off_reference Genotype a chain that the reference walk does not cross, with no line.
     * Without it, leave such a chain out.
     * @param keep_uncrossed Genotype a chain that no called allele crosses, at the parent's
     * ploidy, and retain it, since the linkage model may still move the parent onto an allele
     * that crosses it. Without it, leave such a chain out.
     */
    void place(const Parent& parent, bool off_reference, bool keep_uncrossed,
               const function<void(const Placed& child)>& visit) const;

    /// Count the called parent alleles that cross `child`, capped at `cap`.
    ///
    /// An allele crossing a chain more than once, as in a cycle or tandem duplication, counts
    /// once, since the caller assumes ploidy 1 or 2. Such alleles are tallied in
    /// `DescentCounters::child_multi_crossing`.
    int child_ploidy(const vector<ChildCrossingFinder>& walks, const vector<int>& genotype,
                     const net_handle_t& child, int cap) const;

private:
    /// Return true if the reference walk of `parent`, if it has one, crosses `child`.
    bool reference_crosses(const Parent& parent, const vector<ChildCrossingFinder>& walks,
                           const net_handle_t& child) const;

    const HandleGraph* graph = nullptr;
    /// Gives the children of each site.
    const EditSnarlDecompositionOverlay* decomposition = nullptr;
    /// Says which child chains a difference block already spells.
    const BlockRecordWriter* blocks = nullptr;
    /// Tallies what descent does.
    DescentCounters* counters = nullptr;
};

}

#endif
