#ifndef VG_CHILD_CROSSING_FINDER_HPP_INCLUDED
#define VG_CHILD_CROSSING_FINDER_HPP_INCLUDED

/** \file
 * A finder of the places where one walk through a site crosses the site's children.
 */

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "edit_snarl_decomposition_overlay.hpp"
#include "../src/handle.hpp"
#include "../src/traversal_finder.hpp"

namespace vg::multipass {

using namespace std;

/// One time a walk crosses a child site: where it enters the child.
struct ChildCrossing {
    /// The step at which the walk enters the child.
    int entry_step;
    /// The bases of the walk before the child starts.
    ///
    /// For a child inside nodes, this counts the bases of the entered node before the child, read
    /// the way the walk reads the node.
    int64_t entry_bases;
};

/**
 * A finder of the places where one walk through a site crosses the site's children.
 *
 * A child is a snarl or an edit of the overlay of sites, given by its net handle. A walk crosses a
 * child each time it passes one of the child's boundaries and later the other.
 * - For a child between two nodes, a visit to a boundary node passes the side of it that faces
 *   the child, whichever way the node is read. So the walk crosses the child when it visits one
 *   boundary node and then the other, and enters it at the first of the two visits. A child whose
 *   two bounds are one node, as in a looping chain, is crossed by a second visit to that node.
 * - For a child inside one node, a visit to the node passes both of its boundaries, so every
 *   visit crosses the child once.
 *
 * Make one per walk, so that finding the crossings of each child does not scan the whole walk
 * again. It refers to the overlay, the graph and the walk, which must outlive it.
 */
class ChildCrossingFinder {
public:
    /// Prepare to find where `walk`, a walk in `graph`, crosses children in `sites`.
    ChildCrossingFinder(const EditSnarlDecompositionOverlay& sites, const HandleGraph& graph,
                        const Traversal& walk);

    /// Find every time the walk crosses `child`, in walk order.
    ///
    /// A walk can cross a child more than once, as in a cycle or tandem duplication.
    vector<ChildCrossing> crossings(const net_handle_t& child) const;

private:
    /// Get the bases of the node of step `step` before the first base of `child`, an edit, read
    /// the way the walk reads the node.
    int64_t bases_before_child(const net_handle_t& child, int step) const;

    const EditSnarlDecompositionOverlay* sites;
    const HandleGraph* graph;
    const Traversal* walk;
    /// The steps of the walk that visit each node, in ascending order.
    unordered_map<nid_t, vector<int>> visits_of;
    /// The bases of the walk before each step. It is one longer than the walk.
    vector<int64_t> bases_before;
};

/**
 * Make the crossing mask of `child`, with bit i set where the walk of `walks[i]` crosses it.
 *
 * The mask is indexed by candidate walk, since it is tested against a parent's chosen walks.
 *
 * @param known If not null, set to false when there are more than 64 walks, and true otherwise.
 * @return The mask, or 0 when there are more than 64 walks.
 */
uint64_t crossing_mask(const vector<ChildCrossingFinder>& walks, const net_handle_t& child,
                       bool* known = nullptr);

/**
 * Find how far along the first walk of `genotype` that crosses `child` the walk first enters it,
 * in bases, or 0 when none does.
 *
 * @param travs The walks that `genotype` indexes.
 * @param finders The crossing finder of each walk of `travs`, filled as walks are first used.
 */
size_t entry_bases_along_genotype(const EditSnarlDecompositionOverlay& sites,
                                  const HandleGraph& graph, const vector<Traversal>& travs,
                                  const vector<int>& genotype, const net_handle_t& child,
                                  unordered_map<const Traversal*, ChildCrossingFinder>& finders);

}

#endif
