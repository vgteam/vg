#ifndef VG_ALGORITHMS_WEAKLY_CONNECTED_COMPONENTS_DENSE_HPP_INCLUDED
#define VG_ALGORITHMS_WEAKLY_CONNECTED_COMPONENTS_DENSE_HPP_INCLUDED

/** \file
 * Weakly connected components of a graph whose node IDs are dense.
 */

#include <unordered_set>
#include <vector>

#include "../handle.hpp"

namespace vg {
namespace algorithms {

/**
 * Return the weakly connected components of `graph`, each as the set of its node IDs.
 *
 * Search each component depth-first, starting from the first of its nodes in the graph's own
 * order, and mark the nodes reached in an array indexed by node ID. This gives the components,
 * and inserts the nodes in each, in the same order as
 * `handlealgs::weakly_connected_components()`, and runs much faster on a large graph, since that
 * function marks nodes in a hash set.
 *
 * The array has one entry per ID from the smallest node ID to the largest, so the node IDs should
 * be dense.
 *
 * TODO: Offer this in libhandlegraph alongside `handlealgs::weakly_connected_components()`.
 */
std::vector<std::unordered_set<nid_t>> weakly_connected_components_dense(const HandleGraph& graph);

}
}

#endif
