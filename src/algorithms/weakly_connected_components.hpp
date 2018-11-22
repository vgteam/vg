#ifndef VG_ALGORITHMS_WEAKLY_CONNECTED_COMPONENTS_HPP_INCLUDED
#define VG_ALGORITHMS_WEAKLY_CONNECTED_COMPONENTS_HPP_INCLUDED

/**
 * \file weakly_connected_components.hpp
 *
 * Defines an algorithm for finding weakly connected components in a graph.
 */

#include "../handle.hpp"

#include <unordered_set>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

/// Returns sets of IDs defining components that are connected by any series
/// of nodes and edges, even if it is not a valid bidirected walk. TODO: It
/// might make sense to have a handle-returning version, but the consumers of
/// weakly connected components right now want IDs, and membership in a weakly
/// connected component is orientation-independent.
vector<unordered_set<id_t>> weakly_connected_components(const HandleGraph* graph);

/// Return pairs of weakly connected component ID sets and the handles that are
/// their tips, oriented inward. If a node is both a head and a tail, it will
/// appear in tips in both orientations.
vector<pair<unordered_set<id_t>, vector<handle_t>>> weakly_connected_components_with_tips(const HandleGraph* graph);

}
}

#endif
