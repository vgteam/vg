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

/// Returns sets of handles defining components that are connected by any series
/// of nodes and edges, even if it is not a valid bidirected walk. All handles
/// returned are in the local forward orientation.
vector<unordered_set<handle_t>> weakly_connected_components(const HandleGraph* graph);


}
}

#endif
