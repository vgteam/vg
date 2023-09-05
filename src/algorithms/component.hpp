/** \file
 * Contains algorithms that report per-component info
 */

#ifndef VG_ALGORITHMS_COMPONENT_HPP_INCLUDED
#define VG_ALGORITHMS_COMPONENT_HPP_INCLUDED

#include <structures/rank_pairing_heap.hpp>

#include "handle.hpp"

namespace vg {
namespace algorithms {

using namespace std;

// returns the number of weakly connected components
size_t num_components(const HandleGraph& graph);

// returns the size in number of nodes of each component
vector<size_t> component_sizes(const HandleGraph& graph);

// returns sets of path handles, one set for each component (unless the
// component doesn't have any paths)
vector<unordered_set<path_handle_t>> component_paths(const PathHandleGraph& graph);

// the same semantics as the previous, but multithreaded and more
// memory intensive
vector<unordered_set<path_handle_t>> component_paths_parallel(const PathHandleGraph& graph);
}

}

#endif // VG_ALGORITHMS_COMPONENT_HPP_INCLUDED
