/** \file
 * Contains an algorithm to split paths up by their connected component
 */

#ifndef VG_ALGORITHMS_COMPONENT_PATHS_HPP_INCLUDED
#define VG_ALGORITHMS_COMPONENT_PATHS_HPP_INCLUDED

#include <structures/rank_pairing_heap.hpp>

#include "handle.hpp"

namespace vg {
namespace algorithms {

using namespace std;

// returns sets of path handles, one set for each component (unless the
// component doesn't have any paths)
vector<unordered_set<path_handle_t>> component_paths(const PathHandleGraph& graph);

// the same semantics as the previous, but multithreaded and more
// memory intensive
vector<unordered_set<path_handle_t>> component_paths_parallele(const PathHandleGraph& graph);
}

}

#endif // VG_ALGORITHMS_COMPONENT_PATHS_HPP_INCLUDED
