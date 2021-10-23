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

vector<unordered_set<path_handle_t>> component_paths(const PathHandleGraph& graph);

}

}

#endif // VG_ALGORITHMS_COMPONENT_PATHS_HPP_INCLUDED
