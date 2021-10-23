/** \file
 * Measures the distance between two graph positions along the reference path
 * (approximated by the longest connecting path)
 */

#ifndef VG_ALGORITHMS_REF_PATH_DISTANCE_HPP_INCLUDED
#define VG_ALGORITHMS_REF_PATH_DISTANCE_HPP_INCLUDED

#include <structures/rank_pairing_heap.hpp>

#include "handle.hpp"
#include "position.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// Search the local region around two positions and return the longest distance between
/// them along any paths found during this search. Returns numeric_limits<int64_t>::max()
/// if no shared path is found.
int64_t ref_path_distance(const PathPositionHandleGraph* graph, const pos_t& pos_1, const pos_t& pos_2,
                          const unordered_set<path_handle_t>& ref_paths, int64_t max_search_dist);

}

}

#endif // VG_ALGORITHMS_REF_PATH_DISTANCE_HPP_INCLUDED
