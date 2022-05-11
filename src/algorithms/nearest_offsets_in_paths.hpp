#ifndef VG_ALGORITHMS_FIND_CLOSEST_WITH_PATHS_HPP_INCLUDED
#define VG_ALGORITHMS_FIND_CLOSEST_WITH_PATHS_HPP_INCLUDED

/**
 * \file nearest_offsets_in_paths.hpp
 *
 * Defines algorithm for finding the nearest offset along a path of the
 * the closest position in the graph that overlaps a path
 */

#include "../handle.hpp"

#include <vector>
#include <utility>
#include <unordered_map>
#include <map>

#include "structures/rank_pairing_heap.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// Represents a set of positions and orientations, along a collection of paths.
/// Positions and orientations may or may not be stored sorted.
using path_offset_collection_t = unordered_map<path_handle_t, vector<pair<size_t, bool>>>;
    
/// Return, for the nearest position in a path to the given position,
/// subject to the given max search distance, a mapping from path name to
/// all positions on each path where that pos_t occurs.
/// Stops search when path(s) are ancountered.
///
/// If path_filter is set, ignores paths for which it returns false.
path_offset_collection_t nearest_offsets_in_paths(const PathPositionHandleGraph* graph,
                                                  const pos_t& pos, int64_t max_search,
                                                  const std::function<bool(const path_handle_t&)>* path_filter = nullptr);
    
/// Wrapper for the above to support some earlier code. Only looks for paths
/// that directly touch the position, and returns the paths by name.
map<string, vector<pair<size_t, bool>>> offsets_in_paths(const PathPositionHandleGraph* graph, const pos_t& pos);

/// A "simple" model for path position getting for debugging
path_offset_collection_t simple_offsets_in_paths(const PathPositionHandleGraph* graph, pos_t pos);
    

}
}

#endif
