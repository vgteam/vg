#ifndef VG_ALGORITHMS_NEXT_PATH_POSITION_HPP_INCLUDED
#define VG_ALGORITHMS_NEXT_PATH_POSITION_HPP_INCLUDED

/**
 * \file next_path_position.hpp
 *
 * Defines algorithm for finding the closest position on a path downstream from a handle.
 */

#include "../handle.hpp"
#include "../position.hpp"

#include <utility>

namespace vg {
namespace algorithms {

using namespace std;
    
    /// Get the nearest position that is on a path, if the distance between it
    /// and the current position is <= max_search.
    /// Will search over multiple nodes.
    /// Returns a pair of the position and the distance to it.
    /// If no such position is found, returns a position with node ID 0.
    pair<pos_t, int64_t> next_path_position(const PathHandleGraph& graph, pos_t pos, int64_t max_search);

}
}

#endif
