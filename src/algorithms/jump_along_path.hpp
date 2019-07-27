#ifndef VG_ALGORITHMS_JUMP_ALONG_PATH_HPP_INCLUDED
#define VG_ALGORITHMS_JUMP_ALONG_PATH_HPP_INCLUDED

/**
 * \file jump_along_path.hpp
 *
 * Defines algorithm for jumping a given number of bases away from a position
 * using paths.
 */

#include "../handle.hpp"

#include "structures/rank_pairing_heap.hpp"

namespace vg {
namespace algorithms {

using namespace std;
    
    /// returns a vector of positions that are found by jumping a fixed oriented distance
    /// along path(s) from the given position. if the position is not on a path, searches
    /// from the position to a path and adds/subtracts the search distance to the jump
    /// depending on the search direction. returns an empty vector if there is no path within
    /// the max search distance or if the jump distance goes past the end of the path
    vector<pos_t> jump_along_closest_path(const PathPositionHandleGraph* graph,
                                          const pos_t& pos, int64_t jump_dist,
                                          size_t max_search_dist);

}
}

#endif
