#ifndef VG_ALGORITHMS_FIND_CLOSEST_WITH_PATHS_HPP_INCLUDED
#define VG_ALGORITHMS_FIND_CLOSEST_WITH_PATHS_HPP_INCLUDED

/**
 * \file find_closest_with_paths.hpp
 *
 * Defines algorithm for finding nearby nodes that have paths on then in a PathHandleGraph.
 */

#include "../handle.hpp"

#include <vector>
#include <utility>

#include "structures/rank_pairing_heap.hpp"

namespace vg {
namespace algorithms {

using namespace std;
    
    /// Return a vector of pairs of handles that occur on the same relative
    /// strand as the start handle, the distance from the right or left end
    /// of the start handle needed to reach them, and whether they were reached
    /// going right (true) or left (false) from the start. The handles are the
    /// closest one(s) to the start handle that touch any paths.
    ///
    /// Distances are always positive. If specified, right_extra_dist and
    /// left_extra_dist are added to the search distances, to allow for
    /// searching from a particular point on the starting node. 
    ///
    /// Search does not exceed max_search_dist bases.
    ///
    /// Will only ever return an empty vector or a 1-element vector.
    vector<tuple<handle_t, int64_t, bool>> find_closest_with_paths(const PathHandleGraph& graph, handle_t start, size_t offset, int64_t max_search_dist);
    

}
}

#endif
