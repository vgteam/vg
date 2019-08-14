#ifndef VG_ALGORITHMS_DIJKSTRA_HPP_INCLUDED
#define VG_ALGORITHMS_DIJKSTRA_HPP_INCLUDED

/**
 * \file dijkstra.hpp
 *
 * Definitions of Dijkstra's Algorithm utilities for traversing the graph closest-first.
 */

#include <unordered_map>
#include <vg/vg.pb.h>

#include "../position.hpp"
#include "../hash_map.hpp"
#include "../handle.hpp"

namespace vg {
namespace algorithms {
    
    /// Walk out from the given handle and visit all reachable handles
    /// (including the start) in the graph once, in closest-first order,
    /// accounting for sequence lengths. Walks right unless traverse_leftward
    /// is specified, in which case it walks left. Distances are measured
    /// between the outgoing side of the start node and the incoming side of
    /// the target.
    ///
    /// When the shortest distance to a handle has been determined, calls
    /// reached_callback with that handle and the distance to it. Calls to
    /// reached_callback will occur in ascending order of distance. The
    /// reached_callback function must return true to continue the search, and
    /// false to abort it early.
    ///
    /// Returns true if the search terminated normally, and false if it was
    /// aborted.
    bool dijkstra(const HandleGraph* g, handle_t start,
                  function<bool(const handle_t&, size_t)> reached_callback,
                  bool traverse_leftward = false);
                  
    /// Same as the single-start version, except allows starting from multiple
    /// handles, all at distance 0.
    bool dijkstra(const HandleGraph* g, const unordered_set<handle_t>& starts,
                  function<bool(const handle_t&, size_t)> reached_callback,
                  bool traverse_leftward = false);
                                                      
}
}

#endif
