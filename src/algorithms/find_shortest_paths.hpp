#ifndef VG_ALGORITHMS_FIND_SHORTEST_PATHS_HPP_INCLUDED
#define VG_ALGORITHMS_FIND_SHORTEST_PATHS_HPP_INCLUDED

/**
 * \file find_shortest_paths.hpp
 *
 * Definitions for the find_shortest_paths algorithm.
 */

#include <unordered_map>

#include "../position.hpp"
#include "../vg.pb.h"
#include "../hash_map.hpp"
#include "../handle.hpp"

namespace vg {
namespace algorithms {
    
    /// Finds the length of the shortest oriented path from the given handle
    /// leftward to all reachable oriented nodes on a directed walk. Uses
    /// Dijkstra's Algorithm.
    unordered_map<handle_t, size_t>  find_shortest_paths(const HandleGraph* g, handle_t start);
                                                      
}
}

#endif
