/**
 * \file find_shortest_paths.cpp
 *
 * Implementation for the find_shortest_paths algorithm.
 */
 
#include "find_shortest_paths.hpp"
#include "dijkstra.hpp"
#include <structures/updateable_priority_queue.hpp>

namespace vg {
namespace algorithms {

using namespace structures;

unordered_map<handle_t, size_t>  find_shortest_paths(const HandleGraph* g, handle_t start,
                                                     bool traverse_leftward) {

    // This is the minimum distance to each handle
    unordered_map<handle_t, size_t> distances;
    
    dijkstra(g, start, [&](const handle_t& current, size_t distance) {
        // Record handle's distance
        distances[current] = distance;
        
        // Always keep going
        return true;
    }, traverse_leftward);
    
    return distances;

}

    
}
}
