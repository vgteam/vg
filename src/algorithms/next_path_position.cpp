/**
 * \file next_path_position.cpp
 *
 * Contains implementation of next_path_position function that finds
 * the next downstream position on a path in a PathHandleGraph.
 */

#include "next_path_position.hpp"

#include "find_closest_with_paths.hpp"


//#define debug_algorithms

namespace vg {
namespace algorithms {

using namespace std;

pair<pos_t, int64_t> next_path_position(const PathHandleGraph& graph, pos_t pos, int64_t max_search) {
    // We need the closest on-a-path position to the passed in position, and the approximate rightward offset to it.
    // Offset will be negative if you have to go left instead.
    
    // Get a handle to where we start
    handle_t h = graph.get_handle(id(pos), is_rev(pos));
    
    // Find the closest on-a-path node
    vector<tuple<handle_t, int64_t, bool>> closest = algorithms::find_closest_with_paths(graph,
        h, offset(pos), max_search);
    
    if (!closest.empty()) {
        // We found something.
        // Unpack it.
        auto& found = get<0>(closest.front());
        auto& dist = get<1>(closest.front());
        auto& arrived_at_end = get<2>(closest.front());
        
        // Output the position. Make sure to get the offset we end up at in the
        // node's forward strand for the end of the node we arrive at.
        return make_pair(make_pos_t(graph.get_id(found),
                                    graph.get_is_reverse(found),
                                    arrived_at_end ? graph.get_length(found) : 0),
                         dist);
    } else {
        // We aren't connected to anything on a path within the requested distance limit
        return make_pair(make_pos_t(0,false,0), numeric_limits<int64_t>::max());
    }
}
    
    
}
}
