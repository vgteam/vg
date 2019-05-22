/**
 * \file dijkstra.cpp
 *
 * Implementation of Dijkstra's Algorithm over the bidirected graph.
 */
 
#include "find_shortest_paths.hpp"
#include <structures/updateable_priority_queue.hpp>

namespace vg {
namespace algorithms {

using namespace structures;

bool dijkstra(const HandleGraph* g, handle_t start,
              function<bool(const handle_t&, size_t)> reached_callback,
              bool traverse_leftward) {

    // We keep a priority queue so we can visit the handle with the shortest
    // distance next. We put handles in here whenever we see them with shorter
    // distances (since STL priority queue can't update), so we also need to
    // make sure nodes coming out haven't been visited already.
    using Record = pair<size_t, handle_t>;
    
    // We filter out handles that have already been visited.
    unordered_set<handle_t> visited;
    
    // We need a custom ordering for the queue
    struct IsFirstGreater {
        inline bool operator()(const Record& a, const Record& b) {
            return a.first > b.first;
        }
    };
    
    // We use a filtered priority queue for auto-Dijkstra
    UpdateablePriorityQueue<Record, handle_t, vector<Record>, IsFirstGreater> queue([](const Record& item) {
        return item.second;
    });
    
    // We keep a current handle
    handle_t current = start;
    size_t distance = 0;
    queue.push(make_pair(distance, start));
    
    while (!queue.empty()) {
        // While there are things in the queue, get the first.
        tie(distance, current) = queue.top();
        queue.pop();
        
#ifdef debug_vg_algorithms
        cerr << "Visit " << g->get_id(current) << " " << g->get_is_reverse(current) << " at distance " << distance << endl;
#endif    


        // Emit the handle as being at the given distance
        if (!reached_callback(current, distance)) {
            // The user told us to stop. Return that we stopped early.
            return false;
        }
        
        // Remember that we made it here.
        visited.emplace(current);
        
        if (current != start) {
            // Up the distance with the node's length. We don't do this for the
            // start handle because we want to count distance from the *end* of
            // the start handle unless directed otherwise.
            distance += g->get_length(current);
        }
            
        g->follow_edges(current, traverse_leftward, [&](const handle_t& next) {
            // For each handle to the right of here
            
            if (!visited.count(next)) {
                // New shortest distance. Will never happen after the handle comes out of the queue because of Dijkstra.
                queue.push(make_pair(distance, next));
                
#ifdef debug_vg_algorithms
                cerr << "\tNew best path to " << g->get_id(next) << " " << g->get_is_reverse(next)
                    << " at distance " << distance << endl;
#endif
                
            }
        });
    }

    // If we made it here, we finished the entire graph.
    return true;

}

    
}
}
