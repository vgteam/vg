/**
 * \file find_closest_with_paths.cpp
 *
 * Contains implementation of find_closest_with_paths function that finds
 * nearby nodes with paths in a PathHandleGraph.
 */

#include "find_closest_with_paths.hpp"

#include <unordered_set>
#include <queue>

//#define debug_algorithms

namespace vg {
namespace algorithms {

using namespace std;

vector<tuple<handle_t, int64_t, bool>> find_closest_with_paths(const PathHandleGraph& graph, handle_t start, size_t offset, size_t max_search_dist) {
    
    
#ifdef debug_algorithms
    cerr << "looking for closest postition on a path starting from " << graph.get_id(start) << (graph.get_is_reverse(start) ? "-" : "+") << ":" << offset << " with max search dist " << max_search_dist << endl;
#endif
    
    // Holds all the handles with paths, their unsigned distances from the start handle, and whether they were found going left.
    vector<tuple<handle_t, int64_t, bool>> to_return;
    
    // Use greater so that we traverse in ascending order of distance
    structures::RankPairingHeap<pair<handle_t, bool>, int64_t, greater<int64_t>> queue;
    
    // add in the initial traversals in both directions from the start position
    // distances are measured to the "incoming" side of the node: the end of the node
    // if searching left and the beginning of the node if searching right
    queue.push_or_reprioritize(make_pair(start, true), offset - graph.get_length(start));
    queue.push_or_reprioritize(make_pair(start, false), -offset);
    
    while (!queue.empty() && to_return.empty()) {
        // get the queue that has the next shortest path
        auto trav = queue.top();
        queue.pop();
        
        // unpack this record
        handle_t here = trav.first.first;
        bool search_left = trav.first.second;
        int64_t dist = trav.second;
        
#ifdef debug_algorithms
        cerr << "traversing " << graph.get_id(here) << (graph.get_is_reverse(here) ? "-" : "+")
        << " in " << (search_left ? "leftward" : "rightward") << " direction at distance " << dist << endl;
#endif
        
        graph.for_each_step_on_handle(here, [&](const step_handle_t& step) {
            // For each path visit that occurs on this node
            
            // Record that fact.
            to_return.emplace_back(here, dist, search_left);
        });
        
        if (!to_return.empty()) {
            // we found the closest, we're done
            break;
        }
        
        int64_t dist_thru = dist + graph.get_length(here);
        
        if (dist_thru <= max_search_dist) {
            // we can cross the node within our budget of search distance, enqueue
            // the next nodes in the search direction
            graph.follow_edges(here, search_left, [&](const handle_t& next) {
                
#ifdef debug_algorithms
                cerr << "\tfollowing edge to " << graph.get_id(next) << (graph.get_is_reverse(next) ? "-" : "+")
                << " at dist " << dist_thru << endl;
#endif
                
                queue.push_or_reprioritize(make_pair(next, search_left), dist_thru);
            });
        }
    }
    
    return std::move(to_return);
}
    
    
}
}
