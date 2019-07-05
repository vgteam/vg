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

vector<tuple<handle_t, size_t, bool>> find_closest_with_paths(const PathHandleGraph& graph, handle_t start, size_t max_search_dist,
    size_t right_extra_dist, size_t left_extra_dist) {
    
    // Holds all the handles with paths, their unsigned distances from the start handle, and whether they were found going left.
    vector<tuple<handle_t, size_t, bool>> to_return;
    
    // a local struct for traversals that are ordered by search distance and keep
    // track of whether we're searching to the left or right
    struct Traversal {
        Traversal() {}
        Traversal(int64_t dist, handle_t handle, bool search_left) :
        dist(dist), handle(handle), search_left(search_left) {}
        handle_t handle;
        int64_t dist;
        bool search_left;
        inline bool operator<(const Traversal& other) const {
            return dist > other.dist; // opposite order so priority queue selects minimum
        }
    };
    
    // Check if the node for the given handle occurs on any paths.
    // Remember we reached it at the given distance, and whether we reached it searching left or not.
    // Record the handle and its signed distance if it has paths.
    // Return true if any paths touching the handle were found.
    function<bool(handle_t,int64_t,bool)> look_for_paths = [&](handle_t handle, int64_t search_dist, bool search_left) {
        
        bool found_path = false;
        
        int64_t trav_id = graph.get_id(handle);
        bool trav_is_rev = graph.get_is_reverse(handle);
        
#ifdef debug_algorithms
        cerr << "checking for paths for " << trav_id << (trav_is_rev ? "-" : "+") << " at search dist "
            << search_dist << " from searching " << (search_left ? "leftwards" : "rightwards") << endl;
#endif
        
        // Track if we find any paths on this node
        bool found_anything = false;
        
        graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            // For each path visit that occurs on this node
            
            // Record that fact.
            to_return.emplace_back(handle, search_dist, search_left);
            found_anything = true;
            
            // Don't look at any more
            return false;
        });
        
        return found_anything;
    };
    
    // TODO: This is *NOT* a Dijkstra traversal! We should maybe use a UpdateablePriorityQueue.
    priority_queue<Traversal> queue;
    unordered_set<handle_t> traversed;
    
    handle_t handle = start; 
    
    // add in the initial traversals in both directions from the start position
    queue.emplace(left_extra_dist, handle, true);
    queue.emplace(right_extra_dist, handle, false);
    
    // Start by checking the start node
    traversed.insert(handle);
    bool found_path = look_for_paths(handle, 0, false);
    
    while (!queue.empty() && !found_path) {
        // get the queue that has the next shortest path
        
        Traversal trav = queue.top();
        queue.pop();
        
#ifdef debug_algorithms
        cerr << "traversing " << graph.get_id(trav.handle) << (graph.get_is_reverse(trav.handle) ? "-" : "+")
            << " in " << (trav.search_left ? "leftward" : "rightward") << " direction at distance " << trav.dist << endl;
#endif
        
        function<bool(const handle_t& next)> check_next = [&](const handle_t& next) {
#ifdef debug_algorithms
            cerr << "\tfollowing edge to " << graph.get_id(next) << (graph.get_is_reverse(next) ? "-" : "+")
                << " at dist " << trav.dist << endl;
#endif
            
            if (!traversed.count(next)) {
                found_path = look_for_paths(next, trav.dist, trav.search_left);
                
                int64_t dist_thru = trav.dist + graph.get_length(next);
                
                if (dist_thru <= (int64_t) max_search_dist) {
                    queue.emplace(dist_thru, next, trav.search_left);
                }
                traversed.emplace(next);
            }
            return !found_path;
        };
        
        graph.follow_edges(trav.handle, trav.search_left, check_next);
    }
    
    return std::move(to_return);
}
    
    
}
}
