#ifndef VG_ALGORITHMS_A_STAR_HPP_INCLUDED
#define VG_ALGORITHMS_A_STAR_HPP_INCLUDED

/**
 * \file a_star.hpp
 *
 * Defines an implementation of the A* heuristic-guided search algorithm.
 */

#include "../handle.hpp"
#include "topological_sort.hpp"

#include <unordered_map>
#include <vector>

//#define debug_a_star

namespace vg {
namespace algorithms {

using namespace std;

    /// Implements the A* heuristic-guided search algorithm. Returns the path from pos_1 to
    /// pos_2 that is either minimal or maximal length, according to the parameters. Allows
    /// an extremal distance beyond which the algorithm will cease looking for paths (this
    /// should be a large value when looking for minimal paths and a small value when looking
    /// for maximum paths). If there is no path between the positions, or none within the
    /// extremal length, an empty vector will be returned.
    template<class DistHeuristic>
    vector<handle_t> a_star(const HandleGraph* graph,
                            const pos_t& pos_1, const pos_t& pos_2,
                            const DistHeuristic& dist_heuristic,
                            bool find_min = true,
                            int64_t extremal_distance = numeric_limits<int64_t>::max());
    
    
    
    
    
    
    
    /*
     * Implementation of template functions
     */
    
    template<class DistHeuristic>
    vector<handle_t> a_star(const HandleGraph* graph,
                            const pos_t& pos_1, const pos_t& pos_2,
                            const DistHeuristic& dist_heuristic,
                            bool find_min,
                            int64_t extremal_distance) {
        
        struct AStarPathEnd {
            AStarPathEnd(handle_t handle, int64_t distance, int64_t heuristic_estimate,
                         handle_t predecessor, bool find_min)
            : handle(handle), distance(distance), heuristic_estimate(heuristic_estimate), predecessor(predecessor), find_min(find_min) {
            }
            
            handle_t handle;            // current node
            int64_t distance;           // distance up to this point in traversed path
            int64_t heuristic_estimate; // estimate of the total distance to the target
            handle_t predecessor;       // the node we came from (for backtracking)
            // TODO: there must be a design where I don't need to copy this in every struct
            bool find_min;
            
            // reverse the ordering for a min to match the behavior of a priority queue (which selects
            // the max value)
            inline bool operator<(const AStarPathEnd& other) const {
                return ((find_min && heuristic_estimate > other.heuristic_estimate)
                        || (!find_min && other.heuristic_estimate < heuristic_estimate));
            }
        };
        
        // TODO: this function won't handle cycles as written
        // it needs to be able to move past the end node and circle back around to the start
        // node without breaking the predecessor (i.e. have multiple predecessors)
        
        // TODO: handle the same-node unreachable case
        
        // TODO: handle a sentinel for being unreachable
        
        // when we traverse a node,
        unordered_map<handle_t, handle_t> best_predecessor;
        
        
        handle_t start = graph->get_handle(id(pos_1), is_rev(pos_1));
        handle_t end = graph->get_handle(id(pos_2), is_rev(pos_2));
        
        priority_queue<AStarPathEnd> queue;
        // set negative distance so the search starts from the beginning of the node
        queue.emplace(start,
                      -offset(pos_1),
                      dist_heuristic(make_pos_t(id(pos_1), is_rev(pos_1), 0), pos_2),
                      as_handle(int64_t(0)), // dummy handle, no real predecessor
                      find_min);
        
        while (!queue.empty()) {
            auto path_end = queue.top();
            queue.pop();
            
#ifdef debug_a_star
            cerr << "at node " << graph->get_id(path_end.handle) << (graph->get_is_reverse(path_end.handle) ? "-" : "+") << " with heuristic distance " << path_end.heuristic_estimate << " and prefix distance " << path_end.distance << " from predecessor " << graph->get_id(path_end.predecessor) << (graph->get_is_reverse(path_end.predecessor) ? "-" : "+") << endl;
#endif
            
            if ((find_min && path_end.distance + int64_t(offset(pos_2)) > extremal_distance)
                || (!find_min && path_end.distance + int64_t(offset(pos_2)) < extremal_distance)) {
                // we've crossed over the distance beyond which we're not interested
                
#ifdef debug_a_star
                cerr << "\twe've crossed the extremal distance of " << (find_min ? extremal_distance - offset(pos_2) : extremal_distance + offset(pos_2)) << endl;
#endif
                break;
            }
            
            if (best_predecessor.count(path_end.handle)) {
                continue;
            }
            if (path_end.handle != start) {
                // TODO: prohibits cycles involving start node
                best_predecessor[path_end.handle] = path_end.predecessor;
            }
            
            if (path_end.handle == end) {
                // we've found the end, reconstruct the path and return it
                vector<handle_t> path(1, path_end.handle);
                while (best_predecessor.count(path.back())) {
                    path.emplace_back(best_predecessor[path.back()]);
                }
                // put the path in forward order
                reverse(path.begin(), path.end());
                return path;
            }
            
            int64_t distance = graph->get_length(path_end.handle) + path_end.distance;
            
#ifdef debug_a_star
            cerr << "\tdistance through node is " << distance << endl;
#endif
            
            auto enqueue_next = [&](const handle_t& next) {
                
                int64_t remaining = dist_heuristic(make_pos_t(graph->get_id(next), graph->get_is_reverse(next), 0),
                                                   pos_2);
#ifdef debug_a_star
                cerr << "\ttraversing to " << graph->get_id(next) << (graph->get_is_reverse(next) ? "-" : "+") << " with guess for remaining distance at " << remaining << endl;
#endif
                
                queue.emplace(next, distance, distance + remaining, path_end.handle, find_min);
            };
            
            // enqueue walking forward
            graph->follow_edges(path_end.handle, false, enqueue_next);
        }
        
        // we made it through the loop without finding a path
        return vector<handle_t>();
    }
}
}

#endif
