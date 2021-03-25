#ifndef VG_ALGORITHMS_A_STAR_HPP_INCLUDED
#define VG_ALGORITHMS_A_STAR_HPP_INCLUDED

/**
 * \file a_star.hpp
 *
 * Defines an implementation of the A* heuristic-guided search algorithm.
 */

#include "../handle.hpp"

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
                            int64_t extremal_distance = numeric_limits<int64_t>::max(),
                            bool monotonic_heuristic = true);
    
    // TODO: non-monotonic heuristics in cyclic graphs cannot identify unreachable positions
    
    
    
    
    
    /*
     * Implementation of template functions
     */
    
    template<class DistHeuristic>
    vector<handle_t> a_star(const HandleGraph* graph,
                            const pos_t& pos_1, const pos_t& pos_2,
                            const DistHeuristic& dist_heuristic,
                            bool find_min,
                            int64_t extremal_distance,
                            bool monotonic_heuristic) {
        
#ifdef debug_a_star
        cerr << "doing A* search from " << pos_1 << " to " << pos_2 << " looking for " << (find_min ? "min" : "max") << " distance path with breakout distance " << extremal_distance << endl;
#endif
        
        /*
         * The last node in a search path, pointing back to the previous node for traceback purposes
         */
        struct AStarPathEnd {
            AStarPathEnd(handle_t handle, int64_t predecessor_idx)
                : handle(handle), predecessor_idx(predecessor_idx) { }
            
            handle_t handle;            // current node
            int64_t predecessor_idx;    // the index of the predecessor node in the search history or -1 if none
        };
        
        /*
         * A step in the search, including the node, the traveled distance, and the heuristic estimate of
         * the length of the whole path
         */
        struct AStarSearchRecord {
            
            AStarSearchRecord(const AStarPathEnd& path_end, int64_t distance, int64_t heuristic_estimate)
                : path_end(path_end), distance(distance), heuristic_estimate(heuristic_estimate) { }
            
            AStarPathEnd path_end;      // the end of the current path
            int64_t distance;           // distance up to this point in traversed path
            int64_t heuristic_estimate; // estimate of the total distance to the target
        };
        
        /*
         * A stateful comparator that prioritizes search steps based on their heuristic distance
         */
        struct AStarCompare {
            AStarCompare(bool do_min) : do_min(do_min) { }
            
            // reverse the ordering for a min to match the behavior of a priority queue (which selects
            // the max value)
            inline bool operator()(const AStarSearchRecord& a, const AStarSearchRecord& b) const {
                return ((do_min && a.heuristic_estimate > b.heuristic_estimate)
                        || (!do_min && a.heuristic_estimate < b.heuristic_estimate));
            }
            
            bool do_min;
        };
        
        
        // TODO: this function won't handle cycles as written
        // it needs to be able to move past the end node and circle back around to the start
        // node without breaking the predecessor (i.e. have multiple predecessors)
        
        // TODO: handle a sentinel for being unreachable
        
        // handle the same node reachable case as an edge case
        if (id(pos_1) == id(pos_2) && is_rev(pos_1) == is_rev(pos_2) && offset(pos_1) <= offset(pos_2)
            && (find_min ? offset(pos_2) - offset(pos_1) <= extremal_distance : offset(pos_2) - offset(pos_1) >= extremal_distance)) {
            return vector<handle_t>(1, graph->get_handle(id(pos_1), is_rev(pos_1)));
            
        }
        
        // the paths we've walked so far
        vector<AStarPathEnd> path_search_history;
        
        // the nodes we've decided we don't need to revisit
        unordered_set<handle_t> closed_nodes;
        
        handle_t start = graph->get_handle(id(pos_1), is_rev(pos_1));
        handle_t end = graph->get_handle(id(pos_2), is_rev(pos_2));
        
        // init the priority queue, ordered by either min or max
        AStarCompare compare(find_min);
        priority_queue<AStarSearchRecord, vector<AStarSearchRecord>, AStarCompare> queue(compare);
        
        // set negative distance so the search starts from the beginning of the node
        queue.emplace(AStarPathEnd(start, -1),
                      -offset(pos_1),
                      dist_heuristic(make_pos_t(id(pos_1), is_rev(pos_1), 0), pos_2));
        
        while (!queue.empty()) {
            auto search_record = queue.top();
            queue.pop();
            
            // have we marked this node as one we don't need to traverse again?
            if (closed_nodes.count(search_record.path_end.handle)) {
                continue;
            }
            
            // create a record of this search step in our history
            int64_t current_idx = path_search_history.size();
            path_search_history.emplace_back(search_record.path_end);
            
#ifdef debug_a_star
            cerr << "at node " << graph->get_id(search_record.path_end.handle) << (graph->get_is_reverse(search_record.path_end.handle) ? "-" : "+") << " with heuristic distance " << search_record.heuristic_estimate << " and prefix distance " << search_record.distance << " from predecessor " << graph->get_id(path_search_history[search_record.path_end.predecessor_idx].handle) << (graph->get_is_reverse(path_search_history[search_record.path_end.predecessor_idx].handle) ? "-" : "+") << endl;
#endif
            
            if ((find_min && search_record.distance + int64_t(offset(pos_2)) > extremal_distance)
                || (!find_min && search_record.distance + int64_t(offset(pos_2)) < extremal_distance)) {
                // we've crossed over the distance beyond which we're not interested
                
#ifdef debug_a_star
                cerr << "\twe've crossed the extremal distance of " << (find_min ? extremal_distance - offset(pos_2) : extremal_distance + offset(pos_2)) << endl;
#endif
                break;
            }
            
            // we never need to return here if we're using a monotonic heuristic
            if (monotonic_heuristic) {
                closed_nodes.insert(search_record.path_end.handle);
            }
            
            if (path_search_history.size() > 1 && // in first step, we only got to this loop if offsets are unreachable
                search_record.path_end.handle == end) {
                // we've found the end, reconstruct the path and return it
                
#ifdef debug_a_star
                cerr << "\tfound target" << endl;
                cerr << "path search history:" << endl;
                for (size_t i = 0; i < path_search_history.size(); i++) {
                    cerr << "\t" << i << ": " << graph->get_id(path_search_history[i].handle) << (graph->get_is_reverse(path_search_history[i].handle) ? "-" : "+") << " -> " << path_search_history[i].predecessor_idx << endl;
                }
#endif
                
                vector<handle_t> path;
                // walk the path backwards through the searh history
                int64_t idx = path_search_history.size() - 1;
                while (idx >= 0) {
                    auto& path_step = path_search_history[idx];
                    path.emplace_back(path_step.handle);
                    idx = path_step.predecessor_idx;
                }
                // put the path in forward order
                reverse(path.begin(), path.end());
                return path;
            }
            
            int64_t distance = graph->get_length(search_record.path_end.handle) + search_record.distance;
            
#ifdef debug_a_star
            cerr << "\tdistance through node is " << distance << endl;
#endif
            
            auto enqueue_next = [&](const handle_t& next) {
                
                int64_t remaining = dist_heuristic(make_pos_t(graph->get_id(next), graph->get_is_reverse(next), 0),
                                                   pos_2);
#ifdef debug_a_star
                cerr << "\ttraversing to " << graph->get_id(next) << (graph->get_is_reverse(next) ? "-" : "+") << " with guess for remaining distance at " << remaining << endl;
#endif
                
                queue.emplace(AStarPathEnd(next, current_idx), distance, distance + remaining);
            };
            
            // enqueue walking forward
            graph->follow_edges(search_record.path_end.handle, false, enqueue_next);
        }
        
        // we made it through the loop without finding a path
        return vector<handle_t>();
    }
}
}

#endif
