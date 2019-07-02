/**
 * \file jump_along_path.cpp
 *
 * Contains implementation of jump_along_closest_path algorithm
 */

#include "jump_along_path.hpp"

//#define debug_algorithms

namespace vg {
namespace algorithms {

using namespace std;

    vector<pos_t> jump_along_closest_path(const PathPositionHandleGraph* graph,
                                          const pos_t& pos, int64_t jump_dist,
                                          size_t max_search_dist) {
#ifdef debug_algorithms
        cerr << "jumping " << jump_dist << " from " << pos << " with a max search dist of " << max_search_dist << endl;
#endif
        
        // intialize the return value
        vector<pos_t> to_return;
        
        function<bool(handle_t,int64_t,bool)> look_for_jumpable_paths = [&](handle_t handle, int64_t search_dist, bool search_left) {
            
            bool found_jumpable_path = false;
            
#ifdef debug_algorithms
            cerr << "checking for jumpable paths for " << graph->get_id(handle) << (graph->get_is_reverse(handle) ? "-" : "+") << " at search dist " << search_dist << " from searching " << (search_left ? "leftwards" : "rightwards") << endl;
#endif
            
            for (const step_handle_t& step : graph->steps_of_handle(handle)) {
                bool path_rev_strand = (graph->get_is_reverse(graph->get_handle_of_step(step))
                                        != graph->get_is_reverse(handle));
                
                int64_t dist_to_jump;
                if (path_rev_strand && search_left) {
                    dist_to_jump = -jump_dist - search_dist;
                }
                else if (path_rev_strand) {
                    dist_to_jump = -jump_dist + search_dist + graph->get_length(handle);
                }
                else if (search_left) {
                    dist_to_jump = jump_dist + search_dist + graph->get_length(handle);
                }
                else {
                    dist_to_jump = jump_dist - search_dist;
                }
                
                int64_t target_path_pos = graph->get_position_of_step(step) + dist_to_jump;
                if (target_path_pos >= 0
                    && target_path_pos < graph->get_path_length(graph->get_path_handle_of_step(step))) {
                    
                    step_handle_t jump_step = graph->get_step_at_position(graph->get_path_handle_of_step(step),
                                                                          target_path_pos);
                    
                    
                    size_t step_offset = graph->get_position_of_step(jump_step);
                    handle_t jump_handle = graph->get_handle_of_step(jump_step);
                    
                    to_return.emplace_back(graph->get_id(jump_handle),
                                           graph->get_is_reverse(jump_handle) != path_rev_strand,
                                           path_rev_strand
                                              ? graph->get_position_of_step(jump_step) + graph->get_length(jump_handle) - target_path_pos
                                              : target_path_pos - graph->get_position_of_step(jump_step));
                    
                    found_jumpable_path = true;
#ifdef debug_algorithms
                    cerr << "found jump position " << to_return.back() << " at forward path offset " << target_path_pos << endl;
#endif
                }
            }
            
            return found_jumpable_path;
        };
        
        // records of (handle, searching left) prioritized by distance
        // all distances are measured to the left end of the node
        structures::RankPairingHeap<pair<handle_t, bool>, size_t> queue;
        
        handle_t handle = graph->get_handle(id(pos), is_rev(pos));
        
        // add in the initial traversals in both directions from the start position
        queue.push_or_reprioritize(make_pair(handle, true), offset(pos));
        queue.push_or_reprioritize(make_pair(handle, false), graph->get_length(handle) - offset(pos));
        
        bool found_jumpable_path = false;
        while (!queue.empty() && !found_jumpable_path) {
            // get the queue that has the next shortest path
            
            auto here = queue.top();
            queue.pop();
            
#ifdef debug_algorithms
            cerr << "traversing " << graph->get_id(here.first.first) << (graph->get_is_reverse(here.first.first) ? "-" : "+") << " in " << (here.first.second ? "leftward" : "rightward") << " direction at distance " << here.second << endl;
#endif
            
            graph->follow_edges(here.first.first, here.first.second, [&](const handle_t& next) {
#ifdef debug_algorithms
                cerr << "\tfollowing edge to " << graph->get_id(next) << (graph->get_is_reverse(next) ? "-" : "+") << " at dist " << here.second << endl;
#endif
                
                found_jumpable_path = look_for_jumpable_paths(next, here.second, here.first.second);
                
                int64_t dist_thru = here.second + graph->get_length(next);
                if (dist_thru <= (int64_t) max_search_dist) {
                    queue.push_or_reprioritize(make_pair(next, here.first.second), dist_thru);
                }
                
                return !found_jumpable_path;
            });
        }
        
        return to_return;
        
    }
}
}
