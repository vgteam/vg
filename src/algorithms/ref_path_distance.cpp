/** \file
 * Implementsthe reference path distance function
 */

#include "ref_path_distance.hpp"

namespace vg {
namespace algorithms {

using namespace std;

int64_t ref_path_distance(const PathPositionHandleGraph* graph, const pos_t& pos_1, const pos_t& pos_2,
                          int64_t min_search_dist, int64_t max_search_dist) {
    
    // to record the nearest position on the strands of paths for each of the
    // two positions
    unordered_map<pair<path_handle_t, bool>, int64_t> nearby_paths_1, nearby_paths_2;
    
    // dijkstra priority queue, records of (node, going left, from pos_1)
    structures::RankPairingHeap<tuple<handle_t, bool, bool>, int64_t, greater<int64_t>> queue;
    
    // intialize in both directions from both positions
    handle_t handle_1 = graph->get_handle(id(pos_1), is_rev(pos_1));
    handle_t handle_2 = graph->get_handle(id(pos_2), is_rev(pos_2));
    queue.push_or_reprioritize(make_tuple(handle_1, true, true),
                               offset(pos_1) - graph->get_length(handle_1));
    queue.push_or_reprioritize(make_tuple(handle_1, false, true), -offset(pos_1));
    queue.push_or_reprioritize(make_tuple(handle_2, true, true),
                               offset(pos_2) - graph->get_length(handle_2));
    queue.push_or_reprioritize(make_tuple(handle_2, false, true), -offset(pos_2));
    
    bool found_shared = false;
    while (!queue.empty()) {
        
        auto top = queue.top();
        queue.pop();
        
        handle_t handle;
        bool go_left, from_pos_1;
        tie(handle, go_left, from_pos_1) = top.first;
        int64_t dist = top.second;
        
        // after the min search distance, we can quit as soon as we find a
        // shared path
        if (found_shared && dist > min_search_dist) {
            break;
        }
        
        decltype(nearby_paths_1)* nearby_paths;
        decltype(nearby_paths_1)* other_nearby_paths;
        if (from_pos_1) {
            nearby_paths = &nearby_paths_1;
            other_nearby_paths = &nearby_paths_2;
        }
        else {
            nearby_paths = &nearby_paths_2;
            other_nearby_paths = &nearby_paths_1;
        }
        
        int64_t dist_thru = dist + graph->get_length(handle);
        
        graph->for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            pair<path_handle_t, bool> oriented_path(graph->get_path_handle_of_step(step),
                                                    graph->get_handle_of_step(step) != handle);
            
            if (!nearby_paths->count(oriented_path)) {
                
                int64_t path_offset;
                if (oriented_path.second && go_left) {
                    // traverse left, reverse strand of path
                    path_offset = (graph->get_path_length(graph->get_path_handle_of_step(step))
                                   - graph->get_position_of_step(step) + dist);
                }
                else if (oriented_path.second) {
                    // traverse right, reverse strand of path
                    path_offset = (graph->get_path_length(graph->get_path_handle_of_step(step))
                                   - graph->get_position_of_step(step) - dist_thru);
                }
                else if (go_left) {
                    // traverse left, forward strand of path
                    path_offset = (graph->get_position_of_step(step) + dist_thru);
                }
                else {
                    // traverse right, forward strand of path
                    path_offset = graph->get_position_of_step(step) - dist;
                }
                                   
                (*nearby_paths)[oriented_path] = path_offset;
                found_shared = found_shared || other_nearby_paths->count(oriented_path);
            }
        });
        
        // only queue up the next if we're still within the max distance
        if (dist_thru <= max_search_dist) {
            graph->follow_edges(handle, go_left, [&](const handle_t& next) {
                queue.push_or_reprioritize(make_tuple(next, go_left, from_pos_1), dist_thru);
            });
        }
    }
    
    // check all shared paths that we found
    int64_t approx_ref_dist = numeric_limits<int64_t>::max();
    for (const auto& path_record_1 : nearby_paths_1) {
        auto it = nearby_paths_2.find(path_record_1.first);
        if (it != nearby_paths_2.end()) {
            int64_t dist = it->second - path_record_1.second;
            // we'll assume that the reference distance is the longest distance along
            // a shared path (because others may correspond to spliced transcripts)
            if (approx_ref_dist == numeric_limits<int64_t>::max() ||
                abs(dist) > abs(approx_ref_dist)) {
                approx_ref_dist = dist;
            }
        }
    }
    
    return -1;
}
    
}
}
