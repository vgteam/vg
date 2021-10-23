/** \file
 * Implements the reference path distance function
 */

#include "ref_path_distance.hpp"

//#define debug_ref_path_distance

namespace vg {
namespace algorithms {

using namespace std;

int64_t ref_path_distance(const PathPositionHandleGraph* graph, const pos_t& pos_1, const pos_t& pos_2,
                          const unordered_set<path_handle_t>& ref_paths, int64_t max_search_dist) {
    
    // to record the nearest position on the strands of paths for each of the
    // two positions
    unordered_map<pair<path_handle_t, bool>, int64_t> nearby_paths_1, nearby_paths_2;
    
    // dijkstra priority queue, records of (node, from pos_1)
    structures::RankPairingHeap<tuple<handle_t, bool>, int64_t, greater<int64_t>> queue;
    
    // intialize in both directions from both positions
    handle_t handle_1 = graph->get_handle(id(pos_1), is_rev(pos_1));
    handle_t handle_2 = graph->get_handle(id(pos_2), is_rev(pos_2));
    
    // we only walk forward (toward each other) so that we know the positions can reach each other
    queue.push_or_reprioritize(make_tuple(handle_1, true), -offset(pos_1));
    queue.push_or_reprioritize(make_tuple(handle_2, false),
                               offset(pos_2) - graph->get_length(handle_2));
    
    vector<pair<path_handle_t, bool>> shared_refs;
    while (!queue.empty() && shared_refs.empty()) {
        
        auto top = queue.top();
        handle_t handle;
        bool from_pos_1;
        tie(handle, from_pos_1) = top.first;
        queue.pop();
        
#ifdef debug_ref_path_distance
        cerr << "[ref_path_distance] dequeue " << graph->get_id(handle) << (graph->get_is_reverse(handle) ? "-" : "+") << ", left? " << !from_pos_1 << " from pos " << (from_pos_1 ? 1 : 2) << endl;
#endif
        
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
        
        graph->for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            pair<path_handle_t, bool> oriented_path(graph->get_path_handle_of_step(step),
                                                    graph->get_handle_of_step(step) != handle);
#ifdef debug_ref_path_distance
            cerr << "[ref_path_distance] on path " << graph->get_path_name(oriented_path.first) << ", on rev? " << oriented_path.second << ", path offset " << graph->get_position_of_step(step) << endl;
#endif
            if (!nearby_paths->count(oriented_path)) {
                
                int64_t path_offset;
                if (oriented_path.second && !from_pos_1) {
                    // traverse left, reverse strand of path
                    path_offset = (graph->get_path_length(graph->get_path_handle_of_step(step))
                                   - graph->get_position_of_step(step));
                }
                else if (oriented_path.second) {
                    // traverse right, reverse strand of path
                    path_offset = (graph->get_path_length(graph->get_path_handle_of_step(step))
                                   - graph->get_position_of_step(step) - graph->get_length(handle));
                }
                else if (!from_pos_1) {
                    // traverse left, forward strand of path
                    path_offset = graph->get_position_of_step(step) + graph->get_length(handle);
                }
                else {
                    // traverse right, forward strand of path
                    path_offset = graph->get_position_of_step(step);
                }
                
                // we include the offset on the node if we actually started on the path
                if (from_pos_1 && handle == handle_1) {
                    path_offset += offset(pos_1);
                }
                else if (!from_pos_1 && handle == handle_2) {
                    path_offset += offset(pos_2);
                }
                
#ifdef debug_ref_path_distance
                cerr << "[ref_path_distance] first encounter, recording path offset of " << path_offset << endl;
#endif
                                   
                (*nearby_paths)[oriented_path] = path_offset;
                if (ref_paths.count(oriented_path.first) && other_nearby_paths->count(oriented_path)) {
                    shared_refs.emplace_back(oriented_path);
                }
            }
        });
                    
        
        if (shared_refs.empty()) {
            // only queue up the next if we're still within the max distance
            int64_t dist_thru = top.second + graph->get_length(handle);
            if (dist_thru <= max_search_dist) {
                graph->follow_edges(handle, !from_pos_1, [&](const handle_t& next) {
                    queue.push_or_reprioritize(make_tuple(next, from_pos_1), dist_thru);
                });
            }
        }
    }
    
#ifdef debug_ref_path_distance
    cerr << "[ref_path_distance] choosing max absolute distance on shared paths" << endl;
#endif
    
    int64_t approx_ref_dist = numeric_limits<int64_t>::max();
    if (!shared_refs.empty()) {
        // we found a labeled reference, measure distance using that
        for (const auto& ref : shared_refs) {
            int64_t dist = nearby_paths_2[ref] - nearby_paths_1[ref];
            if (approx_ref_dist == numeric_limits<int64_t>::max() || dist > approx_ref_dist) {
                approx_ref_dist = dist;
            }
        }
    }
    else {
        // try among the remaining paths
        for (const auto& path_record_1 : nearby_paths_1) {
            auto it = nearby_paths_2.find(path_record_1.first);
            if (it != nearby_paths_2.end()) {
                int64_t dist = it->second - path_record_1.second;
                
#ifdef debug_ref_path_distance
                cerr << "[ref_path_distance] distance on path " << graph->get_path_name(path_record_1.first.first) << " is " << dist << endl;
#endif
                if (approx_ref_dist == numeric_limits<int64_t>::max() || dist > approx_ref_dist) {
                    approx_ref_dist = dist;
                }
            }
        }
    }
    
#ifdef debug_ref_path_distance
    cerr << "[ref_path_distance] max absolute distance is " << approx_ref_dist << endl;
#endif
    
    return approx_ref_dist;
}
    
}
}
