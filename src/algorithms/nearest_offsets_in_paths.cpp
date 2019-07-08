/**
 * \file nearest_offsets_in_paths.cpp
 *
 * Contains implementation of nearest_offsets_in_paths function
 */

#include "nearest_offsets_in_paths.hpp"

//#define debug

namespace vg {
namespace algorithms {

using namespace std;

unordered_map<path_handle_t, vector<pair<size_t, bool>>> nearest_offsets_in_paths(const PathPositionHandleGraph* graph,
                                                                                  const pos_t& pos,
                                                                                  int64_t max_search) {
    
    // init the return value
    unordered_map<path_handle_t, vector<pair<size_t, bool>>> return_val;
    
    // use greater so that we traverse in ascending order of distance
    structures::RankPairingHeap<pair<handle_t, bool>, int64_t, greater<int64_t>> queue;
    
    // add in the initial traversals in both directions from the start position
    // distances are measured to the left side of the node
    handle_t start = graph->get_handle(id(pos), is_rev(pos));
    queue.push_or_reprioritize(make_pair(start, false), -offset(pos));
    queue.push_or_reprioritize(make_pair(graph->flip(start), true), offset(pos) - graph->get_length(start));
    
    while (!queue.empty()) {
        // get the queue that has the next shortest path
        auto trav = queue.top();
        queue.pop();
        
        // unpack this record
        handle_t here = trav.first.first;
        bool search_left = trav.first.second;
        int64_t dist = trav.second;
        
#ifdef debug_algorithms
        cerr << "traversing " << graph->get_id(here) << (graph->get_is_reverse(here) ? "-" : "+")
        << " in " << (search_left ? "leftward" : "rightward") << " direction at distance " << dist << endl;
#endif
        
        for (const step_handle_t& step : graph->steps_of_handle(here)) {
            // For each path visit that occurs on this node
#ifdef debug
            cerr << "handle is on step at path offset " << graph->get_position_of_step(step) << endl;
#endif
            
            // flip the handle back to the orientation it started in
            handle_t oriented = search_left ? graph->flip(here) : here;
            
            // the orientation of the position relative to the forward strand of the path
            bool rev_on_path = (oriented != graph->get_handle_of_step(step));
            
            // the offset of this step on the forward strand
            int64_t path_offset = graph->get_position_of_step(step);
            
            if (rev_on_path == search_left) {
                path_offset -= dist;
            }
            else {
                path_offset += graph->get_length(oriented) + dist;
            }
            
#ifdef debug
            cerr << "after adding search distance and node offset, " << path_offset << " on strand " << rev_on_path << endl;
#endif
            
            path_handle_t path_handle = graph->get_path_handle_of_step(step);
            
            // handle possible under/overflow from the search distance
            path_offset = max<int64_t>(min<int64_t>(path_offset, graph->get_path_length(path_handle)), 0);
            
            // add in the search distance and add the result to the output
            return_val[path_handle].emplace_back(path_offset, rev_on_path);
        }
        
        if (!return_val.empty()) {
            // we found the closest, we're done
            break;
        }
        
        int64_t dist_thru = dist + graph->get_length(here);
        
        if (dist_thru <= max_search) {
            
            // we can cross the node within our budget of search distance, enqueue
            // the next nodes in the search direction
            graph->follow_edges(here, false, [&](const handle_t& next) {
                
#ifdef debug_algorithms
                cerr << "\tfollowing edge to " << graph->get_id(next) << (graph->get_is_reverse(next) ? "-" : "+")
                << " at dist " << dist_thru << endl;
#endif
                
                queue.push_or_reprioritize(make_pair(next, search_left), dist_thru);
            });
        }
    }
    
    return return_val;
}

map<string, vector<pair<size_t, bool>>> offsets_in_paths(const PathPositionHandleGraph* graph, const pos_t& pos) {
    auto offsets = nearest_offsets_in_paths(graph, pos, -1);
    map<string, vector<pair<size_t, bool>>> named_offsets;
    for (pair<const path_handle_t, vector<pair<size_t, bool>>>& offset : offsets) {
        named_offsets[graph->get_path_name(offset.first)] = move(offset.second);
    }
    return named_offsets;
}
    
}
}
