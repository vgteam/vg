/**
 * \file nearest_offsets_in_paths.cpp
 *
 * Contains implementation of nearest_offsets_in_paths function
 */

#include "nearest_offsets_in_paths.hpp"
#include "../crash.hpp"

// #define debug
// #define debug_algorithms

namespace vg {
namespace algorithms {

using namespace std;

path_offset_collection_t nearest_offsets_in_paths(const PathPositionHandleGraph* graph,  const pos_t& pos, int64_t max_search,
                                                  const std::unordered_set<PathSense>& desired_senses,
                                                  const std::function<bool(const path_handle_t&)>* path_filter,
                                                  bool subtract_traversal_dist, pair<size_t, bool>* traversal_dist) {
    
    // init the return value
    // This is a map from path handle, to vector of offset and orientation pairs
    path_offset_collection_t return_val;
    
    // use greater so that we traverse in ascending order of distance; the direction and
    // handle in the priority make every priority distinct, so the pop order doesn't
    // depend on the heap's pointer-ordered internals
    typedef tuple<int64_t, bool, uint64_t> priority_t;
    structures::RankPairingHeap<pair<handle_t, bool>, priority_t, greater<priority_t>> queue;
    
    // add in the initial traversals in both directions from the start position
    // distances are measured to the left side of the node
    handle_t start = graph->get_handle(id(pos), is_rev(pos));
    handle_t start_flip = graph->flip(start);
    queue.push_or_reprioritize(make_pair(start, false), priority_t(-offset(pos), false, handlegraph::as_integer(start)));
    queue.push_or_reprioritize(make_pair(start_flip, true), priority_t(offset(pos) - graph->get_length(start), true, handlegraph::as_integer(start_flip)));
    
    while (!queue.empty()) {
        // get the queue that has the next shortest path
        auto trav = queue.top();
        queue.pop();
        
        // unpack this record
        handle_t here = trav.first.first;
        bool search_left = trav.first.second;
        int64_t dist = get<0>(trav.second);
        
#ifdef debug_algorithms
        cerr << "traversing " << graph->get_id(here) << (graph->get_is_reverse(here) ? "-" : "+")
        << " in " << (search_left ? "leftward" : "rightward") << " direction at distance " << dist << endl;
#endif
       
        graph->for_each_step_of_sense(here, desired_senses, [&](const step_handle_t& step) {
            // For each path visit that occurs on this node
#ifdef debug
            cerr << "handle is on step at path offset " << graph->get_position_of_step(step) << endl;
#endif

            path_handle_t path_handle = graph->get_path_handle_of_step(step);

            if (path_filter && !(*path_filter)(path_handle)) {
                // We are to ignore this path
#ifdef debug
                cerr << "handle is on ignored path " << graph->get_path_name(path_handle) << endl;
#endif
                return;
            }

            // flip the handle back to the orientation it started in
            handle_t oriented = search_left ? graph->flip(here) : here;
            
            // the orientation of the position relative to the forward strand of the path
            bool rev_on_path = (oriented != graph->get_handle_of_step(step));
            
            // the offset of this step on the forward strand
            int64_t path_offset = graph->get_position_of_step(step);
            
            if (rev_on_path != search_left) {
                path_offset += graph->get_length(oriented);
            }

            if (subtract_traversal_dist) {
                path_offset += (rev_on_path != search_left ? dist : -dist);
            }
            
#ifdef debug
            cerr << "after adding search distance and node offset, " << path_offset << " on strand " << rev_on_path << endl;
#endif
            
            // handle possible under/overflow from the search distance
            path_offset = max<int64_t>(min<int64_t>(path_offset, graph->get_path_length(path_handle)), 0);
            
            // add in the search distance and add the result to the output
            return_val[path_handle].emplace_back(path_offset, rev_on_path);
        });
        
        if (!return_val.empty()) {
            // we found the closest, we're done
            if (traversal_dist) {
                traversal_dist->first = max<int64_t>(dist, 0);
                traversal_dist->second = search_left;
            }
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
                
                queue.push_or_reprioritize(make_pair(next, search_left), priority_t(dist_thru, search_left, handlegraph::as_integer(next)));
            });
        }
    }
    
    return return_val;
}

map<string, vector<pair<size_t, bool>>> offsets_in_paths(const PathPositionHandleGraph* graph, const pos_t& pos) {
    auto offsets = nearest_offsets_in_paths(graph, pos, -1);
    map<string, vector<pair<size_t, bool>>> named_offsets;
    for (pair<const path_handle_t, vector<pair<size_t, bool>>>& offset : offsets) {
        named_offsets[graph->get_path_name(offset.first)] = std::move(offset.second);
    }
    return named_offsets;
}

path_offset_collection_t simple_offsets_in_paths(const PathPositionHandleGraph* graph, pos_t pos) {
    path_offset_collection_t positions;
    handle_t handle = graph->get_handle(id(pos), is_rev(pos));
    size_t handle_length = graph->get_length(handle);
    for (const step_handle_t& step : graph->steps_of_handle(handle)) {
        // the orientation of the position relative to the forward strand of the path
        bool rev_path = graph->get_is_reverse(graph->get_handle_of_step(step));
        // the offset of this step on the forward strand
        int64_t path_offset = graph->get_position_of_step(step);
        auto& pos_in_path = positions[graph->get_path_handle_of_step(step)];
        // Make sure to interpret the pos_t offset on the correct strand.
        size_t node_forward_strand_offset = is_rev(pos) ? (handle_length - offset(pos) - 1) : offset(pos);
        // Normalize to a forward strand offset.
        size_t off = path_offset + (rev_path ?
                                    (handle_length - node_forward_strand_offset - 1) :
                                    node_forward_strand_offset);
        pos_in_path.push_back(make_pair(off, rev_path));
    }
    return positions;
}
    
}
}
