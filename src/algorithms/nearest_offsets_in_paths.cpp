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
    
    unordered_map<path_handle_t, vector<pair<size_t, bool>>> return_val;
    
    pair<pos_t, int64_t> pz = algorithms::next_path_position(*graph, pos, max_search);
    
#ifdef debug
    cerr << "got next path position " << pz.first << " at search distance " << pz.second << endl;
#endif
    
    auto& path_pos = pz.first;
    auto& diff = pz.second;
    if (id(path_pos)) {
        handle_t handle_on_path = graph->get_handle(id(path_pos), is_rev(path_pos));
        for (const step_handle_t& step : graph->steps_of_handle(handle_on_path)) {
            
#ifdef debug
            cerr << "position is on step at path offset " << graph->get_position_of_step(step) << endl;
#endif
            
            path_handle_t path_handle = graph->get_path_handle_of_step(step);
            
            // the offset of this step on the forward strand of the path plus the search distance
            int64_t path_offset = graph->get_position_of_step(step) + diff;
                        
            // handle the offset on the node in a path-strand appropriate manner
            bool rev_on_path = (handle_on_path != graph->get_handle_of_step(step));
            if (rev_on_path) {
                path_offset += graph->get_length(handle_on_path) - offset(path_pos);
            }
            else {
                path_offset += offset(path_pos);
            }
            
#ifdef debug
            cerr << "after adding search distance and node offset, " << path_offset << " on strand " << rev_on_path << endl;
#endif
            
            // handle possible under/overflow from the search distance
            path_offset = max<int64_t>(min<int64_t>(path_offset, graph->get_path_length(path_handle)), 0);
            
            // add in the search distance and add the result to the output
            return_val[path_handle].emplace_back(path_offset, rev_on_path);
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
