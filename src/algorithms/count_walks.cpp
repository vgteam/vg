#include "count_walks.hpp"

namespace vg {
namespace algorithms {

using namespace std;
    
    

     tuple<vector<handle_t>, unordered_map<handle_t, size_t>, bool> count_walks_through_nodes(const HandleGraph* graph) {
        
        
        tuple<vector<handle_t>, unordered_map<handle_t, size_t>, bool> to_return; 

        vector<handle_t>& sinks = get<0>(to_return);
        unordered_map<handle_t, size_t>& count = get<1>(to_return);
        bool& overflowed = get<2>(to_return);
        
        count.reserve(graph->get_node_count());
        
        // identify sources and sinks
        graph->for_each_handle([&](const handle_t& handle) {
            bool is_source = true, is_sink = true;
            graph->follow_edges(handle, true, [&](const handle_t& prev) {
                is_source = false;
                return false;
            });
            graph->follow_edges(handle, false, [&](const handle_t& next) {
                is_sink = false;
                return false;
            });
            
            // base case for dynamic programming
            if (is_source) {
                count[handle] = 1;
            }
            if (is_sink) {
                sinks.emplace_back(handle);
            }
        });
        
        // count walks by dynamic programming
        for (const handle_t& handle : lazier_topological_order(graph)) {
            size_t count_here = count[handle];
            graph->follow_edges(handle, false, [&](const handle_t& next) {
                size_t& count_next = count[next];
                if (numeric_limits<size_t>::max() - count_here < count_next) {
                    overflowed = true;
                }
                else {
                    count_next += count_here;
                }
            });
        }
        return to_return;
    }
    size_t count_walks(const HandleGraph* graph){

        tuple<vector<handle_t>, unordered_map<handle_t, size_t>, bool>  to_receive = count_walks_through_nodes(graph);

        vector<handle_t>& sinks = get<0>(to_receive);
        unordered_map<handle_t, size_t>& count = get<1>(to_receive);
        bool& overflowed = get<2>(to_receive); 

        if (overflowed) {
                return numeric_limits<size_t>::max();
            }
 
        // total up the walks at the sinks
        size_t total_count = 0;
        for (handle_t& sink : sinks) {
            total_count += count[sink];
        }     
        return total_count;
    }
}
}
