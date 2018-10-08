#include "count_walks.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    size_t count_walks(const HandleGraph* graph) {
        
        vector<handle_t> sinks;
        unordered_map<handle_t, size_t> count;
        count.reserve(graph->node_size());
        
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
        bool overflowed = false;
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
            
            if (overflowed) {
                return numeric_limits<size_t>::max();
            }
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
