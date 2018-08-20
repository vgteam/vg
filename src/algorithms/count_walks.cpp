#include "count_walks.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    size_t count_walks(const HandleGraph* graph) {
        
        vector<handle_t> sinks;
        unordered_map<handle_t, size_t> count;
        count.reserve(graph->node_size());
        
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
            
            if (is_source) {
                count[handle] = 1;
            }
            if (is_sink) {
                sinks.emplace_back(handle);
            }
        });
        
        for (const handle_t& handle : lazier_topological_order(graph)) {
            size_t count_here = count[handle];
            graph->follow_edges(handle, false, [&](const handle_t& next) {
                count[next] += count_here;
            });
        }
        
        size_t total_count = 0;
        for (handle_t& sink : sinks) {
            total_count += count[sink];
        }
        
        return total_count;
    }
}
}
