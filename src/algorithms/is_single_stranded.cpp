#include "is_directed_acyclic.hpp"

#include <unordered_map>

namespace vg {
namespace algorithms {

using namespace std;

bool is_single_stranded(const HandleGraph* graph) {
    
    bool single_stranded = true;
    
    function<bool(const handle_t&)> check_edges = [&](const handle_t& handle) {
        
        function<bool(const handle_t&)> check_edge = [&](const handle_t& next) {
            single_stranded = (graph->get_is_reverse(handle) == graph->get_is_reverse(next));
            return single_stranded;
        };
        
        graph->follow_edges(handle, false, check_edge);
        if (single_stranded) {
            graph->follow_edges(handle, true, check_edge);
        }
        return single_stranded;
    };
    
    graph->for_each_handle(check_edges);
    
    return single_stranded;
}

}
}
