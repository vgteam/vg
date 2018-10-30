#include "weakly_connected_components.hpp"

namespace vg {
namespace algorithms {

using namespace std;

vector<unordered_set<id_t>> weakly_connected_components(const HandleGraph* graph) {
    vector<unordered_set<id_t>> to_return;
    
    // This only holds locally forward handles
    unordered_set<handle_t> traversed;
    
    graph->for_each_handle([&](const handle_t& handle) {
        
        // Only think about it in the forward orientation
        auto forward = graph->forward(handle);
        
        if (traversed.count(forward)) {
            // Already have this node, so don't start a search from it.
            return;
        }
        
        // The stack only holds locally forward handles
        vector<handle_t> stack{forward};
        to_return.emplace_back();
        while (!stack.empty()) {
            handle_t here = stack.back();
            stack.pop_back();
            
            traversed.insert(here);
            to_return.back().insert(graph->get_id(here));
            
            // We have a function to handle all connected handles
            auto handle_other = [&](const handle_t& other) {
                // Again, make it forward
                auto other_forward = graph->forward(other);
                
                if (!traversed.count(other_forward)) {
                    stack.push_back(other_forward);
                }
            };
            
            // Look at edges in both directions
            graph->follow_edges(here, false, handle_other);
            graph->follow_edges(here, true, handle_other);
            
        }
    });
    return to_return;
}

vector<pair<unordered_set<id_t>, vector<handle_t>>> weakly_connected_components_with_tips(const HandleGraph* graph) {
    // TODO: deduplicate with above
    
    vector<pair<unordered_set<id_t>, vector<handle_t>>> to_return;
    
    // This only holds locally forward handles
    unordered_set<handle_t> traversed;
    
    graph->for_each_handle([&](const handle_t& handle) {
        
        // Only think about it in the forward orientation
        auto forward = graph->forward(handle);
        
        if (traversed.count(forward)) {
            // Already have this node, so don't start a search from it.
            return;
        }
        
        // The stack only holds locally forward handles
        vector<handle_t> stack{forward};
        to_return.emplace_back();
        while (!stack.empty()) {
            handle_t here = stack.back();
            stack.pop_back();
            
            traversed.insert(here);
            to_return.back().first.insert(graph->get_id(here));
            
            // We have a counter for the number of edges we saw.
            // If it is 0 after traversing edges we know we have a tip.
            size_t edge_counter = 0;
            
            // We have a function to handle all connected handles
            auto handle_other = [&](const handle_t& other) {
                // Again, make it forward
                auto other_forward = graph->forward(other);
                
                if (!traversed.count(other_forward)) {
                    stack.push_back(other_forward);
                }
                
                edge_counter++;
            };
            
            // Look at edges in both directions
            graph->follow_edges(here, false, handle_other);
            if (edge_counter == 0) {
                // This is a tail node. Put it in reverse as a tip.
                to_return.back().second.push_back(graph->flip(here));
            }
            
            edge_counter = 0;
            graph->follow_edges(here, true, handle_other);
            if (edge_counter == 0) {
                // This is a head node. Put it as a tip.
                to_return.back().second.push_back(here);
            }
            
        }
    });
    return to_return;
}

}
}
