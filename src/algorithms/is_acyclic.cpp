#include "is_acyclic.hpp"

#include <unordered_map>

namespace vg {
namespace algorithms {

using namespace std;

bool is_acyclic(const HandleGraph* graph) {
    
    // the existence of reversing cycles is equivalent to whether a single stranded
    // orientation exists
    if (single_stranded_orientation(graph).size() < graph->get_node_count()) {
        return false;
    }
    // the existence of non-reversing cycles is checked by the directed acyclic algorithm
    return is_directed_acyclic(graph);
}

bool is_directed_acyclic(const HandleGraph* graph) {
    // We track the in and out degrees of all nodes. We then clean up degrees
    // entries from tips until either we've cleaned up all the nodes or there
    // are only directed cycles left.

    // Build the degrees map
    unordered_map<id_t, pair<int64_t, int64_t>> degrees;
    // And also the stack of tips to start at
    vector<handle_t> stack;
    graph->for_each_handle([&](const handle_t& here) {
        size_t start_degree = graph->get_degree(here, true);
        size_t end_degree = graph->get_degree(here, false);
    
        degrees[graph->get_id(here)] = make_pair(start_degree, end_degree);
        
        if (start_degree == 0) {
            // Tip looking forward
            stack.push_back(here);
        }
        if (end_degree == 0) {
            // Tip looking backward
            stack.push_back(graph->flip(here));
        }
        
    });
    
    while (!stack.empty()) {
        handle_t here = stack.back();
        stack.pop_back();
        
        auto iter = degrees.find(graph->get_id(here));
        if (iter == degrees.end()) {
            // Already processed
            continue;
        }
        
        degrees.erase(iter);
        
        graph->follow_edges(here, false, [&](const handle_t& next) {
            auto next_iter = degrees.find(graph->get_id(next));
            if (next_iter != degrees.end()) {
                // We have a node next that we haven't finished yet
                
                // Reduce its degree on the appropriate side.
                int64_t& in_degree = graph->get_is_reverse(next) ? next_iter->second.second : next_iter->second.first;
                in_degree--;
                if (in_degree == 0) {
                    // This is a new tip in this orientation
                    stack.push_back(next);
                }
            }
        });
    }
    
    // If we clean up the whole graph, it must have been directed-acyclic.
    return degrees.empty();
}

}
}
