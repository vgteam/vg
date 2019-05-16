#include "is_single_stranded.hpp"

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
    
    vector<handle_t> single_stranded_orientation(const HandleGraph* graph) {
        
        // the return value
        vector<handle_t> orientation;
        orientation.reserve(graph->get_node_count());
        
        // keep track of which nodes have already been oriented and which orientation
        unordered_map<id_t, bool> recorded_orientation;
        
        // keep track of whether we've encountered a node in two orientations
        bool failed = false;
        
        // DFS through the graph
        graph->for_each_handle([&](const handle_t& handle) {
            if (recorded_orientation.count(graph->get_id(handle))) {
                return true;
            }
            
            // initialize the stack
            vector<handle_t> stack(1, handle);
            
            // record the orientation of the seed for the traversal
            orientation.push_back(handle);
            recorded_orientation[graph->get_id(handle)] = graph->get_is_reverse(handle);
            
            function<bool(const handle_t&)> walk_edges = [&](const handle_t& next) {
                auto iter = recorded_orientation.find(graph->get_id(next));
                if (iter != recorded_orientation.end()) {
                    // we've been here before, but make sure we're encountering it in the same orientation
                    failed = (iter->second != graph->get_is_reverse(next));
                }
                else {
                    // add to the DFS stack
                    stack.push_back(next);
                    
                    // record the orientation we encountered it in
                    orientation.push_back(next);
                    recorded_orientation[graph->get_id(next)] = graph->get_is_reverse(next);
                }
                // continue if we haven't failed
                return !failed;
            };
            
            // DFS
            while (!stack.empty() && !failed) {
                handle_t here = stack.back();
                stack.pop_back();
                
                graph->follow_edges(here, true, walk_edges);
                if (!failed) {
                    graph->follow_edges(here, false, walk_edges);
                }
            }
            
            // continue if there's any more to do and we haven't failed
            return orientation.size() < graph->get_node_count() && !failed;
        });
        
        // if we failed, we return an empty vector as a sentinel
        if (failed) {
            orientation.clear();
        }
        
        return orientation;
    }

    
    unordered_set<id_t> make_single_stranded(MutableHandleGraph* graph) {
        
        auto orientations = single_stranded_orientation(graph);
        
        if (orientations.size() != graph->get_node_count()) {
            // we got the sentinel for an un-single-strandable graph
            cerr << "error:[algorithms] attempted to apply single-stranded orientation to non-single stranded graph" << endl;
            exit(1);
        }
        
        return apply_orientations(graph, orientations);
    }
}
}
