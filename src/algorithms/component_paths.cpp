/** \file
 * Implements the connected component paths function
 */

#include <queue>

#include "component_paths.hpp"
#include "sdsl/bit_vectors.hpp"


namespace vg {
namespace algorithms {

using namespace std;

vector<unordered_set<path_handle_t>> component_paths(const PathHandleGraph& graph) {
    
    vector<unordered_set<path_handle_t>> component_path_sets;
    
    // we will use this as a bit-flag for when a node id has already been
    // traverseed
    id_t min_id = graph.min_node_id();
    sdsl::bit_vector enqueued(graph.max_node_id() - min_id + 1, 0);
    
    graph.for_each_handle([&](const handle_t& handle) {
        
        if (enqueued[graph.get_id(handle) - min_id]) {
            // we've already traversed this node
            return;
        }
        
#ifdef debug_component_paths
        cerr << "starting new component on node " << graph.get_id(handle) << endl;
#endif
        
        // a node that hasn't been traversed means a new component
        component_path_sets.emplace_back();
        unordered_set<path_handle_t>& component_path_set = component_path_sets.back();
        
        // init a BFS queue
        std::queue<handle_t> queue;
        
        // function to call on each subsequent handle we navigate to
        function<bool(const handle_t&)> record_paths_and_enqueue = [&](const handle_t& here) {
            
#ifdef debug_component_paths
            cerr << "traverse to handle on node " << graph.get_id(here) << endl;
#endif
            
            // don't queue up the same node twice
            if (!enqueued[graph.get_id(here) - min_id]) {
                
                // add the paths of the new node
                graph.for_each_step_on_handle(here, [&](const step_handle_t& step) {
                    component_path_set.insert(graph.get_path_handle_of_step(step));
#ifdef debug_component_paths
                    cerr << "found path " << graph.get_path_name(graph.get_path_handle_of_step(step)) << endl;
#endif
                });
                
                // and add it to the queue
                queue.push(here);
                enqueued[graph.get_id(here) - min_id] = 1;
            }
            return true;
        };
        
        // queue up the first node
        record_paths_and_enqueue(handle);
        
        // do the BFS traversal
        while (!queue.empty()) {
            handle_t handle = queue.front();
            queue.pop();
            
            // traverse in both directions
            graph.follow_edges(handle, false, record_paths_and_enqueue);
            graph.follow_edges(handle, true, record_paths_and_enqueue);
        }
    });
    
    return component_path_sets;
}
    
}
}
