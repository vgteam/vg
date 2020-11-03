#include "path_component_index.hpp"

#include <queue>
#include "sdsl/bit_vectors.hpp"

//#define debug_component_index

namespace vg {

    PathComponentIndex::PathComponentIndex() {
        // Nothing to do
    }
    
    PathComponentIndex::PathComponentIndex(const PathHandleGraph* graph) {
        
        // we will use this as a bit-flag for when a node id has already been
        // traverseed
        id_t min_id = graph->min_node_id();
        sdsl::bit_vector enqueued(graph->max_node_id() - min_id + 1, 0);
        
        graph->for_each_handle([&](const handle_t& handle) {
            
            if (enqueued[graph->get_id(handle) - min_id]) {
                // we've already traversed this node
                return;
            }
            
#ifdef debug_component_index
            cerr << "starting new component on node " << graph->get_id(handle) << endl;
#endif
            
            // a node that hasn't been traversed means a new component
            component_path_sets.emplace_back();
            unordered_set<path_handle_t>& component_path_set = component_path_sets.back();
            
            // init a BFS queue
            std::queue<handle_t> queue;
            
            // function to call on each subsequent handle we navigate to
            function<bool(const handle_t&)> record_paths_and_enqueue = [&](const handle_t& here) {
                
#ifdef debug_component_index
                cerr << "traverse to handle on node " << graph->get_id(here) << endl;
#endif
                
                // don't queue up the same node twice
                if (!enqueued[graph->get_id(here) - min_id]) {
                    
                    // add the paths of the new node
                    graph->for_each_step_on_handle(here, [&](const step_handle_t& step) {
                        component_path_set.insert(graph->get_path_handle_of_step(step));
#ifdef debug_component_index
                        cerr << "found path " << graph->get_path_name(graph->get_path_handle_of_step(step)) << endl;
#endif
                    });
                    
                    // and add it to the queue
                    queue.push(here);
                    enqueued[graph->get_id(here) - min_id] = 1;
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
                graph->follow_edges(handle, false, record_paths_and_enqueue);
                graph->follow_edges(handle, true, record_paths_and_enqueue);
            }
        });
        
        // make it so we can index into this with the path rank directly
        component_path_set_of_path.reserve(graph->get_path_count());
        
        // index from the paths to their component set
        for (size_t i = 0; i < component_path_sets.size(); i++) {
            for (const path_handle_t& path : component_path_sets[i]) {
                if (component_path_set_of_path.count(path)) {
                    cerr << "warning:[PathComponentIndex] Graph contains path " << graph->get_path_name(path) << " that spans multiple connected components. This path must follow edges that are not included in the graph. The PathComponentIndex may not be semantically meaningful for this graph." << endl;
                    continue;
                }
                component_path_set_of_path[path] = i;
            }
        }
    }
    
    bool PathComponentIndex::paths_on_same_component(const path_handle_t& path_1,
                                                     const path_handle_t& path_2) const {
        
        return component_path_sets.at(component_path_set_of_path.at(path_1)).count(path_2);
    }
}
