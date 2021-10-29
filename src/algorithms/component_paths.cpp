/** \file
 * Implements the connected component paths function
 */

#include <queue>
#include <thread>
#include <atomic>
#include <mutex>

#include "../cluster.hpp"
#include "component_paths.hpp"
#include "utility.hpp"
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

vector<unordered_set<path_handle_t>> component_paths(const PathHandleGraph& graph) {
    
    vector<unordered_set<path_handle_t>> return_val;
    
    nid_t min_id = graph.min_node_id();
    nid_t id_range = graph.max_node_id() - min_id + 1;
    size_t num_nodes = graph.get_node_count();
    double id_density = double(num_nodes) / double(id_range);
    int thread_count = get_thread_count();
    if (id_density < 1.0 / double(thread_count)) {
        // the coverage of the node ID space is too diffuse for this algorithm
        // to work well
        return_val = component_path(graph);
    }
    else {
        // we use a technique to generate a high-entropy shuffle in O(1) space
        // and time per
        
        int which_prime = 0;
        while (spaced_primes[which_prime] <= id_range + 1) {
            ++which_prime;
        }
        
        size_t prime = spaced_primes[which_prime];
        size_t primitive_root = primitive_roots_of_unity[which_prime];
        size_t permuted = 1;
        
        atomic<size_t> num_visited(0);
        
        vector<int> visitor(num_nodes, -1);
        vector<int> thread_current_id()
        
        mutex get_seed_mutex;
        vector<thread> workers;
        for (int i = 0; i < thread_count; ++i) {
            workers.emplace_back([&]() {
                while (num_visited.load() < num_nodes) {
                    
                    get_seed_mutex.lock();
                    
                    while (true) {
                        id_t node_id = min_node_id + permuted - 1;
                        //
                        permuted = (permuted * primitive_root) % prime;
                    }
                    
                    get_seed_mutex.unlock();
                    
                    
                }
            });
        }
    }
}
    
}
}
