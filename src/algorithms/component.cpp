/** \file
 * Implements per-component algorithms
 */

#include <queue>
#include <thread>
#include <atomic>
#include <mutex>
#include <structures/union_find.hpp>

#include "../cluster.hpp"
#include "component.hpp"
#include "utility.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"

//#define debug_component_paths
//#define debug_parallel_component_paths

namespace vg {
namespace algorithms {

using namespace std;

// internal generic BFS implementation
void traverse_components(const HandleGraph& graph,
                         function<void(void)>& on_new_comp,
                         function<void(handle_t)>& on_node) {
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
        on_new_comp();
        
        // init a BFS queue
        deque<handle_t> queue;
        
        // function to call on each subsequent handle we navigate to
        function<bool(const handle_t&)> on_node_and_enqueue = [&](const handle_t& here) {
            
#ifdef debug_component_paths
            cerr << "traverse to handle on node " << graph.get_id(here) << endl;
#endif
            
            // don't queue up the same node twice
            if (!enqueued[graph.get_id(here) - min_id]) {
                
                // add the paths of the new node
                on_node(here);
                
                // and add it to the queue
                queue.push_back(here);
                enqueued[graph.get_id(here) - min_id] = 1;
            }
            return true;
        };
        
        // queue up the first node
        on_node_and_enqueue(handle);
        
        // do the BFS traversal
        while (!queue.empty()) {
            handle_t handle = queue.front();
            queue.pop_front();
            
            // traverse in both directions
            graph.follow_edges(handle, false, on_node_and_enqueue);
            graph.follow_edges(handle, true, on_node_and_enqueue);
        }
    });
}

vector<size_t> component_sizes(const HandleGraph& graph) {
    
    vector<size_t> comp_sizes;
    
    // add a new size for a new component
    function<void(void)> on_new_comp = [&]() {
        comp_sizes.push_back(0);
    };
    
    // add the paths of the new node
    function<void(handle_t)> on_node = [&](const handle_t here) {
        comp_sizes.back()++;
    };
    
    return comp_sizes;
}


vector<unordered_set<path_handle_t>> component_paths(const PathHandleGraph& graph) {
    
    vector<unordered_set<path_handle_t>> component_path_sets;
    
    // add a new set for a new component
    function<void(void)> on_new_comp = [&]() {
        component_path_sets.emplace_back();
    };
    
    // add the paths of the new node
    function<void(handle_t)> on_node = [&](const handle_t here) {
        graph.for_each_step_on_handle(here, [&](const step_handle_t& step) {
            component_path_sets.back().insert(graph.get_path_handle_of_step(step));
            
#ifdef debug_component_paths
            cerr << "found path " << graph.get_path_name(graph.get_path_handle_of_step(step)) << endl;
#endif
        });
    };
    
    traverse_components(graph, on_new_comp, on_node);
    
    return component_path_sets;
}

template<typename Int1, typename Int2>
void reallocate_atomic_int_vector(vector<atomic<Int1>>*& vec1, vector<atomic<Int2>>*& vec2) {
    
    if (!vec1) {
        // the vector has already been reallocated
        return;
    }
    
    // allocate the second vector
    vec2 = new vector<atomic<Int2>>(vec1->size());
    
    // TODO: these threads will actually be fighting for processor time
    // with the spin-locks holding the main threads busy while they wait...
    
    // parallel copy in contiguous blocks
    int thread_count = get_thread_count();
    static const int64_t block_size = 4096;
    atomic<int64_t> next(0);
    vector<thread> workers;
    for (int i = 0; i < thread_count; ++i) {
        workers.emplace_back([&]() {
            while (true) {
                int64_t begin = block_size * (next++);
                if (begin >= vec2->size()) {
                    // we're past the end of the vector
                    break;
                }
                for (int64_t j = begin, n = min<int64_t>(begin + block_size, vec2->size()); j < n; ++j) {
                    (*vec2)[j].store((*vec1)[j].load());
                }
            }
        });
    }
    
    // barrier sync
    for (auto& worker : workers) {
        worker.join();
    }
    
    // free the first vector
    delete vec1;
    vec1 = nullptr;
};

vector<unordered_set<path_handle_t>> component_paths_parallel(const PathHandleGraph& graph) {

#ifdef debug_parallel_component_paths
    cerr << "computing component paths in parallel" << endl;
#endif
    
    // get all paths
    vector<path_handle_t> paths;
    paths.reserve(graph.get_path_count());
    graph.for_each_path_handle([&](const path_handle_t& path) {
        paths.emplace_back(path);
    });
    
    // sort in descending order by step count
    stable_sort(paths.begin(), paths.end(), [&](path_handle_t a, path_handle_t b) {
        return graph.get_step_count(a) > graph.get_step_count(b);
    });
    
    int thread_count = get_thread_count();
    
    // the number of threads that haven't exited
    atomic<int> threads_active(thread_count);
    
    // a system that lets one thread freeze the others at checkpoints while it does
    // some job
    atomic<int> frozen(0);
    atomic<int64_t> num_frozen(0);
    
    // checkpoint to wait if any other thread wants us to freeze
    auto check_freeze = [&]() {
        if (frozen.load()) {
            ++num_frozen;
            while (frozen.load()) {
                // spin lock
            }
            --num_frozen;
        }
    };
    // wait until all other threads have reached a freeze check point
    // and then execute a function
    auto freeze_and_execute = [&](const function<void(void)>& exec) {
        // continue trying to freeze until we actually get to be the thread
        // the freezes the other threads
        bool was_frozen = true;
        while (was_frozen) {
            was_frozen = frozen.fetch_or(1);
            if (was_frozen) {
                // we weren't the thread who switched the freeze on, freeze
                // ourselves and wait
                check_freeze();
            }
            else {
                while (num_frozen.load() < threads_active.load() - 1) {
                    // spin lock waiting for the other threads to reach a checkpoint
                }
                // execute the function
                exec();
                // unfreeze the rest of the threads
                frozen.fetch_and(0);
                // leave the loop
            }
        }
    };
    
    // we'll be using the ID space as vector indices, calculat evalues we'll use for that
    nid_t min_id = graph.min_node_id();
    nid_t id_range = graph.max_node_id() - min_id + 1;
    
    // we'll try to accomplish the job with the minimum int size possible to
    // keep the memory use down
    //   note: this has to be done on the heap because the deleted copy assignment
    //   and copy construction for atomic ints means vectors can never be resized
    vector<atomic<uint8_t>>* id_vec_8 = new vector<atomic<uint8_t>>(id_range);
    vector<atomic<uint16_t>>* id_vec_16 = nullptr;
    vector<atomic<uint32_t>>* id_vec_32 = nullptr;
    // in parallel initialize with sentinels, which will be replaced by search IDs
    static const size_t block_size = 4096;
    atomic<size_t> block_idx(0);
    vector<thread> initializers;
    for (int i = 0; i < thread_count; ++i) {
        initializers.emplace_back([&]() {
            while (true) {
                size_t begin = block_size * (block_idx++);
                if (begin >= id_vec_8->size()) {
                    break;
                }
                for (size_t j = begin, end = min<size_t>(begin + block_size, id_vec_8->size()); j < end; ++j) {
                    (*id_vec_8)[j].store(0);
                }
            }
        });
    }
    
    // barrier sync
    for (auto& initializer : initializers) {
        initializer.join();
    }
    
    // this value keeps track of which one of these we're actually using, taking
    // the values 0, 1, or 2
    uint8_t which_vec = 0;
    // the last search ID we can accommodate for each vector
    const uint32_t max_search_id[3] = {
        numeric_limits<uint8_t>::max(),
        numeric_limits<uint16_t>::max(),
        numeric_limits<uint32_t>::max()
    };
    
    
    // for each search ID, the other search IDs it encountered adjacent to its search
    vector<unordered_set<int32_t>> neighbors(max_search_id[0] + 1);
    // for each search ID, the path handles it fouund while traversing
    vector<unordered_set<path_handle_t>> search_path_sets(max_search_id[0] + 1);
    
    // define accessors that hide the ugliness of checking which vector we're using:
    
    // perform atomic load on whichever vector we're currently using
    auto load = [&](int64_t i) {
        uint32_t loaded;
        switch (which_vec) {
            case 0:
                loaded = (*id_vec_8)[i].load();
                break;
            case 1:
                loaded = (*id_vec_16)[i].load();
                break;
            default:
                loaded = (*id_vec_32)[i].load();
                break;
        }
        return loaded;
    };
    // perform atomic compare-exchange on whichever vector we're currently using
    auto compare_exchange = [&](int64_t i, uint32_t& expected, uint32_t desired) {
        bool exchanged;
        switch (which_vec) {
            case 0:
            {
                uint8_t expected_8 = expected;
                uint8_t desired_8 = desired;
                exchanged = (*id_vec_8)[i].compare_exchange_strong(expected_8, desired_8);
                if (!exchanged) {
                    expected = expected_8;
                }
                break;
            }
            case 1:
            {
                uint16_t expected_16 = expected;
                uint16_t desired_16 = desired;
                exchanged = (*id_vec_16)[i].compare_exchange_strong(expected_16, desired_16);
                if (!exchanged) {
                    expected = expected_16;
                }
                break;
            }
            default:
            {
                exchanged = (*id_vec_32)[i].compare_exchange_strong(expected, desired);
                break;
            }
        }
        return exchanged;
    };
    
    
    // to keep track of the index of the path we will use to seed a BFS next
    atomic<int64_t> next_path(0);
    // to keep track of the ID of the next seeded search
    atomic<uint32_t> next_search_id(1);
    // initialize the swarm of workers
    vector<thread> workers;
    for (int i = 0; i < thread_count; ++i) {
        workers.emplace_back([&,i]() {
            while (true) {
                
                int64_t path_idx = next_path++;
                
                if (path_idx >= paths.size()) {
                    // all of the paths have been explored, we can exit
                    break;
                }
                
#ifdef debug_parallel_component_paths
                cerr << ("worker " + to_string(i) + " got path idx " + to_string(path_idx) + ": " + graph.get_path_name(paths[path_idx]) + " with step count " + to_string(graph.get_step_count(paths[path_idx])) + "\n");
#endif
                
                path_handle_t path = paths[path_idx];
                if (graph.get_step_count(path) == 0) {
                    // skip an empty path
                    check_freeze();
                    continue;
                }
                
                // seed a BFS search off of the first node in this path
                handle_t seed = graph.get_handle_of_step(graph.path_begin(path));
                
                if (load(graph.get_id(seed) - min_id) != 0) {
                    // another thread has already traversed over this node, no need
                    // to start a search here
#ifdef debug_parallel_component_paths
                    cerr << ("worker " + to_string(i) + " skipping seed " + to_string(graph.get_id(seed)) + ", which was previously visited by " + to_string(load(graph.get_id(seed) - min_id)) + "\n");
#endif
                    check_freeze();
                    continue;
                }
                
                
                
                // we're going to initiate a BFS from the seed, assign a new search ID
                uint32_t search_id = next_search_id++;
                
                // TODO: add finer-grain reallocations so that neighbors and
                // search_path_sets don't need to get so large to guarantee that
                // we don't index past them with a search ID
                if (search_id > max_search_id[which_vec]) {
                    // we need to move up to the next int size in order to acommodate
                    // this search ID, so demand that the other threads stop writing so
                    // we can switch to a larger bit width
                    freeze_and_execute([&]() {
                        // check to make sure another thread didn't already move us over
                        // to the next vector while we were waiting to freeze
                        if (search_id <= max_search_id[which_vec]) {
                            return;
                        }
                        
                        ++which_vec;
                        neighbors.resize(max_search_id[which_vec] + 1);
                        search_path_sets.resize(max_search_id[which_vec] + 1);
                        if (which_vec == 1) {
                            reallocate_atomic_int_vector(id_vec_8, id_vec_16);
                        }
                        else if (which_vec == 2) {
                            reallocate_atomic_int_vector(id_vec_16, id_vec_32);
                        }
                        else {
                            cerr << "error: parallel component paths algorithm ran out of 32-bit search IDs\n";
                            exit(1);
                        }
                    });
                }
                
#ifdef debug_parallel_component_paths
                cerr << ("worker " + to_string(i) + " starting search on seed " + to_string(graph.get_id(seed)) + " with search ID " + to_string(search_id) + "\n");
#endif
                
                // FIFO queue for BFS
                deque<handle_t> queue;
                
                // function to call on each subsequent handle we navigate to
                function<bool(const handle_t&)> record_paths_and_enqueue = [&](const handle_t& here) {
                    int64_t idx = graph.get_id(here) - min_id;
                    uint32_t visit_id = 0;
                    bool exchanged = compare_exchange(idx, visit_id, search_id);
                    if (exchanged) {
                        // we found the unvisited sentinel and replaced it with our search ID
                        
                        // add the paths of the new node
                        graph.for_each_step_on_handle(here, [&](const step_handle_t& step) {
                            search_path_sets[search_id].insert(graph.get_path_handle_of_step(step));
                        });
                        
                        // and add it to the queue
                        queue.push_back(here);
                    }
                    else if (visit_id != search_id) {
                        // we are adjacent to nodes explored by a different search, record the
                        // neighbor
                        neighbors[search_id].insert(visit_id);
                    }
                    return true;
                };
                
                // init the queue
                record_paths_and_enqueue(seed);
                
                while (!queue.empty()) {
                    // set a checkpoint in case a thread is trying to reallocate
                    check_freeze();
                    
                    // de-queue a node
                    handle_t handle = queue.front();
                    queue.pop_front();
                    
                    // traverse in both directions
                    graph.follow_edges(handle, false, record_paths_and_enqueue);
                    graph.follow_edges(handle, true, record_paths_and_enqueue);
                }
                
            }
            // keep track of the fact this thread is exiting
            --threads_active;
#ifdef debug_parallel_component_paths
            cerr << ("worker " + to_string(i) + " exiting\n");
#endif
        });
    }
    
    // barrier sync
    for (auto& worker : workers) {
        worker.join();
    }
    
    // find equivalence classes of search IDs on the same component
    size_t num_search_ids = next_search_id.load();
    structures::UnionFind union_find(num_search_ids, false);
    for (size_t i = 0; i < num_search_ids; ++i) {
        for (auto j : neighbors[i]) {
            union_find.union_groups(i, j);
        }
    }
    // agglomerate the sets of paths ecountered by all of the searches
    // in each equivalence class
    vector<unordered_set<path_handle_t>> return_val;
    for (const auto& component_search_ids : union_find.all_groups()) {
        bool added_new_set = false;
        for (size_t i : component_search_ids) {
            for (path_handle_t path : search_path_sets[i]) {
                if (!added_new_set) {
                    return_val.emplace_back();
                    added_new_set = true;
                }
                return_val.back().insert(path);
            }
        }
    }
    
    delete id_vec_8;
    delete id_vec_16;
    delete id_vec_32;
    
    return return_val;
}
    
}
}
