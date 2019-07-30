#include "expand_context.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    void expand_context_by_steps(const HandleGraph* source, MutableHandleGraph* subgraph,
                                 int64_t dist, bool expand_forward, bool expand_backward) {
        
        // let's make an O(1) lookup for the current edges in case the implementation doesn't provide
        // one (avoids O(N^2) behavior using has_edge)
        unordered_set<edge_t> already_included_edges;
        subgraph->for_each_edge([&](const edge_t& edge) {
            already_included_edges.insert(edge);
        });
        
        // a fifo queue for BFS
        queue<pair<handle_t, int64_t>> bfs_queue;
        // all the edges we encounter along the way
        unordered_set<edge_t> seen_edges;
        
        // initialize the queue with the current subgraph
        subgraph->for_each_handle([&](const handle_t& handle) {
            bfs_queue.emplace(source->get_handle(subgraph->get_id(handle)), 0);
        });
        
        // BFS outward
        while (!bfs_queue.empty()) {
            pair<handle_t, int64_t> here = bfs_queue.front();
            bfs_queue.pop();
            
            int64_t dist_thru = here.second + 1;
            
            if (dist_thru <= dist) {
                if (expand_forward) {
                    source->follow_edges(here.first, false, [&](const handle_t& next) {
                        seen_edges.insert(source->edge_handle(here.first, next));
                        if (!subgraph->has_node(source->get_id(next))) {
                            subgraph->create_handle(source->get_sequence(source->forward(next)),
                                                    source->get_id(next));
                            bfs_queue.emplace(next, dist_thru);
                        }
                    });
                }
                
                if (expand_backward) {
                    source->follow_edges(here.first, true, [&](const handle_t& prev) {
                        seen_edges.insert(source->edge_handle(prev, here.first));
                        if (!subgraph->has_node(source->get_id(prev))) {
                            subgraph->create_handle(source->get_sequence(source->forward(prev)),
                                                    source->get_id(prev));
                            bfs_queue.emplace(prev, dist_thru);
                        }
                    });
                }
            }
        }
        
        // add in the edges we saw
        for (const edge_t& edge : seen_edges) {
            edge_t insertable = subgraph->edge_handle(subgraph->get_handle(source->get_id(edge.first),
                                                                           source->get_is_reverse(edge.first)),
                                                      subgraph->get_handle(source->get_id(edge.second),
                                                                           source->get_is_reverse(edge.second)));
            if (!already_included_edges.count(insertable)) {
                subgraph->create_edge(insertable);
            }
        }
    }
    
    void expand_context_by_length(const HandleGraph* source, MutableHandleGraph* subgraph,
                                  int64_t dist, bool expand_forward, bool expand_backward) {
        
        // let's make an O(1) lookup for the current edges in case the implementation doesn't provide
        // one (avoids O(N^2) behavior using has_edge)
        unordered_set<edge_t> already_included_edges;
        subgraph->for_each_edge([&](const edge_t& edge) {
            already_included_edges.insert(edge);
        });
        
        // extra bool indicates whether the position indicated is at the
        // beginning of the node (beginning = true, end = false)
        structures::RankPairingHeap<pair<handle_t, bool>, int64_t, greater<int64_t>> dijk_queue;
        // all the edges we encounter along the way
        unordered_set<edge_t> seen_edges;
        
        // initialize the queue at the ends pointing out of the node in the direction
        // of our search
        subgraph->for_each_handle([&](const handle_t& handle) {
            handle_t src_handle = source->get_handle(subgraph->get_id(handle));
            if (expand_forward) {
                dijk_queue.push_or_reprioritize(make_pair(src_handle, false), 0);
            }
            if (expand_backward) {
                dijk_queue.push_or_reprioritize(make_pair(src_handle, true), 0);
            }
        });
        
        
        while (!dijk_queue.empty()) {
            // the next closest untraversed node end
            pair<pair<handle_t, bool>, int64_t> here = dijk_queue.top();
            dijk_queue.pop();
            
            if (here.second > dist) {
                break;
            }
            
            if ((here.first.second && expand_forward)
                || (!here.first.second && expand_backward)) {
                // cross the node (traverse the sequence length)
                int64_t dist_across = here.second + source->get_length(here.first.first);
                if (dist_across <= dist) {
                    dijk_queue.push_or_reprioritize(make_pair(here.first.first, !here.first.second), dist_across);
                }
            }
            if ((here.first.second && expand_backward)
                || (!here.first.second && expand_forward)) {
                // cross an edge (no added distance)
                source->follow_edges(here.first.first, here.first.second, [&](const handle_t& next) {
                    dijk_queue.push_or_reprioritize(make_pair(next, !here.first.second), here.second);
                    
                    // the edge handle will be in a different order depending on whether we're
                    // going left or right
                    if (here.first.second) {
                        seen_edges.insert(source->edge_handle(next, here.first.first));
                    }
                    else {
                        seen_edges.insert(source->edge_handle(here.first.first, next));
                    }
                    
                    if (!subgraph->has_node(source->get_id(next))) {
                        subgraph->create_handle(source->get_sequence(source->forward(next)),
                                                source->get_id(next));
                    }
                });
            }
        }
        
        // add in the edges we saw
        for (const edge_t& edge : seen_edges) {
            edge_t insertable = subgraph->edge_handle(subgraph->get_handle(source->get_id(edge.first),
                                                                           source->get_is_reverse(edge.first)),
                                                      subgraph->get_handle(source->get_id(edge.second),
                                                                           source->get_is_reverse(edge.second)));
            if (!already_included_edges.count(insertable)) {
                subgraph->create_edge(insertable);
            }
        }
    }
    
    void expand_context(const HandleGraph* source, MutableHandleGraph* subgraph,
                        int64_t dist, bool use_steps, bool expand_forward,
                        bool expand_backward) {
        
        if (use_steps) {
            expand_context_by_steps(source, subgraph, dist, expand_forward, expand_backward);
        }
        else {
            expand_context_by_length(source, subgraph, dist, expand_forward, expand_backward);
        }
    }
    
    void expand_context_with_paths(const PathHandleGraph* source,
                                   MutablePathMutableHandleGraph* subgraph,
                                   int64_t dist, bool use_steps, bool expand_forward,
                                   bool expand_backward) {
        
        // get the topology of the subgraph we want
        expand_context(source, subgraph, dist, use_steps, expand_forward, expand_backward);
        
        // find all steps on all nodes
        unordered_set<step_handle_t> seen_steps;
        subgraph->for_each_handle([&](const handle_t& handle) {
            handle_t src_handle = source->get_handle(subgraph->get_id(handle));
            source->for_each_step_on_handle(src_handle, [&](const step_handle_t& step) {
                seen_steps.insert(step);
            });
        });
        
        // keep track of how many segments we've added for each path
        unordered_map<path_handle_t, int> segment_count;
        
        while (!seen_steps.empty()) {
            // choose a segment arbitrarily
            step_handle_t segment_seed = *seen_steps.begin();
            path_handle_t path_handle = source->get_path_handle_of_step(segment_seed);
            
            
            // walk backwards until we're out of the subgraph or we loop around to
            // the same step
            auto step = segment_seed;
            step_handle_t prev;
            bool first_iter = true;
            while (seen_steps.count(step) && (first_iter || step != segment_seed)) {
                prev = step;
                step = source->get_previous_step(step);
                first_iter = false;
            }
            
            if (step == segment_seed) {
                // we walked all the way around a circular path
                assert(source->get_is_circular(path_handle));
                
                // copy the whole path, making the same step the "first"
                path_handle_t segment_handle = subgraph->create_path_handle(source->get_path_name(path_handle), true);
                source->for_each_step_in_path(path_handle, [&](const step_handle_t& step) {
                    handle_t handle = source->get_handle_of_step(step);
                    subgraph->append_step(segment_handle,
                                          subgraph->get_handle(source->get_id(handle),
                                                               source->get_is_reverse(handle)));
                    seen_steps.erase(step);
                });
            }
            else {
                // we walked off the subgraph, so this is the start of a segment
                
                // get the path name (and possibly disambiguate segments from the same path)
                string path_name = source->get_path_name(path_handle);
                if (segment_count[path_handle]) {
                    path_name += "-" + to_string(segment_count[path_handle]);
                }
                
                if (subgraph->has_path(path_name)) {
                    cerr << "error: failed to create a unique path name for segment of " << source->get_path_name(path_handle) << ", there is already path named " << path_name << endl;
                    exit(1);
                }
                
                // make a new path for this segment
                path_handle_t segment_handle = subgraph->create_path_handle(path_name);
                
                for (step = prev;                          // start back on the subgraph
                     seen_steps.count(step);               // stop once we leave the subgraph
                     step = source->get_next_step(step)) {
                    
                    handle_t src_handle = source->get_handle_of_step(step);
                    subgraph->append_step(segment_handle,
                                          subgraph->get_handle(source->get_id(src_handle),
                                                               source->get_is_reverse(src_handle)));
                    
                    // keep track of which steps we've already added
                    seen_steps.erase(step);
                }
            }
            
            // make note of the fact that we've made a segment from this source path
            segment_count[path_handle]++;
        }
    }
}
}
