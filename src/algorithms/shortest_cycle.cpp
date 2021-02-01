#include "shortest_cycle.hpp"

namespace vg {
namespace algorithms {

    /// An implementation of Bellman-Ford with Yen's ordering improvement applied
    /// to a layout ideally has a small feedback arc set
    size_t bellman_ford_shortest_cycle_length(const HandleGraph* graph,
                                              const handle_t& source,
                                              const vector<handle_t>& layout,
                                              const unordered_map<handle_t, size_t>& handle_index,
                                              const vector<pair<size_t, size_t>>& feedback_edges) {
        
        // init a dynamic programming structure
        vector<size_t> dp_length(layout.size(), numeric_limits<size_t>::max());
        
        // base case
        size_t source_idx = handle_index.at(source);
        dp_length[source_idx] = 0;
        
        size_t cycle_length = numeric_limits<size_t>::max();
        
        // use dynamic programming over an implicitly dagified graph
        bool any_changed = true;
        for (int64_t i = 0; i < feedback_edges.size() + 1 && any_changed; i++) {
            any_changed = false;
            // iterate over forward edges
            for (const handle_t& handle : layout) {
                size_t idx_from = handle_index.at(handle);
                if (dp_length[idx_from] == numeric_limits<size_t>::max()) {
                    // this node hasn't been reached yet (this saves us checking for overflow)
                    continue;
                }
                size_t dist_thru = dp_length[idx_from] + graph->get_length(handle);
                graph->follow_edges(handle, false, [&](const handle_t& next) {
                    size_t idx_to = handle_index.at(next);
                    if (idx_from < idx_to) {
                        if (idx_to == source_idx) {
                            if (dist_thru < cycle_length) {
                                cycle_length = dist_thru;
                                any_changed = true;
                            }
                        }
                        else {
                            if (dist_thru < dp_length[idx_to]) {
                                dp_length[idx_to] = dist_thru;
                                any_changed = true;
                            }
                        }
                    }
                });
            }
            
            // iterate over feedback edges
            for (const pair<size_t, size_t>& feedback_edge : feedback_edges) {
                if (dp_length[feedback_edge.first] == numeric_limits<size_t>::max()) {
                    // this node hasn't been reached yet (this saves us checking for overflow)
                    continue;
                }
                size_t dist_thru = dp_length[feedback_edge.first] + graph->get_length(layout[feedback_edge.first]);
                if (feedback_edge.second == source_idx) {
                    if (dist_thru < cycle_length) {
                        cycle_length = dist_thru;
                        any_changed = true;
                    }
                }
                else {
                    if (dist_thru < dp_length[feedback_edge.second]) {
                        dp_length[feedback_edge.second] = dist_thru;
                        any_changed = true;
                    }
                }
            }
        }
        
        return cycle_length;
    }
    
    /// Simple Dijkstra implementation that computes shortest cycle
    size_t dijkstra_shortest_cycle_length(const HandleGraph* graph, const handle_t& source) {

        // distance from start of source to incoming side of the handle
        unordered_map<handle_t, size_t> distance_to;
        
        // init the queue
        structures::RankPairingHeap<handle_t, size_t, greater<size_t>> queue;
        queue.push_or_reprioritize(source, 0);
        
        // Dijkstra traversal over entire graph
        while (!queue.empty()) {
            pair<handle_t, size_t> here = queue.top();
            queue.pop();
            
            distance_to[here.first] = here.second;
            
            size_t dist_thru = here.second + graph->get_length(here.first);
            graph->follow_edges(here.first, false, [&](const handle_t& next) {
                queue.push_or_reprioritize(next, dist_thru);
            });
        }
        
        // walk one step in the other direction to complete the cycle
        size_t cycle_length = numeric_limits<size_t>::max();
        graph->follow_edges(source, true, [&](const handle_t& prev) {
            auto iter = distance_to.find(prev);
            if (iter != distance_to.end()) {
                cycle_length = min(cycle_length, iter->second + graph->get_length(prev));
            }
        });
        return cycle_length;
    }
    
    size_t shortest_cycle_length_internal(const HandleGraph* graph,
                                          const handle_t& source,
                                          const vector<handle_t>& layout,
                                          const unordered_map<handle_t, size_t>& handle_index,
                                          const vector<pair<size_t, size_t>>& feedback_edges) {
        
        size_t log_node_size = 0;
        {
            size_t log_counter = layout.size();
            while (log_counter) {
                log_node_size++;
                log_counter /= 2;
            }
        }
        
        // the Bellman-Ford implementation has run time proportional to the number of feedback
        // arcs, so it has an advantage over Dijkstra if it is dominated by log |V|
        if (feedback_edges.size() < log_node_size) {
            return bellman_ford_shortest_cycle_length(graph, source, layout, handle_index, feedback_edges);
        }
        else {
            return dijkstra_shortest_cycle_length(graph, source);
        }
    }
    
    size_t shortest_cycle_length(const HandleGraph* graph, const handle_t& source) {
        
        // compute a small FAS layout
        vector<handle_t> layout = handlealgs::eades_algorithm(graph);
        
        // identify each handle with its index in the layout
        unordered_map<handle_t, size_t> handle_index;
        for (size_t i = 0; i < layout.size(); i++) {
            handle_index[layout[i]] = i;
        }
        
        // collect the backward facing edges
        vector<pair<size_t, size_t>> feedback_edges;
        for (const handle_t& handle : layout) {
            size_t idx_from = handle_index[handle];
            graph->follow_edges(handle, false, [&](const handle_t& next) {
                size_t idx_to = handle_index[next];
                if (idx_from >= idx_to) {
                    feedback_edges.emplace_back(idx_from, idx_to);
                }
            });
        }
        
        return shortest_cycle_length_internal(graph, source, layout, handle_index, feedback_edges);
    }

    size_t shortest_cycle_length(const HandleGraph* graph) {
        
        // compute a small FAS layout
        vector<handle_t> layout = handlealgs::eades_algorithm(graph);
        
        // identify each handle with its index in the layout
        unordered_map<handle_t, size_t> handle_index;
        for (size_t i = 0; i < layout.size(); i++) {
            handle_index[layout[i]] = i;
        }
        
        // collect the backward facing edges
        vector<pair<size_t, size_t>> feedback_edges;
        for (const handle_t& handle : layout) {
            size_t idx_from = handle_index[handle];
            graph->follow_edges(handle, false, [&](const handle_t& next) {
                size_t idx_to = handle_index[next];
                if (idx_from >= idx_to) {
                    feedback_edges.emplace_back(idx_from, idx_to);
                }
            });
        }
        
        size_t min_cycle_length = numeric_limits<size_t>::max();
        
        // TODO: it shouldn't be necessary to do this on all nodes
        for (const handle_t& handle : layout) {
            size_t cycle_length = shortest_cycle_length_internal(graph,
                                                                 handle,
                                                                 layout,
                                                                 handle_index,
                                                                 feedback_edges);
            min_cycle_length = min(min_cycle_length, cycle_length);
        }
        
        return min_cycle_length;
    }
}
}

