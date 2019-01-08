#include "dagify.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    unordered_map<id_t, id_t> dagify(const HandleGraph* graph, MutableHandleGraph* into,
                                     size_t min_preserved_path_length) {
        
        // initialize the translator from the dagified graph back to the original graph
        unordered_map<id_t, id_t> trans;
        
        // generate a canonical orientation across the graph
        vector<handle_t> orientation = single_stranded_orientation(graph);
        
        // mark the ones that whose canonical orientation is reversed
        unordered_set<id_t> reversed_nodes;
        for (size_t i = 0; i < orientation.size(); i++) {
            if (graph->get_is_reverse(orientation[i])) {
                reversed_nodes.insert(graph->get_id(orientation[i]));
            }
        }
        
        // find the strongly connected components of the original graph
        vector<unordered_set<id_t>> strong_components = strongly_connected_components(graph);
        
        // create the nodes of the dagified graph
        unordered_map<id_t, size_t> component_of;
        unordered_map<handle_t, vector<handle_t>> injector;
        for (size_t i = 0; i < strong_components.size(); i++) {
            
            // keep track of which nodes are in which component
            auto& component = strong_components[i];
            for (id_t node_id : component) {
                component_of[node_id] = i;
            }
            
            // figure out how many times we need to copy this SCC
            SubHandleGraph subgraph(graph);
            for (const id_t& node_id : component) {
                subgraph.add_handle(graph->get_handle(node_id));
            }
            
            // get a layout with a low FAS
            vector<handle_t> layout = eades_algorithm(&subgraph);
            
            // make sure the layout matches the canonical orientation of the graph
            if (graph->get_is_reverse(layout.front()) != reversed_nodes.count(graph->get_id(layout.front()))) {
                for (int64_t i = 0, j = layout.size() - 1; i < j; i++, j--) {
                    auto tmp = layout[i];
                    layout[i] = subgraph.flip(layout[j]);
                    layout[j] = subgraph.flip(tmp);
                }
            }
            
            // record the ordering of the layout
            unordered_map<handle_t, size_t> ordering;
            for (size_t i = 0; i < layout.size(); i++) {
                ordering[layout[i]] = i;
            }
            
            // mark the edges as either forward or backward relative to the layout
            vector<vector<size_t>> forward_edges(layout.size());
            vector<pair<size_t, size_t>> backward_edges;
            // also collect the nodes that are the head of a backward-facing edge
            unordered_set<handle_t> backward_edge_heads;
            subgraph.for_each_edge([&](const edge_t& edge) {
                // get the indices of the edge in the layout, according to the canonical
                // orientation
                size_t i, j;
                auto iter = ordering.find(edge.first);
                if (iter != ordering.end()) {
                    i = iter->second;
                    j = ordering[edge.second];
                }
                else {
                    i = ordering[subgraph.flip(edge.second)];
                    j = ordering[subgraph.flip(edge.first)];
                }
                
                // classify the edge as forward or backward
                if (i < j) {
                    forward_edges[i].push_back(j);
                }
                else {
                    backward_edges.emplace_back(i, j);
                    backward_edge_heads.insert(layout[i]);
                }
                
                // always keep going
                return true;
            });
            
            // check for each node whether we've duplicated the component enough times
            // to preserve its cycles
            for (const handle_t& handle : backward_edge_heads) {
                
                vector<int64_t> distances(layout.size(), numeric_limits<int64_t>::max());
                // start with negative distance so that we are measuring from the end of the node
                distances[ordering[handle]] = -subgraph.get_length(handle);
                
                // keep going until the distance to the next copy unit is too large
                int64_t min_relaxed_dist = -subgraph.get_length(handle);
                for (size_t copy_num = 0; min_relaxed_dist <= int64_t(min_preserved_path_length); copy_num++) {
                    
                    // do we need a new copy of this SCC to preserve paths?
                    if (copy_num == injector[layout.front()].size()) {
                        // we haven't added this copy of the connected component yet
                        
                        // add the nodes
                        for (const handle_t& original_handle : layout) {
                            // create the node in the same foward orientation as the original
                            handle_t new_handle = into->create_handle(graph->get_sequence(graph->forward(original_handle)));
                            // use the handle locally in the same orientation as it is in the layout
                            if (graph->get_is_reverse(original_handle)) {
                                new_handle = into->flip(new_handle);
                            }
                            trans[into->get_id(new_handle)] = graph->get_id(original_handle);
                            injector[handle].push_back(new_handle);
                        }
                        
                        // add the forward edges within this copy
                        for (size_t i = 0; i < forward_edges.size(); i++) {
                            handle_t from = injector[layout[i]].back();
                            for (const size_t& j : forward_edges[i]) {
                                into->create_edge(from, injector[layout[j]].back());
                            }
                        }
                        
                        if (copy_num > 0) {
                            // add the backward edges between the copies
                            for (const pair<size_t, size_t>& bwd_edge : backward_edges) {
                                const auto& from_copies = injector[layout[bwd_edge.first]];
                                into->create_edge(from_copies[from_copies.size() - 2],
                                                  injector[layout[bwd_edge.second]].back());
                            }
                        }
                    }
                    
                    // find the shortest path to the nodes, staying within this copy
                    // of the SCC
                    for (size_t i = 0; i < distances.size(); i++) {
                        // skip infinity to avoid overflow
                        if (distances[i] == numeric_limits<int64_t>::max()) {
                            continue;
                        }
                        
                        int64_t dist_thru = distances[i] + subgraph.get_length(layout[i]);
                        for (const size_t& j : forward_edges[i]) {
                            distances[j] = min(distances[j], dist_thru);
                        }
                    }
                    
                    // now find the minimum distance to nodes in the next copy of the SCC (which
                    // may not yet be created in the graph)
                    min_relaxed_dist = numeric_limits<int64_t>::max();
                    vector<int64_t> next_distances(distances.size(), numeric_limits<int64_t>::max());
                    for (const pair<size_t, size_t>& bwd_edge : backward_edges) {
                        if (distances[bwd_edge.first] == numeric_limits<int64_t>::max()) {
                            continue;
                        }
                        
                        int64_t dist_thru = distances[bwd_edge.first] + subgraph.get_length(layout[bwd_edge.first]);
                        if (dist_thru < next_distances[bwd_edge.second]) {
                            next_distances[bwd_edge.second] = dist_thru;
                            min_relaxed_dist = min(min_relaxed_dist, dist_thru);
                        }
                    }
                    
                    distances = move(next_distances);
                }
            }
        }
        
//        // compute the edges between the connected components (i.e. the condensation graph)
//        vector<unordered_set<size_t>> condensation_graph(strong_components.size());
//        for (const handle_t& handle : orientation) {
//            size_t comp_here = component_of[graph->get_id(handle)];
//            graph->follow_edges(handle, false, [&](const handle_t& next) {
//                size_t comp_next = component_of[graph->get_id(next)];
//                if (comp_here != comp_next) {
//                    condensation_graph[comp_here].insert(comp_next);
//                }
//            });
//        }
//
//        // use Kahn's algorithm to compute the topological order of the condensation graph
//        vector<size_t> topological_index(condensation_graph.size(), 0);
//
//        // compute in degrees
//        vector<size_t> in_degree(condensation_graph.size(), 0);
//        for (const unordered_set<size_t>& adj_list : condensation_graph) {
//            for (const size_t& edge : adj_list) {
//                in_degree[edge]++;
//            }
//        }
//        // init stack
//        vector<size_t> stack;
//        for (size_t i = 0; i < in_degree.size(); i++) {
//            if (in_degree[i] == 0) {
//                stack.push_back(i);
//            }
//        }
//        // compute topological order
//        for (size_t next_top_idx = 0; next_top_idx < condensation_graph.size(); next_top_idx++) {
//            size_t here = stack.back();
//            stack.pop_back();
//            topological_index[here] = next_top_idx;
//            for (const size_t& next : condensation_graph[here]) {
//                in_degree[next]--;
//                if (in_degree[next] == 0) {
//                    stack.push_back(next);
//                }
//            }
//        }
        
        // add edges between the strongly connected components
        graph->for_each_edge([&](const edge_t& canonical_edge) {
            // put the edge in the order of the orientation we've imposed on the graph
            edge_t edge = (graph->get_is_reverse(canonical_edge.first) != reversed_nodes.count(graph->get_id(canonical_edge.first)) ?
                           make_pair(graph->flip(canonical_edge.second), graph->flip(canonical_edge.first)) :
                           canonical_edge);

            if (component_of[graph->get_id(edge.first)] != component_of[graph->get_id(edge.second)]) {
                // this edge is between SCCs
                
                // connect the last copy of the first node to all copies of the second
                const handle_t& from = injector[edge.first].back();
                for (const handle_t& to : injector[edge.second]) {
                    into->create_edge(from, to);
                }
            }
            // always keep going
            return true;
        });
        
        // return the ID translator
        return trans;
    }
}
}

