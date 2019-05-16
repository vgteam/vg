/**
 * \file dagify.cpp
 *
 * Defines an algorithm to convert a single-stranded graph into a DAG
 */

//#define debug_dagify

#include "dagify.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    unordered_map<id_t, id_t> dagify(const HandleGraph* graph, MutableHandleGraph* into,
                                     size_t min_preserved_path_length) {
        
        // initialize the translator from the dagified graph back to the original graph
        unordered_map<id_t, id_t> translator;
        
        // generate a canonical orientation across the graph
        vector<handle_t> orientation = single_stranded_orientation(graph);
        
        if (orientation.size() < graph->get_node_count()) {
            cerr << "error:[dagify] Dagify algorithm only valid on graphs with a single stranded orientation, consider using split_strands first" << endl;
            exit(1);
        }
        
#ifdef debug_dagify
        cerr << "got canonical orientation:" << endl;
        for (const handle_t& h : orientation) {
            cerr << "\t" << graph->get_id(h) << (graph->get_is_reverse(h) ? "-" : "+") << endl;
        }
#endif
        
        // mark the ones that whose canonical orientation is reversed
        unordered_set<id_t> reversed_nodes;
        for (size_t i = 0; i < orientation.size(); i++) {
            if (graph->get_is_reverse(orientation[i])) {
                reversed_nodes.insert(graph->get_id(orientation[i]));
            }
        }
        
        // find the strongly connected components of the original graph
        vector<unordered_set<id_t>> strong_components = strongly_connected_components(graph);
        
#ifdef debug_dagify
        cerr << "got strongly connected components:" << endl;
        for (size_t i = 0; i < strong_components.size(); i++) {
            cerr << "\tcomponent " << i << endl;
            for (id_t nid : strong_components[i]) {
                cerr << "\t\t" << nid << endl;
            }
        }
#endif
        
        // duplicate strongly connected components into the dagified graph in such a way
        // that paths are preserved
        
        // a tracker for which SCC a node belongs to
        unordered_map<id_t, size_t> component_of;
        // a map from a node in the original graph to all its copies (in order) in the
        // dagified graph
        unordered_map<handle_t, vector<handle_t>> injector;
        for (size_t i = 0; i < strong_components.size(); i++) {
            
#ifdef debug_dagify
            cerr << "handling component " << i << endl;
#endif
            
            // keep track of which nodes are in which component (for later)
            auto& component = strong_components[i];
            for (id_t node_id : component) {
                component_of[node_id] = i;
            }
            
            // figure out how many times we need to copy this SCC
            
            // wrap the SCC in a handle graph
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
                if (layout.size() % 2) {
                    layout[layout.size() / 2] = subgraph.flip(layout[layout.size() / 2]);
                }
            }
            
#ifdef debug_dagify
            cerr << "layout for component:" << endl;
            for (const handle_t& h : layout) {
                cerr << "\t" << graph->get_id(h) << (graph->get_is_reverse(h) ? "-" : "+") << endl;
            }
#endif
            
            // record the ordering of the layout so we can identify backward edges
            unordered_map<handle_t, size_t> ordering;
            for (size_t i = 0; i < layout.size(); i++) {
                ordering[layout[i]] = i;
            }
            
            // mark the edges as either forward or backward relative to the layout
            vector<vector<size_t>> forward_edges(layout.size());
            vector<pair<size_t, size_t>> backward_edges;
            subgraph.for_each_edge([&](const edge_t& edge) {
                // get the indices of the edge in the layout, making sure to match
                // the canonical orientation
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
                }
                
                // always keep going
                return true;
            });
            
            // check for each node whether we've duplicated the component enough times
            // to preserve its cycles
            
            // dynamic progamming structures that represent distances within the current
            // copy of the SCC and the next copy
            vector<int64_t> distances(layout.size(), numeric_limits<int64_t>::max());
            vector<int64_t> next_distances(layout.size(), numeric_limits<int64_t>::max());
            
            // init the distances so that we are measuring from the end of the heads of
            // backward edges (which cross to the next copy of the SCC)
            for (const pair<size_t, size_t>& bwd_edge : backward_edges) {
                handle_t handle = layout[bwd_edge.first];
                distances[ordering[handle]] = -subgraph.get_length(handle);
            }
            
            // init the tracker that we use for the bail-out condition
            int64_t min_relaxed_dist = -1;
            
            // add copies until the minimum distance to the new copy is longer than the distance we're
            // trying to preserve
            for (size_t copy_num = 0; min_relaxed_dist < int64_t(min_preserved_path_length); copy_num++) {
                
#ifdef debug_dagify
                cerr << "copy number " << copy_num << endl;
#endif
                
                // do we need a new copy of this SCC to preserve paths?
                if (copy_num == injector[layout.front()].size()) {
                    // we haven't added this copy of the connected component yet
                    
#ifdef debug_dagify
                    cerr << "adding nodes for this copy" << endl;
#endif
                    
                    // add the nodes
                    for (const handle_t& original_handle : layout) {
                        // create the node in the same foward orientation as the original
                        handle_t new_handle = into->create_handle(graph->get_sequence(graph->forward(original_handle)));
                        // use the handle locally in the same orientation as it is in the layout
                        if (graph->get_is_reverse(original_handle)) {
                            new_handle = into->flip(new_handle);
                        }
#ifdef debug_dagify
                        cerr << "\t" << graph->get_id(original_handle) << " duplicated to " << into->get_id(new_handle) << endl;
#endif
                        
                        // record the translation between the graphs
                        translator[into->get_id(new_handle)] = graph->get_id(original_handle);
                        injector[original_handle].push_back(new_handle);
                    }
                    
                    // add the forward edges within this copy
                    for (size_t i = 0; i < forward_edges.size(); i++) {
                        handle_t from = injector[layout[i]].back();
                        for (const size_t& j : forward_edges[i]) {
                            into->create_edge(from, injector[layout[j]].back());
#ifdef debug_dagify
                            cerr << "\t\tfwd edge " << into->get_id(from) << (into->get_is_reverse(from) ? "-" : "+") << " -> " << into->get_id(injector[layout[j]].back()) << (into->get_is_reverse(injector[layout[j]].back()) ? "-" : "+") << endl;
#endif
                        }
                    }
                    
                    // is there a previous copy?
                    if (copy_num > 0) {
                        // add the backward edges between the copies
                        for (const pair<size_t, size_t>& bwd_edge : backward_edges) {
                            const auto& from_copies = injector[layout[bwd_edge.first]];
                            into->create_edge(from_copies[from_copies.size() - 2],
                                              injector[layout[bwd_edge.second]].back());
#ifdef debug_dagify
                            cerr << "\t\tbwd edge " << into->get_id(from_copies[from_copies.size() - 2]) << (into->get_is_reverse(from_copies[from_copies.size() - 2]) ? "-" : "+") << " -> " << into->get_id(injector[layout[bwd_edge.second]].back()) << (into->get_is_reverse(injector[layout[bwd_edge.second]].back()) ? "-" : "+") << endl;
#endif
                        }
                    }
                }
                
                // we've finished adding the copy of the SCC that corresponds to this iteration
                // now we will do the dynamic programming to bound the distance to the next SCC
                
                // find the shortest path to the nodes, staying within this copy of the SCC
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
                for (const pair<size_t, size_t>& bwd_edge : backward_edges) {
                    // skip infinity to avoid overflow
                    if (distances[bwd_edge.first] == numeric_limits<int64_t>::max()) {
                        continue;
                    }
                    
                    int64_t dist_thru = distances[bwd_edge.first] + subgraph.get_length(layout[bwd_edge.first]);
                    if (dist_thru < next_distances[bwd_edge.second]) {
                        next_distances[bwd_edge.second] = dist_thru;
                        // keep track of the shortest distance to the next copy
                        min_relaxed_dist = min(min_relaxed_dist, dist_thru);
                    }
                }
                
#ifdef debug_dagify
                cerr << "distances within component" << endl;
                for (size_t i = 0; i < distances.size(); i++) {
                    cerr << "\t" << graph->get_id(layout[i]) << (graph->get_is_reverse(layout[i]) ? "-" : "+") << " " << distances[i] << endl;
                }
                cerr << "distances to next component" << endl;
                for (size_t i = 0; i < next_distances.size(); i++) {
                    cerr << "\t" << graph->get_id(layout[i]) << (graph->get_is_reverse(layout[i]) ? "-" : "+") << " " << next_distances[i] << endl;
                }
#endif
                
                // initialize the DP structures for the next iteration
                distances = move(next_distances);
                next_distances.assign(distances.size(), numeric_limits<int64_t>::max());
            }
        }
        
#ifdef debug_dagify
        cerr << "adding edges between SCCs" << endl;
#endif
        
        // add edges between the strongly connected components
        graph->for_each_edge([&](const edge_t& canonical_edge) {
            if (component_of[graph->get_id(canonical_edge.first)] != component_of[graph->get_id(canonical_edge.second)]) {
                // this edge is between SCCs
                
                // put the edge in the order of the orientation we've imposed on the graph so
                // we can index into the look structures we created
                edge_t edge = (graph->get_is_reverse(canonical_edge.first) != reversed_nodes.count(graph->get_id(canonical_edge.first)) ?
                               edge_t(graph->flip(canonical_edge.second), graph->flip(canonical_edge.first)) :
                               canonical_edge);
                
                // connect the last copy of the first node to all copies of the second
                const handle_t& from = injector[edge.first].back();
                for (const handle_t& to : injector[edge.second]) {
                    into->create_edge(from, to);
#ifdef debug_dagify
                    cerr << "\t" << into->get_id(from) << (into->get_is_reverse(from) ? "-" : "+") << " -> " << into->get_id(to) << (into->get_is_reverse(to) ? "-" : "+") << endl;
#endif
                }
            }
            
            // always keep going
            return true;
        });
        
        // return the ID translator
        return translator;
    }
}
}

