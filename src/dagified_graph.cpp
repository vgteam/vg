/**
 * \file dagified_graph.cpp: contains the implementation of DagifiedGraph
 */


#include "dagified_graph.hpp"

//#define debug_dagify

namespace vg {

using namespace std;

    DagifiedGraph::DagifiedGraph(const HandleGraph* graph, size_t min_preserved_path_length) : graph(graph) {
        
#ifdef debug_dagify
        cerr << "constructing dagified graph" << endl;
#endif
        
        // find the numeric range of handles in the underlying graph (needed for later bookkeeping)
        uint64_t max_handle = std::numeric_limits<uint64_t>::min();
        graph->for_each_handle([&](const handle_t& handle) {
            for (handle_t h : {handle, graph->flip(handle)}) {
                uint64_t integer_handle = handlegraph::as_integer(h);
                min_handle = min(integer_handle, min_handle);
                max_handle = max(integer_handle, max_handle);
            }
        });
        handle_val_range = max_handle - min_handle + 1;
        
#ifdef debug_dagify
        cerr << "graph has handle range " << handle_val_range << ", min handle " << min_handle << ", min ID " << graph->min_node_id() << ", and ID range " << (graph->max_node_id() - graph->min_node_id() + 1) << endl;
        cerr << "preserving walks up to length " << min_preserved_path_length << endl;
#endif
        
        // now we begin the dagify algorithm
        
        vector<vector<handle_t>> strong_components;
        {
            // get a low-FAS layout with a canonical orientation for each handle
            vector<handle_t> layout = handlealgs::eades_algorithm(graph);
            
            // invert the mapping for the layout
            layout_order.reserve(layout.size());
            for (size_t i = 0; i < layout.size(); ++i) {
                layout_order[layout[i]] = i;
            }
            
            
            // identify the SCCs and build the reverse SCC mapping
            // TODO: annoying that we have to work with this return type for strongly_connected_components
            scc_of_handle.resize(layout.size());
            size_t scc_idx = 0;
            for (const unordered_set<id_t>& scc : handlealgs::strongly_connected_components(graph)) {
                // init new component
                strong_components.emplace_back();
                auto& component = strong_components.back();
                component.reserve(scc.size());
                
                // build the reverse mapping for this SCC
                for (const id_t& node_id : scc) {
                    handle_t handle = graph->get_handle(node_id);
                    auto iter = layout_order.find(handle);
                    if (iter == layout_order.end()) {
                        iter = layout_order.find(graph->flip(handle));
                    }
                    scc_of_handle[iter->second] = scc_idx;
                }
                ++scc_idx;
            }
            
            // build the SCC components in layout order
            for (const handle_t& handle : layout) {
                strong_components[scc_of_handle[layout_order[handle]]].push_back(handle);
            }
            
            // let the layout fall out of scope
        }
#ifdef debug_dagify
        cerr << "identified " << strong_components.size() << " strongly connected components out of " << graph->get_node_count() << " nodes" << endl;
#endif
        
        // identify how many times each SCC needs to be duplicated
        scc_copy_count.resize(strong_components.size());
        for (size_t scc_idx = 0; scc_idx < strong_components.size(); ++scc_idx) {
#ifdef debug_dagify
            cerr << "BEGIN NEW SCC " << scc_idx << endl;
#endif
            const vector<handle_t>& component = strong_components[scc_idx];
            
            // record the ordering of the layout so we can build adjacency lists
            unordered_map<handle_t, size_t> ordering;
            for (size_t i = 0; i < component.size(); i++) {
                ordering[component[i]] = i;
            }
            
            // mark the edges as either forward or backward relative to the layout
            vector<vector<size_t>> forward_edges(component.size());
            vector<pair<size_t, size_t>> backward_edges;
            for (size_t i = 0; i < component.size(); ++i) {
                graph->follow_edges(component[i], false, [&](const handle_t& next) {
                    if (scc_of_handle[layout_order[next]] == scc_idx) {
                        // this edge is internal to a strongly connected component
                        size_t j = ordering[next];
                        if (i < j) {
                            // non feedback arc
                            forward_edges[i].push_back(j);
                        }
                        else {
                            // feedback arc
                            backward_edges.emplace_back(i, j);
                        }
                    }
                });
            }
#ifdef debug_dagify
            cerr << "feedforward graph:" << endl;
            for (size_t i = 0; i < component.size(); ++i) {
                cerr << graph->get_id(component[i]) << ":";
                for (auto j : forward_edges[i]) {
                    cerr << " " << graph->get_id(component[j]);
                }
                cerr << endl;
            }
            cerr << "feedback edges:" << endl;
            for (auto edge : backward_edges) {
                cerr << graph->get_id(component[edge.first]) << " -> " << graph->get_id(component[edge.second]) << endl;
            }
#endif
            
            // check for each node whether we've duplicated the component enough times
            // to preserve its cycles
            
            // dynamic progamming structures that represent distances within the current
            // copy of the SCC and the next copy
            vector<int64_t> distances(component.size(), numeric_limits<int64_t>::max());
            
            // init the distances so that we are measuring from the end of the heads of
            // backward edges (which cross to the next copy of the SCC)
            for (const pair<size_t, size_t>& bwd_edge : backward_edges) {
                distances[bwd_edge.first] = -graph->get_length(component[bwd_edge.first]);
            }
            
            // init the tracker that we use for the bail-out condition
            int64_t min_relaxed_dist = -1;
            
            // keep track of how many times we've implicitly copied
            uint64_t copy_num = 0;
            for (; min_relaxed_dist < int64_t(min_preserved_path_length); copy_num++) {
                
#ifdef debug_dagify
                cerr << "making " << copy_num << "-th copy of SCC with incoming min relaxed distance " << min_relaxed_dist << endl;
#endif
                
                // the distances in the next copy unit
                vector<int64_t> next_distances(component.size(), numeric_limits<int64_t>::max());
                
                // find the shortest path to the nodes, staying within this copy of the SCC
                for (size_t i = 0; i < distances.size(); i++) {
                    // skip infinity to avoid overflow
                    if (distances[i] == numeric_limits<int64_t>::max()) {
                        continue;
                    }
                    
                    int64_t dist_thru = distances[i] + graph->get_length(component[i]);
                    for (const size_t& j : forward_edges[i]) {
                        distances[j] = min(distances[j], dist_thru);
                    }
                }
                
                // now find the minimum distance to nodes in the next copy of the SCC
                min_relaxed_dist = numeric_limits<int64_t>::max();
                for (const pair<size_t, size_t>& bwd_edge : backward_edges) {
                    // skip infinity to avoid overflow
                    if (distances[bwd_edge.first] == numeric_limits<int64_t>::max()) {
                        continue;
                    }
                    
                    int64_t dist_thru = distances[bwd_edge.first] +  graph->get_length(component[bwd_edge.first]);
                    if (dist_thru < next_distances[bwd_edge.second]) {
                        next_distances[bwd_edge.second] = dist_thru;
                        // keep track of the shortest distance to the next copy
                        min_relaxed_dist = min(min_relaxed_dist, dist_thru);
                    }
                }
                
#ifdef debug_dagify
                cerr << "distances within component" << endl;
                for (size_t i = 0; i < distances.size(); i++) {
                    cerr << "\t" << graph->get_id(component[i]) << (graph->get_is_reverse(component[i]) ? "-" : "+") << " ";
                    if (distances[i] != numeric_limits<int64_t>::max()) {
                        cerr << distances[i];
                    }
                    else {
                        cerr << ".";
                    }
                    cerr << endl;
                }
                cerr << "distances to next component" << endl;
                for (size_t i = 0; i < next_distances.size(); i++) {
                    cerr << "\t" << graph->get_id(component[i]) << (graph->get_is_reverse(component[i]) ? "-" : "+") << " ";
                    if (distances[i] != numeric_limits<int64_t>::max()) {
                        cerr << distances[i];
                    }
                    else {
                        cerr << ".";
                    }
                    cerr << endl;
                }
#endif
                
                // initialize the DP structures for the next iteration
                distances = move(next_distances);
            }
            
            // now we know that the copy count needs to be, so we can record the information we need
            // from this component
            
            // record the copy count
            scc_copy_count[scc_idx] = copy_num;
            
            // add the number of nodes to the total count
            node_count += component.size() * copy_num;
            
            // find the maximum projected node ID in this component
            id_t comp_max_id = numeric_limits<id_t>::min();
            for (const handle_t& handle : component) {
                comp_max_id = max(comp_max_id, graph->get_id(handle));
            }
            max_id = max<id_t>(max_id, comp_max_id + (copy_num - 1) * (graph->max_node_id() - graph->min_node_id() + 1));
        }
    }
    
    bool DagifiedGraph::has_node(id_t node_id) const {
        id_t original_id = get_underlying_id(node_id);
        bool exists = graph->has_node(original_id);
        if (exists) {
            // was the node duplicated enough times to have created this node ID?
            exists = scc_copy_of_node_id(node_id) < scc_copy_count.at(scc_of_handle.at(layout_order_of_handle(graph->get_handle(original_id))));
        }
        return exists;
    }
    
    handle_t DagifiedGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return nth_copy_of_handle(graph->get_handle(get_underlying_id(node_id), is_reverse), scc_copy_of_node_id(node_id));
    }
    
    id_t DagifiedGraph::get_id(const handle_t& handle) const {
        return (graph->get_id(get_underlying_handle(handle))
                + scc_copy_of_handle(handle) * (graph->max_node_id() - graph->min_node_id() + 1));
    }
    
    bool DagifiedGraph::get_is_reverse(const handle_t& handle) const {
        return graph->get_is_reverse(get_underlying_handle(handle));
    }
    
    handle_t DagifiedGraph::flip(const handle_t& handle) const {
        return nth_copy_of_handle(graph->flip(get_underlying_handle(handle)), scc_copy_of_handle(handle));
    }
    
    size_t DagifiedGraph::get_length(const handle_t& handle) const {
        return graph->get_length(get_underlying_handle(handle));
    }
    
    string DagifiedGraph::get_sequence(const handle_t& handle) const {
        return graph->get_sequence(get_underlying_handle(handle));
    }
    
    bool DagifiedGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                          const function<bool(const handle_t&)>& iteratee) const {
        
        // this is the complicated part where we have to induce an edge structure that is a DAG
        handle_t underlying = get_underlying_handle(handle);
        uint64_t scc_copy = scc_copy_of_handle(handle);
        
        // does the handle match the canonical orientation of the layout we computed?
        bool matches_layout = layout_order.count(underlying);
        
        // are we traversing this node forward along the canonical orientation?
        bool canonical_fwd = matches_layout != go_left;
        
        size_t layout_index = layout_order_of_handle(underlying);
        uint64_t scc_id = scc_of_handle.at(layout_index);
        
        return graph->follow_edges(underlying, go_left, [&](const handle_t& next_underlying) {
            
            size_t next_layout_index = layout_order_of_handle(next_underlying);
            uint64_t next_scc_id = scc_of_handle.at(next_layout_index);
            
            bool keep_going = true;
            if (next_scc_id != scc_id) {
                // this is over an edge that's between two strongly connected component
                uint64_t next_scc_count = scc_copy_count.at(next_scc_id);
                if (canonical_fwd) {
                    // we are traveling in the canonically forward direction
                    if (scc_copy + 1 == scc_copy_count.at(scc_id)) {
                        // only the last copy of a handle is allowed to extend to the next strongly
                        // connected component, and it connects to all copies in the next one
                        for (size_t i = 0; i < next_scc_count && keep_going; ++i) {
                            keep_going = iteratee(nth_copy_of_handle(next_underlying, i));
                        }
                    }
                }
                else {
                    // we are going in the reverse direction of the canonical orientation, so we
                    // can only connect to the last copy of the next handle
                    keep_going = iteratee(nth_copy_of_handle(next_underlying, next_scc_count - 1));
                }
            }
            else {
                // this edge is internal to a strongly connected component
                if (canonical_fwd) {
                    // we are traversing in the direction of the canonical orientation
                    if (next_layout_index > layout_index) {
                        // we are not taking a reversing edge, so we stay in the same copy unit
                        keep_going = iteratee(nth_copy_of_handle(next_underlying, scc_copy));
                    }
                    else if (scc_copy + 1 < scc_copy_count.at(scc_id)) {
                        // we are taking a reversing edge, and there are still copy units ahead
                        keep_going = iteratee(nth_copy_of_handle(next_underlying, scc_copy + 1));
                    }
                }
                else {
                    // we are traversing against the direction of the canonical orientation
                    if (next_layout_index < layout_index) {
                        // we are moving backwards over a forward edge, stay in the same copy unit
                        keep_going = iteratee(nth_copy_of_handle(next_underlying, scc_copy));
                    }
                    else if (scc_copy > 0) {
                        // we are moving backwards over a reversing edge, and there are still copy units before this one
                        keep_going = iteratee(nth_copy_of_handle(next_underlying, scc_copy - 1));
                    }
                }
            }
            
            return keep_going;
        });
    }
    
    bool DagifiedGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                             bool parallel) const {
        return graph->for_each_handle([&](const handle_t& underlying) {
            // iterate over however many copies of the handle there are
            size_t copy_count = scc_copy_count.at(scc_of_handle.at(layout_order_of_handle(underlying)));
            bool keep_going = true;
            for (size_t i = 0; i < copy_count && keep_going; ++i) {
                keep_going = iteratee(nth_copy_of_handle(underlying, i));
            }
            return keep_going;
        }, parallel);
    }
    
    size_t DagifiedGraph::get_node_count() const {
        return node_count;
    }
    
    id_t DagifiedGraph::min_node_id() const {
        // duplicated handles only increase in ID, so the original minimum doesn't change
        return graph->min_node_id();
    }
    
    id_t DagifiedGraph::max_node_id() const {
        return max_id;
    }

    handle_t DagifiedGraph::get_underlying_handle(const handle_t& handle) const {
        return handlegraph::as_handle(((uint64_t(handlegraph::as_integer(handle)) - min_handle) % handle_val_range) + min_handle);
    }

    uint64_t DagifiedGraph::scc_copy_of_handle(const handle_t& handle) const {
        return (uint64_t(handlegraph::as_integer(handle)) - min_handle) / handle_val_range;
    }

    uint64_t DagifiedGraph::scc_copy_of_node_id(const id_t& node_id) const {
        return (node_id - graph->min_node_id()) / (graph->max_node_id() - graph->min_node_id() + 1);
    }

    size_t DagifiedGraph::layout_order_of_handle(const handle_t& handle) const {
        auto iter = layout_order.find(handle);
        if (iter == layout_order.end()) {
            iter = layout_order.find(graph->flip(handle));
        }
        return iter->second;
    }


    id_t DagifiedGraph::get_underlying_id(const id_t& node_id) const {
        return ((node_id - graph->min_node_id()) % (graph->max_node_id() - graph->min_node_id() + 1)) + graph->min_node_id();
    }

    handle_t DagifiedGraph::nth_copy_of_handle(const handle_t& handle, const uint64_t& n) const {
        return handlegraph::as_handle(uint64_t(handlegraph::as_integer(handle)) + n * handle_val_range);
    }
}

