/**
 * \file extract_containing_graph.cpp
 *
 * Implementation for the extract_containing_graph algorithm.
 */
 
#include "extract_containing_graph.hpp"

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {

using namespace structures;

void extract_containing_graph(const HandleGraph* source,
                              MutableHandleGraph* into,
                              const vector<pos_t>& positions,
                              const vector<size_t>& forward_search_lengths,
                              const vector<size_t>& backward_search_lengths,
                              size_t reversing_walk_length) {
    
    if (forward_search_lengths.size() != backward_search_lengths.size()
        || forward_search_lengths.size() != positions.size()) {
        cerr << "error:[extract_containing_graph] subgraph extraction search lengths do not match seed positions" << endl;
        exit(1);
    }
    
    if (into->get_node_count()) {
        cerr << "error:[extract_containing_graph] must extract into an empty graph" << endl;
        exit(1);
    }
    
#ifdef debug_vg_algorithms
    cerr << "[extract_containing_graph] extracting containing graph from the following points:" << endl;
    for (size_t i = 0; i < positions.size(); i ++) {
        cerr << "\t" << positions[i] << ", forward dist " << forward_search_lengths[i] << ", backward dist " << backward_search_lengths[i] << endl;
    }
#endif
    
    
    // computing search distances relative to this maximum allows us to keep the searches
    // from all of the seed nodes in the same priority queue so that we only need to do
    // one Dijkstra traversal
    int64_t max_search_length = max(*std::max_element(forward_search_lengths.begin(), forward_search_lengths.end()),
                                    *std::max_element(backward_search_lengths.begin(), backward_search_lengths.end()));
    
    
#ifdef debug_vg_algorithms
    cerr << "[extract_containing_graph] artificial offset calculated to be " << max_search_length << endl;
#endif
    // for keeping track of all the edges we cross to add later
    //
    // we use spp because the order of the edges affects some tie-breaking behavior
    // later on that we want to remain system independent (i.e. no dependence on the
    // system stdlib)
    spp::sparse_hash_set<edge_t> observed_edges;
    
    // initialize the queue, opposite order so priority queue selects minimum
    // priority represent distance from starting pos to the left side of this node
    RankPairingHeap<handle_t, int64_t, greater<int64_t>> queue;
    
    for (size_t i = 0; i < positions.size(); i++) {
        
        const pos_t& pos = positions[i];
        handle_t source_handle = source->get_handle(id(pos), false);
        
        // add all of the initial nodes to the graph
        if (!into->has_node(id(pos))) {
            into->create_handle(source->get_sequence(source_handle), id(pos));
        }
        
        // compute the modified search lengths
        int64_t dist_forward = -offset(pos) + max_search_length - forward_search_lengths[i];
        int64_t dist_backward = offset(pos) - source->get_length(source_handle) + max_search_length - backward_search_lengths[i];
        
        // add a traversal for each direction
        queue.push_or_reprioritize(is_rev(pos) ? source->flip(source_handle) : source_handle, dist_forward);
        queue.push_or_reprioritize(is_rev(pos) ? source_handle : source->flip(source_handle), dist_backward);
#ifdef debug_vg_algorithms
        cerr << "[extract_containing_graph] init enqueue " << id(pos) << " " << is_rev(pos) << ": " << dist_forward << endl;
        cerr << "[extract_containing_graph] init enqueue " << id(pos) << " " << !is_rev(pos) << ": " << dist_backward << endl;
#endif
    }
    
    while (!queue.empty()) {
        // get the next shortest distance traversal from either the init
        pair<handle_t, int64_t> trav = queue.top();
        queue.pop();
        
        
#ifdef debug_vg_algorithms
        cerr << "[extract_containing_graph] dequeue " << source->get_id(trav.first) << " " << source->get_is_reverse(trav.first) << ": " << trav.second << endl;
#endif
        
        // make sure the node is in the graph
        if (!into->has_node(source->get_id(trav.first))) {
#ifdef debug_vg_algorithms
            cerr << "[extract_containing_graph] adding node to subgraph" << endl;
#endif
            into->create_handle(source->get_sequence(source->forward(trav.first)), source->get_id(trav.first));
        }
        
        int64_t dist_thru = trav.second + source->get_length(trav.first);
        
        if (dist_thru < max_search_length) {
            // we can add more nodes along same path without going over the max length
            
            // look locally right from this position
            source->follow_edges(trav.first, false, [&](const handle_t& next) {
                // record the edge
                observed_edges.insert(source->edge_handle(trav.first, next));
                // add it to the queue
                queue.push_or_reprioritize(next, dist_thru);
#ifdef debug_vg_algorithms
                cerr << "[extract_containing_graph] traverse and (possibly) enqueue " << source->get_id(next) << " " << source->get_is_reverse(next) << ": " << dist_thru << endl;
#endif
            });
        }
        
        // if we're allowing reversing walks and this isn't one of the starting positions...
        if (reversing_walk_length > 0 && trav.second > 0) {
            
            // choose a distance that will let the reversing walk only as far as the minimum of the max
            // search length or the reverse length
            int64_t synthetic_dist = max_search_length > reversing_walk_length ? max_search_length - reversing_walk_length : 0;
            
            handle_t flipped = source->flip(trav.first);
            // look locally right from this position
            source->follow_edges(flipped, false, [&](const handle_t& next) {
                // record the edge
                observed_edges.insert(source->edge_handle(flipped, next));
                
                // add it to the queue
                queue.push_or_reprioritize(next, max_search_length - reversing_walk_length);
#ifdef debug_vg_algorithms
                cerr << "[extract_containing_graph] reverse walk and (possibly) enqueue " << source->get_id(next) << " " << source->get_is_reverse(next) << ": " << max_search_length - reversing_walk_length << endl;
#endif
            });
        }
    }
    
    // add the edges to the graph
    for (const edge_t& edge : observed_edges) {
#ifdef debug_vg_algorithms
        cerr << "[extract_containing_graph] adding edge " << source->get_id(edge.first) << " " << source->get_is_reverse(edge.first) << " -> " << source->get_id(edge.second) << " " << source->get_is_reverse(edge.second) << endl;
#endif
        into->create_edge(into->get_handle(source->get_id(edge.first), source->get_is_reverse(edge.first)),
                          into->get_handle(source->get_id(edge.second), source->get_is_reverse(edge.second)));
    }
}

void extract_containing_graph(const HandleGraph* source, MutableHandleGraph* into, const vector<pos_t>& positions,
                              size_t max_dist, size_t reversing_walk_length) {
    
    extract_containing_graph(source, into, positions, vector<size_t>(positions.size(), max_dist),
                             reversing_walk_length);
}

void extract_containing_graph(const HandleGraph* source, MutableHandleGraph* into, const vector<pos_t>& positions,
                              const vector<size_t>& position_max_dist, size_t reversing_walk_length) {
    
    extract_containing_graph(source, into, positions, position_max_dist, position_max_dist,
                             reversing_walk_length);
}

}
}
