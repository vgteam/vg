/**
 * \file extract_containing_graph.cpp
 *
 * Implementation for the extract_containing_graph algorithm.
 */
 
#include "extract_containing_graph.hpp"
#include <structures/updateable_priority_queue.hpp>

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {

using namespace structures;

void extract_containing_graph(const HandleGraph* source, Graph& g, const vector<pos_t>& positions,
                              const vector<size_t>& forward_search_lengths,
                              const vector<size_t>& backward_search_lengths) {
    
    if (forward_search_lengths.size() != backward_search_lengths.size()
        || forward_search_lengths.size() != positions.size()) {
        cerr << "error:[extract_containing_graph] subgraph extraction search lengths do not match seed positions" << endl;
        assert(false);
    }
    
    if (g.node_size() || g.edge_size()) {
        cerr << "error:[extract_containing_graph] must extract into an empty graph" << endl;
        assert(false);
    }
    
#ifdef debug_vg_algorithms
    cerr << "[extract_containing_graph] extracting containing graph from the following points:" << endl;
    for (size_t i = 0; i < positions.size(); i ++) {
        cerr << "\t" << positions[i] << ", forward dist " << forward_search_lengths[i] << ", backward dist " << backward_search_lengths[i] << endl;
    }
#endif
    
    // TODO: these structs are duplicative with extract_connecting_graph
    
    // a local struct that packages a handle with its distance from the first position
    struct Traversal {
        Traversal(handle_t handle, int64_t dist) : handle(handle), dist(dist) {}
        int64_t dist; // distance from pos to the right side of this node
        handle_t handle; // Oriented node traversal
        inline bool operator<(const Traversal& other) const {
            return dist > other.dist; // opposite order so priority queue selects minimum
        }
    };
    
    // a graph index that we will maintain as we extract the subgraph
    unordered_map<id_t, Node*> graph;
    
    size_t max_search_length = max(*std::max_element(forward_search_lengths.begin(), forward_search_lengths.end()),
                                   *std::max_element(backward_search_lengths.begin(), backward_search_lengths.end()));
    
    
    // initialize the queue
    UpdateablePriorityQueue<Traversal, handle_t> queue([](const Traversal& item) {
        return item.handle;
    });
    for (size_t i = 0; i < positions.size(); i++) {
        const pos_t& pos = positions[i];
        // add all of the initial nodes to the graph
        if (!graph.count(id(pos))) {
            Node* node = g.add_node();
            // TODO: this might require more get_handle calls than we want
            auto handle = source->get_handle(id(pos), false);
            node->set_sequence(source->get_sequence(handle));
            node->set_id(id(pos));
            graph[id(pos)] = node;
        }
        
        // adding this extra distance allows us to keep the searches from all of the seed nodes in
        // the same priority queue so that we only need to do one Dijkstra traversal
        
        // add a traversal for each direction
        size_t dist_forward = graph[id(pos)]->sequence().size() - offset(pos) + max_search_length - forward_search_lengths[i];
        size_t dist_backward = offset(pos) + max_search_length - backward_search_lengths[i];
        if (dist_forward < max_search_length) {
            queue.emplace(source->get_handle(id(pos), is_rev(pos)), dist_forward);
        }
        if (dist_backward < max_search_length) {
            queue.emplace(source->get_handle(id(pos), !is_rev(pos)), dist_backward);
        }
    }
    
    unordered_set<pair<handle_t, handle_t>> observed_edges;
    
    while (!queue.empty()) {
        // get the next shortest distance traversal from either the init
        Traversal trav = queue.top();
        queue.pop();
        
        source->follow_edges(trav.handle, false, [&](const handle_t& next) {
            // Look locally right from this position
            
            // Get the ID of where we're going.
            auto next_id = source->get_id(next);
            
            // record the edge
            observed_edges.insert(source->edge_handle(trav.handle, next));
            
            // make sure the node is in the graph
            if (!graph.count(next_id)) {
                Node* node = g.add_node();
                node->set_sequence(source->get_sequence(source->forward(next)));
                node->set_id(next_id);
                graph[next_id] = node;
            }
            
            // distance to the end of this node
            int64_t dist_thru = trav.dist + graph[next_id]->sequence().size();
            if (dist_thru < max_search_length) {
                // we can add more nodes along same path without going over the max length
                queue.emplace(next, dist_thru);
            }
        });
    }
    
    // add the edges to the graph
    for (const pair<handle_t, handle_t>& edge : observed_edges) {
        Edge* e = g.add_edge();
        e->set_from(source->get_id(edge.first));
        e->set_from_start(source->get_is_reverse(edge.first));
        e->set_to(source->get_id(edge.second));
        e->set_to_end(source->get_is_reverse(edge.second));
    }
}

void extract_containing_graph(const HandleGraph* source, Graph& g, const vector<pos_t>& positions, size_t max_dist) {
    
    // make a dummy vector for all positions at the same distance
    vector<size_t> dists(positions.size(), max_dist);
    return extract_containing_graph(source, g, positions, dists);
}

void extract_containing_graph(const HandleGraph* source, Graph& g, const vector<pos_t>& positions,
                              const vector<size_t>& position_max_dist) {
    
    return extract_containing_graph(source, g, positions, position_max_dist, position_max_dist);
}

}
}
