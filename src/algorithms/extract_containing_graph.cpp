/**
 * \file extract_containing_graph.cpp
 *
 * Implementation for the extract_containing_graph algorithm.
 */
 
#include "extract_containing_graph.hpp"

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {

    void extract_containing_graph_internal(Graph& g, const vector<pos_t>& positions,
                                           const vector<size_t>& forward_search_lengths,
                                           const vector<size_t>& backward_search_lengths,
                                           function<vector<Edge>(id_t)> edges_on_start,
                                           function<vector<Edge>(id_t)> edges_on_end,
                                           function<string(id_t)> node_sequence) {
        
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
        
        // a local struct that packages a node traversal with its distance from the first position
        // note: we don't use NodeTraversals because we may be using an XG, in which case we don't
        // have Node*'s
        struct Traversal {
            Traversal(id_t id, bool rev, int64_t dist) : id(id), rev(rev), dist(dist) {}
            int64_t dist; // distance to the right side of this node
            id_t id;      // node ID
            bool rev;     // strand
            inline bool operator<(const Traversal& other) const {
                return dist > other.dist; // opposite order so priority queue selects minimum
            }
        };
        
        // represent an edge in a canonical ((from, from_start), (to, to_end)) format so that a hash can
        // identify edges as the same or different
        auto canonical_form_edge = [](const id_t& id, const bool& left_side, const id_t& to_id, const bool& reversing) {
            if (id < to_id) {
                return make_pair(make_pair(id, left_side), make_pair(to_id, left_side != reversing));
            }
            else if (to_id < id) {
                return make_pair(make_pair(to_id, left_side == reversing), make_pair(id, !left_side));
            }
            else {
                return make_pair(make_pair(id, reversing && left_side), make_pair(id, reversing && !left_side));
            }
        };
        
        // a graph index that we will maintain as we extract the subgraph
        unordered_map<id_t, Node*> graph;
        
        size_t max_search_length = max(*std::max_element(forward_search_lengths.begin(), forward_search_lengths.end()),
                                       *std::max_element(backward_search_lengths.begin(), backward_search_lengths.end()));
        
        
        // initialize the queue
        priority_queue<Traversal> queue;
        for (size_t i = 0; i < positions.size(); i++) {
            const pos_t& pos = positions[i];
            // add all of the initial nodes to the graph
            if (!graph.count(id(pos))) {
                Node* node = g.add_node();
                node->set_sequence(node_sequence(id(pos)));
                node->set_id(id(pos));
                graph[id(pos)] = node;
            }
            
            // adding this extra distance allows us to keep the searches from all of the seed nodes in
            // the same priority queue so that we only need to do one Dijkstra traversal
            
            // add a traversal for each direction
            size_t dist_forward = graph[id(pos)]->sequence().size() - offset(pos) + max_search_length - forward_search_lengths[i];
            size_t dist_backward = offset(pos) + max_search_length - backward_search_lengths[i];
            if (dist_forward < max_search_length) {
                queue.emplace(id(pos), is_rev(pos), dist_forward);
            }
            if (dist_backward < max_search_length) {
                queue.emplace(id(pos), !is_rev(pos), dist_backward);
            }
        }
        
        unordered_set<pair<id_t, bool>> traversed;
        unordered_set<pair<pair<id_t, bool>, pair<id_t, bool>>> observed_edges;
        
        while (!queue.empty()) {
            // get the next shortest distance traversal from either the init
            Traversal trav = queue.top();
            queue.pop();
            
            // make sure we haven't traversed this node already
            if (traversed.count(make_pair(trav.id, trav.rev))) {
                continue;
            }
            // mark the node as traversed
            traversed.emplace(trav.id, trav.rev);
            
            // which side are we traversing out of?
            for (Edge& edge : (trav.rev ? edges_on_start(trav.id) : edges_on_end(trav.id))) {
                
                // get the orientation and id of the other side of the edge
                id_t next_id;
                bool next_rev;
                if (edge.from() == trav.id && edge.from_start() == trav.rev) {
                    next_id = edge.to();
                    next_rev = edge.to_end();
                }
                else {
                    next_id = edge.from();
                    next_rev = !edge.from_start();
                }
                
                // record the edge
                observed_edges.insert(canonical_form_edge(trav.id, trav.rev, next_id, next_rev != trav.rev));
                
                // make sure the node is in the graph
                if (!graph.count(next_id)) {
                    Node* node = g.add_node();
                    node->set_sequence(node_sequence(next_id));
                    node->set_id(next_id);
                    graph[next_id] = node;
                }
                
                // distance to the end of this node
                int64_t dist_thru = trav.dist + graph[next_id]->sequence().size();
                if (!traversed.count(make_pair(next_id, next_rev)) && dist_thru < max_search_length) {
                    // we can add more nodes along same path without going over the max length
                    queue.emplace(next_id, next_rev, dist_thru);
                }
            }
        }
        
        // add the edges to the graph
        for (const pair<pair<id_t, bool>, pair<id_t, bool>>& edge : observed_edges) {
            Edge* e = g.add_edge();
            e->set_from(edge.first.first);
            e->set_to(edge.second.first);
            e->set_from_start(edge.first.second);
            e->set_to_end(edge.second.second);
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
    
    void extract_containing_graph(const HandleGraph* source, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_forward_max_dist,
                                  const vector<size_t>& position_backward_max_dist) {
        
        return extract_containing_graph_internal(g, positions, position_forward_max_dist, position_backward_max_dist,
                                                 [&](id_t id) {
                                                    // Get edges on start
                                                    vector<Edge> to_return;
                                                    auto here = source->get_handle(id, false);
                                                    source->follow_edges(here, true, [&](const handle_t& left) {
                                                        to_return.push_back(xg::make_edge(source->get_id(left),
                                                            source->get_is_reverse(left), id, false));
                                                    });
                                                    return to_return;
                                                },
                                                [&](id_t id) {
                                                    // Get edges on end
                                                    vector<Edge> to_return;
                                                    auto here = source->get_handle(id, false);
                                                    source->follow_edges(here, false, [&](const handle_t& right) {
                                                        to_return.push_back(xg::make_edge(id, false, source->get_id(right),
                                                            source->get_is_reverse(right)));
                                                    });
                                                    return to_return;
                                                },
                                                [&](id_t id) {
                                                    // Get sequence
                                                    return source->get_sequence(source->get_handle(id, false));
                                                });
    
    } 

}
}
