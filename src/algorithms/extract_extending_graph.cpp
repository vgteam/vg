/**
 * \file extract_extending_graph.cpp
 *
 * Implementation for the extract_extending_graph algorithm.
 */
 
#include "extract_extending_graph.hpp"

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {

    unordered_map<id_t, id_t> extract_extending_graph_internal(Graph& g, int64_t max_dist, pos_t pos, bool backward,
                                                               bool preserve_cycles_on_src,
                                                               function<vector<Edge>(id_t)> edges_on_start,
                                                               function<vector<Edge>(id_t)> edges_on_end,
                                                               function<string(id_t)> node_sequence) {
        
        if (g.node_size() || g.edge_size()) {
            cerr << "error:[extract_extending_graph] must extract into an empty graph" << endl;
            assert(false);
        }
        
#ifdef debug_vg_algorithms
        cerr << "[extract_extending_graph] extracting exending graph from " << pos << " in " << (backward ? "backward" : "forward") << " direction with max search dist " << max_dist << endl;
#endif
        
        // TODO: these structs are duplicative with extract_connecting_graph
        
        // a local struct that packages a node traversal with its distance from the first position
        // note: we don't use NodeTraversals because we may be using an XG, in which case we don't
        // have Node*'s
        struct Traversal {
            Traversal(id_t id, bool rev, int64_t dist) : id(id), rev(rev), dist(dist) {}
            int64_t dist; // distance from pos to the right side of this node
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
        
        // a map from node ids in the extracted graph to the node ids in the original graph
        unordered_map<id_t, id_t> id_trans;
        
        // a graph index that we will maintain as we extract the subgraph
        unordered_map<id_t, Node*> graph;
        Node* src_node = g.add_node();
        src_node->set_sequence(node_sequence(id(pos)));
        src_node->set_id(id(pos));
        graph[id(pos)] = src_node;
        id_trans[id(pos)] = id(pos);
        
        // initialize the queue
        priority_queue<Traversal> queue;
        if (backward) {
            int64_t dist = offset(pos);
            if (dist < max_dist) {
                queue.emplace(id(pos), !is_rev(pos), dist);
            }
        }
        else {
            int64_t dist = g.node(0).sequence().size() - offset(pos);
            if (dist < max_dist) {
                queue.emplace(id(pos), is_rev(pos), dist);
            }
        }
        
        id_t max_id = id(pos);
        bool cycled_to_source = false;
        unordered_set<pair<id_t, bool>> traversed;
        unordered_set<pair<pair<id_t, bool>, pair<id_t, bool>>> observed_edges;
        
        while (!queue.empty()) {
            // get the next shortest distance traversal from either the init
            Traversal trav = queue.top();
            queue.pop();
            
#ifdef debug_vg_algorithms
            cerr << "[extract_extending_graph] traversing " << trav.id << (trav.rev ? "-" : "+") << " at dist " << trav.dist << endl;
#endif
            
            // make sure we haven't traversed this node already
            if (traversed.count(make_pair(trav.id, trav.rev))) {
                continue;
            }
            // mark the node as traversed
            traversed.emplace(trav.id, trav.rev);
            
            // which side are we traversing out of?
            for (Edge& edge : (trav.rev ? edges_on_start(trav.id) : edges_on_end(trav.id))) {
                
#ifdef debug_vg_algorithms
                cerr << "[extract_extending_graph] checking edge " << pb2json(edge) << endl;
#endif
                
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
                
                // do we ever return to the source node?
                cycled_to_source = cycled_to_source || next_id == id(pos);
                
                // record the edge
                observed_edges.insert(canonical_form_edge(trav.id, trav.rev, next_id, next_rev != trav.rev));
                
                // make sure the node is in the graph
                if (!graph.count(next_id)) {
                    Node* node = g.add_node();
                    node->set_sequence(node_sequence(next_id));
                    node->set_id(next_id);
                    graph[next_id] = node;
                    id_trans[next_id] = next_id;
                    max_id = std::max(next_id, max_id);
                }
                
                // distance to the end of this node
                int64_t dist_thru = trav.dist + graph[next_id]->sequence().size();
                if (!traversed.count(make_pair(next_id, next_rev)) && dist_thru < max_dist) {
                    // we can add more nodes along same path without going over the max length
                    queue.emplace(next_id, next_rev, dist_thru);
#ifdef debug_vg_algorithms
                    cerr << "[extract_extending_graph] enqueuing " << next_id << (next_rev ? "-" : "+") << " at dist " << dist_thru << endl;
#endif
                }
            }
        }
        
        vector<pair<pair<id_t, bool>, pair<id_t, bool>>> src_edges;
        
        // add the edges to the graph
        for (const pair<pair<id_t, bool>, pair<id_t, bool>>& edge : observed_edges) {
            bool add_edge;
            if (edge.first.first == id(pos) || edge.second.first == id(pos)) {
                // record the edges that touch the source in case we need to duplicate them
                src_edges.push_back(edge);
                // does the edge only touch the side of the source node that's not going to be cut?
                add_edge = !(edge.first.first == id(pos) && edge.first.second != (backward != is_rev(pos))) &&
                           !(edge.second.first == id(pos) && edge.second.second == (backward != is_rev(pos)));
#ifdef debug_vg_algorithms
                cerr << "[extract_extending_graph] " << (add_edge ? "" : "not ") << "adding edge " << edge.first.first << (edge.first.second ? "-" : "+") << " -> " << edge.second.first << (edge.second.second ? "-" : "+") << " that touches source node" << endl;
#endif
            }
            else {
                // the edge doesn't touch the source node, so it will certainly survive the cut
                add_edge = true;
#ifdef debug_vg_algorithms
                cerr << "[extract_extending_graph] " << "adding edge " << edge.first.first << (edge.first.second ? "-" : "+") << " -> " << edge.second.first << (edge.second.second ? "-" : "+") << " that does not touch source node" << endl;
#endif
            }
            
            if (add_edge) {
                Edge* e = g.add_edge();
                e->set_from(edge.first.first);
                e->set_to(edge.second.first);
                e->set_from_start(edge.first.second);
                e->set_to_end(edge.second.second);
            }
        }
        
        if (cycled_to_source && preserve_cycles_on_src) {
            // we need to duplicate the source node to ensure that cycles on it will be preserved
            // after we cut it
            
            // choose an ID that hasn't been used yet
            id_t dup_id = max_id + 1;
            
            // duplicate the node
            Node* dup_node = g.add_node();
            dup_node->set_id(dup_id);
            dup_node->set_sequence(graph[id(pos)]->sequence());
            
            // record the ID translation
            id_trans[dup_id] = id(pos);
            
            for (const pair<pair<id_t, bool>, pair<id_t, bool>>& edge : src_edges) {
                Edge* e = g.add_edge();
                e->set_from_start(edge.first.second);
                e->set_to_end(edge.second.second);
                
                if (edge.first.first == id(pos) && edge.second.first == id(pos)) {
                    // this edge is a self loop, always make a copy on the duplicated node
                    e->set_from(dup_id);
                    e->set_to(dup_id);
                    
                    // if one of the ends of the edge is on the non-cut side of the source, also make a copy
                    // between the original and the duplicated node
                    if (edge.first.first == id(pos) && edge.first.second == (backward != is_rev(pos))) {
                        e = g.add_edge();
                        e->set_from(id(pos));
                        e->set_to(dup_id);
                        e->set_from_start(edge.first.second);
                        e->set_to_end(edge.second.second);
                    }
                    else if (edge.second.first == id(pos) && edge.second.second != (backward != is_rev(pos))) {
                        e = g.add_edge();
                        e->set_from(dup_id);
                        e->set_to(id(pos));
                        e->set_from_start(edge.first.second);
                        e->set_to_end(edge.second.second);
                    }
                }
                else if (edge.first.first == id(pos)) {
                    // add a copy of the edge from the other node to the duplicated node
                    e->set_from(dup_id);
                    e->set_to(edge.second.first);
                }
                else {
                    // add a copy of the edge from the other node to the duplicated node
                    e->set_from(edge.first.first);
                    e->set_to(dup_id);
                }
            }
        }
        
        // cut the source node at the starting position
        if (is_rev(pos) && backward) {
            src_node->set_sequence(src_node->sequence().substr(src_node->sequence().size() - offset(pos), offset(pos)));
        }
        else if (is_rev(pos)) {
            src_node->set_sequence(src_node->sequence().substr(0, src_node->sequence().size() - offset(pos)));
        }
        else if (backward) {
            src_node->set_sequence(src_node->sequence().substr(0, offset(pos)));
        }
        else {
            src_node->set_sequence(src_node->sequence().substr(offset(pos), src_node->sequence().size() - offset(pos)));
        }
        
        return id_trans;
    }
    
    unordered_map<id_t, id_t> extract_extending_graph(const HandleGraph* source, Graph& g, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src) {
        
        return extract_extending_graph_internal(g, max_dist, pos, backward, preserve_cycles_on_src,
                                                [&](id_t id) {
                                                    // Get edges on start
                                                    vector<Edge> to_return;
                                                    auto here = source->get_handle(id, false);
                                                    source->follow_edges(here, true, [&](const handle_t& left) -> bool {
                                                        to_return.push_back(xg::make_edge(source->get_id(left),
                                                            source->get_is_reverse(left), id, false));
                                                    });
                                                    return to_return;
                                                },
                                                [&](id_t id) {
                                                    // Get edges on start
                                                    vector<Edge> to_return;
                                                    auto here = source->get_handle(id, false);
                                                    source->follow_edges(here, false, [&](const handle_t& right) -> bool {
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
    
    unordered_map<id_t, id_t> extract_extending_graph(VG& vg, Graph& g, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src){
        
        // Just view the vg as a handle graph                                              
        return extract_extending_graph(&vg, g, max_dist, pos, backward, preserve_cycles_on_src);
                                                      
    }
    
    unordered_map<id_t, id_t> extract_extending_graph(xg::XG& xg_index, Graph& g, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src,
                                                      LRUCache<id_t, Node>* node_cache,
                                                      LRUCache<id_t, vector<Edge>>* edge_cache) {
                                                      
                                                      
        // Just view the xg as a handle graph and ignore the caches.                                          
        return extract_extending_graph(&xg_index, g, max_dist, pos, backward, preserve_cycles_on_src);
    }

}
}
