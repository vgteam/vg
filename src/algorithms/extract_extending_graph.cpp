/**
 * \file extract_extending_graph.cpp
 *
 * Implementation for the extract_extending_graph algorithm.
 */
 
#include "extract_extending_graph.hpp"
#include <structures/updateable_priority_queue.hpp>

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {

using namespace structures;

unordered_map<id_t, id_t> extract_extending_graph(const HandleGraph* source, Graph& g, int64_t max_dist, pos_t pos,
                                                  bool backward, bool preserve_cycles_on_src) {
    
    if (g.node_size() || g.edge_size()) {
        cerr << "error:[extract_extending_graph] must extract into an empty graph" << endl;
        assert(false);
    }
    
#ifdef debug_vg_algorithms
    cerr << "[extract_extending_graph] extracting exending graph from " << pos << " in " << (backward ? "backward" : "forward") << " direction with max search dist " << max_dist << endl;
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
    
    // a map from node ids in the extracted graph to the node ids in the original graph
    unordered_map<id_t, id_t> id_trans;
    
    // a graph index that we will maintain as we extract the subgraph
    unordered_map<id_t, Node*> graph;
    Node* src_node = g.add_node();
    src_node->set_sequence(source->get_sequence(source->get_handle(id(pos))));
    src_node->set_id(id(pos));
    graph[id(pos)] = src_node;
    id_trans[id(pos)] = id(pos);
    
    // Keep a handle to where we start from.
    // If the queue stays empty and we do no searching, the handle will be
    // uninitialized.
    handle_t start_handle;
    
    // initialize the queue for Dijkstra traversal.
    UpdateablePriorityQueue<Traversal, handle_t> queue([](const Traversal& item) {
        return item.handle;
    });
    
    if (backward) {
        int64_t dist = offset(pos);
        if (dist < max_dist) {
            // We need to leave the node to get enough sequence
            start_handle = source->get_handle(id(pos), !is_rev(pos));
            queue.emplace(start_handle, dist);
        }
    }
    else {
        int64_t dist = g.node(0).sequence().size() - offset(pos);
        if (dist < max_dist) {
            // We need to leave the node to get enough sequence
            start_handle = source->get_handle(id(pos), is_rev(pos));
            queue.emplace(start_handle, dist);
        }
    }
    
    id_t max_id = id(pos);
    bool cycled_to_source = false;
    unordered_set<pair<handle_t, handle_t>> observed_edges;
    
    while (!queue.empty()) {
        // get the next shortest distance traversal from either the init
        Traversal trav = queue.top();
        queue.pop();
        
#ifdef debug_vg_algorithms
        cerr << "[extract_extending_graph] traversing " << source->get_id(trav.handle)
            << (source->get_is_reverse(trav.handle) ? "-" : "+") << " at dist " << trav.dist << endl;
#endif
        
        // Now go out the right local side
        source->follow_edges(trav.handle, false, [&](const handle_t& next) {
            // For each next handle...
            
#ifdef debug_vg_algorithms
            cerr << "[extract_extending_graph] checking edge " << source->get_id(trav.handle)
            << (source->get_is_reverse(trav.handle) ? "-" : "+") << " -> " << source->get_id(next)
            << (source->get_is_reverse(next) ? "-" : "+") << endl;
#endif
            
            // Get the ID of where we're going.
            auto next_id = source->get_id(next);
            
            // do we ever return to the source node?
            cycled_to_source = cycled_to_source || next_id == id(pos);
            
            // record the edge
            observed_edges.insert(source->edge_handle(trav.handle, next));
            
            // make sure the node is in the graph
            if (!graph.count(next_id)) {
                Node* node = g.add_node();
                // Grab the sequence in its local forward orientation.
                node->set_sequence(source->get_sequence(source->forward(next)));
                node->set_id(next_id);
                graph[next_id] = node;
                id_trans[next_id] = next_id;
                max_id = std::max(next_id, max_id);
            }
            
            // distance to the end of this node
            int64_t dist_thru = trav.dist + graph[next_id]->sequence().size();
            if (dist_thru < max_dist) {
                // we can add more nodes along same path without going over the max length
                queue.emplace(next, dist_thru);
#ifdef debug_vg_algorithms
                cerr << "[extract_extending_graph] enqueuing " << source->get_id(next)
                    << (source->get_is_reverse(next) ? "-" : "+") << " at dist " << dist_thru << endl;
#endif
            }
        });
    }
    
    vector<pair<handle_t, handle_t>> src_edges;
    
    // add the edges to the graph
    for (const pair<handle_t, handle_t>& edge : observed_edges) {
        // This loop will never run unless we actually ran the search, so it is safe to use start_handle here.
        bool add_edge;
        if (source->forward(edge.first) == source->forward(start_handle) ||
            source->forward(edge.second) == source->forward(start_handle)) {
            // record the edges that touch the source in case we need to duplicate them
            src_edges.push_back(edge);
            // does the edge only touch the side of the source node that's not going to be cut?
            // That means it must touch only the right side of our starting handle.
            // We know it touches one side, so make sure it *doesn't* touch the left side.
            // Which means we need the reverse of our starting handle to not
            // be the left handle, and our handle to not be the right
            // handle.
            add_edge = !(edge.first == source->flip(start_handle) || edge.second == start_handle);
#ifdef debug_vg_algorithms
            cerr << "[extract_extending_graph] " << (add_edge ? "" : "not ") << "adding edge " << source->get_id(edge.first) << (source->get_is_reverse(edge.first) ? "-" : "+") << " -> " << source->get_id(edge.second) << (source->get_is_reverse(edge.second) ? "-" : "+") << " that touches source node" << endl;
#endif
        }
        else {
            // the edge doesn't touch the source node, so it will certainly survive the cut
            add_edge = true;
#ifdef debug_vg_algorithms
            cerr << "[extract_extending_graph] " << "adding edge " << source->get_id(edge.first) << (source->get_is_reverse(edge.first) ? "-" : "+") << " -> " << source->get_id(edge.second) << (source->get_is_reverse(edge.second) ? "-" : "+") << " that does not touch source node" << endl;
#endif
        }
        
        if (add_edge) {
            Edge* e = g.add_edge();
            e->set_from(source->get_id(edge.first));
            e->set_from_start(source->get_is_reverse(edge.first));
            e->set_to(source->get_id(edge.second));
            e->set_to_end(source->get_is_reverse(edge.second));
        }
    }
    
    if (cycled_to_source && preserve_cycles_on_src) {
#ifdef debug_vg_algorithms
        cerr << "[extract_extending_graph] we cycled to the source and we want to preserve cycles, duplicating source node" << endl;
#endif
        
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
        
        for (const pair<handle_t, handle_t>& edge : src_edges) {
            Edge* e = g.add_edge();
            e->set_from_start(source->get_is_reverse(edge.first));
            e->set_to_end(source->get_is_reverse(edge.second));
            
            // Load the IDs for this edge
            auto from_id = source->get_id(edge.first);
            auto to_id = source->get_id(edge.second);
            
#ifdef debug_vg_algorithms
            cerr << "[extract_extending_graph] " << "adding edge " << source->get_id(edge.first) << (source->get_is_reverse(edge.first) ? "-" : "+") << " -> " << source->get_id(edge.second) << (source->get_is_reverse(edge.second) ? "-" : "+") << " to source node copy" << endl;
#endif
            
            if (from_id == id(pos) && to_id == id(pos)) {
                // this edge is a self loop, always make a copy on the duplicated node
                e->set_from(dup_id);
                e->set_to(dup_id);
                
                // if one of the ends of the edge is on the non-cut side of the source, also make a copy
                // between the original and the duplicated node
                if (from_id == id(pos) && source->get_is_reverse(edge.first) == (backward != is_rev(pos))) {
                    e = g.add_edge();
                    e->set_from(id(pos));
                    e->set_to(dup_id);
                    e->set_from_start(source->get_is_reverse(edge.first));
                    e->set_to_end(source->get_is_reverse(edge.second));
                }
                else if (to_id == id(pos) && source->get_is_reverse(edge.second) != (backward != is_rev(pos))) {
                    e = g.add_edge();
                    e->set_from(dup_id);
                    e->set_to(id(pos));
                    e->set_from_start(source->get_is_reverse(edge.first));
                    e->set_to_end(source->get_is_reverse(edge.second));
                }
            }
            else if (from_id == id(pos)) {
                // add a copy of the edge from the other node to the duplicated node
                e->set_from(dup_id);
                e->set_to(to_id);
            }
            else {
                // add a copy of the edge from the other node to the duplicated node
                e->set_from(from_id);
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

}
}
