/**
 * \file extract_extending_graph.cpp
 *
 * Implementation for the extract_extending_graph algorithm.
 */
 
#include "extract_extending_graph.hpp"
#include <structures/updateable_priority_queue.hpp>
#include "../crash.hpp"

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {

using namespace structures;

unordered_map<id_t, id_t> extract_extending_graph(const HandleGraph* source, DeletableHandleGraph* into, int64_t max_dist, pos_t pos,
                                                  bool backward, bool preserve_cycles_on_src_node) {
    
    if (into->get_node_count()) {
        cerr << "error:[extract_extending_graph] must extract into an empty graph" << endl;
        assert(false);
    }
    
#ifdef debug_vg_algorithms
    cerr << "[extract_extending_graph] extracting exending graph from " << pos << " in " << (backward ? "backward" : "forward") << " direction with max search dist " << max_dist << endl;
#endif
    
    // TODO: struct is duplicative with extract_connecting_graph
    
    // a local struct that packages a handle with its distance from the first position
    struct Traversal {
        Traversal(handle_t handle, int64_t dist) : handle(handle), dist(dist) {}
        int64_t dist; // distance from pos to the right side of this node
        handle_t handle; // Oriented node traversal
        inline bool operator<(const Traversal& other) const {
            return dist > other.dist; // opposite order so priority queue selects minimum
        }
    };
    
    // A map from node ids in the extracted graph to the node ids in the original graph.
    // The IDs are always the same except for when we cut/duplicate the source node.
    unordered_map<id_t, id_t> id_trans;
    
    // a graph index that we will maintain as we extract the subgraph
    handle_t search_origin = source->get_handle(id(pos));
    handle_t src_node = into->create_handle(source->get_sequence(search_origin), source->get_id(search_origin));
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
            start_handle = is_rev(pos) ? search_origin : source->flip(search_origin);
            queue.emplace(start_handle, dist);
        }
    }
    else {
        int64_t dist = source->get_length(search_origin) - offset(pos);
        if (dist < max_dist) {
            // We need to leave the node to get enough sequence
            start_handle = is_rev(pos) ? source->flip(search_origin) : search_origin;
            queue.emplace(start_handle, dist);
        }
    }
    
    id_t max_id = id(pos);
    bool cycled_to_source = false;
    unordered_set<edge_t> observed_edges;
    unordered_set<id_t> observed_nodes{id(pos)};
    
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
            id_t next_id = source->get_id(next);
            
            // do we ever return to the source node?
            cycled_to_source = cycled_to_source || next_id == id(pos);
            
            // record the edge
            observed_edges.insert(source->edge_handle(trav.handle, next));
            
            // make sure the node is in the graph
            if (!observed_nodes.count(next_id)) {
                // Grab the sequence in its local forward orientation.
                into->create_handle(source->get_sequence(source->forward(next)), next_id);
                id_trans[next_id] = next_id;
                observed_nodes.insert(next_id);
                max_id = std::max(next_id, max_id);
            }
            
            // distance to the end of this node
            int64_t dist_thru = trav.dist + source->get_length(next);
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
    
    vector<edge_t> src_edges;
    
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
            into->create_edge(into->get_handle(source->get_id(edge.first), source->get_is_reverse(edge.first)), into->get_handle(source->get_id(edge.second), source->get_is_reverse(edge.second)));
        }
    }
    
    if (cycled_to_source && preserve_cycles_on_src_node) {
#ifdef debug_vg_algorithms
        cerr << "[extract_extending_graph] we cycled to the source and we want to preserve cycles, duplicating source node" << endl;
#endif
        
        // we need to duplicate the source node to ensure that cycles on it will be preserved
        // after we cut it
        
        // choose an ID that hasn't been used yet
        id_t dup_id = max_id + 1;
        
        // duplicate the node and record the ID translation
        handle_t dup_node = into->create_handle(source->get_sequence(search_origin), dup_id);
        id_trans[dup_id] = id(pos);
        
        for (const edge_t& edge : src_edges) {
            
#ifdef debug_vg_algorithms
            cerr << "[extract_extending_graph] " << "adding edge " << source->get_id(edge.first) << (source->get_is_reverse(edge.first) ? "-" : "+") << " -> " << source->get_id(edge.second) << (source->get_is_reverse(edge.second) ? "-" : "+") << " to source node copy" << endl;
#endif
            
            if (source->get_id(edge.first) == id(pos) && source->get_id(edge.second) == id(pos)) {
                // this edge is a self loop, always make a copy on the duplicated node
                into->create_edge(source->get_is_reverse(edge.first) ? into->flip(dup_node) : dup_node,
                                  source->get_is_reverse(edge.second) ? into->flip(dup_node) : dup_node);
                
                // if one of the ends of the edge is on the non-cut side of the source, also make a copy
                // between the original and the duplicated node
                if (source->get_is_reverse(edge.first) == (backward != is_rev(pos))) {
                    into->create_edge(source->get_is_reverse(edge.first) ? into->flip(src_node) : src_node,
                                      source->get_is_reverse(edge.second) ? into->flip(dup_node) : dup_node);
                }
                else if (source->get_is_reverse(edge.second) != (backward != is_rev(pos))) {
                    into->create_edge(source->get_is_reverse(edge.first) ? into->flip(dup_node) : dup_node,
                                      source->get_is_reverse(edge.second) ? into->flip(src_node) : src_node);
                }
            }
            else if (source->get_id(edge.first) == id(pos)) {
                // add a copy of the edge from the other node to the duplicated node
                into->create_edge(source->get_is_reverse(edge.first) ? into->flip(dup_node) : dup_node,
                                  into->get_handle(source->get_id(edge.second), source->get_is_reverse(edge.second)));
            }
            else {
                // add a copy of the edge from the other node to the duplicated node
                into->create_edge(into->get_handle(source->get_id(edge.first), source->get_is_reverse(edge.first)),
                                  source->get_is_reverse(edge.second) ? into->flip(dup_node) : dup_node);
            }
        }
    }
    
    
    // cut the source node at the starting position and record a new translation
    size_t forward_offset = is_rev(pos) ? source->get_length(search_origin) - offset(pos) : offset(pos);
    pair<handle_t, handle_t> halves = into->divide_handle(src_node, forward_offset);
    id_trans.erase(into->get_id(src_node));
    if (is_rev(pos) == backward) {
        into->destroy_handle(halves.first);
        id_trans[into->get_id(halves.second)] = id(pos);
    }
    else {
        into->destroy_handle(halves.second);
        id_trans[into->get_id(halves.first)] = id(pos);
    }
    
    crash_unless(id_trans.size() == into->get_node_count());
    return id_trans;
}

}
}
