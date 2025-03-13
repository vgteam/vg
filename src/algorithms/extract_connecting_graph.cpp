/**
 * \file extract_connecting_graph.cpp
 *
 * Implementation for the extract_connecting_graph algorithm.
 */
 
#include "extract_connecting_graph.hpp"
#include <structures/updateable_priority_queue.hpp>
#include "../crash.hpp"

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {

using namespace structures;

unordered_map<id_t, id_t> extract_connecting_graph(const HandleGraph* source,
                                                   DeletableHandleGraph* into,
                                                   int64_t max_len,
                                                   pos_t pos_1, pos_t pos_2,
                                                   bool strict_max_len) {
#ifdef debug_vg_algorithms
    cerr << "[extract_connecting_graph] max len: " << max_len << ", pos 1: " << pos_1 << ", pos 2: " << pos_2 << endl;
#endif
    
    if (into->get_node_count()) {
        cerr << "error:[extract_connecting_graph] must extract into an empty graph" << endl;
        exit(1);
    }
    
    // a local struct that packages a handle with its distance from the first position
    struct Traversal {
        Traversal(handle_t handle, int64_t dist) : handle(handle), dist(dist) {}
        int64_t dist; // distance from pos to the right side of this node
        handle_t handle; // Oriented node traversal
        inline bool operator<(const Traversal& other) const {
            return dist > other.dist; // opposite order so priority queue selects minimum
        }
    };
    
    // local enum to keep track of the cases where the positions are on the same node
    enum colocation_t {SeparateNodes, SharedNodeReachable, SharedNodeUnreachable, SharedNodeReverse};
    
    // record whether the positions are on the same node, and if so their relationship to each other
    colocation_t colocation;
    if (id(pos_1) == id(pos_2)) {
        if (is_rev(pos_1) == is_rev(pos_2)) {
            if (offset(pos_1) <= offset(pos_2)) {
                colocation = SharedNodeReachable;
            }
            else {
                colocation = SharedNodeUnreachable;
            }
        }
        else {
            colocation = SharedNodeReverse;
        }
    }
    else {
        colocation = SeparateNodes;
    }
    
    handle_t source_handle_1 = source->get_handle(id(pos_1), is_rev(pos_1));
    handle_t source_handle_2 = source->get_handle(id(pos_2), is_rev(pos_2));
    
    // for finding the largest node id in the subgraph
    id_t max_id = max(id(pos_1), id(pos_2));
    
    // a translator for node ids in 'into' to node ids in 'source'
    unordered_map<id_t, id_t> id_trans;
    
    // the edges we have encountered in the traversal
    unordered_set<edge_t> observed_edges;
    
    // create nodes for the source positions in the new graph
    handle_t into_handle_1 = into->create_handle(source->get_sequence(source->forward(source_handle_1)), id(pos_1));
    id_trans[id(pos_1)] = id(pos_1);
    handle_t into_handle_2;
    if (id(pos_1) != id(pos_2)) {
        into_handle_2 = into->create_handle(source->get_sequence(source->forward(source_handle_2)), id(pos_2));
        id_trans[id(pos_2)] = id(pos_2);
    }
    else {
        into_handle_2 = into_handle_1;
    }
    
    // make the handles match the orientation of the positions
    if (is_rev(pos_1)) {
        into_handle_1 = into->flip(into_handle_1);
    }
    if (is_rev(pos_2)) {
        into_handle_2 = into->flip(into_handle_2);
    }
    
    // keep track of whether we find a path or not
    bool found_target = false;
    
    // mark final position for skipping so that we won't look for additional traversals
    unordered_set<handle_t> skip_handles{source_handle_1, source_handle_2};
    
    // initialize the queue
    UpdateablePriorityQueue<Traversal, handle_t> queue([](const Traversal& item) {
        return item.handle;
    });
    
    // the distance to the ends of the starting nodes
    int64_t first_traversal_length = source->get_length(source_handle_1) - offset(pos_1);
    int64_t last_traversal_length = offset(pos_2);
    
    // the max length of the part of a path preceding the final node in each direction
    int64_t forward_max_len = max_len - last_traversal_length;
    int64_t backward_max_len = max_len - first_traversal_length;
    
    // STEP 1: FORWARD SEARCH (TO EXTRACT SUBGRAPH)
    // separately handle (common) edge case that both positions are on the same node
    // and the second is reachable from the first
    // TODO: is there a more elegant way to do this than as a special case? currently I need to
    // do it because the objects in the queue are implicitly at the "end" of a node traversal,
    // but here we would want to stop before reaching the end of the first node
    if (colocation == SharedNodeReachable) {
        found_target = offset(pos_2) - offset(pos_1) <= max_len;
#ifdef debug_vg_algorithms
        cerr << "FORWARD SEARCH: positions are on same node, skipping forward search and identifying them as " << (found_target ? "" : "not ") << "reachable" << endl;
#endif
    }
    else {
        // search through graph to find the target, or to find cycles involving this node
        
#ifdef debug_vg_algorithms
        cerr << "FORWARD SEARCH: beginning search with forward max len " << forward_max_len << " and first traversal length " << first_traversal_length << endl;
#endif
        
        // if we can reach the end of this node, init the queue with it
        if (first_traversal_length <= forward_max_len) {
            queue.emplace(source_handle_1, first_traversal_length);
        }
        
        // search along a Dijkstra tree
        while (!queue.empty()) {
            // get the next closest node to the starting position
            Traversal trav = queue.top();
            queue.pop();
            
#ifdef debug_vg_algorithms
            cerr << "FORWARD SEARCH: traversing node " << source->get_id(trav.handle) << " in "
                << (source->get_is_reverse(trav.handle) ? "reverse" : "forward")
                << " orientation at distance " << trav.dist << endl;
#endif
            
            source->follow_edges(trav.handle, false, [&](const handle_t& next) {
                // get the orientation and id of the other side of the edge
                
                id_t next_id = source->get_id(next);
                bool next_rev = source->get_is_reverse(next);
                
#ifdef debug_vg_algorithms
                cerr << "FORWARD SEARCH: got edge "
                    << source->get_id(trav.handle) << " " << source->get_is_reverse(trav.handle)
                    << " -> " << next_id << " " << next_rev << endl;
#endif
                found_target = found_target || (next_id == id(pos_2) && next_rev == is_rev(pos_2));
                max_id = max(max_id, next_id);
                
                // make sure the node is in
                if (!id_trans.count(next_id)) {
                    // Make a node with the forward orientation sequence
                    into->create_handle(source->get_sequence(source->forward(next)), next_id);
                    id_trans[next_id] = next_id;
                }
                
                // distance to the end of this node
                int64_t dist_thru = trav.dist + source->get_length(next);
                if (!skip_handles.count(next) && dist_thru <= forward_max_len) {
                    // we can add more nodes along same path without going over the max length
                    // and we do not want to skip the target node
                    queue.emplace(next, dist_thru);
#ifdef debug_vg_algorithms
                    cerr << "FORWARD SEARCH: distance " << dist_thru << " is under maximum, adding to queue" << endl;
#endif
                }
                
                observed_edges.insert(source->edge_handle(trav.handle, next));
            });
        }
    }
    
    // there is no path between the nodes under the maximum distance, leave g empty and return
    // an empty translator
    if (!found_target) {
        into->clear();
        id_trans.clear();
        return id_trans;
    }
    
    // add the edges we saw
    for (const edge_t& edge : observed_edges) {
        into->create_edge(into->get_handle(source->get_id(edge.first), source->get_is_reverse(edge.first)),
                          into->get_handle(source->get_id(edge.second), source->get_is_reverse(edge.second)));
    }
    
#ifdef debug_vg_algorithms
    cerr << "state of graph after forward and backward search:" << endl;
    into->for_each_handle([&](const handle_t& handle) {
        cerr << "node " << into->get_id(handle) << " " << into->get_sequence(handle) << endl;
        cerr << "\tleft" << endl;
        into->follow_edges(handle, true, [&](const handle_t& next) {
            cerr << "\t\t" << into->get_id(next) << (into->get_is_reverse(next) ? "-" : "+") << endl;
        });
        cerr << "\tright" << endl;
        into->follow_edges(handle, false, [&](const handle_t& next) {
            cerr << "\t\t" << into->get_id(next) << (into->get_is_reverse(next) ? "-" : "+") << endl;
        });
    });
#endif
    
    // STEP 2: DUPLICATING NODES
    // if we're trying to detect terminal cycles, duplicate out the node so that the cyclic paths
    // survive the node cutting step
    
    auto duplicate_node = [&](const handle_t& handle, bool preserve_left_edges, bool preserve_right_edges) {
        unordered_set<edge_t> add_edges;
        
        // Do the duplication
        handle_t dup_handle = into->create_handle(into->get_sequence(into->forward(handle)));
        if (into->get_is_reverse(handle)) {
            // Make sure that we don't forget our orientation relative to the
            // local forward.
            dup_handle = into->flip(dup_handle);
        }
        // Record the translation
        id_trans[into->get_id(dup_handle)] = id_trans[into->get_id(handle)];
        
        if (preserve_right_edges) {
            // collect rightward edges
            into->follow_edges(handle, false, [&](const handle_t& next) {
                add_edges.insert(into->edge_handle(dup_handle, next));
                
                if (into->get_id(next) == into->get_id(handle)) {
                    handle_t dup_next = next == handle ? dup_handle : into->flip(dup_handle);
                    // this is a loop, also add a connection between the original node
                    add_edges.insert(into->edge_handle(handle, dup_next));
                    // and with the duplicated node onto itself
                    add_edges.insert(into->edge_handle(dup_handle, dup_next));
                }
                
            });
        }
        
        if (preserve_left_edges) {
            // collect leftward edges
            into->follow_edges(handle, true, [&](const handle_t& prev) {
                add_edges.insert(into->edge_handle(prev, dup_handle));
                
                if (into->get_id(prev) == into->get_id(handle)) {
                    handle_t dup_prev = prev == handle ? dup_handle : into->flip(dup_handle);
                    // this is a loop, also add a connection between the original node
                    add_edges.insert(into->edge_handle(dup_prev, handle));
                    // and with the duplicated node onto itself
                    add_edges.insert(into->edge_handle(dup_prev, dup_handle));
                }
            });
        }
        
        // now add the edges in
        for (const edge_t& edge : add_edges) {
            into->create_edge(edge);
        }
        
        return dup_handle;
    };
    
    // STEP 3: CUTTING NODES
    // now cut the two end nodes at the designated positions and remove the edges on the cut side
    // to make the end positions tips in the graph
    //
    // We need to guarantee that, if two separate end nodes came from one
    // original graph node, we assign the left one the lower ID.
    
    handle_t cut_handle_1, cut_handle_2;
    
    switch (colocation) {
        case SeparateNodes:
        {
            // split the node, update the IDs, and clean up the other side
            cut_handle_1 = into->truncate_handle(into_handle_1, true, offset(pos_1));
            id_trans.erase(id(pos_1));
            id_trans[into->get_id(cut_handle_1)] = id(pos_1);
            
            // repeat for the second position
            cut_handle_2 = into->truncate_handle(into_handle_2, false, offset(pos_2));
            id_trans.erase(id(pos_2));
            id_trans[into->get_id(cut_handle_2)] = id(pos_2);
            
            break;
        }
        case SharedNodeReachable:
        {
            // Split the node, update the IDs, and clean up the two ends.
            // We know the two handles are in the same orientation on the node,
            // or they wouldn't be reachable. So we know the two offsets are
            // along the same strand.
            crash_unless(is_rev(pos_2) == is_rev(pos_1));
            cut_handle_1 = into->truncate_handle(into->truncate_handle(into_handle_2, false, offset(pos_2)), true, offset(pos_1));
            id_trans.erase(id(pos_1));
            id_trans[into->get_id(cut_handle_1)] = id(pos_1);
            // We have one shared end node
            cut_handle_2 = cut_handle_1;
            break;
        }
        case SharedNodeUnreachable:
        case SharedNodeReverse:
        {
            // make a new node that will preserve the edges on the lefthand side
            handle_t dup_node = duplicate_node(into_handle_2, true, false);
            cut_handle_2 = into->truncate_handle(dup_node, false, offset(pos_2));
            id_trans[into->get_id(cut_handle_2)] = id(pos_2);

            // cut the original node and preserve its righthand side edges
            cut_handle_1 = into->truncate_handle(into_handle_1, true, offset(pos_1));
            id_trans.erase(id(pos_1));
            id_trans[into->get_id(cut_handle_1)] = id(pos_1);
            
            if (into->get_id(cut_handle_2) < into->get_id(cut_handle_1)) {
                // We assume that cut_handle_1 will get the lower ID. Make sure that's always true.
                throw std::runtime_error("Graph assigned end node a lower ID than start node. Caller will not be able to identify them properly.");
            }
            
            break;
        }
    }
    
#ifdef debug_vg_algorithms
    cerr << "state of graph after cutting nodes:" << endl;
    into->for_each_handle([&](const handle_t& handle) {
        cerr << "node " << into->get_id(handle) << " " << into->get_sequence(handle) << endl;
        cerr << "\tleft" << endl;
        into->follow_edges(handle, true, [&](const handle_t& next) {
            cerr << "\t\t" << into->get_id(next) << (into->get_is_reverse(next) ? "-" : "+") << endl;
        });
        cerr << "\tright" << endl;
        into->follow_edges(handle, false, [&](const handle_t& next) {
            cerr << "\t\t" << into->get_id(next) << (into->get_is_reverse(next) ? "-" : "+") << endl;
        });
    });
#endif
    
    // STEP 4: PRUNING
    // the graph now contains all the paths we've indicated and the end positions are tips, we now
    // provide three options for pruning away any unnecessary nodes and edges we've added in the
    // process of searching for the subgraph that has this guarantee
    
    // these will be filled by the pruning algorithms
    unordered_set<handle_t> nodes_to_erase;
    unordered_set<edge_t> edges_to_erase;
    
    if (strict_max_len) {
        // OPTION 1: PRUNE TO PATHS UNDER MAX LENGTH
        // some nodes in the current graph may not be on paths, or the paths that they are on may be
        // above the maximum distance, so we do a forward-backward distance search to check
        
        // compute the minimum distance from the two start point s
        unordered_map<handle_t, size_t> forward_dist = handlealgs::find_shortest_paths(into, cut_handle_1, false);
        unordered_map<handle_t, size_t> reverse_dist = handlealgs::find_shortest_paths(into, cut_handle_2, true);
        
        // now we have the lengths of the shortest path remaining in graph to and from each node
        // with these, we can compute the shortest path that uses each node and edge to see if it
        // should be included in the final graph
        
        // forward-backward test for whether a node is used in a sufficiently short walk
        auto should_remove_node = [&](const handle_t& handle) {
            handle_t flipped = into->flip(handle);
            
            bool remove_node = true;
            if (forward_dist.count(handle) && reverse_dist.count(handle)) {
                size_t min_dist_thru_node = (forward_dist[handle]
                                             + into->get_length(handle)
                                             + reverse_dist[handle]);
                if (min_dist_thru_node <= max_len) {
                    remove_node = false;
                }
            }
            
            if (forward_dist.count(flipped) && reverse_dist.count(flipped)) {
                size_t min_dist_thru_node = (forward_dist[flipped]
                                             + into->get_length(flipped)
                                             + reverse_dist[flipped]);
                if (min_dist_thru_node <= max_len) {
                    remove_node = false;
                }
            }
            
            return remove_node;
        };
        
        // forward-backward test for whether an edge is used in a sufficiently short walk
        auto should_remove_edge = [&](const handle_t& prev, const handle_t& next) {
            handle_t flipped_prev = into->flip(next);
            handle_t flipped_next = into->flip(prev);
            
            bool remove_edge = true;
            if (forward_dist.count(prev) && reverse_dist.count(next)) {
                size_t min_dist_over_edge = (forward_dist[prev]
                                             + into->get_length(prev)
                                             + into->get_length(next)
                                             + reverse_dist[next]);
                if (min_dist_over_edge <= max_len) {
                    remove_edge = false;
                }
            }
            
            if (forward_dist.count(flipped_prev) && reverse_dist.count(flipped_next)) {
                size_t min_dist_over_edge = (forward_dist[flipped_prev]
                                             + into->get_length(flipped_prev)
                                             + into->get_length(flipped_next)
                                             + reverse_dist[flipped_next]);
                if (min_dist_over_edge <= max_len) {
                    remove_edge = false;
                }
            }
            
            return remove_edge;
        };
        
        // apply the tests to each node/edge and collect the results
        into->for_each_handle([&](const handle_t& handle) {
            if (should_remove_node(handle)) {
                nodes_to_erase.insert(handle);
            }
            else {
                into->follow_edges(handle, false, [&](const handle_t& next) {
                    if (should_remove_edge(handle, next)) {
                        edges_to_erase.insert(into->edge_handle(handle, next));
                    }
                });
                
                into->follow_edges(handle, true, [&](const handle_t& prev) {
                    if (should_remove_edge(prev, handle)) {
                        edges_to_erase.insert(into->edge_handle(prev, handle));
                    }
                });
            }
        });
    }
    else {
        // OPTION 2: PRUNE TO PATHS
        // some nodes in the current graph may not be on paths between the cut points,
        // so we do a reverse traversal (the initial traversal establishes forward reachability)
        
        unordered_set<handle_t> reverse_reachable{cut_handle_2};
        reverse_reachable.reserve(into->get_node_count());
        vector<handle_t> stack(1, cut_handle_2);
        
#ifdef debug_vg_algorithms
        cerr << "REACHABILILTY PRUNE: beginning reachability test backward from ";
        for (auto handle : stack) {
            cerr << into->get_id(handle) << (into->get_is_reverse(handle) ? "-" : "+") << " ";
        }
        cerr << endl;
#endif
        
        while (!stack.empty()) {
            handle_t handle = stack.back();
            stack.pop_back();
            
#ifdef debug_vg_algorithms
            cerr << "REACHABILILTY PRUNE: traversing node " << into->get_id(handle) << " in " << (into->get_is_reverse(handle) ? "reverse" : "forward") << endl;
#endif
            
            into->follow_edges(handle, true, [&](const handle_t& next) {
                if (!reverse_reachable.count(next)) {
                    stack.emplace_back(next);
                    reverse_reachable.insert(next);
                };
            });
        }
        
        // now we know which nodes are reachable from both ends, to be on a path between the end positions,
        // a node or edge must be reachable from both directions
        
        // apply the tests to each node and collect the results
        into->for_each_handle([&](const handle_t& handle) {
            if (!reverse_reachable.count(handle) && !reverse_reachable.count(into->flip(handle))) {
                nodes_to_erase.insert(handle);
            }
        });
    }
    
    // if we did pruning, actually destroy the nodes
    for (const handle_t& handle : nodes_to_erase) {
        nid_t id = into->get_id(handle);
        into->destroy_handle(handle);
        id_trans.erase(id);
    }
    
    // and the edges
    for (const edge_t& edge : edges_to_erase) {
        if (!nodes_to_erase.count(edge.first)
            && !nodes_to_erase.count(into->flip(edge.first))
            && !nodes_to_erase.count(edge.second)
            && !nodes_to_erase.count(into->flip(edge.second))) {
            into->destroy_edge(edge);
        }
    }
    
#ifdef debug_vg_algorithms
    cerr << "state of graph after pruning:" << endl;
    into->for_each_handle([&](const handle_t& handle) {
        cerr << "node " << into->get_id(handle) << " " << into->get_sequence(handle) << endl;
        cerr << "\tleft" << endl;
        into->follow_edges(handle, true, [&](const handle_t& next) {
            cerr << "\t\t" << into->get_id(next) << (into->get_is_reverse(next) ? "-" : "+") << endl;
        });
        cerr << "\tright" << endl;
        into->follow_edges(handle, false, [&](const handle_t& next) {
            cerr << "\t\t" << into->get_id(next) << (into->get_is_reverse(next) ? "-" : "+") << endl;
        });
    });
#endif
    
    // TODO: it's not enough to return the translator because there's also the issue of the positions
    // on the first node being offset (however this information is fully contained in the arguments of
    // the function, which are obviously available in the environment that calls it)
    crash_unless(id_trans.size() == into->get_node_count());
    return id_trans;
}

}
}
