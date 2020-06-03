/**
 * \file extract_connecting_graph.cpp
 *
 * Implementation for the extract_connecting_graph algorithm.
 */
 
#include "extract_connecting_graph.hpp"
#include <structures/updateable_priority_queue.hpp>

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {

using namespace structures;

unordered_map<id_t, id_t> extract_connecting_graph(const HandleGraph* source,
                                                   DeletableHandleGraph* into,
                                                   int64_t max_len,
                                                   pos_t pos_1, pos_t pos_2,
                                                   bool detect_terminal_cycles,
                                                   bool only_paths,
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
    
    unordered_set<handle_t> skip_handles{source_handle_1};
    // mark final position for skipping so that we won't look for additional traversals unless that's
    // the only way to find terminal cycles
    if (!(colocation == SharedNodeReverse && detect_terminal_cycles)) {
        skip_handles.insert(source_handle_2);
    }
    
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
    
    // STEP 2: BACKWARD SEARCH (TO EXTRACT CYCLES ON THE FINAL NODE)
    // the forward search doesn't traverse through the second position, so we need to traverse
    // backwards from last position too if we're detecting cycles
    // also we cannot find any new nodes/edges that will pass future distance filters if
    // both forward and backward traversals are starting along the same edges, or if all paths
    // are already cyclical, so we exclude those cases to simplify some case checking in the loop
    if (detect_terminal_cycles &&
        (colocation == SeparateNodes || colocation == SharedNodeReachable)) {
        
        
#ifdef debug_vg_algorithms
        cerr << "BACKWARD SEARCH: beginning search with backward max len " << backward_max_len << " and last traversal length " << last_traversal_length << endl;
#endif
        
        // initialize the queue going backward from the last position if it's reachable
        queue.clear();
        if (last_traversal_length <= backward_max_len) {
            queue.emplace(source->flip(source_handle_2), last_traversal_length);
        }
        
        // reset the traversal list to skip and add the two reverse traversals
        skip_handles.clear();
        skip_handles.insert(source->flip(source_handle_1));
        skip_handles.insert(source->flip(source_handle_2));
        
        // search along a Dijkstra tree
        while (!queue.empty()) {
            // get the next closest node to the starting position
            Traversal trav = queue.top();
            queue.pop();
            
#ifdef debug_vg_algorithms
            cerr << "BACKWARD SEARCH: traversing node " << source->get_id(trav.handle)
                << " in " << (source->get_is_reverse(trav.handle) ? "reverse" : "forward") 
                << " orientation at distance " << trav.dist << endl;
#endif
            
            source->follow_edges(trav.handle, false, [&](const handle_t& next) {
                // get the orientation and id of the other side of the edge
                id_t next_id = source->get_id(next);
                bool next_rev = source->get_is_reverse(next);
                
#ifdef debug_vg_algorithms
                cerr << "BACKWARD SEARCH: got edge "
                    << source->get_id(trav.handle) << " " << source->get_is_reverse(trav.handle)
                    << " -> " << next_id << " " << next_rev << endl;
#endif

                max_id = max(max_id, next_id);
                
                // make sure the node is in the graph
                if (!id_trans.count(next_id)) {
                    into->create_handle(source->get_sequence(source->forward(next)), next_id);
                    id_trans[next_id] = next_id;
                }
                
                // distance to the end of this node
                int64_t dist_thru = trav.dist + source->get_length(next);
                if (!skip_handles.count(next) && dist_thru <= forward_max_len) {
                    // we can add more nodes along same path without going over the max length
                    // and we have not reached the target node yet
                    queue.emplace(next, dist_thru);
#ifdef debug_vg_algorithms
                    cerr << "BACKWARD SEARCH: distance " << dist_thru << " is under maximum, adding to queue" << endl;
#endif
                }
                
                observed_edges.insert(source->edge_handle(trav.handle, next));
            });
        }
    }
    
    // add the edges we saw
    for (const edge_t& edge : observed_edges) {
        into->create_edge(edge);
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
    
    // STEP 3: DUPLICATING NODES
    // if we're trying to detect terminal cycles, duplicate out the node so that the cyclic paths
    // survive the node cutting step
    
    auto duplicate_node = [&](const handle_t& handle, bool preserve_left_edges, bool preserve_right_edges) {
        unordered_set<edge_t> add_edges;
        
        handle_t dup_handle = into->create_handle(into->get_sequence(into->forward(handle)));
        
        // record the translation
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
    
    id_t duplicate_node_1 = 0, duplicate_node_2 = 0;
    
    if (detect_terminal_cycles) {
        // if there are edges traversed in both directions from the boundary position's nodes, then
        // they might be in a cycle
        bool has_left_edges_1 = false, has_left_edges_2 = false, has_right_edges_1 = false, has_right_edges_2 = false;
        into->follow_edges(into_handle_1, true, [&](const handle_t& ignored) {
            has_left_edges_1 = true;
            return false;
        });
        into->follow_edges(into_handle_1, false, [&](const handle_t& ignored) {
            has_right_edges_1 = true;
            return false;
        });
        into->follow_edges(into_handle_2, true, [&](const handle_t& ignored) {
            has_left_edges_2 = true;
            return false;
        });
        into->follow_edges(into_handle_2, false, [&](const handle_t& ignored) {
            has_right_edges_2 = true;
            return false;
        });
        
        bool possibly_in_cycle_1 = has_left_edges_1 && has_right_edges_1;
        bool possibly_in_cycle_2 = has_left_edges_2 && has_right_edges_2;
        
        // logic changes depending on colocation of positions on same node
        switch (colocation) {
            case SeparateNodes:
            {
                // the two positions are on separate nodes, so we can duplicate cycles independently
                
                if (possibly_in_cycle_1) {
                    duplicate_node(into_handle_1, true, true);
                }
                
                if (possibly_in_cycle_2) {
                    duplicate_node(into_handle_2, true, true);
                }
                break;
            }
            case SharedNodeReachable:
            {
                // one position is reachable from the next within the same node
                
                if (possibly_in_cycle_1) {
                    
                    // later, we're going to trim this node to it's middle portion between the two positions
                    // so now that we want to preserve cycles, we need to make two new nodes that will hold
                    // the prefix and suffix of the node so that the edges have somewhere to attach to
                    
                    
                    duplicate_node_1 = into->get_id(duplicate_node(into_handle_1, false, true));
                    duplicate_node_2 = into->get_id(duplicate_node(into_handle_2, true, false));
                    
                    // we also need a simple duplicate node in case there are any cycles that pass all the way
                    // through the node
                    
                    duplicate_node(into_handle_2, true, true);
                }
                break;
            }
            case SharedNodeUnreachable:
            case SharedNodeReverse:
            {
                // all paths between these positions are cyclical, but we still duplicate the node
                // so that any cycles that pass all the way through the node are there to be accepted
                // or rejected by the distance filter even after we cut the node up
                
                if (possibly_in_cycle_1) {
                    duplicate_node(into_handle_1, true, true);
                }
                break;
            }
        }
    }
    
#ifdef debug_vg_algorithms
    cerr << "state of graph after duplicating nodes to preserve cycles:" << endl;
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
    
    // STEP 4: CUTTING NODES
    // now cut the two end nodes at the designated positions and remove the edges on the cut side
    // to make the end positions tips in the graph
    
    handle_t cut_handle_1, cut_handle_2;
    
    switch (colocation) {
        case SeparateNodes:
        {
            // split the node, update the IDs, and clean up the other side
            auto halves_1 = into->divide_handle(into_handle_1, offset(pos_1));
            id_trans.erase(id(pos_1));
            id_trans[into->get_id(halves_1.second)] = id(pos_1);
            into->destroy_handle(halves_1.first);
            cut_handle_1 = halves_1.second;
            
            // repeat for the second position
            auto halves_2 = into->divide_handle(into_handle_2, offset(pos_2));
            id_trans.erase(id(pos_2));
            id_trans[into->get_id(halves_2.first)] = id(pos_2);
            into->destroy_handle(halves_2.second);
            cut_handle_2 = halves_2.first;
            
            break;
        }
        case SharedNodeReachable:
        {
            // split the node, update the IDs, and clean up the two ends
            auto thirds = into->divide_handle(into_handle_2, vector<size_t>{offset(pos_1), offset(pos_2)});
            id_trans.erase(id(pos_1));
            id_trans[into->get_id(thirds[1])] = id(pos_1);
            into->destroy_handle(thirds.front());
            into->destroy_handle(thirds.back());
            cut_handle_1 = thirds[1];
            cut_handle_2 = thirds[1];
            
            // if we created duplicate nodes to hold the right and left side edges in cycles, cut
            // those as well, update the IDs, and then clean up the other side
            if (duplicate_node_1) {
                handle_t dup_handle = into->get_handle(duplicate_node_1, is_rev(pos_1));
                auto halves = into->divide_handle(dup_handle, offset(pos_1));
                id_trans.erase(duplicate_node_1);
                id_trans[into->get_id(halves.second)] = id(pos_1);
                duplicate_node_1 = into->get_id(halves.second);
                into->destroy_handle(halves.first);
            }
            
            if (duplicate_node_2) {
                handle_t dup_handle = into->get_handle(duplicate_node_2, is_rev(pos_2));
                auto halves = into->divide_handle(dup_handle, offset(pos_2));
                id_trans.erase(duplicate_node_2);
                id_trans[into->get_id(halves.first)] = id(pos_2);
                duplicate_node_2 = into->get_id(halves.first);
                into->destroy_handle(halves.second);
            }
            
            break;
        }
        case SharedNodeUnreachable:
        case SharedNodeReverse:
        {
            // make a new node that will preserve the edges on the righthand side
            handle_t dup_node = duplicate_node(into_handle_1, false, true);
            auto halves_1 = into->divide_handle(dup_node, offset(pos_1));
            id_trans[into->get_id(halves_1.second)] = id(pos_1);
            into->destroy_handle(halves_1.first);
            cut_handle_1 = halves_1.second;
            
            // cut the original node and preserve its lefthand side edges
            auto halves_2 = into->divide_handle(into_handle_2, offset(pos_2));
            id_trans.erase(id(pos_2));
            id_trans[into->get_id(halves_2.first)] = id(pos_2);
            into->destroy_handle(halves_2.second);
            cut_handle_2 = halves_2.first;
            
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
    
    // STEP 5: PRUNING
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
        unordered_map<handle_t, size_t> forward_dist = find_shortest_paths(into, cut_handle_1, false);
        unordered_map<handle_t, size_t> reverse_dist = find_shortest_paths(into, cut_handle_2, true);

        // also consider minimum distances from alternate start points if we have them
        if (duplicate_node_1) {
            auto alt_forward_dist = find_shortest_paths(into, into->get_handle(duplicate_node_1, is_rev(pos_1)), false);
            for (const auto& alt_dist : alt_forward_dist) {
                auto iter = forward_dist.find(alt_dist.first);
                if (iter != forward_dist.end()) {
                    iter->second = min(iter->second, alt_dist.second);
                }
                else {
                    forward_dist.insert(alt_dist);
                }
            }
        }
        
        // also consider minimum distances from alternate start points if we have them
        if (duplicate_node_2) {
            auto alt_reverse_dist = find_shortest_paths(into, into->get_handle(duplicate_node_2, is_rev(pos_2)), true);
            for (const auto& alt_dist : alt_reverse_dist) {
                auto iter = reverse_dist.find(alt_dist.first);
                if (iter != reverse_dist.end()) {
                    iter->second = min(iter->second, alt_dist.second);
                }
                else {
                    reverse_dist.insert(alt_dist);
                }
            }
        }
        
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
    else if (only_paths) {
        // OPTION 2: PRUNE TO PATHS
        // some nodes in the current graph may not be on paths, so we do a forward-backward
        // reachability search to check
        
        auto identify_reachable = [&](handle_t start, bool search_leftward) {
            
            vector<handle_t> stack(1, start);
            unordered_set<handle_t> reachable{start};
            
#ifdef debug_vg_algorithms
            cerr << "REACHABILILTY PRUNE: beginning reachability test " << (search_leftward ? "backward" : "forward") << " from ";
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
                
                into->follow_edges(handle, search_leftward, [&](const handle_t& next) {
                    if (!reachable.count(next)) {
                        stack.emplace_back(next);
                        reachable.insert(next);
                    };
                });
            }
            
            return reachable;
        };
        
        unordered_set<handle_t> forward_reachable = identify_reachable(cut_handle_1, false);
        unordered_set<handle_t> reverse_reachable = identify_reachable(cut_handle_2, true);
        
        // now we know which nodes are reachable from both ends, to be on a path between the end positions,
        // a node or edge must be reachable from both directions
        
        // is node reachable from both positions?
        auto should_remove_node = [&](const handle_t& handle) {
            handle_t flipped = into->flip(handle);
            return !((forward_reachable.count(handle) && reverse_reachable.count(handle))
                     || (forward_reachable.count(flipped) && reverse_reachable.count(flipped)));
        };
        
        // is edge reachable from both positions?
        auto should_remove_edge = [&](const handle_t& prev, const handle_t& next) {
            handle_t flipped_prev = into->flip(next);
            handle_t flipped_next = into->flip(prev);
            
            return !((forward_reachable.count(prev) && reverse_reachable.count(next))
                     || (forward_reachable.count(flipped_prev) && reverse_reachable.count(flipped_next)));
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
    
    // if we did pruning, actually destroy the nodes
    for (const handle_t& handle : nodes_to_erase) {
        into->destroy_handle(handle);
    }
    
    // and the edges
    for (const edge_t& edge : edges_to_erase) {
        if (!nodes_to_erase.count(into->forward(edge.first))
            && !nodes_to_erase.count(into->forward(edge.second))) {
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
    return id_trans;
}

}
}
