//
//  vg_algorithms.cpp
//
//

#include "vg_algorithms.hpp"

//#define debug_vg_algorithms

namespace vg {
namespace algorithms {
    unordered_map<id_t, id_t> extract_connecting_graph_internal(Graph& g, int64_t max_len,
                                                                pos_t pos_1, pos_t pos_2,
                                                                bool include_terminal_positions,
                                                                bool detect_terminal_cycles,
                                                                bool no_additional_tips,
                                                                bool only_paths,
                                                                bool strict_max_len,
                                                                function<vector<Edge>(id_t)> edges_on_start,
                                                                function<vector<Edge>(id_t)> edges_on_end,
                                                                function<string(id_t)> node_sequence) {
#ifdef debug_vg_algorithms
        cerr << "[extract_connecting_graph] max len: " << max_len << ", pos 1: " << pos_1 << ", pos 2: " << pos_2 << endl;
#endif
        
        if (g.node_size() || g.edge_size()) {
            cerr << "error:[extract_connecting_graph] must extract into an empty graph" << endl;
            exit(1);
        }
        
        // a local struct for Nodes that maintains edge lists
        struct LocalNode {
            LocalNode() {}
            LocalNode(string sequence) : sequence(sequence) {}
            string sequence;
            // edges are stored as (node id, is reversing?)
            vector<pair<id_t, bool>> edges_left;
            vector<pair<id_t, bool>> edges_right;
        };
        
        // a local struct that packages a node traversal with its distance from the first position
        // note: we don't use NodeTraversals because we may be using an XG, in which case we don't
        // have Node*'s
        struct Traversal {
            Traversal(id_t id, bool rev, int64_t dist) : id(id), rev(rev), dist(dist) {}
            int64_t dist; // distance from pos_1 to the right side of this node
            id_t id;      // node ID
            bool rev;     // strand
            inline bool operator<(const Traversal& other) const {
                return dist > other.dist; // opposite order so priority queue selects minimum
            }
        };
        
        // local enum to keep track of the cases where the positions are on the same node
        enum colocation_t {SeparateNodes, SharedNodeReachable, SharedNodeUnreachable, SharedNodeReverse};
        
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
        
        
        // functions to extract the part of the node string past the first and second positions:
        // get sequence to the right
        auto trimmed_seq_right = [&](const string& seq, int64_t offset, bool rev) {
            if (rev) {
                return seq.substr(0, seq.size() - offset - 1 + include_terminal_positions);
            }
            else {
                return seq.substr(offset + 1 - include_terminal_positions,
                                  seq.size() - offset - 1 + include_terminal_positions);
            }
        };
        // get sequence to the left
        auto trimmed_seq_left = [&](const string& seq, int64_t offset, bool rev) {
            if (rev) {
                return seq.substr(seq.size() - offset - include_terminal_positions,
                                  offset + include_terminal_positions);
            }
            else {
                return seq.substr(0, offset + include_terminal_positions);
            }
        };
        
        // record whether the positions are on the same node, and if so their relationship to each other
        colocation_t colocation;
        if (id(pos_1) == id(pos_2)) {
            if (is_rev(pos_1) == is_rev(pos_2)) {
                if (offset(pos_1) < offset(pos_2) + include_terminal_positions) {
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
        
        // for finding the largest node id in the subgraph
        id_t max_id = max(id(pos_1), id(pos_2));
        
        // a translator for node ids in g to node ids in the original graph
        unordered_map<id_t, id_t> id_trans;
        
        // the edges we have encountered in the traversal
        unordered_set<pair<pair<id_t, bool>, pair<id_t, bool>>> observed_edges;
        
        // the representation of the graph we're going to build up before storing in g (allows easier
        // subsetting operations than Graph, XG, or VG objects)
        unordered_map<id_t, LocalNode> graph;
        graph[id(pos_1)] = LocalNode(node_sequence(id(pos_1)));
        graph[id(pos_2)] = LocalNode(node_sequence(id(pos_2)));
        
        // keep track of whether we find a path or not
        bool found_target = false;
        
        unordered_set<pair<id_t, bool>> queued_traversals{make_pair(id(pos_1), is_rev(pos_1))};
        // mark final position as "queued" so that we won't look for additional traversals unless that's
        // the only way to find terminal cycles
        if (!(colocation == SharedNodeReverse && detect_terminal_cycles)) {
            queued_traversals.insert(make_pair(id(pos_2), is_rev(pos_2)));
        }
        // initialize the queue
        priority_queue<Traversal> queue;
        
        // the distance to the ends of the starting nodes
        int64_t first_traversal_length = graph[id(pos_1)].sequence.size() - offset(pos_1);
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
            found_target = (offset(pos_2) - offset(pos_1) <= max_len);
        }
        else {
            // search through graph to find the target, or to find cycles involving this node
            
#ifdef debug_vg_algorithms
            cerr << "FORWARD SEARCH: beginning search with forward max len " << forward_max_len << " and first traversal length " << first_traversal_length << endl;
#endif
            
            // if we can reach the end of this node, init the queue with it
            if (first_traversal_length <= forward_max_len) {
                queue.emplace(id(pos_1), is_rev(pos_1), first_traversal_length);
            }
            
            // search along a Dijkstra tree
            while (!queue.empty()) {
                // get the next closest node to the starting position
                Traversal trav = queue.top();
                queue.pop();
                
#ifdef debug_vg_algorithms
                cerr << "FORWARD SEARCH: traversing node " << trav.id << " in " << (trav.rev ? "reverse" : "forward") << " orientation at distance " << trav.dist << endl;
#endif
                
                // which side are we traversing out of?
                auto& edges_out = trav.rev ? graph[trav.id].edges_left : graph[trav.id].edges_right;
                for (Edge& edge : (trav.rev ? edges_on_start(trav.id) : edges_on_end(trav.id))) {
                    // get the orientation and id of the other side of the edge
#ifdef debug_vg_algorithms
                    cerr << "FORWARD SEARCH: got edge " << pb2json(edge) << endl;
#endif
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
                    
                    found_target = found_target || (next_id == id(pos_2) && next_rev == is_rev(pos_2));
                    max_id = max(max_id, next_id);
                    
                    // make sure the node is in
                    if (!graph.count(next_id)) {
                        graph[next_id] = LocalNode(node_sequence(next_id));
                    }
                    
                    // distance to the end of this node
                    int64_t dist_thru = trav.dist + graph[next_id].sequence.size();
                    if (!queued_traversals.count(make_pair(next_id, next_rev)) && dist_thru <= forward_max_len) {
                        // we can add more nodes along same path without going over the max length
                        // and we have not reached the target node yet
                        queue.emplace(next_id, next_rev, dist_thru);
                        queued_traversals.emplace(next_id, next_rev);
#ifdef debug_vg_algorithms
                        cerr << "FORWARD SEARCH: distance " << dist_thru << " is under maximum, adding to queue" << endl;
#endif
                    }
                    
                    bool reversing = (trav.rev != next_rev);
                    auto canonical_edge = canonical_form_edge(trav.id, trav.rev, next_id, reversing);
                    if (!observed_edges.count(canonical_edge)) {
                        // what side does this edge enter on the next node?
                        auto& edges_in = next_rev ? graph[next_id].edges_right : graph[next_id].edges_left;
                        // is the edge reversing?
                        // add this edge to the edge list on the current node
                        edges_out.push_back(make_pair(next_id, reversing));
                        // add to other node, but if it is a self-loop to the same side don't add it twice
                        if (!(trav.id == next_id && reversing) ) {
                            edges_in.push_back(make_pair(trav.id, reversing));
                        }
                        observed_edges.insert(canonical_edge);
                    }
                }
            }
        }
        
        // there is no path between the nodes under the maximum distance, leave g empty and return
        // an empty translator
        if (!found_target) {
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
            if (last_traversal_length <= backward_max_len) {
                queue.emplace(id(pos_2), !is_rev(pos_1), last_traversal_length);
            }
            
            // reset the queued traversal list and add the two reverse traversals
            queued_traversals.clear();
            queued_traversals.insert(make_pair(id(pos_2), !is_rev(pos_2)));
            queued_traversals.insert(make_pair(id(pos_1), !is_rev(pos_1)));
            
            // search along a Dijkstra tree
            while (!queue.empty()) {
                // get the next closest node to the starting position
                Traversal trav = queue.top();
                queue.pop();
                
#ifdef debug_vg_algorithms
                cerr << "BACKWARD SEARCH: traversing node " << trav.id << " in " << (trav.rev ? "reverse" : "forward") << " orientation at distance " << trav.dist << endl;
#endif
                
                for (Edge& edge : (trav.rev ? edges_on_start(trav.id) : edges_on_end(trav.id))) {
                    // get the orientation and id of the other side of the edge
#ifdef debug_vg_algorithms
                    cerr << "BACKWARD SEARCH: got edge " << pb2json(edge) << endl;
#endif
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
                    max_id = max(max_id, next_id);
                    
                    // make sure the node is in the graph
                    
                    if (!graph.count(next_id)) {
                        graph[next_id] = LocalNode(node_sequence(next_id));
                    }
                    
                    
                    // distance to the end of this node
                    int64_t dist_thru = trav.dist + graph[next_id].sequence.size();
                    if (!queued_traversals.count(make_pair(next_id, next_rev)) && dist_thru <= forward_max_len) {
                        // we can add more nodes along same path without going over the max length
                        // and we have not reached the target node yet
                        queue.emplace(next_id, next_rev, dist_thru);
                        queued_traversals.emplace(next_id, next_rev);
#ifdef debug_vg_algorithms
                        cerr << "BACKWARD SEARCH: distance " << dist_thru << " is under maximum, adding to queue" << endl;
#endif
                    }
                    
                    // is the edge reversing?
                    bool reversing = (trav.rev != next_rev);
                    auto canonical_edge = canonical_form_edge(trav.id, trav.rev, next_id, reversing);
                    if (!observed_edges.count(canonical_edge)) {
                        // which side are we traversing out of?
                        auto& edges_out = trav.rev ? graph[trav.id].edges_left : graph[trav.id].edges_right;
                        // what side does this edge enter on the next node?
                        auto& edges_in = next_rev ? graph[next_id].edges_right : graph[next_id].edges_left;
                        // add this edge to the edge list on the current node
                        edges_out.push_back(make_pair(next_id, reversing));
                        // add to other node, but if it is a self-loop to the same side don't add it twice
                        if (!(trav.id == next_id && reversing) ) {
                            edges_in.push_back(make_pair(trav.id, reversing));
                        }
                        observed_edges.insert(canonical_edge);
                    }
                }
            }
        }
        
#ifdef debug_vg_algorithms
        cerr << "state of graph after forward and backward search:" << endl;
        for (const pair<id_t, LocalNode>& node_record : graph) {
            cerr << node_record.first << " " << node_record.second.sequence << ": L( ";
            for (pair<id_t, bool> edge : node_record.second.edges_left) {
                cerr << edge.first << (edge.second ? "+" : "-") << " ";
            }
            cerr << " ) R( ";
            for (pair<id_t, bool> edge : node_record.second.edges_right) {
                cerr << edge.first << (edge.second ? "-" : "+") << " ";
            }
            cerr << ")" << endl;
        }
#endif
        
        // if we have to add new nodes, any id this large or larger will not have conflicts
        id_t next_id = max_id + 1;
        
        id_t duplicate_node_1 = 0, duplicate_node_2 = 0;
        
        // STEP 3: DUPLICATING NODES
        // if we're trying to detect terminal cycles, duplicate out the node so that the cyclic paths
        // survive the node cutting step
        if (detect_terminal_cycles) {
            // if there are edges traversed in both directions from the boundary position's nodes, then
            // they must be in cycles
            bool in_cycle_1 = !(graph[id(pos_1)].edges_left.empty() || graph[id(pos_1)].edges_right.empty());
            bool in_cycle_2 = !(graph[id(pos_2)].edges_left.empty() || graph[id(pos_2)].edges_right.empty());
            
            // logic changes depending on colocation of positions on same node
            switch (colocation) {
                case SeparateNodes:
                {
                    // the two positions are on separate nodes, so we can duplicate cycles independently
                    
                    if (in_cycle_1) {
                        LocalNode& node_1 = graph[id(pos_1)];
                        graph[next_id] = LocalNode(node_1.sequence);
                        LocalNode& new_node = graph[next_id];
                        
                        bool add_looping_connection = false;
                        
                        for (pair<id_t, bool>& edge : node_1.edges_right) {
                            if (edge.first == id(pos_1) && edge.second) {
                                // this is a reversing self loop, move it onto the cyclic node
                                edge.first = next_id;
                                // add a backwards reference to the original node
                                new_node.edges_right.emplace_back(id(pos_1), edge.second);
                                // add a copy of the edge onto the cyclic node
                                new_node.edges_right.emplace_back(next_id, edge.second);
                            }
                            else if (edge.first == id(pos_1)) {
                                // this is a non reversing self loop
                                
                                // mark that we need to make a connection between the old and new node, but
                                // don't add it yet so we don't mess up the iteration
                                add_looping_connection = true;
                                // make another edge for the nonreversing self loop on only the cyclic node
                                new_node.edges_right.emplace_back(next_id, edge.second);
                                new_node.edges_left.emplace_back(next_id, edge.second);
                            }
                            else {
                                // copy the edge and add a backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = edge.second ? next_node.edges_right : next_node.edges_left;
                                edges_backward.emplace_back(next_id, edge.second);
                                new_node.edges_right.emplace_back(edge.first, edge.second);
                            }
                        }
                        
                        for (pair<id_t, bool>& edge : node_1.edges_left) {
                            if (edge.first == id(pos_1) && edge.second) {
                                // this is a reversing self loop, move it onto the cyclic node
                                edge.first = next_id;
                                // add a backwards reference to the original node
                                new_node.edges_left.emplace_back(id(pos_1), edge.second);
                                // add a copy of the edge onto the cyclic node
                                new_node.edges_left.emplace_back(next_id, edge.second);
                            }
                            else if (edge.first != id(pos_1)) {
                                // copy the edge and add a backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = edge.second ? next_node.edges_left : next_node.edges_right;
                                edges_backward.emplace_back(next_id, edge.second);
                                new_node.edges_left.emplace_back(edge.first, edge.second);
                            }
                        }
                        
                        if (add_looping_connection) {
                            auto& new_incoming_edges = is_rev(pos_1) ? new_node.edges_right : new_node.edges_left;
                            auto& old_outgoing_edges = is_rev(pos_1) ? node_1.edges_left : node_1.edges_right;
                            new_incoming_edges.emplace_back(id(pos_1), false);
                            old_outgoing_edges.emplace_back(next_id, false);
                        }
                        
                        // record the translation
                        id_trans[next_id] = id(pos_1);
                        next_id++;
                    }
                    
                    if (in_cycle_2) {
                        LocalNode& node_2 = graph[id(pos_2)];
                        graph[next_id] = LocalNode(node_2.sequence);
                        LocalNode& new_node = graph[next_id];
                        
                        bool add_looping_connection = false;
                        
                        for (pair<id_t, bool>& edge : node_2.edges_right) {
                            if (edge.first == id(pos_2) && edge.second) {
                                // this is a reversing self loop, move it onto the cyclic node
                                edge.first = next_id;
                                // add a backwards reference to the original node
                                new_node.edges_right.emplace_back(id(pos_2), edge.second);
                                // add a copy of the edge onto the cyclic node
                                new_node.edges_right.emplace_back(next_id, edge.second);
                            }
                            else if (edge.first == id(pos_2)) {
                                // this is a non reversing self loop
                                
                                // mark that we need to make a connection between the old and new node, but
                                // don't add it yet so we don't mess up the iteration
                                add_looping_connection = true;
                                // make another edge for the nonreversing self loop on only the cyclic node
                                new_node.edges_right.emplace_back(next_id, edge.second);
                                new_node.edges_left.emplace_back(next_id, edge.second);
                            }
                            else {
                                // copy the edge and add a backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = edge.second ? next_node.edges_right : next_node.edges_left;
                                edges_backward.emplace_back(next_id, edge.second);
                                new_node.edges_right.emplace_back(edge.first, edge.second);
                            }
                        }
                        
                        for (pair<id_t, bool>& edge : node_2.edges_left) {
                            if (edge.first == id(pos_2) && edge.second) {
                                // this is a reversing self loop, move it onto the cyclic node
                                edge.first = next_id;
                                // add a backwards reference to the original node
                                new_node.edges_left.emplace_back(id(pos_1), edge.second);
                                // add a copy of the edge onto the cyclic node
                                new_node.edges_left.emplace_back(next_id, edge.second);
                            }
                            else if (edge.first != id(pos_2)) {
                                // copy the edge and add a backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = edge.second ? next_node.edges_left : next_node.edges_right;
                                edges_backward.emplace_back(next_id, edge.second);
                                new_node.edges_left.emplace_back(edge.first, edge.second);
                            }
                        }
                        
                        if (add_looping_connection) {
                            auto& new_outgoing_edges = is_rev(pos_2) ? new_node.edges_left : new_node.edges_right;
                            auto& old_incoming_edges = is_rev(pos_2) ? node_2.edges_right : node_2.edges_left;
                            new_outgoing_edges.emplace_back(id(pos_2), false);
                            old_incoming_edges.emplace_back(next_id, false);
                        }
                        
                        // record the translation
                        id_trans[next_id] = id(pos_2);
                        next_id++;
                    }
                    break;
                }
                case SharedNodeReachable:
                {
                    // one position is reachable from the next within the same node
                    
                    if (in_cycle_1) {
                        // later, we're going to trim this node to it's middle portion between the two positions
                        // so now that we want to preserve cycles, we need to make two new nodes that will hold
                        // the prefix and suffix of the node so that the edges have somewhere to attach to
                        
                        LocalNode& node = graph[id(pos_1)];
                        
                        // make a node for the righthand side of the traversal
                        id_t righthand_id = next_id;
                        graph[righthand_id] = LocalNode(trimmed_seq_right(node.sequence, offset(pos_1), is_rev(pos_1)));
                        LocalNode& righthand_node = graph[righthand_id];
                        
                        // move over the edges going out of the side that the traversal leaves
                        auto& righthand_new_edges = is_rev(pos_1) ? righthand_node.edges_left : righthand_node.edges_right;
                        righthand_new_edges = std::move(is_rev(pos_1) ? node.edges_left : node.edges_right);
                        
                        // update the edge references for the new node
                        for (pair<id_t, bool>& edge : righthand_new_edges) {
                            if (edge.first == id(pos_1) && edge.second) {
                                // if this is a reversing self loop, update it to the new node (the lefthand
                                // node hasn't been made yet, so let the ID on the edges pointing to it stay for
                                // the moment)
                                edge.first = righthand_id;
                            }
                            else {
                                // update the backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_right : next_node.edges_left;
                                auto edge_iter = std::find(edges_backward.begin(), edges_backward.end(),
                                                           make_pair(id(pos_1), edge.second));
                                (*edge_iter).first = righthand_id;
                            }
                        }
                        // record the translation
                        id_trans[righthand_id] = id(pos_1);
                        next_id++;
                        
                        // make a node for the lefthand side of the traversal
                        id_t lefthand_id = next_id;
                        graph[lefthand_id] = LocalNode(trimmed_seq_left(node.sequence, offset(pos_2), is_rev(pos_2)));
                        LocalNode& lefthand_node = graph[lefthand_id];
                        
                        // move over the edges going out of the side that the traversal leaves
                        auto& lefthand_new_edges = is_rev(pos_1) ? lefthand_node.edges_right : lefthand_node.edges_left;
                        lefthand_new_edges = std::move(is_rev(pos_1) ? node.edges_right : node.edges_left);
                        
                        // update the edge references for the new node
                        for (pair<id_t, bool>& edge : lefthand_new_edges) {
                            if (edge.first == id(pos_1)) {
                                // if this is a reversing self loop (non-reversing have already been updated
                                // to point to the righthand node), update it to the new node
                                edge.first = lefthand_id;
                            }
                            
                            if (!(edge.first == lefthand_id && edge.second)) {
                                // update the backward reference unless this is a reversing self loop
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_left : next_node.edges_right;
                                auto edge_iter = std::find(edges_backward.begin(), edges_backward.end(),
                                                           make_pair(id(pos_1), edge.second));
                                (*edge_iter).first = lefthand_id;
                            }
                        }
                        // record the translation
                        id_trans[lefthand_id] = id(pos_1);
                        next_id++;
                        
                        // now we have nodes to hold the edges, but we haven't preserved cycles that go through
                        // the node itself yet. to do that we need to duplicate it
                        
                        graph[next_id] = LocalNode(node.sequence);
                        LocalNode& cycle_node = graph[next_id];
                        
                        bool add_looping_connection = false;
                        
                        auto& edges_out = is_rev(pos_1) ? cycle_node.edges_left : cycle_node.edges_right;
                        auto& edges_in = is_rev(pos_1) ? cycle_node.edges_right : cycle_node.edges_left;
                        for (pair<id_t, bool>& edge : righthand_new_edges) {
                            if (edge.first == righthand_id) {
                                // this must be a reversing self loop, move it onto the cyclic node
                                edge.first = next_id;
                                // add a backwards reference to the righthand node
                                edges_out.emplace_back(righthand_id, edge.second);
                                // add a copy of the edge onto the cyclic node
                                edges_out.emplace_back(next_id, edge.second);
                            }
                            else if (edge.first == lefthand_id) {
                                // there is a non reversing self loop, mark it now but wait to add until later
                                // so we don't mess up the iteration
                                add_looping_connection = true;
                            }
                            else {
                                // copy the edge and add a backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_right : next_node.edges_left;
                                edges_backward.emplace_back(next_id, edge.second);
                                edges_out.emplace_back(edge.first, edge.second);
                            }
                        }
                        
                        for (pair<id_t, bool>& edge : lefthand_new_edges) {
                            if (edge.first == lefthand_id) {
                                // this must be a reversing self loop, move it onto the cyclic node
                                edge.first = next_id;
                                // add a backwards reference to the lefthand node
                                edges_in.emplace_back(lefthand_id, edge.second);
                                // add a copy of the edge onto the cyclic node
                                edges_in.emplace_back(next_id, edge.second);
                            }
                            else if (edge.first != righthand_id) {
                                // copy the edge and add a backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_left : next_node.edges_right;
                                edges_backward.emplace_back(next_id, edge.second);
                                edges_in.emplace_back(edge.first, edge.second);
                            }
                        }
                        
                        if (add_looping_connection) {
                            // add an edge from the righthand side to the cyclic node
                            righthand_new_edges.emplace_back(next_id, false);
                            edges_in.emplace_back(righthand_id, false);
                            
                            // add an edge from the lefthand side to the cyclic node
                            lefthand_new_edges.emplace_back(next_id, false);
                            edges_out.emplace_back(lefthand_id, false);
                            
                            // add an edge for the nonreversing self loop on only the cyclic node
                            edges_out.emplace_back(next_id, false);
                            edges_in.emplace_back(next_id, false);
                        }
                        
                        
                        id_trans[next_id] = id(pos_1);
                        next_id++;
                        
                        
                        duplicate_node_1 = righthand_id;
                        duplicate_node_2 = lefthand_id;
                    }
                    break;
                }
                case SharedNodeUnreachable:
                {
                    // all paths between these positions are cyclical, but we still duplicate the node
                    // so that any cycles that pass all the way through the node are there to be accepted
                    // or rejected by the distance filter
                    
                    LocalNode& node = graph[id(pos_1)];
                    
                    graph[next_id] = LocalNode(node.sequence);
                    LocalNode& cycle_node = graph[next_id];
                    
                    auto& new_outgoing_edges = is_rev(pos_1) ? cycle_node.edges_left : cycle_node.edges_right;
                    auto& new_incoming_edges = is_rev(pos_1) ? cycle_node.edges_right : cycle_node.edges_left;
                    
                    auto& old_outgoing_edges = is_rev(pos_1) ? node.edges_left : node.edges_right;
                    auto& old_incoming_edges = is_rev(pos_1) ? node.edges_right : node.edges_left;
                    
                    bool add_looping_connection = false;
                    
                    for (pair<id_t, bool>& edge : old_outgoing_edges) {
                        if (edge.first == id(pos_1) && edge.second) {
                            // this is a reversing self loop, move it onto the cyclic node
                            edge.first = next_id;
                            // add a backwards reference to the original node
                            new_outgoing_edges.emplace_back(id(pos_1), edge.second);
                            // add a copy of the edge onto the cyclic node
                            new_outgoing_edges.emplace_back(next_id, edge.second);
                        }
                        else if (edge.first == id(pos_1)) {
                            // this is a non reversing self loop, mark it but don't add
                            // the edge yet so we don't mess up the iteration
                            add_looping_connection = true;
                        }
                        else {
                            // copy the edge and add a backward reference
                            LocalNode& next_node = graph[edge.first];
                            auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_right : next_node.edges_left;
                            edges_backward.emplace_back(next_id, edge.second);
                            new_outgoing_edges.emplace_back(edge.first, edge.second);
                        }
                    }
                    
                    for (pair<id_t, bool>& edge : old_incoming_edges) {
                        if (edge.first == id(pos_1) && edge.second) {
                            // this is a reversing self loop, move it onto the cyclic node
                            edge.first = next_id;
                            // add a backwards reference to the original node
                            new_incoming_edges.emplace_back(id(pos_1), edge.second);
                            // add a copy of the edge onto the cyclic node
                            new_incoming_edges.emplace_back(next_id, edge.second);
                        }
                        else if (edge.first != id(pos_1)) {
                            // copy the edge and add a backward reference
                            LocalNode& next_node = graph[edge.first];
                            auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_left : next_node.edges_right;
                            edges_backward.emplace_back(next_id, edge.second);
                            new_incoming_edges.emplace_back(edge.first, edge.second);
                        }
                    }
                    
                    if (add_looping_connection) {
                        // add edge outward to new node
                        old_outgoing_edges.emplace_back(next_id, false);
                        new_incoming_edges.emplace_back(id(pos_1), false);
                        // edge edge inward from new node
                        old_incoming_edges.emplace_back(next_id, false);
                        new_outgoing_edges.emplace_back(id(pos_1), false);
                        // add cycle from new node to itself
                        new_outgoing_edges.emplace_back(next_id, false);
                        new_incoming_edges.emplace_back(next_id, false);
                    }
                    
                    id_trans[next_id] = id(pos_1);
                    next_id++;
                    
                    break;
                }
                case SharedNodeReverse:
                {
                    if (in_cycle_1) {
                        LocalNode& node = graph[id(pos_1)];
                        
                        graph[next_id] = LocalNode(node.sequence);
                        LocalNode& cycle_node = graph[next_id];
                        
                        auto& new_outgoing_edges = is_rev(pos_1) ? cycle_node.edges_left : cycle_node.edges_right;
                        auto& new_incoming_edges = is_rev(pos_1) ? cycle_node.edges_right : cycle_node.edges_left;
                        
                        auto& old_outgoing_edges = is_rev(pos_1) ? node.edges_left : node.edges_right;
                        auto& old_incoming_edges = is_rev(pos_1) ? node.edges_right : node.edges_left;
                        
                        bool add_reversing_connection = false;
                        bool add_looping_connection = false;
                        
                        for (pair<id_t, bool>& edge : old_outgoing_edges) {
                            if (edge.first == id(pos_1) && edge.second) {
                                // this is a reversing self loop
                                // indicate that we need to add a reversing edge between these
                                // but don't add it yet so we don't screw up iteration
                                add_reversing_connection = true;
                            }
                            else if (edge.first == id(pos_1)) {
                                // this is a non reversing self loop
                                // indicate that we need to add a reversing edge between these
                                // but don't add it yet so we don't screw up iteration
                                add_looping_connection = true;
                            }
                            else {
                                // copy the edge and add a backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_right : next_node.edges_left;
                                edges_backward.emplace_back(next_id, edge.second);
                                new_outgoing_edges.emplace_back(edge.first, edge.second);
                            }
                        }
                            
                        
                        for (pair<id_t, bool>& edge : old_incoming_edges) {
                            if (edge.first == id(pos_1) && edge.second) {
                                // this is a reversing self loop, add a copy of the edge onto the cyclic node
                                new_incoming_edges.emplace_back(next_id, edge.second);
                            }
                            else if (edge.first != id(pos_1)) {
                                // copy the edge and add a backward reference
                                LocalNode& next_node = graph[edge.first];
                                auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_left : next_node.edges_right;
                                edges_backward.emplace_back(next_id, edge.second);
                                new_incoming_edges.emplace_back(edge.first, edge.second);
                            }
                        }
                        
                        // preserve cycles involving a reversing self loop
                        if (add_reversing_connection) {
                            // add a connection to the new cyclic node
                            old_outgoing_edges.emplace_back(next_id, true);
                            new_outgoing_edges.emplace_back(id(pos_1), true);
                            // add a self loop on the cyclic node
                            new_outgoing_edges.emplace_back(next_id, true);
                        }
                        
                        if (add_looping_connection) {
                            // add a connection to the new cyclic node
                            old_outgoing_edges.emplace_back(next_id, false);
                            new_incoming_edges.emplace_back(id(pos_1), false);
                            // add a self loop on the cyclic node
                            new_outgoing_edges.emplace_back(next_id, false);
                            new_incoming_edges.emplace_back(next_id, false);
                        }
                        
                        id_trans[next_id] = id(pos_1);
                        next_id++;
                    }
                    
                    break;
                }
            }
        }
        
#ifdef debug_vg_algorithms
        cerr << "state of graph after duplicating nodes to preserve cycles:" << endl;
        for (const pair<id_t, LocalNode>& node_record : graph) {
            cerr << node_record.first << " " << node_record.second.sequence << ": L( ";
            for (pair<id_t, bool> edge : node_record.second.edges_left) {
                cerr << edge.first << (edge.second ? "+" : "-") << " ";
            }
            cerr << " ) R( ";
            for (pair<id_t, bool> edge : node_record.second.edges_right) {
                cerr << edge.first << (edge.second ? "-" : "+") << " ";
            }
            cerr << ")" << endl;
        }
#endif
        
        // STEP 4: CUTTING NODES
        // now cut the two end nodes at the designated positions and remove the edges on the cut side
        // to make the end positions tips in the graph
        
        switch (colocation) {
            case SeparateNodes:
            {
                LocalNode& node_1 = graph[id(pos_1)];
                LocalNode& node_2 = graph[id(pos_2)];
                auto& outward_edges_1 = is_rev(pos_1) ? node_1.edges_right : node_1.edges_left;
                auto& outward_edges_2 = is_rev(pos_2) ? node_2.edges_left : node_2.edges_right;
                // delete all backward edge references
                for (pair<id_t, bool>& edge : outward_edges_1) {
                    if (!(edge.first == id(pos_1) && edge.second)) {
                        LocalNode& next_node = graph[edge.first];
                        auto& backward_edges = is_rev(pos_1) != edge.second ? next_node.edges_left : next_node.edges_right;
                        backward_edges.erase(std::find(backward_edges.begin(), backward_edges.end(),
                                                       make_pair(id(pos_1), edge.second)));
                    }
                }
                for (pair<id_t, bool>& edge : outward_edges_2) {
                    if (!(edge.first == id(pos_2) && edge.second)) {
                        LocalNode& next_node = graph[edge.first];
                        auto& backward_edges = is_rev(pos_2) != edge.second ? next_node.edges_right : next_node.edges_left;
                        backward_edges.erase(std::find(backward_edges.begin(), backward_edges.end(),
                                                       make_pair(id(pos_2), edge.second)));
                    }
                }
                // clear the original edge lists
                outward_edges_1.clear();
                outward_edges_2.clear();
                
                // cut the node sequence
                node_1.sequence = trimmed_seq_right(node_1.sequence, offset(pos_1), is_rev(pos_1));
                node_2.sequence = trimmed_seq_left(node_2.sequence, offset(pos_2), is_rev(pos_2));
                break;
            }
            case SharedNodeReachable:
            {
                LocalNode& node = graph[id(pos_1)];
                // delete all backward edge references in both directions
                for (pair<id_t, bool>& edge : node.edges_right) {
                    if (!(edge.first == id(pos_1) && edge.second)) {
                        LocalNode& next_node = graph[edge.first];
                        auto& backward_edges = is_rev(pos_1) != edge.second ? next_node.edges_left : next_node.edges_right;
                        backward_edges.erase(std::find(backward_edges.begin(), backward_edges.end(),
                                                       make_pair(id(pos_1), edge.second)));
                    }
                }
                for (pair<id_t, bool>& edge : node.edges_left) {
                    if (!(edge.first == id(pos_2) && edge.second)) {
                        LocalNode& next_node = graph[edge.first];
                        auto& backward_edges = is_rev(pos_2) != edge.second ? next_node.edges_right : next_node.edges_left;
                        backward_edges.erase(std::find(backward_edges.begin(), backward_edges.end(),
                                                       make_pair(id(pos_2), edge.second)));
                    }
                }
                
                // clear the node's own edge lists
                node.edges_right.clear();
                node.edges_left.clear();
                
                // cut the node sequence
                if (is_rev(pos_1)) {
                    node.sequence = node.sequence.substr(node.sequence.size() - offset(pos_2) - include_terminal_positions,
                                                         offset(pos_2) - offset(pos_1) - 1 + 2 * include_terminal_positions);
                }
                else {
                    node.sequence = node.sequence.substr(offset(pos_1) + 1 - include_terminal_positions,
                                                         offset(pos_2) - offset(pos_1) - 1 + 2 * include_terminal_positions);
                }
                break;
            }
            case SharedNodeUnreachable:
            {
                LocalNode& node = graph[id(pos_1)];
                graph[next_id] = LocalNode(node.sequence);
                LocalNode& new_node = graph[next_id];
                
                // move the edges from one side onto the new node
                new_node.edges_right = std::move(node.edges_right);
                node.edges_right.clear();
                
                // relabel the edges pointing back into this side
                for (pair<id_t, bool>& edge : new_node.edges_right) {
                    LocalNode& next_node = graph[edge.first];
                    for (pair<id_t, bool>& edge_backward : edge.second ? next_node.edges_right : next_node.edges_left) {
                        if (edge_backward.first == id(pos_1)) {
                            edge_backward.first = next_id;
                            break;
                        }
                    }
                }
                
                // cut the sequences of the two nodes according to the search positions and switch
                // the pointer for one of the positions onto the new node
                if (is_rev(pos_1)) {
                    id_trans[next_id] = id(pos_2);
                    get_id(pos_2) = next_id;
                    node.sequence = trimmed_seq_right(node.sequence, offset(pos_1), is_rev(pos_1));
                    new_node.sequence = trimmed_seq_left(new_node.sequence, offset(pos_2), is_rev(pos_2));
                }
                else {
                    id_trans[next_id] = id(pos_1);
                    get_id(pos_1) = next_id;
                    new_node.sequence = trimmed_seq_right(new_node.sequence, offset(pos_1), is_rev(pos_1));
                    node.sequence = trimmed_seq_left(node.sequence, offset(pos_2), is_rev(pos_2));
                }
                
                next_id++;
                
                break;
            }
            case SharedNodeReverse:
            {
                LocalNode& node = graph[id(pos_1)];
                // delete all backward edge references
                auto& incoming_edges = is_rev(pos_1) ? node.edges_right : node.edges_left;
                for (pair<id_t, bool>& edge : incoming_edges) {
                    if (!(edge.first == id(pos_1) && edge.second)) {
                        LocalNode& next_node = graph[edge.first];
                        auto& backward_edges = is_rev(pos_1) != edge.second ? next_node.edges_left : next_node.edges_right;
                        backward_edges.erase(std::find(backward_edges.begin(), backward_edges.end(),
                                                       make_pair(id(pos_1), edge.second)));
                    }
                }
                incoming_edges.clear();
                
                // now make a new node to be the sink
                
                graph[next_id] = LocalNode(node.sequence);
                LocalNode& new_node = graph[next_id];
                
                auto& old_outgoing_edges = is_rev(pos_1) ? node.edges_left : node.edges_right;
                auto& new_outgoing_edges = is_rev(pos_1) ? new_node.edges_left : new_node.edges_right;
                
                for (pair<id_t, bool>& edge : old_outgoing_edges) {
                    if (edge.first == id(pos_1)) {
                        // since we clear all edges from the other side, this must be a reversing self loop
                        // change it to being a connecting edge between the two nodes
                        edge.first = next_id;
                        new_outgoing_edges.emplace_back(id(pos_1), edge.second);
                    }
                    else {
                        // copy the edge and add a backwards reference
                        LocalNode& next_node = graph[edge.first];
                        auto& backward_edges = is_rev(pos_1) != edge.second ? next_node.edges_right : next_node.edges_left;
                        new_outgoing_edges.emplace_back(edge.first, edge.second);
                        backward_edges.emplace_back(next_id, edge.second);
                    }
                }
                
                // record the node translation
                id_trans[next_id] = id(pos_1);
                // mark the duplicated node as the new sink
                get_id(pos_2) = next_id;
                next_id++;
                                
                // trim the sequences
                node.sequence = trimmed_seq_right(node.sequence, offset(pos_1), is_rev(pos_1));
                new_node.sequence = trimmed_seq_left(new_node.sequence, offset(pos_2), is_rev(pos_2));
                
                break;
            }
        }
        
#ifdef debug_vg_algorithms
        cerr << "state of graph after cutting nodes:" << endl;
        for (const pair<id_t, LocalNode>& node_record : graph) {
            cerr << node_record.first << " " << node_record.second.sequence << ": L( ";
            for (pair<id_t, bool> edge : node_record.second.edges_left) {
                cerr << edge.first << (edge.second ? "+" : "-") << " ";
            }
            cerr << " ) R( ";
            for (pair<id_t, bool> edge : node_record.second.edges_right) {
                cerr << edge.first << (edge.second ? "-" : "+") << " ";
            }
            cerr << ")" << endl;
        }
#endif
        
        // STEP 5: PRUNING
        // the graph now contains all the paths we've indicated and the end positions are tips, we now
        // provide three options for pruning away any unnecessary nodes and edges we've added in the
        // process of searching for the subgraph that has this guarantee
        
        if (strict_max_len) {
            // OPTION 1: PRUNE TO PATHS UNDER MAX LENGTH
            // some nodes in the current graph may not be on paths, or the paths that they are on may be
            // above the maximum distance, so we do a forward-backward distance search to check
            
            unordered_map<pair<id_t, bool>, int64_t> forward_trav_dist;
            unordered_map<pair<id_t, bool>, int64_t> reverse_trav_dist;
            
            // re-initialize the queue in the forward direction
            queue.emplace(id(pos_1), is_rev(pos_1), graph[id(pos_1)].sequence.size());
            
            // reset the queued traversal list and the first traversal
            queued_traversals.clear();
            queued_traversals.insert(make_pair(id(pos_1), is_rev(pos_1)));
            
            // if we duplicated the start node, add that too
            if (duplicate_node_1) {
                queue.emplace(duplicate_node_1, is_rev(pos_1), graph[duplicate_node_1].sequence.size());
                queued_traversals.insert(make_pair(duplicate_node_1, is_rev(pos_1)));
            }
            
            while (!queue.empty()) {
                // get the next closest node traversal
                Traversal trav = queue.top();
                queue.pop();
                forward_trav_dist[make_pair(trav.id, trav.rev)] = trav.dist;
                
#ifdef debug_vg_algorithms
                cerr << "FORWARD PRUNE: traversing node " << trav.id << " in " << (trav.rev ? "reverse" : "forward") << " orientation at distance " << trav.dist << endl;
#endif
                
                // the edges in the direction of this traversal
                auto& edges_out = trav.rev ? graph[trav.id].edges_left : graph[trav.id].edges_right;
                for (const pair<id_t, bool>& edge : edges_out) {
                    
                    // the distance to the opposite side of this next node
                    int64_t dist_thru = trav.dist + graph[edge.first].sequence.size();
                    
                    // queue up the node traversal if it hasn't been seen before
                    pair<id_t, bool> next_trav = make_pair(edge.first, edge.second != trav.rev);
                    if (!queued_traversals.count(next_trav)) {
                        queue.emplace(next_trav.first, next_trav.second, dist_thru);
                        queued_traversals.insert(next_trav);
                    }
                }
            }
            
            // re-initialize the queue
            queue.emplace(id(pos_2), !is_rev(pos_2), 0);
            
            // reset the queued traversal list and add the reverse traversal
            queued_traversals.clear();
            queued_traversals.insert(make_pair(id(pos_2), !is_rev(pos_2)));
            
            // if we duplicated the end node, add that too
            if (duplicate_node_2) {
                queue.emplace(duplicate_node_2, !is_rev(pos_2), 0);
                queued_traversals.insert(make_pair(duplicate_node_2, !is_rev(pos_2)));
            }
            
            while (!queue.empty()) {
                // get the next closest node traversal
                Traversal trav = queue.top();
                queue.pop();
                reverse_trav_dist[make_pair(trav.id, trav.rev)] = trav.dist;
                
#ifdef debug_vg_algorithms
                cerr << "BACKWARD PRUNE: traversing node " << trav.id << " in " << (trav.rev ? "reverse" : "forward") << " orientation at distance " << trav.dist << endl;
#endif
                
                // the distance to the opposite side of this next node
                int64_t dist_thru = trav.dist + graph[trav.id].sequence.size();
                
                // the edges in the direction of this traversal
                auto& edges_out = trav.rev ? graph[trav.id].edges_left : graph[trav.id].edges_right;
                
                for (const pair<id_t, bool>& edge : edges_out) {
                    // queue up the node traversal if it hasn't been seen before
                    pair<id_t, bool> next_trav = make_pair(edge.first, edge.second != trav.rev);
                    if (!queued_traversals.count(next_trav)) {
                        queue.emplace(next_trav.first, next_trav.second, dist_thru);
                        queued_traversals.insert(next_trav);
                    }
                }
            }
            
            // now we have the lengths of the shortest path remaining in graph to and from each node
            // with these, we can compute the shortest path that uses each node and edge to see if it
            // should be included in the final graph
            
            vector<unordered_map<id_t, LocalNode>::iterator> to_erase;
            for (auto iter = graph.begin(); iter != graph.end(); iter++) {
                bool erase_node = true;
                id_t node_id = (*iter).first;
                // did a short enough path use one or the other traversal directions?
                if (forward_trav_dist.count(make_pair(node_id, true)) &&
                    reverse_trav_dist.count(make_pair(node_id, false))) {
                    if (forward_trav_dist[make_pair(node_id, true)]
                        + reverse_trav_dist[make_pair(node_id, false)] <= max_len) {
                        erase_node = false;
                    }
                }
                if (forward_trav_dist.count(make_pair(node_id, false)) &&
                    reverse_trav_dist.count(make_pair(node_id, true))) {
                    if (forward_trav_dist[make_pair(node_id, false)]
                        + reverse_trav_dist[make_pair(node_id, true)] <= max_len) {
                        erase_node = false;
                    }
                }
                
                if (erase_node) {
                    // the shortest path using this node is too long
                    to_erase.push_back(iter);
                }
                else {
                    LocalNode& node = (*iter).second;
                    // find which edges are traversed on sufficiently short paths
                    auto new_right_end = std::remove_if(node.edges_right.begin(), node.edges_right.end(),
                                                        [&](const pair<id_t, bool>& edge) {
                                                            bool erase_edge = true;
                                                            if (forward_trav_dist.count(make_pair(node_id, false))
                                                                && reverse_trav_dist.count(make_pair(edge.first, !edge.second))) {
                                                                if (forward_trav_dist[make_pair(node_id, false)]
                                                                    + reverse_trav_dist[make_pair(edge.first, !edge.second)]
                                                                    + graph[edge.first].sequence.size() <= max_len) {
                                                                    erase_edge = false;
                                                                }
                                                            }
                                                            if (forward_trav_dist.count(make_pair(edge.first, !edge.second))
                                                                && reverse_trav_dist.count(make_pair(node_id, false))) {
                                                                if (forward_trav_dist[make_pair(edge.first, !edge.second)]
                                                                    + reverse_trav_dist[make_pair(node_id, false)]
                                                                    + graph[node_id].sequence.size() <= max_len) {
                                                                    erase_edge = false;
                                                                }
                                                            }
                                                            return erase_edge;
                                                        });
                    auto new_left_end = std::remove_if(node.edges_left.begin(), node.edges_left.end(),
                                                       [&](const pair<id_t, bool>& edge) {
                                                           bool erase_edge = true;
                                                           if (forward_trav_dist.count(make_pair(node_id, true))
                                                               && reverse_trav_dist.count(make_pair(edge.first, edge.second))) {
                                                               if (forward_trav_dist[make_pair(node_id, true)]
                                                                   + reverse_trav_dist[make_pair(edge.first, edge.second)]
                                                                   + graph[edge.first].sequence.size() <= max_len) {
                                                                   erase_edge = false;
                                                               }
                                                           }
                                                           if (forward_trav_dist.count(make_pair(edge.first, edge.second))
                                                               && reverse_trav_dist.count(make_pair(node_id, true))) {
                                                               if (forward_trav_dist[make_pair(edge.first, edge.second)]
                                                                   + reverse_trav_dist[make_pair(node_id, true)]
                                                                   + graph[node_id].sequence.size() <= max_len) {
                                                                   erase_edge = false;
                                                               }
                                                           }
                                                           return erase_edge;
                                                       });
                    // remove the edges that only occurred on path that were too long
                    node.edges_right.resize(new_right_end - node.edges_right.begin());
                    node.edges_left.resize(new_left_end - node.edges_left.begin());
                }
            }
            
            // remove the nodes
            for (auto& iter : to_erase) {
                // if we're removing one of the duplicated nodes, remove it from the ID translator
                if (id_trans.count((*iter).first)) {
                    id_trans.erase((*iter).first);
                }
                graph.erase(iter);
            }
        }
        else if (only_paths) {
            // OPTION 2: PRUNE TO PATHS
            // some nodes in the current graph may not be on paths, so we do a forward-backward
            // reachability search to check
            
            list<pair<id_t, bool>> stack;
            
            unordered_set<pair<id_t, bool>> forward_reachable;
            unordered_set<pair<id_t, bool>> reverse_reachable;
            
            // initialize the stack in the forward direction
            stack.emplace_back(id(pos_1), is_rev(pos_1));
            forward_reachable.emplace(id(pos_1), is_rev(pos_1));
            
            // if we duplicated the start node, add that too
            if (duplicate_node_1) {
                stack.emplace_back(duplicate_node_1, is_rev(pos_1));
                forward_reachable.emplace(duplicate_node_1, is_rev(pos_1));
            }
            
            while (!stack.empty()) {
                // get the next closest node traversal
                pair<id_t, bool> trav = stack.back();
                stack.pop_back();
                
                // the edges in the direction of this traversal
                auto& edges_out = trav.second ? graph[trav.first].edges_left : graph[trav.first].edges_right;
                for (const pair<id_t, bool>& edge : edges_out) {
                    
                    // queue up the node traversal if it hasn't been seen before
                    pair<id_t, bool> next_trav = make_pair(edge.first, edge.second != trav.second);
                    if (!forward_reachable.count(next_trav)) {
                        stack.emplace_back(next_trav);
                        forward_reachable.insert(next_trav);
                    }
                }
            }
            
            // re-initialize the stack in the reverse direction
            stack.emplace_back(id(pos_2), !is_rev(pos_2));
            reverse_reachable.emplace(id(pos_2), !is_rev(pos_2));
            
            // if we duplicated the second end node, add that too
            if (duplicate_node_2) {
                stack.emplace_back(duplicate_node_2, !is_rev(pos_2));
                reverse_reachable.emplace(duplicate_node_2, !is_rev(pos_2));
            }
            
            while (!stack.empty()) {
                // get the next closest node traversal
                pair<id_t, bool> trav = stack.back();
                stack.pop_back();
                
                // the edges in the direction of this traversal
                auto& edges_out = trav.second ? graph[trav.first].edges_left : graph[trav.first].edges_right;
                for (const pair<id_t, bool>& edge : edges_out) {
                    
                    // queue up the node traversal if it hasn't been seen before
                    pair<id_t, bool> next_trav = make_pair(edge.first, edge.second != trav.second);
                    if (!reverse_reachable.count(next_trav)) {
                        stack.emplace_back(next_trav);
                        reverse_reachable.insert(next_trav);
                    }
                }
            }
            
            // now we know which nodes are reachable from both ends, to be on a path between the end positions,
            // a node or edge must be reachable from both directions
            
            vector<unordered_map<id_t, LocalNode>::iterator> to_erase;
            for (auto iter = graph.begin(); iter != graph.end(); iter++) {
                id_t node_id = (*iter).first;
                // did a path use one or the other traversal directions?
                if (!(forward_reachable.count(make_pair(node_id, true)) &&
                      reverse_reachable.count(make_pair(node_id, false))) &&
                    !(forward_reachable.count(make_pair(node_id, false)) &&
                      reverse_reachable.count(make_pair(node_id, true)))) {
                        
                    to_erase.push_back(iter);
                }
                else {
                    LocalNode& node = (*iter).second;
                    // find which edges are also on traversed paths
                    auto new_right_end = std::remove_if(node.edges_right.begin(), node.edges_right.end(),
                                                        [&](const pair<id_t, bool>& edge) {
                                                            return !(forward_reachable.count(make_pair(node_id, false)) &&
                                                                     reverse_reachable.count(make_pair(edge.first, !edge.second))) &&
                                                                   !(forward_reachable.count(make_pair(edge.first, !edge.second))
                                                                     && reverse_reachable.count(make_pair(node_id, false)));
                                                        });
                    auto new_left_end = std::remove_if(node.edges_left.begin(), node.edges_left.end(),
                                                       [&](const pair<id_t, bool>& edge) {
                                                           return !(forward_reachable.count(make_pair(node_id, true)) &&
                                                                    reverse_reachable.count(make_pair(edge.first, edge.second))) &&
                                                                  !(forward_reachable.count(make_pair(edge.first, edge.second))
                                                                    && reverse_reachable.count(make_pair(node_id, true)));
                                                       });
                    // remove the edges that only occurred on path that were too long
                    node.edges_right.resize(new_right_end - node.edges_right.begin());
                    node.edges_left.resize(new_left_end - node.edges_left.begin());
                }
            }
            
            // remove the nodes
            for (auto& iter : to_erase) {
                // if we're removing one of the duplicated nodes, remove it from the ID translator
                if (id_trans.count((*iter).first)) {
                    id_trans.erase((*iter).first);
                }
                graph.erase(iter);
            }
        }
        else if (no_additional_tips) {
            // OPTION 3: PRUNE ADDITIONAL TIPS
            // all cycles to the original (non-duplicated) nodes are dangling tips. there may also be tips that
            // resulted from paths we explored until hitting the max search length in the graph extraction step.
            // next we remove all tips (except if the tip is a node with our end position on it)
            
            if (no_additional_tips) {
                unordered_map<id_t, int64_t> left_degree;
                unordered_map<id_t, int64_t> right_degree;
                
                for (const pair<id_t, LocalNode>& node_record : graph) {
                    left_degree[node_record.first] = node_record.second.edges_left.size();
                    right_degree[node_record.first] = node_record.second.edges_right.size();
                }
                
                // remove nodes from the graph if they are tips or only connect to tips
                list<id_t> to_check;
                for (const pair<id_t, int64_t>& left_deg_record : left_degree) {
                    // check every node in the graph once
                    to_check.push_front(left_deg_record.first);
#ifdef debug_vg_algorithms
                    cerr << "TIP REMOVAL: initializing queue with node " << left_deg_record.first << endl;
#endif
                    while (!to_check.empty()) {
                        id_t node_id = to_check.back();
                        to_check.pop_back();
                        if (node_id == id(pos_1) || node_id == id(pos_2) ||
                            node_id == duplicate_node_1 || node_id == duplicate_node_2 ||
                            !graph.count(node_id)) {
                            // the end nodes get a free pass on begin tips, or we may
                            // have already pruned this node
                            continue;
                        }
                        if (left_degree[node_id] == 0) {
#ifdef debug_vg_algorithms
                            cerr << "TIP REMOVAL: node " << node_id << " is a left tip" << endl;
#endif
                            if (id_trans.count(node_id)) {
                                id_trans.erase(node_id);
                            }
                            for (pair<id_t, bool>& edge : graph[node_id].edges_right) {
                                if (edge.second) {
                                    right_degree[edge.first]--;
                                }
                                else {
                                    left_degree[edge.first]--;
                                }
                                to_check.push_front(edge.first);
                            }
                            graph.erase(node_id);
                        }
                        else if (right_degree[node_id] == 0) {
#ifdef debug_vg_algorithms
                            cerr << "TIP REMOVAL: node " << node_id << " is a right tip" << endl;
#endif
                            if (id_trans.count(node_id)) {
                                id_trans.erase(node_id);
                            }
                            for (pair<id_t, bool>& edge : graph[node_id].edges_left) {
                                if (edge.second) {
                                    left_degree[edge.first]--;
                                }
                                else {
                                    right_degree[edge.first]--;
                                }
                                to_check.push_front(edge.first);
                            }
                            graph.erase(node_id);
                        }
                    }
                }
                
                // remove edges that point to nodes we removed
                for (auto&& node_record : graph) {
                    LocalNode& node = node_record.second;
                    auto new_left_end = std::remove_if(node.edges_left.begin(), node.edges_left.end(),
                                                       [&](const pair<id_t, bool>& edge) {return !graph.count(edge.first);});
                    auto new_right_end = std::remove_if(node.edges_right.begin(), node.edges_right.end(),
                                                        [&](const pair<id_t, bool>& edge) {return !graph.count(edge.first);});
                    node.edges_left.resize(new_left_end - node.edges_left.begin());
                    node.edges_right.resize(new_right_end - node.edges_right.begin());
                    
                }
            }
        }
        
            
#ifdef debug_vg_algorithms
        cerr << "state of graph after pruning:" << endl;
        for (const pair<id_t, LocalNode>& node_record : graph) {
            cerr << node_record.first << " " << node_record.second.sequence << ": L( ";
            for (pair<id_t, bool> edge : node_record.second.edges_left) {
                cerr << edge.first << (edge.second ? "+" : "-") << " ";
            }
            cerr << " ) R( ";
            for (pair<id_t, bool> edge : node_record.second.edges_right) {
                cerr << edge.first << (edge.second ? "-" : "+") << " ";
            }
            cerr << ")" << endl;
        }
#endif
        
        // STEP 6: TRANSLATION TO PROTOBUF
        // transfer the local graph we've been building to g
        
        // add all remaining nodes that do not have recorded translations to the ID translator
        for (auto& node_record : graph) {
            if (!id_trans.count(node_record.first)) {
                id_trans[node_record.first] = node_record.first;
            }
        }
        
        for (const pair<id_t, LocalNode>& node_record : graph) {
            // add in each node
            Node* node = g.add_node();
            node->set_id(node_record.first);
            node->set_sequence(node_record.second.sequence);
            
            // add each incoming edge
            for (const pair<id_t, bool>& edge : node_record.second.edges_left) {
                // break symmetry on the edge to avoid adding it from both edge lists
                if (edge.first > node_record.first || (edge.first == node_record.first && edge.second)) {
                    Edge* pb_edge = g.add_edge();
                    pb_edge->set_from(node_record.first);
                    pb_edge->set_to(edge.first);
                    pb_edge->set_from_start(true);
                    pb_edge->set_to_end(!edge.second);
                }
            }
            for (const pair<id_t, bool>& edge : node_record.second.edges_right) {
                // break symmetry on the edge to avoid adding it from both edge lists
                if (edge.first >= node_record.first) {
                    Edge* pb_edge = g.add_edge();
                    pb_edge->set_from(node_record.first);
                    pb_edge->set_to(edge.first);
                    pb_edge->set_from_start(false);
                    pb_edge->set_to_end(edge.second);
                }
            }
        }
        
        // TODO: it's not enough to return the translator because there's also the issue of the positions
        // on the first node being offset (however this information is fully contained in the arguments of
        // the function, which are obviously available in the environment that calls it)
        return id_trans;
    }
    
    // wrapper for connecting graph algorithm using VG to access graph elements
    unordered_map<id_t, id_t> extract_connecting_graph(VG& vg, Graph& g, int64_t max_len,
                                                       pos_t pos_1, pos_t pos_2,
                                                       bool include_terminal_positions,
                                                       bool detect_terminal_cycles,
                                                       bool no_additional_tips,
                                                       bool only_paths,
                                                       bool strict_max_len){
        
        return extract_connecting_graph_internal(g, max_len, pos_1, pos_2,
                                                 include_terminal_positions,
                                                 detect_terminal_cycles,
                                                 no_additional_tips,
                                                 only_paths,
                                                 strict_max_len,
                                                 [&](id_t id) {
                                                     vector<Edge> to_return;
                                                     for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                         if ((edge->from() == id && edge->from_start())
                                                             || (edge->to() == id && !edge->to_end())) {
                                                             to_return.push_back(*edge);
                                                         }
                                                     }
                                                     return to_return;
                                                 },
                                                 [&](id_t id) {
                                                     vector<Edge> to_return;
                                                     for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                         if ((edge->from() == id && !edge->from_start())
                                                             || (edge->to() == id && edge->to_end())) {
                                                             to_return.push_back(*edge);
                                                         }
                                                     }
                                                     return to_return;
                                                 },
                                                 [&](id_t id) {
                                                     return vg.get_node(id)->sequence();
                                                 });
    }
    
    // wrapper for connecting graph algorithm using XG to access graph elements
    unordered_map<id_t, id_t> extract_connecting_graph(xg::XG& xg_index, Graph& g, int64_t max_len,
                                                       pos_t pos_1, pos_t pos_2,
                                                       bool include_terminal_positions,
                                                       bool detect_terminal_cycles,
                                                       bool no_additional_tips,
                                                       bool only_paths,
                                                       bool strict_max_len,
                                                       LRUCache<id_t, Node>* node_cache,
                                                       LRUCache<id_t, vector<Edge>>* edge_cache) {
        if (node_cache && edge_cache) {
            return extract_connecting_graph_internal(g, max_len, pos_1, pos_2,
                                                     include_terminal_positions,
                                                     detect_terminal_cycles,
                                                     no_additional_tips,
                                                     only_paths,
                                                     strict_max_len,
                                                     [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (node_cache) {
            return extract_connecting_graph_internal(g, max_len, pos_1, pos_2,
                                                     include_terminal_positions,
                                                     detect_terminal_cycles,
                                                     no_additional_tips,
                                                     only_paths,
                                                     strict_max_len,
                                                     [&](id_t id) {return xg_index.edges_on_start(id);},
                                                     [&](id_t id) {return xg_index.edges_on_end(id);},
                                                     [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (edge_cache) {
            return extract_connecting_graph_internal(g, max_len, pos_1, pos_2,
                                                     include_terminal_positions,
                                                     detect_terminal_cycles,
                                                     no_additional_tips,
                                                     only_paths,
                                                     strict_max_len,
                                                     [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_index.node_sequence(id);});
        }
        else {
            return extract_connecting_graph_internal(g, max_len, pos_1, pos_2,
                                                     include_terminal_positions,
                                                     detect_terminal_cycles,
                                                     no_additional_tips,
                                                     only_paths,
                                                     strict_max_len,
                                                     [&](id_t id) {return xg_index.edges_on_start(id);},
                                                     [&](id_t id) {return xg_index.edges_on_end(id);},
                                                     [&](id_t id) {return xg_index.node_sequence(id);});
        }
    }
    
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
    
    void extract_containing_graph(VG& vg, Graph& g, const vector<pos_t>& positions, size_t max_dist) {
        
        // make a dummy vector for all positions at the same distance
        vector<size_t> dists(positions.size(), max_dist);
        return extract_containing_graph_internal(g, positions, dists, dists,
                                                 [&](id_t id) {
                                                     vector<Edge> to_return;
                                                     for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                         if ((edge->from() == id && edge->from_start())
                                                             || (edge->to() == id && !edge->to_end())) {
                                                             to_return.push_back(*edge);
                                                         }
                                                     }
                                                     return to_return;
                                                 },
                                                 [&](id_t id) {
                                                     vector<Edge> to_return;
                                                     for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                         if ((edge->from() == id && !edge->from_start())
                                                             || (edge->to() == id && edge->to_end())) {
                                                             to_return.push_back(*edge);
                                                         }
                                                     }
                                                     return to_return;
                                                 },
                                                 [&](id_t id) {
                                                     return vg.get_node(id)->sequence();
                                                 });
    }
    
    void extract_containing_graph(VG& vg, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_max_dist) {
        
        return extract_containing_graph_internal(g, positions, position_max_dist, position_max_dist,
                                                 [&](id_t id) {
                                                     vector<Edge> to_return;
                                                     for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                         if ((edge->from() == id && edge->from_start())
                                                             || (edge->to() == id && !edge->to_end())) {
                                                             to_return.push_back(*edge);
                                                         }
                                                     }
                                                     return to_return;
                                                 },
                                                 [&](id_t id) {
                                                     vector<Edge> to_return;
                                                     for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                         if ((edge->from() == id && !edge->from_start())
                                                             || (edge->to() == id && edge->to_end())) {
                                                             to_return.push_back(*edge);
                                                         }
                                                     }
                                                     return to_return;
                                                 },
                                                 [&](id_t id) {
                                                     return vg.get_node(id)->sequence();
                                                 });
    }
    
    void extract_containing_graph(VG& vg, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_forward_max_dist,
                                  const vector<size_t>& position_backward_max_dist) {
        
        return extract_containing_graph_internal(g, positions, position_forward_max_dist, position_backward_max_dist,
                                                 [&](id_t id) {
                                                     vector<Edge> to_return;
                                                     for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                         if ((edge->from() == id && edge->from_start())
                                                             || (edge->to() == id && !edge->to_end())) {
                                                             to_return.push_back(*edge);
                                                         }
                                                     }
                                                     return to_return;
                                                 },
                                                 [&](id_t id) {
                                                     vector<Edge> to_return;
                                                     for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                         if ((edge->from() == id && !edge->from_start())
                                                             || (edge->to() == id && edge->to_end())) {
                                                             to_return.push_back(*edge);
                                                         }
                                                     }
                                                     return to_return;
                                                 },
                                                 [&](id_t id) {
                                                     return vg.get_node(id)->sequence();
                                                 });
    }
    
    void extract_containing_graph(xg::XG& xg_index, Graph& g, const vector<pos_t>& positions, size_t max_dist,
                                  LRUCache<id_t, Node>* node_cache, LRUCache<id_t, vector<Edge>>* edge_cache) {
        
        // make a dummy vector for all positions at the same distance
        vector<size_t> dists(positions.size(), max_dist);
        if (node_cache && edge_cache) {
            return extract_containing_graph_internal(g, positions, dists, dists,
                                                     [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (node_cache) {
            return extract_containing_graph_internal(g, positions, dists, dists,
                                                     [&](id_t id) {return xg_index.edges_on_start(id);},
                                                     [&](id_t id) {return xg_index.edges_on_end(id);},
                                                     [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (edge_cache) {
            return extract_containing_graph_internal(g, positions, dists, dists,
                                                     [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_index.node_sequence(id);});
        }
        else {
            return extract_containing_graph_internal(g, positions, dists, dists,
                                                     [&](id_t id) {return xg_index.edges_on_start(id);},
                                                     [&](id_t id) {return xg_index.edges_on_end(id);},
                                                     [&](id_t id) {return xg_index.node_sequence(id);});
        }
    }
    
    void extract_containing_graph(xg::XG& xg_index, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_max_dist,
                                  LRUCache<id_t, Node>* node_cache, LRUCache<id_t, vector<Edge>>* edge_cache) {
        
        if (node_cache && edge_cache) {
            return extract_containing_graph_internal(g, positions, position_max_dist, position_max_dist,
                                                     [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (node_cache) {
            return extract_containing_graph_internal(g, positions, position_max_dist, position_max_dist,
                                                     [&](id_t id) {return xg_index.edges_on_start(id);},
                                                     [&](id_t id) {return xg_index.edges_on_end(id);},
                                                     [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (edge_cache) {
            return extract_containing_graph_internal(g, positions, position_max_dist, position_max_dist,
                                                     [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_index.node_sequence(id);});
        }
        else {
            return extract_containing_graph_internal(g, positions, position_max_dist, position_max_dist,
                                                     [&](id_t id) {return xg_index.edges_on_start(id);},
                                                     [&](id_t id) {return xg_index.edges_on_end(id);},
                                                     [&](id_t id) {return xg_index.node_sequence(id);});
        }
    }
    
    void extract_containing_graph(xg::XG& xg_index, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_forward_max_dist,
                                  const vector<size_t>& position_backward_max_dist,
                                  LRUCache<id_t, Node>* node_cache, LRUCache<id_t, vector<Edge>>* edge_cache) {
        
        if (node_cache && edge_cache) {
            return extract_containing_graph_internal(g, positions, position_forward_max_dist, position_backward_max_dist,
                                                     [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (node_cache) {
            return extract_containing_graph_internal(g, positions, position_forward_max_dist, position_backward_max_dist,
                                                     [&](id_t id) {return xg_index.edges_on_start(id);},
                                                     [&](id_t id) {return xg_index.edges_on_end(id);},
                                                     [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (edge_cache) {
            return extract_containing_graph_internal(g, positions, position_forward_max_dist, position_backward_max_dist,
                                                     [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                     [&](id_t id) {return xg_index.node_sequence(id);});
        }
        else {
            return extract_containing_graph_internal(g, positions, position_forward_max_dist, position_backward_max_dist,
                                                     [&](id_t id) {return xg_index.edges_on_start(id);},
                                                     [&](id_t id) {return xg_index.edges_on_end(id);},
                                                     [&](id_t id) {return xg_index.node_sequence(id);});
        }
    }
    
    unordered_map<id_t, id_t> extract_extending_graph_internal(Graph& g, int64_t max_dist, pos_t pos, bool backward,
                                                               bool preserve_cycles_on_src,
                                                               function<vector<Edge>(id_t)> edges_on_start,
                                                               function<vector<Edge>(id_t)> edges_on_end,
                                                               function<string(id_t)> node_sequence) {
        
        if (g.node_size() || g.edge_size()) {
            cerr << "error:[extract_extending_graph] must extract into an empty graph" << endl;
            assert(false);
        }
        
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
            }
            else {
                // the edge doesn't touch the source node, so it will certainly survive the cut
                add_edge = true;
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
    
    unordered_map<id_t, id_t> extract_extending_graph(VG& vg, Graph& g, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src){
        return extract_extending_graph_internal(g, max_dist, pos, backward, preserve_cycles_on_src,
                                                [&](id_t id) {
                                                    vector<Edge> to_return;
                                                    for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                        if ((edge->from() == id && edge->from_start())
                                                            || (edge->to() == id && !edge->to_end())) {
                                                            to_return.push_back(*edge);
                                                        }
                                                    }
                                                    return to_return;
                                                },
                                                [&](id_t id) {
                                                    vector<Edge> to_return;
                                                    for (Edge* edge : vg.edges_of(vg.get_node(id))) {
                                                        if ((edge->from() == id && !edge->from_start())
                                                            || (edge->to() == id && edge->to_end())) {
                                                            to_return.push_back(*edge);
                                                        }
                                                    }
                                                    return to_return;
                                                },
                                                [&](id_t id) {
                                                    return vg.get_node(id)->sequence();
                                                });
    }
    
    unordered_map<id_t, id_t> extract_extending_graph(xg::XG& xg_index, Graph& g, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src,
                                                      LRUCache<id_t, Node>* node_cache,
                                                      LRUCache<id_t, vector<Edge>>* edge_cache) {
        
        if (node_cache && edge_cache) {
            return extract_extending_graph_internal(g, max_dist, pos, backward, preserve_cycles_on_src,
                                                    [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                    [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                    [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        if (node_cache) {
            return extract_extending_graph_internal(g, max_dist, pos, backward, preserve_cycles_on_src,
                                                    [&](id_t id) {return xg_index.edges_on_start(id);},
                                                    [&](id_t id) {return xg_index.edges_on_end(id);},
                                                    [&](id_t id) {return xg_cached_node_sequence(id, &xg_index, *node_cache);});
        }
        else if (edge_cache) {
            return extract_extending_graph_internal(g, max_dist, pos, backward, preserve_cycles_on_src,
                                                    [&](id_t id) {return xg_cached_edges_on_start(id, &xg_index, *edge_cache);},
                                                    [&](id_t id) {return xg_cached_edges_on_end(id, &xg_index, *edge_cache);},
                                                    [&](id_t id) {return xg_index.node_sequence(id);});
        }
        else {
            return extract_extending_graph_internal(g, max_dist, pos, backward, preserve_cycles_on_src,
                                                    [&](id_t id) {return xg_index.edges_on_start(id);},
                                                    [&](id_t id) {return xg_index.edges_on_end(id);},
                                                    [&](id_t id) {return xg_index.node_sequence(id);});
        }
    }
}
}


