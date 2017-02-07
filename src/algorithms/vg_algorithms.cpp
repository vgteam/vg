//
//  vg_algorithms.cpp
//
//

#include "vg_algorithms.hpp"

namespace vg {
namespace algorithms {
    unordered_map<id_t, id_t> extract_connecting_graph_internal(Graph& g, int64_t max_len,
                                                                const pos_t& pos_1, const pos_t& pos_2,
                                                                function<vector<Edge>(id_t)> edges_on_start,
                                                                function<vector<Edge>(id_t)> edges_on_end,
                                                                function<string(id_t)> node_sequence) {
        if (g.node_size() || g.edge_size()) {
            cerr << "error:[extract_connecting_graph] must extract into an empty graph" << endl;
            exit(1);
        }
        
        // a local struct for Nodes that maintains edge lists
        struct TemporaryNode {
            TemporaryNode() {}
            TemporaryNode(string sequence) : sequence(sequence) {}
            string sequence;
            // edges are stored as (node id, is reversing?)
            vector<pair<id_t, bool>> edges_left;
            vector<pair<id_t, bool>> edges_right;
        };
        
        // a local struct that packages a node traversal with its distance from the first position
        struct Traversal {
            Traversal(id_t id, bool rev, int64_t dist) : id(id), rev(rev), dist(dist) {}
            int64_t dist; // distance from pos_1 to the right side of this node
            id_t id;      // node ID
            bool rev;     // strand
            inline bool operator<(const Traversal& other) const {
                return dist > other.dist; // opposite order so priority queue selects minimum
            }
        };
        
        // for finding the largest node id in the subgraph
        id_t max_id = max(id(pos_1), id(pos_2));
        
        // the representation of the graph we're going to build up before storing in g (allows easier
        // subsetting operations than Graph or VG objects)
        unordered_map<id_t, TemporaryNode> temp_graph;
        temp_graph[id(pos_1)] = TemporaryNode(node_sequence(id(pos_1)));
        temp_graph[id(pos_2)] = TemporaryNode(node_sequence(id(pos_2)));
        
        bool found_target = false;
        
        // mark both positions as "queued" so that we won't look for additional traversals
        unordered_set<pair<id_t, bool>> queued_traversals{make_pair(id(pos_1), is_rev(pos_1)), make_pair(id(pos_2), is_rev(pos_2))};
        // initialize the queue
        priority_queue<Traversal> queue;
        
        // the distance to the end of the starting node
        int64_t first_traversal_length =  temp_graph[id(pos_1)].sequence.size() - offset(pos_1);
        
        // the max length of the part of a path preceding the final node in each direction
        int64_t forward_max_len = max_len - offset(pos_2);
        int64_t backward_max_len = max_len - first_traversal_length;
        
        // separately handle (common) edge case that both positions are on the same node
        // and the second is reachable from the first
        // TODO: is there a more elegant way to do this than as a special case? currently I need to
        // do it because the objects in the queue are implicitly at the "end" of a node traversal,
        // but here we would want to "stop" before reaching the end of the first node
        if (id(pos_1) == id(pos_2) && is_rev(pos_1) == is_rev(pos_2) && offset(pos_1) < offset(pos_2)) {
            found_target = (offset(pos_2) - offset(pos_1) < max_len);
        }
        else {
            
            // if we can reach the end of this node, init the queue with it
            if (first_traversal_length < forward_max_len) {
                queue.emplace(id(pos_1), is_rev(pos_1), first_traversal_length);
            }
            
            // search along a Dijkstra tree for shortest path
            while (!queue.empty()) {
                // get the next closest node to the starting position
                Traversal trav = queue.top();
                queue.pop();
                
                // which side are we traversing out of?
                auto& edges_out = trav.rev ? temp_graph[trav.id].edges_left : temp_graph[trav.id].edges_right;
                for (Edge& edge : (trav.rev ? edges_on_start(trav.id) : edges_on_end(trav.id))) {
                    // get the orientation and id of the other side of the edge
                    id_t next_id;
                    bool next_rev;
                    if (edge.from() == trav.id && edge.from_start() == trav.rev) {
                        next_id = edge.to();
                        next_rev = edge.to_end() != trav.rev;
                    }
                    else {
                        next_id = edge.from();
                        next_rev = edge.from_start() != trav.rev;
                    }
                    found_target = found_target || (next_id == id(pos_2) && next_rev == is_rev(pos_2));
                    max_id = max(max_id, next_id);
                    
                    // make sure the node is in
                    if (!temp_graph.count(next_id)) {
                        temp_graph[next_id] = TemporaryNode(node_sequence(next_id));
                    }
                    
                    // distance to the end of this node
                    int64_t dist_thru = trav.dist + temp_graph[next_id].sequence.size();
                    if (!queued_traversals.count(make_pair(next_id, next_rev)) && dist_thru < forward_max_len) {
                        // we can add more nodes along same path without going over the max length
                        // and we have not reached the target node yet
                        queue.emplace(next_id, next_rev, dist_thru);
                        queued_traversals.emplace(next_id, next_rev);
                    }
                    
                    // what side does this edge enter on the next node?
                    auto& edges_in = next_rev ? temp_graph[next_id].edges_right : temp_graph[next_id].edges_left;
                    // is the edge reversing?
                    bool reversing = (trav.rev != next_rev);
                    // add this edge to the edge list on the current node
                    edges_out.push_back(make_pair(next_id, reversing));
                    // add to other node, but if it is a self-loop to the same side don't add it twice
                    if (!(trav.id == next_id && reversing) ) {
                        edges_in.push_back(make_pair(trav.id, reversing));
                    }
                }
            }
        }
        // if we have to add new nodes, any id this large or larger will not have conflicts
        id_t next_id = max_id + 1;
        
        // a translator for node ids in g to node ids in the original graph
        unordered_map<id_t, id_t> id_trans;
        
        // leave g empty if we can't connect these positions
        if (!found_target) {
            return id_trans;
        }
        
        // now search backwards from the second position to remove anything that isn't on
        // some path between the two positions
        
        // edges stored as ((from, from_start), (to, to_end))
        unordered_set<pair<pair<id_t, bool>, pair<id_t, bool>>> reachable_edges;
        
        // put an edge from a given side of a node into canonical form
        auto canonical_form_edge = [](const pair<id_t, bool>& edge, id_t id, bool left) {
            if (edge.first < id || (left && edge.first == id)) {
                return make_pair(make_pair(edge.first, edge.second == left), make_pair(id, !left));
            }
            else {
                return make_pair(make_pair(id, left), make_pair(edge.first, edge.second != left));
            }
        };
        
        // re-initialize the queue
        queue.emplace(id(pos_2), !is_rev(pos_2), offset(pos_2));
        
        // reset the queued traversal list and add the two reverse traversals
        queued_traversals.clear();
        queued_traversals.insert(make_pair(id(pos_2), !is_rev(pos_2)));
        // again, mark the final node as "queued" so the search won't traverse through it
        queued_traversals.insert(make_pair(id(pos_1), !is_rev(pos_1)));
        
        while (!queue.empty()) {
            // get the next closest node traversal
            Traversal trav = queue.top();
            queue.pop();
            
            // the edges in the direction of this traversal
            auto& edges_out = trav.rev ? temp_graph[trav.id].edges_left : temp_graph[trav.id].edges_right;
            for (const pair<id_t, bool>& edge : edges_out) {
                
                // the distance to the opposite side of this next node
                int64_t dist_thru = trav.dist + temp_graph[edge.first].sequence.size();
                
                // have to represent edges canonically for hash to correctly identify duplicates
                // note: this edge is reachable even if we don't add the node to the stack because
                // the distance cutoff for the stack is to the opposite side of the node
                reachable_edges.insert(canonical_form_edge(edge, trav.id, trav.rev));
                
                // queue up the node traversal if it hasn't been seen before and is under the max distance
                pair<id_t, bool> next_trav = make_pair(edge.first, edge.second != trav.rev);
                if (!queued_traversals.count(next_trav) && dist_thru < backward_max_len) {
                    queue.emplace(next_trav.first, next_trav.second, dist_thru);
                    queued_traversals.insert(next_trav);
                }
            }
        }
        
        // remove any edges we didn't traverse and record which nodes to erase
        vector<unordered_map<id_t, TemporaryNode>::iterator> to_erase;
        for (auto iter = temp_graph.begin(); iter != temp_graph.end(); iter++) {
            id_t id = (*iter).first;
            if (queued_traversals.count(make_pair(id, true)) || queued_traversals.count(make_pair(id, false))) {
                // this node is on some path between the two positions
                TemporaryNode& node = (*iter).second;
                // move all edges we don't want to the end
                auto new_end = std::remove_if(node.edges_left.begin(), node.edges_left.end(),
                                              [&](const pair<id_t, bool>& edge) {
                                                  return !reachable_edges.count(canonical_form_edge(edge, id, true)); });
                // trim the vector down to only the ones we ant
                node.edges_left.erase(new_end, node.edges_left.end());
                // repeat for the other side
                new_end = std::remove_if(node.edges_right.begin(), node.edges_right.end(),
                                         [&](const pair<id_t, bool>& edge) {
                                             return !reachable_edges.count(canonical_form_edge(edge, id, false)); });
                node.edges_right.erase(new_end, node.edges_right.end());
                
            }
            else {
                // this node was found in the initial traversal, but is not actually on a path
                // between the two positions
                to_erase.push_back(iter);
            }
        }
        
        // remove the nodes
        for (auto& iter : to_erase) {
            temp_graph.erase(iter);
        }
        
        // the temporary graph now only consists of only the original nodes that we are going to keep,
        // so record them in the node id translator
        for (auto& node_record : temp_graph) {
            id_trans[node_record.first] = node_record.first;
        }
        
        // if there are edges traversed in both directions from the boundary position's nodes, then
        // they must be in cycles
        bool in_cycle_1 = !(temp_graph[id(pos_1)].edges_left.empty() || temp_graph[id(pos_1)].edges_right.empty());
        bool in_cycle_2 = !(temp_graph[id(pos_2)].edges_left.empty() || temp_graph[id(pos_2)].edges_right.empty());
        
        // functions to extract the part of the node string past the first and second positions:
        
        // get sequence to the right
        auto trimmed_seq_right = [](const string& seq, int64_t offset, bool rev) {
            if (rev) {
                return seq.substr(0, seq.size() - offset - 1);
            }
            else {
                return seq.substr(offset + 1, seq.size() - offset - 1);
            }
        };
        // get sequence to the left
        auto trimmed_seq_left = [](const string& seq, int64_t offset, bool rev) {
            if (rev) {
                return seq.substr(seq.size() - offset, offset - 1);
            }
            else {
                // can't let the length go negative because it's going to be cast to size_t
                return offset ? seq.substr(0, offset - 1) : string();
            }
        };
        
        // TODO: when we duplicate the source node, how is it going to be decided which is the
        // the source and which the sink?
        
        // cut the nodes containing the start positions, duplicating or trimming them if necessary
        if (id(pos_1) == id(pos_2)) {
            // the two positions are on the same node, so we have to be careful not to mess up
            // the situation for the second position while trimming or duplicating for the first
            if (is_rev(pos_1) == is_rev(pos_2) && offset(pos_1) < offset(pos_2)) {
                // this should have forestalled the search, so the graph is just the one node and we can trim it
                // from both sides
                TemporaryNode& node = temp_graph[id(pos_1)];
                node.sequence = is_rev(pos_1) ? node.sequence.substr(node.sequence.size() - offset(pos_2),  offset(pos_2) - offset(pos_1) - 1)
                : node.sequence.substr(offset(pos_1) + 1, offset(pos_2) - offset(pos_1) - 1);
            }
            else if (is_rev(pos_1) == is_rev(pos_2)) {
                // we found the target on the same strand of the same node, but by starting searching
                // in the opposite direction. this means the node is definitely in a cycle, but the
                // interior of the node was never traversed so we can cut it out
                
                TemporaryNode& node = temp_graph[id(pos_1)];
                temp_graph[next_id] = TemporaryNode(node.sequence);
                TemporaryNode& new_node = temp_graph[next_id];
                
                // move the edges from one side onto the new node
                new_node.edges_right = std::move(node.edges_right);
                node.edges_right.clear();
                
                // relabel the edges pointing back into this side
                for (pair<id_t, bool>& edge : new_node.edges_right) {
                    TemporaryNode& next_node = temp_graph[edge.first];
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
                    node.sequence = trimmed_seq_right(node.sequence, offset(pos_1), is_rev(pos_1));
                    new_node.sequence = trimmed_seq_left(new_node.sequence, offset(pos_2), is_rev(pos_2));
                }
                else {
                    id_trans[next_id] = id(pos_1);
                    new_node.sequence = trimmed_seq_right(new_node.sequence, offset(pos_1), is_rev(pos_1));
                    node.sequence = trimmed_seq_left(node.sequence, offset(pos_2), is_rev(pos_2));
                }
                
                next_id++;
            }
            else {
                // because they are on opposite strands, the traversals leave out of the same side of the
                // node, which may also be in a cycle
                
                TemporaryNode& node = temp_graph[id(pos_1)];
                
                if (in_cycle_1) {
                    // create a new node that will retain all of the edges to preserve cycles
                    temp_graph[next_id] = TemporaryNode(node.sequence);
                    TemporaryNode& new_node = temp_graph[next_id];
                    
                    // deep copy all edges
                    new_node.edges_left = node.edges_left;
                    new_node.edges_right = node.edges_right;
                    
                    // add edges back from all connected nodes
                    for (pair<id_t, bool>& edge : new_node.edges_left) {
                        TemporaryNode& next_node = temp_graph[edge.first];
                        auto& edges_backward = edge.second ? next_node.edges_left : next_node.edges_right;
                        edges_backward.emplace_back(next_id, edge.second);
                    }
                    for (pair<id_t, bool>& edge : new_node.edges_right) {
                        TemporaryNode& next_node = temp_graph[edge.first];
                        auto& edges_backward = edge.second ? next_node.edges_right : next_node.edges_left;
                        edges_backward.emplace_back(next_id, edge.second);
                    }
                    
                    id_trans[next_id] = id(pos_1);
                    next_id++;
                }
                
                // make a new node so we can trim the sequence in two places
                temp_graph[next_id] = TemporaryNode(node.sequence);
                TemporaryNode& new_node = temp_graph[next_id];
                
                if (is_rev(pos_1)) {
                    // deep copy the edge list
                    new_node.edges_left = node.edges_left;
                    // add edges backward to the new node
                    for (pair<id_t, bool>& edge : new_node.edges_left) {
                        TemporaryNode& next_node = temp_graph[edge.first];
                        auto& edges_backward = edge.second ? next_node.edges_left : next_node.edges_right;
                        edges_backward.emplace_back(next_id, edge.second);
                    }
                    // remove the node from all backwards edge lists on the other side
                    for (pair<id_t, bool>& edge : node.edges_right) {
                        TemporaryNode& next_node = temp_graph[edge.first];
                        auto& edges_backward = edge.second ? next_node.edges_right : next_node.edges_left;
                        // find the edge backward in the other node's edge list
                        for (auto iter = edges_backward.begin(); iter != edges_backward.end(); iter++) {
                            if ((*iter).first == id(pos_1)) {
                                // erase it
                                edges_backward.erase(iter);
                                break;
                            }
                        }
                    }
                    // remove the edges in from the other side on the old node
                    node.edges_right.clear();
                }
                else {
                    // deep copy the edge list
                    new_node.edges_right = node.edges_right;
                    // add edges backward to the new node
                    for (pair<id_t, bool>& edge : new_node.edges_right) {
                        TemporaryNode& next_node = temp_graph[edge.first];
                        auto& edges_backward = edge.second ? next_node.edges_right : next_node.edges_left;
                        edges_backward.emplace_back(next_id, edge.second);
                    }
                    // remove the node from all backwards edge lists on the other side
                    for (pair<id_t, bool>& edge : node.edges_left) {
                        TemporaryNode& next_node = temp_graph[edge.first];
                        auto& edges_backward = edge.second ? next_node.edges_left : next_node.edges_right;
                        // find the edge backward in the other node's edge list
                        for (auto iter = edges_backward.begin(); iter != edges_backward.end(); iter++) {
                            if ((*iter).first == id(pos_1)) {
                                // erase it
                                edges_backward.erase(iter);
                                break;
                            }
                        }
                    }
                    // remove the edges in from the other side on the old node
                    node.edges_left.clear();
                }
                
                //  cut the sequences of the two nodes according to the search positions
                node.sequence = trimmed_seq_right(node.sequence, offset(pos_1), is_rev(pos_1));
                new_node.sequence = trimmed_seq_left(new_node.sequence, offset(pos_2), is_rev(pos_2));
                // record the translation
                id_trans[next_id] = id(pos_2);
            }
        }
        else {
            // the two positions are on separate nodes, so we can handle the two sides independently
            
            TemporaryNode& node_1 = temp_graph[id(pos_1)];
            // duplicate the node if it's in a cycle
            if (in_cycle_1) {
                temp_graph[next_id] = TemporaryNode(node_1.sequence);
                TemporaryNode& new_node = temp_graph[next_id];
                // copy the edges going out of the side that the traversal leaves
                auto& new_edges = is_rev(pos_1) ? new_node.edges_left : new_node.edges_right;
                for (pair<id_t, bool>& edge : is_rev(pos_1) ? node_1.edges_left : node_1.edges_right) {
                    TemporaryNode& next_node = temp_graph[edge.first];
                    auto& edges_backward = is_rev(pos_1) != edge.second ? next_node.edges_right : next_node.edges_left;
                    new_edges.push_back(edge);
                    edges_backward.emplace_back(next_id, edge.second);
                }
                // record the translation
                id_trans[next_id] = id(pos_1);
                // set up to cut the new node instead the previous one
                node_1 = new_node;
                next_id++;
            }
            // cut the sequence
            node_1.sequence = trimmed_seq_right(node_1.sequence, offset(pos_1), is_rev(pos_1));
            
            TemporaryNode& node_2 = temp_graph[id(pos_2)];
            // duplicate the node if it's in a cycle
            if (in_cycle_2) {
                temp_graph[next_id] = TemporaryNode(node_2.sequence);
                TemporaryNode& new_node = temp_graph[next_id];
                // copy the edges going out of the side that the traversal leaves
                auto& new_edges = is_rev(pos_2) ? new_node.edges_right : new_node.edges_left;
                for (pair<id_t, bool>& edge : is_rev(pos_2) ? node_2.edges_right : node_2.edges_left) {
                    TemporaryNode& next_node = temp_graph[edge.first];
                    auto& edges_backward = is_rev(pos_2) != edge.second ? next_node.edges_left : next_node.edges_right;
                    new_edges.push_back(edge);
                    edges_backward.emplace_back(next_id, edge.second);
                }
                // record the translation
                id_trans[next_id] = id(pos_2);
                // set up to cut the new node instead the previous one
                node_2 = new_node;
                next_id++;
            }
            // cut the sequence
            node_2.sequence = trimmed_seq_left(node_2.sequence, offset(pos_2), is_rev(pos_2));
        }
        
        // transfer to g
        for (const pair<id_t, TemporaryNode>& node_record : temp_graph) {
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
        // the function, which are obviously available in the environment that calls this)
        return id_trans;
    }
    
    // wrapper for connecting graph algorithm using VG to access graph elements
    unordered_map<id_t, id_t> extract_connecting_graph(VG& vg, Graph& g, int64_t max_len,
                                                       const pos_t& pos_1, const pos_t& pos_2){
        
        return extract_connecting_graph_internal(g, max_len, pos_1, pos_2,
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
                                                       const pos_t& pos_1, const pos_t& pos_2) {
        return extract_connecting_graph_internal(g, max_len, pos_1, pos_2,
                                                 [&](id_t id) {return xg_index.edges_on_start(id);},
                                                 [&](id_t id) {return xg_index.edges_on_end(id);},
                                                 [&](id_t id) {return xg_index.node_sequence(id);});
    }
}
}


