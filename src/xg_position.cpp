#include "xg_position.hpp"

namespace vg {

Node xg_node(id_t id, const xg::XG* xgidx) {
    Node node;
    node.set_id(id);
    node.set_sequence(xgidx->get_sequence(xgidx->get_handle(id)));
    return node;
}

string xg_node_sequence(id_t id, const xg::XG* xgidx) {
    return xgidx->get_sequence(xgidx->get_handle(id));
}

size_t xg_node_length(id_t id, const xg::XG* xgidx) {
    return xgidx->get_length(xgidx->get_handle(id));
}

int64_t xg_node_start(id_t id, const xg::XG* xgidx) {
    return xgidx->node_start(id);
}

char xg_pos_char(pos_t pos, const xg::XG* xgidx) {
    return xgidx->pos_char(id(pos), is_rev(pos), offset(pos));
}

map<pos_t, char> xg_next_pos_chars(pos_t pos, const xg::XG* xgidx) {

    map<pos_t, char> nexts;
    // See if the node is cached (did we just visit it?)
    Node node = xg_node(id(pos), xgidx);
    // if we are still in the node, return the next position and character
    if (offset(pos) < node.sequence().size()-1) {
        ++get_offset(pos);
        //nexts[pos] = xg_pos_char(pos, xgidx, node_cache);
        char c = node.sequence().at(offset(pos));
        if (is_rev(pos)) c = reverse_complement(c);
        nexts[pos] = c;
    } else {
        handlegraph::handle_t handle = xgidx->get_handle(id(pos), is_rev(pos));
        xgidx->follow_edges(handle, false, [&](const handlegraph::handle_t& next) {
                pos_t p = make_pos_t(xgidx->get_id(next), xgidx->get_is_reverse(next), 0);
                nexts[p] = xg_pos_char(p, xgidx);
            });
    }
    return nexts;
}

set<pos_t> xg_next_pos(pos_t pos, bool whole_node, const xg::XG* xgidx) {
    set<pos_t> nexts;
    // See if the node is cached (did we just visit it?)
    Node node = xg_node(id(pos), xgidx);
    // if we are still in the node, return the next position and character
    if (!whole_node && offset(pos) < node.sequence().size()-1) {
        ++get_offset(pos);
        nexts.insert(pos);
    } else {
        handlegraph::handle_t handle = xgidx->get_handle(id(pos), is_rev(pos));
        xgidx->follow_edges(handle, false, [&](const handlegraph::handle_t& next) {
                nexts.insert(make_pos_t(xgidx->get_id(next), xgidx->get_is_reverse(next), 0));
            });
    }
    return nexts;
}

int64_t xg_distance(pos_t pos1, pos_t pos2, int64_t maximum, const xg::XG* xgidx) {
    //cerr << "distance from " << pos1 << " to " << pos2 << endl;
    if (pos1 == pos2) return 0;
    int64_t adj = (offset(pos1) == xg_node_length(id(pos1), xgidx) ? 0 : 1);
    set<pos_t> seen;
    set<pos_t> nexts = xg_next_pos(pos1, false, xgidx);
    int64_t distance = 0;
    while (!nexts.empty()) {
        set<pos_t> todo;
        for (auto& next : nexts) {
            if (!seen.count(next)) {
                seen.insert(next);
                if (next == pos2) {
                    return distance+adj;
                }
                // handle the edge case that we are looking for the position after the end of this node
                if (make_pos_t(id(next), is_rev(next), offset(next)+1) == pos2) {
                    return distance+adj+1;
                }
                for (auto& x : xg_next_pos(next, false, xgidx)) {
                    todo.insert(x);
                }
            }
        }
        if (distance == maximum) {
            break;
        }
        nexts = todo;
        ++distance;
    }
    return numeric_limits<int64_t>::max();
}

set<pos_t> xg_positions_bp_from(pos_t pos, int64_t distance, bool rev, const xg::XG* xgidx) {
    // handle base case
    //size_t xg_node_length(id_t id, xg::XG* xgidx);
    if (rev) {
        pos = reverse(pos, xg_node_length(id(pos), xgidx));
    }
    set<pos_t> positions;
    if (distance == 0) {
        positions.insert(pos);
        //return positions;
    } else {
        set<pos_t> seen;
        set<pos_t> nexts = xg_next_pos(pos, false, xgidx);
        int64_t walked = 0;
        while (!nexts.empty()) {
            if (walked+1 == distance) {
                for (auto& next : nexts) {
                    positions.insert(next);
                }
                break;
            }
            set<pos_t> todo;
            for (auto& next : nexts) {
                if (!seen.count(next)) {
                    seen.insert(next);
                    for (auto& x : xg_next_pos(next, false, xgidx)) {
                        todo.insert(x);
                    }
                }
            }
            nexts = todo;
            ++walked;
        }
    }
    if (rev) {
        set<pos_t> rev_pos;
        for (auto& p : positions) {
            rev_pos.insert(reverse(p, xg_node_length(id(p), xgidx)));
        }
        return rev_pos;
    } else {
        return positions;
    }
}

map<string, vector<pair<size_t, bool> > > xg_alignment_path_offsets(const xg::XG* xgidx, const Alignment& aln, bool just_min,
    bool nearby, size_t search_limit) {
    
    if (nearby && search_limit == 0) {
        // Fill in the search limit
        search_limit = aln.sequence().size();
    }
    
#ifdef debug
    cerr << "Searching for path positions for " << aln.name();
    if (nearby) {
        cerr << " within " << search_limit << " bp";
    } else {
        cerr << " that are actually touched";
    }
    cerr << endl;
#endif
    
    map<string, vector<pair<size_t, bool> > > offsets;
    for (auto& mapping : aln.path().mapping()) {
    
        // How many bases does this Mapping cover over?
        size_t mapping_width = mapping_from_length(mapping);
        
        if (mapping_width == 0 && !nearby) {
            // Just skip over this mapping; it touches no bases.
            continue;
        }
    
        // We may have to consider both the starts and ends of mappings
        vector<bool> end = {false};
        if (just_min && !nearby) {
            // We want the min actually touched position along each path. It
            // could come from the Mapping start or the Mapping end.
            end.push_back(true);
        }
        
        for (auto look_at_end : end) {
            // For the start and the end of the Mapping, as needed
            
            if (mapping_width == 1 && look_at_end) {
                // The end is the same as the start so we don't need to look it
                // up separately.
                continue;
            }
       
            // Find the position of this end of this mapping
            pos_t mapping_pos = make_pos_t(mapping.position());
            if (look_at_end) {
                // Look at the end of the mapping instead of the start
                get_offset(mapping_pos) += mapping_width - 1;
            }
       
            // Find the positions for this end of this Mapping
            auto pos_offs = (nearby ?
                             xgidx->nearest_offsets_in_paths(mapping_pos, search_limit)
                             : xgidx->offsets_in_paths(mapping_pos));
            for (auto& p : pos_offs) {
                // For each path, splice the list of path positions for this
                // Mapping onto the end of the list of positions we found in that
                // path
                auto& v = offsets[xgidx->get_path_name(p.first)];
                auto& y = p.second;
                v.reserve(v.size() + distance(y.begin(),y.end()));
                v.insert(v.end(),y.begin(),y.end());
                
#ifdef debug
                cerr << "\tFound hit on path " << p.first << endl; 
#endif
            }
        }
    }
    if (!nearby && offsets.empty()) {
        // find the nearest if we couldn't find any before
        return xg_alignment_path_offsets(xgidx, aln, just_min, true, search_limit);
    }
    if (just_min) {
        // We need the minimum position for each path
        for (auto& p : offsets) {
            auto& v = p.second;
            auto m = *min_element(v.begin(), v.end(),
                                  [](const pair<size_t, bool>& a,
                                     const pair<size_t, bool>& b)
                                  { return a.first < b.first; });
            v.clear();
            v.push_back(m);
        }
    }
    return offsets;
}

void xg_annotate_with_initial_path_positions(const xg::XG* xgidx, Alignment& aln, size_t search_limit) {
    if (!aln.refpos_size()) {
        auto init_path_positions = xg_alignment_path_offsets(xgidx, aln, true, false, search_limit);
        for (const pair<string, vector<pair<size_t, bool> > >& pos_record : init_path_positions) {
            for (auto& pos : pos_record.second) {
                Position* refpos = aln.add_refpos();
                refpos->set_name(pos_record.first);
                refpos->set_offset(pos.first);
                refpos->set_is_reverse(pos.second);
            }
        }
    }
}

void xg_neighborhood(const xg::XG& xgidx,int64_t id, size_t dist, Graph& g, bool use_steps) {
    if (xgidx.has_node(id)) {
        Node& node = *g.add_node();
        node.set_id(id);
        node.set_sequence(xgidx.get_sequence(xgidx.get_handle(id)));
        xg_expand_context(xgidx, g, dist, true, use_steps);
    }
}

void xg_expand_context(const xg::XG& xgidx, Graph& g, size_t dist, bool add_paths, bool use_steps,
                        bool expand_forward, bool expand_backward,
                        int64_t until_node) {
    if (use_steps) {
        xg_expand_context_by_steps(xgidx, g, dist, add_paths, expand_forward, expand_backward, until_node);
    } else {
        xg_expand_context_by_length(xgidx, g, dist, add_paths, expand_forward, expand_backward, until_node);
    }
}

void xg_expand_context_by_steps(const xg::XG& xgidx, Graph& g, size_t steps, bool add_paths,
                                bool expand_forward, bool expand_backward,
                                int64_t until_node) {
    map<int64_t, Node*> nodes;
    map<pair<side_t, side_t>, Edge*> edges;
    set<int64_t> to_visit;
    // start with the nodes in the graph
    for (size_t i = 0; i < g.node_size(); ++i) {
        to_visit.insert(g.node(i).id());
        // handles the single-node case: we should still get the paths
        Node* np = g.mutable_node(i);
        nodes[np->id()] = np;
    }
    for (size_t i = 0; i < g.edge_size(); ++i) {
        auto& edge = g.edge(i);
        to_visit.insert(edge.from());
        to_visit.insert(edge.to());
        edges[make_pair(make_side(edge.from(), edge.from_start()),
                        make_side(edge.to(), edge.to_end()))] = g.mutable_edge(i);
    }
    // and expand
    for (size_t i = 0; i < steps; ++i) {
        set<int64_t> to_visit_next;
        for (auto id : to_visit) {
            // build out the graph
            // if we have nodes we haven't seeen
            if (nodes.find(id) == nodes.end()) {
                Node* np = g.add_node();
                nodes[id] = np;
                *np = node(id);
            }
            vector<Edge> edges_todo;
            if (expand_forward && expand_backward) {
                edges_todo = edges_of(id);
            } else if (expand_forward) {
                edges_todo = edges_from(id);
            } else if (expand_backward) {
                edges_todo = edges_to(id);
            } else {
                cerr << "[xg] error: Requested neither forward no backward context expansion" << endl;
                exit(1);
            }
            for (auto& edge : edges_todo) {
                auto sides = make_pair(make_side(edge.from(), edge.from_start()),
                                       make_side(edge.to(), edge.to_end()));
                if (edges.find(sides) == edges.end()) {
                    Edge* ep = g.add_edge(); *ep = edge;
                    edges[sides] = ep;
                }
                if (edge.from() == id) {
                    to_visit_next.insert(edge.to());
                } else {
                    to_visit_next.insert(edge.from());
                }
            }
            if (until_node != 0 && nodes.find(until_node) != nodes.end()) {
                break;
            }
        }
        to_visit = to_visit_next;
    }
    // then add connected nodes that we have edges to but didn't pull in yet.
    // These are the nodes reached on the last step; we won't follow their edges
    // to new noded.
    set<int64_t> last_step_nodes;
    for (auto& e : edges) {
        auto& edge = e.second;
        // get missing nodes
        int64_t f = edge->from();
        if (nodes.find(f) == nodes.end()) {
            Node* np = g.add_node();
            nodes[f] = np;
            *np = node(f);
            last_step_nodes.insert(f);
        }
        int64_t t = edge->to();
        if (nodes.find(t) == nodes.end()) {
            Node* np = g.add_node();
            nodes[t] = np;
            *np = node(t);
            last_step_nodes.insert(t);
        }
    }
    // We do need to find edges that connect the nodes we just grabbed on the
    // last step. Otherwise we'll produce something that isn't really a useful
    // subgraph, because there might be edges connecting the nodes you have that
    // you don't see.
    for(auto& n : last_step_nodes) {
        for (auto& edge : edges_from(n)) {
            if(last_step_nodes.count(edge.to())) {
                // This edge connects two nodes that were added on the last
                // step, and so wouldn't have been found by the main loop.
                
                // Determine if it's been seen already (somehow).
                // TODO: it probably shouldn't have been, unless it's a self loop or something.
                auto sides = make_pair(make_side(edge.from(), edge.from_start()),
                                       make_side(edge.to(), edge.to_end()));
                if (edges.find(sides) == edges.end()) {
                    // If it isn't there already, add it to the graph
                    Edge* ep = g.add_edge(); *ep = edge;
                    edges[sides] = ep;
                }
            }
        }       
    }
    // Edges between the last step nodes and other nodes will have already been
    // pulled in, on the step when those other nodes were processed by the main
    // loop.
    if (add_paths) {
        add_paths_to_graph(nodes, g);
    }
}

void xg_expand_context_by_length(const xg::XG& xgidx, Graph& g, size_t length, bool add_paths,
                                 bool expand_forward, bool expand_backward,
                                 int64_t until_node) {

    // map node_id --> min-distance-to-left-side, min-distance-to-right-side
    // these distances include the length of the node in the table. 
    map<int64_t, pair<int64_t, int64_t> > node_table;
    // nodes and edges in graph, so we don't duplicate when we add to protobuf
    map<int64_t, Node*> nodes;
    map<pair<side_t, side_t>, Edge*> edges;
    // bfs queue (id, enter-on-left-size)
    queue<int64_t> to_visit;

    // add starting graph with distance 0
    for (size_t i = 0; i < g.node_size(); ++i) {
        Node* np = g.mutable_node(i);
        node_table[np->id()] = pair<int64_t, int64_t>(0, 0);
        nodes[np->id()] = np;
        to_visit.push(np->id());
    }

    // add starting edges
    for (size_t i = 0; i < g.edge_size(); ++i) {
        auto& edge = g.edge(i);
        edges[make_pair(make_side(edge.from(), edge.from_start()),
                        make_side(edge.to(), edge.to_end()))] = g.mutable_edge(i);
    }

    // expand outward breadth-first
    while (!to_visit.empty() && (until_node == 0 || nodes.find(until_node) == nodes.end())) {
        int64_t id = to_visit.front();
        to_visit.pop();
        pair<int64_t, int64_t> dists = node_table[id];
        if (dists.first < length || dists.second < length) {
            vector<Edge> edges_todo;
            if (expand_forward && expand_backward) {
                edges_todo = edges_of(id);
            } else if (expand_forward) {
                edges_todo = edges_from(id);
            } else if (expand_backward) {
                edges_todo = edges_to(id);
            } else {
                cerr << "[xg] error: Requested neither forward no backward context expansion" << endl;
                exit(1);
            }
            for (auto& edge : edges_todo) {
                // update distance table with other end of edge
                function<void(int64_t, bool, bool)> lambda = [&](
                    int64_t other, bool from_start, bool to_end) {

                    int64_t dist = !from_start ? dists.first : dists.second;
                    Node other_node = node(other);
                    int64_t other_dist = dist + other_node.sequence().size();
                    if (dist < length) {
                        auto it = node_table.find(other);
                        bool updated = false;
                        if (it == node_table.end()) {
                            auto entry = make_pair(numeric_limits<int64_t>::max(),
                                                   numeric_limits<int64_t>::max());
                            it = node_table.insert(make_pair(other, entry)).first;
                            updated = true;
                        }
                        if (!to_end && other_dist < it->second.first) {
                            updated = true;
                            node_table[other].first = other_dist;
                        } else if (to_end && other_dist < it->second.second) {
                            updated = true;
                            node_table[other].second = other_dist;
                        }
                        // create the other node
                        if (nodes.find(other) == nodes.end()) {
                            Node* np = g.add_node();
                            nodes[other] = np;
                            *np = other_node;
                        }
                        // create all links back to graph, so as not to break paths
                        for (auto& other_edge : edges_of(other)) {
                            auto sides = make_pair(make_side(other_edge.from(),
                                                             other_edge.from_start()),
                                                   make_side(other_edge.to(),
                                                             other_edge.to_end()));
                            int64_t other_from = other_edge.from() == other ? other_edge.to() :
                                other_edge.from();
                            if (nodes.find(other_from) != nodes.end() &&
                                edges.find(sides) == edges.end()) {
                                Edge* ep = g.add_edge(); *ep = other_edge;
                                edges[sides] = ep;
                            }
                        }
                        // revisit the other node
                        if (updated) {
                            // this may be overly conservative (bumping any updated node)
                            to_visit.push(other);
                        }
                    }
                };
                // we can actually do two updates if we have a self loop, hence no else below
                if (edge.from() == id) {
                    lambda(edge.to(), edge.from_start(), edge.to_end());
                }
                if (edge.to() == id) {
                    lambda(edge.from(), !edge.to_end(), !edge.from_start());
                }
            }
        }
    }

    if (add_paths) {
        xg_add_paths_to_graph(xgidx, nodes, g);
    }
}
    
// if the graph ids partially ordered, this works no prob
// otherwise... owch
// the paths become disordered due to traversal of the node ids in order
void xg_add_paths_to_graph(const xg::XG& xgidx, map<int64_t, Node*>& nodes, Graph& g) {
    // map from path name to (map from mapping rank to mapping)
    map<string, map<size_t, Mapping>> paths;
    // mappings without 
    map<string, vector<Mapping>> unplaced;
    // use:
    //size_t node_position_in_path(int64_t id, const string& name) const;

    // pick up current paths in the graph
    for (size_t i = 0; i < g.path_size(); ++i) {
        auto& path = g.path(i);
        for (size_t j = 0; j < path.mapping_size(); ++j) {
            auto& m = path.mapping(j);
            if (m.rank()) {
                paths[path.name()][m.rank()] = m;
            } else {
                unplaced[path.name()].push_back(m);
            }                
        }
    }
    // do the same for the mappings in the list of nodes
    for (auto& n : nodes) {
        auto& id = n.first;
        auto& node = n.second;
        for (auto& n : node_mappings(id)) {
            auto& name = n.first;
            for (auto& m : n.second) {
                if (m.rank()) {
                    paths[name][m.rank()] = m;
                } else {
                    unplaced[name].push_back(m);
                }
            }
        }
    }
    // rebuild graph's paths
    // NB: mapping ranks allow us to remove this bit
    // only adding what we haven't seen before
    g.clear_path();
    for (auto& p : paths) {
        auto& name = p.first;
        auto& mappings = p.second;
        Path* path = g.add_path();
        path->set_name(name);
        for (auto& n : mappings) {
            *path->add_mapping() = n.second;
        }
        if (unplaced.find(name) != unplaced.end()) {
            auto& unp = unplaced[name];
            for (auto& m : unp) {
                *path->add_mapping() = m;
            }
        }
    }
}

void xg_get_id_range(const xg::XG& xgidx, int64_t id1, int64_t id2, Graph& g) {
    id1 = max(min_id, id1);
    id2 = min(max_id, id2);
    for (auto i = id1; i <= id2; ++i) {
        if (xgidx.has_node(i)) {
            Node& node = *g.add_node();
            node.set_id(i);
            node.set_sequence(xgidx.get_sequence(xgidx.get_handle(i)));
        }
    }
}

}
