#include "xg_position.hpp"

//#define debug

namespace vg {

Node xg_node(id_t id, const XG* xgidx) {
    return xgidx->node(id);
}

vector<Edge> xg_edges_on_start(id_t id, const XG* xgidx) {
    vector<Edge> all_edges = xgidx->edges_of(id);
    auto new_end = std::remove_if(all_edges.begin(), all_edges.end(),
                                  [&](const Edge& edge) {
                                      return (edge.from() == id && edge.from_start()) ||
                                             (edge.to() == id && !edge.to_end());
                                  });
    all_edges.resize(new_end - all_edges.begin());
    return all_edges;
}

vector<Edge> xg_edges_on_end(id_t id, const XG* xgidx) {
    vector<Edge> all_edges = xgidx->edges_of(id);
    auto new_end = std::remove_if(all_edges.begin(), all_edges.end(),
                                  [&](const Edge& edge) {
                                      return (edge.from() == id && !edge.from_start()) ||
                                             (edge.to() == id && edge.to_end());
                                  });
    all_edges.resize(new_end - all_edges.begin());
    return all_edges;
}

string xg_node_sequence(id_t id, const XG* xgidx) {
    return xgidx->node(id).sequence();
}

size_t xg_node_length(id_t id, const XG* xgidx) {
    return xgidx->node_length(id);
}

int64_t xg_node_start(id_t id, const XG* xgidx) {
    return xgidx->node_start(id);
}

char xg_pos_char(pos_t pos, const XG* xgidx) {
    return xgidx->pos_char(id(pos), is_rev(pos), offset(pos));
}

map<pos_t, char> xg_next_pos_chars(pos_t pos, const XG* xgidx) {

    map<pos_t, char> nexts;
    // See if the node is cached (did we just visit it?)
    Node node = xgidx->node(id(pos));
    // if we are still in the node, return the next position and character
    if (offset(pos) < node.sequence().size()-1) {
        ++get_offset(pos);
        //nexts[pos] = xg_pos_char(pos, xgidx, node_cache);
        char c = node.sequence().at(offset(pos));
        if (is_rev(pos)) c = reverse_complement(c);
        nexts[pos] = c;
    } else {
        // helper
        auto is_inverting = [](const Edge& e) {
            return !(e.from_start() == e.to_end())
            && (e.from_start() || e.to_end());
        };
        // check our cache
        vector<Edge> edges = xgidx->edges_of(node.id());
        // look at the next positions we could reach
        if (!is_rev(pos)) {
            // we are on the forward strand, the next things from this node come off the end
            for (auto& edge : edges) {
                if((edge.to() == id(pos) && edge.to_end()) || (edge.from() == id(pos) && !edge.from_start())) {
                    id_t nid = (edge.from() == id(pos) ?
                                edge.to()
                                : edge.from());
                    pos_t p = make_pos_t(nid, is_inverting(edge), 0);
                    nexts[p] = xg_pos_char(p, xgidx);
                }
            }
        } else {
            // we are on the reverse strand, the next things from this node come off the start
            for (auto& edge : edges) {
                if((edge.to() == id(pos) && !edge.to_end()) || (edge.from() == id(pos) && edge.from_start())) {
                    id_t nid = (edge.to() == id(pos) ?
                                edge.from()
                                : edge.to());
                    pos_t p = make_pos_t(nid, !is_inverting(edge), 0);
                    nexts[p] = xg_pos_char(p, xgidx);
                }
            }
        }
    }
    return nexts;
}

set<pos_t> xg_next_pos(pos_t pos, bool whole_node, const XG* xgidx) {
    set<pos_t> nexts;
    // See if the node is cached (did we just visit it?)
    Node node = xgidx->node(id(pos));
    // if we are still in the node, return the next position and character
    if (!whole_node && offset(pos) < node.sequence().size()-1) {
        ++get_offset(pos);
        nexts.insert(pos);
    } else {
        // helper
        auto is_inverting = [](const Edge& e) {
            return !(e.from_start() == e.to_end())
            && (e.from_start() || e.to_end());
        };
        // check our cache
        auto edges = xgidx->edges_of(id(pos));
        // look at the next positions we could reach
        if (!is_rev(pos)) {
            // we are on the forward strand, the next things from this node come off the end
            for (auto& edge : edges) {
                if((edge.to() == id(pos) && edge.to_end()) || (edge.from() == id(pos) && !edge.from_start())) {
                    id_t nid = (edge.from() == id(pos) ?
                                edge.to()
                                : edge.from());
                    nexts.insert(make_pos_t(nid, is_inverting(edge), 0));
                }
            }
        } else {
            // we are on the reverse strand, the next things from this node come off the start
            for (auto& edge : edges) {
                if((edge.to() == id(pos) && !edge.to_end()) || (edge.from() == id(pos) && edge.from_start())) {
                    id_t nid = (edge.to() == id(pos) ?
                                edge.from()
                                : edge.to());
                    nexts.insert(make_pos_t(nid, !is_inverting(edge), 0));
                }
            }
        }
    }
    return nexts;
}

int64_t xg_distance(pos_t pos1, pos_t pos2, int64_t maximum, const XG* xgidx) {
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

set<pos_t> xg_positions_bp_from(pos_t pos, int64_t distance, bool rev, const XG* xgidx) {
    // handle base case
    //size_t xg_node_length(id_t id, XG* xgidx);
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

map<string, vector<pair<size_t, bool> > > xg_alignment_path_offsets(const XG* xgidx, const Alignment& aln, bool just_min,
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
        
        // Find the position of this end of this mapping
        pos_t mapping_pos = make_pos_t(mapping.position());
        
        // Find the positions for this end of this Mapping
        auto pos_offs = algorithms::nearest_offsets_in_paths(xgidx, mapping_pos,
                                                             nearby ? search_limit : -1);
        
        for (auto look_at_end : end) {
            // For the start and the end of the Mapping, as needed
            
            for (auto& p : pos_offs) {
                // For each path, splice the list of path positions for this
                // Mapping onto the end of the list of positions we found in that
                // path
                auto& v = offsets[xgidx->get_path_name(p.first)];
                
                for (pair<size_t, bool>& y : p.second) {
                    v.emplace_back(y.second ? y.first - mapping_width : y.first,
                                   y.second);
                }
                
#ifdef debug
                cerr << "\tFound hit on path " << xgidx->get_path_name(p.first) << endl;
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

void xg_annotate_with_initial_path_positions(const XG* xgidx, Alignment& aln, size_t search_limit) {
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

}
