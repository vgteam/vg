#include "xg_position.hpp"

namespace vg {

Node xg_node(id_t id, xg::XG* xgidx) {
    return xgidx->node(id);
}

vector<Edge> xg_edges_on_start(id_t id, xg::XG* xgidx) {
    vector<Edge> all_edges = xgidx->edges_of(id);
    auto new_end = std::remove_if(all_edges.begin(), all_edges.end(),
                                  [&](const Edge& edge) {
                                      return (edge.from() == id && edge.from_start()) ||
                                             (edge.to() == id && !edge.to_end());
                                  });
    all_edges.resize(new_end - all_edges.begin());
    return all_edges;
}

vector<Edge> xg_edges_on_end(id_t id, xg::XG* xgidx) {
    vector<Edge> all_edges = xgidx->edges_of(id);
    auto new_end = std::remove_if(all_edges.begin(), all_edges.end(),
                                  [&](const Edge& edge) {
                                      return (edge.from() == id && !edge.from_start()) ||
                                             (edge.to() == id && edge.to_end());
                                  });
    all_edges.resize(new_end - all_edges.begin());
    return all_edges;
}

string xg_node_sequence(id_t id, xg::XG* xgidx) {
    return xgidx->node(id).sequence();
}

size_t xg_node_length(id_t id, xg::XG* xgidx) {
    return xgidx->node_length(id);
}

int64_t xg_node_start(id_t id, xg::XG* xgidx) {
    return xgidx->node_start(id);
}

char xg_pos_char(pos_t pos, xg::XG* xgidx) {
    return xgidx->pos_char(id(pos), is_rev(pos), offset(pos));
}

map<pos_t, char> xg_next_pos_chars(pos_t pos, xg::XG* xgidx) {

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

set<pos_t> xg_next_pos(pos_t pos, bool whole_node, xg::XG* xgidx) {
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

int64_t xg_distance(pos_t pos1, pos_t pos2, int64_t maximum, xg::XG* xgidx) {
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

set<pos_t> xg_positions_bp_from(pos_t pos, int64_t distance, bool rev, xg::XG* xgidx) {
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


}
