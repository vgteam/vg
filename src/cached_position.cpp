#include "cached_position.hpp"

namespace vg {

Node xg_cached_node(id_t id, XG* xgidx, LRUCache<id_t, Node>& node_cache) {
    pair<Node, bool> cached = node_cache.retrieve(id);
    if(!cached.second) {
        cached.first = xgidx->node(id);
        node_cache.put(id, cached.first);
    }
    Node& node = cached.first;
    return node;
}

vector<Edge> xg_cached_edges_of(id_t id, XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache) {
    pair<vector<Edge>, bool> cached = edge_cache.retrieve(id);
    if(!cached.second) {
        for (auto& edge : xgidx->edges_of(id)) {
            cached.first.push_back(edge);
        }
        edge_cache.put(id, cached.first);
    }
    return cached.first;
}
    
vector<Edge> xg_cached_edges_on_start(id_t id, XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache) {
    vector<Edge> all_edges = xg_cached_edges_of(id, xgidx, edge_cache);
    auto new_end = std::remove_if(all_edges.begin(), all_edges.end(),
                                  [&](const Edge& edge) {
                                      return (edge.from() == id && edge.from_start()) ||
                                             (edge.to() == id && !edge.to_end());
                                  });
    all_edges.resize(new_end - all_edges.begin());
    return all_edges;
}

vector<Edge> xg_cached_edges_on_end(id_t id, XG* xgidx, LRUCache<id_t, vector<Edge> >& edge_cache) {
    vector<Edge> all_edges = xg_cached_edges_of(id, xgidx, edge_cache);
    auto new_end = std::remove_if(all_edges.begin(), all_edges.end(),
                                  [&](const Edge& edge) {
                                      return (edge.from() == id && !edge.from_start()) ||
                                             (edge.to() == id && edge.to_end());
                                  });
    all_edges.resize(new_end - all_edges.begin());
    return all_edges;
}

string xg_cached_node_sequence(id_t id, XG* xgidx, LRUCache<id_t, Node>& node_cache) {
    pair<Node, bool> cached = node_cache.retrieve(id);
    if(!cached.second) {
        cached.first = xgidx->node(id);
        node_cache.put(id, cached.first);
    }
    Node& node = cached.first;
    return node.sequence();
}

size_t xg_cached_node_length(id_t id, XG* xgidx, LRUCache<id_t, Node>& node_cache) {
    pair<Node, bool> cached = node_cache.retrieve(id);
    if(!cached.second) {
        cached.first = xgidx->node(id);
        node_cache.put(id, cached.first);
    }
    Node& node = cached.first;
    return node.sequence().size();
}

int64_t xg_cached_node_start(id_t id, XG* xgidx, LRUCache<id_t, int64_t>& node_start_cache) {
    pair<int64_t, bool> cached = node_start_cache.retrieve(id);
    if(!cached.second) {
        cached.first = (int64_t)xgidx->node_start(id);
        node_start_cache.put(id, cached.first);
    }
    return cached.first;
}

char xg_cached_pos_char(pos_t pos, XG* xgidx, LRUCache<id_t, Node>& node_cache) {
    pair<Node, bool> cached = node_cache.retrieve(id(pos));
    if(!cached.second) {
        // If it's not in the cache, put it in
        cached.first = xgidx->node(id(pos));
        node_cache.put(id(pos), cached.first);
    }
    Node& node = cached.first;
    if (is_rev(pos)) {
        /*
        cerr << "reversed... " << endl;
        cerr << "rev pos " << offset(reverse(pos, node.sequence().size())) << endl;
        cerr << "seq is " << node.sequence() << " and got " <<
            reverse_complement(node.sequence()[offset(reverse(pos, node.sequence().size()))-1]) << endl;
        */
        return reverse_complement(node.sequence()[offset(reverse(pos, node.sequence().size()))-1]);
    } else {
        /*
        cerr << "forward... " << endl;
        cerr << "seq is " << node.sequence() << " and got " << node.sequence().at(offset(pos)) << endl;
        */
        return node.sequence().at(offset(pos));
    }
}

map<pos_t, char> xg_cached_next_pos_chars(pos_t pos, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache) {

    map<pos_t, char> nexts;
    // See if the node is cached (did we just visit it?)
    pair<Node, bool> cached = node_cache.retrieve(id(pos));
    if(!cached.second) {
        // If it's not in the cache, put it in
        cached.first = xgidx->node(id(pos));
        node_cache.put(id(pos), cached.first);
    }
    Node& node = cached.first;
    // if we are still in the node, return the next position and character
    if (offset(pos) < node.sequence().size()-1) {
        ++get_offset(pos);
        nexts[pos] = xg_cached_pos_char(pos, xgidx, node_cache);
    } else {
        // helper
        auto is_inverting = [](const Edge& e) {
            return !(e.from_start() == e.to_end())
            && (e.from_start() || e.to_end());
        };
        // check our cache
        pair<vector<Edge>, bool> cached = edge_cache.retrieve(id(pos));
        if(!cached.second) {
            // If it's not in the cache, put it in
            for (auto& edge : xgidx->edges_of(id(pos))) {
                cached.first.push_back(edge);
            }
            edge_cache.put(id(pos), cached.first);
        }
        auto& edges = cached.first;
        // look at the next positions we could reach
        if (!is_rev(pos)) {
            // we are on the forward strand, the next things from this node come off the end
            for (auto& edge : edges) {
                if((edge.to() == id(pos) && edge.to_end()) || (edge.from() == id(pos) && !edge.from_start())) {
                    id_t nid = (edge.from() == id(pos) ?
                                edge.to()
                                : edge.from());
                    pos_t p = make_pos_t(nid, is_inverting(edge), 0);
                    nexts[p] = xg_cached_pos_char(p, xgidx, node_cache);
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
                    nexts[p] = xg_cached_pos_char(p, xgidx, node_cache);
                }
            }
        }
    }
    return nexts;
}

set<pos_t> xg_cached_next_pos(pos_t pos, bool whole_node, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache) {
    set<pos_t> nexts;
    // See if the node is cached (did we just visit it?)
    pair<Node, bool> cached = node_cache.retrieve(id(pos));
    if(!cached.second) {
        // If it's not in the cache, put it in
        cached.first = xgidx->node(id(pos));
        node_cache.put(id(pos), cached.first);
    }
    Node& node = cached.first;
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
        pair<vector<Edge>, bool> cached = edge_cache.retrieve(id(pos));
        if(!cached.second) {
            // If it's not in the cache, put it in
            for (auto& edge : xgidx->edges_of(id(pos))) {
                cached.first.push_back(edge);
            }
            edge_cache.put(id(pos), cached.first);
        }
        auto& edges = cached.first;
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

int64_t xg_cached_distance(pos_t pos1, pos_t pos2, int64_t maximum, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache) {
    //cerr << "distance from " << pos1 << " to " << pos2 << endl;
    if (pos1 == pos2) return 0;
    int64_t adj = (offset(pos1) == xg_cached_node_length(id(pos1), xgidx, node_cache) ? 0 : 1);
    set<pos_t> seen;
    set<pos_t> nexts = xg_cached_next_pos(pos1, false, xgidx, node_cache, edge_cache);
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
                for (auto& x : xg_cached_next_pos(next, false, xgidx, node_cache, edge_cache)) {
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

set<pos_t> xg_cached_positions_bp_from(pos_t pos, int64_t distance, bool rev, XG* xgidx, LRUCache<id_t, Node>& node_cache, LRUCache<id_t, vector<Edge> >& edge_cache) {
    // handle base case
    //size_t xg_cached_node_length(id_t id, XG* xgidx, LRUCache<id_t, Node>& node_cache);
    if (rev) {
        pos = reverse(pos, xg_cached_node_length(id(pos), xgidx, node_cache));
    }
    set<pos_t> positions;
    if (distance == 0) {
        positions.insert(pos);
        //return positions;
    } else {
        set<pos_t> seen;
        set<pos_t> nexts = xg_cached_next_pos(pos, false, xgidx, node_cache, edge_cache);
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
                    for (auto& x : xg_cached_next_pos(next, false, xgidx, node_cache, edge_cache)) {
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
            rev_pos.insert(reverse(p, xg_cached_node_length(id(p), xgidx, node_cache)));
        }
        return rev_pos;
    } else {
        return positions;
    }
}


}
