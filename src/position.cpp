#include "position.hpp"

namespace vg {

pos_t make_pos_t(const Position& pos) {
    return make_tuple(pos.node_id(), pos.is_reverse(), pos.offset());
}

pos_t make_pos_t(const position_t& pos) {
    return make_tuple(pos.node_id(), pos.is_reverse(), pos.offset());
}
    
pos_t make_pos_t(gcsa::node_type node) {
    return make_tuple(gcsa::Node::id(node), gcsa::Node::rc(node), gcsa::Node::offset(node));
}
    
Position make_position(const pos_t& pos) {
    Position p;
    p.set_node_id(id(pos));
    p.set_is_reverse(is_rev(pos));
    p.set_offset(offset(pos));
    return p;
}

Position make_position(id_t id, bool is_rev, offset_t off) {
    Position p;
    p.set_node_id(id);
    p.set_is_reverse(is_rev);
    p.set_offset(off);
    return p;
}
    
Position make_position(gcsa::node_type node) {
    Position p;
    p.set_node_id(gcsa::Node::id(node));
    p.set_is_reverse(gcsa::Node::rc(node));
    p.set_offset(gcsa::Node::offset(node));
    return p;
}

Position reverse(const Position& pos, size_t node_length) {
    auto p = pos;
    p.set_offset(node_length - pos.offset());
    p.set_is_reverse(!pos.is_reverse());
    return p;
}

pair<int64_t, int64_t> min_oriented_distances(const unordered_map<path_handle_t, vector<pair<size_t, bool> > >& path_offsets1,
                                              const unordered_map<path_handle_t, vector<pair<size_t, bool> > >& path_offsets2) {
    int64_t distance_same = std::numeric_limits<int64_t>::max();
    int64_t distance_diff = std::numeric_limits<int64_t>::max();
    for (auto& path : path_offsets1) {
        auto& name = path.first;
        auto f = path_offsets2.find(name);
        if (f != path_offsets2.end()) {
            auto& pos1 = path.second;
            auto& pos2 = f->second;
            for (auto& p1 : pos1) {
                for (auto& p2 : pos2) {
                    int64_t proposal = abs((int64_t)p1.first - (int64_t)p2.first);
                    if (p1.second == p2.second) {
                        distance_same = min(distance_same, proposal);
                    } else if (p1.second != p2.second) {
                        distance_diff = min(distance_diff, proposal);
                    }
                }
            }
        }
    }
    return make_pair(distance_same, distance_diff);
}

string debug_string(const position_t& pos) {
    string to_return = "{id: " + to_string(pos.node_id()) + ", off: " + to_string(pos.offset()) + ", rev: " + (pos.is_reverse() ? "T" : "F") + "}";
    return to_return;
}

void from_proto_position(const Position& from, position_t& to) {
    to.set_node_id(from.node_id());
    to.set_offset(from.offset());
    to.set_is_reverse(from.is_reverse());
}

}
