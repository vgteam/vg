#include "position.hpp"

namespace vg {

pos_t make_pos_t(const Position& pos) {
    return make_tuple(pos.node_id(), pos.is_reverse(), pos.offset());
}

pos_t make_pos_t(id_t id, bool is_rev, off_t off) {
    return make_tuple(id, is_rev, off);
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

Position make_position(id_t id, bool is_rev, off_t off) {
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

bool is_empty(const pos_t& pos) {
    return id(pos) == 0;
}

id_t id(const pos_t& pos) {
    return get<0>(pos);
}

bool is_rev(const pos_t& pos) {
    return get<1>(pos);
}

off_t offset(const pos_t& pos) {
    return get<2>(pos);
}

id_t& get_id(pos_t& pos) {
    return get<0>(pos);
}

bool& get_is_rev(pos_t& pos) {
    return get<1>(pos);
}

off_t& get_offset(pos_t& pos) {
    return get<2>(pos);
}

pos_t reverse(const pos_t& pos, size_t node_length) {
    pos_t rev = pos;
    // swap the offset onto the other strand
    get_offset(rev) = node_length - offset(rev);
    // invert the position
    get_is_rev(rev) = !is_rev(rev);
    return rev;
}

Position reverse(const Position& pos, size_t node_length) {
    auto p = pos;
    p.set_offset(node_length - pos.offset());
    p.set_is_reverse(!pos.is_reverse());
    return p;
}

ostream& operator<<(ostream& out, const pos_t& pos) {
    return out << id(pos) << (is_rev(pos) ? "-" : "+") << offset(pos);
}

pair<int64_t, int64_t> min_oriented_distances(const map<string, vector<pair<size_t, bool> > >& path_offsets1,
                                              const map<string, vector<pair<size_t, bool> > >& path_offsets2) {
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

}
