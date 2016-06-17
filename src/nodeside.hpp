#ifndef VG_NODESIDE_HPP
#define VG_NODESIDE_HPP

#include <ostream>
#include <utility>

#include "vg.pb.h"
#include "types.hpp"
#include "hash_map.hpp"

namespace vg {

using namespace std;

// Represents one side of a Node, identified by ID, for the purposes of
// indexing edges.
class NodeSide {
public:
    id_t node;
    bool is_end;

    // We need this to be a converting constructor so we can represent the empty and deleted item keys in a pair_hash_map.
    inline NodeSide(id_t node, bool is_end = false): node(node), is_end(is_end) {
        // Nothing to do
    }

    inline NodeSide(): NodeSide(0, false) {
        // Nothing to do
    }

    inline bool operator==(const NodeSide& other) const {
        return node == other.node && is_end == other.is_end;
    }

    inline bool operator!=(const NodeSide& other) const {
        return node != other.node || is_end != other.is_end;
    }

    inline bool operator<(const NodeSide& other) const {
        return node < other.node || (node == other.node && is_end < other.is_end);
    }

    // Make an edge into a canonically ordered pair of NodeSides
    static inline pair<NodeSide, NodeSide> pair_from_edge(Edge* e) {
        return minmax(NodeSide(e->from(), !e->from_start()), NodeSide(e->to(), e->to_end()));
    }

    // Make an edge into a canonically ordered pair of NodeSides
    static inline pair<NodeSide, NodeSide> pair_from_edge(const Edge& e) {
        return minmax(NodeSide(e.from(), !e.from_start()), NodeSide(e.to(), e.to_end()));
    }

    // Make a canonically ordered pair of NodeSides from an edge off of the
    // start of a node, to another node in the given relative orientation.
    static inline pair<NodeSide, NodeSide> pair_from_start_edge(id_t start_id, const pair<id_t, bool>& oriented_other) {
        // If it's in the same relative orientation, we go to its end.
        return minmax(NodeSide(start_id, false), NodeSide(oriented_other.first, !oriented_other.second));
    }

    // Make a canonically ordered pair of NodeSides from an edge off of the
    // end of a node, to another node in the given relative orientation.
    static inline pair<NodeSide, NodeSide> pair_from_end_edge(id_t end_id, const pair<id_t, bool>& oriented_other) {
        // If it's in the same relative orientation, we go to its start.
        return minmax(NodeSide(end_id, true), NodeSide(oriented_other.first, oriented_other.second));
    }

    // reverse complement the node side
    inline NodeSide flip(void) const {
        return NodeSide(node, !is_end);
    }

};

// helpers to be more clear
inline NodeSide node_start(id_t id) {
    return NodeSide(id, false);
}

inline NodeSide node_end(id_t id) {
    return NodeSide(id, true);
}

inline ostream& operator<<(ostream& out, const NodeSide& nodeside) {
    return out << nodeside.node << " " << (nodeside.is_end ? "end" : "start");
}

}

namespace std {
// We need to implement a hash function for these if we want to be able to use them in keys.
template <> struct hash<vg::NodeSide>
{
    // Produce a hash of a NodeSide
    size_t operator()(const vg::NodeSide& item) const
    {
        // Hash it just as we would a pair.
        return hash<pair<id_t, bool>>()(make_pair(item.node, item.is_end));
    }
};
}




#endif
