/**
 * \file scanner.cpp
 * Implementations for traversing Protpbuf object trees
 */

#include "scanner.hpp"

namespace vg {

using namespace std;

// Specializations have to be defined in dependency order to avoid complaints of specialization after instantiation.

template<>
bool PositionIDScanner<Mapping>::scan(const Mapping& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee) {
    
    if (msg.position().node_id() != 0) {
        // Just enumerate the position we have
        return pos_iteratee(msg.position());
    } else {
        // Skip this empty Position
        return true;
    }
}

template<>
bool PositionIDScanner<Path>::scan(const Path& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee) {
    
    // If we don't see any real positions or IDs in this path, we have to emit a 0 node ID sentinel.
    bool path_is_empty = true;
    
    // We have to wrap the iteratees to do this
    auto record_pos = [&](const Position& pos) -> bool {
        path_is_empty = false;
        return pos_iteratee(pos);
    };
    auto record_id = [&](const id_t& id) -> bool {
        path_is_empty = false;
        return id_iteratee(id);
    };
    
    bool keep_going = true;
    for (size_t i = 0; keep_going && i < msg.mapping_size(); i++) {
        // Scan over all the mappings
        keep_going &= PositionIDScanner<Mapping>::scan(msg.mapping(i), record_pos, record_id);
    }
    
    if (keep_going && path_is_empty) {
        // Visit the sentinel zero node ID
        keep_going &= id_iteratee(0);
    }
    
    return keep_going;
}

template<>
bool PositionIDScanner<Alignment>::scan(const Alignment& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee) {
    
    // Visit the Path
    return PositionIDScanner<Path>::scan(msg.path(), pos_iteratee, id_iteratee);
}

template<>
bool PositionIDScanner<Node>::scan(const Node& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee) {

    if (msg.id() != 0) {
        // Just announce the node's ID
        return id_iteratee(msg.id());
    }
    return true;

}

template<>
bool PositionIDScanner<Edge>::scan(const Edge& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee) {

    bool keep_going = true;
    
    // Make sure to filter out zero node IDs

    if (msg.from() != 0) {
        keep_going &= id_iteratee(msg.from());
    }
    
    if (keep_going && msg.to() != 0) {
        keep_going &= id_iteratee(msg.to());
    }

    return keep_going;

}

template<>
bool PositionIDScanner<Graph>::scan(const Graph& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee) {
    
    // If we don't see any real positions or IDs in this graph, we have to emit a 0 node ID sentinel.
    bool graph_is_empty = true;
    
    // We have to wrap the iteratees to do this
    auto record_pos = [&](const Position& pos) -> bool {
        graph_is_empty = false;
        return pos_iteratee(pos);
    };
    auto record_id = [&](const id_t& id) -> bool {
        graph_is_empty = false;
        return id_iteratee(id);
    };
    
    // Note that it's OK if we catch a 0 node ID sentinel from a contained
    // path. Then 0 has already been emitted and we don't need to emit it again
    // for the graph as a whole, even if the graph has no nodes/edges.
    
    bool keep_going = true;
    for (size_t i = 0; keep_going && i < msg.node_size(); i++) {
        // Scan over all the nodes
        keep_going &= PositionIDScanner<Node>::scan(msg.node(i), record_pos, record_id);
    }
    for (size_t i = 0; keep_going && i < msg.edge_size(); i++) {
        // Scan over all the edges
        keep_going &= PositionIDScanner<Edge>::scan(msg.edge(i), record_pos, record_id);
    }
    for (size_t i = 0; keep_going && i < msg.path_size(); i++) {
        // Scan over all the paths
        keep_going &= PositionIDScanner<Path>::scan(msg.path(i), record_pos, record_id);
    }

    if (keep_going && graph_is_empty) {
        // Visit the sentinel zero node ID
        keep_going &= id_iteratee(0);
    }
    
    return keep_going;
    
}
    
    
}
