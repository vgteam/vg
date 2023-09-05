#ifndef VG_SCANNER_HPP_INCLUDED
#define VG_SCANNER_HPP_INCLUDED

/**
 * \file scanner.hpp
 * Define some "scanners" that can traverse a tree of VG
 * Protobuf objects and iterate over items found in the tree.
 */

#include <vg/vg.pb.h>
#include "types.hpp"
#include <functional>

namespace vg {

using namespace std;

/** 
 * We define a PositionIDScanner that scans a VG Protobuf message tree for whole Position objects and node IDs.
 * Each of the two is visited with its own iteratee function.
 * May emit the same Position or id multiple times.
 * Will never emit an empty Position.
 * Will only emit the 0 node ID if a Graph or Path contains no nonzero node IDs.
 */
template<typename Message>
struct PositionIDScanner {
    /// Scan over the Position objects and non-Position-wrapped node IDs in
    /// this message and all its children. Returns false if an iteratee
    /// returned false and asked to stop.
    static bool scan(const Message& msg, const function<bool(const Position&)>& pos_iteratee,
        const function<bool(const id_t&)>& id_iteratee);
};

/**
 * We define an IDScanner which scans over all node ID references in a tree of VG Protobuf objects.
 * May emit the same ID multiple times.
 * Will only emit the 0 node ID if a Graph or Path contains no nonzero node IDs.
 */
template<typename Message>
struct IDScanner {
    /// Scan over the node IDs in this message and all its children.
    /// Returns false if an iteratee returned false and asked to stop.
    static bool scan(const Message& msg, const function<bool(const id_t&)>& iteratee);
};

/**
 * We define a PositionScanner which scans over all Position objects and node
 * IDs, wrapped as Positions, in a tree of VG Protobuf objects.
 * Will only emit the empty Position if a Graph or Path contains no nonzero node IDs.
 */
template<typename Message>
struct WrappingPositionScanner {
    /// Scan over the node IDs in this message and all its children.
    /// Returns false if an iteratee returned false and asked to stop.
    static bool scan(const Message& msg, const function<bool(const Position&)>& iteratee);
};

/////////////
// Template Specializations
/////////////

// Declare specializations of the above that we will implement

template<>
bool PositionIDScanner<Mapping>::scan(const Mapping& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee);
    
template<>
bool PositionIDScanner<Path>::scan(const Path& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee);
    
template<>
bool PositionIDScanner<Alignment>::scan(const Alignment& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee);
    
template<>
bool PositionIDScanner<Node>::scan(const Node& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee);
    
template<>
bool PositionIDScanner<Edge>::scan(const Edge& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee);
    
template<>
bool PositionIDScanner<Graph>::scan(const Graph& msg, const function<bool(const Position&)>& pos_iteratee,
    const function<bool(const id_t&)>& id_iteratee);

/////////////
// Template Implementations
/////////////


// We implement everything in terms of the PositionIDScanner.

template<typename Message>
bool IDScanner<Message>::scan(const Message& msg, const function<bool(const id_t&)>& iteratee) {
    // Get the node ID form the position and iterate over that
    return PositionIDScanner<Message>::scan(msg, [&](const Position& pos) {
        return iteratee(pos.node_id());
    }, iteratee);
}

template<typename Message>
bool WrappingPositionScanner<Message>::scan(const Message& msg, const function<bool(const Position&)>& iteratee) {
    // Wrap the node ID and iterate over that
    return PositionIDScanner<Message>::scan(msg, iteratee, [&](const id_t& id) {
        Position wrapped;
        wrapped.set_node_id(id);
        return iteratee(wrapped);
    });
}

}


#endif
