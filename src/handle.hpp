#ifndef VG_HANDLE_HPP_INCLUDED
#define VG_HANDLE_HPP_INCLUDED

/** \file 
 * Defines a handle type that can refer to oriented nodes of, and be used to
 * traverse, any backing graph implementation. Not just an ID or a pos_t because
 * XG (and maybe other implementations) provide more efficient local traversal
 * mechanisms if you can skip ID lookups.
 */

#include "types.hpp"
#include <functional>
#include <cstdint>

namespace vg {

using namespace std;

/// A handle is 8 (assuming id_t is still int64_t) opaque bytes.
/// A handle refers to an oriented node.
/// Two handles are equal iff they refer to the same orientation of the same node in the same graph.
/// Handles have no ordering, but can be hashed.
struct handle_t {
    char data[sizeof(id_t)];
};

// XG is going to store node index in its g-vector and node orientation in there
// VG is going to store ID and orientation
// Other implementations can store other things (or maybe int indexes into tables)

/// View a handle as an integer
inline int64_t& as_integer(handle_t& handle) {
    return reinterpret_cast<int64_t&>(handle);
}

/// View a const handle as a const integer
inline const int64_t& as_integer(const handle_t& handle) {
    return reinterpret_cast<const int64_t&>(handle);
}

/// View an integer as a handle
inline handle_t& as_handle(int64_t& value) {
    return reinterpret_cast<handle_t&>(value);
}

/// View a const integer as a const handle
inline const handle_t& as_handle(const int64_t& value) {
    return reinterpret_cast<const handle_t&>(value);
}

/// Define equality on handles
inline bool operator==(const handle_t& a, const handle_t& b) {
    return as_integer(a) == as_integer(b);
}

/// Define inequality on handles
inline bool operator!=(const handle_t& a, const handle_t& b) {
    return as_integer(a) != as_integer(b);
}

}

// This needs to be outside the vg namespace

namespace std {

/**
 * Define hashes for handles.
 */
template<> struct hash<vg::handle_t> {
public:
    inline size_t operator()(const vg::handle_t& handle) const {
        return std::hash<int64_t>()(vg::as_integer(handle));
    }
};

}

namespace vg {

using namespace std;

/**
 * This is the interface that a graph that uses handles needs to support.
 */
class HandleGraph {

    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse) const = 0;
    
    /// Get the ID from a handle
    virtual id_t get_id(const handle_t& handle) const = 0;
    
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const = 0;
    
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const = 0;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual string get_sequence(const handle_t& handle) const = 0;
    
    /// Loop over all the handles to next nodes. Passes them to a callback which
    /// returns false to stop iterating and true to continue.
    virtual void get_next(const handle_t& handle, const function<bool(const handle_t&)>& iteratee) = 0;
    
    /// Loop over all the handles to previous nodes. Passes them to a callback which
    /// returns false to stop iterating and true to continue.
    virtual void get_prev(const handle_t& handle, const function<bool(const handle_t&)>& iteratee) = 0;
    
};

}

#endif
