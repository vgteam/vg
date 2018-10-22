#ifndef VG_HANDLE_HPP_INCLUDED
#define VG_HANDLE_HPP_INCLUDED

/** \file 
 * Defines a handle type that can refer to oriented nodes of, and be used to
 * traverse, any backing graph implementation. Not just an ID or a pos_t because
 * XG (and maybe other implementations) provide more efficient local traversal
 * mechanisms if you can skip ID lookups.
 */

#include "types.hpp"
#include "vg.pb.h"
#include "hash_map.hpp"

#include <functional>
#include <cstdint>
#include <vector>

namespace vg {

using namespace std;

/// A handle is 8 (assuming id_t is still int64_t) opaque bytes.
/// A handle refers to an oriented node.
/// Two handles are equal iff they refer to the same orientation of the same node.
/// Only handles in the same graph may be compared.
/// Handles have no ordering, but can be hashed.
/// Bit usage should be right-justified in machine endian order (i.e. a graph
/// implementation must use the lowest bits of the lowest byte first, and must
/// only use other bits of other bytes when the number of nodes/handles is
/// sufficiently large as to require itin the chosen representation).
/// This allows for overlays that can wrap their backing graphs' handles without using any more space.
/// TODO: Note that this precludes us using pointers as handles unless we make them wider.
struct handle_t { char data[sizeof(id_t)];
};

typedef pair<handle_t, handle_t> edge_t;

// XG is going to store node index in its g-vector and node orientation in there
// VG is going to store ID and orientation
// Other implementations can store other things (or maybe int indexes into tables)

/// View a handle as an integer
inline uint64_t& as_integer(handle_t& handle) {
    return reinterpret_cast<uint64_t&>(handle);
}

/// View a const handle as a const integer
inline const uint64_t& as_integer(const handle_t& handle) {
    return reinterpret_cast<const uint64_t&>(handle);
}

/// View an integer as a handle
inline handle_t& as_handle(uint64_t& value) {
    return reinterpret_cast<handle_t&>(value);
}

/// View a const integer as a const handle
inline const handle_t& as_handle(const uint64_t& value) {
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

/// Define a way to pack an integer and an orientation bit into a handle_t. XG
/// and VG both ought to use these functions instead of doing fiddly bit
/// packing themselves.
struct EasyHandlePacking {

    /// Extract the packed integer
    inline static uint64_t unpack_number(const handle_t& handle) {
        return as_integer(handle) >> 1;
    }
    
    /// Extract the packed bit
    inline static bool unpack_bit(const handle_t& handle) {
        return as_integer(handle) & 1;
    }
    
    /// Pack up an integer and a bit into a handle
    inline static handle_t pack(const uint64_t& number, const bool& bit) {
        // Make sure the number doesn't use all the bits
        assert(number < (0x1ULL << 63));
        
        return as_handle((number << 1) | (bit ? 1 : 0));
    }
    
    /// Toggle the packed bit and return a new handle
    inline static handle_t toggle_bit(const handle_t& handle) {
        return as_handle(as_integer(handle) ^ 1);
    }

};

/**
 * Define hashes for handles.
 */
template<>
struct wang_hash<handle_t> {
    size_t operator()(const vg::handle_t& handle) const {
        return wang_hash<std::int64_t>()(as_integer(handle));
    }
};
    
struct path_handle_t {
    char data[sizeof(int64_t)];
};
    
/// View a path handle as an integer
inline int64_t& as_integer(path_handle_t& path_handle) {
    return reinterpret_cast<int64_t&>(path_handle);
}

/// View a const path handle as a const integer
inline const int64_t& as_integer(const path_handle_t& path_handle) {
    return reinterpret_cast<const int64_t&>(path_handle);
}

/// View an integer as a path handle
inline path_handle_t& as_path_handle(int64_t& value) {
    return reinterpret_cast<path_handle_t&>(value);
}

/// View a const integer as a const path handle
inline const path_handle_t& as_path_handle(const int64_t& value) {
    return reinterpret_cast<const path_handle_t&>(value);
}

/// Define equality on path handles
inline bool operator==(const path_handle_t& a, const path_handle_t& b) {
    return as_integer(a) == as_integer(b);
}

/// Define inequality on path handles
inline bool operator!=(const path_handle_t& a, const path_handle_t& b) {
    return as_integer(a) != as_integer(b);
}

struct occurrence_handle_t {
    char data[2 * sizeof(int64_t)];
};
    
/// View an occurrence handle as an integer
inline int64_t* as_integers(occurrence_handle_t& occurrence_handle) {
    return reinterpret_cast<int64_t*>(&occurrence_handle);
}

/// View a const occurrence handle as a const integer
inline const int64_t* as_integers(const occurrence_handle_t& occurrence_handle) {
    return reinterpret_cast<const int64_t*>(&occurrence_handle);
}

/// Define equality on occurrence handles
inline bool operator==(const occurrence_handle_t& a, const occurrence_handle_t& b) {
    return as_integers(a)[0] == as_integers(b)[0] && as_integers(a)[1] == as_integers(b)[1];
}

/// Define inequality on occurrence handles
inline bool operator!=(const occurrence_handle_t& a, const occurrence_handle_t& b) {
    return !(a == b);
}
    
}   // namespace vg

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

/**
 * Define hashes for path handles.
 */
template<> struct hash<vg::path_handle_t> {
public:
    inline size_t operator()(const vg::path_handle_t& path_handle) const {
        return std::hash<int64_t>()(vg::as_integer(path_handle));
    }
};

}

namespace vg {

using namespace std;

/**
 * This is the interface that a graph that uses handles needs to support.
 * It is also the interface that users should code against.
 */
class HandleGraph {
public:

    ////////////////////////////////////////////////////////////////////////////
    // Interface that needs to be implemented
    ////////////////////////////////////////////////////////////////////////////
    
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse = false) const = 0;
    
    /// Get the ID from a handle
    virtual id_t get_id(const handle_t& handle) const = 0;
    
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const = 0;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    virtual handle_t flip(const handle_t& handle) const = 0;
    
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const = 0;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual string get_sequence(const handle_t& handle) const = 0;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    virtual bool follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const = 0;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    virtual void for_each_handle(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const = 0;
    
    /// Return the number of nodes in the graph
    /// TODO: can't be node_count because XG has a field named node_count.
    virtual size_t node_size() const = 0;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual id_t min_node_id() const = 0;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual id_t max_node_id() const = 0;
    
    ////////////////////////////////////////////////////////////////////////////
    // Interface that needs to be using'd
    ////////////////////////////////////////////////////////////////////////////
    
    /// Loop over all the handles to next/previous (right/left) nodes. Works
    /// with a callback that just takes all the handles and returns void.
    /// MUST be pulled into implementing classes with `using` in order to work!
    template <typename T>
    auto follow_edges(const handle_t& handle, bool go_left, T&& iteratee) const
    -> typename std::enable_if<std::is_void<decltype(iteratee(get_handle(0, false)))>::value>::type {
        // Implementation only for void-returning iteratees
        // We ought to just overload on the std::function but that's not allowed until C++14.
        // See <https://stackoverflow.com/q/13811180>
        
        // We also can't use result_of<T(handle_t)>::type to sniff the return
        // type out because that ::type would not exist (since that's what you
        // get for a void apparently?) and we couldn't check if it's bool or
        // void.
        
        // So we do this nonsense thing with a trailing return type (to get the
        // actual arg into scope) and a decltype (which is allowed to resolve to
        // void) and is_void (which is allowed to take void) and a fake
        // get_handle call (which is the shortest handle_t-typed expression I
        // could think of).
        
        // Make a wrapper that puts a bool return type on.
        function<bool(const handle_t&)> lambda = [&](const handle_t& found) {
            iteratee(found);
            return true;
        };
        
        // Use that
        follow_edges(handle, go_left, lambda);
        
        // During development I managed to get earlier versions of this template to build infinitely recursive functions.
        static_assert(!std::is_void<decltype(lambda(get_handle(0, false)))>::value, "can't take our own lambda");
    }
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Works with void-returning iteratees.
    /// MUST be pulled into implementing classes with `using` in order to work!
    template <typename T>
    auto for_each_handle(T&& iteratee, bool parallel = false) const
    -> typename std::enable_if<std::is_void<decltype(iteratee(get_handle(0, false)))>::value>::type {
        // Make a wrapper that puts a bool return type on.
        function<bool(const handle_t&)> lambda = [&](const handle_t& found) {
            iteratee(found);
            return true;
        };
        
        // Use that
        for_each_handle(lambda, parallel);
    }
    
    /// Get a handle from a Visit Protobuf object.
    /// Must be using'd to avoid shadowing.
    handle_t get_handle(const Visit& visit) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Additional optional interface with a default implementation
    ////////////////////////////////////////////////////////////////////////////
    
    /// Get the number of edges on the right (go_left = false) or left (go_left
    /// = true) side of the given handle. The default implementation is O(n) in
    /// the number of edges returned, but graph implementations that track this
    /// information more efficiently can override this method.
    virtual size_t get_degree(const handle_t& handle, bool go_left) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Concrete utility methods
    ////////////////////////////////////////////////////////////////////////////
    
    /// Get a Protobuf Visit from a handle.
    Visit to_visit(const handle_t& handle) const;
    
    /// Get the locally forward version of a handle
    handle_t forward(const handle_t& handle) const;
    
    /// A pair of handles can be used as an edge. When so used, the handles have a
    /// canonical order and orientation.
    edge_t edge_handle(const handle_t& left, const handle_t& right) const;
    
    /// Such a pair can be viewed from either inward end handle and produce the
    /// outward handle you would arrive at.
    handle_t traverse_edge_handle(const edge_t& edge, const handle_t& left) const;
};
    
/**
 * This is the interface for a handle graph that stores embedded paths.
 */
class PathHandleGraph : virtual public HandleGraph {
public:
    
    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface that needs to be implemented
    ////////////////////////////////////////////////////////////////////////////
    
    /// Look up the path handle for the given path name
    virtual path_handle_t get_path_handle(const string& path_name) const = 0;
    
    /// Look up the name of a path from a handle to it
    virtual string get_path_name(const path_handle_t& path_handle) const = 0;
    
    /// Returns the number of node occurrences in the path
    virtual size_t get_occurrence_count(const path_handle_t& path_handle) const = 0;
    
    /// Returns the number of paths stored in the graph
    virtual size_t get_path_count() const = 0;
    
    /// Execute a function on each path in the graph
    virtual void for_each_path_handle(const function<void(const path_handle_t&)>& iteratee) const = 0;
    
    /// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
    virtual handle_t get_occurrence(const occurrence_handle_t& occurrence_handle) const = 0;
    
    /// Get a handle to the first occurrence in a path
    virtual occurrence_handle_t get_first_occurrence(const path_handle_t& path_handle) const = 0;
    
    /// Get a handle to the last occurrence in a path
    virtual occurrence_handle_t get_last_occurrence(const path_handle_t& path_handle) const = 0;
    
    /// Returns true if the occurrence is not the last occurence on the path, else false
    virtual bool has_next_occurrence(const occurrence_handle_t& occurrence_handle) const = 0;
    
    /// Returns true if the occurrence is not the first occurence on the path, else false
    virtual bool has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const = 0;
    
    /// Returns a handle to the next occurrence on the path
    virtual occurrence_handle_t get_next_occurrence(const occurrence_handle_t& occurrence_handle) const = 0;
    
    /// Returns a handle to the previous occurrence on the path
    virtual occurrence_handle_t get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const = 0;
    
    /// Returns a handle to the path that an occurrence is on
    virtual path_handle_t get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const = 0;
    
    /// Returns the 0-based ordinal rank of a occurrence on a path
    virtual size_t get_ordinal_rank_of_occurrence(const occurrence_handle_t& occurrence_handle) const = 0;
};

/**
 * This is the interface for a handle graph that supports modification.
 */
class MutableHandleGraph : virtual public HandleGraph {
public:
    /*
     * Note: All operations may invalidate path handles and occurrence handles.
     */
    
    /// Create a new node with the given sequence and return the handle.
    virtual handle_t create_handle(const string& sequence) = 0;

    /// Create a new node with the given id and sequence, then return the handle.
    virtual handle_t create_handle(const string& sequence, const id_t& id) = 0;
    
    /// Remove the node belonging to the given handle and all of its edges.
    /// Does not update any stored paths.
    virtual void destroy_handle(const handle_t& handle) = 0;
    
    /// Create an edge connecting the given handles in the given order and orientations.
    /// Ignores existing edges.
    virtual void create_edge(const handle_t& left, const handle_t& right) = 0;
    
    /// Convenient wrapper for create_edge.
    inline void create_edge(const edge_t& edge) {
        create_edge(edge.first, edge.second);
    }
    
    /// Remove the edge connecting the given handles in the given order and orientations.
    /// Ignores nonexistent edges.
    /// Does not update any stored paths.
    virtual void destroy_edge(const handle_t& left, const handle_t& right) = 0;
    
    /// Convenient wrapper for destroy_edge.
    inline void destroy_edge(const edge_t& edge) {
        destroy_edge(edge.first, edge.second);
    }
    
    /// Remove all nodes and edges. Does not update any stored paths.
    virtual void clear() = 0;
    
    /// Swap the nodes corresponding to the given handles, in the ordering used
    /// by for_each_handle when looping over the graph. Other handles to the
    /// nodes being swapped must not be invalidated. If a swap is made while
    /// for_each_handle is running, it affects the order of the handles
    /// traversed during the current traversal (so swapping an already seen
    /// handle to a later handle's position will make the seen handle be visited
    /// again and the later handle not be visited at all).
    virtual void swap_handles(const handle_t& a, const handle_t& b) = 0;
    
    /// Alter the node that the given handle corresponds to so the orientation
    /// indicated by the handle becomes the node's local forward orientation.
    /// Rewrites all edges pointing to the node and the node's sequence to
    /// reflect this. Invalidates all handles to the node (including the one
    /// passed). Returns a new, valid handle to the node in its new forward
    /// orientation. Note that it is possible for the node's ID to change.
    /// Does not update any stored paths. May change the ordering of the underlying
    /// graph.
    virtual handle_t apply_orientation(const handle_t& handle) = 0;
    
    /// Split a handle's underlying node at the given offsets in the handle's
    /// orientation. Returns all of the handles to the parts. Other handles to
    /// the node being split may be invalidated. The split pieces stay in the
    /// same local forward orientation as the original node, but the returned
    /// handles come in the order and orientation appropriate for the handle
    /// passed in.
    /// Updates stored paths.
    virtual vector<handle_t> divide_handle(const handle_t& handle, const vector<size_t>& offsets) = 0;
    
    /// Specialization of divide_handle for a single division point
    inline pair<handle_t, handle_t> divide_handle(const handle_t& handle, size_t offset) {
        auto parts = divide_handle(handle, vector<size_t>{offset});
        return make_pair(parts.front(), parts.back());
    }
}; 

}

#endif
