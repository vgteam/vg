//
//  packed_graph.hpp
//  
//  Contains a implementation of a sequence graph based on bit-packed integer
//  vectors.
//

#ifndef VG_PACKED_GRAPH_HPP_INCLUDED
#define VG_PACKED_GRAPH_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <vector>
#include <utility>
#include <functional>
#include "handle.hpp"
#include "packed_structs.hpp"


namespace vg {
    
using namespace std;


class PackedGraph : public DeletableHandleGraph {
        
public:
    PackedGraph();
    ~PackedGraph();
    
    // Method to check if a node exists by ID
    virtual bool has_node(id_t node_id) const;
    
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    virtual id_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    virtual handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward orientation.
    virtual string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    virtual bool follow_edges(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    virtual void for_each_handle(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    /// TODO: can't be node_count because XG has a field named node_count.
    virtual size_t node_size(void) const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual id_t min_node_id(void) const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual id_t max_node_id(void) const;
    
    
/**
 * This is the interface for a handle graph that stores embedded paths.
 */
    
//    ////////////////////////////////////////////////////////////////////////////
//    // Path handle interface that needs to be implemented
//    ////////////////////////////////////////////////////////////////////////////
//
//    /// Determine if a path name exists and is legal to get a path handle for.
//    bool has_path(const std::string& path_name) const;
//
//    /// Look up the path handle for the given path name.
//    /// The path with that name must exist.
//    path_handle_t get_path_handle(const std::string& path_name) const;
//
//    /// Look up the name of a path from a handle to it
//    std::string get_path_name(const path_handle_t& path_handle) const;
//
//    /// Returns the number of node occurrences in the path
//    size_t get_occurrence_count(const path_handle_t& path_handle) const;
//
//    /// Returns the number of paths stored in the graph
//    size_t get_path_count() const;
//
//    /// Execute a function on each path in the graph
//    // TODO: allow stopping early?
//    void for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const;
//
//    /// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
//    handle_t get_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Get a handle to the first occurrence in a path.
//    /// The path MUST be nonempty.
//    occurrence_handle_t get_first_occurrence(const path_handle_t& path_handle) const;
//
//    /// Get a handle to the last occurrence in a path
//    /// The path MUST be nonempty.
//    occurrence_handle_t get_last_occurrence(const path_handle_t& path_handle) const;
//
//    /// Returns true if the occurrence is not the last occurence on the path, else false
//    bool has_next_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns true if the occurrence is not the first occurence on the path, else false
//    bool has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns a handle to the next occurrence on the path
//    occurrence_handle_t get_next_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns a handle to the previous occurrence on the path
//    occurrence_handle_t get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns a handle to the path that an occurrence is on
//    path_handle_t get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    /// Returns the 0-based ordinal rank of a occurrence on a path
//    size_t get_ordinal_rank_of_occurrence(const occurrence_handle_t& occurrence_handle) const;
//
//    ////////////////////////////////////////////////////////////////////////////
//    // Additional optional interface with a default implementation
//    ////////////////////////////////////////////////////////////////////////////
//
//    /// Returns true if the given path is empty, and false otherwise
//    bool is_empty(const path_handle_t& path_handle) const;


/**
 * This is the interface for a handle graph that supports modification.
 */
    /*
     * Note: All operations may invalidate path handles and occurrence handles.
     */
    
    /// Create a new node with the given sequence and return the handle.
    virtual handle_t create_handle(const std::string& sequence);

    /// Create a new node with the given id and sequence, then return the handle.
    virtual handle_t create_handle(const std::string& sequence, const id_t& id);
    
    /// Remove the node belonging to the given handle and all of its edges.
    /// Does not update any stored paths.
    /// Invalidates the destroyed handle.
    /// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
    /// May **NOT** be called during parallel for_each_handle iteration.
    /// May **NOT** be called on the node from which edges are being followed during follow_edges.
    virtual void destroy_handle(const handle_t& handle);
    
    /// Create an edge connecting the given handles in the given order and orientations.
    /// Ignores existing edges.
    virtual void create_edge(const handle_t& left, const handle_t& right);
    
    /// Remove the edge connecting the given handles in the given order and orientations.
    /// Ignores nonexistent edges.
    /// Does not update any stored paths.
    virtual void destroy_edge(const handle_t& left, const handle_t& right);
    
    /// Remove all nodes and edges. Does not update any stored paths.
    virtual void clear(void);
    
    /// TODO: swap_handles is actually not yet supported
    
    /// Swap the nodes corresponding to the given handles, in the ordering used
    /// by for_each_handle when looping over the graph. Other handles to the
    /// nodes being swapped must not be invalidated. If a swap is made while
    /// for_each_handle is running, it affects the order of the handles
    /// traversed during the current traversal (so swapping an already seen
    /// handle to a later handle's position will make the seen handle be visited
    /// again and the later handle not be visited at all).
    virtual void swap_handles(const handle_t& a, const handle_t& b);
    
    /// Alter the node that the given handle corresponds to so the orientation
    /// indicated by the handle becomes the node's local forward orientation.
    /// Rewrites all edges pointing to the node and the node's sequence to
    /// reflect this. Invalidates all handles to the node (including the one
    /// passed). Returns a new, valid handle to the node in its new forward
    /// orientation. Note that it is possible for the node's ID to change.
    /// Does not update any stored paths. May change the ordering of the underlying
    /// graph.
    virtual handle_t apply_orientation(const handle_t& handle);
    
    /// Split a handle's underlying node at the given offsets in the handle's
    /// orientation. Returns all of the handles to the parts. Other handles to
    /// the node being split may be invalidated. The split pieces stay in the
    /// same local forward orientation as the original node, but the returned
    /// handles come in the order and orientation appropriate for the handle
    /// passed in.
    /// Updates stored paths.
    virtual vector<handle_t> divide_handle(const handle_t& handle, const std::vector<size_t>& offsets);

/**
 * This is the interface for a handle graph with embedded paths where the paths can be modified.
 * Note that if the *graph* can also be modified, the implementation will also
 * need to inherit from MutableHandleGraph, via the combination
 * MutablePathMutableHandleGraph interface.
 * TODO: This is a very limited interface at the moment. It will probably need to be extended.
 */
    
//    /**
//     * Destroy the given path. Invalidates handles to the path and its node occurrences.
//     */
//    void destroy_path(const path_handle_t& path);
//
//    /**
//     * Create a path with the given name. The caller must ensure that no path
//     * with the given name exists already, or the behavior is undefined.
//     * Returns a handle to the created empty path. Handles to other paths must
//     * remain valid.
//     */
//    path_handle_t create_path_handle(const std::string& name);
//
//    /**
//     * Append a visit to a node to the given path. Returns a handle to the new
//     * final occurrence on the path which is appended. Handles to prior
//     * occurrences on the path, and to other paths, must remain valid.
//     */
//    occurrence_handle_t append_occurrence(const path_handle_t& path, const handle_t& to_append);

/// These are the backing data structures that we use to fulfill the above functions

private:
    // TODO: delete this later, very duplicative
    inline std::string reverse_complement(const std::string& seq) const;
    
    id_t max_id = 0;
    id_t min_id = std::numeric_limits<id_t>::max();
    
    size_t new_node_record(id_t node_id);
    inline uint8_t encode_nucleotide(const char& nt) const;
    inline char decode_nucleotide(const uint64_t& val) const;
    inline uint64_t complement_encoded_nucleotide(const uint64_t& val) const;
    inline size_t graph_iv_index(const handle_t& handle) const;
    inline uint64_t graph_index_to_seq_len_index(const size_t& graph_index) const;
    inline uint64_t graph_index_to_seq_start_index(const size_t& graph_index) const;
    inline const uint64_t& encode_edge_target(const handle_t& handle) const;
    inline const handle_t& decode_edge_target(const uint64_t& val) const;
    inline const uint64_t get_next_edge_index(const uint64_t& edge_index) const;
    inline const uint64_t get_edge_target(const uint64_t& edge_index) const;
    inline void set_edge_target(const uint64_t& edge_index, const handle_t& handle);
    void remove_edge_reference(const handle_t& on, const handle_t& to);
    void defragment(void);
    
    const static size_t PAGE_WIDTH = 64;
    
    // TODO: some of these offsets are a little silly and
    // only are around as legacy. I should remove them once
    // the factoring stabilizes
    
    const static size_t GRAPH_RECORD_SIZE = 2;
    const static size_t SEQ_START_RECORD_SIZE = 1;
    const static size_t SEQ_LENGTH_RECORD_SIZE = 1;
    
    const static size_t GRAPH_START_EDGES_OFFSET = 0;
    const static size_t GRAPH_END_EDGES_OFFSET = 1;
    
    const static size_t EDGE_RECORD_SIZE = 2;
    
    const static size_t EDGE_TRAV_OFFSET = 0;
    const static size_t EDGE_NEXT_OFFSET = 1;
    
    // defragment when the deleted records are this fraction of the whole
    const static double defrag_factor;
    
    /// Encodes the topology of the graph. Consists of fixed width records that represent
    /// offsets in edge_lists_iv.
    /// {start edge list index, end edge list index}
    PagedVector graph_iv;
    
    /// Encodes the start of a node's sequence in seq_iv. Matches the order of graph_iv.
    PagedVector seq_start_iv;
    
    /// Encodes the length of a node's sequence in seq_iv. Matches the order of graph_iv.
    PackedVector seq_length_iv;

    /// Encodes a series of edges lists of nodes.
    /// {ID|orientation (bit-packed), next edge index}
    PagedVector edge_lists_iv;
    
    /// Encodes the 1-based offset of an ID in graph_iv in units of GRAPH_RECORD_SIZE.
    /// If no node with that ID exists, contains a 0.
    PackedDeque id_to_graph_iv;

    /// Encodes all of the sequences of all nodes and all paths in the graph.
    /// The node sequences occur in the same order as in graph_iv;
    PackedVector seq_iv;

//    /// Same length as seq_iv. 0's indicate that a base is still touched by some
//    /// node or some path. 1's indicate that all nodes or paths that touch this
//    /// base have been deleted.
//    PackedVector dead_bv;
//
//    /// Encodes a self-balancing binary tree as integers. Consists of fixed-width
//    /// records that have the following structure:
//    /// {interval start, members index, parent index, left child index, right child index}
//    /// Interval start variable indicates the start of a range in seq_iv (corresponding to
//    /// a node, unless the node has been deleted), members index indicates the 1-based index
//    /// of the first path membership record corresponding to this interval in
//    /// path_membership_value_iv, and parent/left child/right child index indicates the
//    /// topology of a binary tree search structure for these intervals. The indexes are 1-based
//    /// with 0 indicating that the neighbor does not exist.
//    SuccinctSplayTree path_membership_range_iv;
//
//    /// Encodes a series of linked lists. Consists of fixed-width records that have
//    /// the following structure:
//    /// {path id, rank, next index}
//    /// Path ID indicates which path the node occurs on, rank indicates the ordinal
//    /// position of this occurrence in the path, and next index indicates the 1-based
//    /// index of the next occurrence of this node in this vector (or 0 if there is none)
//    PackedVector path_membership_value_iv;
//
//    /// Encodes the embedded paths of the graph. Each path is represented as three vectors
//    /// starts, lengths, orientations
//    /// The values in starts correspond to the 0-based indexes of an interval in seq_iv.
//    /// The values in lengths are simply the length.
//    /// The strand of this interval is given by the corresponding bit in orientations, with 1
//    /// indicating reverse strand.
//    std::vector<path_t> paths;
//
//    size_t dead_bases = 0;
    size_t deleted_node_records = 0;
    size_t deleted_edge_records = 0;
};
    
    
    
inline uint8_t PackedGraph::encode_nucleotide(const char& nt) const {
    if (nt == 'a' || nt == 'A') {
        return 0;
    }
    else if (nt == 'c' || nt == 'C') {
        return 1;
    }
    else if (nt == 'g' || nt == 'G') {
        return 2;
    }
    else if (nt == 't' || nt == 'T') {
        return 3;
    }
    else {
        // all others, but probably N's
        return 4;
    }
}
    
inline uint64_t PackedGraph::complement_encoded_nucleotide(const uint64_t& val) const {
    return val == 4 ? 4 : 3 - val;
}
    
inline char PackedGraph::decode_nucleotide(const uint64_t& val) const {
    static const char* alphabet = "ACGTN";
    return alphabet[val];
}
    
inline size_t PackedGraph::graph_iv_index(const handle_t& handle) const {
    return (id_to_graph_iv.get(get_id(handle) - min_id) - 1) * GRAPH_RECORD_SIZE;
}

inline uint64_t PackedGraph::graph_index_to_seq_len_index(const size_t& graph_index) const {
    return (graph_index * SEQ_LENGTH_RECORD_SIZE) / GRAPH_RECORD_SIZE;
}

inline uint64_t PackedGraph::graph_index_to_seq_start_index(const size_t& graph_index) const {
    return (graph_index * SEQ_START_RECORD_SIZE) / GRAPH_RECORD_SIZE;
}
    
inline const uint64_t& PackedGraph::encode_edge_target(const handle_t& handle) const {
    return reinterpret_cast<const uint64_t&>(handle);
}
    
inline const handle_t& PackedGraph::decode_edge_target(const uint64_t& val) const {
    return reinterpret_cast<const handle_t&>(val);
}

inline const uint64_t PackedGraph::get_next_edge_index(const uint64_t& edge_index) const {
    return edge_lists_iv.get((edge_index - 1) * EDGE_RECORD_SIZE + EDGE_NEXT_OFFSET);
}

inline const uint64_t PackedGraph::get_edge_target(const uint64_t& edge_index) const {
    return edge_lists_iv.get((edge_index - 1) * EDGE_RECORD_SIZE + EDGE_TRAV_OFFSET);
}
    
inline void PackedGraph::set_edge_target(const uint64_t& edge_index, const handle_t& handle) {
    edge_lists_iv.set((edge_index - 1) * EDGE_RECORD_SIZE + EDGE_TRAV_OFFSET, encode_edge_target(handle));
}

inline std::string PackedGraph::reverse_complement(const std::string& seq) const {
    std::string rev_comp(seq.size(), 'A');
    for (size_t i = 0; i < seq.size(); i++) {
        char nt = seq.at(i);
        if (nt == 'a' || nt == 'A') {
            rev_comp[rev_comp.size() - i - 1] = 'T';
        }
        else if (nt == 'c' || nt == 'C') {
            rev_comp[rev_comp.size() - i - 1] = 'G';
        }
        else if (nt == 'g' || nt == 'G') {
            rev_comp[rev_comp.size() - i - 1] = 'C';
        }
        else if (nt == 't' || nt == 'T') {
            rev_comp[rev_comp.size() - i - 1] = 'A';
        }
        else {
            rev_comp[rev_comp.size() - i - 1] = 'N';
        }
    }
    return rev_comp;
}

} // end dankness

#endif /* dgraph_hpp */
