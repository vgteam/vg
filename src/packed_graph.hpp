//
//  packed_graph.hpp
//  
//  Contains a implementation of a sequence graph based on bit-packed integer
//  vectors.
//

#ifndef VG_PACKED_GRAPH_HPP_INCLUDED
#define VG_PACKED_GRAPH_HPP_INCLUDED

#include <utility>

#include "handle.hpp"
#include "packed_structs.hpp"
#include "hash_map.hpp"
#include "utility.hpp"


namespace vg {
    
using namespace std;


class PackedGraph : public MutablePathDeletableHandleGraph {
        
public:
    PackedGraph();
    ~PackedGraph();
    
    /// Construct from a stream
    PackedGraph(istream& in);
    
    /// Output to a stream
    void serialize(ostream& out) const;
    
    /// Load contents from a stream and replace current contents
    void deserialize(istream& in);
    
    ////////////////////////////////////////////////////////////////////////////
    // HandleGraph interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Method to check if a node exists by ID
    bool has_node(id_t node_id) const;
    
    /// Look up the handle for the node with the given ID in the given orientation
    handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    id_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward orientation.
    string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    /// TODO: can't be node_count because XG has a field named node_count.
    size_t node_size(void) const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t min_node_id(void) const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t max_node_id(void) const;

    /// Create a new node with the given sequence and return the handle.
    handle_t create_handle(const std::string& sequence);

    /// Create a new node with the given id and sequence, then return the handle.
    handle_t create_handle(const std::string& sequence, const id_t& id);
    
    /// Remove the node belonging to the given handle and all of its edges.
    /// Does not update any stored paths.
    /// Invalidates the destroyed handle.
    /// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
    /// May **NOT** be called during parallel for_each_handle iteration.
    /// May **NOT** be called on the node from which edges are being followed during follow_edges.
    void destroy_handle(const handle_t& handle);
    
    /// Create an edge connecting the given handles in the given order and orientations.
    /// Ignores existing edges.
    void create_edge(const handle_t& left, const handle_t& right);
    
    /// Remove the edge connecting the given handles in the given order and orientations.
    /// Ignores nonexistent edges.
    /// Does not update any stored paths.
    void destroy_edge(const handle_t& left, const handle_t& right);
    
    /// Remove all nodes and edges. Does not update any stored paths.
    void clear(void);
    
    /// TODO: swap_handles is actually not yet supported
    
    /// Swap the nodes corresponding to the given handles, in the ordering used
    /// by for_each_handle when looping over the graph. Other handles to the
    /// nodes being swapped must not be invalidated. If a swap is made while
    /// for_each_handle is running, it affects the order of the handles
    /// traversed during the current traversal (so swapping an already seen
    /// handle to a later handle's position will make the seen handle be visited
    /// again and the later handle not be visited at all).
    void swap_handles(const handle_t& a, const handle_t& b);
    
    /// Alter the node that the given handle corresponds to so the orientation
    /// indicated by the handle becomes the node's local forward orientation.
    /// Rewrites all edges pointing to the node and the node's sequence to
    /// reflect this. Invalidates all handles to the node (including the one
    /// passed). Returns a new, valid handle to the node in its new forward
    /// orientation. Note that it is possible for the node's ID to change.
    /// Does not update any stored paths. May change the ordering of the underlying
    /// graph.
    handle_t apply_orientation(const handle_t& handle);
    
    /// Split a handle's underlying node at the given offsets in the handle's
    /// orientation. Returns all of the handles to the parts. Other handles to
    /// the node being split may be invalidated. The split pieces stay in the
    /// same local forward orientation as the original node, but the returned
    /// handles come in the order and orientation appropriate for the handle
    /// passed in.
    /// Updates stored paths.
    vector<handle_t> divide_handle(const handle_t& handle, const std::vector<size_t>& offsets);
    
    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface
    ////////////////////////////////////////////////////////////////////////////

    /// Determine if a path name exists and is legal to get a path handle for.
    bool has_path(const std::string& path_name) const;

    /// Look up the path handle for the given path name.
    /// The path with that name must exist.
    path_handle_t get_path_handle(const std::string& path_name) const;

    /// Look up the name of a path from a handle to it
    string get_path_name(const path_handle_t& path_handle) const;

    /// Returns the number of node occurrences in the path
    size_t get_occurrence_count(const path_handle_t& path_handle) const;

    /// Returns the number of paths stored in the graph
    size_t get_path_count() const;

    /// Execute a function on each path in the graph
    bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;

    /// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
    handle_t get_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Get a handle to the first occurrence in a path.
    /// The path MUST be nonempty.
    occurrence_handle_t get_first_occurrence(const path_handle_t& path_handle) const;

    /// Get a handle to the last occurrence in a path
    /// The path MUST be nonempty.
    occurrence_handle_t get_last_occurrence(const path_handle_t& path_handle) const;

    /// Returns true if the occurrence is not the last occurence on the path, else false
    bool has_next_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Returns true if the occurrence is not the first occurence on the path, else false
    bool has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Returns a handle to the next occurrence on the path
    occurrence_handle_t get_next_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Returns a handle to the previous occurrence on the path
    occurrence_handle_t get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Returns a handle to the path that an occurrence is on
    path_handle_t get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const;
    
    /// Calls the given function for each occurrence of the given handle on a path.
    bool for_each_occurrence_on_handle_impl(const handle_t& handle,
                                                    const function<bool(const occurrence_handle_t&)>& iteratee) const;
    
    /**
     * Destroy the given path. Invalidates handles to the path and its node occurrences.
     */
    void destroy_path(const path_handle_t& path);

    /**
     * Create a path with the given name. The caller must ensure that no path
     * with the given name exists already, or the behavior is undefined.
     * Returns a handle to the created empty path. Handles to other paths must
     * remain valid.
     */
    path_handle_t create_path_handle(const string& name);

    /**
     * Append a visit to a node to the given path. Returns a handle to the new
     * final occurrence on the path which is appended. Handles to prior
     * occurrences on the path, and to other paths, must remain valid.
     */
    occurrence_handle_t append_occurrence(const path_handle_t& path, const handle_t& to_append);

    ////////////////////////////////////////////////////////////////////////////
    /// Specialized PackedGraph methods
    ////////////////////////////////////////////////////////////////////////////
    
    /// Attempt to compress data into less memory, possibly using more memory temporarily
    /// (especially useful before serializing).
    void compactify(void);
    
private:
    
    /// The maximum ID in the graph
    id_t max_id = 0;
    /// The minimum ID in the graph
    id_t min_id = std::numeric_limits<id_t>::max();
    
    /// Initialize all of the data corresponding with a new node and return
    /// it's 1-based offset
    size_t new_node_record(id_t node_id);
    
    /// Find and edge on given handle, to a given handle, and remove it from the edge list
    void remove_edge_reference(const handle_t& on, const handle_t& to);
    
    /// Check if have orphaned enough records in the graph's various linked lists to
    /// warrant reallocating and defragmenting them. If so, do it. Optionally, defragment
    /// even if we have not deleted many things.
    void defragment(bool force = false);
    
    /// Defragment when the orphaned records are this fraction of the whole.
    const static double defrag_factor;
    
    /// We use a standard page width for all page-compressed vectors
    const static size_t PAGE_WIDTH = 128;
    
    // TODO: some of these offsets are a little silly and only are around as legacy.
    // They could be removed once the factoring stabilizes, but optimization will also
    // probably handle it.
    
    /// Encodes the topology of the graph. Consists of fixed width records that represent
    /// offsets in edge_lists_iv.
    /// {start edge list index, end edge list index}
    PagedVector graph_iv;
    const static size_t GRAPH_RECORD_SIZE = 2;
    const static size_t GRAPH_START_EDGES_OFFSET = 0;
    const static size_t GRAPH_END_EDGES_OFFSET = 1;
    
    /// Encodes the start of a node's sequence in seq_iv. Matches the order of graph_iv.
    PagedVector seq_start_iv;
    const static size_t SEQ_START_RECORD_SIZE = 1;
    
    /// Encodes the length of a node's sequence in seq_iv. Matches the order of graph_iv.
    PackedVector seq_length_iv;
    const static size_t SEQ_LENGTH_RECORD_SIZE = 1;

    /// Encodes a series of edges lists of nodes.
    /// {ID|orientation (bit-packed), next edge index}
    PagedVector edge_lists_iv;
    const static size_t EDGE_RECORD_SIZE = 2;
    const static size_t EDGE_TRAV_OFFSET = 0;
    const static size_t EDGE_NEXT_OFFSET = 1;
    
    // TODO: template out the deque and back id_to_graph_iv with a paged vector? might
    // provide better compression now that it can handle 0's gracefully. unsure how the
    // wrapping around would act with pages though...
    
    /// Encodes the 1-based offset of an ID in graph_iv in units of GRAPH_RECORD_SIZE.
    /// If no node with that ID exists, contains a 0. The index of a given ID is
    /// computed by (ID - min ID).
    PackedDeque id_to_graph_iv;

    /// Encodes all of the sequences of all nodes in the graph.
    PackedVector seq_iv;
    
    /// Encodes the membership of a node in all paths. In the same order as graph_iv.
    /// Consists of 1-based offset to the corresponding heads of linked lists in
    /// path_membership_value_iv, which contains the actual pointers into the paths.
    PagedVector path_membership_node_iv;
    const static size_t NODE_MEMBER_RECORD_SIZE = 1;
    
    /// Encodes a series of linked lists of the memberships within paths. Consists of
    /// fixed width records of the following form. Path IDs are 0-based indexes, the
    /// other two indexes are 1-based and expressed in units of PATH_RECORD_SIZE and
    /// MEMBERSHIP_RECORD_SIZE respectively.
    /// {path ID, index in path, next membership record index}
    PagedVector path_membership_value_iv;
    const static size_t MEMBERSHIP_RECORD_SIZE = 3;
    const static size_t MEMBERSHIP_PATH_OFFSET = 0;
    const static size_t MEMBERSHIP_OCCURRENCE_OFFSET = 1;
    const static size_t MEMBERSHIP_NEXT_OFFSET = 2;
    
    /*
     * A struct to package the data associated with a path through the graph.
     */
    struct PackedPath {
        PackedPath(const string& name) : name(name), occurrences_iv(PAGE_WIDTH) {}
        
        /// The path's name
        string name;
        
        /// Marks whether this path has been deleted
        bool is_deleted = false;
        
        // TODO: split up occurrences_iv into adjacencies and travs separately
        
        /// Linked list records that encode the oriented nodes of the path. Indexes are
        /// 1-based, with 0 used as a sentinel to indicate none further.
        /// {ID|orientation (bit-packed), prev occurrence index, next occurrence index}
        PagedVector occurrences_iv;
        
        /// 1-based index of the head of the linked list in occurrences_iv.
        size_t head = 0;
        
        /// 1-based index of the tail of the linked list in occurrences_iv.
        size_t tail = 0;
    };
    const static size_t PATH_RECORD_SIZE = 3;
    const static size_t PATH_TRAV_OFFSET = 0;
    const static size_t PATH_PREV_OFFSET = 1;
    const static size_t PATH_NEXT_OFFSET = 2;
    
    /// Map from path names to index in the paths vector.
    string_hash_map<string, int64_t> path_id;
    
    /// Vector of the embedded paths in the graph
    vector<PackedPath> paths;
    
    ///////////////////////////////////////////////////////////////////////
    /// Convenience functions to translate between encodings in the vectors
    ///////////////////////////////////////////////////////////////////////
    
    inline uint8_t encode_nucleotide(const char& nt) const;
    inline char decode_nucleotide(const uint64_t& val) const;
    inline uint64_t complement_encoded_nucleotide(const uint64_t& val) const;
    
    inline size_t graph_iv_index(const handle_t& handle) const;
    
    inline uint64_t graph_index_to_seq_len_index(const size_t& graph_index) const;
    inline uint64_t graph_index_to_seq_start_index(const size_t& graph_index) const;
    inline uint64_t graph_index_to_node_member_index(const size_t& graph_index) const;
    
    inline const uint64_t& encode_traversal(const handle_t& handle) const;
    inline const handle_t& decode_traversal(const uint64_t& val) const;
    
    inline uint64_t get_next_edge_index(const uint64_t& edge_index) const;
    inline uint64_t get_edge_target(const uint64_t& edge_index) const;
    inline void set_edge_target(const uint64_t& edge_index, const handle_t& handle);
    
    inline uint64_t get_next_membership(const uint64_t& membership_index) const;
    inline uint64_t get_membership_occurrence(const uint64_t& membership_index) const;
    inline uint64_t get_membership_path(const uint64_t& membership_index) const;
    inline void set_next_membership(const uint64_t& membership_index, const uint64_t& next);
    
    inline uint64_t get_occurrence_trav(const PackedPath& path, const uint64_t& occurrence_index) const;
    inline uint64_t get_occurrence_prev(const PackedPath& path, const uint64_t& occurrence_index) const;
    inline uint64_t get_occurrence_next(const PackedPath& path, const uint64_t& occurrence_index) const;
    inline void set_occurrence_trav(PackedPath& path, const uint64_t& occurrence_index, const uint64_t& trav);
    inline void set_occurrence_prev(PackedPath& path, const uint64_t& occurrence_index, const uint64_t& prev_index);
    inline void set_occurrence_next(PackedPath& path, const uint64_t& occurrence_index, const uint64_t& next_index);
    
    uint64_t deleted_node_records = 0;
    uint64_t deleted_edge_records = 0;
    uint64_t deleted_membership_records = 0;
};
    
    
    
inline uint8_t PackedGraph::encode_nucleotide(const char& nt) const {
    
    uint8_t encoded;
    switch (nt) {
        case 'a':
        case 'A':
            encoded = 0;
            break;
            
        case 'c':
        case 'C':
            encoded = 1;
            break;
            
        case 'g':
        case 'G':
            encoded = 2;
            break;
            
        case 't':
        case 'T':
            encoded = 3;
            break;
            
        default:
            // all others, but probably N's
            encoded = 4;
            break;
    }
    
    return encoded;
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

inline uint64_t PackedGraph::graph_index_to_node_member_index(const size_t& graph_index) const {
    return (graph_index * NODE_MEMBER_RECORD_SIZE) / GRAPH_RECORD_SIZE;
}
    
inline const uint64_t& PackedGraph::encode_traversal(const handle_t& handle) const {
    return reinterpret_cast<const uint64_t&>(handle);
}
    
inline const handle_t& PackedGraph::decode_traversal(const uint64_t& val) const {
    return reinterpret_cast<const handle_t&>(val);
}

inline uint64_t PackedGraph::get_next_edge_index(const uint64_t& edge_index) const {
    return edge_lists_iv.get((edge_index - 1) * EDGE_RECORD_SIZE + EDGE_NEXT_OFFSET);
}

inline uint64_t PackedGraph::get_edge_target(const uint64_t& edge_index) const {
    return edge_lists_iv.get((edge_index - 1) * EDGE_RECORD_SIZE + EDGE_TRAV_OFFSET);
}
    
inline void PackedGraph::set_edge_target(const uint64_t& edge_index, const handle_t& handle) {
    edge_lists_iv.set((edge_index - 1) * EDGE_RECORD_SIZE + EDGE_TRAV_OFFSET, encode_traversal(handle));
}
    
inline uint64_t PackedGraph::get_next_membership(const uint64_t& membership_index) const {
    return path_membership_value_iv.get((membership_index - 1) * MEMBERSHIP_RECORD_SIZE + MEMBERSHIP_NEXT_OFFSET);
}
    
inline uint64_t PackedGraph::get_membership_occurrence(const uint64_t& membership_index) const {
    return path_membership_value_iv.get((membership_index - 1) * MEMBERSHIP_RECORD_SIZE + MEMBERSHIP_OCCURRENCE_OFFSET);
}

inline uint64_t PackedGraph::get_membership_path(const uint64_t& membership_index) const {
    return path_membership_value_iv.get((membership_index - 1) * MEMBERSHIP_RECORD_SIZE + MEMBERSHIP_PATH_OFFSET);
}

inline void PackedGraph::set_next_membership(const uint64_t& membership_index, const uint64_t& next) {
    path_membership_value_iv.set((membership_index - 1) * MEMBERSHIP_RECORD_SIZE + MEMBERSHIP_NEXT_OFFSET, next);
}

inline uint64_t PackedGraph::get_occurrence_trav(const PackedPath& path, const uint64_t& occurrence_index) const {
    return path.occurrences_iv.get((occurrence_index - 1) * PATH_RECORD_SIZE + PATH_TRAV_OFFSET);
}

inline uint64_t PackedGraph::get_occurrence_prev(const PackedPath& path, const uint64_t& occurrence_index) const {
    return path.occurrences_iv.get((occurrence_index - 1) * PATH_RECORD_SIZE + PATH_PREV_OFFSET);
}

inline uint64_t PackedGraph::get_occurrence_next(const PackedPath& path, const uint64_t& occurrence_index) const {
    return path.occurrences_iv.get((occurrence_index - 1) * PATH_RECORD_SIZE + PATH_NEXT_OFFSET);
}

inline void PackedGraph::set_occurrence_trav(PackedPath& path, const uint64_t& occurrence_index, const uint64_t& trav) {
    path.occurrences_iv.set((occurrence_index - 1) * PATH_RECORD_SIZE + PATH_TRAV_OFFSET, trav);
}

inline void PackedGraph::set_occurrence_prev(PackedPath& path, const uint64_t& occurrence_index, const uint64_t& prev_index) {
    path.occurrences_iv.set((occurrence_index - 1) * PATH_RECORD_SIZE + PATH_PREV_OFFSET, prev_index);
}

inline void PackedGraph::set_occurrence_next(PackedPath& path, const uint64_t& occurrence_index, const uint64_t& next_index) {
    path.occurrences_iv.set((occurrence_index - 1) * PATH_RECORD_SIZE + PATH_NEXT_OFFSET, next_index);
}

} // end dankness

#endif /* dgraph_hpp */
