//
//  hash_graph.hpp
//  
//  Contains a implementation of a sequence graph based on memory efficient
//  hash tables
//

#ifndef VG_HASH_GRAPH_HPP_INCLUDED
#define VG_HASH_GRAPH_HPP_INCLUDED

#include "handle.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "endianness.hpp"

namespace vg {
    
using namespace std;

class HashGraph : public MutablePathDeletableHandleGraph {
        
public:
    
    HashGraph();
    ~HashGraph();
    
    /// Deserialize from a stream of data
    HashGraph(istream& in);
    
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
    
    /// Efficiently get the number of edges attached to one side of a handle.
    size_t get_degree(const handle_t& handle, bool go_left) const;
    
    /// Returns one base of a handle's sequence, in the orientation of the
    /// handle.
    char get_base(const handle_t& handle, size_t index) const;
    
    /// Returns a substring of a handle's sequence, in the orientation of the
    /// handle. If the indicated substring would extend beyond the end of the
    /// handle's sequence, the return value is truncated to the sequence's end.
    string get_subsequence(const handle_t& handle, size_t index, size_t size) const;

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
    
    /// Adjust the representation of the graph in memory to improve performance.
    /// Optionally, allow the node IDs to be reassigned to further improve
    /// performance.
    /// Note: Ideally, this method is called one time once there is expected to be
    /// few graph modifications in the future.
    void optimize(bool allow_id_reassignment = true);
    
    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Returns the number of paths stored in the graph
    size_t get_path_count() const;

    /// Determine if a path name exists and is legal to get a path handle for.
    bool has_path(const std::string& path_name) const;

    /// Look up the path handle for the given path name.
    /// The path with that name must exist.
    path_handle_t get_path_handle(const std::string& path_name) const;

    /// Look up the name of a path from a handle to it
    string get_path_name(const path_handle_t& path_handle) const;
    
    /// Look up whether a path is circular
    bool get_is_circular(const path_handle_t& path_handle) const;

    /// Returns the number of node steps in the path
    size_t get_step_count(const path_handle_t& path_handle) const;

    /// Get a node handle (node ID and orientation) from a handle to a step on a path
    handle_t get_handle_of_step(const step_handle_t& step_handle) const;

    /// Returns a handle to the path that a step is on
    path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
    
    /// Get a handle to the first step, or in a circular path to an arbitrary step
    /// considered "first". If the path is empty, returns the past-the-last step
    /// returned by path_end.
    step_handle_t path_begin(const path_handle_t& path_handle) const;
    
    /// Get a handle to a fictitious position past the end of a path. This position is
    /// return by get_next_step for the final step in a path in a non-circular path.
    /// Note that get_next_step will *NEVER* return this value for a circular path.
    step_handle_t path_end(const path_handle_t& path_handle) const;
    
    /// Returns a handle to the next step on the path. If the given step is the final step
    /// of a non-circular path, returns the past-the-last step that is also returned by
    /// path_end. In a circular path, the "last" step will loop around to the "first" (i.e.
    /// the one returned by path_begin).
    /// Note: to iterate over each step one time, even in a circular path, consider
    /// for_each_step_in_path.
    step_handle_t get_next_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the previous step on the path. If the given step is the first
    /// step of a non-circular path, this method has undefined behavior. In a circular path,
    /// it will loop around from the "first" step (i.e. the one returned by path_begin) to
    /// the "last" step.
    /// Note: to iterate over each step one time, even in a circular path, consider
    /// for_each_step_in_path.
    step_handle_t get_previous_step(const step_handle_t& step_handle) const;
    
    /// Execute a function on each path in the graph
    bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;
    
    /// Calls a function with all steps of a node on paths.
    bool for_each_step_on_handle_impl(const handle_t& handle,
                                      const function<bool(const step_handle_t&)>& iteratee) const;
    
    /**
     * Destroy the given path. Invalidates handles to the path and its node steps.
     */
    void destroy_path(const path_handle_t& path);

    /**
     * Create a path with the given name. The caller must ensure that no path
     * with the given name exists already, or the behavior is undefined.
     * Returns a handle to the created empty path. Handles to other paths must
     * remain valid.
     */
    path_handle_t create_path_handle(const string& name, bool is_circular = false);

    /**
     * Append a visit to a node to the given path. Returns a handle to the new
     * final step on the path which is appended. If the path is cirular, the new
     * step is placed between the steps considered "last" and "first" by the
     * method path_begin. Handles to prior steps on the path, and to other paths,
     * must remain valid.
     */
    step_handle_t append_step(const path_handle_t& path, const handle_t& to_append);

    /**
     * Prepend a visit to a node to the given path. Returns a handle to the new
     * first step on the path which is appended. If the path is cirular, the new
     * step is placed between the steps considered "last" and "first" by the
     * method path_begin. Handles to later steps on the path, and to other paths,
     * must remain valid.
     */
    step_handle_t prepend_step(const path_handle_t& path, const handle_t& to_prepend);
    
    /**
     * Delete a segment of a path and rewrite it as some other sequence of steps. Returns a pair
     * of step_handle_t's that indicate the range of the new segment in the path. The segment to
     * delete should be designated by the first and the past-the-last step handle.  If the step
     * that is returned by path_begin is deleted, path_begin will now return the first step from
     * the new segment or, in the case that the new segment is empty, segment_end.
     */
    pair<step_handle_t, step_handle_t> rewrite_segment(const step_handle_t& segment_begin,
                                                       const step_handle_t& segment_end,
                                                       const std::vector<handle_t>& new_segment);
    
    /**
     * Make a path circular or non-circular. If the path is becoming circular, the
     * last step is joined to the first step. If the path is becoming linear, the
     * step considered "last" is unjoined from the step considered "first" according
     * to the method path_begin.
     */
    void set_circularity(const path_handle_t& path, bool circular);
    
    ////////////////////////////////////////////////////////////////////////////
    // Non-handle methods
    ////////////////////////////////////////////////////////////////////////////
    
    /// Write the graph to an out stream.
    void serialize(ostream& out) const;
    
    /// Read the graph (in the format written by serialize()) from an in stream.
    void deserialize(istream& in);
    
private:
    
    
    /*
     * A linked list record representing a single node in an embedded path
     */
    struct path_mapping_t {
        path_mapping_t(const handle_t& handle,
                       const int64_t& path_id) : handle(handle), path_id(path_id) {}
        
        handle_t handle;
        int64_t path_id;
        struct path_mapping_t* prev = nullptr;
        struct path_mapping_t* next = nullptr;
    };
    
    /*
     * A node object with the sequence and its edge lists
     */
    struct node_t {
        node_t() {}
        node_t(const string& sequence) : sequence(sequence) {}
        string sequence;
        /// Adjacency list from the left side of the node
        vector<handle_t> left_edges;
        /// Adjacency list from the right side of the node
        vector<handle_t> right_edges;
        /// The occurrences of this node on paths;
        vector<path_mapping_t*> occurrences;
        
        /// Write the node to an out stream.
        void serialize(ostream& out) const;
        /// Read the node (in the format written by serialize()) from an in stream.
        void deserialize(istream& in);
    };
    
    /*
     * A simple linked list implementation of an embedded path
     */
    class path_t {
    public:
        
        path_t();
        path_t(const string& name, const int64_t& path_id, bool is_circular = false);
        
        /// Move constructor
        path_t(path_t&& other);
        
        /// Move assignment
        path_t& operator=(path_t&& other);
        
        /// Copy constructor
        path_t(const path_t& other);
        
        /// Copy assignment
        path_t& operator=(const path_t& other);
        
        ~path_t();
        
        /// Add a node to the end of the path
        path_mapping_t* push_back(const handle_t& handle);
        
        /// Add a node to the front of the path
        path_mapping_t* push_front(const handle_t& handle);
        
        /// Remove the mapping from the path and free its memory
        void remove(path_mapping_t* mapping);
        
        /// Insert a new node into the middle of the path. If the the provided node
        /// is null, inserts at the end.
        path_mapping_t* insert_before(const handle_t& handle, path_mapping_t* mapping);
        
        /// Write the path to an out stream.
        void serialize(ostream& out) const;
        
        /// Read the path (in the format written by serialize()) from an in stream.
        void deserialize(istream& in);
        
        path_mapping_t* head = nullptr;
        path_mapping_t* tail = nullptr;
        size_t count = 0;
        int64_t path_id = 0;
        string name;
        bool is_circular = false;
    };
    
    /// The maximum ID in the graph
    id_t max_id = 0;
    /// The minimum ID in the graph
    id_t min_id = numeric_limits<id_t>::max();
    
    /// Encodes the graph topology
    hash_map<id_t, node_t> graph;
    
    /// Maps path names to path IDs
    string_hash_map<string, int64_t> path_id;
    
    /// Maps path IDs to the actual paths
    hash_map<int64_t, path_t> paths;
    
    /// The next path ID we will assign to a new path
    int64_t next_path_id = 1;
};
    
    

} // end dankness

#endif
