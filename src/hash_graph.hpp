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

namespace vg {
    
using namespace std;

class HashGraph : public MutablePathDeletableHandleGraph {
        
public:
    
    HashGraph();
    ~HashGraph();
    
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
    
    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface
    ////////////////////////////////////////////////////////////////////////////

    /// Determine if a path name exists and is legal to get a path handle for.
    virtual bool has_path(const std::string& path_name) const;

    /// Look up the path handle for the given path name.
    /// The path with that name must exist.
    virtual path_handle_t get_path_handle(const std::string& path_name) const;

    /// Look up the name of a path from a handle to it
    virtual string get_path_name(const path_handle_t& path_handle) const;

    /// Returns the number of node occurrences in the path
    virtual size_t get_occurrence_count(const path_handle_t& path_handle) const;

    /// Returns the number of paths stored in the graph
    virtual size_t get_path_count() const;

    /// Execute a function on each path in the graph
    // TODO: allow stopping early?
    virtual void for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const;

    /// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
    virtual handle_t get_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Get a handle to the first occurrence in a path.
    /// The path MUST be nonempty.
    virtual occurrence_handle_t get_first_occurrence(const path_handle_t& path_handle) const;

    /// Get a handle to the last occurrence in a path
    /// The path MUST be nonempty.
    virtual occurrence_handle_t get_last_occurrence(const path_handle_t& path_handle) const;

    /// Returns true if the occurrence is not the last occurence on the path, else false
    virtual bool has_next_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Returns true if the occurrence is not the first occurence on the path, else false
    virtual bool has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Returns a handle to the next occurrence on the path
    virtual occurrence_handle_t get_next_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Returns a handle to the previous occurrence on the path
    virtual occurrence_handle_t get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const;

    /// Returns a handle to the path that an occurrence is on
    virtual path_handle_t get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const;
    
    /// Returns a vector of all occurrences of a node on paths. Optionally restricts to
    /// occurrences that match the handle in orientation.
    virtual vector<occurrence_handle_t> occurrences_of_handle(const handle_t& handle,
                                                              bool match_orientation = false) const;
    
    /**
     * Destroy the given path. Invalidates handles to the path and its node occurrences.
     */
    virtual void destroy_path(const path_handle_t& path);

    /**
     * Create a path with the given name. The caller must ensure that no path
     * with the given name exists already, or the behavior is undefined.
     * Returns a handle to the created empty path. Handles to other paths must
     * remain valid.
     */
    virtual path_handle_t create_path_handle(const string& name);

    /**
     * Append a visit to a node to the given path. Returns a handle to the new
     * final occurrence on the path which is appended. Handles to prior
     * occurrences on the path, and to other paths, must remain valid.
     */
    virtual occurrence_handle_t append_occurrence(const path_handle_t& path, const handle_t& to_append);

private:
    
    /*
     * A node object with the sequence and its edge lists
     */
    struct node_t {
        node_t() {}
        node_t(const string& sequence) : sequence(sequence) {}
        string sequence;
        vector<handle_t> left_edges;
        vector<handle_t> right_edges;
    };
    
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
     * A simple linked list implementation of an embedded path
     */
    struct path_t {
        path_t() {}
        path_t(const string& name, const int64_t& path_id) : name(name), path_id(path_id) {}
        ~path_t() {
            for (path_mapping_t* mapping = head; mapping != nullptr;) {
                path_mapping_t* next = mapping->next;
                delete mapping;
                mapping = next;
            }
        }
        
        /// Add a node to the end of the path
        path_mapping_t* push_back(const handle_t& handle) {
            path_mapping_t* pushing = new path_mapping_t(handle, path_id);
            if (!tail) {
                tail = head = pushing;
            }
            else {
                tail->next = pushing;
                pushing->prev = tail;
                tail = pushing;
            }
            count++;
            return pushing;
        }
        
        /// Add a node to the front of the path
        path_mapping_t* push_front(const handle_t& handle) {
            path_mapping_t* pushing = new path_mapping_t(handle, path_id);
            if (!head) {
                tail = head = pushing;
            }
            else {
                head->prev = pushing;
                pushing->next = head;
                head = pushing;
            }
            count++;
            return pushing;
        }
        
        /// Remove the mapping from the path and free its memory
        void remove(path_mapping_t* mapping) {
            if (mapping == head) {
                head = mapping->next;
            }
            if (mapping == tail) {
                tail = mapping->prev;
            }
            if (mapping->next) {
                mapping->next->prev = mapping->prev;
            }
            if (mapping->prev) {
                mapping->prev->next = mapping->next;
            }
            count--;
            delete mapping;
        }
        
        /// Insert a new node into the middle of the path. If the the profvided node
        /// is null, inserts at the beginning.
        path_mapping_t* insert_after(const handle_t& handle, path_mapping_t* mapping) {
            path_mapping_t* inserting = new path_mapping_t(handle, path_id);
            if (mapping) {
                inserting->next = mapping->next;
                if (mapping->next) {
                    mapping->next->prev = inserting;
                }
                mapping->next = inserting;
                inserting->prev = mapping;
                
                if (mapping == tail) {
                    tail = inserting;
                }
            }
            else {
                if (head) {
                    inserting->next = head;
                    head->prev = inserting;
                    head = inserting;
                }
                else {
                    head = tail = inserting;
                }
            }
            
            count++;
            return inserting;
        }
        
        
        
        path_mapping_t* head = nullptr;
        path_mapping_t* tail = nullptr;
        size_t count = 0;
        int64_t path_id = 0;
        string name;
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
    
    /// Maps node IDs to their occurrences on paths
    hash_map<id_t, vector<path_mapping_t*>> occurrences;
    
    /// The next path ID we will assign to a new path
    int64_t next_path_id = 1;
};
    
    

} // end dankness

#endif
