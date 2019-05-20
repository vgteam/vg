#ifndef VG_EXTRA_NODE_GRAPH_HPP_INCLUDED
#define VG_EXTRA_NODE_GRAPH_HPP_INCLUDED

/**
 * \file extra_node_graph.hpp
 *
 * Provides a HandleGraph implementation that can add an extra node and edges
 * to/from it on top of a backing HandleGraph.
 *
 */


#include "handle.hpp"

#include <handlegraph/util.hpp>


namespace vg {

using namespace handlegraph;

/**
 * Present a HandleGraph that is a backing HandleGraph with an additional node
 * and edges to/from it added.
 */
class ExtraNodeGraph : public HandleGraph {

public:
    /**
     * Make a new ExtraNodeGraph. The backing graph must not be modified
     * while the overlay exists.
     *
     * Creates a new handle with the given sequence, with edges from all the
     * edges_in handles to it, and edges from it to all the edges_out handles.
     *
     * Self loops on the new node are not supported.
     */
    ExtraNodeGraph(const HandleGraph* backing, const string& sequence, const vector<handle_t>& edges_in, const vector<handle_t>& edges_out);
    
    /// Expose the handle to the new extra node
    handle_t get_created_handle() const;
    
    ////////////////////////////////////////////////////////////////////////////
    // Handle-based interface
    ////////////////////////////////////////////////////////////////////////////

    /// Check if a node exists by ID
    virtual bool has_node(id_t node_id) const;
    
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse) const;
    
    /// Get the ID from a handle
    virtual id_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    virtual handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee returns false.
    virtual bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph
    virtual size_t get_node_count() const;
    
    /// Return the smallest ID in the graph.
    virtual id_t min_node_id() const;
    
    /// Return the largest ID in the graph.
    virtual id_t max_node_id() const;
    
    /// Compute the degree of one side of a handle in O(1) time, if the backing
    /// graph also provides this facility in O(1) time. Takes O(n) time
    /// otherwise in the returned degree.
    virtual size_t get_degree(const handle_t& handle, bool go_left) const;
    
    ////////////////////////////////////////////////////////////////////////////
    // (Future) Overlay Interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Convert a backing graph handle to our handle to the same node
    inline handle_t from_backing(const handle_t& backing_handle) const {
        return as_handle(as_integer(backing_handle) + 2);
    }
    
protected:

    // TODO: a lot of this code can be unified with SourceSinkOverlay
    
    /// What backing graph do we overlay?
    const HandleGraph* backing;
    
    /// What are the handles in the backing graph that read into our new node forward?
    unordered_set<handle_t> in_from;
    /// And where do we go after our new node forward?
    unordered_set<handle_t> out_to;
    
    /// What is our projected node ID?
    id_t added_id;
    // And sequence?
    string sequence;
    
    // We reserve the 2 low numbers of the handles for our new node, and shift everything else up.
    const handle_t added_fwd = as_handle(0);
    const handle_t added_rev = as_handle(1);
    
    /// Convert our handle to a backing graph node into a backing graph handle to the same node
    inline handle_t to_backing(const handle_t& our_handle) const {
        return as_handle(as_integer(our_handle) - 2);
    }
    
    /// Determine if a handle points to an overlay-added node or not
    inline bool is_ours(const handle_t& our_handle) const {
        return ((uint64_t) as_integer(our_handle)) < 2;
    }
    
};


}

#endif
