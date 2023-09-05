#ifndef VG_SOURCE_SINK_OVERLAY_HPP_INCLUDED
#define VG_SOURCE_SINK_OVERLAY_HPP_INCLUDED

/**
 * \file source_sink_overlay.hpp
 *
 * Provides SourceSinkOverlay, a HandleGraph implementation that joins all the
 * heads and tails of a backing graph to single source and sink nodes. 
 *
 */


#include "handle.hpp"

#include <handlegraph/util.hpp>


namespace vg {

using namespace handlegraph;

/**
 * Present a HandleGraph that is a backing HandleGraph with all its head nodes
 * connected to a single source node, and all its tail nodes connected to a
 * single sink node.
 */
class SourceSinkOverlay : public ExpandingOverlayGraph {

public:
    /**
     * Make a new SourceSinkOverlay. The backing graph must not be modified
     * while the overlay exists.
     *
     * The overlay will project a source node consisting of '#' characters, and
     * a sink node consisting of '$' characters. The lengths of the nodes may
     * be specified, and default to 1024, the max length that GCSA2 supports.
     * The IDs of the nodes will be autodetected from the backing graph's max
     * ID if not specified (or given as 0). If either is specified, both must
     * be specified.
     *
     * Also breaks into disconnected components with no tips, unless
     * break_disconnected is false. When breaking into such a component, we
     * choose an arbitrary node, link the source node to its start, and link
     * everything that also went to its start to the sink node.
     */
    SourceSinkOverlay(const HandleGraph* backing, size_t length = 1024, id_t source_id = 0, id_t sink_id = 0,
        bool break_disconnected = true);
    
    /// Expose the handle to the synthetic source
    handle_t get_source_handle() const;
    
    /// Expose the handle to the synthetic sink
    handle_t get_sink_handle() const;
    
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
    // Overlay interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Get the handle in the underlying graph that corresponds to the handle in
    /// the overlay.
    /// Throws an error if called on either the source or sink node
    virtual handle_t get_underlying_handle(const handle_t& handle) const;
    
protected:
    
    /// How long are the projected nodes?
    size_t node_length;
    
    /// What backing graph do we overlay?
    const HandleGraph* backing;
    
    /// What is our projected source node ID?
    id_t source_id;
    /// What is our projected sink node ID?
    id_t sink_id;
    
    /// We keep a set of backing graph head handles, in backing graph handle
    /// space. This also includes anything else we need to hook up to our
    /// source node to break into tipless components.
    unordered_set<handle_t> backing_heads;
    /// And similarly for the tails. These handles read out of their components.
    unordered_set<handle_t> backing_tails;
    
    // We reserve the 4 low numbers of the handles for our new source and sink, and shift everything else up.
    // These could have been static, but I couldn't figure out how to properly initialize them using constexpr functions.
    const handle_t source_fwd = as_handle(0);
    const handle_t source_rev = as_handle(1);
    const handle_t sink_fwd = as_handle(2);
    const handle_t sink_rev = as_handle(3);
    
    /// Convert a backing graph handle to our handle to the same node
    inline handle_t from_backing(const handle_t& backing_handle) const {
        return as_handle(as_integer(backing_handle) + 4);
    }
    
    /// Convert our handle to a backing graph node into a backing graph handle to the same node
    inline handle_t to_backing(const handle_t& our_handle) const {
        return as_handle(as_integer(our_handle) - 4);
    }
    
    /// Determine if a handle points to an overlay-added node or not
    inline bool is_ours(const handle_t& our_handle) const {
        return ((uint64_t) as_integer(our_handle)) < 4;
    }
    
};


}

#endif
