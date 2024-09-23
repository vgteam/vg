#ifndef VG_SUBPATH_OVERLAY_HPP_INCLUDED
#define VG_SUBPATH_OVERLAY_HPP_INCLUDED

/**
 * \file subpath_overlay.hpp
 *
 * Provides SubpathOverlay, a PathHandleGraph overlay that presents an subinterval
 * of a path as a linearized graph.
 *
 */


#include "handle.hpp"

#include <handlegraph/util.hpp>


namespace vg {

using namespace handlegraph;

/**
 * Overlay that presents a linear graph corresponding to a subinterval of a path
 */
class SubpathOverlay : public ExpandingOverlayGraph {

public:
    /**
     * Construct new SubpathOverlay between two steps. The 'end' step is past-the-last.
     */
    SubpathOverlay(const PathHandleGraph* backing,
                   const step_handle_t& begin, const step_handle_t& end);
    SubpathOverlay() = default;

    ~SubpathOverlay() = default;

    ////////////////////////////////////////////////////////////////////////////
    // Handle-based interface
    ////////////////////////////////////////////////////////////////////////////

    /// Method to check if a node exists by ID
    bool has_node(nid_t node_id) const;
   
    /// Look up the handle for the node with the given ID in the given orientation
    handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    nid_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    std::string get_sequence(const handle_t& handle) const;
    
    /// Return the number of nodes in the graph
    size_t get_node_count() const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    nid_t min_node_id() const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    nid_t max_node_id() const;

    ///////////////////////////////////
    /// ExpandingOverlayGraph interface
    ///////////////////////////////////
    
    /**
     * Returns the handle in the underlying graph that corresponds to a handle in the
     * overlay
     */
    handle_t get_underlying_handle(const handle_t& handle) const;
    
protected:
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined. Returns true if we finished and false if we 
    /// stopped early.
    bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;

    /// the backing graph
    const HandleGraph* backing_graph = nullptr;

    std::vector<handle_t> subpath_handles;
};

}

#endif // VG_SUBPATH_OVERLAY_HPP_INCLUDED
