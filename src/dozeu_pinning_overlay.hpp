#ifndef VG_DOZEU_PINNING_OVERLAY_HPP_INCLUDED
#define VG_DOZEU_PINNING_OVERLAY_HPP_INCLUDED

/** \file
 * dozeu_pinning_overlay.hpp: defines a handle graph implementation that handles nodes
 * with no sequence in a way that plays well with dozeu pinned alignment
 */

#include "handle.hpp"

namespace vg {

using namespace std;

/*
 * An overlay that 1) masks source/sink nodes that have no sequence
 * associated with them and 2) duplicates any nodes necessary to
 * ensure that their neighbors are still sourc/sink nodes. This
 * transformation is invalid on graphs that have source-to-sink path
 * that consists entirely of null nodes. Assumes graph is single-
 * stranded.
 */
class DozeuPinningOverlay : public ExpandingOverlayGraph {
public:
    
    /// Initialize with the graph we want to do pinned alignment to. Boolean flag indicates
    /// whether sources or sinks should be preserved (for pinning right or left respectively).
    DozeuPinningOverlay(const HandleGraph* graph, bool preserve_sinks);
    
    /// Default constructor -- not actually functional
    DozeuPinningOverlay() = default;
    
    /// Default destructor
    ~DozeuPinningOverlay() = default;
    
    /// Returns true if any node has been duplicated in the overlay, in which case it may
    /// be necessary to translate between the IDs of this overlay and the underlying graph.
    bool performed_duplications() const;
    
    //////////////////////////
    /// HandleGraph interface
    //////////////////////////
    
    // Method to check if a node exists by ID
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
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    string get_sequence(const handle_t& handle) const;
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined.
    bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// Return the number of nodes in the graph.
    size_t get_node_count() const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t min_node_id() const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    id_t max_node_id() const;
    
    ///////////////////////////////////
    /// ExpandingOverlayGraph interface
    ///////////////////////////////////
    
    /// Returns the handle in the underlying graph that corresponds to a handle in the
    /// overlay
    handle_t get_underlying_handle(const handle_t& handle) const;
    
private:
    
    bool is_a_duplicate_handle(const handle_t& handle) const;
    
    bool is_a_duplicate_id(const id_t& node_id) const;
    
    id_t get_underlying_id(const id_t& node_id) const;
    
    handle_t get_duplicate_handle(const handle_t& handle) const;
    
    /// The graph we're masking empty nodes in
    const HandleGraph* graph = nullptr;
    
    /// Are we duplicating to preserve sources or sinks?
    bool preserve_sinks;
    
    /// The total number of null nodes
    size_t num_null_nodes = 0;
    
    /// The minimum numeric value of a handle
    uint64_t max_handle = numeric_limits<uint64_t>::min();
    
    /// The length of the interval between the max and min numeric values of handles
    uint64_t handle_val_range = 0;
    
    unordered_set<handle_t> duplicated_handles;
};
}

#endif // VG_DOZEU_PINNING_OVERLAY_HPP_INCLUDED
