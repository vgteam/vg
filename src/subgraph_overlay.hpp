#ifndef VG_SUBGRAPH_OVERLAY_HPP_INCLUDED
#define VG_SUBGRAPH_OVERLAY_HPP_INCLUDED

/**
 * \file subgraph_overlay.hpp
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
 * Present a HandleGraph that is a backing HandleGraph but restricted
 * to a subset of nodes.  It won't give handles to nodes not in the 
 * subset, but it's not bulletproof: handles from outside the subset
 * won't undergo any special checks.  
 */
class SubgraphOverlay : virtual public HandleGraph {

public:
    /**
     * Make a new PathSubgraphOverlay. The backing graph must not be modified
     * while the overlay exists.
     *
     */
    SubgraphOverlay(const HandleGraph* backing, const unordered_set<nid_t>* node_subset);

    virtual ~SubgraphOverlay();

    ////////////////////////////////////////////////////////////////////////////
    // Handle-based interface
    ////////////////////////////////////////////////////////////////////////////

    /// Method to check if a node exists by ID
    virtual bool has_node(nid_t node_id) const;
   
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
    
    /// Get the ID from a handle
    virtual nid_t get_id(const handle_t& handle) const;
    
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const;
    
    /// Invert the orientation of a handle (potentially without getting its ID)
    virtual handle_t flip(const handle_t& handle) const;
    
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const;
    
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual std::string get_sequence(const handle_t& handle) const;
    
    /// Return the number of nodes in the graph
    virtual size_t get_node_count() const;
    
    /// Return the smallest ID in the graph, or some smaller number if the
    /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual nid_t min_node_id() const;
    
    /// Return the largest ID in the graph, or some larger number if the
    /// largest ID is unavailable. Return value is unspecified if the graph is empty.
    virtual nid_t max_node_id() const;

protected:
    
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;
    
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee
    /// returns false. Can be told to run in parallel, in which case stopping
    /// after a false return value is on a best-effort basis and iteration
    /// order is not defined. Returns true if we finished and false if we 
    /// stopped early.
    virtual bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;

protected:

    /// the backing graph
    const HandleGraph* backing_graph;

    /// the node subset. note, we don't own this so its up to client to keep in scope,
    /// just like the backing graph
    const unordered_set<nid_t>* node_subset;
    
    /// keep min_node_id() and max_node_id() constant
    nid_t min_node = 0;
    nid_t max_node = 0;
};

/**
 * Present a PathHandleGraph that is a backing HandleGraph but restricted
 * to a subset of nodes.
 *
 * Warning: we don't yet have a subgraph interface.  So we only consider paths
 * from the backing graph that are fully contained in the subgraph.
 */
class PathSubgraphOverlay : virtual public SubgraphOverlay, virtual public PathHandleGraph  {

public:
    /**
     * Make a new PathSubgraphOverlay. The backing graph must not be modified
     * while the overlay exists.
     *
     */
    PathSubgraphOverlay(const PathHandleGraph* backing, const unordered_set<nid_t>* node_subset);

    virtual ~PathSubgraphOverlay();

    ////////////////////////////////////////////////////////////////////////////
    // Path handle interface
    ////////////////////////////////////////////////////////////////////////////
    
    /// Returns the number of paths stored in the graph
    virtual size_t get_path_count() const;
    
    /// Determine if a path name exists and is legal to get a path handle for.
    virtual bool has_path(const std::string& path_name) const;
    
    /// Look up the path handle for the given path name.
    /// The path with that name must exist.
    virtual path_handle_t get_path_handle(const std::string& path_name) const;
    
    /// Look up the name of a path from a handle to it
    virtual std::string get_path_name(const path_handle_t& path_handle) const;
    
    /// Look up whether a path is circular
    virtual bool get_is_circular(const path_handle_t& path_handle) const;
    
    /// Returns the number of node steps in the path
    virtual size_t get_step_count(const path_handle_t& path_handle) const;
    
    /// Get a node handle (node ID and orientation) from a handle to an step on a path
    virtual handle_t get_handle_of_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the path that an step is on
    virtual path_handle_t get_path_handle_of_step(const step_handle_t& step_handle) const;
    
    /// Get a handle to the first step, which will be an arbitrary step in a circular path
    /// that we consider "first" based on our construction of the path. If the path is empty,
    /// then the implementation must return the same value as path_end().
    virtual step_handle_t path_begin(const path_handle_t& path_handle) const;
    
    /// Get a handle to a fictitious position past the end of a path. This position is
    /// returned by get_next_step for the final step in a path in a non-circular path.
    /// Note: get_next_step will *NEVER* return this value for a circular path.
    virtual step_handle_t path_end(const path_handle_t& path_handle) const;
    
    /// Get a handle to the last step, which will be an arbitrary step in a circular path that
    /// we consider "last" based on our construction of the path. If the path is empty
    /// then the implementation must return the same value as path_front_end().
    virtual step_handle_t path_back(const path_handle_t& path_handle) const;
    
    /// Get a handle to a fictitious position before the beginning of a path. This position is
    /// return by get_previous_step for the first step in a path in a non-circular path.
    /// Note: get_previous_step will *NEVER* return this value for a circular path.
    virtual step_handle_t path_front_end(const path_handle_t& path_handle) const;

    /// Returns true if the step is not the last step in a non-circular path.
    virtual bool has_next_step(const step_handle_t& step_handle) const;

    /// Returns true if the step is not the first step in a non-circular path.
    virtual bool has_previous_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the next step on the path. If the given step is the final step
    /// of a non-circular path, this method has undefined behavior. In a circular path,
    /// the "last" step will loop around to the "first" step.
    virtual step_handle_t get_next_step(const step_handle_t& step_handle) const;
    
    /// Returns a handle to the previous step on the path. If the given step is the first
    /// step of a non-circular path, this method has undefined behavior. In a circular path,
    /// it will loop around from the "first" step (i.e. the one returned by path_begin) to
    /// the "last" step.
    virtual step_handle_t get_previous_step(const step_handle_t& step_handle) const;

protected:    

    /// Execute a function on each path in the graph. If it returns false, stop
    /// iteration. Returns true if we finished and false if we stopped early.
    virtual bool for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const;

    /// Execute a function on each step of a handle in any path. If it
    /// returns false, stop iteration. Returns true if we finished and false if
    /// we stopped early.
    virtual bool for_each_step_on_handle_impl(const handle_t& handle,
                                              const std::function<bool(const step_handle_t&)>& iteratee) const;


protected:

    /// the backing path graph, just to not have to bother with dynamic cast
    const PathHandleGraph* backing_path_graph;

    /// the subset of paths from the backing graph that are entirely contained within our subgraph
    unordered_set<path_handle_t> path_subset;
};

}

#endif
