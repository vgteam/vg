#ifndef VG_PATH_SUBGRAPH_HPP_INCLUDED
#define VG_PATH_SUBGRAPH_HPP_INCLUDED

/** \file
 * path_subgraph.hpp: represents a subgraph defined by a path in another graph.
 */

#include "handle.hpp"
#include <handlegraph/expanding_overlay_graph.hpp>
#include <unordered_set>
#include <vector>

namespace vg {

using namespace std;

    /**
     * A HandleGraph implementation that represents a subgraph of another HandleGraph, defined by a Path.
     * The leading and trailing nodes are cut according to the Path.
     * Supports translation of other Paths from this graph into the base graph.
     *
     * Nodes are numbered 1 to n along the path; multiple visits on the path to the same backing node will become distinct.
     */
    class PathSubgraph : public handlegraph::ExpandingOverlayGraph {
    public:
        
        /// Create a PathSubgraph describing the subgraph of the given graph
        /// defined by the given path. The path must not be empty. The path
        /// can take part of the start and end nodes, but all mappings must be
        /// perfect matches, and all adjacent mappings must properly cross a real
        /// edge.
        PathSubgraph(const HandleGraph* base, const Path& path);
        
        /// Get a topological order very easily, since the path defines one.
        vector<handle_t> get_topological_order() const;
        
        //////////////////////////
        /// HandleGraph interface
        //////////////////////////
        
        /// Method to check if a node exists by ID
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
        
        /// Get the sequence of a node, presented in the handle's local forward
        /// orientation.
        virtual string get_sequence(const handle_t& handle) const;
        
    protected:
        /// Loop over all the handles to next/previous (right/left) nodes. Passes
        /// them to a callback which returns false to stop iterating and true to
        /// continue. Returns true if we finished and false if we stopped early.
        virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
        
        /// Loop over all the nodes in the graph in their local forward
        /// orientations, in their internal stored order. Stop if the iteratee
        /// returns false. Can be told to run in parallel, in which case stopping
        /// after a false return value is on a best-effort basis and iteration
        /// order is not defined.
        virtual bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
        
    public:
        /// Return the number of nodes in the graph
        virtual size_t get_node_count() const;
        
        /// Return the smallest ID in the graph, or some smaller number if the
        /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t min_node_id() const;
        
        /// Return the largest ID in the graph, or some larger number if the
        /// largest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t max_node_id() const;
        
        //////////////////////////
        /// ExpandingOverlayGraph interface
        //////////////////////////
        
        /// Get the handle in the backing graph that the given handle in this graph represents.
        virtual handle_t get_underlying_handle(const handle_t& handle) const;

        //////////////////////////
        /// Additional Interface
        //////////////////////////

        /// Translate a Path against us to a Path against the base graph
        Path translate_down(const Path& path_against_subgraph) const;
        
    private:
        const HandleGraph* super = nullptr;
        Path defining_path;
    };
}

#endif
