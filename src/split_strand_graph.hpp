#ifndef VG_SPLIT_STRAND_GRAPH_HPP_INCLUDED
#define VG_SPLIT_STRAND_GRAPH_HPP_INCLUDED

/** \file
 * split_strand_graph.hpp: defines a handle graph overlay that duplicates nodes
 * and edges so that both the forward and reverse strand of the underlying graph
 * are now on the forward strand
 */

#include "handle.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

    /**
     * A HandleGraph implementation that overlays some other handle graph and splits
     * the two strands of its nodes into separate nodes
     */
    class StrandSplitGraph : public ExpandingOverlayGraph {
    public:
        
        /// Initialize as the reverse version of another graph, optionally also
        /// complementing
        StrandSplitGraph(const HandleGraph* graph);
        
        /// Default constructor -- not actually functional
        StrandSplitGraph() = default;
        
        /// Default destructor
        ~StrandSplitGraph() = default;
        
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
        bool follow_edges_impl(const handle_t& handle, bool go_left,
                                       const function<bool(const handle_t&)>& iteratee) const;
        
        /// Loop over all the nodes in the graph in their local forward
        /// orientations, in their internal stored order. Stop if the iteratee
        /// returns false. Can be told to run in parallel, in which case stopping
        /// after a false return value is on a best-effort basis and iteration
        /// order is not defined.
        bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                          bool parallel = false) const;
        
        /// Return the number of nodes in the graph
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
        
        /**
         * Returns the handle in the underlying graph that corresponds to a handle in the
         * overlay
         */
        handle_t get_underlying_handle(const handle_t& handle) const;

        ///////////////////////////////////
        /// Extra methods
        ///////////////////////////////////
        
        /**
         * Returns true if any nodes in the overlay correspond to the given node in the underlying graph.
         */
        bool has_overlay_node_for(const nid_t& underlying_id) const;
        
        /**
         * Returns the handle in the overlay graph that corresponds to a handle
         * and orientation in the underlying graph. Reverse versions of
         * underlying graph nodes become the locally-forward overlay node that
         * represents them.
         */
        handle_t get_overlay_handle(const handle_t& underlying_handle) const;
        
    private:
        /// The underlying graph we're making splitting
        const HandleGraph* graph = nullptr;
    };
}

#endif
