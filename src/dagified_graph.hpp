#ifndef VG_DAGIFIED_GRAPH_HPP_INCLUDED
#define VG_DAGIFIED_GRAPH_HPP_INCLUDED

/** \file
 * dagified_graph.hpp: defines a handle graph overlay implementation transforms cyclic graphs into DAGs
 */

#include "handle.hpp"

namespace vg {

using namespace std;

    /**
     * A HandleGraph implementation that wraps some other handle graph and converts it into a
     * DAG, preserving all paths up a a minimum length.
     */
    class DagifiedGraph : public ExpandingOverlayGraph {
    public:
        
        /// Expand a single-stranded graph into a DAG, preserving all walks up to the minimum length.
        /// If max duplications is provided, limits the number of times any node is copied.
        DagifiedGraph(const HandleGraph* graph, size_t min_preserved_path_length, 
                      size_t max_num_duplications = std::numeric_limits<size_t>::max());
        
        /// Default constructor -- not actually functional
        DagifiedGraph() = default;
        
        /// Default destructor
        ~DagifiedGraph() = default;
        
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
        
    protected:
        
        /*
         * Helper methods
         */
        
        /// Helper function to identify which ordinal copy of a strongly connected component
        /// in the underlying graph the handle belongs to
        uint64_t scc_copy_of_handle(const handle_t& handle) const;
        
        /// Helper function to identify which ordinal copy of a strongly connected component
        /// in the underlying graph the node ID belongs to
        uint64_t scc_copy_of_node_id(const id_t& node_id) const;
        
        /// Helper function that returns the layout order of a handle from the underlying graph
        size_t layout_order_of_handle(const handle_t& handle) const;
        
        /// Helper function to identify the node ID in the original graph that corresponds to a
        /// node in the dagified graph
        id_t get_underlying_id(const id_t& node_id) const;
        
        /// Helper function, returns the n-th copy in the dagified graph of a handle in the underlying graph
        handle_t nth_copy_of_handle(const handle_t& handle, const uint64_t& n) const;
                
        /*
         * Member variables
         */
                
        /// The underlying graph we're dagifiying
        const HandleGraph* graph = nullptr;
        
        /// Map from a canonical orientation of underlying handles (not necessarily forward!)
        /// to the ordinal position of the node in a low-FAS layout
        unordered_map<handle_t, size_t> layout_order;
        
        /// The ID of the strongly connected component that the handle at each layout position
        /// belongs to
        vector<uint64_t> scc_of_handle;
        
        /// The number of times each strongly connected component is duplicated in the
        /// dagified graph
        vector<uint64_t> scc_copy_count;
        
        /// The minimum value of a handle in the underlying graph, interpreted as an integer
        uint64_t min_handle = std::numeric_limits<uint64_t>::max();
        
        /// The width of the range of values that handles in the underlying graph take
        uint64_t handle_val_range = 0;
        
        /// The number of nodes including duplicated nodes, computed during construction
        size_t node_count = 0;
        
        /// The maximum ID of the graph, computed during construction
        id_t max_id = std::numeric_limits<uint64_t>::min();
        
    };
}

#endif
