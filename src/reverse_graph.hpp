#ifndef VG_REVERSE_GRAPH_HPP_INCLUDED
#define VG_REVERSE_GRAPH_HPP_INCLUDED

/** \file
 * reverse_graph.hpp: defines a handle graph implementation that reverses the sequences
 * of some other graph
 */

#include "handle.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

    /**
     * A HandleGraph implementation that wraps some other handle graph and reverses
     * and optionally complements the sequences.
     */
    class ReverseGraph : public ExpandingOverlayGraph {
    public:
        
        /// Initialize as the reverse version of another graph, optionally also
        /// complementing
        ReverseGraph(const HandleGraph* forward_graph, bool complement);
        
        /// Default constructor -- not actually functional
        ReverseGraph() = default;
        
        /// Default destructor
        ~ReverseGraph() = default;
        
        //////////////////////////
        /// HandleGraph interface
        //////////////////////////
        
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
        
        /// Get the sequence of a node, presented in the handle's local forward
        /// orientation.
        virtual string get_sequence(const handle_t& handle) const;
        
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
        
        /// Return the number of nodes in the graph
        virtual size_t get_node_count() const;
        
        /// Return the smallest ID in the graph, or some smaller number if the
        /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t min_node_id() const;
        
        /// Return the largest ID in the graph, or some larger number if the
        /// largest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t max_node_id() const;
        
        ////////////////////////////////////////////////////////////////////////////
        /// (Future) Overlay Interface
        ////////////////////////////////////////////////////////////////////////////
        
        /// Convert a backing graph handle to our handle to the same node
        inline handle_t from_backing(const handle_t& backing_handle) const {
            return backing_handle;
        }
        
    protected:
        
        ///////////////////////////////////
        /// ExpandingOverlayGraph interface
        ///////////////////////////////////
        
        /**
         * Returns the handle in the underlying graph that corresponds to a handle in the
         * overlay
         */
        virtual handle_t get_underlying_handle(const handle_t& handle) const;
        
        /// The forward version of the graph we're making backwards
        const HandleGraph* forward_graph = nullptr;
        
        /// Complement the sequences?
        bool complement = false;
    };
}

#endif
