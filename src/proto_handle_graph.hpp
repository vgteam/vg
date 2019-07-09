#ifndef VG_PROTO_GRAPH_HPP_INCLUDED
#define VG_PROTO_GRAPH_HPP_INCLUDED

/** \file
 * proto_handle_graph.hpp: a handle graph wrapper for an unindexed Protobuf graph
 */

#include "handle.hpp"

namespace vg {

using namespace std;

    /**
     * THIS IS NOT AN EFFICIENT GRAPH IMPLMENTATION. It is a HandleGraph wrapper for an
     * unindexed Protobuf graph, which is primarily intended to provide a shim that eases
     * the transition away from using Protobuf. In particular, it does not maintain
     * adjacency lists, so edge and node lookup are O(N). However, it does provide the
     * ability to look up nodes and edges by ordinal index in addition to the standard
     * HandleGraph interface.
     */
    class ProtoHandleGraph : public HandleGraph {
    public:
        
        /// Initialize as a wrapper for a Protobuf graph
        ProtoHandleGraph(const Graph* graph);
        
        /// Default constructor -- not actually functional
        ProtoHandleGraph() = default;
        
        /// Default destructor
        ~ProtoHandleGraph() = default;
        
        //////////////////////////////////////////
        /// Non-HandleGraph methods from Protobuf
        //////////////////////////////////////////

        /// Get a handle to the i-th node
        handle_t get_handle_by_index(const size_t& i) const;
        
        /// Returns the number of edges in the graph
        size_t edge_size() const;
        
        /// Get a handle to the i-th edge
        edge_t get_edge_by_index(const size_t& i) const;
        
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
        /// TODO: can't be node_count because XG has a field named node_count.
        virtual size_t get_node_count() const;
        
        /// Return the smallest ID in the graph, or some smaller number if the
        /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t min_node_id() const;
        
        /// Return the largest ID in the graph, or some larger number if the
        /// largest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t max_node_id() const;
        
        /// Loop over all edges in their canonical orientation (as returned by edge_handle) and
        /// execute an iteratee on each one. Can stop early by returning false from the iteratee.
        /// Early stopping may not be immediate if executing in parallel.
        /// TODO: Promote to HandleGraph interface with a default implementation.
        virtual bool for_each_edge(const function<bool(const edge_t&)>& iteratee, bool parallel = false) const;
        
    private:
        /// The Protobuf graph we're wrapping
        const Graph* graph = nullptr;
    };
}

#endif
