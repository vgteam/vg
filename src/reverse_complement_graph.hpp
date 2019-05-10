#ifndef VG_REVERSE_COMPLEMENT_GRAPH_HPP_INCLUDED
#define VG_REVERSE_COMPLEMENT_GRAPH_HPP_INCLUDED

/** \file
 * reverse_complement_graph.hpp: defines a handle graph implementation that
 * reverses and complements the sequences of some other graph
 */

#include <unordered_map>
#include "handle.hpp"
#include "backwards_graph.hpp"

namespace vg {

using namespace std;

    /**
     * A HandleGraph implementation that wraps some other handle graph and reverses
     * but does *NOT* complement the sequences.
     *
     * See also: BackwardsGraph
     *
     * TODO: we inherit publically from BackwardsGraph, because inheriting it
     * privately and publically inheriting HandleGraph again produces a
     * warning. DO NOT use this as a BackwardsGraph; it does a different thing.
     */
    class ReverseComplementGraph : public BackwardsGraph {
    public:
        
        /// Initialize as the reverse complement version of another graph
        ReverseComplementGraph(const HandleGraph* forward_graph);
        
        /// Default constructor -- not actually functional
        ReverseComplementGraph() = default;
        
        /// Default destructor
        ~ReverseComplementGraph() = default;
        
        // The only thing we need to provide from the handle graph interface
        // that we don't get from BackwardsGraph is a complementing
        // get_sequence.
        
        /// Get the sequence of a node, presented in the handle's local forward
        /// orientation.
        virtual string get_sequence(const handle_t& handle) const;
    };
}

#endif
