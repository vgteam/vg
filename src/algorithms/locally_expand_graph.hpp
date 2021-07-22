#ifndef VG_ALGORITHMS_LOCALLY_EXPAND_GRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_LOCALLY_EXPAND_GRAPH_HPP_INCLUDED

/**
 * \file locally_expand_graph.hpp
 *
 * Definitions for the locally_expand_graph algorithm.
 */

#include "../handle.hpp"

namespace vg {
namespace algorithms {
    
/// Add to a subgraph (with same node ID space as parent) by walking forward from a given
/// node and adding all walks up to a maximum distance away. The handle provided graph
/// should be from the subgraph, not the parent graph.
void locally_expand_graph(const HandleGraph& parent, MutableHandleGraph& subgraph,
                          handle_t from, int64_t max_dist);
                                                      
}
}

#endif
