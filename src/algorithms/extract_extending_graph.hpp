#ifndef VG_ALGORITHMS_EXTRACT_EXTENDING_GRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_EXTRACT_EXTENDING_GRAPH_HPP_INCLUDED

/**
 * \file extract_extending_graphs.hpp
 *
 * Definitions for the extract_extending_graph algorithm.
 */

#include <unordered_map>
#include <vg/vg.pb.h>

#include "../position.hpp"
#include "../hash_map.hpp"
#include "../handle.hpp"

namespace vg {
namespace algorithms {
    
    /// Fills graph 'into' with the subgraph of the handle graph 'source' that extends in one direction from
    /// a given position, up to a maximum distance. The node containing the position will be "cut" so that only
    /// the portion that is forward in the search direction remains. Node IDs may be changed in the extracted
    /// graph, but they can be translated back to node IDs in the original graph with the returned map, although
    /// that translation procedure MUST handle the node that pos is on specially, as it may be cut.
    /// translate_node_ids from path.hpp can do this as long as you pass along what part of the node was removed.
    /// The node containing the source position may optionally be duplicated to preserve cycles on it after its
    /// cut, but no other nodes will will duplicated.
    ///
    /// Args:
    ///  source                  graph to extract subgraph from
    ///  into                    graph to extract into
    ///  max_dist                include all nodes and edges that can be reached in at most this distance
    ///  pos                     extend from this position
    ///  backward                extend in this direction
    ///  preserve_cycles_on_src  if necessary, duplicate starting node to preserve cycles after cutting it
    unordered_map<id_t, id_t> extract_extending_graph(const HandleGraph* source, DeletableHandleGraph* into, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src_node);
                                                      
}
}

#endif
