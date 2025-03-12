#ifndef VG_ALGORITHMS_EXTRACT_CONNECTING_GRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_EXTRACT_CONNECTING_GRAPH_HPP_INCLUDED

/**
 * \file extract_connecting_graph.cpp
 *
 * Implementation for the extract_connecting_graph algorithm.
 */

#include <unordered_map>
#include <vg/vg.pb.h>

#include "../position.hpp"
#include "../handle.hpp"
#include "../hash_map.hpp"

namespace vg {
namespace algorithms {
    
    /// Fills a DeletableHandleGraph with the subgraph of a HandleGraph that connects two positions. The nodes
    /// that contain the two positions will be 'cut' at the position and will be tips in the returned graph. The
    /// algorithm guarantees that 'into' contains all walks between pos_1 and pos_2 under the maximum length
    /// except walks that include a cycle involving either position. If no walk between the two positions under
    /// the maximum length exists, 'into' will be left empty. An error is thrown if 'into' is not empty when
    /// passed to function.
    ///
    /// If pos_1 and pos_2 face each other on the same node, the intervening
    /// portion of the node is produced in into. If they are on the same node
    /// but do not face each other, portions of the original node will exist as
    /// distinct nodes in into, and the one correspondign to pos_1 will have
    /// the lower node ID.
    ///
    /// Node local forward orientations are not changed.
    ///
    /// Args:
    ///  source                     graph to extract subgraph from
    ///  into                       graph to extract into
    ///  max_len                    guarantee finding walks along which pos_1 and pos_2 are this distance apart
    ///  pos_1                      start position, subgraph walks begin from here in same orientation
    ///  pos_2                      end position, subgraph walks end here in the same orientation
    ///  strict_max_len             only extract nodes and edges if they fall on some walk between pos_1 and pos_2
    ///                             that is under the maximum length (implies only_walks = true)
    ///
    /// Returns: a map from node ids in the extracted graph to the node ids in
    /// the original graph. The map and the graph will have the same number of
    /// entries.
    unordered_map<id_t, id_t> extract_connecting_graph(const HandleGraph* source,
                                                       DeletableHandleGraph* into,
                                                       int64_t max_len,
                                                       pos_t pos_1, pos_t pos_2,
                                                       bool strict_max_len = false);

}
}

#endif
