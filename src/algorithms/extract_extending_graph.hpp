#ifndef VG_EXTRACT_EXTENDING_GRAPH_HPP
#define VG_EXTRACT_EXTENDING_GRAPH_HPP

/**
 * \file extract_estending_graphs.hpp
 *
 * Definitions for the extract_extending_graph algorithm.
 */

#include <unordered_map>

#include "../position.hpp"
#include "../cached_position.hpp"
#include "../xg.hpp"
#include "../vg.hpp"
#include "../vg.pb.h"
#include "../hash_map.hpp"
#include "../handle.hpp"

namespace vg {
namespace algorithms {
    
    /// Fills graph g with the subgraph of the VG graph vg that extends in one direction from a given
    /// position, up to a maximum distance. The node containing the position will be "cut" so that only
    /// the portion that is forward in the search direction remains.
    ///
    /// Args:
    ///  vg                      graph to extract subgraph from
    ///  g                       graph to extract into
    ///  max_dist                include all nodes and edges that can be reached in at most this distance
    ///  pos                     extend from this position
    ///  backward                extend in this direction
    ///  preserve_cycles_on_src  if necessary, duplicate starting node to preserve cycles after cutting it
    unordered_map<id_t, id_t> extract_extending_graph(VG& vg, Graph& g, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src);
    
    /// Same semantics as previous, but accesses graph through an XG instead of a VG. Optionally uses
    /// an LRUCache to speed up Node queries.
    unordered_map<id_t, id_t> extract_extending_graph(xg::XG& xg_index, Graph& g, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src,
                                                      LRUCache<id_t, Node>* node_cache = nullptr,
                                                      LRUCache<id_t, vector<Edge>>* edge_cache = nullptr);
    
    /// Same semantics as previous, but accesses graph through a handle graph interface instead of a VG.
    unordered_map<id_t, id_t> extract_extending_graph(const HandleGraph* source, Graph& g, int64_t max_dist, pos_t pos,
                                                      bool backward, bool preserve_cycles_on_src);
    
}
}

#endif
