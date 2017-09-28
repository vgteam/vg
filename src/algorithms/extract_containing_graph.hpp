#ifndef VG_EXTRACT_CONTAINING_GRAPH_HPP
#define VG_EXTRACT_CONTAINING_GRAPH_HPP

/**
 * \file extract_containing_graphs.hpp
 *
 * Definitions for the extract_containing_graph algorithm.
 */

#include <vector>

#include "../position.hpp"
#include "../cached_position.hpp"
#include "../xg.hpp"
#include "../vg.hpp"
#include "../vg.pb.h"
#include "../handle.hpp"
#include "../hash_map.hpp"

namespace vg {
namespace algorithms {

    /// Fills graph g with the subgraph of the VG graph vg that contains all of the positions in
    /// the positions vector and all other nodes and edges that can be reached within a maximum
    /// distance from any of these positions. Node IDs in the subgraph are retained from the full graph.
    ///
    /// Args:
    ///  vg                         graph to extract subgraph from
    ///  g                          graph to extract into
    ///  positions                  search outward from these positions
    ///  max_dist                   include all nodes and edges that can be reached in at most this distance
    void extract_containing_graph(VG& vg, Graph& g, const vector<pos_t>& positions, size_t max_dist);
    
    /// Same semantics as previous, but accesses graph through an XG instead of a VG. Optionally uses
    /// an LRUCache to speed up Node queries.
    void extract_containing_graph(xg::XG& xg_index, Graph& g, const vector<pos_t>& positions, size_t max_dist,
                                  LRUCache<id_t, Node>* node_cache = nullptr,
                                  LRUCache<id_t, vector<Edge>>* edge_cache = nullptr);
    
    /// Same semantics as previous except that there is a separate maximum distance for different
    /// positions in the graph. Each distance is associated with the position with the same index. Throws
    /// an error if the position and distance vectors are not the same length.
    void extract_containing_graph(VG& vg, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_max_dist);
    
    /// Same semantics as previous, but accesses graph through an XG instead of a VG. Optionally uses
    /// an LRUCache to speed up Node queries.
    void extract_containing_graph(xg::XG& xg_index, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_max_dist,
                                  LRUCache<id_t, Node>* node_cache = nullptr,
                                  LRUCache<id_t, vector<Edge>>* edge_cache = nullptr);
    
    /// Same semantics as previous except that there is a separate maximum distance for different
    /// positions in the graph and for each search direction. Each distance is associated with the
    /// position with the same index. The forward distance is in the same orientation as the position,
    /// and the backward distance is in the reverse orientation of the position. Throws an error if
    /// the position and distance vectors are not the same length.
    void extract_containing_graph(VG& vg, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_forward_max_dist,
                                  const vector<size_t>& position_backward_max_dist);
    
    /// Same semantics as previous, but accesses graph through an XG instead of a VG. Optionally uses
    /// an LRUCache to speed up Node queries.
    void extract_containing_graph(xg::XG& xg_index, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_forward_max_dist,
                                  const vector<size_t>& position_backward_max_dist,
                                  LRUCache<id_t, Node>* node_cache = nullptr,
                                  LRUCache<id_t, vector<Edge>>* edge_cache = nullptr);
    
    /// Extract the containing graph from a handle graph.               
    void extract_containing_graph(const HandleGraph* source, Graph& g, const vector<pos_t>& positions,
                                  const vector<size_t>& position_forward_max_dist,
                                  const vector<size_t>& position_backward_max_dist);

}
}

#endif
