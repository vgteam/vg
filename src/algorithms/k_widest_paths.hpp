#ifndef VG_ALGORITHMS_K_WIDEST_PATHS_HPP_INCLUDED
#define VG_ALGORITHMS_K_WIDEST_PATHS_HPP_INCLUDED

/**
 * \file k_widest_paths.hpp
 *
 * Yen's algorithm to find the K widest 
 */

#include <vector>

#include "../position.hpp"
#include "../handle.hpp"

namespace vg {
namespace algorithms {

/// This Dijkstra is the same underlying algorithm as the one in dijkstra.hpp
/// but the interface is different enough that I opted to make it a seprate
/// thing rather than add loads of optional arguments.   The key differences
/// are these generalizations:
///  -- looks for the "widest" path (maximum minimum weight) instead of shortest
///  -- counts node and edge weights (via callbakcs)
///  -- returns the path as well as the score
///  -- option for ignoring certain nodes and edges in search (required by Yen's algorithm) 
pair<double, vector<handle_t>> widest_dijkstra(const HandleGraph* g, handle_t source, handle_t sink,
                                               function<double(const handle_t&)> node_weight_callback,
                                               function<double(const edge_t&)> edge_weight_callback,
                                               function<bool(const handle_t&)> is_node_ignored_callback,
                                               function<bool(const edge_t&)> is_edge_ignored_callbback);

/// Find the k widest paths
vector<pair<double, vector<handle_t>>> yens_k_widest_paths(const HandleGraph* g, handle_t source, handle_t sink,
                                                           size_t K,
                                                           function<double(const handle_t&)> node_weight_callback,
                                                           function<double(const edge_t&)> edge_weight_callback);

}
}

#endif
