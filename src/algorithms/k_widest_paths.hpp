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
///  -- when avg_flow is set to true, use the average instead of minimum flow.
///     the guarantee of getting the true maximum is lost, due to the lack of optimal substructure.
///     (a short crappy path can be a better prefix than a long better-supported path)
pair<double, vector<handle_t>> widest_dijkstra(const HandleGraph* g, handle_t source, handle_t sink,
                                               function<double(const handle_t&)> node_weight_callback,
                                               function<double(const edge_t&)> edge_weight_callback,
                                               function<bool(const handle_t&)> is_node_ignored_callback,
                                               function<bool(const edge_t&)> is_edge_ignored_callbback,
                                               bool avg_flow = false);


/// Find the k widest paths
/// Cutoff: stop searching once best path has score < cutoff
/// Avg_flow: search using maximum average flow instead of width (min flow).  Unlike min flow, the solution
///           is not guarateed to be optimal.
vector<pair<double, vector<handle_t>>> yens_k_widest_paths(const HandleGraph* g, handle_t source, handle_t sink,
                                                           size_t K,
                                                           function<double(const handle_t&)> node_weight_callback,
                                                           function<double(const edge_t&)> edge_weight_callback,
                                                           double cutoff = numeric_limits<double>::lowest(),
                                                           bool avg_flow = false);

}
}

#endif
