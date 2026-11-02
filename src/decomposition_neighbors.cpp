/** \file
 * Stepping from a traversal in a snarl decomposition to each neighboring net
 * graph element.
 */

#include "decomposition_neighbors.hpp"

#include <algorithm>
#include <vector>

namespace vg {

using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

bool for_each_neighbor_impl(const SnarlDecomposition& decomposition, const net_handle_t& here,
                              const HandleGraph* graph, bool go_left,
                              const std::function<bool(const net_handle_t&)>& iteratee) {
    // The canonical traversal of each neighboring net graph element, in the
    // order that follow_net_edges() first gives a traversal of it.
    std::vector<net_handle_t> neighbors;
    // The traversal to give for each element, at the element's index in
    // neighbors.
    std::vector<net_handle_t> chosen;
    decomposition.follow_net_edges(here, graph, go_left, [&](const net_handle_t& traversal) {
        net_handle_t neighbor = decomposition.canonical(traversal);
        // A traversal crosses its element when it goes from the start to the
        // end or from the end to the start.
        bool crosses = decomposition.starts_at(traversal) != decomposition.ends_at(traversal) &&
            !decomposition.starts_at_tip(traversal) && !decomposition.ends_at_tip(traversal);
        for (size_t i = 0; i < neighbors.size(); i++) {
            if (neighbors[i] == neighbor) {
                if (crosses) {
                    chosen[i] = traversal;
                }
                return;
            }
        }
        neighbors.push_back(neighbor);
        chosen.push_back(traversal);
    });
    return std::ranges::all_of(chosen, iteratee);
}

}
