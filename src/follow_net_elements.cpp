/** \file
 * Walking a snarl decomposition one net graph element at a time.
 */

#include "follow_net_elements.hpp"

#include <vector>

namespace vg {

using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

bool follow_net_elements_impl(const SnarlDecomposition& decomposition, const net_handle_t& here,
                              const HandleGraph* graph, bool go_left,
                              const std::function<bool(const net_handle_t&)>& iteratee) {
    // The canonical traversal of each element, in the order that
    // follow_net_edges() first gives a traversal of the element.
    std::vector<net_handle_t> elements;
    // The traversal to give for each element, at the element's index in
    // elements.
    std::vector<net_handle_t> chosen;
    decomposition.follow_net_edges(here, graph, go_left, [&](const net_handle_t& traversal) {
        net_handle_t element = decomposition.canonical(traversal);
        bool crosses = decomposition.starts_at(traversal) != decomposition.ends_at(traversal) &&
            !decomposition.starts_at_tip(traversal) && !decomposition.ends_at_tip(traversal);
        for (size_t i = 0; i < elements.size(); i++) {
            if (elements[i] == element) {
                if (crosses) {
                    chosen[i] = traversal;
                }
                return;
            }
        }
        elements.push_back(element);
        chosen.push_back(traversal);
    });
    for (const net_handle_t& traversal : chosen) {
        if (!iteratee(traversal)) {
            return false;
        }
    }
    return true;
}

}
