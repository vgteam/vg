/** \file
 * Implements weakly_connected_components_dense.
 */

#include "weakly_connected_components_dense.hpp"

namespace vg {
namespace algorithms {

/**
 * Add to `component` every node that a depth-first search from `start` reaches.
 *
 * `traversed` marks the nodes reached by this and earlier searches, indexed by node ID minus
 * `min_id`. The search does not enter nodes already marked.
 */
static void collect_weak_component(const HandleGraph& graph, handle_t start, nid_t min_id,
                                   std::vector<bool>& traversed,
                                   std::unordered_set<nid_t>& component) {
    // The stack only holds locally forward handles
    std::vector<handle_t> stack{graph.forward(start)};
    while (!stack.empty()) {
        handle_t here = stack.back();
        stack.pop_back();
        traversed[graph.get_id(here) - min_id] = true;
        component.insert(graph.get_id(here));
        auto handle_other = [&](const handle_t& other) {
            handle_t other_forward = graph.forward(other);
            if (!traversed[graph.get_id(other_forward) - min_id]) {
                stack.push_back(other_forward);
            }
        };
        // Look at edges in both directions
        graph.follow_edges(here, false, handle_other);
        graph.follow_edges(here, true, handle_other);
    }
}

std::vector<std::unordered_set<nid_t>> weakly_connected_components_dense(const HandleGraph& graph) {
    std::vector<std::unordered_set<nid_t>> to_return;
    const nid_t min_id = graph.min_node_id();
    const nid_t max_id = graph.max_node_id();
    std::vector<bool> traversed(max_id >= min_id ? (size_t)(max_id - min_id + 1) : 0, false);

    graph.for_each_handle([&](const handle_t& handle) {
        if (traversed[graph.get_id(handle) - min_id]) {
            // Already have this node, so don't start a search from it.
            return;
        }
        to_return.emplace_back();
        collect_weak_component(graph, handle, min_id, traversed, to_return.back());
    });
    return to_return;
}

}
}
