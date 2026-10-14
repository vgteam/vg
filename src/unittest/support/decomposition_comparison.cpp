/**
 * \file decomposition_comparison.cpp
 * Describes a SnarlDecomposition's snarl tree in names that any
 * implementation shares, so that tests can compare implementations.
 */

#include "decomposition_comparison.hpp"

#include <algorithm>
#include <format>
#include <functional>

namespace vg::unittest {

/// Write a graph handle as its node ID, then + if forward or - if reverse.
static std::string handle_text(const HandleGraph& graph, const handle_t& handle) {
    return std::format("{}{}", graph.get_id(handle), graph.get_is_reverse(handle) ? "-" : "+");
}

std::string element_name(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& element) {
    if (decomposition.is_node(element)) {
        return std::format("node {}", graph.get_id(decomposition.get_handle(element, &graph)));
    }
    handle_t enter = decomposition.get_handle(decomposition.get_bound(element, false, true), &graph);
    handle_t leave = decomposition.get_handle(decomposition.get_bound(element, true, false), &graph);
    std::string bound_pair = std::min(handle_text(graph, enter) + "/" + handle_text(graph, leave),
                                      handle_text(graph, graph.flip(leave)) + "/" +
                                      handle_text(graph, graph.flip(enter)));
    return std::format("{} {}", decomposition.is_snarl(element) ? "snarl" : "chain", bound_pair);
}

/// Return true if a decomposition element is a snarl with no children.
static bool is_childless_snarl(const SnarlDecomposition& decomposition, const net_handle_t& element) {
    return decomposition.is_snarl(element) && decomposition.for_each_child(element, [](const net_handle_t&) {
        return false;
    });
}

/**
 * @brief Add an element, and everything under it, to a tree shape.
 *
 * @param name The element's name in the shape.
 * @param element A start-to-end traversal of the element, or the root.
 */
static void add_to_shape(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const std::string& name, const net_handle_t& element, TreeShape& shape) {
    std::vector<std::string>& children = shape[name];
    // Name each child, and add the children with children of their own.
    std::function<void(const net_handle_t&)> visit = [&](const net_handle_t& child) {
        if (decomposition.is_root(child)) {
            // A grouping inside the root counts as part of the root.
            decomposition.for_each_child(child, visit);
            return;
        }
        if (is_childless_snarl(decomposition, child)) {
            return;
        }
        std::string child_name = element_name(decomposition, graph, child);
        children.push_back(child_name);
        if (!decomposition.is_node(child)) {
            add_to_shape(decomposition, graph, child_name, decomposition.canonical(child), shape);
        }
    };
    if (decomposition.is_chain(element)) {
        // List the children along a start-to-end traversal of the chain.
        net_handle_t forward = decomposition.get_parent(decomposition.get_bound(element, false, true));
        decomposition.for_each_child(forward, visit);
        std::vector<std::string> reversed(children.rbegin(), children.rend());
        bool cyclic = decomposition.get_handle(decomposition.get_bound(element, false, true), &graph) ==
            decomposition.get_handle(decomposition.get_bound(element, true, false), &graph);
        if (cyclic && !reversed.empty()) {
            // Start the reverse order at the bound node too.
            std::rotate(reversed.begin(), reversed.end() - 1, reversed.end());
        }
        children = std::min(children, reversed);
    } else {
        decomposition.for_each_child(element, visit);
        std::ranges::sort(children);
    }
}

TreeShape tree_shape(const SnarlDecomposition& decomposition, const HandleGraph& graph) {
    TreeShape shape;
    add_to_shape(decomposition, graph, "root", decomposition.get_root(), shape);
    return shape;
}

}
