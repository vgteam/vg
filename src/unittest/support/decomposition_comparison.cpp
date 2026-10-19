/**
 * \file decomposition_comparison.cpp
 * Names for the snarls, chains and nodes of a snarl tree that every
 * SnarlDecomposition implementation shares.
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

/// Return true if a traversal is of a snarl with no children.
static bool is_childless_snarl(const SnarlDecomposition& decomposition, const net_handle_t& element) {
    return decomposition.is_snarl(element) && decomposition.for_each_child(element, [](const net_handle_t&) {
        return false;
    });
}

/**
 * Add a snarl, chain or root, and everything under it, to a tree shape.
 *
 * @param name The name of the element in the shape.
 * @param element A traversal of a snarl or chain, or the root.
 */
static void add_to_shape(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const std::string& name, const net_handle_t& element, TreeShape& shape) {
    std::vector<std::string>& children = shape[name];
    // Name each child, and add each snarl and chain child to the shape.
    std::function<void(const net_handle_t&)> visit = [&](const net_handle_t& child) {
        if (decomposition.is_root(child)) {
            // Replace the child with its own children.
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
        // List the children in the order a start-to-end traversal of the chain
        // visits them, and also in the reverse order.
        net_handle_t forward = decomposition.get_parent(decomposition.get_bound(element, false, true));
        decomposition.for_each_child(forward, visit);
        std::vector<std::string> reversed(children.rbegin(), children.rend());
        if (bool looping = decomposition.get_handle(decomposition.get_bound(element, false, true), &graph) ==
                decomposition.get_handle(decomposition.get_bound(element, true, false), &graph);
            looping && !reversed.empty()) {
            // Start the reverse order at the boundary node too. A trivial chain
            // also has one handle for both bounds, and its one child stays put.
            std::ranges::rotate(reversed, reversed.end() - 1);
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
