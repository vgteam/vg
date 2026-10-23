/**
 * \file decomposition_comparison.cpp
 * Names for the net graph elements of a snarl tree that every
 * SnarlDecomposition implementation shares.
 */

#include "decomposition_comparison.hpp"

#include <algorithm>
#include <functional>
#include <string>

namespace vg::unittest {

/// Write a graph handle as its node ID, then + if forward or - if reverse.
static std::string handle_text(const HandleGraph& graph, const handle_t& handle) {
    return std::to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+");
}

/**
 * Write the bound pair of a snarl or chain both ways: from its start bound to
 * its end bound, and from its end bound flipped to its start bound flipped.
 */
static std::pair<std::string, std::string> bound_pairs(const SnarlDecomposition& decomposition,
                                                       const HandleGraph& graph, const net_handle_t& element) {
    handle_t enter = decomposition.get_handle(decomposition.get_bound(element, false, true), &graph);
    handle_t leave = decomposition.get_handle(decomposition.get_bound(element, true, false), &graph);
    return {handle_text(graph, enter) + "/" + handle_text(graph, leave),
            handle_text(graph, graph.flip(leave)) + "/" + handle_text(graph, graph.flip(enter))};
}

std::string element_name(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& element) {
    if (decomposition.is_node(element)) {
        return "node " + std::to_string(graph.get_id(decomposition.get_handle(element, &graph)));
    }
    auto [forward, flipped] = bound_pairs(decomposition, graph, element);
    return (decomposition.is_snarl(element) ? "snarl " : "chain ") + std::min(forward, flipped);
}

/**
 * Write where a traversal starts and ends as two letters: S, E or T.
 *
 * The letters are relative to the bound pair in the element_name() of the
 * element, so that implementations that orient the element differently give
 * the same letters. A node's letters are relative to its forward orientation.
 */
static std::string traversal_letters(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                     const net_handle_t& traversal) {
    bool swap;
    if (decomposition.is_node(traversal)) {
        swap = graph.get_is_reverse(decomposition.get_handle(decomposition.get_bound(traversal, false, true), &graph));
    } else {
        auto [forward, flipped] = bound_pairs(decomposition, graph, traversal);
        swap = flipped < forward;
    }
    auto letter = [&](SnarlDecomposition::endpoint_t endpoint) {
        if (endpoint == SnarlDecomposition::TIP) {
            return 'T';
        }
        return (endpoint == SnarlDecomposition::START) != swap ? 'S' : 'E';
    };
    return {letter(decomposition.starts_at(traversal)), letter(decomposition.ends_at(traversal))};
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
            std::vector<std::string>& traversals = shape[child_name + " traversals"];
            decomposition.for_each_traversal(child, [&](const net_handle_t& traversal) {
                traversals.push_back(traversal_letters(decomposition, graph, traversal));
            });
            std::ranges::sort(traversals);
            std::vector<std::string>& tippy = shape[child_name + " tippy children"];
            decomposition.for_each_tippy_child(child, [&](const net_handle_t& traversal) {
                tippy.push_back(element_name(decomposition, graph, traversal) + " " +
                                traversal_letters(decomposition, graph, traversal));
            });
            std::ranges::sort(tippy);
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
