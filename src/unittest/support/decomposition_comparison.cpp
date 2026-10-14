/**
 * \file decomposition_comparison.cpp
 * Implements TraversalLabeler and require_same_tree().
 */

#include "decomposition_comparison.hpp"

#include <algorithm>
#include <format>
#include <map>
#include <set>
#include <vector>

#include "../catch.hpp"

namespace vg::unittest {

using endpoint_t = SnarlDecomposition::endpoint_t;

/// @brief Describe a graph handle as its node ID and + if forward or - if
/// reverse.
static std::string handle_text(const HandleGraph& graph, const handle_t& handle) {
    return std::format("{}{}", graph.get_id(handle), graph.get_is_reverse(handle) ? "-" : "+");
}

/// @brief Get the letter a node label uses for an endpoint.
static char endpoint_letter(endpoint_t endpoint) {
    switch (endpoint) {
    case SnarlDecomposition::START:
        return 'S';
    case SnarlDecomposition::END:
        return 'E';
    default:
        return 'T';
    }
}

TraversalLabeler::TraversalLabeler(const SnarlDecomposition& decomposition, const HandleGraph& graph) :
    decomposition(decomposition), graph(graph) {
}

std::string TraversalLabeler::side(const net_handle_t& net, endpoint_t endpoint, bool entering) const {
    if (endpoint == SnarlDecomposition::TIP) {
        return "tip";
    }
    bool at_start = (endpoint == SnarlDecomposition::START);
    // The start bound facing in and the end bound facing out both read along
    // the element.
    handle_t along = decomposition.get_handle(decomposition.get_bound(net, !at_start, at_start), &graph);
    return handle_text(graph, entering == at_start ? along : graph.flip(along));
}

std::string TraversalLabeler::bounds_label(const net_handle_t& net) const {
    handle_t start = decomposition.get_handle(decomposition.get_bound(net, false, true), &graph);
    handle_t end = decomposition.get_handle(decomposition.get_bound(net, true, false), &graph);
    return std::min(handle_text(graph, start) + "/" + handle_text(graph, end),
                    handle_text(graph, graph.flip(end)) + "/" + handle_text(graph, graph.flip(start)));
}

std::string TraversalLabeler::element_label(const net_handle_t& net) const {
    if (decomposition.is_node(net)) {
        return std::format("node {}", graph.get_id(decomposition.get_handle(net, &graph)));
    }
    if (decomposition.is_snarl(net) || decomposition.is_chain(net)) {
        return std::format("{} {}", decomposition.is_snarl(net) ? "snarl" : "chain", bounds_label(net));
    }
    return label(net);
}

std::string TraversalLabeler::label(const net_handle_t& net) const {
    if (decomposition.is_root(net)) {
        return "root";
    }
    if (decomposition.is_node(net)) {
        return std::format("node {} {}{}", handle_text(graph, decomposition.get_handle(net, &graph)),
                           endpoint_letter(decomposition.starts_at(net)), endpoint_letter(decomposition.ends_at(net)));
    }
    if (decomposition.is_sentinel(net)) {
        return std::format("sentinel {} of {}", handle_text(graph, decomposition.get_handle(net, &graph)),
                           bounds_label(decomposition.get_parent(net)));
    }
    return std::format("{} {}>{}", decomposition.is_snarl(net) ? "snarl" : "chain",
                       side(net, decomposition.starts_at(net), true),
                       side(net, decomposition.ends_at(net), false));
}

/// @brief Traversals of one decomposition, by label.
using LabelledTraversals = std::map<std::string, net_handle_t, std::less<>>;

/**
 * @brief Collect, by label, the realizable traversals under a net graph
 * element that start and end at bounds, and the bound sentinels of the
 * snarls under it.
 *
 * @param net The element, listed as its parent lists it.
 * @param found The collection to add each traversal to, under its label.
 */
static void collect_traversals(const SnarlDecomposition& decomposition, const TraversalLabeler& labeler,
                               const net_handle_t& net, LabelledTraversals& found) {
    if (!decomposition.is_root(net)) {
        decomposition.for_each_traversal(net, [&decomposition, &labeler, &found](const net_handle_t& traversal) {
            if (!decomposition.starts_at_tip(traversal) && !decomposition.ends_at_tip(traversal)) {
                found[labeler.label(traversal)] = traversal;
            }
        });
    }
    if (decomposition.is_snarl(net)) {
        for (bool get_end : {false, true}) {
            for (bool face_in : {false, true}) {
                net_handle_t bound = decomposition.get_bound(net, get_end, face_in);
                found[labeler.label(bound)] = bound;
            }
        }
    }
    if (!decomposition.is_node(net)) {
        decomposition.for_each_child(net, [&decomposition, &labeler, &found](const net_handle_t& child) {
            collect_traversals(decomposition, labeler, child, found);
        });
    }
}

/// @brief List the labels of a collection's traversals.
static std::vector<std::string> labels_of(const LabelledTraversals& traversals) {
    std::vector<std::string> labels;
    for (const auto& [label, traversal] : traversals) {
        labels.push_back(label);
    }
    return labels;
}

/**
 * @brief Label the children of a traversal.
 *
 * @details For an end-to-start traversal of a chain, works out the children
 * from those of the start-to-end traversal, read in reverse.
 *
 * TODO: A distance index lists a chain's children start-to-end for every
 * traversal of the chain, where the interface wants an end-to-start
 * traversal's children listed end-to-start, so we cannot ask it for them.
 *
 * @return For a start-to-end or end-to-start traversal of a chain, the
 * children's traversal labels in the order the traversal reads them.
 * Otherwise, the children's element labels, sorted.
 */
static std::vector<std::string> child_labels(const SnarlDecomposition& decomposition, const TraversalLabeler& labeler,
                                             const net_handle_t& traversal) {
    std::vector<std::string> labels;
    if (decomposition.is_node(traversal) || decomposition.is_sentinel(traversal)) {
        return labels;
    }
    bool forward = decomposition.starts_at(traversal) == SnarlDecomposition::START &&
        decomposition.ends_at(traversal) == SnarlDecomposition::END;
    if (bool backward = decomposition.starts_at(traversal) == SnarlDecomposition::END &&
                        decomposition.ends_at(traversal) == SnarlDecomposition::START;
        decomposition.is_chain(traversal) && (forward || backward)) {
        std::vector<net_handle_t> children;
        decomposition.for_each_child(forward ? traversal : decomposition.flip(traversal),
                                     [&children](const net_handle_t& child) {
            children.push_back(child);
        });
        if (backward) {
            std::ranges::reverse(children);
            for (net_handle_t& child : children) {
                child = decomposition.flip(child);
            }
        }
        for (const net_handle_t& child : children) {
            labels.push_back(labeler.label(child));
        }
        return labels;
    }
    decomposition.for_each_child(traversal, [&labeler, &labels](const net_handle_t& child) {
        labels.push_back(labeler.element_label(child));
    });
    std::ranges::sort(labels);
    return labels;
}

/**
 * @brief Label the traversals reached by following net edges from a
 * traversal.
 *
 * @details TODO: Going left from a chain node to the neighbouring chain node,
 * a distance index reads the node it reaches the way the traversal reads, or
 * the way it travels, depending on how it packed the nodes into records. So
 * a chain node reached that way is labelled by its element label.
 *
 * @return The labels, sorted.
 */
static std::vector<std::string> neighbour_labels(const SnarlDecomposition& decomposition,
                                                 const TraversalLabeler& labeler, const HandleGraph& graph,
                                                 const net_handle_t& traversal, bool go_left) {
    std::vector<std::string> labels;
    decomposition.follow_net_edges(traversal, &graph, go_left,
                                   [&decomposition, &labeler, &labels, &traversal, go_left](const net_handle_t& next) {
        bool node_to_node = go_left && decomposition.is_node(traversal) && decomposition.is_node(next);
        labels.push_back(node_to_node ? labeler.element_label(next) : labeler.label(next));
    });
    std::ranges::sort(labels);
    return labels;
}

/**
 * @brief Return true if a parent in one decomposition matches a parent in
 * another.
 *
 * @param must_match_way_round True if the parent must be read the same way
 * in both, and false if it may be read either way round.
 */
static bool parents_match(const TraversalLabeler& expected_labeler, const net_handle_t& expected_parent,
                          const SnarlDecomposition& actual, const TraversalLabeler& actual_labeler,
                          const net_handle_t& actual_parent, bool must_match_way_round) {
    std::string wanted = expected_labeler.label(expected_parent);
    return actual_labeler.label(actual_parent) == wanted ||
        (!must_match_way_round && actual_labeler.label(actual.flip(actual_parent)) == wanted);
}

void require_same_tree(const SnarlDecomposition& expected, const SnarlDecomposition& actual,
                       const HandleGraph& graph) {
    TraversalLabeler expected_labeler(expected, graph);
    TraversalLabeler actual_labeler(actual, graph);
    LabelledTraversals expected_traversals;
    LabelledTraversals actual_traversals;
    collect_traversals(expected, expected_labeler, expected.get_root(), expected_traversals);
    collect_traversals(actual, actual_labeler, actual.get_root(), actual_traversals);
    REQUIRE(labels_of(actual_traversals) == labels_of(expected_traversals));

    for (const auto& [label, expected_traversal] : expected_traversals) {
        INFO("Comparing at " << label);
        const net_handle_t& actual_traversal = actual_traversals.at(label);
        REQUIRE(child_labels(actual, actual_labeler, actual_traversal) ==
                child_labels(expected, expected_labeler, expected_traversal));
        if (!expected.is_sentinel(expected_traversal)) {
            net_handle_t expected_parent = expected.get_parent(expected_traversal);
            bool directed = expected.starts_at(expected_traversal) != expected.ends_at(expected_traversal) &&
                !expected.starts_at_tip(expected_traversal) && !expected.ends_at_tip(expected_traversal);
            REQUIRE(parents_match(expected_labeler, expected_parent, actual, actual_labeler,
                                  actual.get_parent(actual_traversal),
                                  directed && expected.is_chain(expected_parent)));
        }
        for (bool go_left : {false, true}) {
            INFO("Going " << (go_left ? "left" : "right"));
            REQUIRE(neighbour_labels(actual, actual_labeler, graph, actual_traversal, go_left) ==
                    neighbour_labels(expected, expected_labeler, graph, expected_traversal, go_left));
        }
    }
}

}
