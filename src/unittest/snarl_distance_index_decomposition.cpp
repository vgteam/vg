/// \file snarl_distance_index_decomposition.cpp
///
/// Unit tests that check SnarlDistanceIndex against the
/// handlegraph::SnarlDecomposition interface.
///
/// Most tests build a small graph for one rule of the interface. The last
/// checks that indexes built without stored distances give the same
/// traversals as one built with them.

#include "catch.hpp"
#include "../snarl_distance_index.hpp"
#include "../integrated_snarl_finder.hpp"
#include "support/snarl_decomposition_fuzzer.hpp"
#include "support/decomposition_lineup.hpp"
#include "support/randomness.hpp"

#include <bdsg/hash_graph.hpp>

#include <algorithm>
#include <vector>

namespace vg::unittest {

using handlegraph::net_handle_t;
using endpoint_t = SnarlDecomposition::endpoint_t;

/// Build a distance index over a graph, from the IntegratedSnarlFinder's
/// decomposition.
static void index_graph(const HandleGraph& graph, SnarlDistanceIndex& index) {
    IntegratedSnarlFinder finder(graph);
    fill_in_distance_index(&index, &graph, &finder);
}

/// Build a distance index over a graph from a decomposition written as
/// snarl finder events.
static void index_events(const HandleGraph& graph, const std::vector<DecompositionEvent>& events,
                         SnarlDistanceIndex& index) {
    ReplaySnarlFinder finder(&graph, events);
    fill_in_distance_index(&index, &graph, &finder);
}

/// List the children of a net graph element or the root, in the order
/// for_each_child() gives them.
static std::vector<net_handle_t> children_of(const SnarlDistanceIndex& index, const net_handle_t& parent) {
    std::vector<net_handle_t> children;
    index.for_each_child(parent, [&children](const net_handle_t& child) {
        children.push_back(child);
    });
    return children;
}

/// Return true if for_each_traversal() gives a traversal of a net graph
/// element that starts at one endpoint and ends at another.
static bool has_traversal(const SnarlDistanceIndex& index, const net_handle_t& net, endpoint_t start,
                          endpoint_t end) {
    bool found = false;
    index.for_each_traversal(net, [&](const net_handle_t& traversal) {
        found |= index.starts_at(traversal) == start && index.ends_at(traversal) == end;
    });
    return found;
}

/// List the traversals that follow_net_edges() gives from a traversal.
static std::vector<net_handle_t> adjacent_traversals(const SnarlDistanceIndex& index, const HandleGraph& graph,
                                                    const net_handle_t& here, bool go_left) {
    std::vector<net_handle_t> adjacent;
    index.follow_net_edges(here, &graph, go_left, [&adjacent](const net_handle_t& next) {
        adjacent.push_back(next);
    });
    return adjacent;
}

/// Get the start-to-end traversal of the chain that holds a node.
static net_handle_t chain_of_node(const SnarlDistanceIndex& index, nid_t id) {
    net_handle_t chain = index.get_parent(index.get_node_net_handle(id));
    return index.starts_at(chain) == SnarlDecomposition::START ? chain : index.flip(chain);
}

TEST_CASE("SnarlDistanceIndex for_each_child() of an END_START chain lists the end node first",
          "[snarl_distance_interface]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    SnarlDistanceIndex index;
    index_graph(graph, index);

    net_handle_t backward = index.flip(chain_of_node(index, graph.get_id(n1)));
    std::vector<net_handle_t> children = children_of(index, backward);
    REQUIRE(children.size() == 3);
    REQUIRE(index.is_node(children.front()));
    REQUIRE(index.get_handle(children.front(), &graph) == graph.flip(n4));
    REQUIRE(index.get_handle(children.back(), &graph) == graph.flip(n1));
}

TEST_CASE("SnarlDistanceIndex for_each_traversal() gives START_TIP for a snarl with a walk from its start to a tip",
          "[snarl_distance_interface]") {
    // The snarl from 1 to 4 holds nodes 2, 3 and 5. An edge joins the end of 3
    // to node 5, and the end of 5 is a tip. The test replays a decomposition
    // written by hand, because a snarl finder may make node 5 an end of the
    // root chain.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A", 1);
    handle_t n2 = graph.create_handle("C", 2);
    handle_t n3 = graph.create_handle("G", 3);
    handle_t n4 = graph.create_handle("T", 4);
    handle_t n5 = graph.create_handle("A", 5);
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    graph.create_edge(n3, n5);
    using enum DecompositionEventType;
    SnarlDistanceIndex index;
    index_events(graph, {
        {BEGIN_CHAIN, 1, false}, {BEGIN_SNARL, 1, false},
        {BEGIN_CHAIN, 2, false}, {END_CHAIN, 2, false},
        {BEGIN_CHAIN, 3, false}, {END_CHAIN, 3, false},
        {BEGIN_CHAIN, 5, false}, {END_CHAIN, 5, false},
        {END_SNARL, 4, false}, {END_CHAIN, 4, false}
    }, index);

    net_handle_t snarl = index.get_parent(index.get_parent(index.get_node_net_handle(graph.get_id(n2))));
    REQUIRE(index.is_snarl(snarl));
    REQUIRE(has_traversal(index, snarl, SnarlDecomposition::START, SnarlDecomposition::TIP));
}

TEST_CASE("SnarlDistanceIndex for_each_traversal() gives START_TIP for a chain with a walk from its start to a tip in one of its snarls",
          "[snarl_distance_interface]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4. In the snarl, an
    // edge joins the end of 1 to node 3, and the end of 3 is a tip. The test
    // replays a decomposition written by hand, because a snarl finder may make
    // node 3 an end of the root chain.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A", 1);
    handle_t n2 = graph.create_handle("C", 2);
    handle_t n3 = graph.create_handle("G", 3);
    handle_t n4 = graph.create_handle("T", 4);
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n1, n4);
    using enum DecompositionEventType;
    SnarlDistanceIndex index;
    index_events(graph, {
        {BEGIN_CHAIN, 1, false}, {BEGIN_SNARL, 1, false},
        {BEGIN_CHAIN, 2, false}, {END_CHAIN, 2, false},
        {BEGIN_CHAIN, 3, false}, {END_CHAIN, 3, false},
        {END_SNARL, 4, false}, {END_CHAIN, 4, false}
    }, index);

    net_handle_t chain = chain_of_node(index, graph.get_id(n1));
    REQUIRE(index.get_handle(index.get_bound(chain, false, true), &graph) == n1);
    REQUIRE(has_traversal(index, chain, SnarlDecomposition::START, SnarlDecomposition::TIP));
}

TEST_CASE("SnarlDistanceIndex for_each_traversal() gives START_TIP for a snarl with a walk from its start to a tip in a child chain",
          "[snarl_distance_interface]") {
    // The snarl from 1 to 7 holds the chain of node 2, the snarl from 2 to 5,
    // and node 5. In the snarl from 2 to 5, an edge joins the end of 3 to node
    // 6, and the end of 6 is a tip. The test replays a decomposition written
    // by hand, because a snarl finder may make node 6 an end of the root
    // chain.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A", 1);
    handle_t n2 = graph.create_handle("C", 2);
    handle_t n3 = graph.create_handle("G", 3);
    handle_t n4 = graph.create_handle("T", 4);
    handle_t n5 = graph.create_handle("A", 5);
    handle_t n6 = graph.create_handle("C", 6);
    handle_t n7 = graph.create_handle("G", 7);
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n5);
    graph.create_edge(n4, n5);
    graph.create_edge(n3, n6);
    graph.create_edge(n5, n7);
    graph.create_edge(n1, n7);
    using enum DecompositionEventType;
    std::vector<DecompositionEvent> events {
        {BEGIN_CHAIN, 1, false}, {BEGIN_SNARL, 1, false},
        {BEGIN_CHAIN, 2, false}, {BEGIN_SNARL, 2, false},
        {BEGIN_CHAIN, 3, false}, {END_CHAIN, 3, false},
        {BEGIN_CHAIN, 4, false}, {END_CHAIN, 4, false},
        {BEGIN_CHAIN, 6, false}, {END_CHAIN, 6, false},
        {END_SNARL, 5, false}, {END_CHAIN, 5, false},
        {END_SNARL, 7, false}, {END_CHAIN, 7, false}
    };
    SnarlDistanceIndex index;
    index_events(graph, events, index);

    net_handle_t inner_chain = chain_of_node(index, graph.get_id(n2));
    REQUIRE(index.is_chain(inner_chain));
    REQUIRE(children_of(index, inner_chain).size() == 3);
    net_handle_t outer_snarl = index.get_parent(inner_chain);
    REQUIRE(index.is_snarl(outer_snarl));
    REQUIRE(has_traversal(index, outer_snarl, SnarlDecomposition::START, SnarlDecomposition::TIP));
}

TEST_CASE("SnarlDistanceIndex get_bound() of a looping chain gives a START_END or END_START node traversal",
          "[snarl_distance_interface]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4, and an edge
    // from the end of 4 to the start of 1 makes the chain a looping chain.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n1);
    SnarlDistanceIndex index;
    index_graph(graph, index);

    net_handle_t chain = chain_of_node(index, graph.get_id(n1));
    // The interface gives a node only start-to-end and end-to-start
    // traversals.
    for (bool get_end : {false, true}) {
        for (bool face_in : {false, true}) {
            net_handle_t bound = index.get_bound(chain, get_end, face_in);
            INFO("get_end " << get_end << ", face_in " << face_in);
            REQUIRE(index.is_node(bound));
            REQUIRE(index.starts_at(bound) != index.ends_at(bound));
        }
    }
}

TEST_CASE("SnarlDistanceIndex for_each_traversal() gives START_START for a chain with a self-loop in a childless snarl",
          "[snarl_distance_interface]") {
    // The chain is node 1, the snarl from 1 to 4, node 4, a snarl with no
    // children from 4 to 5, and node 5. An edge from the end of 4 back into
    // the end of 4 lets a walk from the chain's start cross the snarl from 1
    // to 4, turn around, and leave by the chain's start.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("A");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n4, graph.flip(n4));
    SnarlDistanceIndex index;
    index_graph(graph, index);

    net_handle_t chain = chain_of_node(index, graph.get_id(n1));
    REQUIRE(has_traversal(index, chain, SnarlDecomposition::START, SnarlDecomposition::START));
}

TEST_CASE("SnarlDistanceIndex follow_net_edges() going right undoes going left from a chain node",
          "[snarl_distance_interface]") {
    // The chain is nodes 1, 2 and 3 in a row. An edge from the start of 2
    // back into the start of 2 makes the index store 2 in a different record
    // from 1.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(graph.flip(n2), n2);
    SnarlDistanceIndex index;
    index_graph(graph, index);

    net_handle_t chain = chain_of_node(index, graph.get_id(n1));
    for (const net_handle_t& child : children_of(index, chain)) {
        for (const net_handle_t& here : {child, index.flip(child)}) {
            for (const net_handle_t& left : adjacent_traversals(index, graph, here, true)) {
                INFO("going left from " << index.net_handle_as_string(here) << " gives "
                     << index.net_handle_as_string(left));
                std::vector<net_handle_t> right = adjacent_traversals(index, graph, left, false);
                REQUIRE(std::find(right.begin(), right.end(), here) != right.end());
            }
        }
    }
}


TEST_CASE("SnarlDistanceIndex get_parent_traversal() of a trivial chain's bounds gives the trivial chain START_END",
          "[snarl_distance_interface]") {
    // The snarl from 1 to 4 holds nodes 2 and 3, each a chain by itself.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    SnarlDistanceIndex index;
    index_graph(graph, index);

    net_handle_t trivial_chain = index.get_parent(index.get_node_net_handle(graph.get_id(n2)));
    REQUIRE(index.is_chain(trivial_chain));
    net_handle_t forward = index.get_parent_traversal(index.get_bound(trivial_chain, false, true),
                                                      index.get_bound(trivial_chain, true, false));
    REQUIRE(index.starts_at(forward) == SnarlDecomposition::START);
    REQUIRE(index.ends_at(forward) == SnarlDecomposition::END);
    REQUIRE(index.canonical(forward) == index.canonical(trivial_chain));
}


TEST_CASE("SnarlDistanceIndex get_parent() of a sentinel gives its snarl START_END",
          "[snarl_distance_interface]") {
    // The snarl from 1 to 4 holds nodes 2 and 3.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    SnarlDistanceIndex index;
    index_graph(graph, index);

    net_handle_t snarl = index.get_parent(index.get_parent(index.get_node_net_handle(graph.get_id(n2))));
    REQUIRE(index.is_snarl(snarl));
    for (bool get_end : {false, true}) {
        for (bool face_in : {false, true}) {
            INFO("get_end " << get_end << ", face_in " << face_in);
            net_handle_t parent = index.get_parent(index.get_bound(snarl, get_end, face_in));
            REQUIRE(index.starts_at(parent) == SnarlDecomposition::START);
            REQUIRE(index.ends_at(parent) == SnarlDecomposition::END);
        }
    }
}


TEST_CASE("SnarlDistanceIndex follow_net_edges() from a chain node gives each traversal of the next snarl from its start",
          "[snarl_distance_interface]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4. In the snarl,
    // an edge from the end of 2 into the end of 3 lets a walk from 1 turn
    // around, so the snarl has START_END and START_START traversals.
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("A");
    handle_t n2 = graph.create_handle("C");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    graph.create_edge(n2, graph.flip(n3));
    SnarlDistanceIndex index;
    index_graph(graph, index);

    net_handle_t chain = chain_of_node(index, graph.get_id(n1));
    REQUIRE(index.get_handle(index.get_bound(chain, false, true), &graph) == n1);
    std::vector<net_handle_t> adjacent = adjacent_traversals(index, graph, index.get_bound(chain, false, true), false);
    REQUIRE(adjacent.size() == 2);
    for (const net_handle_t& next : adjacent) {
        REQUIRE(index.is_snarl(next));
        REQUIRE(index.starts_at(next) == SnarlDecomposition::START);
    }
    REQUIRE(index.ends_at(adjacent[0]) != index.ends_at(adjacent[1]));
}


/// List the start and end endpoints of each realizable traversal of a net
/// graph element, in the order for_each_traversal() gives them.
static std::vector<std::pair<SnarlDecomposition::endpoint_t, SnarlDecomposition::endpoint_t>>
traversal_kinds(const SnarlDistanceIndex& index, const net_handle_t& net) {
    std::vector<std::pair<SnarlDecomposition::endpoint_t, SnarlDecomposition::endpoint_t>> kinds;
    index.for_each_traversal(net, [&](const net_handle_t& traversal) {
        kinds.emplace_back(index.starts_at(traversal), index.ends_at(traversal));
    });
    return kinds;
}

/**
 * Require that two distance indexes of one decomposition agree on a net graph
 * element or the root, and on everything under it.
 *
 * expected_item and actual_item are the same element, or the root, in each
 * index.
 *
 * @post Each snarl, chain and root has the same number of children in both
 * indexes, in the same order.
 * @post Each node has the same graph handle in both indexes.
 * @post Each net graph element has the same realizable traversals in both
 * indexes.
 */
static void require_same_traversals(const SnarlDistanceIndex& expected, const net_handle_t& expected_item,
                                    const SnarlDistanceIndex& actual, const net_handle_t& actual_item,
                                    const HandleGraph& graph) {
    if (!expected.is_root(expected_item)) {
        INFO("Element " << expected.net_handle_as_string(expected_item));
        REQUIRE(traversal_kinds(actual, actual_item) == traversal_kinds(expected, expected_item));
    }
    if (expected.is_node(expected_item)) {
        REQUIRE(actual.get_handle(actual_item, &graph) == expected.get_handle(expected_item, &graph));
        return;
    }
    std::vector<net_handle_t> expected_children = children_of(expected, expected_item);
    std::vector<net_handle_t> actual_children = children_of(actual, actual_item);
    REQUIRE(actual_children.size() == expected_children.size());
    for (size_t i = 0; i < expected_children.size(); i++) {
        require_same_traversals(expected, expected_children[i], actual, actual_children[i], graph);
    }
}

/// Require that distance indexes built from one decomposition without stored
/// distances, with a size limit that leaves snarls oversized, and with
/// distances only for top-level chains, give the same traversals as one built
/// with distances.
static void require_traversals_without_distances(const HandleGraph& graph,
                                                 const std::vector<DecompositionEvent>& events) {
    SnarlDistanceIndex with_distances;
    ReplaySnarlFinder finder(&graph, events);
    fill_in_distance_index(&with_distances, &graph, &finder);
    for (auto [size_limit, only_top_level] : {std::pair<size_t, bool>(0, false), std::pair<size_t, bool>(1, false),
                                              std::pair<size_t, bool>(50000, true)}) {
        INFO("Size limit " << size_limit << (only_top_level ? ", top-level chain distances only" : ""));
        SnarlDistanceIndex index;
        ReplaySnarlFinder replay(&graph, events);
        fill_in_distance_index(&index, &graph, &replay, size_limit, only_top_level);
        require_same_traversals(with_distances, with_distances.get_root(), index, index.get_root(), graph);
    }
}

TEST_CASE("SnarlDistanceIndex for_each_traversal() gives the same traversals with or without stored distances",
          "[snarl_distance_interface]") {
    SECTION("The same traversals for snarls with tips and turnarounds") {
        // The chain is node 1, the snarl from 1 to 4, node 4, a snarl with no
        // children from 4 to 5, and node 5. In the snarl from 1 to 4, node 6
        // hangs off the end of 3 and the end of 6 is a tip. An edge from the
        // end of 2 into the end of 3 lets a walk turn around. An edge from the
        // end of 4 back into the end of 4 is a self-loop.
        bdsg::HashGraph graph;
        handle_t n1 = graph.create_handle("A", 1);
        handle_t n2 = graph.create_handle("C", 2);
        handle_t n3 = graph.create_handle("G", 3);
        handle_t n4 = graph.create_handle("T", 4);
        handle_t n5 = graph.create_handle("A", 5);
        handle_t n6 = graph.create_handle("C", 6);
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n3);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n4);
        graph.create_edge(n3, n6);
        graph.create_edge(n2, graph.flip(n3));
        graph.create_edge(n4, n5);
        graph.create_edge(n4, graph.flip(n4));
        using enum DecompositionEventType;
        require_traversals_without_distances(graph, {
            {BEGIN_CHAIN, 1, false}, {BEGIN_SNARL, 1, false},
            {BEGIN_CHAIN, 2, false}, {END_CHAIN, 2, false},
            {BEGIN_CHAIN, 3, false}, {END_CHAIN, 3, false},
            {BEGIN_CHAIN, 6, false}, {END_CHAIN, 6, false},
            {END_SNARL, 4, false}, {BEGIN_SNARL, 4, false}, {END_SNARL, 5, false},
            {END_CHAIN, 5, false}
        });
    }

    SECTION("The same traversals for a snarl with a tip in a child chain") {
        // The snarl from 1 to 7 holds the chain of node 2, the snarl from 2 to
        // 5, and node 5. In the snarl from 2 to 5, node 6 hangs off the end of
        // 3, and the end of 6 is a tip.
        bdsg::HashGraph graph;
        handle_t n1 = graph.create_handle("A", 1);
        handle_t n2 = graph.create_handle("C", 2);
        handle_t n3 = graph.create_handle("G", 3);
        handle_t n4 = graph.create_handle("T", 4);
        handle_t n5 = graph.create_handle("A", 5);
        handle_t n6 = graph.create_handle("C", 6);
        handle_t n7 = graph.create_handle("G", 7);
        graph.create_edge(n1, n2);
        graph.create_edge(n2, n3);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n5);
        graph.create_edge(n4, n5);
        graph.create_edge(n3, n6);
        graph.create_edge(n5, n7);
        graph.create_edge(n1, n7);
        using enum DecompositionEventType;
        require_traversals_without_distances(graph, {
            {BEGIN_CHAIN, 1, false}, {BEGIN_SNARL, 1, false},
            {BEGIN_CHAIN, 2, false}, {BEGIN_SNARL, 2, false},
            {BEGIN_CHAIN, 3, false}, {END_CHAIN, 3, false},
            {BEGIN_CHAIN, 4, false}, {END_CHAIN, 4, false},
            {BEGIN_CHAIN, 6, false}, {END_CHAIN, 6, false},
            {END_SNARL, 5, false}, {END_CHAIN, 5, false},
            {END_SNARL, 7, false}, {END_CHAIN, 7, false}
        });
    }

    SECTION("The same traversals for a snarl with a walk between tips in two child chains") {
        // The snarl from 1 to 10 holds two chains. The first is node 2, the
        // snarl from 2 to 4, and node 4. In it, node 3 enters the start of 4,
        // and the start of 3 is a tip. The second is node 5, the snarl from 5
        // to 7, and node 7. In it, the end of 5 enters node 8, and the end of
        // 8 is a tip. An edge from the end of 4 to the start of 5 makes a walk
        // from the tip of 3 to the tip of 8.
        bdsg::HashGraph graph;
        handle_t n1 = graph.create_handle("A", 1);
        handle_t n2 = graph.create_handle("C", 2);
        handle_t n3 = graph.create_handle("G", 3);
        handle_t n4 = graph.create_handle("T", 4);
        handle_t n5 = graph.create_handle("A", 5);
        handle_t n7 = graph.create_handle("C", 7);
        handle_t n8 = graph.create_handle("G", 8);
        handle_t n10 = graph.create_handle("T", 10);
        graph.create_edge(n1, n2);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n4);
        graph.create_edge(n4, n10);
        graph.create_edge(n1, n5);
        graph.create_edge(n5, n7);
        graph.create_edge(n5, n8);
        graph.create_edge(n7, n10);
        graph.create_edge(n4, n5);
        using enum DecompositionEventType;
        require_traversals_without_distances(graph, {
            {BEGIN_CHAIN, 1, false}, {BEGIN_SNARL, 1, false},
            {BEGIN_CHAIN, 2, false}, {BEGIN_SNARL, 2, false},
            {BEGIN_CHAIN, 3, false}, {END_CHAIN, 3, false},
            {END_SNARL, 4, false}, {END_CHAIN, 4, false},
            {BEGIN_CHAIN, 5, false}, {BEGIN_SNARL, 5, false},
            {BEGIN_CHAIN, 8, false}, {END_CHAIN, 8, false},
            {END_SNARL, 7, false}, {END_CHAIN, 7, false},
            {END_SNARL, 10, false}, {END_CHAIN, 10, false}
        });
    }

    SECTION("The same traversals for a looping chain") {
        // The chain is node 1, the snarl from 1 to 4, and node 4, and an edge
        // from the end of 4 to the start of 1 closes the loop.
        bdsg::HashGraph graph;
        handle_t n1 = graph.create_handle("A");
        handle_t n2 = graph.create_handle("C");
        handle_t n3 = graph.create_handle("G");
        handle_t n4 = graph.create_handle("T");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n3);
        graph.create_edge(n2, n4);
        graph.create_edge(n3, n4);
        graph.create_edge(n4, n1);
        IntegratedSnarlFinder finder(graph);
        require_traversals_without_distances(graph, capture_events(finder, graph));
    }

    SECTION("The same traversals for random graphs") {
        std::default_random_engine generator(test_seed_source());
        for (size_t repeat = 0; repeat < 100; repeat++) {
            bdsg::HashGraph graph = random_flipped_graph(generator);
            require_traversals_without_distances(graph, fuzzed_decomposition(graph, 0.5, generator));
        }
    }
}

}
