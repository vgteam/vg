/// \file snarl_distance_index_decomposition.cpp
///
/// Unit tests that check SnarlDistanceIndex against the rules of the
/// handlegraph::SnarlDecomposition interface, each on the smallest graph that
/// shows a rule being broken.

#include "catch.hpp"
#include "../snarl_distance_index.hpp"
#include "../integrated_snarl_finder.hpp"
#include "support/snarl_decomposition_fuzzer.hpp"

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

/// List the children of a net graph element, as for_each_child() gives them.
static std::vector<net_handle_t> children_of(const SnarlDistanceIndex& index, const net_handle_t& parent) {
    std::vector<net_handle_t> children;
    index.for_each_child(parent, [&children](const net_handle_t& child) {
        children.push_back(child);
    });
    return children;
}

/// Return true if for_each_traversal() lists a traversal from start to end of
/// the element.
static bool has_traversal(const SnarlDistanceIndex& index, const net_handle_t& item, endpoint_t start,
                          endpoint_t end) {
    bool found = false;
    index.for_each_traversal(item, [&](const net_handle_t& traversal) {
        found |= index.starts_at(traversal) == start && index.ends_at(traversal) == end;
    });
    return found;
}

/// List the traversals that follow_net_edges() gives from a traversal.
static std::vector<net_handle_t> neighbours(const SnarlDistanceIndex& index, const HandleGraph& graph,
                                           const net_handle_t& here, bool go_left) {
    std::vector<net_handle_t> reached;
    index.follow_net_edges(here, &graph, go_left, [&reached](const net_handle_t& next) {
        reached.push_back(next);
    });
    return reached;
}

/// Get the start-to-end traversal of the chain that holds a node.
static net_handle_t chain_of_node(const SnarlDistanceIndex& index, nid_t id) {
    net_handle_t chain = index.get_parent(index.get_node_net_handle(id));
    return index.starts_at(chain) == SnarlDecomposition::START ? chain : index.flip(chain);
}

TEST_CASE("A distance index lists an end-to-start chain's children end-to-start",
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

TEST_CASE("A distance index gives a snarl a start-to-tip traversal when its start reaches a tip",
          "[snarl_distance_interface]") {
    // The snarl from 1 to 4 holds nodes 2, 3 and 5, and node 5 hangs off 3
    // with nothing after it. The test replays a decomposition written by
    // hand, because a snarl finder would end the root chain at node 5.
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

TEST_CASE("A distance index gives a chain a start-to-tip traversal when its start reaches a tip",
          "[snarl_distance_interface][!mayfail]") {
    // The chain is node 1, the snarl from 1 to 4, and node 4. In the snarl,
    // node 3 hangs off 1 with nothing after it. The test replays a
    // decomposition written by hand, because a snarl finder would end the
    // root chain at node 3.
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

TEST_CASE("A distance index gives a snarl a start-to-tip traversal through a child chain with a tip",
          "[snarl_distance_interface][!mayfail]") {
    // The snarl from 1 to 7 holds the chain of 2, the snarl from 2 to 5, and
    // 5. In that inner snarl, node 6 hangs off 3 with nothing after it. The
    // test replays a decomposition written by hand, because a snarl finder
    // would end the root chain at node 6.
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

TEST_CASE("A distance index gives a cyclic chain an end bound that is a node traversal",
          "[snarl_distance_interface]") {
    // The chain goes from 1 through the snarl from 1 to 4 and back to 1.
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
    // The interface allows a node only start-to-end and end-to-start
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

TEST_CASE("A distance index lets a chain turn around at a self-loop between two of its nodes",
          "[snarl_distance_interface][!mayfail]") {
    // The chain is node 1, the snarl from 1 to 4, node 4, and node 5. An
    // edge from the end of 4 back into the end of 4 lets a walk from the
    // chain's start cross the snarl, turn around, and leave by the start.
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

TEST_CASE("Going left from a chain node gives the traversals that give it going right",
          "[snarl_distance_interface]") {
    // The chain is nodes 1, 2 and 3 in a row. A self-loop on the start side
    // of 2 makes the index store 2 in a different record from 1.
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
            for (const net_handle_t& left : neighbours(index, graph, here, true)) {
                INFO("going left from " << index.net_handle_as_string(here) << " reaches "
                     << index.net_handle_as_string(left));
                std::vector<net_handle_t> right = neighbours(index, graph, left, false);
                REQUIRE(std::find(right.begin(), right.end(), here) != right.end());
            }
        }
    }
}

}
