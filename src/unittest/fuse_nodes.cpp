/// \file fuse_nodes.cpp
///
/// Unit tests for the fuse_nodes algorithm.

#include "catch.hpp"
#include "../handle.hpp"
#include "../algorithms/fuse_nodes.hpp"

#include "bdsg/hash_graph.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// Spell out the sequence a path visits, so tests can assert on what a path
/// means rather than on which nodes happen to carry it.
static string spell_path(const PathHandleGraph& graph, const string& name) {
    string sequence;
    path_handle_t path = graph.get_path_handle(name);
    for (handle_t h : graph.scan_path(path)) {
        sequence += graph.get_sequence(h);
    }
    return sequence;
}

TEST_CASE("can_fuse_nodes accepts a clean boundary and rejects the rest", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAA");
    handle_t n2 = graph.create_handle("CC");
    handle_t n3 = graph.create_handle("GGGG");
    handle_t n4 = graph.create_handle("TT");
    graph.create_edge(n1, n2);
    graph.create_edge(n3, n4);

    SECTION("nothing on the inner sides of the pair") {
        REQUIRE(algorithms::can_fuse_nodes(graph, n2, n3));
    }
    SECTION("the left node has something on its right") {
        REQUIRE(!algorithms::can_fuse_nodes(graph, n1, n3));
    }
    SECTION("the right node has something on its left") {
        REQUIRE(!algorithms::can_fuse_nodes(graph, n2, n4));
    }
    SECTION("a node cannot be fused with itself") {
        REQUIRE(!algorithms::can_fuse_nodes(graph, n2, n2));
    }
    SECTION("an edge between the pair would be swallowed, so it is refused") {
        graph.create_edge(n2, n3);
        REQUIRE(!algorithms::can_fuse_nodes(graph, n2, n3));
    }
}

TEST_CASE("fuse_nodes concatenates sequence, re-hangs edges, and keeps paths", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAA");
    handle_t n2 = graph.create_handle("CC");
    handle_t n3 = graph.create_handle("GGGG");
    handle_t n4 = graph.create_handle("TT");
    graph.create_edge(n1, n2);
    graph.create_edge(n3, n4);

    // Three paths over the seam: one ends on it, one starts on it, and one
    // crosses it backwards, an orientation the rewrite has to preserve.
    path_handle_t left_path = graph.create_path_handle("left");
    graph.append_step(left_path, n1);
    graph.append_step(left_path, n2);
    path_handle_t right_path = graph.create_path_handle("right");
    graph.append_step(right_path, n3);
    graph.append_step(right_path, n4);
    path_handle_t back_path = graph.create_path_handle("backwards");
    graph.append_step(back_path, graph.flip(n2));
    graph.append_step(back_path, graph.flip(n1));

    REQUIRE(algorithms::can_fuse_nodes(graph, n2, n3));
    handle_t fused = algorithms::fuse_nodes(&graph, n2, n3);

    REQUIRE(graph.get_sequence(fused) == "CCGGGG");
    REQUIRE(graph.get_node_count() == 3);

    size_t edge_count = 0;
    graph.for_each_edge([&](const edge_t&) { edge_count++; });
    REQUIRE(edge_count == 2);
    REQUIRE(graph.has_edge(n1, fused));
    REQUIRE(graph.has_edge(fused, n4));

    REQUIRE(spell_path(graph, "left") == "AAACCGGGG");
    REQUIRE(spell_path(graph, "right") == "CCGGGGTT");
    REQUIRE(spell_path(graph, "backwards") == "CCCCGGTTT");
    REQUIRE(graph.get_step_count(graph.get_path_handle("left")) == 2);
}

TEST_CASE("fuse_nodes turns a cycle through the pair into a self-loop", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAA");
    handle_t n2 = graph.create_handle("CC");
    handle_t n3 = graph.create_handle("GGGG");
    graph.create_edge(n1, n2);
    // Back edge closing a cycle over the pair. Neither end of it sits on the
    // inner side of the seam, so the pair is still fusable.
    graph.create_edge(n3, n2);

    REQUIRE(algorithms::can_fuse_nodes(graph, n2, n3));
    handle_t fused = algorithms::fuse_nodes(&graph, n2, n3);

    REQUIRE(graph.get_sequence(fused) == "CCGGGG");
    REQUIRE(graph.get_node_count() == 2);
    REQUIRE(graph.has_edge(n1, fused));
    REQUIRE(graph.has_edge(fused, fused));

    size_t edge_count = 0;
    graph.for_each_edge([&](const edge_t&) { edge_count++; });
    REQUIRE(edge_count == 2);
}

}
}
