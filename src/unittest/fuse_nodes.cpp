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

TEST_CASE("can_fuse_nodes accepts a left node with nothing on its right and a right node with nothing on its left", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t a = graph.create_handle("TT");
    handle_t left = graph.create_handle("AAA");
    handle_t right = graph.create_handle("CC");
    handle_t b = graph.create_handle("GGGG");
    // Edges on the outer sides don't matter.
    graph.create_edge(a, left);
    graph.create_edge(right, b);
    REQUIRE(algorithms::can_fuse_nodes(graph, left, right));
}

TEST_CASE("can_fuse_nodes rejects a left node with an edge on its right side", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t left = graph.create_handle("AAA");
    handle_t right = graph.create_handle("CC");
    handle_t other = graph.create_handle("GGGG");
    graph.create_edge(left, other);
    REQUIRE(!algorithms::can_fuse_nodes(graph, left, right));
}

TEST_CASE("can_fuse_nodes rejects a right node with an edge on its left side", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t left = graph.create_handle("AAA");
    handle_t right = graph.create_handle("CC");
    handle_t other = graph.create_handle("GGGG");
    graph.create_edge(other, right);
    REQUIRE(!algorithms::can_fuse_nodes(graph, left, right));
}

TEST_CASE("can_fuse_nodes rejects fusing a node with itself", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t node = graph.create_handle("AAA");
    REQUIRE(!algorithms::can_fuse_nodes(graph, node, node));
}

TEST_CASE("can_fuse_nodes rejects two nodes already joined by an edge", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t left = graph.create_handle("AAA");
    handle_t right = graph.create_handle("CC");
    graph.create_edge(left, right);
    REQUIRE(!algorithms::can_fuse_nodes(graph, left, right));
}

TEST_CASE("fuse_nodes concatenates the sequences and keeps outer edges and paths", "[fuse_nodes][algorithms]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAA");
    handle_t n2 = graph.create_handle("CC");
    handle_t n3 = graph.create_handle("GGGG");
    handle_t n4 = graph.create_handle("TT");
    graph.create_edge(n1, n2);
    graph.create_edge(n3, n4);

    // One path ends on the fused pair, one starts on it, and one visits it in reverse;
    // fusing must keep all three, including the reverse one and its orientation.
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
    // This edge closes a cycle, but it doesn't touch the sides being fused,
    // so fusing is still allowed.
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
