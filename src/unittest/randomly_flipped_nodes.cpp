#include "catch.hpp"
#include "../handle.hpp"
#include "../utility.hpp"
#include <bdsg/hash_graph.hpp>

#include "support/randomly_flipped_nodes.hpp"
#include "support/randomness.hpp"
#include "support/random_graph.hpp"

#include <set>
#include <random>

namespace vg {
namespace unittest {

/// Get the canonicalized set of edge sequence pairs from a graph.
/// Each edge is represented as a pair of sequences (left_seq, right_seq) read
/// in the orientation of the edge. To canonicalize, we compare each pair
/// against its reverse complement (RC(right_seq), RC(left_seq)) and keep the
/// lexicographically smaller one.
///
/// This doesn't fully constrain the graph, but if this doesn't match what it's
/// supposed to, it can tell us that the graph smells off and is wrong.
static set<pair<string, string>> canonical_edge_pairs(const HandleGraph& graph) {
    set<pair<string, string>> result;
    graph.for_each_edge([&](const edge_t& edge) {
        string left_seq = graph.get_sequence(edge.first);
        string right_seq = graph.get_sequence(edge.second);

        // The reverse complement pair: RC(right) on the left, RC(left) on the right
        string rc_right = reverse_complement(right_seq);
        string rc_left = reverse_complement(left_seq);

        pair<string, string> forward_pair = {left_seq, right_seq};
        pair<string, string> rc_pair = {rc_right, rc_left};

        // Use the lexicographically smaller one as canonical
        if (rc_pair < forward_pair) {
            result.insert(rc_pair);
        } else {
            result.insert(forward_pair);
        }
        return true;
    });
    return result;
}

/// Make sure that observed and expected graphs are not obviously not
/// isomorphic.
static void validate_graph(const HandleGraph& observed, const HandleGraph& expected, const set<pair<string, string>>& expected_edges) {
    REQUIRE(observed.get_node_count() == expected.get_node_count());
    REQUIRE(observed.get_edge_count() == expected.get_edge_count());

    auto observed_edges = canonical_edge_pairs(observed);
    REQUIRE(observed_edges == expected_edges);
}

TEST_CASE("randomly_flipped_nodes preserves graph structure on a simple linear graph", "[randomly_flipped_nodes]") {
    bdsg::HashGraph graph;
    std::string stick_sequence = "GGACTGACTCGCATGTCGAGCGACTCGCGCGAGCTATCGTAGTACGCGAGTCATATTATATTATCACG";
    size_t node_length = 3;
    handle_t prev_handle;
    for (size_t i = 0; i < stick_sequence.size(); i += node_length) {
        handle_t h = graph.create_handle(stick_sequence.substr(i, node_length));
        if (i > 0) {
            graph.create_edge(prev_handle, h);
        }
        prev_handle = h;
    }

    auto original_edges = canonical_edge_pairs(graph);

    SECTION("flipping no nodes preserves edges exactly") {
        default_random_engine gen(test_seed_source());
        auto flipped = randomly_flipped_nodes(graph, 0.0, gen);
        validate_graph(flipped, graph, original_edges);
    }

    SECTION("flipping all nodes preserves canonical edge pairs") {
        default_random_engine gen(test_seed_source());
        auto flipped = randomly_flipped_nodes(graph, 1.0, gen);
        validate_graph(flipped, graph, original_edges);
    }

    SECTION("flipping 50% of nodes preserves canonical edge pairs") {
        default_random_engine gen(test_seed_source());
        auto flipped = randomly_flipped_nodes(graph, 0.5, gen);
        validate_graph(flipped, graph, original_edges);
    }
}

TEST_CASE("randomly_flipped_nodes preserves structure on graph with reversing edges", "[randomly_flipped_nodes]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GATT", 1);
    handle_t h2 = graph.create_handle("ACA", 2);
    handle_t h3 = graph.create_handle("CGAT", 3);
    handle_t h4 = graph.create_handle("TCGAA", 4);

    // Forward edges
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    // Reversing edge: 4 fwd -> 3 rev
    graph.create_edge(h4, graph.flip(h3));

    auto original_edges = canonical_edge_pairs(graph);

    default_random_engine gen(test_seed_source());
    for (int i = 0; i < 10; i++) {
        auto flipped = randomly_flipped_nodes(graph, 0.5, gen);
        validate_graph(flipped, graph, original_edges);
    }
}

TEST_CASE("randomly_flipped_nodes preserves structure on graph with self-loops", "[randomly_flipped_nodes]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("ACGT", 1);
    handle_t h2 = graph.create_handle("TTCC", 2);

    graph.create_edge(h1, h2);
    // Self-loop on h1: fwd -> fwd
    graph.create_edge(h1, h1);
    // Inverting self-loop on h2: fwd -> rev
    graph.create_edge(h2, graph.flip(h2));

    auto original_edges = canonical_edge_pairs(graph);

    default_random_engine gen(test_seed_source());
    for (int i = 0; i < 10; i++) {
        auto flipped = randomly_flipped_nodes(graph, 0.5, gen);
        validate_graph(flipped, graph, original_edges);
    }
}

TEST_CASE("randomly_flipped_nodes preserves structure on random graphs", "[randomly_flipped_nodes]") {
    for (int trial = 0; trial < 5; trial++) {
        bdsg::HashGraph graph;
        random_graph(100, 10, 10, &graph);

        auto original_edges = canonical_edge_pairs(graph);

        default_random_engine gen(test_seed_source());
        for (int i = 0; i < 5; i++) {
            auto flipped = randomly_flipped_nodes(graph, 0.5, gen);
            validate_graph(flipped, graph, original_edges);
        }
    }
}

TEST_CASE("randomly_flipped_nodes preserves node IDs", "[randomly_flipped_nodes]") {
    bdsg::HashGraph graph;
    graph.create_handle("AAA", 5);
    graph.create_handle("CCC", 10);
    graph.create_handle("GGG", 15);
    graph.create_edge(graph.get_handle(5), graph.get_handle(10));
    graph.create_edge(graph.get_handle(10), graph.get_handle(15));

    default_random_engine gen(test_seed_source());
    auto flipped = randomly_flipped_nodes(graph, 0.5, gen);

    REQUIRE(flipped.has_node(5));
    REQUIRE(flipped.has_node(10));
    REQUIRE(flipped.has_node(15));
}

TEST_CASE("randomly_flipped_nodes actually flips node sequences", "[randomly_flipped_nodes]") {
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("AAAC", 1);  // RC = GTTT

    default_random_engine gen(test_seed_source());
    // Guarantee a flip
    auto flipped = randomly_flipped_nodes(graph, 1.0, gen);

    // The forward sequence should be the RC of the original
    REQUIRE(flipped.get_sequence(flipped.get_handle(1)) == "GTTT");
}

} // namespace unittest
} // namespace vg
