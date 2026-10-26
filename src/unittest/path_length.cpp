/// \file path_length.cpp
///
/// Unit tests for handlegraph::algorithms::get_path_length.

#include "catch.hpp"
#include "../handle.hpp"

#include "bdsg/hash_graph.hpp"
#include "bdsg/overlays/path_position_overlays.hpp"

#include "handlegraph/algorithms/path_length.hpp"

namespace vg {
namespace unittest {

using handlegraph::algorithms::get_path_length;

TEST_CASE("get_path_length adds up visited node lengths", "[algorithms][path_length]") {
    bdsg::HashGraph graph;
    handle_t n1 = graph.create_handle("AAAA");
    handle_t n2 = graph.create_handle("CC");
    handle_t n3 = graph.create_handle("GGGGG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    path_handle_t path = graph.create_path_handle("path");
    graph.append_step(path, n1);
    graph.append_step(path, n2);
    graph.append_step(path, n3);

    SECTION("a graph without positions is scanned step by step") {
        // HashGraph is not a PathPositionHandleGraph, so this takes the scan_path branch.
        REQUIRE(get_path_length(&graph, path) == 11);
    }

    SECTION("a positional graph reports the same length") {
        bdsg::PositionOverlay positional(&graph);
        REQUIRE(get_path_length(&positional, positional.get_path_handle("path")) == 11);
    }

    SECTION("an empty path has no length") {
        path_handle_t empty = graph.create_path_handle("empty");
        REQUIRE(get_path_length(&graph, empty) == 0);
    }

    SECTION("a step's orientation does not change how much it contributes") {
        handle_t n4 = graph.create_handle("TTT");
        path_handle_t flipped = graph.create_path_handle("flipped");
        graph.append_step(flipped, graph.flip(n4));
        REQUIRE(get_path_length(&graph, flipped) == 3);
    }
}

}
}
