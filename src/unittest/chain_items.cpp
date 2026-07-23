/// \file chain_items.cpp
///  
/// unit tests for chaining logic

#include <iostream>
#include "../algorithms/chain_items.hpp"
#include "../integrated_snarl_finder.hpp"
#include "bdsg/hash_graph.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

/// Turn inline test data of read start, graph handle and offset, length, and score into Anchor objects.
static vector<algorithms::Anchor> make_anchors(const vector<tuple<size_t, handle_t, size_t, size_t, int>>& test_data, const HandleGraph& graph) {
    vector<algorithms::Anchor> to_score;
    for (auto& item : test_data) {
        pos_t graph_pos = make_pos_t(graph.get_id(get<1>(item)), graph.get_is_reverse(get<1>(item)), get<2>(item));
        to_score.emplace_back(get<0>(item), graph_pos, get<3>(item), 0, 0, get<4>(item));
    }
    
    // Sort by read interval as is required
    std::sort(to_score.begin(), to_score.end(), [](const algorithms::Anchor& a, const algorithms::Anchor& b) {
        return (a.read_start() < b.read_start()) ||
            (a.read_start() == b.read_start() && a.length() < b.length());
    });
    
    return to_score;
}

/// Make a disconnected graph of fixed-length nodes
static HashGraph make_disconnected_graph(size_t nodes, size_t length = 32) {
    HashGraph graph;
    
    // What node sequence should we use for everything?
    string seq(length, 'A');
    
    for (size_t i = 0; i < nodes; i++) {
        // Make all the nodes
        graph.create_handle(seq, (nid_t) (i + 1));
    }
    
    return graph;
}

/// Make a long graph of fixed-length nodes
static HashGraph make_long_graph(size_t nodes, size_t length = 32) {
    HashGraph graph = make_disconnected_graph(nodes, length);
    
    for (size_t i = 1; i < nodes; i++) {
        // Link them up
        graph.create_edge(graph.get_handle((nid_t) i, false), graph.get_handle((nid_t) (i + 1), false));
    }
    
    return graph;
}

/// Get a vector of all handles in the graph by node ID.
static vector<handle_t> get_handles(const HandleGraph& graph) {
    vector<handle_t> handles;
    // We just assume the graph has small dense node IDs and pack them in there.
    handles.resize(graph.max_node_id() + 1);
    graph.for_each_handle([&](const handle_t& h) {
        handles[graph.get_id(h)] = h;
    });
    return handles;
}

TEST_CASE("find_best_chain chains two extensions abutting in read and graph correctly", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    
    // Set up extensions
    auto to_score = make_anchors({{1, h[1], 1, 9, 9},
                                  {10, h[1], 10, 9, 9}}, graph);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, distance_index, graph);
    REQUIRE(result.chain_score == (9 + 9));
    REQUIRE(result.anchors == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions abutting in read with a gap in graph correctly", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    
    // Set up extensions
    auto to_score = make_anchors({{1, h[1], 1, 9, 9},
                                  {10, h[1], 11, 9, 9}}, graph);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, distance_index, graph);
    // TODO: why is this gap free under the current scoring?
    REQUIRE(result.chain_score == (9 + 9));
    REQUIRE(result.anchors == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions abutting in graph with a gap in read correctly", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    // Set up extensions
    auto to_score = make_anchors({{1, h[1], 1, 9, 9},
                                  {11, h[1], 10, 9, 9}}, graph);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, distance_index, graph);
    // TODO: why is this gap free under the current scoring?
    REQUIRE(result.chain_score == (9 + 9));
    REQUIRE(result.anchors == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain is willing to leave the main diagonal if the items suggest it", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(10, 10);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    // Set up extensions.
    // We're going to have to pay for at least 2 gaps so we need to make sure that doing that is worth it.
    auto to_score = make_anchors({{10, h[1], 0, 10, 10}, // First one on main diagonal
                                  {41, h[4], 0, 10, 10}, // Middle one that is further in the read than the graph
                                  {61, h[6], 0, 10, 10}, // Another middle one that is further in the read than the graph
                                  {100, h[10], 0, 10, 10}}, graph); // Last one on main diagonal
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, distance_index, graph);
    // We should take all of the items in order and not be scared off by the indels.
    REQUIRE(result.anchors == std::vector<size_t>{0, 1, 2, 3});
}

TEST_CASE("find_best_chain gets edge case possibility for X", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_disconnected_graph(5, 10);
    // 1-3-5 diagonal
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(3, false));
    graph.create_edge(graph.get_handle(3, false), graph.get_handle(5, false));
    // 2-3-4 diagonal
    graph.create_edge(graph.get_handle(2, false), graph.get_handle(3, false));
    graph.create_edge(graph.get_handle(3, false), graph.get_handle(4, false));
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    // One anchor on each node
    // read start, graph handle and offset, length, and score
    auto to_score = make_anchors({{1, h[1], 0, 5, 5},
                                  {1, h[2], 0, 5, 5},
                                  {11, h[3], 0, 5, 5},
                                  {21, h[4], 0, 5, 5},
                                  {21, h[5], 0, 5, 5}}, graph);
    
    /// Actually run the chaining and test
    auto result = algorithms::find_best_chains(to_score, distance_index, graph, algorithms::ChainScoringScheme(), 2);
    // We should see all possible paths
    REQUIRE(result.subchains.size() == 5);
    REQUIRE(result.connections.size() == 4);
}

TEST_CASE("X with different length chains", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_disconnected_graph(7, 10);
    // 1-3-5-6 diagonal
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(3, false));
    graph.create_edge(graph.get_handle(3, false), graph.get_handle(5, false));
    graph.create_edge(graph.get_handle(5, false), graph.get_handle(6, false));
    // 2-3-4-7 diagonal
    graph.create_edge(graph.get_handle(2, false), graph.get_handle(3, false));
    graph.create_edge(graph.get_handle(3, false), graph.get_handle(4, false));
    graph.create_edge(graph.get_handle(4, false), graph.get_handle(7, false));
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    // One anchor on each node
    // read start, graph handle and offset, length, and score
    auto to_score = make_anchors({{1, h[1], 0, 10, 10},
                                  {1, h[2], 0, 10, 10},
                                  {11, h[3], 0, 10, 10},
                                  {21, h[4], 0, 10, 10},
                                  {21, h[5], 0, 10, 10},
                                  {31, h[6], 0, 10, 10}}, graph);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chains(to_score, distance_index, graph, algorithms::ChainScoringScheme(), 2);
    // We should see all possible paths
    REQUIRE(result.subchains.size() == 5);
    REQUIRE(result.connections.size() == 4);
}

}

}

