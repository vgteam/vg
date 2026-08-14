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
/// Assume things are in the correct order; they are given seed numbers in order of appearance
static vector<algorithms::Anchor> make_anchors(const vector<tuple<size_t, handle_t, size_t, size_t, int>>& test_data, const HandleGraph& graph) {
    vector<algorithms::Anchor> to_score;
    for (size_t i = 0; i < test_data.size(); i++) {
        const auto& item = test_data[i];
        pos_t graph_pos = make_pos_t(graph.get_id(get<1>(item)), graph.get_is_reverse(get<1>(item)), get<2>(item));
        to_score.emplace_back(get<0>(item), graph_pos, get<3>(item), 0, 0, get<4>(item), i);
    }
    
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

static algorithms::SparseAnchorChain run_ziptree_iterator_best_chain(const HashGraph& graph, 
                                                                     const vector<algorithms::Anchor>& anchors,
                                                                     size_t read_length) {
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    // Convert these into Seed type
    vector<SnarlDistanceIndexClusterer::Seed> seeds;
    for (const auto& pos : anchors) {
        ZipCode zipcode;
        zipcode.fill_in_zipcode_from_pos(distance_index, pos.graph_start());
        zipcode.fill_in_full_decoder();
        seeds.push_back({pos.graph_start(), 0, zipcode});
    }

    // Next, make a ZipCodeForest for the graph/seeds
    ZipCodeForest zip_forest;
    zip_forest.fill_in_forest(seeds, distance_index);

    // Make iterator for only the first tree
    // Seriously this is for test cases, only one tree at once
    return algorithms::find_best_chain(anchors, distance_index, graph, read_length,
                                       algorithms::zip_tree_transition_iterator(seeds,
                                                                                zip_forest.trees.front(),
                                                                                std::numeric_limits<size_t>::max(),
                                                                                std::numeric_limits<size_t>::max()
                                                                               )
                                       );
}

static vector<algorithms::SubchainGroup> run_ziptree_iterator_multi_chain(const HashGraph& graph, 
                                                                          const vector<algorithms::Anchor>& anchors,
                                                                          size_t read_length,
                                                                          size_t num_chains) {
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    // Convert these into Seed type
    vector<SnarlDistanceIndexClusterer::Seed> seeds;
    for (const auto& pos : anchors) {
        ZipCode zipcode;
        zipcode.fill_in_zipcode_from_pos(distance_index, pos.graph_start());
        zipcode.fill_in_full_decoder();
        seeds.push_back({pos.graph_start(), 0, zipcode});
    }

    // Next, make a ZipCodeForest for the graph/seeds
    ZipCodeForest zip_forest;
    zip_forest.fill_in_forest(seeds, distance_index);

    // Make iterator for only the first tree
    // Seriously this is for test cases, only one tree at once
    return algorithms::find_best_chains(anchors, distance_index, graph, read_length,
                                        algorithms::zip_tree_transition_iterator(seeds,
                                                                                zip_forest.trees.front(),
                                                                                std::numeric_limits<size_t>::max(),
                                                                                std::numeric_limits<size_t>::max()
                                                                                ),
                                        algorithms::ChainScoringScheme(), num_chains
                                        );
}

TEST_CASE("find_best_chain chains two extensions abutting in read and graph correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    // Set up extensions
    auto to_score = make_anchors({{1, h[1], 1, 9, 9},
                                  {10, h[1], 10, 9, 9}}, graph);
    
    // Actually run the chaining and test
    auto result = run_ziptree_iterator_best_chain(graph, to_score, 20);
    REQUIRE(result.chain_score == (9 + 9));
    REQUIRE(result.anchors == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions abutting in read with a gap in graph correctly", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    // Set up extensions
    auto to_score = make_anchors({{1, h[1], 1, 9, 9},
                                  {10, h[1], 11, 9, 9}}, graph);
    
    // Actually run the chaining and test
    auto result = run_ziptree_iterator_best_chain(graph, to_score, 20);
    // TODO: why is this gap free under the current scoring?
    REQUIRE(result.chain_score == (9 + 9));
    REQUIRE(result.anchors == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions abutting in graph with a gap in read correctly", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    // Set up extensions
    auto to_score = make_anchors({{1, h[1], 1, 9, 9},
                                  {11, h[1], 10, 9, 9}}, graph);
    
    // Actually run the chaining and test
    auto result = run_ziptree_iterator_best_chain(graph, to_score, 20);
    // TODO: why is this gap free under the current scoring?
    REQUIRE(result.chain_score == (9 + 9));
    REQUIRE(result.anchors == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain is willing to leave the main diagonal if the items suggest it", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(10, 10);
    auto h = get_handles(graph);
    // Set up extensions.
    // We're going to have to pay for at least 2 gaps so we need to make sure that doing that is worth it.
    auto to_score = make_anchors({{10, h[1], 0, 10, 10}, // First one on main diagonal
                                  {41, h[4], 0, 10, 10}, // Middle one that is further in the read than the graph
                                  {61, h[6], 0, 10, 10}, // Another middle one that is further in the read than the graph
                                  {100, h[10], 0, 10, 10}}, graph); // Last one on main diagonal
    
    // Actually run the chaining and test
    auto result = run_ziptree_iterator_best_chain(graph, to_score, 120);
    // We should take all of the items in order and not be scared off by the indels.
    REQUIRE(result.anchors == std::vector<size_t>{0, 1, 2, 3});
}

TEST_CASE("Simple X case", "[chain_items]") {
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

    // One anchor on each node
    // read start, graph handle and offset, length, and score
    auto to_score = make_anchors({{1, h[1], 0, 5, 5},
                                  {1, h[2], 0, 5, 5},
                                  {11, h[3], 0, 5, 5},
                                  {21, h[4], 0, 5, 5},
                                  {21, h[5], 0, 5, 5},
                                  {31, h[6], 0, 10, 10},
                                  {31, h[7], 0, 10, 10}}, graph);
    
    /// Actually run the chaining and test
    auto result = run_ziptree_iterator_multi_chain(graph, to_score, 40, 2);
    // We should see all possible paths
    REQUIRE(result.front().subchains.size() == 5);
    REQUIRE(result.front().connections.size() == 5);
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

    // One anchor on each node
    // read start, graph handle and offset, length, and score
    auto to_score = make_anchors({{1, h[1], 0, 10, 10},
                                  {1, h[2], 0, 10, 10},
                                  {11, h[3], 0, 10, 10},
                                  {21, h[4], 0, 10, 10},
                                  {21, h[5], 0, 10, 10},
                                  {31, h[6], 0, 10, 10}}, graph);
    
    // Actually run the chaining and test
    auto result = run_ziptree_iterator_multi_chain(graph, to_score, 40, 2);
    // We should see all possible paths
    REQUIRE(result.front().subchains.size() == 5);
    REQUIRE(result.front().connections.size() == 5);
}

}

}

