/// \file chain_items.cpp
///  
/// unit tests for chaining logic

#include <iostream>
#include "../algorithms/chain_items.hpp"
#include "../integrated_snarl_finder.hpp"
#include "bdsg/hash_graph.hpp"
#include "minimizer_mapper.hpp"
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

static pair<int, vector<size_t>> run_ziptree_iterator(const HashGraph& graph,
                                                      const vector<algorithms::Anchor>& anchors) {
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    // Convert these into Seed type
    vector<SnarlDistanceIndexClusterer::Seed> seeds;
    vector<vg::MinimizerMapper::Minimizer> minimizers;
    for (const auto& pos : anchors) {
        ZipCode zipcode;
        zipcode.fill_in_zipcode_from_pos(distance_index, pos.graph_start());
        zipcode.fill_in_full_decoder();
        seeds.push_back({pos.graph_start(), 0, zipcode});

        minimizers.emplace_back();
        minimizers.back().value.offset = pos.read_start();
        minimizers.back().value.is_reverse = false;
    }

    VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

    // Next, make a ZipCodeForest for the graph/seeds
    ZipCodeForest zip_forest;
    zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());

    // Make iterator for only the first tree
    // Seriously this is for test cases, only one tree at once
    return algorithms::find_best_chain(anchors, distance_index, graph, 6, 1,
                                       algorithms::zip_tree_transition_iterator(seeds,
                                                                                zip_forest.trees.front(),
                                                                                std::numeric_limits<size_t>::max(),
                                                                                std::numeric_limits<size_t>::max()
                                                                               )
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
    auto result = run_ziptree_iterator(graph, to_score);
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions abutting in read with a gap in graph correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    // Set up extensions
    auto to_score = make_anchors({{1, h[1], 1, 9, 9},
                                  {10, h[1], 11, 9, 9}}, graph);
    
    // Actually run the chaining and test
    auto result = run_ziptree_iterator(graph, to_score);
    // TODO: why is this gap free under the current scoring?
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions abutting in graph with a gap in read correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    // Set up extensions
    auto to_score = make_anchors({{1, h[1], 1, 9, 9},
                                  {11, h[1], 10, 9, 9}}, graph);
    
    // Actually run the chaining and test
    auto result = run_ziptree_iterator(graph, to_score);
    // TODO: why is this gap free under the current scoring?
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain is willing to leave the main diagonal if the items suggest it", "[chain_items][find_best_chain]") {
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
    auto result = run_ziptree_iterator(graph, to_score);
    // We should take all of the items in order and not be scared off by the indels.
    REQUIRE(result.second == std::vector<size_t>{0, 1, 2, 3});
}

}

}

