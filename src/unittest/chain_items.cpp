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

/// Turn inline test data of start and end positions and scores into GaplessExtensions
static vector<GaplessExtension> fake_extensions(const vector<tuple<size_t, size_t, int>>& test_data) {
    vector<GaplessExtension> to_score;
    for (auto& item : test_data) {
        to_score.emplace_back();
        to_score.back().read_interval = std::make_pair(get<0>(item), get<1>(item));
        to_score.back().score = get<2>(item);
    }
    
    // Sort by read interval as is required
    std::sort(to_score.begin(), to_score.end(), [](const GaplessExtension& a, const GaplessExtension& b) {
        return (a.read_interval.first < b.read_interval.first) ||
            (a.read_interval.first == b.read_interval.first && a.read_interval.second < b.read_interval.second);
    });
    
    return to_score;
}

/// Turn inline test data of start and end positions in read, path and start offset in graph, and score into GaplessExtensions
static vector<GaplessExtension> fake_extensions(const vector<tuple<size_t, size_t, vector<handle_t>, size_t, int>>& test_data) {
    vector<GaplessExtension> to_score;
    for (auto& item : test_data) {
        to_score.emplace_back();
        to_score.back().read_interval = std::make_pair(get<0>(item), get<1>(item));
        to_score.back().path = get<2>(item);
        to_score.back().offset = get<3>(item);
        to_score.back().score = get<4>(item);
    }
    
    // Sort by read interval as is required
    std::sort(to_score.begin(), to_score.end(), [](const GaplessExtension& a, const GaplessExtension& b) {
        return (a.read_interval.first < b.read_interval.first) ||
            (a.read_interval.first == b.read_interval.first && a.read_interval.second < b.read_interval.second);
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

TEST_CASE("ChainingSpace::get_graph_overlap works when extensions intersect but don't overlap for chaining purposes", "[chain_items][get_graph_overlap]") {
    // Set up graph fixture
    HashGraph graph = make_disconnected_graph(5, 4);
    auto h = get_handles(graph);
    
    graph.create_edge(h[1], h[2]);
    
    graph.create_edge(h[5], h[2]);
    
    graph.create_edge(h[2], h[3]);
    graph.create_edge(h[3], h[4]);
    
    // Set up the aligner fixture
    Aligner scoring;
    // Set up the chaining space fixture
    algorithms::ChainingSpace<GaplessExtension> space(scoring, nullptr, &graph);
    
    
    // Set up extensions
    auto extensions = fake_extensions({{1, 10, {h[1], h[2], h[3]}, 1, 9},
                                       {1, 10, {h[5], h[2], h[3]}, 1, 9}});
                                       
    // Check overlap
    REQUIRE(space.get_graph_overlap(extensions[0], extensions[1]) == 0);
    REQUIRE(space.get_graph_overlap(extensions[1], extensions[0]) == 0);
}

TEST_CASE("ChainingSpace::get_graph_overlap works when everything is on one giant handle", "[chain_items][get_graph_overlap]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1, 18);
    auto h = get_handles(graph);
    
    // Set up the aligner fixture
    Aligner scoring;
    // Set up the chaining space fixture
    algorithms::ChainingSpace<GaplessExtension> space(scoring, nullptr, &graph);
    
    // Set up extensions
    auto extensions = fake_extensions({{0, 5, {h[1]}, 0, 5},
                                       {2, 7, {h[1]}, 2, 5},
                                       {4, 9, {h[1]}, 4, 5},
                                       {6, 11, {h[1]}, 6, 5}});
    
    // Check overlap
    REQUIRE(space.get_graph_overlap(extensions[0], extensions[1]) == 3);
    REQUIRE(space.get_graph_overlap(extensions[1], extensions[0]) == 7);
    
    REQUIRE(space.get_graph_overlap(extensions[0], extensions[2]) == 1);
    REQUIRE(space.get_graph_overlap(extensions[2], extensions[0]) == 9);
    
    REQUIRE(space.get_graph_overlap(extensions[0], extensions[3]) == 0);
    REQUIRE(space.get_graph_overlap(extensions[3], extensions[0]) == 11);
}

TEST_CASE("score_best_chain scores no gapless extensions as 0", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<GaplessExtension> to_score;
    REQUIRE(algorithms::score_best_chain(to_score, space) == 0);
}

TEST_CASE("score_best_chain scores a 1-base gapless extension as 1", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 2, 1}});
    REQUIRE(algorithms::score_best_chain(to_score, space) == 1);
}

TEST_CASE("score_best_chain scores two adjacent 1-base gapless extensions as 2", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 2, 1},
                                     {2, 3, 1}});
    REQUIRE(algorithms::score_best_chain(to_score, space) == 2);
}

TEST_CASE("score_best_chain scores one 9-base extension as 9", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9}});
    REQUIRE(algorithms::score_best_chain(to_score, space) == 9);
}

TEST_CASE("score_best_chain scores two abutting 9-base extensions correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {10, 19, 9}});
    REQUIRE(algorithms::score_best_chain(to_score, space) == (9 + 9));
}

TEST_CASE("score_best_chain scores two 9-base extensions separated by a 1-base gap correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {11, 20, 9}});
    // We assume the gap is the same in the graph and the read, and that it could all be matches.
    REQUIRE(algorithms::score_best_chain(to_score, space) == (9 + 9 + 1));
}

TEST_CASE("score_best_chain scores two 9-base extensions separated by a 2-base gap correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {12, 21, 9}});
    // We assume the gap is the same in the graph and the read, and that it could all be matches.
    REQUIRE(algorithms::score_best_chain(to_score, space) == (9 + 9 + 2));
}

TEST_CASE("score_best_chain scores two 9-base extensions with a 1-base overlap correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {9, 18, 9}});
    // The answer should be that we don't allow the overlap.
    REQUIRE(algorithms::score_best_chain(to_score, space) == 9);
}

TEST_CASE("score_best_chain scores two 9-base extensions with a 2-base overlap correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {8, 17, 9}});
    // The answer should be that we don't allow the overlap.
    REQUIRE(algorithms::score_best_chain(to_score, space) == 9);
}

TEST_CASE("find_best_chain chains no gapless extensions as score 0", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<GaplessExtension> to_score;
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == 0);
    REQUIRE(result.second.size() == 0);
}

TEST_CASE("find_best_chain chains a 1-base gapless extension as score 1", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 2, 1}});
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == 1);
    REQUIRE(result.second == std::vector<size_t>{0});
}

TEST_CASE("find_best_chain chains two adjacent 1-base gapless extensions as score 2", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 2, 1},
                                     {2, 3, 1}});
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == 2);
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains one 9-base extension as 9", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9}});
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == 9);
    REQUIRE(result.second == std::vector<size_t>{0});
}

TEST_CASE("find_best_chain chains two abutting 9-base extensions correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {10, 19, 9}});
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two 9-base extensions separated by a 1-base gap correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {11, 20, 9}});
    auto result = algorithms::find_best_chain(to_score, space);
    // We assume the gap is the same in the graph and the read, and that it could all be matches.
    REQUIRE(result.first == (9 + 9 + 1));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two 9-base extensions separated by a 2-base gap correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {12, 21, 9}});
    auto result = algorithms::find_best_chain(to_score, space);
    // We assume the gap is the same in the graph and the read, and that it could all be matches.
    REQUIRE(result.first == (9 + 9 + 2));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two 9-base extensions with a 1-base overlap correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {9, 18, 9}});
    auto result = algorithms::find_best_chain(to_score, space);
    // The answer should be that we don't allow the overlap.
    REQUIRE(result.first == 9);
    REQUIRE(result.second == std::vector<size_t>{0});
}

TEST_CASE("find_best_chain chains two 9-base extensions with a 2-base overlap correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {8, 17, 9}});
    auto result = algorithms::find_best_chain(to_score, space);
    // The answer should be that we don't allow the overlap.
    REQUIRE(result.first == 9);
    REQUIRE(result.second == std::vector<size_t>{0});
}

TEST_CASE("find_best_chain chains two extensions abutting in read and graph correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 10, {h[1]}, 1, 9},
                                     {10, 19, {h[1]}, 10, 9}});
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions abutting in read with a gap in graph correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 10, {h[1]}, 1, 9},
                                     {10, 19, {h[1]}, 11, 9}});
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions abutting in graph with a gap in read correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(1);
    auto h = get_handles(graph);
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 10, {h[1]}, 1, 9},
                                     {11, 20, {h[1]}, 10, 9}});
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions that abut over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(5, 4);
    auto h = get_handles(graph);
    graph.create_edge(h[1], h[2]);
    graph.create_edge(h[2], h[3]);
    graph.create_edge(h[3], h[4]);
    graph.create_edge(h[4], h[5]);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //             ** **** ***
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 10, {h[1], h[2], h[3]}, 1, 9},
                                     {10, 19, {h[3], h[4], h[5]}, 2, 9}});
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions that abut in graph with a gap in read over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(5, 4);
    auto h = get_handles(graph);
    graph.create_edge(h[1], h[2]);
    graph.create_edge(h[2], h[3]);
    graph.create_edge(h[3], h[4]);
    graph.create_edge(h[4], h[5]);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //             ** **** ***
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 10, {h[1], h[2], h[3]}, 1, 9},
                                     {11, 20, {h[3], h[4], h[5]}, 2, 9}});
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions that abut in read with a gap in graph over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(5, 4);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //              * **** ****
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 10, {h[1], h[2], h[3]}, 1, 9},
                                     {10, 19, {h[3], h[4], h[5]}, 3, 9}});
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions that have gaps of different sizes in graph and read over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(5, 4);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //              * **** ****
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 10, {h[1], h[2], h[3]}, 1, 9},
                                     {12, 21, {h[3], h[4], h[5]}, 3, 9}});
    
    // Gap in graph is 1 and gap in read is 2
    REQUIRE(space.get_graph_distance(to_score[0], to_score[1]) == 1);
    REQUIRE(space.get_read_distance(to_score[0], to_score[1]) == 2);
    // So this is going to be 1bp of possible match and 1 bp of indel
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    // So we pay 1 gap open, the difference in gap length.
    REQUIRE(result.first == (9 + 9 + 1 - 6));
    REQUIRE(result.second == std::vector<size_t>{0, 1});
}

TEST_CASE("find_best_chain chains two extensions that overlap the same amount in graph and read over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(5, 4);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //         * **** ****
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 10, {h[1], h[2], h[3]}, 1, 9},
                                     {7, 16, {h[2], h[3], h[4]}, 3, 9}});
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    // The answer should be that we don't allow the overlap.
    REQUIRE(result.first == 9);
    REQUIRE(result.second == std::vector<size_t>{0});
}

TEST_CASE("find_best_chain is willing to leave the main diagonal if the items suggest it", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(10, 10);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions
    auto to_score = fake_extensions({{10, 20, {h[1]}, 0, 10}, // First one on main diagonal
                                     {41, 51, {h[4]}, 0, 10}, // Middle one that is further in the read than the graph
                                     {100, 110, {h[10]}, 0, 10}}); // Last one on main diagonal
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    // We should take all of the items in order and not be scared off by the indels.
    REQUIRE(result.second == std::vector<size_t>{0, 1, 2});
}

TEST_CASE("find_best_chain refuses to let the score go negative", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(10, 1);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions
    auto to_score = fake_extensions({{1, 2, {h[1]}, 0, 1}, // First one on main diagonal
                                     {5, 6, {h[4]}, 0, 1}, // Middle one that is further in the read than the graph
                                     {10, 11, {h[10]}, 0, 1}}); // Last one on main diagonal
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    // We'd have to let the score go negative to take either of the indels, so we just stay with the first best hit.
    REQUIRE(result.second == std::vector<size_t>{0});
}

TEST_CASE("find_best_chain is willing to keep indels split if the items suggest it", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(10, 10);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions
    auto to_score = fake_extensions({{10, 20, {h[1]}, 0, 10}, // First one on main diagonal
                                     {41, 51, {h[4]}, 0, 10}, // Middle one that is further in the read than the graph
                                     {89, 99, {h[9]}, 0, 10}, // Middle one that is further in the graph than the read
                                     {100, 110, {h[10]}, 0, 10}}); // Last one on main diagonal
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    // We should take all of the items in order and not be scared off by the indels.
    REQUIRE(result.second == std::vector<size_t>{0, 1, 2, 3});
}

}

}

