/// \file chain_items.cpp
///  
/// unit tests for chaining logic

#include <iostream>
#include "../algorithms/chain_items.hpp"
#include "../integrated_snarl_finder.hpp"
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

// Define a test chaining space for within-node-only matchings.

struct TestItem {
    pos_t pos;
    size_t source;
    
    inline bool operator==(const TestItem& other) const {
        return (this->pos == other.pos && this->source == other.source);
    }
    
    inline bool operator!=(const TestItem& other) const {
        return !(*this == other);
    }
};
    
struct TestSource {
    size_t start;
    size_t length;
    
    inline bool operator==(const TestSource& other) const {
        return (this->start == other.start && this->length == other.length);
    }
    
    inline bool operator!=(const TestSource& other) const {
        return !(*this == other);
    }
};

}

namespace algorithms {

using unittest::TestItem;
using unittest::TestSource;

template<>
struct ChainingSpace<TestItem, TestSource> : public SourceChainingSpace<TestItem, TestSource> {
    using Item = TestItem;
    using Source = TestSource;
    
    ChainingSpace(const VectorView<Source>& sources,
                  const Aligner& scoring,
                  const SnarlDistanceIndex* distance_index,
                  const HandleGraph* graph) :
        SourceChainingSpace<Item, Source>(sources, scoring, distance_index, graph) {
        
        // Nothing to do!
    }
    
    virtual size_t source_read_start(const Source& source) const {
        return source.start;
    }
    
    virtual size_t source_read_end(const Source& source) const {
        return source.start + source.length;
    }
    
    virtual size_t source_read_length(const Source& source) const {
        return source.length;
    }
    
    pos_t graph_start(const Item& item) const {
        return item.pos;
    }
    
    pos_t graph_end(const Item& item) const {
        pos_t end = item.pos;
        get_offset(end) += this->graph_length(item);
        return end;
    }
    
    size_t graph_path_size(const Item& item) const {
        return 1;
    }
    
    handle_t graph_path_at(const Item& item, size_t index) const {
        return this->graph->get_handle(id(item.pos), is_rev(item.pos));
    }
    
    size_t graph_path_offset(const Item& item) const {
        return offset(item.pos);
    }
    
    
};

}

namespace unittest {

TEST_CASE("A forward-strand fallow region with one missing seed can be filled in", "[chain_items][reseed_fallow_region]") {
    
    HashGraph graph = make_long_graph(10, 10);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    
    vector<TestSource> sources {
        {0, 1},
        {50, 1},
        {90, 5}
    };
    
    vector<TestItem> items {
        {make_pos_t(1, false, 0), 0},
        {make_pos_t(10, false, 0), 2}
    };
    
    auto for_each_pos_for_source_in_subgraph = [&](const TestSource& source, const std::vector<nid_t>& subgraph_ids, const std::function<void(const pos_t&)>& iteratee) -> void {
        if (source.start == 50 && std::find(subgraph_ids.begin(), subgraph_ids.end(), 6) != subgraph_ids.end()) {
            // We've been asked about our no-item source and we have the right node so make an item for it.
            iteratee(make_pos_t(6, false, 0));
        }
    };
    
    Aligner scoring;
    algorithms::ChainingSpace<TestItem, TestSource> space(sources, scoring, &distance_index, &graph);
    
    std::unique_ptr<VectorViewInverse> source_sort_inverse;

    vector<TestItem> reseeded = algorithms::reseed_fallow_region<TestItem, TestSource>(
        items[0],
        items[1],
        space,
        source_sort_inverse,
        for_each_pos_for_source_in_subgraph
    );
                                                                 
    REQUIRE(reseeded.size() == 1);
    REQUIRE(reseeded[0].source == 1);
    REQUIRE(reseeded[0].pos == make_pos_t(6, false, 0));
}

/// Take a table of sources, hit positions, and whether the hit should
/// initially be visible.
///
/// Returns a vector of deduplicated sources, a vector of initially visible
/// hits, a vector of hits that could be recovered, and a function that finds
/// the recoverable hits when asked.
static std::tuple<vector<TestSource>, vector<TestItem>, vector<TestItem>, std::function<void(const TestSource&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>> make_reseed_problem(const vector<tuple<TestSource, pos_t, bool>>& all_items) {
    // Build out sources, initial items, and recoverable items
    vector<TestSource> sources;
    vector<TestItem> items;
    vector<TestItem> recoverable;
    for (size_t item_num = 0; item_num < all_items.size(); item_num++) {
        // Make sure we have all the unique sources.
        if (sources.size() == 0 || std::get<0>(all_items[item_num]) != sources.back()) {
            // New source
            sources.push_back(std::get<0>(all_items[item_num]));
        }
        size_t source_num = sources.size() - 1;
        
        if (std::get<2>(all_items[item_num])) {
            // For each item we should see initially, make it
            items.emplace_back(std::get<1>(all_items[item_num]), source_num);
        } else {
            // If we shouldn't see it initially, store it for later
            recoverable.emplace_back(std::get<1>(all_items[item_num]), source_num);
        }
    }
    
    // We need to capture sources and recoverable by value to ensure that we
    // have access to them after the enclosing function returns.
    auto for_each_pos_for_source_in_subgraph = [sources, recoverable](const TestSource& source, const std::vector<nid_t>& subgraph_ids, const std::function<void(const pos_t&)>& iteratee) -> void {
        // Just do a scan for every query
        for (auto& item : recoverable) {
            if (sources[item.source] == source && subgraph_ids.count(id(item.pos))) {
                iteratee(item.pos);
            }
        }
    };
    
    return std::tie(sources, items, recoverable, for_each_pos_for_source_in_subgraph);
}

TEST_CASE("A fallow region with several seeds can be filled in", "[chain_items][reseed_fallow_region]") {
    
    HashGraph graph = make_long_graph(10, 10);
    auto h = get_handles(graph);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    
    // What items are there, and should they be initially found?
    // We assume this is sorted by source.
    vector<tuple<TestSource, pos_t, bool>> all_items {
        {{0, 1}, make_pos_t(1, false, 0), true},
        {{45, 3}, make_pos_t(5, false, 5), false},
        {{49, 1}, make_pos_t(5, true, 0), false},
        {{50, 1}, make_pos_t(6, false, 0), false},
        {{52, 1}, make_pos_t(6, false, 2), false},
        {{55, 1}, make_pos_t(6, false, 5), false},
        {{55, 1}, make_pos_t(6, false, 6), false},
        {{55, 1}, make_pos_t(6, false, 7), false},
        {{90, 5}, make_pos_t(10, false, 0), true}
    };
    
    auto problem = make_reseed_problem(all_items);
    vector<TestSource>& sources = std::get<0>(problem);
    vector<TestItem>& items = std::get<1>(problem);
    vector<TestItem>& recoverable = std::get<2>(problem);
    auto& for_each_pos_for_source_in_subgraph = std::get<3>(problem);
    
    Aligner scoring;
    algorithms::ChainingSpace<TestItem, TestSource> space(sources, scoring, &distance_index, &graph);
    
    std::unique_ptr<VectorViewInverse> source_sort_inverse;

    vector<TestItem> reseeded = algorithms::reseed_fallow_region<TestItem, TestSource>(
        items.front(),
        items.back(),
        space,
        source_sort_inverse,
        for_each_pos_for_source_in_subgraph
    );
                                                                 
    REQUIRE(reseeded.size() == recoverable.size());
    REQUIRE(reseeded == recoverable);
}

}

}

