/// \file chain_items.cpp
///  
/// unit tests for chaining logic

#include <iostream>
#include "../algorithms/chain_items.hpp"
#include "../integrated_snarl_finder.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

/// Generate an Alignment that all the given extensions could be against the sequence of
static Alignment fake_alignment(const vector<GaplessExtension>& extensions) {
    
    // Work out how many bases we need
    size_t max_length = 0;
    for (auto& extension : extensions) {
        max_length = std::max(max_length, extension.read_interval.second);
    }
    // Generate them
    stringstream s;
    for (size_t i = 0; i < max_length; i++) {
        s << "A";
    }
    
    // Wrap in an Alignment
    Alignment aln;
    aln.set_sequence(s.str());
    return aln;
}

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

TEST_CASE("ChainingSpace::get_graph_overlap works when extensions intersect but don't overlap for chaining purposes", "[chain_items][get_graph_overlap]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAA");
    handle_t h2 = graph.create_handle("AAAA");
    handle_t h3 = graph.create_handle("AAAA");
    handle_t h4 = graph.create_handle("AAAA");
    handle_t h5 = graph.create_handle("AAAA");
    
    graph.create_edge(h1, h2);
    
    graph.create_edge(h5, h2);
    
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    
    // Set up the aligner fixture
    Aligner scoring;
    // Set up the chaining space fixture
    algorithms::ChainingSpace<GaplessExtension> space(scoring, nullptr, &graph);
    
    
    // Set up extensions
    auto extensions = fake_extensions({{1, 10, {h1, h2, h3}, 1, 9},
                                       {1, 10, {h5, h2, h3}, 1, 9}});
                                       
    // Check overlap
    REQUIRE(space.get_graph_overlap(extensions[0], extensions[1]) == 0);
    REQUIRE(space.get_graph_overlap(extensions[1], extensions[0]) == 0);
}

TEST_CASE("ChainingSpace::get_graph_overlap works when everything is on one giant handle", "[chain_items][get_graph_overlap]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAAAAAAAAAAAAAAAA");
    
    // Set up the aligner fixture
    Aligner scoring;
    // Set up the chaining space fixture
    algorithms::ChainingSpace<GaplessExtension> space(scoring, nullptr, &graph);
    
    // Set up extensions
    auto extensions = fake_extensions({{0, 5, {h1}, 0, 5},
                                       {2, 7, {h1}, 2, 5},
                                       {4, 9, {h1}, 4, 5},
                                       {6, 11, {h1}, 6, 5}});
    
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
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == 0);
}

TEST_CASE("score_best_chain scores a 1-base gapless extension as 1", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 2, 1}});
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == 1);
}

TEST_CASE("score_best_chain scores two adjacent 1-base gapless extensions as 2", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 2, 1},
                                     {2, 3, 1}});
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == 2);
}

TEST_CASE("score_best_chain scores one 9-base extension as 9", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9}});
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == 9);
}

TEST_CASE("score_best_chain scores two abutting 9-base extensions correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {10, 19, 9}});
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == (9 + 9));
}

TEST_CASE("score_best_chain scores two 9-base extensions separated by a 1-base gap correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {11, 20, 9}});
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == (9 + 9 - 6));
}

TEST_CASE("score_best_chain scores two 9-base extensions separated by a 2-base gap correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {12, 21, 9}});
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == (9 + 9 - 6 - 1));
}

TEST_CASE("score_best_chain scores two 9-base extensions with a 1-base overlap correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {9, 18, 9}});
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == (9 + 9 - 6));
}

TEST_CASE("score_best_chain scores two 9-base extensions with a 2-base overlap correctly", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {8, 17, 9}});
    auto aln = fake_alignment(to_score);
    REQUIRE(algorithms::score_best_chain(to_score, space) == (9 + 9 - 6 - 1));
}

TEST_CASE("score_best_chain scores correctly when other overlapping extensions exist", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<tuple<size_t, size_t, int>> plan{{0, 1, 1}};
    for (size_t i = 0; i < 35; i++) {
        plan.emplace_back(i + 1, i + 1 + 30, 30);
    }
    auto to_score = fake_extensions(plan);
    auto aln = fake_alignment(to_score);
    // We should get a 1-base extension and two 30-base extensions
    REQUIRE(algorithms::score_best_chain(to_score, space) == (1 + 30 + 30));
}

TEST_CASE("score_best_chain scores a backtrack correctly when other overlapping extensions exist", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<tuple<size_t, size_t, int>> plan{{0, 1, 1}, {28, 28 + 45, 45}};
    for (size_t i = 0; i < 35; i++) {
        plan.emplace_back(i + 1, i + 1 + 30, 30);
    }
    auto to_score = fake_extensions(plan);
    auto aln = fake_alignment(to_score);
    // We should get a 1-base extension, a 30-base extension, a backtrack of 3, and a 45-base extension
    REQUIRE(algorithms::score_best_chain(to_score, space) == (1 + 30 + 45 - 6 - 2));
}

TEST_CASE("score_best_chain scores a backtrack correctly when even more overlapping extensions exist", "[chain_items][score_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<tuple<size_t, size_t, int>> plan{{0, 1, 1}, {28, 28 + 45, 45}};
    for (size_t i = 0; i < 35; i++) {
        plan.emplace_back(i + 1, i + 1 + 30, 30);
    }
    for (size_t i = 3; i < 29; i++) {
        plan.emplace_back(i, i + 15, 15);
    }
    auto to_score = fake_extensions(plan);
    auto aln = fake_alignment(to_score);
    // We should get a 1-base extension, a 30-base extension, a backtrack of 3, and a 45-base extension
    REQUIRE(algorithms::score_best_chain(to_score, space) == (1 + 30 + 45 - 6 - 2));
}


TEST_CASE("find_best_chain chains no gapless extensions as 0", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<GaplessExtension> to_score;
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == 0);
    REQUIRE(result.second.size() == 0);
}

TEST_CASE("find_best_chain chains a 1-base gapless extension as 1", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 2, 1}});
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == 1);
    REQUIRE(result.second.at(0) == 0);
}

TEST_CASE("find_best_chain chains two adjacent 1-base gapless extensions as 2", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 2, 1},
                                     {2, 3, 1}});
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == 2);
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains one 9-base extension as 9", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9}});
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == 9);
    REQUIRE(result.second.at(0) == 0);
}

TEST_CASE("find_best_chain chains two abutting 9-base extensions correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {10, 19, 9}});
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two 9-base extensions separated by a 1-base gap correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {11, 20, 9}});
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two 9-base extensions separated by a 2-base gap correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {12, 21, 9}});
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6 - 1));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two 9-base extensions with a 1-base overlap correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {9, 18, 9}});
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two 9-base extensions with a 2-base overlap correctly", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    auto to_score = fake_extensions({{1, 10, 9},
                                     {8, 17, 9}});
    auto aln = fake_alignment(to_score);
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6 - 1));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains correctly when other overlapping extensions exist", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<tuple<size_t, size_t, int>> plan{{0, 1, 1}};
    for (size_t i = 0; i < 35; i++) {
        plan.emplace_back(i + 1, i + 1 + 30, 30);
    }
    auto to_score = fake_extensions(plan);
    auto aln = fake_alignment(to_score);
    // We should get a 1-base extension and two 30-base extensions
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (1 + 30 + 30));
    REQUIRE(to_score.at(result.second.at(0)).read_interval.first == 0);
    REQUIRE(to_score.at(result.second.at(1)).read_interval.first == 1);
    REQUIRE(to_score.at(result.second.at(2)).read_interval.first == 31);
}

TEST_CASE("find_best_chain chains a backtrack correctly when other overlapping extensions exist", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<tuple<size_t, size_t, int>> plan{{0, 1, 1}, {28, 28 + 45, 45}};
    for (size_t i = 0; i < 35; i++) {
        plan.emplace_back(i + 1, i + 1 + 30, 30);
    }
    auto to_score = fake_extensions(plan);
    auto aln = fake_alignment(to_score);
    // We should get a 1-base extension, a 30-base extension, a backtrack of 3, and a 45-base extension
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (1 + 30 + 45 - 6 - 2));
    REQUIRE(to_score.at(result.second.at(0)).read_interval.first == 0);
    REQUIRE(to_score.at(result.second.at(1)).read_interval.first == 1);
    REQUIRE(to_score.at(result.second.at(2)).read_interval.first == 28);
}

TEST_CASE("find_best_chain chains a backtrack correctly when even more overlapping extensions exist", "[chain_items][find_best_chain]") {
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring);
    vector<tuple<size_t, size_t, int>> plan{{0, 1, 1}, {28, 28 + 45, 45}};
    for (size_t i = 0; i < 35; i++) {
        plan.emplace_back(i + 1, i + 1 + 30, 30);
    }
    for (size_t i = 3; i < 29; i++) {
        plan.emplace_back(i, i + 15, 15);
    }
    auto to_score = fake_extensions(plan);
    auto aln = fake_alignment(to_score);
    // We should get a 1-base extension, a 30-base extension, a backtrack of 3,
    // and a 45-base extension, or the equivalently-scored 1-base extension,
    // 30-base extension, backtrack of 18, 15 base extension, 45 base
    // extension, which is allowed even though it has one extension containing
    // another.
    //
    // TODO: In practice we might want match and gap-extend scores to differ.
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (1 + 30 + 45 - 6 - 2));
    REQUIRE(result.second.size() >= 3);
    REQUIRE(result.second.size() <= 4);
    if (result.second.size() == 3) {
        REQUIRE(to_score.at(result.second.at(0)).read_interval.first == 0);
        REQUIRE(to_score.at(result.second.at(1)).read_interval.first == 1);
        REQUIRE(to_score.at(result.second.at(2)).read_interval.first == 28);
    }
    if (result.second.size() == 4) {
        REQUIRE(to_score.at(result.second.at(0)).read_interval.first == 0);
        REQUIRE(to_score.at(result.second.at(1)).read_interval.first == 1);
        REQUIRE(to_score.at(result.second.at(2)).read_interval.first == 13);
        REQUIRE(to_score.at(result.second.at(3)).read_interval.first == 28);
    }
}

TEST_CASE("find_best_chain chains two extensions abutting in read and graph correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions and alignment
    auto to_score = fake_extensions({{1, 10, {h1}, 1, 9},
                                     {10, 19, {h1}, 10, 9}});
    auto aln = fake_alignment(to_score);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two extensions abutting in read with a gap in graph correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions and alignment
    auto to_score = fake_extensions({{1, 10, {h1}, 1, 9},
                                     {10, 19, {h1}, 11, 9}});
    auto aln = fake_alignment(to_score);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two extensions abutting in graph with a gap in read correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // Set up extensions and alignment
    auto to_score = fake_extensions({{1, 10, {h1}, 1, 9},
                                     {11, 20, {h1}, 10, 9}});
    auto aln = fake_alignment(to_score);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two extensions that abut over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAA");
    handle_t h2 = graph.create_handle("AAAA");
    handle_t h3 = graph.create_handle("AAAA");
    handle_t h4 = graph.create_handle("AAAA");
    handle_t h5 = graph.create_handle("AAAA");
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //             ** **** ***
    
    // Set up extensions and alignment
    auto to_score = fake_extensions({{1, 10, {h1, h2, h3}, 1, 9},
                                     {10, 19, {h3, h4, h5}, 2, 9}});
    auto aln = fake_alignment(to_score);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two extensions that abut in graph with a gap in read over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAA");
    handle_t h2 = graph.create_handle("AAAA");
    handle_t h3 = graph.create_handle("AAAA");
    handle_t h4 = graph.create_handle("AAAA");
    handle_t h5 = graph.create_handle("AAAA");
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //             ** **** ***
    
    // Set up extensions and alignment
    auto to_score = fake_extensions({{1, 10, {h1, h2, h3}, 1, 9},
                                     {11, 20, {h3, h4, h5}, 2, 9}});
    auto aln = fake_alignment(to_score);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two extensions that abut in read with a gap in graph over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAA");
    handle_t h2 = graph.create_handle("AAAA");
    handle_t h3 = graph.create_handle("AAAA");
    handle_t h4 = graph.create_handle("AAAA");
    handle_t h5 = graph.create_handle("AAAA");
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //              * **** ****
    
    // Set up extensions and alignment
    auto to_score = fake_extensions({{1, 10, {h1, h2, h3}, 1, 9},
                                     {10, 19, {h3, h4, h5}, 3, 9}});
    auto aln = fake_alignment(to_score);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two extensions that have gaps of different sizes in graph and read over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAA");
    handle_t h2 = graph.create_handle("AAAA");
    handle_t h3 = graph.create_handle("AAAA");
    handle_t h4 = graph.create_handle("AAAA");
    handle_t h5 = graph.create_handle("AAAA");
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //              * **** ****
    
    // Set up extensions and alignment
    auto to_score = fake_extensions({{1, 10, {h1, h2, h3}, 1, 9},
                                     {12, 21, {h3, h4, h5}, 3, 9}});
    auto aln = fake_alignment(to_score);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    // Gap in grapoh is 1 and gap in read is 2 so we pay 1 gap open
    REQUIRE(result.first == (9 + 9 - 6));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

TEST_CASE("find_best_chain chains two extensions that overlap the same amount in graph and read over multiple nodes correctly", "[chain_items][find_best_chain]") {
    // Set up graph fixture
    HashGraph graph;
    handle_t h1 = graph.create_handle("AAAA");
    handle_t h2 = graph.create_handle("AAAA");
    handle_t h3 = graph.create_handle("AAAA");
    handle_t h4 = graph.create_handle("AAAA");
    handle_t h5 = graph.create_handle("AAAA");
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    Aligner scoring;
    algorithms::ChainingSpace<GaplessExtension> space(scoring, &distance_index, &graph);
    
    // AAAA|AAAA|AAAA|AAAA|AAAA
    //  *** **** **
    //         * **** ****
    
    // Set up extensions and alignment
    auto to_score = fake_extensions({{1, 10, {h1, h2, h3}, 1, 9},
                                     {7, 16, {h2, h3, h4}, 3, 9}});
    auto aln = fake_alignment(to_score);
    
    // Actually run the chaining and test
    auto result = algorithms::find_best_chain(to_score, space);
    // We shouldn't charge for a gap at all here.
    // TODO: we probably want to back out the shared matches though??? But what if some are mismatches???
    REQUIRE(result.first == (9 + 9));
    REQUIRE(result.second.at(0) == 0);
    REQUIRE(result.second.at(1) == 1);
}

}

}

