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

// Make a graph with a generic X structure
static HashGraph make_x_graph() {
    // Set up graph fixture
    HashGraph graph = make_disconnected_graph(9, 10);
    // 1-3-5-6 diagonal
    graph.create_edge(graph.get_handle(2, false), graph.get_handle(4, false));
    graph.create_edge(graph.get_handle(4, false), graph.get_handle(6, false));
    graph.create_edge(graph.get_handle(6, false), graph.get_handle(7, false));
    // 2-3-4-7 diagonal
    graph.create_edge(graph.get_handle(3, false), graph.get_handle(4, false));
    graph.create_edge(graph.get_handle(4, false), graph.get_handle(5, false));
    graph.create_edge(graph.get_handle(5, false), graph.get_handle(8, false));
    // Tying off the ends
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(2, false));
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(3, false));
    graph.create_edge(graph.get_handle(7, false), graph.get_handle(9, false));
    graph.create_edge(graph.get_handle(8, false), graph.get_handle(9, false));
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
    // Provide more weight for start & end nodes to anchor graph
    std::unordered_map<nid_t, size_t> extra_node_weight;
    extra_node_weight[graph.min_node_id()] = 10000000000;
    extra_node_weight[graph.max_node_id()] = 10000000000;
    IntegratedSnarlFinder snarl_finder(graph, extra_node_weight);
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
                                                                          size_t num_chains,
                                                                          size_t max_alts = 2) {
    // Provide more weight for start & end nodes to anchor graph
    std::unordered_map<nid_t, size_t> extra_node_weight;
    extra_node_weight[graph.min_node_id()] = 10000000000;
    extra_node_weight[graph.max_node_id()] = 10000000000;
    IntegratedSnarlFinder snarl_finder(graph, extra_node_weight);
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
                                        algorithms::ChainScoringScheme(), num_chains, max_alts
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
    HashGraph graph = make_x_graph();
    auto h = get_handles(graph);

    // One anchor on each node
    // read start, graph handle and offset, length, and score
    auto to_score = make_anchors({{1, h[2], 0, 5, 5},
                                  {1, h[3], 0, 5, 5},
                                  {11, h[4], 0, 5, 5},
                                  {21, h[5], 0, 5, 5},
                                  {21, h[6], 0, 5, 5},
                                  {31, h[7], 0, 10, 10},
                                  {31, h[8], 0, 10, 10}}, graph);
    
    /// Actually run the chaining and test
    auto result = run_ziptree_iterator_multi_chain(graph, to_score, 41, 2);
    // We should see all possible paths
    REQUIRE(result.front().subchains.size() == 5);
    REQUIRE(result.front().connections.size() == 5);
}

TEST_CASE("X with different length chains", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_x_graph();
    auto h = get_handles(graph);

    // One anchor on each node
    // read start, graph handle and offset, length, and score
    auto to_score = make_anchors({{1, h[2], 0, 10, 5},
                                  {1, h[3], 0, 10, 5},
                                  {11, h[4], 0, 10, 5},
                                  {21, h[5], 0, 10, 5},
                                  {21, h[6], 0, 10, 5},
                                  {31, h[7], 0, 10, 5}}, graph);
    
    // Actually run the chaining and test
    auto result = run_ziptree_iterator_multi_chain(graph, to_score, 41, 2);
    // We should see all possible paths
    REQUIRE(result.front().subchains.size() == 5);
    REQUIRE(result.front().connections.size() == 5);
}

TEST_CASE("X with haplotype paths annotates subchains", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_x_graph();
    auto h = get_handles(graph);

    // One anchor on each node
    // read start, graph handle and offset, length, and score
    auto to_score = make_anchors({{1, h[2], 0, 5, 5},
                                  {1, h[3], 0, 5, 5},
                                  {11, h[4], 0, 5, 5},
                                  {21, h[5], 0, 5, 5},
                                  {21, h[6], 0, 5, 5},
                                  {31, h[7], 0, 10, 10},
                                  {31, h[8], 0, 10, 10}}, graph);

    // Haplotype 0 covers nodes 2 and 6, haplotype 1 covers everything but node
    // 2, so that the paths supported narrow down along the 6-7 subchain, and
    // getting from the subchain at 1 to the subchain at 3 needs a recombination.
    std::unordered_map<nid_t, algorithms::path_flags_t> haplotypes {
        {2, 1}, {3, 2}, {4, 2}, {5, 2}, {6, 3}, {7, 2}, {8, 2}
    };
    for (auto& anchor : to_score) {
        anchor.set_paths(haplotypes.at(id(anchor.graph_start())));
    }

    // Actually run the chaining and test
    auto result = run_ziptree_iterator_multi_chain(graph, to_score, 41, 2);
    REQUIRE(result.size() == 1);
    auto& group = result.front();
    REQUIRE(group.subchains.size() == 5);

    // Index the subchains by the node their first anchor is on
    std::unordered_map<nid_t, size_t> subchain_at;
    for (size_t i = 0; i < group.subchains.size(); i++) {
        subchain_at[id(to_score[group.subchains[i].anchors.front()].graph_start())] = i;
    }
    REQUIRE(subchain_at.size() == group.subchains.size());

    // The 6-7 subchain starts on both haplotypes but only haplotype 1 makes it
    // to the end, without any recombination being forced
    auto& subchain_5 = group.subchains[subchain_at.at(6)];
    REQUIRE(subchain_5.anchors.size() == 2);
    REQUIRE(subchain_5.rec_count == 0);
    REQUIRE(subchain_5.start_paths == 3);
    REQUIRE(subchain_5.end_paths == 2);

    // And connecting the subchain at 1 to the subchain at 3 needs a
    // recombination, while connecting the one at 2 to it does not
    auto& subchain_3 = group.subchains[subchain_at.at(4)];
    REQUIRE((group.subchains[subchain_at.at(2)].end_paths & subchain_3.start_paths) == 0);
    REQUIRE((group.subchains[subchain_at.at(3)].end_paths & subchain_3.start_paths) != 0);
}

TEST_CASE("single internally-recombinant anchor forces a subchain-internal recombination", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(2, 10);
    auto h = get_handles(graph);

    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    auto to_score = make_anchors({{1, h[1], 0, 5, 5},
                                  {11, h[2], 0, 5, 5}}, graph);

    to_score[0].set_paths(1);
    // Enters on haplotype 1 and leaves on haplotype 2
    to_score[1].set_paths(2, 4);

    auto result = run_ziptree_iterator_multi_chain(graph, to_score, 20, 2);
    REQUIRE(result.size() == 1);
    REQUIRE(result.front().subchains.size() == 1);

    auto& subchain = result.front().subchains.front();
    REQUIRE(subchain.anchors.size() == 2);
    REQUIRE(subchain.rec_count == 1);
    REQUIRE(subchain.start_paths == 1);
    REQUIRE(subchain.end_paths == 4);
}

TEST_CASE("subchain with no recombination has equal start and end paths", "[chain_items]") {
    // Set up graph fixture
    HashGraph graph = make_long_graph(2, 10);
    auto h = get_handles(graph);

    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    auto to_score = make_anchors({{1, h[1], 0, 5, 5},
                                  {11, h[2], 0, 5, 5}}, graph);

    for (auto& anchor : to_score) {
        anchor.set_paths(1);
    }

    auto result = run_ziptree_iterator_multi_chain(graph, to_score, 20, 2);
    REQUIRE(result.size() == 1);
    REQUIRE(result.front().subchains.size() == 1);

    auto& subchain = result.front().subchains.front();
    REQUIRE(subchain.anchors.size() == 2);
    REQUIRE(subchain.rec_count == 0);
    REQUIRE(subchain.start_paths == 1);
    REQUIRE(subchain.end_paths == 1);
}

TEST_CASE("recombination_penalty changes which predecessor chain_items_dp picks", "[chain_items]") {
    // Set up graph fixture: two alternative starts converging on one node
    HashGraph graph = make_disconnected_graph(4, 10);
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(2, false));
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(3, false));
    graph.create_edge(graph.get_handle(2, false), graph.get_handle(4, false));
    graph.create_edge(graph.get_handle(3, false), graph.get_handle(4, false));
    auto h = get_handles(graph);

    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    auto to_score = make_anchors({{1, h[2], 0, 5, 4},
                                  {1, h[3], 0, 5, 12},
                                  {11, h[4], 0, 5, 5}}, graph);

    // Index the anchors by the node they are on
    std::unordered_map<nid_t, size_t> anchor_at;
    for (size_t i = 0; i < to_score.size(); i++) {
        anchor_at[id(to_score[i].graph_start())] = i;
    }
    size_t a = anchor_at.at(2);
    size_t b = anchor_at.at(3);
    size_t d = anchor_at.at(4);

    // The anchor at 3 shares a haplotype with the low-scoring anchor at 1 but
    // not with the high-scoring anchor at 2
    to_score[a].set_paths(1);
    to_score[b].set_paths(2);
    to_score[d].set_paths(1);

    // Convert these into Seed type
    vector<SnarlDistanceIndexClusterer::Seed> seeds;
    for (const auto& pos : to_score) {
        ZipCode zipcode;
        zipcode.fill_in_zipcode_from_pos(distance_index, pos.graph_start());
        zipcode.fill_in_full_decoder();
        seeds.push_back({pos.graph_start(), 0, zipcode});
    }

    // Next, make a ZipCodeForest for the graph/seeds
    ZipCodeForest zip_forest;
    zip_forest.fill_in_forest(seeds, distance_index);
    algorithms::transition_iterator zip_it = algorithms::zip_tree_transition_iterator(
        seeds, zip_forest.trees.front(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());

    algorithms::ChainScoringScheme scheme;

    SECTION("without a penalty the higher-scoring recombinant predecessor wins") {
        scheme.recombination_penalty = 0;

        vector<vector<algorithms::TracedScore>> chain_scores;
        algorithms::chain_items_dp(chain_scores, to_score, distance_index, graph, zip_it, 3, scheme);

        REQUIRE(chain_scores[d].front().source == b);
        REQUIRE(chain_scores[d].front().score == 17);
    }

    SECTION("a large penalty makes the consistent predecessor win") {
        scheme.recombination_penalty = 20;

        vector<vector<algorithms::TracedScore>> chain_scores;
        algorithms::chain_items_dp(chain_scores, to_score, distance_index, graph, zip_it, 3, scheme);

        REQUIRE(chain_scores[d].front().source == a);
        REQUIRE(chain_scores[d].front().score == 9);
    }
}

TEST_CASE("consistency_bonus changes which predecessor chain_items_dp picks despite a lower raw score", "[chain_items]") {
    // Set up graph fixture: two alternative starts converging on one node
    HashGraph graph = make_disconnected_graph(4, 10);
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(2, false));
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(3, false));
    graph.create_edge(graph.get_handle(2, false), graph.get_handle(4, false));
    graph.create_edge(graph.get_handle(3, false), graph.get_handle(4, false));
    auto h = get_handles(graph);

    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    auto to_score = make_anchors({{1, h[2], 0, 5, 4},
                                  {1, h[3], 0, 5, 12},
                                  {11, h[4], 0, 5, 5}}, graph);

    // Index the anchors by the node they are on
    std::unordered_map<nid_t, size_t> anchor_at;
    for (size_t i = 0; i < to_score.size(); i++) {
        anchor_at[id(to_score[i].graph_start())] = i;
    }
    size_t a = anchor_at.at(2);
    size_t b = anchor_at.at(3);
    size_t d = anchor_at.at(4);

    to_score[a].set_paths(1);
    to_score[b].set_paths(2);
    to_score[d].set_paths(1);

    // Convert these into Seed type
    vector<SnarlDistanceIndexClusterer::Seed> seeds;
    for (const auto& pos : to_score) {
        ZipCode zipcode;
        zipcode.fill_in_zipcode_from_pos(distance_index, pos.graph_start());
        zipcode.fill_in_full_decoder();
        seeds.push_back({pos.graph_start(), 0, zipcode});
    }

    // Next, make a ZipCodeForest for the graph/seeds
    ZipCodeForest zip_forest;
    zip_forest.fill_in_forest(seeds, distance_index);
    algorithms::transition_iterator zip_it = algorithms::zip_tree_transition_iterator(
        seeds, zip_forest.trees.front(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());

    algorithms::ChainScoringScheme scheme;
    scheme.recombination_penalty = 0;

    SECTION("without a bonus the higher-scoring recombinant predecessor wins") {
        scheme.consistency_bonus = 0;

        vector<vector<algorithms::TracedScore>> chain_scores;
        algorithms::chain_items_dp(chain_scores, to_score, distance_index, graph, zip_it, 3, scheme);

        REQUIRE(chain_scores[d].front().source == b);
        REQUIRE(chain_scores[d].front().score == 17);
    }

    SECTION("a bonus larger than the score gap makes the consistent predecessor win") {
        scheme.consistency_bonus = 20;

        vector<vector<algorithms::TracedScore>> chain_scores;
        algorithms::chain_items_dp(chain_scores, to_score, distance_index, graph, zip_it, 3, scheme);

        REQUIRE(chain_scores[d].front().source == a);
        // The bonus only steers the choice; it stays out of the recorded score
        REQUIRE(chain_scores[d].front().score == 9);
    }
}

TEST_CASE("count_total_recombinations counts recombinations between subchains", "[chain_items]") {
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

    auto to_score = make_anchors({{1, h[1], 0, 5, 5},
                                  {1, h[2], 0, 5, 5},
                                  {11, h[3], 0, 5, 5},
                                  {21, h[4], 0, 5, 5},
                                  {21, h[5], 0, 5, 5},
                                  {31, h[6], 0, 10, 10},
                                  {31, h[7], 0, 10, 10}}, graph);

    // Haplotype 0 covers nodes 1 and 5, haplotype 1 covers everything but node 1
    std::unordered_map<nid_t, algorithms::path_flags_t> haplotypes {
        {1, 1}, {2, 2}, {3, 2}, {4, 2}, {5, 3}, {6, 2}, {7, 2}
    };
    for (auto& anchor : to_score) {
        anchor.set_paths(haplotypes.at(id(anchor.graph_start())));
    }

    auto result = run_ziptree_iterator_multi_chain(graph, to_score, 41, 2);
    REQUIRE(result.size() == 1);
    auto& group = result.front();
    REQUIRE(group.subchains.size() == 5);

    // Index the subchains by the node their first anchor is on
    std::unordered_map<nid_t, size_t> subchain_at;
    for (size_t i = 0; i < group.subchains.size(); i++) {
        subchain_at[id(to_score[group.subchains[i].anchors.front()].graph_start())] = i;
    }
    REQUIRE(subchain_at.size() == group.subchains.size());

    // Going from the subchain at 1 to the one at 3 needs a recombination
    REQUIRE(algorithms::count_total_recombinations(group, {subchain_at.at(1), subchain_at.at(3)}) == 1);
    // Going from the subchain at 2 to the one at 3 does not
    REQUIRE(algorithms::count_total_recombinations(group, {subchain_at.at(2), subchain_at.at(3)}) == 0);
}

TEST_CASE("a haplotype-consistent middle subchain still forces a recombination between incompatible ends", "[chain_items]") {
    // A only sees haplotype 1, B sees both 1 and 2 (a shared/generic node),
    // C only sees haplotype 2, D only sees haplotype 1 like A. No single
    // haplotype is consistent with both A and C, so chaining A-B-C must
    // force exactly one recombination overall, even though B itself is not
    // internally recombinant and is compatible with both neighbors
    // individually. Chaining A-B-D instead stays on haplotype 1 throughout,
    // so it must force zero recombinations.
    algorithms::Subchain a({});
    a.rec_count = 0;
    a.start_paths = 1;
    a.end_paths = 1;

    algorithms::Subchain b({});
    b.rec_count = 0;
    b.start_paths = 3; // haplotypes 1 and 2
    b.end_paths = 3;

    algorithms::Subchain c({});
    c.rec_count = 0;
    c.start_paths = 2;
    c.end_paths = 2;

    algorithms::Subchain d({});
    d.rec_count = 0;
    d.start_paths = 1;
    d.end_paths = 1;

    algorithms::SubchainGroup group;
    group.subchains = {a, b, c, d};

    REQUIRE(algorithms::count_total_recombinations(group, {0, 1, 2}) == 1);
    REQUIRE(algorithms::count_total_recombinations(group, {0, 1, 3}) == 0);
}


}

}
