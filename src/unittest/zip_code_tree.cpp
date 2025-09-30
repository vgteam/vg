#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include "../io/json2graph.hpp"
#include "../vg.hpp"
#include "catch.hpp"
#include "bdsg/hash_graph.hpp"
#include "../integrated_snarl_finder.hpp"
#include "random_graph.hpp"
#include "../minimizer_mapper.hpp"
#include <random>
#include <time.h>

//#define print

namespace vg {
namespace unittest {
    ZipCodeForest make_and_validate_forest(const vector<pos_t>& positions, const SnarlDistanceIndex& distance_index,
                                           size_t distance_limit = std::numeric_limits<size_t>::max()) {
        // Convert these into Seed type
        vector<SnarlDistanceIndexClusterer::Seed> seeds;
        for (const auto& pos : positions) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, pos);
            zipcode.fill_in_full_decoder();
            seeds.push_back({pos, 0, zipcode});
        }

        // Next, make a ZipCodeForest for the graph/seeds, and validate it
        ZipCodeForest zip_forest;
        zip_forest.fill_in_forest(seeds, distance_index, distance_limit);
        zip_forest.validate_zip_forest(distance_index, &seeds, distance_limit);

        return zip_forest;
    }
    unordered_map<ZipCodeTree::oriented_seed_t, vector<ZipCodeTree::seed_result_t>> get_reverse_views(
        const ZipCodeForest& zip_forest, size_t distance_limit = std::numeric_limits<size_t>::max()) {
        // For each seed, what seeds and distances do we see in reverse from it?
        unordered_map<ZipCodeTree::oriented_seed_t, vector<ZipCodeTree::seed_result_t>> reverse_views;
        // Follow the the usual iteration process
        for (const auto& zip_tree : zip_forest.trees) {
            for (auto seed_itr = zip_tree.begin(); seed_itr != zip_tree.end(); ++seed_itr) {
                auto dest = *seed_itr;
                
                for (auto& d: dest) {
                    reverse_views[d] = vector<ZipCodeTree::seed_result_t>();
                }
                
                for (auto dist_itr = zip_tree.find_distances(seed_itr, distance_limit);
                     !dist_itr.done(); ++dist_itr) {
                    for (const auto& d: dest) {
                        reverse_views[d].push_back(*dist_itr);
                    }
                }
            }
        }
        return reverse_views;
    }
    TEST_CASE("zip tree one node", "[zip_tree]" ) {
        VG graph;

        // Define the graph structure at the top of the test
        Node* n1 = graph.create_node("GCA");

        // Construct a distance index once the graph is built
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed") {
            // [1+0] (Section starts with the expected ziptree)

            // Define the seed positions at the top of the section
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            // Run make_and_validate_forest()
            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);

            // Finally, run any other spot checks as desired
            REQUIRE(zip_forest.trees.size() == 1);

            ZipCodeTree zip_tree = zip_forest.trees[0];
            REQUIRE(zip_tree.get_tree_size() == 3);
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);
            REQUIRE(zip_tree.get_item_at_index(0).get_value() == std::numeric_limits<size_t>::max());
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(1).get_value() == 0);
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::CHAIN_END);
            REQUIRE(zip_tree.get_item_at_index(2).get_value() == std::numeric_limits<size_t>::max());

            // We see all the seeds in order
            std::vector<ZipCodeTree::oriented_seed_t> seed_indexes = zip_tree.get_all_seeds();
            REQUIRE(seed_indexes.size() == 1);
            REQUIRE(seed_indexes.at(0).seed == 0);

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 1);
                // The only seed can't see any other seeds
                REQUIRE(reverse_views.count({0, false}));
                REQUIRE(reverse_views[{0, false}].empty());
            }
        }
        SECTION("Two seeds") {
            // [1+0 1 1+1]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);

            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];

            REQUIRE(zip_tree.get_tree_size() == 5);

            // Chain start
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);
            REQUIRE(zip_tree.get_item_at_index(0).get_value() == std::numeric_limits<size_t>::max());

            // Seed (either one)
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(1).get_value() == 0 ||
                     zip_tree.get_item_at_index(1).get_value() == 1));

            // Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(2).get_value() == 1);

            // The other seed
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(3).get_value() == 0 ||
                     zip_tree.get_item_at_index(3).get_value() == 1));

            // Chain end
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::CHAIN_END);
            REQUIRE(zip_tree.get_item_at_index(4).get_value() == std::numeric_limits<size_t>::max());

            // We see all the seeds in order
            std::vector<ZipCodeTree::oriented_seed_t> seed_indexes = zip_tree.get_all_seeds();
            REQUIRE(seed_indexes.size() == 2);
            REQUIRE(seed_indexes.at(0).seed == 0);
            REQUIRE(seed_indexes.at(1).seed == 1);

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 2);
                // The first seed can't see any other seeds
                REQUIRE(reverse_views.count({0, false}));
                REQUIRE(reverse_views[{0, false}].empty());
                
                REQUIRE(reverse_views.count({1, false}));
                REQUIRE(reverse_views[{1, false}].size() == 1);
                // The second seed can see the first seed at distance 1
                REQUIRE(reverse_views[{1, false}][0].seed == 0);
                REQUIRE(reverse_views[{1, false}][0].distance == 1);
                REQUIRE(reverse_views[{1, false}][0].is_reversed == false);
            }
        }
        SECTION("Three seeds") {
            // [1+0 1 1+1 1 1+2]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 1);
            positions.emplace_back(1, false, 2);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);

            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];

            REQUIRE(zip_tree.get_tree_size() == 7);


            // Chain start
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);
            REQUIRE(zip_tree.get_item_at_index(0).get_value() == std::numeric_limits<size_t>::max());

            // Seed
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(1).get_value() == 0);

            // Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(2).get_value() == 1);

            // The next seed
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);

            // Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(4).get_value() == 1);

            // The final seed
            REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(5).get_value() == 2);

            // Chain end
            REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);
            REQUIRE(zip_tree.get_item_at_index(6).get_value() == std::numeric_limits<size_t>::max());

            // We see all the seeds in order
            std::vector<ZipCodeTree::oriented_seed_t> seed_indexes = zip_tree.get_all_seeds();
            REQUIRE(seed_indexes.size() == 3);
            REQUIRE(seed_indexes.at(0).seed == 0);
            REQUIRE(seed_indexes.at(1).seed == 1);
            REQUIRE(seed_indexes.at(2).seed == 2);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 0);
            }

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 3);
                // The first seed can't see any other seeds
                REQUIRE(reverse_views.count({0, false}));
                REQUIRE(reverse_views[{0, false}].empty());
                
                REQUIRE(reverse_views.count({1, false}));
                REQUIRE(reverse_views[{1, false}].size() == 1);
                // The second seed can see the first seed at distance 1
                REQUIRE(reverse_views[{1, false}][0].seed == 0);
                REQUIRE(reverse_views[{1, false}][0].distance == 1);
                REQUIRE(reverse_views[{1, false}][0].is_reversed == false);
                
                REQUIRE(reverse_views.count({2, false}));
                REQUIRE(reverse_views[{2, false}].size() == 2);
                // The third seed can see both previous seeds, in reverse order
                REQUIRE(reverse_views[{2, false}][0].seed == 1);
                REQUIRE(reverse_views[{2, false}][0].distance == 1);
                REQUIRE(reverse_views[{2, false}][0].is_reversed == false);
                REQUIRE(reverse_views[{2, false}][1].seed == 0);
                REQUIRE(reverse_views[{2, false}][1].distance == 2);
                REQUIRE(reverse_views[{2, false}][1].is_reversed == false);
            }
        }
    }
    TEST_CASE("zip tree two node chain", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAAGGT");

        Edge* e1 = graph.create_edge(n1, n2);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Three seeds") {
            // [1+0 1 1+1 4 2+2]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 1);
            positions.emplace_back(2, false, 2);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];

            REQUIRE(zip_tree.get_tree_size() == 7);

            // The order should either be 0-1-2, or 2-1-0
            bool is_rev = zip_tree.get_item_at_index(1).get_value() == 2;
            if (is_rev) {
                // Chain start
                REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);
                REQUIRE(zip_tree.get_item_at_index(0).get_value() == std::numeric_limits<size_t>::max());

                // First seed 
                REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 2);
                REQUIRE(zip_tree.get_item_at_index(1).get_is_reversed() == true);

                // Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(2).get_value() == 4);

                // The next seed
                REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);
                REQUIRE(zip_tree.get_item_at_index(3).get_is_reversed() == true);

                // Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(4).get_value() == 1);

                // The last seed
                REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 0);
                REQUIRE(zip_tree.get_item_at_index(5).get_is_reversed() == true);

                // Chain end
                REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);
                REQUIRE(zip_tree.get_item_at_index(6).get_value() == std::numeric_limits<size_t>::max());
            } else {
                // Chain start
                REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);
                REQUIRE(zip_tree.get_item_at_index(0).get_value() == std::numeric_limits<size_t>::max());

                // First seed 
                REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 0);
                REQUIRE(zip_tree.get_item_at_index(1).get_is_reversed() == false);

                // Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(2).get_value() == 1);

                // The next seed
                REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);
                REQUIRE(zip_tree.get_item_at_index(3).get_is_reversed() == false);

                // Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(4).get_value() == 4);

                // The last seed
                REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 2);
                REQUIRE(zip_tree.get_item_at_index(5).get_is_reversed() == false);

                // Chain end
                REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);
                REQUIRE(zip_tree.get_item_at_index(6).get_value() == std::numeric_limits<size_t>::max());
            }
            
            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 0);
            }

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 3);
                // The first seed can't see any other seeds
                REQUIRE(reverse_views.count({0, false}));
                REQUIRE(reverse_views[{0, false}].empty());

                REQUIRE(reverse_views.count({1, false}));
                REQUIRE(reverse_views[{1, false}].size() == 1);
                // The second seed can see the first seed at distance 1
                REQUIRE(reverse_views[{1, false}][0].seed == 0);
                REQUIRE(reverse_views[{1, false}][0].distance == 1);
                REQUIRE(reverse_views[{1, false}][0].is_reversed == false);
                
                REQUIRE(reverse_views.count({2, false}));
                REQUIRE(reverse_views[{2, false}].size() == 2);
                // The third seed can see both previous seeds, in reverse order, at distances 4 and 5.
                REQUIRE(reverse_views[{2, false}][0].seed == 1);
                REQUIRE(reverse_views[{2, false}][0].distance == 4);
                REQUIRE(reverse_views[{2, false}][0].is_reversed == false);
                REQUIRE(reverse_views[{2, false}][1].seed == 0);
                REQUIRE(reverse_views[{2, false}][1].distance == 5);
                REQUIRE(reverse_views[{2, false}][1].is_reversed == false);
            }

            SECTION("Check iterator with distance limit") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest, 2);
                REQUIRE(reverse_views.size() == 3);
                // The first seed can't see any other seeds
                REQUIRE(reverse_views.count({0, false}));
                REQUIRE(reverse_views[{0, false}].empty());

                REQUIRE(reverse_views.count({1, false}));
                REQUIRE(reverse_views[{1, false}].size() == 1);
                // The second seed can see the first seed at distance 1
                REQUIRE(reverse_views[{1, false}][0].seed == 0);
                REQUIRE(reverse_views[{1, false}][0].distance == 1);
                REQUIRE(reverse_views[{1, false}][0].is_reversed == false);

                // The third seed can't see any other seeds
                REQUIRE(reverse_views.count({2, false}));
                REQUIRE(reverse_views[{2, false}].empty());
            }
        }
        SECTION("Two buckets") {
            // [1+2 1 2+0] and [2+6]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            // New tree with distance limit 4
            positions.emplace_back(2, false, 6);
            
            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 4);
            REQUIRE(zip_forest.trees.size() == 2);
        }
    }
    TEST_CASE("zip tree two two node chains", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAAGGT");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("GCAAGGT");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n3, n4);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed on each component") {
            // [3+0] and [1+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 2);
            for (auto& zip_tree : zip_forest.trees) {
                REQUIRE(zip_tree.get_tree_size() == 3);

                // Chain start
                REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

                // First seed 
                REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);

                // Chain end
                REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::CHAIN_END);
            }
                
            SECTION("Count dags") {
                for (auto& zip_tree : zip_forest.trees) {
                    pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                    REQUIRE(dag_non_dag_count.first == 0);
                    REQUIRE(dag_non_dag_count.second == 0);
                }
            }

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 2);
                // Neither seed can see any other seeds
                REQUIRE(reverse_views.count({0, false}));
                REQUIRE(reverse_views[{0, false}].empty());
                REQUIRE(reverse_views.count({1, false}));
                REQUIRE(reverse_views[{1, false}].empty());
            }
        }
        SECTION("Four seeds") {
            // [3+0 5 4+2] and [1+0 5 2+2]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 2);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 2);

            for (auto& zip_tree : zip_forest.trees) {
                REQUIRE(zip_tree.get_tree_size() == 5);

                // Chain start
                REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

                // First seed 
                REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);

                // Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(2).get_value() == 5);

                // The next seed
                REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);

                // Chain end
                REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::CHAIN_END);
            }

            SECTION("Count dags") {
                for (auto& zip_tree : zip_forest.trees) {
                    pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                    REQUIRE(dag_non_dag_count.first == 0);
                    REQUIRE(dag_non_dag_count.second == 0);
                }
            }
        }
        SECTION("Four buckets") {
            // [3+0], [4+5], [1+0], and [2+5]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 5);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 5);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 3);
            REQUIRE(zip_forest.trees.size() == 4);
        }
    }
    TEST_CASE("zip tree simple bubbles in chains", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAAGGT");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("GCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Seeds on chain nodes") {
            // [6+0rev 6 3+0rev 3 1+0rev]
            // Note that the ziptree may also be reversed
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];

            REQUIRE(zip_tree.get_tree_size() == 7);

            // Chain start
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

            // First seed 
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(1).get_is_reversed()) {
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 2);
            } else {
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 0);
            }

            // Distance between them
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(2).get_value() == 3 ||
                    zip_tree.get_item_at_index(2).get_value() == 6));

            // The next seed
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);

            // Distance between them
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(4).get_value() == 3 ||
                    zip_tree.get_item_at_index(4).get_value() == 6));

            // The last seed
            REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(5).get_is_reversed()) {
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 0);
            } else {
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 2);
            }

            // Chain end
            REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);
            
            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 0);
            }

            // We see all the seeds in order
            std::vector<ZipCodeTree::oriented_seed_t> seed_indexes = zip_tree.get_all_seeds();
            REQUIRE(seed_indexes.size() == 3);
            if (seed_indexes.at(0).is_reversed) {
                REQUIRE(seed_indexes.at(0).seed == 2);
                REQUIRE(seed_indexes.at(1).seed == 1);
                REQUIRE(seed_indexes.at(2).seed == 0);
            } else {
                REQUIRE(seed_indexes.at(0).seed == 0);
                REQUIRE(seed_indexes.at(1).seed == 1);
                REQUIRE(seed_indexes.at(2).seed == 2);    
            }

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 3);
                if (seed_indexes.at(0).is_reversed) {
                    // The first seed can't see any other seeds
                    REQUIRE(reverse_views.count({2, true}));
                    REQUIRE(reverse_views[{2, true}].empty());

                    REQUIRE(reverse_views.count({1, true}));
                    REQUIRE(reverse_views[{1, true}].size() == 1);
                    // The second seed can see the first seed at distance 6
                    REQUIRE(reverse_views[{1, true}][0].seed == 2);
                    REQUIRE(reverse_views[{1, true}][0].distance == 6);
                    REQUIRE(reverse_views[{1, true}][0].is_reversed == true);

                    REQUIRE(reverse_views.count({0, true}));
                    REQUIRE(reverse_views[{0, true}].size() == 2);
                    // The third seed can't see both the others at distances 3 and 9
                    REQUIRE(reverse_views[{0, true}][0].seed == 1);
                    REQUIRE(reverse_views[{0, true}][0].distance == 3);
                    REQUIRE(reverse_views[{0, true}][0].is_reversed == true);
                    REQUIRE(reverse_views[{0, true}][1].seed == 2);
                    REQUIRE(reverse_views[{0, true}][1].distance == 9);
                    REQUIRE(reverse_views[{0, true}][1].is_reversed == true);
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }
        }
        SECTION("Seeds on chain nodes one reversed") {
            // [6+0rev 6 3+0rev 2 1-2]
            vector<pos_t> positions;
            positions.emplace_back(1, true, 2);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];;

            REQUIRE(zip_tree.get_tree_size() == 7);

            // Chain start
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

            // First seed 
            // Either seed on 1 going backwards, or seed on 6 going backwards
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(1).get_value() == 0) {
                REQUIRE(zip_tree.get_item_at_index(1).get_is_reversed());
            } else {
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 2);
                REQUIRE(zip_tree.get_item_at_index(1).get_is_reversed());
            }

            // Distance between them
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(2).get_value() == 2 ||
                    zip_tree.get_item_at_index(2).get_value() == 6));

            // The next seed
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);

            // Distance between them
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(4).get_value() == 2 ||
                    zip_tree.get_item_at_index(4).get_value() == 6));

            // The last seed
            REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(5).get_value() == 0) {
                REQUIRE(!zip_tree.get_item_at_index(5).get_is_reversed());
            } else {
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 2);
                REQUIRE(!zip_tree.get_item_at_index(5).get_is_reversed());
            }

            // Chain end
            REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION("One seed on snarl") {
            // [6+0rev 6 (1  6  0  1 [2+1rev]) 3 1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 1);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            ZipCodeTree zip_tree = zip_forest.trees[0];
            REQUIRE(zip_tree.get_tree_size() == 15);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 0);
            }

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 3);
                if (zip_tree.get_item_at_index(1).get_is_reversed()) {
                    // The first seed can't see any other seeds
                    REQUIRE(reverse_views.count({2, true}));
                    REQUIRE(reverse_views[{2, true}].empty());

                    REQUIRE(reverse_views.count({1, true}));
                    REQUIRE(reverse_views[{1, true}].size() == 1);
                    // The second seed can see the first seed at distance 12
                    REQUIRE(reverse_views[{1, true}][0].seed == 2);
                    REQUIRE(reverse_views[{1, true}][0].distance == 12);
                    REQUIRE(reverse_views[{1, true}][0].is_reversed == true);

                    REQUIRE(reverse_views.count({0, true}));
                    REQUIRE(reverse_views[{0, true}].size() == 2);
                    // The third seed can't see both the others at distances 4 and 9
                    REQUIRE(reverse_views[{0, true}][0].seed == 1);
                    REQUIRE(reverse_views[{0, true}][0].distance == 4);
                    REQUIRE(reverse_views[{0, true}][0].is_reversed == true);
                    REQUIRE(reverse_views[{0, true}][1].seed == 2);
                    REQUIRE(reverse_views[{0, true}][1].distance == 9);
                    REQUIRE(reverse_views[{0, true}][1].is_reversed == true);
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }
        }
        SECTION("Reversed chain in snarl") {
            // [(1  1  0  6 [2-1]) 0 1-0]
            vector<pos_t> positions;
            positions.emplace_back(1, true, 0);
            positions.emplace_back(2, true, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            REQUIRE(zip_tree.get_tree_size() == 13);
            // Make sure the distance matrix is OK
            // (Reversed chains used to have inf dist-to-end bugs)
            // C1->start
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);
            // end->start
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(4).get_value() == 0);
            // end->C1
            REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(5).get_value() == 6);

            SECTION("Check iterator") {
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 2);
                if (!zip_forest.trees[0].get_item_at_index(11).get_is_reversed()) {
                    // The first seed can't see any other seeds
                    REQUIRE(reverse_views.count({1, false}));
                    REQUIRE(reverse_views[{1, false}].empty());

                    REQUIRE(reverse_views.count({0, false}));
                    REQUIRE(reverse_views[{0, false}].size() == 1);
                    // The second seed can see the first seed at distance 6
                    REQUIRE(reverse_views[{0, false}][0].seed == 1);
                    REQUIRE(reverse_views[{0, false}][0].distance == 6);
                    REQUIRE(reverse_views[{0, false}][0].is_reversed == false);
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }

            SECTION("Check iterator with distance limit") {
                auto reverse_views = get_reverse_views(zip_forest, 2);
                REQUIRE(reverse_views.size() == 2);
                if (!zip_forest.trees[0].get_item_at_index(11).get_is_reversed()) {
                    // Neither seed can see any other
                    REQUIRE(reverse_views.count({1, false}));
                    REQUIRE(reverse_views[{1, false}].empty());

                    REQUIRE(reverse_views.count({0, false}));
                    REQUIRE(reverse_views[{0, false}].empty());
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }
        }
        SECTION("Three seeds on snarl") {
            // [6+0rev 6 3+0rev 0 (1  3  0  1 [2+4rev 2 2+2rev 1 2+1rev]) 3 1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 1);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(2, false, 4);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            ZipCodeTree zip_tree = zip_forest.trees[0];
            REQUIRE(zip_tree.get_tree_size() == 21);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION("Two children of a snarl") {
            // [6+0rev 0 (2  3  2  inf  3  0  0 [4+0rev][5+1rev 1 5+0rev]) 3 3+0rev 3 1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(5, false, 1);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            ZipCodeTree zip_tree = zip_forest.trees[0];
            REQUIRE(zip_tree.get_tree_size() == 25);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 0);
            }

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 6);
                if (zip_tree.get_item_at_index(1).get_is_reversed()) {
                    // Only checking to make sure we skip the chains correctly
                    
                    // Seed right after the snarl
                    REQUIRE(reverse_views.count({1, true}));
                    // We see all the seeds to the left, not skipping any
                    REQUIRE(reverse_views[{1, true}].size() == 4);

                    // Seed in the first chain
                    REQUIRE(reverse_views.count({3, true}));
                    // We see the other seed in this chain, skip the other chain
                    // and finally see the leftmost seed
                    REQUIRE(reverse_views[{3, true}].size() == 2);
                    REQUIRE(reverse_views[{3, true}][0].seed == 4);
                    REQUIRE(reverse_views[{3, true}][0].distance == 1);
                    REQUIRE(reverse_views[{3, true}][0].is_reversed == true);
                    REQUIRE(reverse_views[{3, true}][1].seed == 5);
                    REQUIRE(reverse_views[{3, true}][1].distance == 3);
                    REQUIRE(reverse_views[{3, true}][1].is_reversed == true);
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }
        }
        SECTION("Only snarls in a chain" ) {
            // [(2  3  2  inf  3  0  0 [4+0rev][5+1rev 1 5+0rev]) 3 (1  7  0  0 [2+0rev])]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(5, false, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            ZipCodeTree zip_tree = zip_forest.trees[0];
            REQUIRE(zip_tree.get_tree_size() == 29);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 2);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION("Seeds on chain nodes bucket") {
            // [6+0rev] and [3+0rev 3 1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 4);
            REQUIRE(zip_forest.trees.size() == 2);
        }
        SECTION("Only snarls in two buckets") {
            // 0: [(2  3  2  inf  3  0  1 [4+0rev][5+1rev])]
            // 1: [(1  7  0  0 [2+0rev])]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 2);
        }
        SECTION("Snarls and nodes in three buckets") {
            // 0: [(2  3  2  inf  3  0  1 [4+0rev][5+1rev])]
            // 1: [(1  7  0  0 [2+0rev])]
            // 2: [1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 1);
            REQUIRE(zip_forest.trees.size() == 3);
        }
        SECTION("Chain in snarl in a separate bucket") {
            // [3+0rev 1 1+2rev] and [2+3rev 2+3rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 3);
            positions.emplace_back(2, false, 3);
            positions.emplace_back(3, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 2);
            REQUIRE(zip_forest.trees[0].get_tree_size() == 5);
            // Seeds on the same position should have no edge
            REQUIRE(zip_forest.trees[1].get_tree_size() == 4);
        }
        SECTION("Chain in snarl in a separate bucket another connected to end (or maybe start)") {
            // 0: [3+0rev 0 (1  7  0  0 [2+0rev]) 1 1+2rev]
            // 1: [2+3rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 3);
            positions.emplace_back(3, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 2);
        }
    }
    TEST_CASE("zip tree simple nested bubbles in chains", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAAGGT");
        Node* n3 = graph.create_node("GCAAGGT");
        Node* n4 = graph.create_node("GCAGCAAGGT");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("GCA");
        Node* n7 = graph.create_node("GGCAGCAAGGTCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n5);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Slice of snarl removed") {
            // 0: [1+0 3 (1  0  0  inf [2+0]) 0 5+0]
            // 1: [(1  6  0  1 [3+6]) 0 4+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 6);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 4);
            REQUIRE(zip_forest.trees.size() == 2);

            ZipCodeTree zip_tree = zip_forest.trees[0];
            REQUIRE(zip_tree.get_tree_size() == 15);
            // The inf edge
            REQUIRE(zip_tree.get_item_at_index(7).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(7).get_value() == std::numeric_limits<size_t>::max());

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 5);
                // Only checking to make sure 5+0 skips the snarl's chain
                if (!zip_tree.get_item_at_index(1).get_is_reversed()) {
                    REQUIRE(reverse_views.count({4, false}));
                    // Go straight to 1+0
                    REQUIRE(reverse_views[{4, false}].size() == 1);
                    REQUIRE(reverse_views[{4, false}][0].seed == 0);
                    REQUIRE(reverse_views[{4, false}][0].distance == 3);
                    REQUIRE(reverse_views[{4, false}][0].is_reversed == false);
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }
        }
    }
    TEST_CASE("zip tree bubble in cyclic snarl", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCA");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GCAAAAAAAAA");
        Node* n6 = graph.create_node("GCA");
        Node* n7 = graph.create_node("GGCAAAAAAAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n6);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n2, n5, true, true);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Two sides of nested snp unordered along read") {
            // 0: [[1+0 3 {1  inf  3  inf  inf  17  inf  0  inf  inf  inf [(2  0  0  inf  3  3  3 [3+0][4+0])]}]
            // 1: [5+5 5+5]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(5, false, 5);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 5);
            positions.emplace_back(3, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 4);
        }
    }
    TEST_CASE("zip tree bubble nested in inversion", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCAAAAAAAAAA");
        Node* n2 = graph.create_node("GCAA");
        Node* n3 = graph.create_node("GCAG");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GCAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n4, false, true);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n2, n5, true, false);
        Edge* e6 = graph.create_edge(n3, n4);
        Edge* e7 = graph.create_edge(n4, n5);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Traverse nested chain forwards but the orientation of the chain is backwards in the snarl tree") {
            // [1+0 5 1+5 7 {1  inf  3  inf  0  inf  inf  7  3  0  inf
            //     [4+0rev 0 (1  4  0  0 [3+0rev]) 4 2+0rev]} 0 5+0 4 5+4]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 5);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(5, false, 4);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
        }
    }
    TEST_CASE("zip tree bubble nested in cyclic snarl", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCAAAAAAAAAA");
        Node* n2 = graph.create_node("GCAA");
        Node* n3 = graph.create_node("GCAG");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GCAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n5);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n5, n5, true, false);
        Edge* e6 = graph.create_edge(n3, n4);
        Edge* e7 = graph.create_edge(n4, n5);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Traverse nested chain forwards but the orientation of the chain is backwards in the snarl tree") {
            // [1+0 5 1+5 7 {1  inf  0  inf  inf  inf  inf  0  inf  3  inf
            //     [2+0 4 (1  0  0  4 [3+0]) 0 4+0]} 0 5+0 4 5+4]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 5);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(5, false, 4);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
        }
    }
    TEST_CASE("zip tree snarl with inversion", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAA");
        Node* n3 = graph.create_node("GCAGGT");
        Node* n4 = graph.create_node("GC");
        Node* n5 = graph.create_node("GCCCCCCCCCCCCCCCCCCCC");

        Edge* e1 = graph.create_edge(n1, n2, false, true);
        Edge* e2 = graph.create_edge(n1, n4);
        Edge* e3 = graph.create_edge(n2, n3, true, false);
        Edge* e4 = graph.create_edge(n3, n4, false, true);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n3, n5, true, false);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Traverse 3 backwards") {
            // [1+0 3 {2  inf  0  inf  12  inf  inf  9  inf  inf  inf  2  inf
            //     2  inf  inf  8  inf  8  5  0  inf [4+0][3-1rev 1 3-0rev]} 0 5+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(3, true, 0);
            positions.emplace_back(3, true, 1);
            positions.emplace_back(5, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            bool chain_is_reversed = distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id()));
            if (chain_is_reversed) {
                cerr << "This test didn't get run because I'm lazy and didn't write it for a reversed chain" << endl;           
            } else {
                // Check some random elements

                // First seed
                REQUIRE(zip_forest.trees[0].get_item_at_index(1).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(1).get_value() == 0);
                // Chain count
                REQUIRE(zip_forest.trees[0].get_item_at_index(4).get_type() == ZipCodeTree::CHAIN_COUNT);
                // Chain start
                REQUIRE(zip_forest.trees[0].get_item_at_index(26).get_type() == ZipCodeTree::CHAIN_START);
                // Second seed (4)
                REQUIRE(zip_forest.trees[0].get_item_at_index(27).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(27).get_value() == 1);
                // Third seed (3-1)
                REQUIRE(zip_forest.trees[0].get_item_at_index(30).get_type() == ZipCodeTree::SEED);
                // Second chain within snarl may be reversed
                if (zip_forest.trees[0].get_item_at_index(30).get_value() == 2) {
                    REQUIRE(!zip_forest.trees[0].get_item_at_index(30).get_is_reversed());
                } else {
                    REQUIRE(zip_forest.trees[0].get_item_at_index(30).get_is_reversed());
                    REQUIRE(zip_forest.trees[0].get_item_at_index(30).get_value() == 3);
                }
            }

            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                // All five seeds go R->L,
                // and the three in the cyclic snarl also go L->R
                REQUIRE(reverse_views.size() == 8);
                if (!chain_is_reversed) {
                    // 4+0 R->L can only leave the snarl; it skips all chains
                    REQUIRE(reverse_views.count({1, false}));
                    REQUIRE(reverse_views[{1, false}].size() == 1);
                    // Edge to 1+0
                    REQUIRE(reverse_views[{1, false}][0].seed == 0);
                    REQUIRE(reverse_views[{1, false}][0].distance == 3);
                    REQUIRE(reverse_views[{1, false}][0].is_reversed == false);

                    // 4+0 L->R can see the other chain & outside the snarl
                    REQUIRE(reverse_views.count({1, true}));
                    REQUIRE(reverse_views[{1, true}].size() == 3);
                    // Edge to 3-0rev
                    REQUIRE(reverse_views[{1, true}][0].seed == 2);
                    REQUIRE(reverse_views[{1, true}][0].distance == 2);
                    REQUIRE(reverse_views[{1, true}][0].is_reversed == true);
                    // Edge to 3-1rev
                    REQUIRE(reverse_views[{1, true}][1].seed == 3);
                    REQUIRE(reverse_views[{1, true}][1].distance == 3);
                    REQUIRE(reverse_views[{1, true}][1].is_reversed == true);
                    // Edge to 5+0rev (yes, rev - we're going L->R here)
                    REQUIRE(reverse_views[{1, true}][2].seed == 4);
                    REQUIRE(reverse_views[{1, true}][2].distance == 8);
                    REQUIRE(reverse_views[{1, true}][2].is_reversed == true);

                    // 3-1 can see the rest of its chain & the other chain
                    // Not rev since we're going L->R
                    REQUIRE(reverse_views.count({3, false}));
                    REQUIRE(reverse_views[{3, false}].size() == 2);
                    // Edge to 3-0
                    REQUIRE(reverse_views[{3, false}][0].seed == 2);
                    REQUIRE(reverse_views[{3, false}][0].distance == 1);
                    REQUIRE(reverse_views[{3, false}][0].is_reversed == false);
                    // Edge to 4+0
                    REQUIRE(reverse_views[{3, false}][1].seed == 1);
                    REQUIRE(reverse_views[{3, false}][1].distance == 3);
                    REQUIRE(reverse_views[{3, false}][1].is_reversed == false);

                    // 5+0 can see all other seeds once
                    REQUIRE(reverse_views.count({4, false}));
                    REQUIRE(reverse_views[{4, false}].size() == 4);
                    // Edge to 3-1 (not rev since going L->R)
                    REQUIRE(reverse_views[{4, false}][0].seed == 3);
                    REQUIRE(reverse_views[{4, false}][0].distance == 5);
                    REQUIRE(reverse_views[{4, false}][0].is_reversed == false);
                    // Edge to 3-0
                    REQUIRE(reverse_views[{4, false}][1].seed == 2);
                    REQUIRE(reverse_views[{4, false}][1].distance == 6);
                    REQUIRE(reverse_views[{4, false}][1].is_reversed == false);
                    // Edge to 4+0
                    REQUIRE(reverse_views[{4, false}][2].seed == 1);
                    REQUIRE(reverse_views[{4, false}][2].distance == 8);
                    REQUIRE(reverse_views[{4, false}][2].is_reversed == false);
                    // Edge to 1+0
                    REQUIRE(reverse_views[{4, false}][3].seed == 0);
                    REQUIRE(reverse_views[{4, false}][3].distance == 11);
                    REQUIRE(reverse_views[{4, false}][3].is_reversed == false);
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }
        }
    }
    TEST_CASE("zip tree non-simple DAG", "[zip_tree]") {
        // bubble between 1 and 3, non-simple dag between 3 and 8 
        // containing node 7 and chain 4-6
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAA");
        Node* n3 = graph.create_node("GCAGGT");
        Node* n4 = graph.create_node("GC");
        Node* n5 = graph.create_node("GC");
        Node* n6 = graph.create_node("GCA");
        Node* n7 = graph.create_node("GCA");
        Node* n8 = graph.create_node("GCAGGGGGGGGGGGGAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n7);
        Edge* e6 = graph.create_edge(n3, n8);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n4, n6);
        Edge* e9 = graph.create_edge(n5, n6);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n7, n8);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Make the zip tree") {
            // [1+0 3 (1  0  0  4 [2+0]) 0 3+0 1 3+1 5 (2  0  1  4  0  6  2
            //     [4+0 2 (1  0  0  2 [5+0]) 0 6+0][7+1]) 0 8+0 2 8+2]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, false, 1);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(8, false, 0);
            positions.emplace_back(8, false, 2);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];

            bool chain_is_reversed = distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id()));
            if (!chain_is_reversed) {
                // Check some random elements

                // First seed
                REQUIRE(zip_forest.trees[0].get_item_at_index(1).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(1).get_value() == 0);
                // Start of snarl
                REQUIRE(zip_forest.trees[0].get_item_at_index(17).get_type() == ZipCodeTree::SNARL_START);
                
                // Parts of distance matrix
                REQUIRE(zip_forest.trees[0].get_item_at_index(20).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_forest.trees[0].get_item_at_index(20).get_value() == 1);
                REQUIRE(zip_forest.trees[0].get_item_at_index(21).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_forest.trees[0].get_item_at_index(21).get_value() == 4);

                REQUIRE(zip_forest.trees[0].get_item_at_index(34).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(34).get_value() == 5);

                REQUIRE(zip_forest.trees[0].get_item_at_index(36).get_type() == ZipCodeTree::SNARL_END);
            }

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 3);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION("Three buckets") {
            // 0: [1+0 3 (1  0  0  4 [2+0]) 0 3+0]
            // 1: [(2  0  1  4  0  6  2 [4+0 2 (1  0  0  2 [5+0]) 0 6+0][7+1]) 0 8+0]
            // 2: [8-0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(8, false, 0);
            positions.emplace_back(8, true, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 3);
            REQUIRE(zip_forest.trees.size() == 3);
        }
    }
    TEST_CASE("zip tree three-chain DAG") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAA");
        Node* n3 = graph.create_node("GCAGGT");
        Node* n4 = graph.create_node("GC");
        Node* n5 = graph.create_node("GC");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n1, n4);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed on each node") {
            // [5+0rev 0 (3  2  6  6  6  4  inf  2  0  0  0
            //     [4+0rev][3+0rev][2+0rev]) 3 1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            
            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            SECTION("Check iterator") {
                auto reverse_views = get_reverse_views(zip_forest);
                REQUIRE(reverse_views.size() == 5);
                if (zip_forest.trees[0].get_item_at_index(1).get_is_reversed()) {
                    // 2+0rev skips 3+0rev's chain but sees the rest on the left
                    REQUIRE(reverse_views.count({1, true}));
                    REQUIRE(reverse_views[{1, true}].size() == 2);
                    // Edge to 4+0rev
                    REQUIRE(reverse_views[{1, true}][0].seed == 3);
                    REQUIRE(reverse_views[{1, true}][0].distance == 4);
                    REQUIRE(reverse_views[{1, true}][0].is_reversed == true);
                    // Edge to 5+0rev
                    REQUIRE(reverse_views[{1, true}][1].seed == 4);
                    REQUIRE(reverse_views[{1, true}][1].distance == 6);
                    REQUIRE(reverse_views[{1, true}][1].is_reversed == true);
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }
        }
    }
    TEST_CASE("zip tree deeply nested bubbles", "[zip_tree]") {
        // top-level chain 1-12-13-16
        // bubble 2-10 containing two bubbles 3-5 and 6-9
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCA");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GAC");
        Node* n6 = graph.create_node("GCA");
        Node* n7 = graph.create_node("GCA");
        Node* n8 = graph.create_node("GCA");
        Node* n9 = graph.create_node("GCA");
        Node* n10 = graph.create_node("GCA");
        Node* n11 = graph.create_node("GCA");
        Node* n12 = graph.create_node("GCA");
        Node* n13 = graph.create_node("GCA");
        Node* n14 = graph.create_node("GCA");
        Node* n15 = graph.create_node("GCA");
        Node* n16 = graph.create_node("GCGGGGGGGGGGGGGGGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n11);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n6);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n10);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n6, n8);
        Edge* e11 = graph.create_edge(n7, n9);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n10, n12);
        Edge* e15 = graph.create_edge(n11, n12);
        Edge* e16 = graph.create_edge(n12, n13);
        Edge* e17 = graph.create_edge(n13, n14);
        Edge* e18 = graph.create_edge(n13, n15);
        Edge* e19 = graph.create_edge(n14, n16);
        Edge* e20 = graph.create_edge(n15, n16);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Make the zip tree with a seed on each node") {
            // [1+0 3 (2  2  0  inf  3  1  1 [11+2][2+0 3 (2  0  0  inf  6  1  3
            //     [6+0 3 (2  1  0  inf  3  2  3 [7+1][8+0]) 2 9+2][3+0 3 (1  0  0  3 [4+0]) 0 5+0])
            //     2 10+2]) 2 12+2 3 13+2 1 (2  2  2  inf  3  1  1 [14+2][15+2]) 2 16+2]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(8, false, 0);
            positions.emplace_back(9, false, 2);
            positions.emplace_back(10, false, 2);
            positions.emplace_back(11, false, 2);
            positions.emplace_back(12, false, 2);
            positions.emplace_back(13, false, 2);
            positions.emplace_back(14, false, 2);
            positions.emplace_back(15, false, 2);
            positions.emplace_back(16, false, 2);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 5);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION("Make the zip tree with a few seeds") {
            // [1+0 3 (1  3  3  3 [(2  0  0  inf  6  9  3 [6+0][3+0 3 5+0])])
            //     5 13+2 1 (1  2  3  1 [15+2])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);
            positions.emplace_back(13, false, 2);
            positions.emplace_back(15, false, 2);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_forest.trees[0].dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 3);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION("3 buckets") {
            // [1+2 1 (1  9  3  3 [10+0])], [13+2], and [16+5]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(10, false, 0);
            positions.emplace_back(13, false, 2);
            positions.emplace_back(16, false, 5);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 4);
            REQUIRE(zip_forest.trees.size() == 3);
        }
        SECTION("Remove empty snarls") {
            // 0: [1+2]
            // 1: [(1  1  0  2 [4+1])]
            // 2: [(1  1  6  3 [6+1 2 (1  1  3  2 [7+1])])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(6, false, 1);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(4, false, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 3);
        }
        SECTION("Chain connected on one end") {
            // 0: [1+2 1 (1  0  3  3 [2+0 2 2+2 1 (1  1  6  3 [6+1 2 (1  1  3  2 [7+1])])])]
            // 1: [(1  1  0  2 [4+1])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(6, false, 1);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(4, false, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 2);
        }
        SECTION("Chain connected on the other end") {
            // 0: [1+2 1 (1  3  3  1 [(1  3  6  2 [(1  1  3  2 [7+1]) 1 9+1]) 0 10+0 2 10+2])]
            // 1: [(1  1  0  2 [4+1])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(10, false, 0);
            positions.emplace_back(10, false, 2);
            positions.emplace_back(9, false, 1);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(4, false, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 2);
        }
        SECTION("One chain removed from a snarl") {
            // 0: [1+2 1 (1  1  3  2 [11+1])]
            // 1: [(2  1  1  inf  3  2  2 [7+1][8+1])]
            // 2: [(1  0  0  3 [4+0])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(8, false, 1);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(11, false, 1);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 3);
        }
    }
    TEST_CASE("zip tree long nested chain", "[zip_tree]") {
        // top-level chain 1-12-13-16
        // bubble 2-10 containing two bubbles 3-5 and 6-9
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCA");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GAC");
        Node* n6 = graph.create_node("GCA");
        Node* n7 = graph.create_node("GCA");
        Node* n8 = graph.create_node("GCA");
        Node* n9 = graph.create_node("GCA");
        Node* n10 = graph.create_node("GCA");
        Node* n11 = graph.create_node("GCA");
        Node* n12 = graph.create_node("GCA");
        Node* n13 = graph.create_node("GCA");
        Node* n14 = graph.create_node("GCA");
        Node* n15 = graph.create_node("GCA");
        Node* n16 = graph.create_node("GCG");
        Node* n17 = graph.create_node("GCA");
        Node* n18 = graph.create_node("GCA");
        Node* n19 = graph.create_node("GCA");
        Node* n20 = graph.create_node("GCA");
        Node* n21 = graph.create_node("GCA");
        Node* n22 = graph.create_node("GCA");
        Node* n23 = graph.create_node("GCAAAAAAAAAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n2, n14);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n7, n9);
        Edge* e12 = graph.create_edge(n8, n10);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n10, n11);
        Edge* e15 = graph.create_edge(n10, n12);
        Edge* e16 = graph.create_edge(n11, n12);
        Edge* e17 = graph.create_edge(n12, n13);
        Edge* e18 = graph.create_edge(n13, n21);
        Edge* e19 = graph.create_edge(n14, n15);
        Edge* e20 = graph.create_edge(n14, n16);
        Edge* e21 = graph.create_edge(n15, n16);
        Edge* e22 = graph.create_edge(n16, n17);
        Edge* e23 = graph.create_edge(n16, n20);
        Edge* e24 = graph.create_edge(n17, n18);
        Edge* e25 = graph.create_edge(n17, n19);
        Edge* e26 = graph.create_edge(n18, n19);
        Edge* e27 = graph.create_edge(n19, n20);
        Edge* e28 = graph.create_edge(n20, n21);
        Edge* e29 = graph.create_edge(n21, n22);
        Edge* e30 = graph.create_edge(n21, n23);
        Edge* e31 = graph.create_edge(n22, n23);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One slice from nodes in the middle of a nested chain") {
            // 0: [1+0 3 2+0 3 (2  0  0  inf  9  3  3 [14+0 3 16+0 3 20+0][3+0 inf 13+0]) 0 21+0]
            // 1: [10+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(10, false, 0);
            positions.emplace_back(13, false, 0);
            positions.emplace_back(21, false, 0);
            positions.emplace_back(14, false, 0);
            positions.emplace_back(16, false, 0);
            positions.emplace_back(20, false, 0);  

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 3);
            REQUIRE(zip_forest.trees.size() == 2);

        }
        SECTION("Two slices from snarls in the middle of a nested chain") {
            // 0: [1+2 1 2+0]
            // 1: [21+0]
            // 2: [(1  0  3  3 [4+0]) 1 6+1 2 (1  0  3  6 [7+0])]
            // 3: [(1  0  0  3 [11+0]) 0 12+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(6, false, 1);
            positions.emplace_back(7, false, 0);
            positions.emplace_back(11, false, 0);
            positions.emplace_back(12, false, 0);
            positions.emplace_back(21, false, 0);  

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 4);
        }
        SECTION("One slice from the start of a chain, connected to the end") {
            // 0: [1+2 1 2+0 3 (1  16  9  3 [12+1 2 13+0]) 0 21+0]
            // 1: [(1  0  3  6 [7+0])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(7, false, 0);
            positions.emplace_back(12, false, 1);
            positions.emplace_back(13, false, 0);
            positions.emplace_back(21, false, 0);    

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 3);
            REQUIRE(zip_forest.trees.size() == 2);
        }
        SECTION("One slice from the end of a chain, connected to the start") {
            // 0: [1+2 1 2+0 3 (2  0  0  inf  9  3  inf [14+0 3 16+0 3 20+0][3+0]) 0 21+0]
            // 1: [(1  0  3  6 [7+0])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(7, false, 0);
            positions.emplace_back(14, false, 0);    
            positions.emplace_back(16, false, 0);    
            positions.emplace_back(20, false, 0);    
            positions.emplace_back(21, false, 0);    

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 3);
            REQUIRE(zip_forest.trees.size() == 2);
        }
    }
    TEST_CASE("zip tree non-dag", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCA");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GAC");
        Node* n6 = graph.create_node("GCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3, false, true);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Make the zip tree with a seed on each node") {
            // [1+0 3 {2  inf  0  inf  6  inf  inf  0  inf  inf  inf  6  inf  6  inf
            //     inf  3  inf  3  inf  3  inf [3+0][2+0]} 0 4+0 3 (1  0  0  3 [5+0]) 0 6+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_forest.trees[0].dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 1);
            }
        }
    }
    TEST_CASE("zip tree nested cyclic non-dag", "[zip_tree]") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCA");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GAC");
        Node* n6 = graph.create_node("AAAAAAAAAAAAAAAGCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n5);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n3);
        Edge* e7 = graph.create_edge(n4, n2);
        Edge* e8 = graph.create_edge(n4, n5);
        Edge* e9 = graph.create_edge(n5, n6);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Make the zip tree with a seed on each node") {
            // [1+0 3 {1  inf  0  inf  inf  3  inf  0  inf  3  inf [2+0 3
            //     {1  inf  0  inf  inf  3  inf  0  inf  3  inf [3+0]} 0 4+0]} 0 5+0 3 6+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_forest.trees[0].dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 2);
            }
        }

    }
    TEST_CASE("zip tree nested inversions", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCA");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("GCA");
        Node* n5 = graph.create_node("GAC");
        Node* n6 = graph.create_node("AAAAAAAAAAAAAAAGCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n4, false, true);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n3, false, true);
        Edge* e5 = graph.create_edge(n2, n5, true, false);
        Edge* e6 = graph.create_edge(n3, n4);
        Edge* e7 = graph.create_edge(n3, n4, true, false);
        Edge* e8 = graph.create_edge(n4, n5);
        Edge* e9 = graph.create_edge(n5, n6);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Go forward through the inversions") {
            // [1+0 3 {1  inf  3  inf  0  inf  inf  9  3  0  inf [4+0rev 0
            //     {1  inf  0  inf  2  inf  inf  3  0  2  inf [3+0 1 3+1]} 3 2+0rev]} 0 5+0 3 6+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, false, 1);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];

            assert(zip_tree.get_tree_size() == 45);

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 2);
            }

            SECTION("Check iterator memorization") {
                auto reverse_views = get_reverse_views(zip_forest);
                // 5+0 should only see 3+1 once due to memorization
                REQUIRE(reverse_views.count({5, false}));
                size_t seen_3 = 0;
                for (auto& view : reverse_views[{5, false}]) {
                    if (view.seed == 3) seen_3++;
                }
                REQUIRE(seen_3 == 1);
            }
        }
        SECTION("Reverse both inversions") {
            // [1+0 3 {1  inf  0  inf  3  inf  inf  9  0  3  inf [4-0 3
            //     {1  inf  0  inf  2  inf  inf  3  0  2  inf [3+0 1 3+1]} 0 2-0]} 0 5+0 3 6+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(4, true, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, false, 1);
            positions.emplace_back(2, true, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];

            SECTION("Count dags") {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_cyclic_snarl_count();
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 2);
            }
        }
    }
    TEST_CASE("zip tree cyclic snarl with overlapping seeds", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCAAAAAAAAAAAAAAAAAAAAAA");
        Node* n2 = graph.create_node("AAAGCA");
        Node* n3 = graph.create_node("GCAAAA");
        Node* n4 = graph.create_node("GCAAAA");
        Node* n5 = graph.create_node("GACAAAAAAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3, false, true);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5, true, false);
        Edge* e6 = graph.create_edge(n4, n5);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Cyclic snarl with seeds on either side") {
            // [5+4rev 4 {3  inf  2  inf  inf  inf  inf  0  inf  inf  inf  inf
            //     inf  inf  inf  inf  inf  inf  inf  inf  inf  inf  8  inf  2
            //     inf  4  inf  inf  12  inf  6  0  8  0  8  inf
            //     [4+4rev 4+4rev 2 4+2rev 4+2rev 2 4+0rev 4+0rev]
            //     [3+0 3+0 2 3+2 3+2 2 3+4 3+4][2+0 2+0 2 2+2 2+2 2 2+4 2+4]}
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(2, false, 4);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(2, false, 4);

            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, false, 2);
            positions.emplace_back(3, false, 4);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, false, 2);
            positions.emplace_back(3, false, 4);

            positions.emplace_back(4, false, 0);
            positions.emplace_back(4, false, 2);
            positions.emplace_back(4, false, 4);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(4, false, 2);
            positions.emplace_back(4, false, 4);
            positions.emplace_back(5, false, 4);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
        }
    }
    TEST_CASE("zip tree duplication", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCAAAAAAAAAAAAAAAAAAAAAA");
        Node* n2 = graph.create_node("AAAGCAAAAAA");
        Node* n3 = graph.create_node("GACAAAAAAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n2);
        Edge* e3 = graph.create_edge(n2, n3);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Cyclic snarl with seeds on either side") {
            // [3+0rev 0 {1  inf  9  inf  inf  9  inf  11  inf  0  inf 
            //     [2+2rev 2+2rev 1 2+1rev 2+1rev 1 2+0rev 2+0rev]} 24 1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 1);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 1);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(3, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            ZipCodeTree zip_tree = zip_forest.trees[0];
            // Make sure that duplicated seeds get collapsed on each other
            REQUIRE(zip_tree.get_tree_size() == 29);
            auto seeds_in_order = zip_tree.get_all_seeds();
            
            SECTION("Check iterator") {
                // For each seed, what seeds and distances do we see in reverse from it?
                auto reverse_views = get_reverse_views(zip_forest);
                // All eight seeds go R->L,
                // and the six in the cyclic snarl also go L->R
                REQUIRE(reverse_views.size() == 14);
                 if (zip_tree.get_item_at_index(1).get_is_reversed()) {
                    // Checking that middle seed can loop around
                    REQUIRE(reverse_views[seeds_in_order[3]].size() == 9);
                    // Edge to 2+2rev
                    REQUIRE(reverse_views[seeds_in_order[3]][0].seed == seeds_in_order[2].seed);
                    REQUIRE(reverse_views[seeds_in_order[3]][0].distance == 1);
                    REQUIRE(reverse_views[seeds_in_order[3]][0].is_reversed == true);
                    REQUIRE(reverse_views[seeds_in_order[3]][1].seed == seeds_in_order[1].seed);
                    REQUIRE(reverse_views[seeds_in_order[3]][1].distance == 1);
                    REQUIRE(reverse_views[seeds_in_order[3]][1].is_reversed == true);
                    // Edge to 2+0rev (loop around)
                    REQUIRE(reverse_views[seeds_in_order[3]][2].seed == seeds_in_order[6].seed);
                    REQUIRE(reverse_views[seeds_in_order[3]][2].distance == 10);
                    REQUIRE(reverse_views[seeds_in_order[3]][2].is_reversed == true);
                    REQUIRE(reverse_views[seeds_in_order[3]][3].seed == seeds_in_order[5].seed);
                    REQUIRE(reverse_views[seeds_in_order[3]][3].distance == 10);
                    REQUIRE(reverse_views[seeds_in_order[3]][3].is_reversed == true);
                    // Edge to 2+1rev (self-loop)
                    REQUIRE(reverse_views[seeds_in_order[3]][4].seed == seeds_in_order[4].seed);
                    REQUIRE(reverse_views[seeds_in_order[3]][4].distance == 11);
                    REQUIRE(reverse_views[seeds_in_order[3]][4].is_reversed == true);
                    REQUIRE(reverse_views[seeds_in_order[3]][5].seed == seeds_in_order[3].seed);
                    REQUIRE(reverse_views[seeds_in_order[3]][5].distance == 11);
                    REQUIRE(reverse_views[seeds_in_order[3]][5].is_reversed == true);
                    // Edge to 3+0rev
                    REQUIRE(reverse_views[seeds_in_order[3]][8].seed == seeds_in_order[0].seed);
                    REQUIRE(reverse_views[seeds_in_order[3]][8].distance == 10);
                    REQUIRE(reverse_views[seeds_in_order[3]][8].is_reversed == true);
                } else {
                    cerr << "Not testing reverse views since I didn't bother writing it" << endl;
                }
            }
        }
    }
    TEST_CASE("zip tree self loops", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCAAAAAAAAAAAAAAAAAAAAAA");
        Node* n2 = graph.create_node("AAAGCAAAAAA");
        Node* n3 = graph.create_node("TT");
        Node* n4 = graph.create_node("GACAAAAAAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n2, true, false);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n2, false, true);
        Edge* e5 = graph.create_edge(n2, n4, true, false);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One position") {
            // [{1  inf  0  0  24  24  24  24  0  24  inf [2+0 0 2-11rev]}]
            vector<pos_t> positions;
            // Same position, but going in either direction
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, true, 11);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            SECTION("Check iterator") {
                auto reverse_views = get_reverse_views(zip_forest);
                // Both seeds have two entries, going R->L and L->R
                REQUIRE(reverse_views.size() == 4);

                // 2+0 sees 2-ll and then itself
                REQUIRE(reverse_views[{0, false}].size() == 1);
                // Edge to 2-11rev taking C1L->C1L self-loop is ignored
                // (due to having distance 0)
                // Edge to self now circling back around
                REQUIRE(reverse_views[{0, false}][0].seed == 0);
                REQUIRE(reverse_views[{0, false}][0].distance == 24);
                REQUIRE(reverse_views[{0, false}][0].is_reversed == false);

                // 2+0rev sees self and then 2-11rev
                REQUIRE(reverse_views[{0, true}].size() == 2);
                // Edge to self taking C1R->C1L self-loop
                REQUIRE(reverse_views[{0, true}][0].seed == 0);
                REQUIRE(reverse_views[{0, true}][0].distance == 24);
                REQUIRE(reverse_views[{0, true}][0].is_reversed == true);
                // Edge to 2-11rev now circling back around
                REQUIRE(reverse_views[{0, true}][1].seed == 1);
                REQUIRE(reverse_views[{0, true}][1].distance == 24);
                REQUIRE(reverse_views[{0, true}][1].is_reversed == true);
            }
        }
        SECTION("Cyclic snarl with seeds on either side") {
            // [4+0rev 0 {1  inf  0  0  22  22  20  24  0  22  inf 
            //     [2+0 0 2-11rev 1 2+1 0 2-10rev 1 2+2 0 2-9rev]} 24 1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 1);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(2, true, 11);
            positions.emplace_back(2, true, 10);
            positions.emplace_back(2, true, 9);
            positions.emplace_back(4, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            // Check self-loop distances
            // c1_left -> c1_left
            REQUIRE(zip_forest.trees[0].get_item_at_index(7).get_value() == 0);
            // c1_left -> c1_right
            REQUIRE(zip_forest.trees[0].get_item_at_index(9).get_value() == 22);
            // c1_right -> c1_right
            REQUIRE(zip_forest.trees[0].get_item_at_index(10).get_value() == 20);
        }
        SECTION("Duplicate seed with reversed in between") {
            // [{1  inf  0  0  24  24  24  24  0  24  inf [2+0 2+0 0 2-11rev]}]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, true, 11);
            positions.emplace_back(2, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            REQUIRE(zip_tree.get_tree_size() == 21);
            // Check items in chain
            REQUIRE(zip_tree.get_item_at_index(14).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(14).get_value() == 0);
            REQUIRE(zip_tree.get_item_at_index(15).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(15).get_value() == 2);
            REQUIRE(zip_tree.get_item_at_index(16).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(16).get_value() == 0);
            REQUIRE(zip_tree.get_item_at_index(17).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(17).get_value() == 1);
        }
    }
    TEST_CASE("zip tree handles complicated nested snarls", "[zip_tree]") {
        // Load an example graph
        VG graph;
        io::json2graph(R"({"node":[{"id": "1","sequence":"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"},
                                   {"id":"2","sequence":"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"},
                                   {"id":"3","sequence":"T"},
                                   {"id":"4","sequence":"T"},
                                   {"id":"5","sequence":"ATATCTATACATATAATACAG"},
                                   {"id":"6","sequence":"AT"},
                                   {"id":"7","sequence":"T"},
                                   {"id":"8","sequence":"A"},
                                   {"id":"9","sequence":"C"},
                                   {"id":"10","sequence":"AT"},
                                   {"id":"11","sequence":"A"},
                                   {"id":"12","sequence":"C"}],
                           "edge":[{"from":"3","to":"10"},
                                   {"from":"4","to":"5"},
                                   {"from":"5","to":"11"},
                                   {"from":"6","to":"7"},
                                   {"from":"7","to":"11"},
                                   {"from":"7","to":"12","to_end":true},
                                   {"from":"7","to":"8"},
                                   {"from":"8","to":"4"},
                                   {"from":"9","to":"10"},
                                   {"from":"11","to":"3"},
                                   {"from":"11","to":"9"},
                                   {"from":"12","from_start":true,"to":"3"},
                                   {"from":"1","to":"6"},
                                   {"from":"10","to":"2"}]})", &graph);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Three seeds") {
            // [6+0 3 (1  3  2  22 [5+1]) 1 10+1]
            vector<pos_t> positions;
            positions.emplace_back(6, false, 0);
            positions.emplace_back(5, false, 1);
            positions.emplace_back(10, false, 1);
            
            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);
        }
    }
    TEST_CASE("Root snarl", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GTGCACA");//8
        Node* n2 = graph.create_node("GTGCACA");
        Node* n3 = graph.create_node("GT");
        Node* n4 = graph.create_node("GATTCTTATAG");//11

        Edge* e1 = graph.create_edge(n1, n3);
        Edge* e2 = graph.create_edge(n1, n4);
        Edge* e3 = graph.create_edge(n3, n2);
        Edge* e4 = graph.create_edge(n3, n4, false, true);
        Edge* e5 = graph.create_edge(n2, n4);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed on each node") {
            // ([3-0rev][4+0][2+0][1+0])
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, true, 0);
            positions.emplace_back(4, false, 0);
            
            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index);
            REQUIRE(zip_forest.trees.size() == 1);

            SECTION("Check iterator") {
                auto reverse_views = get_reverse_views(zip_forest);
                // None of the four seeds can see anything
                REQUIRE(reverse_views.size() == 4);
                for (auto& rv : reverse_views) {
                    REQUIRE(rv.second.size() == 0);
                }
            }
        }
        SECTION("Splice out chain") {
            // ([1+6]) and [1+3]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 3);
            positions.emplace_back(1, false, 6);
            
            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 2);
            REQUIRE(zip_forest.trees.size() == 2);
        }
    }
    TEST_CASE("One nested dag snarl", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("TGTTTAAGGCTCGATCATCCGCTCACAGTCCGTCGTAGACGCATCAGACTTGGTTTCCCAAGC");
        Node* n2 = graph.create_node("G");
        Node* n3 = graph.create_node("A");
        Node* n4 = graph.create_node("CTCGCGG");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("ACCAGGCAGAATCGAGGGATGTTC");
        Node* n7 = graph.create_node("AACAGTGTCCAACACTGG");

        //Inversion
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n7);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n6, n7); 

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed in nested snarl and one outside") {
            // [7+17rev 17 (1  24  1  8 [(1  1  0  0 [5+0rev])])]
            vector<pos_t> positions;
            positions.emplace_back(5, false, 0);
            positions.emplace_back(7, false, 17);
            
            make_and_validate_forest(positions, distance_index, 61);
        }
    }
    TEST_CASE("Components of root", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GTGCACA");//8
        Node* n2 = graph.create_node("GTGAAAAAAAAAAAAAAACACA");
        Node* n3 = graph.create_node("AAAAAAAAAAAAGT");
        Node* n4 = graph.create_node("GATTCTTATAG");//11
        Node* n5 = graph.create_node("GATTCTTATAG");//11

        //Inversion
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n2, false, true);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n3, true, false);
        
        ofstream out ("testGraph.hg");
        graph.serialize(out);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Five buckets") {
            // 0: [1+0 3 1+3 2 1+5 2 {1  inf  0  inf  inf  inf  inf  22  0  inf  inf [2+0]}]
            // 1: [2+7 2 2+9 1 2+10]
            // 2: [3-3rev]
            // 3: [5+0]
            // 4: [4+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 3);
            positions.emplace_back(1, false, 5);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 7);
            positions.emplace_back(2, false, 9);
            positions.emplace_back(2, false, 10);
            positions.emplace_back(3, true, 3);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            
            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 5);
            REQUIRE(zip_forest.trees.size() == 5);
        }
    }
    TEST_CASE("Another non-dag snarl", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GTG");
        Node* n2 = graph.create_node("G");
        Node* n3 = graph.create_node("A");
        Node* n4 = graph.create_node("GAAAAAAAAT");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("G");
        Node* n7 = graph.create_node("GAAAAAAAAAT");
        Node* n8 = graph.create_node("GAT");
        Node* n9 = graph.create_node("GATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n7, false, true);
        Edge* e6 = graph.create_edge(n4, n8, true, false);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n4, n6);
        Edge* e9 = graph.create_edge(n5, n7);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n8, n9);
        
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Multiple seeds in snarl") {
            // [{3  inf  0  inf  23  inf  22  12  inf  11  inf  11  inf  10  inf 
            //     inf  0  inf  inf  inf  inf  inf  24  inf  23  inf  11  inf
            //     inf  23  inf  22  11  10  inf  23  inf [3+0 1 3-0rev]
            //     [(2  0  0  inf  1  1  1 [6-0][5-0])][2+0]}]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, true,  0);
            positions.emplace_back(5, true,  0);
            positions.emplace_back(6, true,  0);
            
            make_and_validate_forest(positions, distance_index);
        }
    }
    TEST_CASE("Remove snarl and then a chain slice", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GTG");
        Node* n2 = graph.create_node("GTG");
        Node* n3 = graph.create_node("AAA");
        Node* n4 = graph.create_node("GAT");
        Node* n5 = graph.create_node("GAAT");
        Node* n6 = graph.create_node("GATAAAAA");
        Node* n7 = graph.create_node("GAT");
        Node* n8 = graph.create_node("GAT");
        Node* n9 = graph.create_node("GAT");
        Node* n10 = graph.create_node("GAT");
        Node* n11 = graph.create_node("GATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n11);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n7, n9);
        Edge* e12 = graph.create_edge(n8, n10);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n10, n11);
        
        //ofstream out ("testGraph.hg");
        //graph.serialize(out);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Node first") {
            // [(1  0  0  3 [2+0 inf 10+0])], [5+0], and [6+4]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 4);
            positions.emplace_back(10, false, 0);
            
            make_and_validate_forest(positions, distance_index, 3);
        }
        SECTION("Snarl first") {
            // [(1  3  0  3 [(1  0  3  3 [3+0]) 10 10+0])], [6+4]
            vector<pos_t> positions;
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 4);
            positions.emplace_back(10, false, 0);
            
            make_and_validate_forest(positions, distance_index, 3);
        }
    }
    TEST_CASE("Remove a child of the top-level chain", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GTGGGGGGG");
        Node* n2 = graph.create_node("GGGGGGGTG");
        Node* n3 = graph.create_node("GGGGGGAAA");
        Node* n4 = graph.create_node("GGGGGGGAT");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
       
        //ofstream out ("testGraph.hg");
        //graph.serialize(out);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One tree on each node") {
            // [(1  7  0  2 [2+7]) 3 3+3] and [4+7]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 7);
            positions.emplace_back(3, false, 3);
            positions.emplace_back(4, false, 7);
            
            make_and_validate_forest(positions, distance_index, 3);
        }
        SECTION("Remove second child of snarl") {
            // [3+8] and [4+5]
            vector<pos_t> positions;
            positions.emplace_back(3, false, 8);
            positions.emplace_back(4, false, 5);
            
            make_and_validate_forest(positions, distance_index, 3);
        }
    }
    TEST_CASE("Remove a child of the top-level snarl", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GTGGGGGGG");
        Node* n2 = graph.create_node("GGGGGGGTG");
        Node* n3 = graph.create_node("GGGGGGAAA");
        Node* n4 = graph.create_node("GGGGGGGAT");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n4, false, true);

        ofstream out ("testGraph.hg");
        graph.serialize(out);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One tree on each node") {
            // [4+5], [1+5], [2+5], and [3+5]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 5);
            positions.emplace_back(2, false, 5);
            positions.emplace_back(3, false, 5);
            positions.emplace_back(4, false, 5);
            
            make_and_validate_forest(positions, distance_index, 3);
        }
        SECTION("Remove second child of snarl") {
            // ([3+8]) and [4+5]
            vector<pos_t> positions;
            positions.emplace_back(3, false, 8);
            positions.emplace_back(4, false, 5);
            
            make_and_validate_forest(positions, distance_index, 3);
        }
        SECTION("Remove first child of snarl") {
            // ([4+0]) and [3+5]
            vector<pos_t> positions;
            positions.emplace_back(3, false, 5);
            positions.emplace_back(4, false, 0);
            
            make_and_validate_forest(positions, distance_index, 3);
        }
        SECTION("Remove one chain") {
            // [4+4]
            vector<pos_t> positions;
            positions.emplace_back(4, false, 4);
            
            make_and_validate_forest(positions, distance_index, 3);
        }
    }
    TEST_CASE("Snp nested in looping snarl", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GTGGGGGGG");
        Node* n2 = graph.create_node("GGGGGGGTG");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("G");
        Node* n5 = graph.create_node("GGGGGGGAT");
        Node* n6 = graph.create_node("GGGGGGGAT");
        Node* n7 = graph.create_node("GGGGGGGATTTTTTTTTTTTTTTTTTTTTT");
        Node* n8 = graph.create_node("GGGGGGGAT");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n5);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n5, n6);
        Edge* e7 = graph.create_edge(n6, n2);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n1, n8);
        Edge* e10 = graph.create_edge(n8, n7);
       
        //ofstream out ("testGraph.hg");
        //graph.serialize(out);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("Snps alone") {
            // [1+0 9 {1  inf  8  inf  inf  26  inf  9  inf  18  inf
            //     [2+8 2+8 1 (2  0  0  inf  1  1  1 [3+0][4+0]) 0 5+0<6/3>]} 0 7+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 8);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(2, false, 8);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(7, false, 0);
            
            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 100);
        }
    }
    TEST_CASE("zipcode tree simple chain with multiple connected components", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("T");
        Node* n8 = graph.create_node("TTTTTTTTT");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n7);
        Edge* e8 = graph.create_edge(n6, n7);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One cluster on the same node plus extra node") {
            // [8+3] and [4+3rev 2 4+1rev 1 4+0rev]
            vector<pos_t> positions;
            positions.emplace_back(4, false, 0);
            positions.emplace_back(4, false, 1);
            positions.emplace_back(4, false, 3);
            positions.emplace_back(8, false, 3);

            make_and_validate_forest(positions, distance_index, 100);
        }
    }
    TEST_CASE("zipcode tree multicomponent chain nested in irregular snarl", "[zip_tree]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCAAAAAAAAAAAAAAAAAAAAAAAAA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("T");
        Node* n8 = graph.create_node("TTTTTTTTT");
        Node* n9 = graph.create_node("TTTTTTTTT");
        Node* n10 = graph.create_node("GCAAAAAAAAAAAAA");
        Node* n11 = graph.create_node("TTT");
        Node* n12 = graph.create_node("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
        Node* n13 = graph.create_node("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n12);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n10);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7, true, false);
        Edge* e9 = graph.create_edge(n7, n8);
        Edge* e10 = graph.create_edge(n7, n9);
        Edge* e11 = graph.create_edge(n8, n9);
        Edge* e12 = graph.create_edge(n9, n11);
        Edge* e13 = graph.create_edge(n10, n11);
        Edge* e14 = graph.create_edge(n10, n10, false, true);
        Edge* e15 = graph.create_edge(n11, n12);
        Edge* e16 = graph.create_edge(n12, n13);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        ofstream out ("testGraph.hg");
        graph.serialize(out);

        SECTION("Cross unreachable chain") {
            // [{1  inf  1  inf  inf  inf  inf  0  inf  3  inf
            //    [{1  inf  0  inf  inf  inf  inf  15  inf  9  inf
            //    [3+0 1 (1  0  0  inf [4+0]) 1 5+1 inf 7+0 1 (1  0  0  9 [8+0]) 0 9+0]}]}]
            vector<pos_t> positions;
            positions.emplace_back(n3->id(), false, 0);
            positions.emplace_back(n4->id(), false, 0);
            positions.emplace_back(n5->id(), false, 1);
            positions.emplace_back(n7->id(), false, 0);
            positions.emplace_back(n8->id(), false, 0);
            positions.emplace_back(n9->id(), false, 0);

            make_and_validate_forest(positions, distance_index, 100);
        }
        SECTION("Cross unreachable chain including snarl that is not start-end reachable") {
            // 0: [{1  inf  1  inf  inf  inf  inf  0  inf  3  inf
            //        [{1  inf  0  inf  inf  inf  inf  15  inf  9  inf 
            //        [3+0 1 (1  0  0  inf [4+0]) 1 5+1 inf 7+0 1 (1  0  0  9 [8+0]) 0 9+0]}]}]
            // 1: [(1  inf  inf  0 [6+0rev])]
            vector<pos_t> positions;
            positions.emplace_back(n3->id(), false, 0);
            positions.emplace_back(n4->id(), false, 0);
            positions.emplace_back(n5->id(), false, 1);
            positions.emplace_back(n6->id(), false, 0);
            positions.emplace_back(n7->id(), false, 0);
            positions.emplace_back(n8->id(), false, 0);
            positions.emplace_back(n9->id(), false, 0);

            make_and_validate_forest(positions, distance_index, 100);
        }
    }
    /*
    This test case will "pass" if you run it, but the forest will be weird
    because the snarl finder gets very confused by the looping chain,
    so the snarl/chain decomposition is quite odd.

         3             <-- the graph looks like this, but the distance index
        / \                thinks that 1/5 are in an irregular snarl,
       2 - 4               with 2 outside; then it makes an irregular snarl
      /     \              with 3 inside and 4 outside in a separate section.
    1 ------- 5
    TEST_CASE("Looping chain zipcode tree", "[zip_tree]") {
        // chain 2rev->2rev
        VG graph;

        Node* n1 = graph.create_node("ACACGTTGC");
        Node* n2 = graph.create_node("TCTCCACCGGCAAGTTTCACTTCACTT");
        Node* n3 = graph.create_node("A");
        Node* n4 = graph.create_node("AT");
        Node* n5 = graph.create_node("CGTGGGG");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n5);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);

        SECTION("One cluster on the same node plus extra node") {
            // 0: [2+0rev 0 (2  inf  9  9  inf  0  inf [5+0rev][1+0rev])]
            // 1: [4+0rev 0 (1  1  0  0 [3+0rev])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);

            make_and_validate_forest(positions, distance_index, 100);
        }
    }
    */
    TEST_CASE("ziptree with inversion inside of duplication", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("AAA");
        Node* n2 = graph.create_node("C");
        Node* n3 = graph.create_node("GAT");
        Node* n4 = graph.create_node("AT");
        Node* n5 = graph.create_node("CCC");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n2, n3, false, true);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n4, true, false);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n2);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed in inversion") {
            // [{1  inf  7  3  1  6  9  6  2  8  inf
            //     [{1  inf  0  3  3  6  9  3  0  3  inf [3+0]}]}]
            vector<pos_t> positions;
            positions.emplace_back(3, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed outside inversion") {
            // [{1  inf  7  9  4  6  3  6  2  11  inf [4+0rev]}]
            vector<pos_t> positions;
            positions.emplace_back(4, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed inside and one seed outside inversion") {
            // [{1  inf  7  9  1  6  9  6  2  8  inf [4+0rev 0
            //     {1  inf  0  3  3  6  9  3  0  3  inf [3+0]}]}]
            vector<pos_t> positions;
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);

            vector<SnarlDistanceIndexClusterer::Seed> seeds;

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed on either side of inversion") {
            // [{1  inf  7  9  0  2  7  6  2  7  inf [4+0rev 4 2+0rev]}]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);

            vector<SnarlDistanceIndexClusterer::Seed> seeds;

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed on each node") {
            // [1+0 3 {1  inf  7  9  0  2  7  6  2  7  inf [4+0rev 0
            //     {1  inf  0  3  3  6  9  3  0  3  inf [3+0]} 1 2+0rev]} 0 5+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);

            vector<SnarlDistanceIndexClusterer::Seed> seeds;

            make_and_validate_forest(positions, distance_index);
        }
    }
    TEST_CASE("ziptree with duplications sharing a start point", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("AA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("CGC");
        Node* n4 = graph.create_node("GT");
        Node* n5 = graph.create_node("AA");
        Node* n6 = graph.create_node("ACAC");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n3, n4);
        Edge* e4 = graph.create_edge(n4, n5);
        Edge* e5 = graph.create_edge(n5, n6);
        // Backtracks to n2
        Edge* e6 = graph.create_edge(n4, n2);
        Edge* e7 = graph.create_edge(n5, n2);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed on each node") {
            // [1+0 2 {2  inf  6  inf  inf  8  inf  0  inf  2  inf  inf  2  inf
            //     2  inf  8  inf  2  inf  4  inf [5+0][2+0 1 3+0 3 4+0]} 0 6+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
    }
    TEST_CASE("ziptree with properly nested duplications", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("AA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("CGC");
        Node* n4 = graph.create_node("GT");
        Node* n5 = graph.create_node("AA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n3, n4);
        Edge* e4 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n2);
        Edge* e7 = graph.create_edge(n3, n3);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed in inner duplication") {
            // [{1  inf  1  inf  inf  0  inf  6  inf  2  inf
            //     [{1  inf  0  inf  inf  3  inf  3  inf  3  inf [3+0]}]}]
            vector<pos_t> positions;
            positions.emplace_back(3, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed in each duplication") {
            // [{1  inf  0  inf  inf  2  inf  6  inf  2  inf
            //     [2+0 1 {1  inf  0  inf  inf  3  inf  3  inf  3  inf [3+0]}]}]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed on either side of inner duplication") {
            // [{1  inf  0  inf  inf  2  inf  6  inf  2  inf [2+0 4 4+0]}]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed on each node") {
            // [1+0 2 {1  inf  0  inf  inf  2  inf  6  inf  2  inf [2+0 1
            //     {1  inf  0  inf  inf  3  inf  3  inf  3  inf [3+0]} 0 4+0]} 0 5+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
    }
    TEST_CASE("ziptree with duplication around an insertion", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("AA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("CGC");
        Node* n4 = graph.create_node("GT");
        Node* n5 = graph.create_node("AA");

        // Main chain
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n3, n5);
        // Insertion
        Edge* e5 = graph.create_edge(n2, n4);
        Edge* e6 = graph.create_edge(n4, n3);
        // Duplication
        Edge* e7 = graph.create_edge(n3, n2);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed inside insertion") {
            // [{1  inf  1  inf  inf  4  inf  4  inf  3  inf [(1  0  0  2 [4+0])]}]
            vector<pos_t> positions;
            positions.emplace_back(4, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed in insertion and one right outside") {
            // [{1  inf  0  inf  inf  3  inf  4  inf  3  inf
            //     [2+0 1 (1  2  0  0 [4-0rev])]}]
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, true, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed in insertion and two right outside") {
            // [{1  inf  1  inf  inf  2  inf  4  inf  1  inf
            //     [(1  2  0  0 [4-0rev]) 0 3+0 2 3+2]}]
            vector<pos_t> positions;
            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, false, 2);
            positions.emplace_back(4, true, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("One seed on each node") {
            // [1+0 2 {1  inf  0  inf  inf  3  inf  4  inf  3  inf
            //     [2+0 1 (1  0  0  2 [4+0]) 0 3+0]} 0 5+0]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
    }
    TEST_CASE("ziptree with wacky cyclic snarl stuff", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("AA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("CGC");
        Node* n4 = graph.create_node("GT");
        Node* n5 = graph.create_node("ACAC");
        Node* n6 = graph.create_node("AA");

        // A regular snarl
        Edge* e1 = graph.create_edge(n2, n3);
        Edge* e2 = graph.create_edge(n2, n4);
        Edge* e3 = graph.create_edge(n3, n5);
        Edge* e4 = graph.create_edge(n4, n5);
        // + a duplication
        Edge* e5 = graph.create_edge(n5, n2);
        // Edges leave snarl
        Edge* e6 = graph.create_edge(n1, n2);
        Edge* e7 = graph.create_edge(n5, n6);
        // Reversion
        Edge* e8 = graph.create_edge(n6, n6);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed on each node") {
            // [2+0rev 0 {2  inf  2  inf  inf  inf  inf  inf  inf  inf  inf  inf
            //     inf  inf  2  inf  0  inf  inf  inf  0  inf [1+0rev][6+0rev]}
            //     4 5+0rev 0 (2  3  2  inf  2  0  0 [3+0rev][4+0rev])]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
    }
    TEST_CASE("Snarl within child chain", "[zip_tree]") {
        VG graph;

        Node* n1 = graph.create_node("AAAAAAAAAAAAA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("CGCTTTTGA");
        Node* n4 = graph.create_node("C");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("AAAAAAAAAAAA");

        // Inner snarl (insertion)
        Edge* e1 = graph.create_edge(n2, n3);
        Edge* e2 = graph.create_edge(n2, n4);
        Edge* e3 = graph.create_edge(n3, n5);
        Edge* e4 = graph.create_edge(n4, n5);
        // Outer snarl (deletion)
        Edge* e5 = graph.create_edge(n1, n2);
        Edge* e6 = graph.create_edge(n1, n6);
        Edge* e7 = graph.create_edge(n5, n6);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION("One seed on each node") {
            // [6+0rev 0 (1  1  0  0 [5+0rev 0 (2  5  1  inf  1  4  0 
            //     [3+4rev][4+0rev]) 1 2+0rev]) 13 1+0rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 4);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);

            make_and_validate_forest(positions, distance_index);
        }
        SECTION("Snip out inner snarl") {
            // [6+0rev 0 (1  3  0  0 [2+0rev]) 1 1+12rev] and [3+4rev]
            vector<pos_t> positions;
            positions.emplace_back(1, false, 12);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 4);
            positions.emplace_back(6, false, 0);

            ZipCodeForest zip_forest = make_and_validate_forest(positions, distance_index, 3);
            REQUIRE(zip_forest.trees.size() == 2);
        }
    }
    TEST_CASE("Random graphs zip tree", "[zip_tree][zip_tree_random]") {
        for (int i = 0; i < 10; i++) {
            // For each random graph
    
            default_random_engine generator(time(NULL));
            uniform_int_distribution<int> variant_count(1, 10);
            uniform_int_distribution<int> chrom_len(10, 200);
            uniform_int_distribution<int> distance_limit(5, 100);
    
            // Make a random graph with three chromosomes of random lengths
            HashGraph graph;
            random_graph({chrom_len(generator),chrom_len(generator),chrom_len(generator)}, 30, variant_count(generator), &graph);
            graph.serialize("testGraph.hg");

            IntegratedSnarlFinder snarl_finder(graph);
            SnarlDistanceIndex distance_index;
            fill_in_distance_index(&distance_index, &graph, &snarl_finder);

            vector<id_t> all_nodes;
            graph.for_each_handle([&](const handle_t& h)->bool{
                id_t id = graph.get_id(h);
                all_nodes.push_back(id);
                return true;
            });

            uniform_int_distribution<int> randPosIndex(0, all_nodes.size()-1);

            // Check k random sets of seeds
            for (size_t k = 0; k < 10 ; k++) {
                vector<pos_t> positions;

                uniform_int_distribution<int> randPosCount(3, 70);
                for (int j = 0; j < randPosCount(generator); j++) {
                    // Check clusters of j random positions

                    id_t nodeID1 = all_nodes[randPosIndex(generator)];
                    handle_t node1 = graph.get_handle(nodeID1);

                    offset_t offset1 = uniform_int_distribution<int>(0,graph.get_length(node1) - 1)(generator);

                    positions.emplace_back(nodeID1,
                                           uniform_int_distribution<int>(0,1)(generator) == 0,
                                           offset1);
                }
                size_t limit = distance_limit(generator);

                make_and_validate_forest(positions, distance_index, limit);
                REQUIRE(true); // Just to count
            }
        }
    }
    /*
    TEST_CASE("Failed zip tree unit test", "[failed]") {
        // Load failed random graph
        HashGraph graph;
        graph.deserialize("testGraph.hg");
        // print with vg view -j testGraph.hg

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        vector<pos_t> positions;
        // add seeds as needed

        make_and_validate_forest(positions, distance_index);
    }
    */
}
}
