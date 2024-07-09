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

    TEST_CASE( "zip tree one node",
                   "[zip_tree]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One seed" ) {
 
            id_t seed_nodes[] = {1};
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            REQUIRE(zip_tree.get_tree_size() == 3);
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(1).get_value() == 0);
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::CHAIN_END);

            // We see all the seeds in order
            std::vector<ZipCodeTree::oriented_seed_t> seed_indexes;
            std::copy(zip_tree.begin(), zip_tree.end(), std::back_inserter(seed_indexes));
            REQUIRE(seed_indexes.size() == 1);
            REQUIRE(seed_indexes.at(0).seed == 0);

            // For each seed, what seeds and distances do we see in reverse from it?
            std::unordered_map<ZipCodeTree::oriented_seed_t, std::vector<ZipCodeTree::seed_result_t>> reverse_views;
            for (auto forward = zip_tree.begin(); forward != zip_tree.end(); ++forward) {
                std::copy(zip_tree.look_back(forward), zip_tree.rend(), std::back_inserter(reverse_views[*forward]));
            }
            REQUIRE(reverse_views.size() == 1);
            // The only seed can't see any other seeds
            REQUIRE(reverse_views.count({0, false}));
            REQUIRE(reverse_views[{0, false}].size() == 0);
        }

        SECTION( "Two seeds" ) {
 
            id_t seed_nodes[] = {1, 1};
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            REQUIRE(zip_tree.get_tree_size() == 5);


            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

            //Seed (either one because they're the same position)
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(1).get_value() == 0 ||
                     zip_tree.get_item_at_index(1).get_value() == 1));

            //Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(2).get_value() == 0);

            //THe other seed
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(3).get_value() == 0 ||
                     zip_tree.get_item_at_index(3).get_value() == 1));

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::CHAIN_END);

            // We see all the seeds in order
            std::vector<ZipCodeTree::oriented_seed_t> seed_indexes;
            std::copy(zip_tree.begin(), zip_tree.end(), std::back_inserter(seed_indexes));
            REQUIRE(seed_indexes.size() == 2);
            REQUIRE(seed_indexes.at(0).seed == 0);
            REQUIRE(seed_indexes.at(1).seed == 1);

            // For each seed, what seeds and distances do we see in reverse from it?
            std::unordered_map<ZipCodeTree::oriented_seed_t, std::vector<ZipCodeTree::seed_result_t>> reverse_views;
            for (auto forward = zip_tree.begin(); forward != zip_tree.end(); ++forward) {
                std::copy(zip_tree.look_back(forward), zip_tree.rend(), std::back_inserter(reverse_views[*forward]));
            }
            REQUIRE(reverse_views.size() == 2);
            // The first seed can't see any other seeds
            REQUIRE(reverse_views.count({0, false}));
            REQUIRE(reverse_views[{0, false}].size() == 0);
            // The second seed can see the first seed at distance 0
            REQUIRE(reverse_views.count({1, false}));
            REQUIRE(reverse_views[{1, false}].size() == 1);
            REQUIRE(reverse_views[{1, false}][0].seed == 0);
            REQUIRE(reverse_views[{1, false}][0].distance == 0);
            REQUIRE(reverse_views[{1, false}][0].is_reverse == false);
        }

        SECTION( "Three seeds" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 2);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            REQUIRE(zip_tree.get_tree_size() == 7);


            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

            //Seed (either one because they're the same position)
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(1).get_value() == 0 ||
                     zip_tree.get_item_at_index(1).get_value() == 1));

            //Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(2).get_value() == 0);

            //THe other seed
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(3).get_value() == 0 ||
                     zip_tree.get_item_at_index(3).get_value() == 1));

            //Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(4).get_value() == 2);

            //The other seed
            REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(5).get_value() == 2);

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);

            // We see all the seeds in order
            std::vector<ZipCodeTree::oriented_seed_t> seed_indexes;
            std::copy(zip_tree.begin(), zip_tree.end(), std::back_inserter(seed_indexes));
            REQUIRE(seed_indexes.size() == 3);
            REQUIRE(seed_indexes.at(0).seed == 0);
            REQUIRE(seed_indexes.at(1).seed == 1);
            REQUIRE(seed_indexes.at(2).seed == 2);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 0);
            }

            // For each seed, what seeds and distances do we see in reverse from it?
            std::unordered_map<ZipCodeTree::oriented_seed_t, std::vector<ZipCodeTree::seed_result_t>> reverse_views;
            for (auto forward = zip_tree.begin(); forward != zip_tree.end(); ++forward) {
                std::copy(zip_tree.look_back(forward), zip_tree.rend(), std::back_inserter(reverse_views[*forward]));
            }
            REQUIRE(reverse_views.size() == 3);
            // The first seed can't see any other seeds
            REQUIRE(reverse_views.count({0, false}));
            REQUIRE(reverse_views[{0, false}].size() == 0);
            // The second seed can see the first seed at distance 0
            REQUIRE(reverse_views.count({1, false}));
            REQUIRE(reverse_views[{1, false}].size() == 1);
            REQUIRE(reverse_views[{1, false}][0].seed == 0);
            REQUIRE(reverse_views[{1, false}][0].distance == 0);
            REQUIRE(reverse_views[{1, false}][0].is_reverse == false);
            // The third seed can see both previous seeds, in reverse order, at distance 2.
            REQUIRE(reverse_views.count({2, false}));
            REQUIRE(reverse_views[{2, false}].size() == 2);
            REQUIRE(reverse_views[{2, false}][0].seed == 1);
            REQUIRE(reverse_views[{2, false}][0].distance == 2);
            REQUIRE(reverse_views[{2, false}][0].is_reverse == false);
            REQUIRE(reverse_views[{2, false}][1].seed == 0);
            REQUIRE(reverse_views[{2, false}][1].distance == 2);
            REQUIRE(reverse_views[{2, false}][1].is_reverse == false);
        }
    }
    TEST_CASE( "zip tree two node chain", "[zip_tree]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAAGGT");

        Edge* e1 = graph.create_edge(n1, n2);


        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "Three seeds" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(1, false, 1);
            positions.emplace_back(2, false, 2);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            REQUIRE(zip_tree.get_tree_size() == 7);

            //The order should either be 0-1-2, or 2-1-0
            bool is_rev = zip_tree.get_item_at_index(1).get_value() == 2;
            if (is_rev) {

                //Chain start
                REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

                //first seed 
                REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 2);
                REQUIRE(zip_tree.get_item_at_index(1).get_is_reversed() == true);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(2).get_value() == 4);

                //The next seed
                REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);
                REQUIRE(zip_tree.get_item_at_index(3).get_is_reversed() == true);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(4).get_value() == 1);

                //The last seed
                REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 0);
                REQUIRE(zip_tree.get_item_at_index(5).get_is_reversed() == true);

                //Chain end
                REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);
            } else {

                //Chain start
                REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

                //first seed 
                REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 0);
                REQUIRE(zip_tree.get_item_at_index(1).get_is_reversed() == false);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(2).get_value() == 1);

                //The next seed
                REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);
                REQUIRE(zip_tree.get_item_at_index(3).get_is_reversed() == false);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(4).get_value() == 4);

                //The last seed
                REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 2);
                REQUIRE(zip_tree.get_item_at_index(5).get_is_reversed() == false);

                //Chain end
                REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);
            }
            
            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 0);
            }

            // For each seed, what seeds and distances do we see in reverse from it?
            std::unordered_map<ZipCodeTree::oriented_seed_t, std::vector<ZipCodeTree::seed_result_t>> reverse_views;
            for (auto forward = zip_tree.begin(); forward != zip_tree.end(); ++forward) {
                std::copy(zip_tree.look_back(forward), zip_tree.rend(), std::back_inserter(reverse_views[*forward]));
            }
            REQUIRE(reverse_views.size() == 3);
            // The first seed can't see any other seeds
            REQUIRE(reverse_views.count({0, false}));
            REQUIRE(reverse_views[{0, false}].size() == 0);
            // The second seed can see the first seed at distance 1
            REQUIRE(reverse_views.count({1, false}));
            REQUIRE(reverse_views[{1, false}].size() == 1);
            REQUIRE(reverse_views[{1, false}][0].seed == 0);
            REQUIRE(reverse_views[{1, false}][0].distance == 1);
            REQUIRE(reverse_views[{1, false}][0].is_reverse == false);
            // The third seed can see both previous seeds, in reverse order, at distances 4 and 5.
            REQUIRE(reverse_views.count({2, false}));
            REQUIRE(reverse_views[{2, false}].size() == 2);
            REQUIRE(reverse_views[{2, false}][0].seed == 1);
            REQUIRE(reverse_views[{2, false}][0].distance == 4);
            REQUIRE(reverse_views[{2, false}][0].is_reverse == false);
            REQUIRE(reverse_views[{2, false}][1].seed == 0);
            REQUIRE(reverse_views[{2, false}][1].distance == 5);
            REQUIRE(reverse_views[{2, false}][1].is_reverse == false);
        }

        SECTION( "Two buckets" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            //New tree with distance limit 4
            positions.emplace_back(2, false, 6);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 4);
            REQUIRE(zip_forest.trees.size() == 2);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }


    }
    TEST_CASE( "zip tree two two node chains", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "One seed on each component" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 2);
            zip_forest.print_self(&seeds, &minimizers);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);

                //The tree should be:
                // [pos1] [pos3]
                REQUIRE(zip_tree.get_tree_size() == 3);

                //Chain start
                REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

                //first seed 
                REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);

                //Chain end
                REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::CHAIN_END);

            }
                
            SECTION( "Count dags" ) {
                for (auto& zip_tree : zip_forest.trees) {
                    pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                    REQUIRE(dag_non_dag_count.first == 0);
                    REQUIRE(dag_non_dag_count.second == 0);
                }
            }
            //TODO: This doesn't work now that it is a forest
            
            // For each seed, what seeds and distances do we see in reverse from it?
            std::unordered_map<ZipCodeTree::oriented_seed_t, std::vector<ZipCodeTree::seed_result_t>> reverse_views;
            for (auto& zip_tree : zip_forest.trees) {
                for (auto forward = zip_tree.begin(); forward != zip_tree.end(); ++forward) {
                    std::copy(zip_tree.look_back(forward), zip_tree.rend(), std::back_inserter(reverse_views[*forward]));
                }
            }
            REQUIRE(reverse_views.size() == 2);
            // Neither seed can see any other seeds
            REQUIRE(reverse_views.count({0, false}));
            REQUIRE(reverse_views[{0, false}].size() == 0);
            REQUIRE(reverse_views.count({1, false}));
            REQUIRE(reverse_views[{1, false}].size() == 0);
        }
        SECTION( "Four seeds" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 2);

            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 2);

            zip_forest.print_self(&seeds, &minimizers);


                //The tree should be:
                // [pos1 5 pos2] [pos3 5 pos4]
                // or
                // [pos2 5 pos1] [ pos3 5 pos4]
                // etc...
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
                REQUIRE(zip_tree.get_tree_size() == 5);

                //Chain start
                REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

                //first seed 
                REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(2).get_value() == 5);

                //The next seed
                REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);

                //Chain end
                REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::CHAIN_END);
            }

            SECTION( "Count dags" ) {
                for (auto& zip_tree : zip_forest.trees) {
                    pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                    REQUIRE(dag_non_dag_count.first == 0);
                    REQUIRE(dag_non_dag_count.second == 0);
                }
            }
            //TODO: This fails now that it is a forest

            // For each seed, what seeds and distances do we see in reverse from it?
            //std::unordered_map<ZipCodeTree::oriented_seed_t, std::vector<ZipCodeTree::seed_result_t>> reverse_views;
            //for (auto forward = zip_tree.begin(); forward != zip_tree.end(); ++forward) {
            //    std::copy(zip_tree.look_back(forward), zip_tree.rend(), std::back_inserter(reverse_views[*forward]));
            //}
            //REQUIRE(reverse_views.size() == 4);
            //// The first seed can't see any other seeds
            //REQUIRE(reverse_views.count({0, false}));
            //REQUIRE(reverse_views[{0, false}].size() == 0);
            //// The second seed can see the first seed at distance 5
            //REQUIRE(reverse_views.count({1, false}));
            //REQUIRE(reverse_views[{1, false}].size() == 1);
            //REQUIRE(reverse_views[{1, false}][0].seed == 0);
            //REQUIRE(reverse_views[{1, false}][0].distance == 5);
            //REQUIRE(reverse_views[{1, false}][0].is_reverse == false);
            //// The third seed can't see any other seeds
            //REQUIRE(reverse_views.count({2, false}));
            //REQUIRE(reverse_views[{2, false}].size() == 0);
            //// The fourth seed can see the third seed at distance 5
            //REQUIRE(reverse_views.count({3, false}));
            //REQUIRE(reverse_views[{3, false}].size() == 1);
            //REQUIRE(reverse_views[{3, false}][0].seed == 2);
            //REQUIRE(reverse_views[{3, false}][0].distance == 5);
            //REQUIRE(reverse_views[{3, false}][0].is_reverse == false);
        }
        SECTION( "Four buckets" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 5);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 5);

            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            REQUIRE(zip_forest.trees.size() == 4);

            zip_forest.print_self(&seeds, &minimizers);
        }
    }
    TEST_CASE( "zip tree simple bubbles in chains", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);

        //graph.to_dot(cerr);

        SECTION( "Seeds on chain nodes" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            //The tree should be:
            // [pos1 3 pos3 6 pos6]
            //or backwards
            REQUIRE(zip_tree.get_tree_size() == 7);

            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

            //first seed 
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(1).get_is_reversed()) {
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 2);
            } else {
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 0);
            }

            //distance between them
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(2).get_value() == 3 ||
                    zip_tree.get_item_at_index(2).get_value() == 6));

            //the next seed
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);

            //distance between them
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(4).get_value() == 3 ||
                    zip_tree.get_item_at_index(4).get_value() == 6));

            //the last seed
            REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(5).get_is_reversed()) {
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 0);
            } else {
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 2);
            }

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);
            
            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 0);
            }

            // TODO: This time we happen to visit the seeds in reverse order.
            // How are we doing querying in a particular direction relative to a particular seed?

            // We see all the seeds in order
            std::vector<ZipCodeTree::oriented_seed_t> seed_indexes;
            std::copy(zip_tree.begin(), zip_tree.end(), std::back_inserter(seed_indexes));
            REQUIRE(seed_indexes.size() == 3);
            if (seed_indexes.at(0).is_reverse) {
                REQUIRE(seed_indexes.at(0).seed == 2);
                REQUIRE(seed_indexes.at(1).seed == 1);
                REQUIRE(seed_indexes.at(2).seed == 0);
            } else {
                REQUIRE(seed_indexes.at(0).seed == 0);
                REQUIRE(seed_indexes.at(1).seed == 1);
                REQUIRE(seed_indexes.at(2).seed == 2);    
            }

            // For each seed, what seeds and distances do we see in reverse from it?
            std::unordered_map<ZipCodeTree::oriented_seed_t, std::vector<ZipCodeTree::seed_result_t>> reverse_views;
            for (auto forward = zip_tree.begin(); forward != zip_tree.end(); ++forward) {
                std::copy(zip_tree.look_back(forward), zip_tree.rend(), std::back_inserter(reverse_views[*forward]));
            }
            REQUIRE(reverse_views.size() == 3);
            if (seed_indexes.at(0).is_reverse) {
                // The first seed can't see any other seeds
                REQUIRE(reverse_views.count({2, true}));
                REQUIRE(reverse_views[{2, true}].size() == 0);
                // The second seed can see the first seed at distance 6
                REQUIRE(reverse_views.count({1, true}));
                REQUIRE(reverse_views[{1, true}].size() == 1);
                REQUIRE(reverse_views[{1, true}][0].seed == 2);
                REQUIRE(reverse_views[{1, true}][0].distance == 6);
                REQUIRE(reverse_views[{1, true}][0].is_reverse == true);
                // The third seed can't see both the others at distances 3 and 9
                REQUIRE(reverse_views.count({0, true}));
                REQUIRE(reverse_views[{0, true}].size() == 2);
                REQUIRE(reverse_views[{0, true}][0].seed == 1);
                REQUIRE(reverse_views[{0, true}][0].distance == 3);
                REQUIRE(reverse_views[{0, true}][0].is_reverse == true);
                REQUIRE(reverse_views[{0, true}][1].seed == 2);
                REQUIRE(reverse_views[{0, true}][1].distance == 9);
                REQUIRE(reverse_views[{0, true}][1].is_reverse == true);
            } else {
                // The first seed can't see any other seeds
                REQUIRE(reverse_views.count({0, false}));
                REQUIRE(reverse_views[{0, false}].size() == 0);
                // The second seed can see the first seed at distance 3
                REQUIRE(reverse_views.count({1, false}));
                REQUIRE(reverse_views[{1, false}].size() == 1);
                REQUIRE(reverse_views[{1, false}][0].seed == 0);
                REQUIRE(reverse_views[{1, false}][0].distance == 3);
                REQUIRE(reverse_views[{1, false}][0].is_reverse == false);
                // The third seed can't see both the others at distances 6 and 9
                REQUIRE(reverse_views.count({2, false}));
                REQUIRE(reverse_views[{2, false}].size() == 2);
                REQUIRE(reverse_views[{2, false}][0].seed == 1);
                REQUIRE(reverse_views[{2, false}][0].distance == 6);
                REQUIRE(reverse_views[{2, false}][0].is_reverse == false);
                REQUIRE(reverse_views[{2, false}][1].seed == 2);
                REQUIRE(reverse_views[{2, false}][1].distance == 9);
                REQUIRE(reverse_views[{2, false}][1].is_reverse == false);
            }
        }
        SECTION( "Seeds on chain nodes one reversed" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, true, 2);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            //The tree should be:
            // [pos1 3 pos3 6 pos6]
            //or backwards
            REQUIRE(zip_tree.get_tree_size() == 7);

            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).get_type() == ZipCodeTree::CHAIN_START);

            //first seed 
            //This is either the first seed on 1 going backwards, or the third seed on 6 going backwards
            REQUIRE(zip_tree.get_item_at_index(1).get_type() == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(1).get_value() == 0) {
                REQUIRE(zip_tree.get_item_at_index(1).get_is_reversed());
            } else {
                REQUIRE(zip_tree.get_item_at_index(1).get_value() == 2);
                REQUIRE(zip_tree.get_item_at_index(1).get_is_reversed());
            }

            //distance between them
            REQUIRE(zip_tree.get_item_at_index(2).get_type() == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(2).get_value() == 2 ||
                    zip_tree.get_item_at_index(2).get_value() == 6));

            //the next seed
            REQUIRE(zip_tree.get_item_at_index(3).get_type() == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(3).get_value() == 1);

            //distance between them
            REQUIRE(zip_tree.get_item_at_index(4).get_type() == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(4).get_value() == 2 ||
                    zip_tree.get_item_at_index(4).get_value() == 6));

            //the last seed
            REQUIRE(zip_tree.get_item_at_index(5).get_type() == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(5).get_value() == 0) {
                REQUIRE(!zip_tree.get_item_at_index(5).get_is_reversed());
            } else {
                REQUIRE(zip_tree.get_item_at_index(5).get_value() == 2);
                REQUIRE(!zip_tree.get_item_at_index(5).get_is_reversed());
            }

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(6).get_type() == ZipCodeTree::CHAIN_END);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION( "One seed on snarl" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 1);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            //The tree should be:
            // [pos1 3 ( 2 [ pos2 ] 6 0 1 ) 0  pos3 6 pos6]
            //or backwards
            REQUIRE(zip_tree.get_tree_size() == 17);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION( "Three seeds on snarl" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 1);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(2, false, 4);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            //The tree should be:
            // [pos1 0 ( 0 [ pos2 x pos2 x pos2 ] 0 0 1 ) 0  pos3 6 pos6]
            //or backwards
            REQUIRE(zip_tree.get_tree_size() == 21);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION( "Two children of a snarl" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(5, false, 1);
            positions.emplace_back(6, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            //The tree should be:
            // [pos1 0  pos3 0 ( 0 [ pos4 ] inf 0 [ pos5 1 pos5 ] 2 3 3 2) 0 pos6]
            //or backwards
            REQUIRE(zip_tree.get_tree_size() == 25);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION( "Only snarls in a chain" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(5, false, 1);

            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            //The tree should be:
            // [( 0 [ pos2 ] 7 0 1) 3 ( 0 [pos4 ] 3 inf [pos5 1 pos5 ] 2 0 3 2 )]
            //or backwards
            REQUIRE(zip_tree.get_tree_size() == 29);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 2);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION( "Seeds on chain nodes bucket" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 4);
            REQUIRE(zip_forest.trees.size() == 2);
            zip_forest.print_self(&seeds, &minimizers);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
        SECTION( "Only snarls in two buckets" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 1);

            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 2);
            REQUIRE(zip_forest.trees.size() == 2);
            zip_forest.print_self(&seeds, &minimizers);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
        SECTION( "Snarls and nodes in three buckets" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 1);

            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 1);
            REQUIRE(zip_forest.trees.size() == 3);
            zip_forest.print_self(&seeds, &minimizers);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
        SECTION( "Chain in snarl in a separate bucket" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 3);
            positions.emplace_back(2, false, 3);
            positions.emplace_back(3, false, 0);

            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 2);
            REQUIRE(zip_forest.trees.size() == 2);
            zip_forest.print_self(&seeds, &minimizers);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
        SECTION( "Chain in snarl in a separate bucket another connected to end (or maybe start)" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 3);
            positions.emplace_back(3, false, 0);

            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 2);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 2);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
    }
    TEST_CASE( "zip tree simple nested bubbles in chains", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);

        //graph.to_dot(cerr);

        SECTION( "Slice of snarl removed" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 6);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);

            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 4);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 2);
            zip_forest.validate_zip_forest(distance_index, &seeds, 4);
        }
    }
    TEST_CASE( "zip tree bubble in cyclic snarl", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);

        //graph.to_dot(cerr);

        SECTION( "Two sides of nested snp unordered along read" ) {
 
            vector<std::pair<pos_t, size_t>> positions;
            positions.emplace_back(make_pos_t(1, false, 0), 0);
            positions.emplace_back(make_pos_t(5, false, 5), 1);
            positions.emplace_back(make_pos_t(4, false, 0), 2);
            positions.emplace_back(make_pos_t(5, false, 5), 3);
            positions.emplace_back(make_pos_t(3, false, 0), 4);

            vector<MinimizerMapper::Minimizer> minimizers;
            

            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (auto pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos.first);
                seeds.push_back({ pos.first, pos.second, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = pos.second;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max(), 4);
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_forest.validate_zip_forest(distance_index, &seeds, 4);
        }
    }
    TEST_CASE( "zip tree snarl with inversion", "[zip_tree]" ) {

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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);


      
        //graph.to_dot(cerr);

        SECTION( "Traverse 3 backwards" ) {
 
            vector<std::pair<pos_t, size_t>> positions;
            positions.emplace_back(make_pos_t(1, false, 0), 0);
            positions.emplace_back(make_pos_t(4, false, 0), 1);
            positions.emplace_back(make_pos_t(3, true, 0), 2);
            positions.emplace_back(make_pos_t(3, true, 1), 3);
            positions.emplace_back(make_pos_t(5, false, 0), 4);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (auto pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos.first);
                seeds.push_back({ pos.first, pos.second, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = pos.second;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);


            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_forest.validate_zip_forest(distance_index, &seeds);

            bool chain_is_reversed = distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id()));
            if (chain_is_reversed) {
                cerr << "This test didn't get run because I'm lazy and didn't write it for a reversed chain" << endl;
                
            } else {
                //For a forward traversal of the chain, the zip tree should be:
                //[1+0/0 3 ( 0 [4+0/1] 18446744073709551615  12 [4+0/1rev] 18446744073709551615  2  2 [3-0/2 1 3-1/3] 5  18446744073709551615  8  8  3) 0 5+0/4]

                //Check some random elements

                //First seed
                REQUIRE(zip_forest.trees[0].get_item_at_index(1).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(1).get_value() == 0);
                //Chain start
                REQUIRE(zip_forest.trees[0].get_item_at_index(5).get_type() == ZipCodeTree::CHAIN_START);
                //Second seed (4)
                REQUIRE(zip_forest.trees[0].get_item_at_index(6).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(6).get_value() == 1);

                //Third seed (4 in the other direction
                REQUIRE(zip_forest.trees[0].get_item_at_index(11).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(6).get_value() == 1);

                //Fourth seed (3-1
                REQUIRE(zip_forest.trees[0].get_item_at_index(17).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(17).get_value() == 2);

            }

        }

    }
    TEST_CASE( "zip tree non-simple DAG", "[zip_tree]" ) {

        //bubble between 1 and 3, non-simple dag between 3 and 8 
        //containing node 7 and chain 4-6
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
      
        //graph.to_dot(cerr);

        SECTION( "Make the zip tree" ) {
 
            vector<std::pair<pos_t, size_t>> positions;
            positions.emplace_back(make_pos_t(1, false, 0), 0);
            positions.emplace_back(make_pos_t(2, false, 0), 1);
            positions.emplace_back(make_pos_t(3, false, 0), 2);
            positions.emplace_back(make_pos_t(3, false, 1), 3);
            positions.emplace_back(make_pos_t(4, false, 0), 4);
            positions.emplace_back(make_pos_t(5, false, 0), 5);
            positions.emplace_back(make_pos_t(6, false, 0), 6);
            positions.emplace_back(make_pos_t(7, false, 1), 7);
            positions.emplace_back(make_pos_t(8, false, 0), 8);
            positions.emplace_back(make_pos_t(8, false, 2), 9);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (auto pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos.first);
                seeds.push_back({ pos.first, pos.second, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = pos.second;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            bool chain_is_reversed = distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id()));
            if (chain_is_reversed) {
                
            } else {
                //For a forward traversal of the chain, the zip tree should be:
                //[1+0/0 3 ( 0 [2+0/0] 4  0  1) 0 3+0/0 1 3+1/0 5 ( 0 [4+0/0 2 ( 0 [5+0/0] 2  0  1) 0 6+0/0] 4  1 [7+1/0] 2  6  0  2) 0 8+0/0 2 8+2/0]
                //Check some random elements

                //First seed
                REQUIRE(zip_forest.trees[0].get_item_at_index(1).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(1).get_value() == 0);
                //Start of cyclic snarl
                REQUIRE(zip_forest.trees[0].get_item_at_index(17).get_type() == ZipCodeTree::SNARL_START);
                REQUIRE(zip_forest.trees[0].get_item_at_index(25).get_type() == ZipCodeTree::SEED);
                REQUIRE(zip_forest.trees[0].get_item_at_index(25).get_value() == 5);

                REQUIRE(zip_forest.trees[0].get_item_at_index(30).get_type() == ZipCodeTree::SNARL_END);

                REQUIRE(zip_forest.trees[0].get_item_at_index(34).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_forest.trees[0].get_item_at_index(34).get_value() == 4);
                REQUIRE(zip_forest.trees[0].get_item_at_index(35).get_type() == ZipCodeTree::EDGE);
                REQUIRE(zip_forest.trees[0].get_item_at_index(35).get_value() == 1);

            }

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 3);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION( "Three buckets" ) {
 
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
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            REQUIRE(zip_forest.trees.size() == 3);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);
        }


    }

    TEST_CASE( "zip tree deeply nested bubbles", "[zip_tree]" ) {
        //top-level chain 1-12-13-16
        //bubble 2-10 containing two bubbles 3-5 and 6-9
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        
        //graph.to_dot(cerr);

        SECTION( "Make the zip tree with a seed on each node" ) {
 
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
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 5);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION( "Make the zip tree with a few seeds" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);
            positions.emplace_back(13, false, 2);
            positions.emplace_back(15, false, 2);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizers);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 3);
                REQUIRE(dag_non_dag_count.second == 0);
            }
        }
        SECTION( "3 buckets" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(10, false, 0);
            positions.emplace_back(13, false, 2);
            positions.emplace_back(16, false, 5);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 4);
            REQUIRE(zip_forest.trees.size() == 3);
            zip_forest.print_self(&seeds, &minimizers);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
        SECTION( "Remove empty snarls" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(6, false, 1);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(4, false, 1);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 2);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 3);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
        SECTION( "Chain connected on one end" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(2, false, 2);
            positions.emplace_back(6, false, 1);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(4, false, 1);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 2);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 2);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
        SECTION( "Chain connected on the other end" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(10, false, 0);
            positions.emplace_back(10, false, 2);
            positions.emplace_back(9, false, 1);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(4, false, 1);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 2);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 2);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
        SECTION( "One chain removed from a snarl" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(8, false, 1);
            positions.emplace_back(7, false, 1);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(11, false, 1);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 2);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 3);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }
        }
    }
    TEST_CASE( "zip tree long nested chain", "[zip_tree]" ) {
        //top-level chain 1-12-13-16
        //bubble 2-10 containing two bubbles 3-5 and 6-9
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        cerr << distance_index.net_handle_as_string(distance_index.get_parent(distance_index.get_node_net_handle(n1->id()))) << endl;
        
        //graph.to_dot(cerr);

        SECTION( "One slice from nodes in the middle of a nested chain" ) {
 
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


            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 2);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }

        }
        SECTION( "Two slices from snarls in the middle of a nested chain" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(6, false, 1);
            positions.emplace_back(7, false, 0);
            positions.emplace_back(11, false, 0);
            positions.emplace_back(12, false, 0);
            positions.emplace_back(21, false, 0);  


            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 2);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 4);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }

        }
        SECTION( "One slice from the start of a chain, connected to the end" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(7, false, 0);
            positions.emplace_back(12, false, 1);
            positions.emplace_back(13, false, 0);
            positions.emplace_back(21, false, 0);    


            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 2);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }

        }
        SECTION( "One slice from the end of a chain, connected to the start" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 2);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(7, false, 0);
            positions.emplace_back(14, false, 0);    
            positions.emplace_back(16, false, 0);    
            positions.emplace_back(20, false, 0);    
            positions.emplace_back(21, false, 0);    


            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            REQUIRE(zip_forest.trees.size() == 2);
            for (auto& zip_tree : zip_forest.trees) {
                zip_tree.validate_zip_tree(distance_index, &seeds);
            }

        }
    }

    TEST_CASE( "zip tree non-dag", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
       
        //graph.to_dot(cerr);

        SECTION( "Make the zip tree with a seed on each node" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = 0;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 1);
            }
        }

    }
    TEST_CASE( "zip tree nested cyclic non-dag", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        
        ofstream out ("testGraph.hg");
        graph.serialize(out);


        //graph.to_dot(cerr);

        SECTION( "Make the zip tree with a seed on each node" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = 0;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 2);
            }
        }

    }
    TEST_CASE( "zip tree nested inversions", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);

        //graph.to_dot(cerr);

        SECTION( "Go forward through the inversions" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(2, false, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, false, 1);
            positions.emplace_back(4, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);


            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (size_t i = 0 ; i < positions.size() ; i++) {
                pos_t pos = positions[i];
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, i, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = i;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            assert(zip_tree.get_tree_size() == 31);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 2);
            }
        }
        SECTION( "Reverse both inversions" ) {
 
            vector<pos_t> positions;
            positions.emplace_back(1, false, 0);
            positions.emplace_back(4, true, 0);
            positions.emplace_back(3, false, 0);
            positions.emplace_back(3, false, 1);
            positions.emplace_back(2, true, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 0);


            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (size_t i = 0 ; i < positions.size() ; i++) {
                pos_t pos = positions[i];
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, i, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = i;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            ZipCodeTree zip_tree = zip_forest.trees[0];
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_tree.validate_zip_tree(distance_index, &seeds);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 0);
                REQUIRE(dag_non_dag_count.second == 2);
            }
        }
    }
    TEST_CASE( "zip tree cyclic snarl with overlapping seeds", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        

        //graph.to_dot(cerr);

        SECTION( "Cyclic snarl with seeds on either side" ) {
 
            vector<pair<pos_t, size_t>> positions;
            positions.emplace_back(make_pos_t(1, false, 0), 0);
            positions.emplace_back(make_pos_t(2, false, 0), 1);
            positions.emplace_back(make_pos_t(2, false, 2), 2);
            positions.emplace_back(make_pos_t(2, false, 4), 3);
            positions.emplace_back(make_pos_t(2, false, 0), 4);
            positions.emplace_back(make_pos_t(2, false, 2), 5);
            positions.emplace_back(make_pos_t(2, false, 4), 6);

            positions.emplace_back(make_pos_t(3, false, 0), 6);
            positions.emplace_back(make_pos_t(3, false, 2), 5);
            positions.emplace_back(make_pos_t(3, false, 4), 4);
            positions.emplace_back(make_pos_t(3, false, 0), 3);
            positions.emplace_back(make_pos_t(3, false, 2), 2);
            positions.emplace_back(make_pos_t(3, false, 4), 1);

            positions.emplace_back(make_pos_t(4, false, 0), 1);
            positions.emplace_back(make_pos_t(4, false, 2), 2);
            positions.emplace_back(make_pos_t(4, false, 4), 3);
            positions.emplace_back(make_pos_t(4, false, 0), 4);
            positions.emplace_back(make_pos_t(4, false, 2), 5);
            positions.emplace_back(make_pos_t(4, false, 4), 6);
            positions.emplace_back(make_pos_t(5, false, 4), 7);
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (auto pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos.first);
                seeds.push_back({ pos.first, pos.second, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = pos.second;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_forest.validate_zip_forest(distance_index, &seeds);

        }
    }
    TEST_CASE( "zip tree duplication", "[zip_tree]" ) {
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
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);
        

        //graph.to_dot(cerr);

        SECTION( "Cyclic snarl with seeds on either side" ) {
 
            vector<pair<pos_t, size_t>> positions;
            positions.emplace_back(make_pos_t(1, false, 0), 0);
            positions.emplace_back(make_pos_t(2, false, 0), 1);
            positions.emplace_back(make_pos_t(2, false, 1), 2);
            positions.emplace_back(make_pos_t(2, false, 2), 3);
            positions.emplace_back(make_pos_t(2, false, 0), 4);
            positions.emplace_back(make_pos_t(2, false, 1), 5);
            positions.emplace_back(make_pos_t(2, false, 2), 6);
            positions.emplace_back(make_pos_t(3, false, 0), 7);

            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (auto pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos.first);
                seeds.push_back({ pos.first, pos.second, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = pos.second;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            REQUIRE(zip_forest.trees.size() == 1);
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_forest.validate_zip_forest(distance_index, &seeds);

        }
    }

    TEST_CASE("zip tree handles complicated nested snarls", "[zip_tree]" ) {
        
        // Load an example graph
        VG graph;
        io::json2graph(R"({"node":[{"id": "1","sequence":"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"},{"id":"2","sequence":"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"},{"id":"63004428","sequence":"T"},{"id":"63004425","sequence":"T"},{"id":"63004426","sequence":"ATATCTATACATATAATACAG"},{"id":"63004421","sequence":"AT"},{"id":"63004422","sequence":"T"},{"id":"63004424","sequence":"A"},{"id":"63004429","sequence":"C"},{"id":"63004430","sequence":"AT"},{"id":"63004427","sequence":"A"},{"id":"63004423","sequence":"C"}],"edge":[{"from":"63004428","to":"63004430"},{"from":"63004425","to":"63004426"},{"from":"63004426","to":"63004427"},{"from":"63004421","to":"63004422"},{"from":"63004422","to":"63004427"},{"from":"63004422","to":"63004423","to_end":true},{"from":"63004422","to":"63004424"},{"from":"63004424","to":"63004425"},{"from":"63004429","to":"63004430"},{"from":"63004427","to":"63004428"},{"from":"63004427","to":"63004429"},{"from":"63004423","from_start":true,"to":"63004428"},{"from":"1","to":"63004421"},{"from":"63004430","to":"2"}]})", &graph);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);
        SnarlDistanceIndexClusterer clusterer(distance_index, &graph);

  

        // I observed:
        // 63004421+0 2 ( 4 [63004426+1] 19  2  1) 2 63004430+1 
        // But we want 63004426+1 to 63004430+1 to be 23 and not 21.

        vector<pos_t> positions;
        positions.emplace_back(63004421, false, 0);
        positions.emplace_back(63004426, false, 1);
        positions.emplace_back(63004430, false, 1);
        
        vector<SnarlDistanceIndexClusterer::Seed> seeds;
        for (pos_t pos : positions) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, pos);
            seeds.push_back({ pos, 0, zipcode});
        }

        VectorView<MinimizerMapper::Minimizer> minimizers;

        ZipCodeForest zip_forest;
        zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
        REQUIRE(zip_forest.trees.size() == 1);
        ZipCodeTree zip_tree = zip_forest.trees[0];
        zip_forest.print_self(&seeds, &minimizers);
        zip_tree.validate_zip_tree(distance_index, &seeds);
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


        vector<pos_t> positions;
        positions.emplace_back(1, false, 0);
        positions.emplace_back(2, false, 0);
        positions.emplace_back(3, true, 0);
        positions.emplace_back(4, false, 0);
        
        vector<SnarlDistanceIndexClusterer::Seed> seeds;
        for (pos_t pos : positions) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, pos);
            seeds.push_back({ pos, 0, zipcode});
        }

        VectorView<MinimizerMapper::Minimizer> minimizers;

        ZipCodeForest zip_forest;
        zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
        REQUIRE(zip_forest.trees.size() == 1);
        ZipCodeTree zip_tree = zip_forest.trees[0];
        zip_forest.print_self(&seeds, &minimizers);
        //TODO: This doesn't actually have the right distances yet, I just want to make sure it won't crash
        //zip_tree.validate_zip_tree(distance_index, &seeds);
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

        vector<pos_t> positions;
        positions.emplace_back(5, false, 0);
        positions.emplace_back(7, false, 17);
        
        vector<SnarlDistanceIndexClusterer::Seed> seeds;
        for (pos_t pos : positions) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, pos);
            seeds.push_back({ pos, 0, zipcode});
        }

        VectorView<MinimizerMapper::Minimizer> minimizers;

        ZipCodeForest zip_forest;
        zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 61);
        zip_forest.print_self(&seeds, &minimizers);
        zip_forest.validate_zip_forest(distance_index, &seeds, 61);
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
        
        vector<SnarlDistanceIndexClusterer::Seed> seeds;
        vector<MinimizerMapper::Minimizer> minimizers;
        for (pos_t pos : positions) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, pos);
            seeds.push_back({ pos, 0, zipcode});

            minimizers.emplace_back();
            minimizers.back().value.offset = 0;
            minimizers.back().value.is_reverse = false;
        }

        VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);



        ZipCodeForest zip_forest;
        zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max(), 5);
        zip_forest.print_self(&seeds, &minimizer_vector);
        REQUIRE(zip_forest.trees.size() == 6);
        for (auto& tree : zip_forest.trees) {
            tree.validate_zip_tree(distance_index, &seeds);
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

        SECTION( "Multiple seeds in snarl" ) {
            vector<std::pair<pos_t, size_t>> positions;
            positions.emplace_back(make_pos_t(2, false, 0), 0);
            positions.emplace_back(make_pos_t(3, false, 0), 1);
            positions.emplace_back(make_pos_t(3, true,  0), 2);
            positions.emplace_back(make_pos_t(5, true,  0), 3);
            positions.emplace_back(make_pos_t(6, true,  0), 4);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;
            for (size_t i = 0 ; i < positions.size(); ++i) {
                auto pos = positions[i];
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos.first);
                seeds.push_back({ pos.first, i, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = pos.second;
                minimizers.back().value.is_reverse = false;
            }

            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, std::numeric_limits<size_t>::max());
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_forest.validate_zip_forest(distance_index, &seeds);
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

        SECTION( "Node first" ) {
            vector<pos_t> positions;
            positions.emplace_back(2, false, 0);
            positions.emplace_back(5, false, 0);
            positions.emplace_back(6, false, 4);
            positions.emplace_back(10, false, 0);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            zip_forest.validate_zip_forest(distance_index, &seeds, 3);
        }
        SECTION( "Snarl first" ) {
            vector<pos_t> positions;
            positions.emplace_back(3, false, 0);
            positions.emplace_back(6, false, 4);
            positions.emplace_back(10, false, 0);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            zip_forest.validate_zip_forest(distance_index, &seeds, 3);
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

        SECTION( "One tree on each node" ) {
            vector<pos_t> positions;
            positions.emplace_back(2, false, 7);
            positions.emplace_back(3, false, 3);
            positions.emplace_back(4, false, 7);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            zip_forest.validate_zip_forest(distance_index, &seeds, 3);
        }
        SECTION( "Remove second child of snarl" ) {
            vector<pos_t> positions;
            positions.emplace_back(3, false, 8);
            positions.emplace_back(4, false, 5);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            zip_forest.validate_zip_forest(distance_index, &seeds, 3);
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
       

        //ofstream out ("testGraph.hg");
        //graph.serialize(out);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION( "One tree on each node" ) {
            vector<pos_t> positions;
            positions.emplace_back(1, false, 5);
            positions.emplace_back(2, false, 5);
            positions.emplace_back(3, false, 5);
            positions.emplace_back(4, false, 5);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            zip_forest.validate_zip_forest(distance_index, &seeds, 3);
        }
        SECTION( "Remove second child of snarl" ) {
            vector<pos_t> positions;
            positions.emplace_back(3, false, 8);
            positions.emplace_back(4, false, 5);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            zip_forest.validate_zip_forest(distance_index, &seeds, 3);
        }
        SECTION( "Remove first child of snarl" ) {
            vector<pos_t> positions;
            positions.emplace_back(3, false, 5);
            positions.emplace_back(4, false, 0);
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            VectorView<MinimizerMapper::Minimizer> minimizers;

            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max(), 3);
            zip_forest.print_self(&seeds, &minimizers);
            zip_forest.validate_zip_forest(distance_index, &seeds, 3);
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

        SECTION( "Snps alone" ) {
            vector<pair<pos_t, size_t>> positions;
            positions.emplace_back(make_pos_t(1, false, 0), 1);
            positions.emplace_back(make_pos_t(2, false, 8), 2);
            positions.emplace_back(make_pos_t(3, false, 0), 3);
            positions.emplace_back(make_pos_t(5, false, 0), 4);
            positions.emplace_back(make_pos_t(2, false, 8), 15);
            positions.emplace_back(make_pos_t(4, false, 5), 16);
            positions.emplace_back(make_pos_t(5, false, 0), 17);
            positions.emplace_back(make_pos_t(7, false, 0), 18);

        distance_index.for_each_child(distance_index.get_root(), [&](net_handle_t child) {
            cerr << distance_index.net_handle_as_string(child) << endl;
        });
            
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            vector<MinimizerMapper::Minimizer> minimizers;

            for (size_t i = 0 ; i < positions.size() ; ++i) {
                auto pos = positions[i];
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos.first);
                seeds.push_back({ pos.first, i, zipcode});

                minimizers.emplace_back();
                minimizers.back().value.offset = pos.second;
                minimizers.back().value.is_reverse = false;
            }
            VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);


            ZipCodeForest zip_forest;
            zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, 100, 100);
            zip_forest.print_self(&seeds, &minimizer_vector);
            zip_forest.validate_zip_forest(distance_index, &seeds, 100);
        }


    }
    /*

    TEST_CASE("Failed unit test", "[failed]") {
        //Load failed random graph
        HashGraph graph;
        graph.deserialize("testGraph.hg");

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        vector<std::pair<pos_t, size_t>> positions;
        positions.emplace_back(make_pos_t(20, false, 7), 0);
        positions.emplace_back(make_pos_t(23, false, 0), 3);
        positions.emplace_back(make_pos_t(13, true, 3), 1);
        positions.emplace_back(make_pos_t(18, false, 0), 8);
        positions.emplace_back(make_pos_t(17, true, 0), 5);
        positions.emplace_back(make_pos_t(19, false, 1), 14);
        positions.emplace_back(make_pos_t(33, false, 0), 15);
        positions.emplace_back(make_pos_t(11, false, 0), 2);
        positions.emplace_back(make_pos_t(10, false, 3), 16);



        vector<SnarlDistanceIndexClusterer::Seed> seeds;
        for (auto pos : positions) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, pos.first);
            seeds.push_back({ pos.first, pos.second, zipcode});
        }
        distance_index.for_each_child(distance_index.get_root(), [&](net_handle_t child) {
            cerr << distance_index.net_handle_as_string(child) << endl;
        });
        VectorView<MinimizerMapper::Minimizer> minimizers;

        ZipCodeForest zip_forest;
        zip_forest.fill_in_forest(seeds, minimizers, distance_index, std::numeric_limits<size_t>::max());
        zip_forest.print_self(&seeds, &minimizers);
        zip_forest.validate_zip_forest(distance_index, &seeds);
    }
    */



    TEST_CASE("Random graphs zip tree", "[zip_tree][zip_tree_random]"){
    
        for (int i = 0; i < 0; i++) {
            // For each random graph
    
            default_random_engine generator(time(NULL));
            uniform_int_distribution<int> variant_count(1, 10);
            uniform_int_distribution<int> chrom_len(10, 200);
            uniform_int_distribution<int> distance_limit(5, 100);
    
            //Make a random graph with three chromosomes of random lengths
            HashGraph graph;
            random_graph({chrom_len(generator),chrom_len(generator),chrom_len(generator)}, 30, variant_count(generator), &graph);
            graph.serialize("testGraph.hg");

            IntegratedSnarlFinder snarl_finder(graph);
            SnarlDistanceIndex distance_index;
            fill_in_distance_index(&distance_index, &graph, &snarl_finder);
            SnarlDistanceIndexClusterer clusterer(distance_index, &graph);

            vector<id_t> all_nodes;
            graph.for_each_handle([&](const handle_t& h)->bool{
                id_t id = graph.get_id(h);
                all_nodes.push_back(id);
                return true;
            });

            uniform_int_distribution<int> randPosIndex(0, all_nodes.size()-1);

            //Check k random sets of seeds
            for (size_t k = 0; k < 10 ; k++) {

                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                vector<MinimizerMapper::Minimizer> minimizers;

                uniform_int_distribution<int> randPosCount(3, 70);
                for (int j = 0; j < randPosCount(generator); j++) {
                    //Check clusters of j random positions

                    id_t nodeID1 = all_nodes[randPosIndex(generator)];
                    handle_t node1 = graph.get_handle(nodeID1);

                    offset_t offset1 = uniform_int_distribution<int>(0,graph.get_length(node1) - 1)(generator);

                    pos_t pos = make_pos_t(nodeID1,
                                           uniform_int_distribution<int>(0,1)(generator) == 0,
                                           offset1 );

                    ZipCode zipcode;
                    zipcode.fill_in_zipcode(distance_index, pos);

                    seeds.push_back({ pos, (size_t)j, zipcode});

                    minimizers.emplace_back();
                    minimizers.back().value.offset = (size_t) j;
                    minimizers.back().value.is_reverse = false;

                }
                size_t limit = distance_limit(generator);

                VectorView<MinimizerMapper::Minimizer> minimizer_vector(minimizers);

                ZipCodeForest zip_forest;
                zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, limit, limit);
                zip_forest.print_self(&seeds, &minimizer_vector);
                zip_forest.validate_zip_forest(distance_index, &seeds, limit);
                REQUIRE(true); //Just to count
            }
        }
    }

}
}
