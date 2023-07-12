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
#include "../zip_code_tree.hpp"
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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            REQUIRE(zip_tree.get_tree_size() == 3);
            REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);
            REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(1).value == 0);
            REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::CHAIN_END);

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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            REQUIRE(zip_tree.get_tree_size() == 5);


            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);

            //Seed (either one because they're the same position)
            REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(1).value == 0 ||
                     zip_tree.get_item_at_index(1).value == 1));

            //Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(2).value == 0);

            //THe other seed
            REQUIRE(zip_tree.get_item_at_index(3).type == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(3).value == 0 ||
                     zip_tree.get_item_at_index(3).value == 1));

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(4).type == ZipCodeTree::CHAIN_END);

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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            REQUIRE(zip_tree.get_tree_size() == 7);


            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);

            //Seed (either one because they're the same position)
            REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(1).value == 0 ||
                     zip_tree.get_item_at_index(1).value == 1));

            //Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(2).value == 0);

            //THe other seed
            REQUIRE(zip_tree.get_item_at_index(3).type == ZipCodeTree::SEED);
            REQUIRE((zip_tree.get_item_at_index(3).value == 0 ||
                     zip_tree.get_item_at_index(3).value == 1));

            //Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(4).type == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(4).value == 2);

            //The other seed
            REQUIRE(zip_tree.get_item_at_index(5).type == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(5).value == 2);

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(6).type == ZipCodeTree::CHAIN_END);

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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            REQUIRE(zip_tree.get_tree_size() == 7);

            //The order should either be 0-1-2, or 2-1-0
            bool is_rev = zip_tree.get_item_at_index(1).value == 2;
            if (is_rev) {

                //Chain start
                REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);

                //first seed 
                REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(1).value == 2);
                REQUIRE(zip_tree.get_item_at_index(1).is_reversed == true);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(2).value == 4);

                //The next seed
                REQUIRE(zip_tree.get_item_at_index(3).type == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(3).value == 1);
                REQUIRE(zip_tree.get_item_at_index(3).is_reversed == true);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(4).type == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(4).value == 1);

                //The last seed
                REQUIRE(zip_tree.get_item_at_index(5).type == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(5).value == 0);
                REQUIRE(zip_tree.get_item_at_index(5).is_reversed == true);

                //Chain end
                REQUIRE(zip_tree.get_item_at_index(6).type == ZipCodeTree::CHAIN_END);
            } else {

                //Chain start
                REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);

                //first seed 
                REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(1).value == 0);
                REQUIRE(zip_tree.get_item_at_index(1).is_reversed == false);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(2).value == 1);

                //The next seed
                REQUIRE(zip_tree.get_item_at_index(3).type == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(3).value == 1);
                REQUIRE(zip_tree.get_item_at_index(3).is_reversed == false);

                //Distance between the seeds
                REQUIRE(zip_tree.get_item_at_index(4).type == ZipCodeTree::EDGE);
                REQUIRE(zip_tree.get_item_at_index(4).value == 4);

                //The last seed
                REQUIRE(zip_tree.get_item_at_index(5).type == ZipCodeTree::SEED);
                REQUIRE(zip_tree.get_item_at_index(5).value == 2);
                REQUIRE(zip_tree.get_item_at_index(5).is_reversed == false);

                //Chain end
                REQUIRE(zip_tree.get_item_at_index(6).type == ZipCodeTree::CHAIN_END);
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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            //The tree should be:
            // [pos1] [pos3]
            REQUIRE(zip_tree.get_tree_size() == 6);

            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);

            //first seed 
            REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::CHAIN_END);

            //Chain start
            REQUIRE(zip_tree.get_item_at_index(3).type == ZipCodeTree::CHAIN_START);

            //The first seed in the new chain
            REQUIRE(zip_tree.get_item_at_index(4).type == ZipCodeTree::SEED);

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(5).type == ZipCodeTree::CHAIN_END);
            
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
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            //The tree should be:
            // [pos1 5 pos2] [pos3 5 pos4]
            // of
            // [pos2 5 pos1] [ pos3 5 pos4]
            // etc...
            REQUIRE(zip_tree.get_tree_size() == 10);

            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);

            //first seed 
            REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);

            //Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(2).value == 5);

            //The next seed
            REQUIRE(zip_tree.get_item_at_index(3).type == ZipCodeTree::SEED);

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(4).type == ZipCodeTree::CHAIN_END);

            //Chain start
            REQUIRE(zip_tree.get_item_at_index(5).type == ZipCodeTree::CHAIN_START);

            //The first seed in the new chain
            REQUIRE(zip_tree.get_item_at_index(6).type == ZipCodeTree::SEED);

            //Distance between the seeds
            REQUIRE(zip_tree.get_item_at_index(7).type == ZipCodeTree::EDGE);
            REQUIRE(zip_tree.get_item_at_index(7).value == 5);

            //The last seed
            REQUIRE(zip_tree.get_item_at_index(8).type == ZipCodeTree::SEED);

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(9).type == ZipCodeTree::CHAIN_END);

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
            REQUIRE(reverse_views.size() == 4);
            // The first seed can't see any other seeds
            REQUIRE(reverse_views.count({0, false}));
            REQUIRE(reverse_views[{0, false}].size() == 0);
            // The second seed can see the first seed at distance 5
            REQUIRE(reverse_views.count({1, false}));
            REQUIRE(reverse_views[{1, false}].size() == 1);
            REQUIRE(reverse_views[{1, false}][0].seed == 0);
            REQUIRE(reverse_views[{1, false}][0].distance == 5);
            REQUIRE(reverse_views[{1, false}][0].is_reverse == false);
            // The third seed can't see any other seeds
            REQUIRE(reverse_views.count({2, false}));
            REQUIRE(reverse_views[{2, false}].size() == 0);
            // The fourth seed can see the third seed at distance 5
            REQUIRE(reverse_views.count({3, false}));
            REQUIRE(reverse_views[{3, false}].size() == 1);
            REQUIRE(reverse_views[{3, false}][0].seed == 2);
            REQUIRE(reverse_views[{3, false}][0].distance == 5);
            REQUIRE(reverse_views[{3, false}][0].is_reverse == false);
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



        ofstream out ("testGraph.hg");
        graph.serialize(out);
        
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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            //The tree should be:
            // [pos1 3 pos3 6 pos6]
            //or backwards
            REQUIRE(zip_tree.get_tree_size() == 7);

            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);

            //first seed 
            REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(1).is_reversed) {
                REQUIRE(zip_tree.get_item_at_index(1).value == 2);
            } else {
                REQUIRE(zip_tree.get_item_at_index(1).value == 0);
            }

            //distance between them
            REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(2).value == 3 ||
                    zip_tree.get_item_at_index(2).value == 6));

            //the next seed
            REQUIRE(zip_tree.get_item_at_index(3).type == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(3).value == 1);

            //distance between them
            REQUIRE(zip_tree.get_item_at_index(4).type == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(4).value == 3 ||
                    zip_tree.get_item_at_index(4).value == 6));

            //the last seed
            REQUIRE(zip_tree.get_item_at_index(5).type == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(5).is_reversed) {
                REQUIRE(zip_tree.get_item_at_index(5).value == 0);
            } else {
                REQUIRE(zip_tree.get_item_at_index(5).value == 2);
            }

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(6).type == ZipCodeTree::CHAIN_END);
            
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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            //The tree should be:
            // [pos1 3 pos3 6 pos6]
            //or backwards
            REQUIRE(zip_tree.get_tree_size() == 7);

            //Chain start
            REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);

            //first seed 
            //This is either the first seed on 1 going backwards, or the third seed on 6 going backwards
            REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(1).value == 0) {
                REQUIRE(zip_tree.get_item_at_index(1).is_reversed);
            } else {
                REQUIRE(zip_tree.get_item_at_index(1).value == 2);
                REQUIRE(zip_tree.get_item_at_index(1).is_reversed);
            }

            //distance between them
            REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(2).value == 3 ||
                    zip_tree.get_item_at_index(2).value == 6));

            //the next seed
            REQUIRE(zip_tree.get_item_at_index(3).type == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(3).value == 1);

            //distance between them
            REQUIRE(zip_tree.get_item_at_index(4).type == ZipCodeTree::EDGE);
            REQUIRE((zip_tree.get_item_at_index(4).value == 3 ||
                    zip_tree.get_item_at_index(4).value == 6));

            //the last seed
            REQUIRE(zip_tree.get_item_at_index(5).type == ZipCodeTree::SEED);
            if (zip_tree.get_item_at_index(5).value == 0) {
                REQUIRE(!zip_tree.get_item_at_index(5).is_reversed);
            } else {
                REQUIRE(zip_tree.get_item_at_index(5).value == 2);
                REQUIRE(!zip_tree.get_item_at_index(5).is_reversed);
            }

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(6).type == ZipCodeTree::CHAIN_END);

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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

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
        SECTION( "Only snarls in a snarl" ) {
 
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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

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
            //all are in the same cluster
            vector<SnarlDistanceIndexClusterer::Seed> seeds;
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 3);
                REQUIRE(dag_non_dag_count.second == 0);
            }
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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

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

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 3);
                REQUIRE(dag_non_dag_count.second == 0);
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
            for (pos_t pos : positions) {
                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);
                seeds.push_back({ pos, 0, zipcode});
            }

            ZipCodeTree zip_tree;
            zip_tree.fill_in_tree(seeds, distance_index);
            zip_tree.print_self();
            zip_tree.validate_zip_tree(distance_index);

            SECTION( "Count dags" ) {
                pair<size_t, size_t> dag_non_dag_count = zip_tree.dag_and_non_dag_snarl_count(seeds, distance_index);
                REQUIRE(dag_non_dag_count.first == 1);
                REQUIRE(dag_non_dag_count.second == 1);
            }
        }

    }

    TEST_CASE("zip tree handles complicated nested snarls") {
        
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

        ZipCodeTree zip_tree;
        zip_tree.fill_in_tree(seeds, distance_index);
        zip_tree.print_self();
        zip_tree.validate_zip_tree(distance_index);
    }
}
}
