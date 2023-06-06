#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include "vg/io/json2pb.h"
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

            ZipCodeTree zip_tree(seeds, distance_index);
            zip_tree.print_self();

            REQUIRE(zip_tree.get_tree_size() == 3);
            REQUIRE(zip_tree.get_item_at_index(0).type == ZipCodeTree::CHAIN_START);
            REQUIRE(zip_tree.get_item_at_index(1).type == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(1).value == 0);
            REQUIRE(zip_tree.get_item_at_index(2).type == ZipCodeTree::CHAIN_END);
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

            ZipCodeTree zip_tree(seeds, distance_index);
            zip_tree.print_self();

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

            ZipCodeTree zip_tree(seeds, distance_index);
            zip_tree.print_self();

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

            //THe other seed
            REQUIRE(zip_tree.get_item_at_index(5).type == ZipCodeTree::SEED);
            REQUIRE(zip_tree.get_item_at_index(5).value == 2);

            //Chain end
            REQUIRE(zip_tree.get_item_at_index(6).type == ZipCodeTree::CHAIN_END);
        }
    }

}
}
