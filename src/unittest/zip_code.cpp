#include "catch.hpp"
#include <stdio.h>
#include <iostream>
#include "../zip_code.hpp"
#include "../integrated_snarl_finder.hpp"

namespace vg{
namespace unittest{
using namespace std;

    TEST_CASE("One node graph", "[zipcode]") {
        VG graph;
 
        Node* n1 = graph.create_node("GCAAACAGATT");
 
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION ("zip code") {
            zip_code_t zip_code;
            zip_code.fill_in_zip_code(distance_index, make_pos_t(n1->id(), 0, false));

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zip_code.zip_code.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the rank of the node (chain) in the root-snarl
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Third value is the length of the node
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 11);

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
    }
    TEST_CASE("Simple chain graph", "[zipcode]") {
        VG graph;
 
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("GCAA");
        Node* n3 = graph.create_node("GCA");
        Node* n4 = graph.create_node("TT");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("GCAAA");


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
        bool chain_is_reversed = distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id()));

        SECTION ("zip code for node on top-level chain") {
            zip_code_t zip_code;
            zip_code.fill_in_zip_code(distance_index, make_pos_t(n1->id(), 0, false));

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zip_code.zip_code.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the node code
            //Third value is the prefix sum of the node

            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n1->id())));

            //Fourth is the node length
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3);

            //Fifth is if the node is reversed
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id())));

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
        SECTION ("zip code for node in simple snarl") {
            zip_code_t zip_code;
            zip_code.fill_in_zip_code(distance_index, make_pos_t(n4->id(), 0, false));

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zip_code.zip_code.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the snarl code

            //1 for a regular snarl
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //prefix sum of the snarl
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (chain_is_reversed ? 5 : 6));

            //length of the snarl
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //node is reversed in the snarl
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(distance_index.get_parent(
                                                distance_index.get_node_net_handle(n4->id()))));

            //Next is the chain code
            //rank of the chain in the snarl
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(distance_index.get_parent(
                                                distance_index.get_node_net_handle(n4->id()))));

            //node length
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 2);

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
    }
    TEST_CASE("Nested snarl graph", "[zipcode]") {
 
        // This graph will have a snarl from 1 to 8, a snarl from 2 to 7,
        // and a snarl from 3 to 5, all nested in each other.
        VG graph;
 
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
 
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n8);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n6);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);

 
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);
        bool chain_is_reversed = distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id()));

        SECTION ("zip code for node on top-level chain") {
            zip_code_t zip_code;
            zip_code.fill_in_zip_code(distance_index, make_pos_t(n1->id(), 0, false));

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zip_code.zip_code.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the node code
            //Third value is the prefix sum of the node

            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n1->id())));

            //Fourth is the node length
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3);

            //Fifth is if the node is reversed
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id())));

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
        SECTION ("zip code for node on in nested chain") {
            zip_code_t zip_code;
            zip_code.fill_in_zip_code(distance_index, make_pos_t(n2->id(), 0, false));

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zip_code.zip_code.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the regular snarl code

            //1 for regular snarl tag
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //Prefix sum of the snarl
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (chain_is_reversed ? 4 : 3));

            //snarl length
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Is the chain is reversed in the snarl 
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(
                                                distance_index.get_parent(distance_index.get_node_net_handle(n2->id()))));
            //Next is the chain code
            //rank in snarl
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(
                                                distance_index.get_parent(distance_index.get_node_net_handle(n2->id()))));
            
            //chain length
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3);

            //Next is the node code
            //Offset of the node in the chain
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n2->id())));

            //length of the node
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //is the node reversed in the parent
            value_and_index = zip_code.zip_code.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n2->id())));

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
    }
}
}
