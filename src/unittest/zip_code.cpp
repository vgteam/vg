#include "catch.hpp"
#include <stdio.h>
#include <iostream>
#include "../zip_code.hpp"
#include "../integrated_snarl_finder.hpp"

namespace vg{
namespace unittest{
using namespace std;

    TEST_CASE("One node zipcode", "[zipcode]") {
        VG graph;
 
        Node* n1 = graph.create_node("GCAAACAGATT");
 
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION ("zip code") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the rank of the node (chain) in the root-snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Third value is the length of the node
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 11+1);

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
        SECTION("decoding code") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));

            net_handle_t chain1 = distance_index.get_parent(distance_index.get_node_net_handle(n1->id()));


            REQUIRE(zipcode.get_length(0) == distance_index.minimum_length(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_NODE);
        }
        SECTION("n1 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            }
        }
        SECTION("Distances within one node") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            REQUIRE(ZipCode::minimum_distance_between(zipcode, make_pos_t(n1->id(), false, 0),
                                                      zipcode, make_pos_t(n1->id(), false, 3),
                                                         distance_index)
                    == 3);
        }
    }
    TEST_CASE("Simple chain zipcode", "[zipcode]") {
        //Snarl 1-3, snarl 3-6
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
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));

            REQUIRE(zipcode.get_max_depth() == 1);

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the node code
            //Third value is the prefix sum of the node

            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n1->id()))+1);

            //Fourth is the node length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3+1);

            //Fifth is if the node is reversed
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id())));

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());

        }
        SECTION ("decoded zip code for node on top-level chain") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            net_handle_t node1 = distance_index.get_node_net_handle(n1->id());
            net_handle_t chain1 = distance_index.get_parent(node1);


            REQUIRE(distance_index.canonical(zipcode.get_net_handle(0, &distance_index)) == 
                    distance_index.canonical(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_CHAIN);


            //Next is the node code
            REQUIRE(zipcode.get_code_type( 1) == ZipCode::NODE);
            REQUIRE(zipcode.get_length( 1) == distance_index.minimum_length(node1));
            REQUIRE(zipcode.get_offset_in_chain(1) == distance_index.get_prefix_sum_value(node1));
            REQUIRE(zipcode.get_is_reversed_in_parent(1) == distance_index.is_reversed_in_parent(node1));

        }
        SECTION ("zip code for node in simple snarl") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));

            REQUIRE(zipcode.get_max_depth() == 2);

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the snarl code

            //1 for a regular snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //prefix sum of the snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (chain_is_reversed ? 5 : 6)+1);

            //length of the snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1+1);

            //node is reversed in the snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            net_handle_t chain4 = distance_index.get_parent(distance_index.get_node_net_handle(n4->id()));
            net_handle_t snarl = distance_index.get_parent(chain4);
            bool is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(chain4)) != 0;
            REQUIRE(value_and_index.first == is_rev);

            //Next is the chain code
            //rank of the chain in the snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(distance_index.get_parent(
                                                distance_index.get_node_net_handle(n4->id()))));

            //node length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 2+1);

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
        SECTION ("decoded zip code for node in simple snarl") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));


            net_handle_t chain4 = distance_index.get_parent(distance_index.get_node_net_handle(n4->id()));
            net_handle_t snarl36 = distance_index.get_parent(chain4); 
            net_handle_t chain1 = distance_index.get_parent(snarl36);


            REQUIRE(distance_index.canonical(zipcode.get_net_handle(0, &distance_index)) == 
                    distance_index.canonical(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_CHAIN);

            //values for the snarl
            REQUIRE(zipcode.get_length(1) == distance_index.minimum_length(snarl36));
            REQUIRE(zipcode.get_offset_in_chain(1) == (chain_is_reversed ? 5 : 6));
            REQUIRE(zipcode.get_code_type(1) == ZipCode::REGULAR_SNARL);
            bool is_rev = distance_index.distance_in_parent(snarl36, distance_index.get_bound(snarl36, false, true),
                                                                   distance_index.flip(chain4)) != 0;

            //values for the chain
            REQUIRE(zipcode.get_length(2) == distance_index.minimum_length(chain4));
            REQUIRE(zipcode.get_rank_in_snarl(2) == distance_index.get_rank_in_parent(chain4));
            REQUIRE(zipcode.get_code_type(2) == ZipCode::CHAIN);
            REQUIRE(zipcode.get_is_reversed_in_parent(2) == is_rev);
        }
        SECTION("Distances") {
            ZipCode zip1;
            zip1.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));

            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                      zip2, make_pos_t(n2->id(), false, 0),
                                                         distance_index)
                    == 3);

            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                      zip3, make_pos_t(n3->id(), false, 0),
                                                         distance_index)
                    == 3);
            REQUIRE(ZipCode::minimum_distance_between(zip3, make_pos_t(n3->id(), true, 2),
                                                         zip1, make_pos_t(n1->id(), true, 2),
                                                         distance_index)
                    == 3);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip4, make_pos_t(n4->id(), false, 0),
                                                         distance_index)
                    == 6);
            REQUIRE(ZipCode::minimum_distance_between(zip5, make_pos_t(n5->id(), false, 0),
                                                        zip4, make_pos_t(n4->id(), false, 0),
                                                         distance_index)
                    == std::numeric_limits<size_t>::max());
            REQUIRE(ZipCode::minimum_distance_between(zip4, make_pos_t(n4->id(), false, 0),
                                                         zip4, make_pos_t(n4->id(), false, 1),
                                                         distance_index)
                    == 1);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip6, make_pos_t(n6->id(), false, 0),
                                                         distance_index)
                    == 7);
        }
        SECTION("n1 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                 ZipCode decoded;
                 decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            }
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
    }
    TEST_CASE("Nested snarl zipcode", "[zipcode]") {
 
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
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));

            REQUIRE(zipcode.get_max_depth() == 1); 

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the node code
            //Third value is the prefix sum of the node


            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n1->id()))+1);

            //Fourth is the node length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3+1);

            //Fifth is if the node is reversed
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id())));

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
        SECTION ("decode zip code for node on top-level chain") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));

            net_handle_t node1 = distance_index.get_node_net_handle(n1->id());
            net_handle_t chain1 = distance_index.get_parent(node1);


            REQUIRE(distance_index.canonical(zipcode.get_net_handle(0, &distance_index)) == 
                    distance_index.canonical(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_CHAIN);


            REQUIRE(zipcode.get_length(1) == distance_index.minimum_length(node1));
            REQUIRE(zipcode.get_offset_in_chain(1) == distance_index.get_prefix_sum_value(node1));
            REQUIRE(zipcode.get_code_type(1) == ZipCode::NODE);
            REQUIRE(zipcode.get_is_reversed_in_parent(1) == distance_index.is_reversed_in_parent(node1));

        }
        SECTION ("zip code for node on in nested chain") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));

            REQUIRE(zipcode.get_max_depth() == 3); 

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the regular snarl code

            //1 for regular snarl tag
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //Prefix sum of the snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (chain_is_reversed ? 4 : 3)+1);

            //snarl length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0+1);

            //Is the chain is reversed in the snarl 
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            net_handle_t chain2 = distance_index.get_parent(distance_index.get_node_net_handle(n2->id()));
            net_handle_t snarl = distance_index.get_parent(chain2);
            bool is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain2))) != 0;
            REQUIRE(value_and_index.first == is_rev);
            //Next is the chain code
            //rank in snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(
                                                distance_index.get_parent(distance_index.get_node_net_handle(n2->id()))));
            
            //chain length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3+1);

            //Next is the node code
            //Offset of the node in the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n2->id()))+1);

            //length of the node
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1+1);

            //is the node reversed in the parent
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n2->id())));

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }
        SECTION ("decode zip code for node on in nested chain") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));

            net_handle_t node2 = distance_index.get_node_net_handle(n2->id());
            net_handle_t chain2 = distance_index.get_parent(node2);
            net_handle_t snarl1 = distance_index.get_parent(chain2);
            net_handle_t chain1 = distance_index.get_parent(snarl1);


            REQUIRE(distance_index.canonical(zipcode.get_net_handle(0, &distance_index)) == 
                    distance_index.canonical(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_CHAIN);

            //Snarl at depth 1
            REQUIRE(zipcode.get_length(1) == 0);
            REQUIRE(zipcode.get_offset_in_chain(1) == (chain_is_reversed ? 4 : 3));
            REQUIRE(zipcode.get_code_type(1) == ZipCode::REGULAR_SNARL);
            bool is_rev = distance_index.distance_in_parent(snarl1, distance_index.get_bound(snarl1, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain2))) != 0;

            //Chain at depth 2
            REQUIRE(zipcode.get_length(2) == 3);
            REQUIRE(zipcode.get_rank_in_snarl(2) == distance_index.get_rank_in_parent(chain2));
            REQUIRE(zipcode.get_code_type(2) == ZipCode::CHAIN);
            REQUIRE(zipcode.get_is_reversed_in_parent(2) == is_rev);

            //Node at depth 3
            REQUIRE(zipcode.get_length(3) == 1);
            REQUIRE(zipcode.get_offset_in_chain(3) == distance_index.get_prefix_sum_value(node2));
            REQUIRE(zipcode.get_code_type(3) == ZipCode::NODE);
            REQUIRE(zipcode.get_is_reversed_in_parent(3) == distance_index.is_reversed_in_parent(node2));

        }
        SECTION ("zip code for more deeply nested node") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            REQUIRE(zipcode.get_max_depth() == 6); 


            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the regular snarl code for snarl 1-8

            //1 for regular snarl tag
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //Prefix sum of the snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (chain_is_reversed ? 4 : 3)+1);

            //snarl length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0+1);

            //Is the chain is reversed in the snarl 
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            net_handle_t chain2 = distance_index.get_parent(distance_index.get_node_net_handle(n2->id()));
            net_handle_t snarl = distance_index.get_parent(chain2);
            bool is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain2))) != 0;
            REQUIRE(value_and_index.first == is_rev);
            //Next is the chain code for chain 2-7
            //rank in snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(
                                                distance_index.get_parent(distance_index.get_node_net_handle(n2->id()))));
            
            //chain length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3+1);

            //Next is the regular snarl code for snarl 2-7
            //1 as tag for regular snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //offset in chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1+1);

            //length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1+1);

            //is_reversed
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            net_handle_t chain3 = distance_index.get_parent(distance_index.get_node_net_handle(n3->id()));
            snarl = distance_index.get_parent(chain3);
            is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain3))) != 0;
            REQUIRE(value_and_index.first == is_rev);

            //Chain code for chain 3-5
            //Rank in parent
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n3->id()))) );

            //length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.minimum_length(distance_index.get_parent(distance_index.get_node_net_handle(n3->id()))) +1);

            //REgular snarl code for snarl 3-5
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1);

            //offset in chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n3->id())) ? 3 : 1)+1);

            //length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0+1);

            //is_reversed
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            net_handle_t chain4 = distance_index.get_parent(distance_index.get_node_net_handle(n4->id()));
            snarl = distance_index.get_parent(chain4);
            is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain4))) != 0;
            REQUIRE(value_and_index.first == is_rev);

            //Chain code for node 4
            //rank in snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(distance_index.get_node_net_handle(n4->id()))) ;

            //length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 4+1) ;


            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());


        }

        SECTION ("decoded zip code for more deeply nested node") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));

            net_handle_t chain4 = distance_index.get_parent(distance_index.get_node_net_handle(n4->id()));
            net_handle_t snarl3 = distance_index.get_parent(chain4);
            net_handle_t chain3 = distance_index.get_parent(snarl3);
            net_handle_t snarl2 = distance_index.get_parent(chain3);
            net_handle_t chain2 = distance_index.get_parent(snarl2);
            net_handle_t snarl1 = distance_index.get_parent(chain2);
            net_handle_t chain1 = distance_index.get_parent(snarl1);


            REQUIRE(distance_index.canonical(zipcode.get_net_handle(0, &distance_index)) == 
                        distance_index.canonical(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_CHAIN);

            //Snarl at depth 1
            REQUIRE(zipcode.get_length(1) == 0);
            REQUIRE(zipcode.get_offset_in_chain(1) == (chain_is_reversed ? 4 : 3));
            REQUIRE(zipcode.get_code_type(1) == ZipCode::REGULAR_SNARL);
            net_handle_t snarl = distance_index.get_parent(chain2);
            bool is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain2))) != 0;


            //Chain at depth 2
            REQUIRE(zipcode.get_is_reversed_in_parent(2) == is_rev);
            REQUIRE(zipcode.get_length(2) == 3);
            REQUIRE(zipcode.get_rank_in_snarl(2) == distance_index.get_rank_in_parent(chain2));
            REQUIRE(zipcode.get_code_type(2) == ZipCode::CHAIN);


            //Snarl at depth 3
            REQUIRE(zipcode.get_length(3) == 1);
            REQUIRE(zipcode.get_offset_in_chain(3) == 1);
            REQUIRE(zipcode.get_code_type(3) == ZipCode::REGULAR_SNARL);
            snarl = distance_index.get_parent(chain3);
            is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain3))) != 0;

            //Chain at depth 4
            REQUIRE(zipcode.get_is_reversed_in_parent(4) == is_rev);
            REQUIRE(zipcode.get_length(4) == distance_index.minimum_length(chain3));
            REQUIRE(zipcode.get_rank_in_snarl(4) == distance_index.get_rank_in_parent(chain3));
            REQUIRE(zipcode.get_code_type(4) == ZipCode::CHAIN);


            //Snarl3 at depth 5
            REQUIRE(zipcode.get_length(5) == 0);
            REQUIRE(zipcode.get_offset_in_chain(5) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n3->id())) ? 3 : 1));
            REQUIRE(zipcode.get_code_type(5) == ZipCode::REGULAR_SNARL);
            snarl = distance_index.get_parent(chain4);
            is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain4))) != 0;

            //node/chain at depth 6
            REQUIRE(zipcode.get_is_reversed_in_parent(6) == is_rev);
            REQUIRE(zipcode.get_length(6) == 4);
            REQUIRE(zipcode.get_rank_in_snarl(6) == distance_index.get_rank_in_parent(chain4));
            REQUIRE(zipcode.get_code_type(6) == ZipCode::CHAIN);

        }
        SECTION("Distances") {
            ZipCode zip1;
            zip1.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            ZipCode zip7;
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            ZipCode zip8;
            zip8.fill_in_zipcode(distance_index, make_pos_t(n8->id(), 0, false));

            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip2, make_pos_t(n2->id(), false, 0),
                                                         distance_index)
                    == 3);

            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip6, make_pos_t(n6->id(), false, 0),
                                                         distance_index)
                    == 4);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip7, make_pos_t(n7->id(), false, 0),
                                                         distance_index)
                    == 5);
            REQUIRE(ZipCode::minimum_distance_between(zip2, make_pos_t(n2->id(), false, 0),
                                                         zip7, make_pos_t(n7->id(), false, 0),
                                                         distance_index)
                    == 2);
            REQUIRE(ZipCode::minimum_distance_between(zip4, make_pos_t(n4->id(), false, 0),
                                                         zip8, make_pos_t(n8->id(), false, 0),
                                                         distance_index)
                    == 8);
            REQUIRE(ZipCode::minimum_distance_between(zip4, make_pos_t(n4->id(), false, 0),
                                                         zip6, make_pos_t(n6->id(), false, 0),
                                                         distance_index)
                    == std::numeric_limits<size_t>::max());
            REQUIRE(ZipCode::minimum_distance_between(zip4, make_pos_t(n4->id(), false, 0),
                                                         zip8, make_pos_t(n8->id(), true, 0),
                                                         distance_index)
                    == std::numeric_limits<size_t>::max());
            REQUIRE(ZipCode::minimum_distance_between(zip5, make_pos_t(n5->id(), false, 0),
                                                         zip6, make_pos_t(n6->id(), false, 0),
                                                         distance_index)
                    == std::numeric_limits<size_t>::max());
            REQUIRE(ZipCode::minimum_distance_between(zip7, make_pos_t(n7->id(), true, 0),
                                                         zip2, make_pos_t(n2->id(), true, 0),
                                                         distance_index)
                    == 2);
        }
        SECTION("Distance is greater than") {
            ZipCode zip1;
            zip1.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            ZipCode zip7;
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            ZipCode zip8;
            zip8.fill_in_zipcode(distance_index, make_pos_t(n8->id(), 0, false));


            REQUIRE(!ZipCode::is_farther_than(zip1, zip2, 0));
            REQUIRE(!ZipCode::is_farther_than(zip2, zip7, 0));
        }
        SECTION("n1 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n7 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n8 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n8->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
    }
    TEST_CASE("Irregular snarl zipcode", "[zipcode]") {
 
        VG graph;
 
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
 
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n4);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n3, false, true);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n4, n6);
        Edge* e9 = graph.create_edge(n5, n7);
        Edge* e10 = graph.create_edge(n6, n7);

 
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);
        bool chain_is_reversed = distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id()));

        graph.serialize_to_file("test_graph.hg");

        SECTION ("zip code for node in irregular snarl") { 
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));

            REQUIRE(zipcode.get_max_depth() == 2); 


            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Irregular snarl code for snarl 1-4
            //0 as tag for irregular snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            net_handle_t irregular_snarl = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(n2->id())));

            //Snarl prefix sum
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            net_handle_t bound = distance_index.get_node_from_sentinel(distance_index.get_bound(irregular_snarl, false, true));
            REQUIRE(value_and_index.first == SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(bound),
                                                                     distance_index.minimum_length(bound))+1);

            //Snarl length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.minimum_length(irregular_snarl)+1);

            //Snarl record offset
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_record_offset(irregular_snarl));

            //Distance to snarl start
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 0 : 1));

            //Distance to snarl end
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 1 : 0));

            //Node 3 as a chain
            //Rank in snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n3->id()))));

            //Length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1+1);

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());
        }
        SECTION ("decode zip code for node in irregular snarl") { 
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));

            net_handle_t chain3 = distance_index.get_parent(distance_index.get_node_net_handle(n3->id()));
            net_handle_t snarl1 = distance_index.get_parent(chain3);
            net_handle_t chain1 = distance_index.get_parent(snarl1);


            REQUIRE(distance_index.canonical(zipcode.get_net_handle(0, &distance_index)) == 
                    distance_index.canonical(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_CHAIN);

            //Snarl1 at depth 1
            REQUIRE(zipcode.get_offset_in_chain(1, &distance_index) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 6 : 3));
            REQUIRE(zipcode.get_code_type(1) == ZipCode::IRREGULAR_SNARL);

            //chain3 at depth 3
            REQUIRE(zipcode.get_length(2) == 1);
            REQUIRE(zipcode.get_rank_in_snarl(2) == distance_index.get_rank_in_parent(chain3));
            REQUIRE(zipcode.get_code_type(2) == ZipCode::CHAIN);
            REQUIRE(zipcode.get_distance_to_snarl_end(2) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 1 : 0));
            REQUIRE(zipcode.get_distance_to_snarl_start(2) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 0 : 1));
        }
        SECTION("Distances") {
            ZipCode zip1;
            zip1.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            ZipCode zip7;
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));


            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip2, make_pos_t(n2->id(), false, 0),
                                                         distance_index)
                    == 3);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip3, make_pos_t(n3->id(), false, 0),
                                                         distance_index)
                    == 4);
            REQUIRE(ZipCode::minimum_distance_between(zip3, make_pos_t(n3->id(), false, 0),
                                                         zip1, make_pos_t(n1->id(), true, 0),
                                                         distance_index)
                    == 3);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip4, make_pos_t(n4->id(), false, 0),
                                                         distance_index)
                    == 3);

            //Shouldn't take the loop in the chain
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 1),
                                                         zip1, make_pos_t(n1->id(), false, 0),
                                                         distance_index)
                    == std::numeric_limits<size_t>::max());
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 1),
                                                         zip2, make_pos_t(n2->id(), true, 0),
                                                         distance_index)
                    == 5);
            REQUIRE(ZipCode::minimum_distance_between(zip3, make_pos_t(n3->id(), false, 0),
                                                         zip4, make_pos_t(n4->id(), false, 0),
                                                         distance_index)
                    == 1);
            REQUIRE(ZipCode::minimum_distance_between(zip2, make_pos_t(n2->id(), false, 0),
                                                         zip2, make_pos_t(n2->id(), true, 0),
                                                         distance_index)
                    == 3);
            REQUIRE(ZipCode::minimum_distance_between(zip2, make_pos_t(n2->id(), false, 0),
                                                         zip2, make_pos_t(n2->id(), true, 0),
                                                         distance_index)
                    == 3);
            REQUIRE(ZipCode::minimum_distance_between(zip3, make_pos_t(n3->id(), false, 0),
                                                         zip2, make_pos_t(n2->id(), true, 0),
                                                         distance_index)
                    == 2);
            REQUIRE(ZipCode::minimum_distance_between(zip3, make_pos_t(n3->id(), true, 0),
                                                         zip2, make_pos_t(n2->id(), true, 0),
                                                         distance_index)
                    == 1);
            REQUIRE(ZipCode::minimum_distance_between(zip4, make_pos_t(n4->id(), false, 1),
                                                         zip4, make_pos_t(n4->id(), false, 0),
                                                         distance_index)
                    == std::numeric_limits<size_t>::max());
        }
        SECTION("n1 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n7 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
    }

    TEST_CASE("Top-level snarl zipcode", "[zipcode][bug]") {
 
        VG graph;
 
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
 
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n2, true, false);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n3, n5, false, true);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n6, n7);

 
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION ("zip code for node in top-level snarl") { 
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));

            REQUIRE(zipcode.get_max_depth() == 1); 


            //0 to indicate that it's a top-level snarl
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 0);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_connected_component_number(distance_index.get_node_net_handle(n1->id())));

            //Next is node 1 as a chain
            //rank in snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n1->id()))));
            //length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3+1);
        }
        SECTION ("decoded zip code for node in top-level snarl") { 
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));


            net_handle_t chain1 = distance_index.get_parent(distance_index.get_node_net_handle(n1->id()));
            net_handle_t root_snarl = distance_index.get_parent(chain1);


            //Root snarl
            REQUIRE(distance_index.canonical(zipcode.get_net_handle(0, &distance_index)) == 
                    distance_index.canonical(distance_index.get_parent(chain1)));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_SNARL);

            //Chain1 at depth 1
            REQUIRE(zipcode.get_length(1) == 3);
            REQUIRE(zipcode.get_rank_in_snarl(1) == distance_index.get_rank_in_parent(chain1));
            REQUIRE(zipcode.get_code_type(1) == ZipCode::CHAIN);
        }
        SECTION ("zip code for node in chain in top-level snarl") { 
            net_handle_t node1 = distance_index.get_node_net_handle(n3->id());
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));

            REQUIRE(zipcode.get_max_depth() == 2); 


            //0 to indicate that it's a top-level snarl
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 0);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_connected_component_number(distance_index.get_node_net_handle(n1->id())));

            //Next is chain 2-3
            //rank in snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_rank_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n3->id()))));
            //length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 2+1);

            //Node 3
            //rank in snarl
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n3->id())) ? 0 : 1)+1);
            //length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 1+1);
        }
        SECTION ("decode zip code for node in chain in top-level snarl") { 
            net_handle_t node3 = distance_index.get_node_net_handle(n3->id());
            net_handle_t chain2 = distance_index.get_parent(node3);
            net_handle_t root_snarl = distance_index.get_parent(chain2);

            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));


            //Root snarl
            REQUIRE(zipcode.get_distance_index_address(0) == distance_index.get_connected_component_number(node3));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_SNARL);

            //chain2 at depth 1
            REQUIRE(zipcode.get_length(1) == 2);
            REQUIRE(zipcode.get_rank_in_snarl(1) == distance_index.get_rank_in_parent(chain2));
            REQUIRE(zipcode.get_code_type(1) == ZipCode::CHAIN);

            //node3 at depth 2
            REQUIRE(zipcode.get_length(2) == 1);
            REQUIRE(zipcode.get_offset_in_chain(2) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n3->id())) ? 0 : 1));
            REQUIRE(zipcode.get_code_type(2) == ZipCode::NODE);
            REQUIRE(zipcode.get_is_reversed_in_parent(2) == distance_index.is_reversed_in_parent(node3));
        }
        SECTION("Distances") {
            ZipCode zip1;
            zip1.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            ZipCode zip7;
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));

            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip2, make_pos_t(n2->id(), false, 0),
                                                         distance_index)
                    == 3);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), true, 0),
                                                         zip2, make_pos_t(n2->id(), false, 0),
                                                         distance_index)
                    == 3);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip3, make_pos_t(n3->id(), false, 0),
                                                         distance_index)
                    == 4);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip3, make_pos_t(n3->id(), true, 0),
                                                         distance_index)
                    == 8);
            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                         zip6, make_pos_t(n6->id(), false, 0),
                                                         distance_index)
                    == std::numeric_limits<size_t>::max());
            REQUIRE(ZipCode::minimum_distance_between(zip6, make_pos_t(n6->id(), false, 0),
                                                         zip7, make_pos_t(n7->id(), false, 0),
                                                         distance_index)
                    == 1);
        }
        SECTION("n1 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n7 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
    }
    TEST_CASE("Top-level chain zipcode", "[zipcode][bug]") {
 
        VG graph;
 
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("TGCGT");
        Node* n7 = graph.create_node("G");
 
        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);

 
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION ("zip code for node on top-level chain") {
            net_handle_t node1 = distance_index.get_node_net_handle(n1->id());
            net_handle_t parent = distance_index.get_parent(node1);
            net_handle_t grandparent = distance_index.get_parent(parent);
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));

            REQUIRE(zipcode.get_max_depth() == 1);

            //1st value is 1 to indicate that it's a chain
            pair<size_t, size_t> value_and_index = zipcode.zipcode.get_value_and_next_index(0);
            REQUIRE(value_and_index.first == 1);

            //Second value is the connected component number of the chain
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 0);

            //Next is the node code
            //Third value is the prefix sum of the node

            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n1->id()))+1);

            //Fourth is the node length
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == 3+1);

            //Fifth is if the node is reversed
            value_and_index = zipcode.zipcode.get_value_and_next_index(value_and_index.second);
            REQUIRE(value_and_index.first == distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id())));

            //That's it
            REQUIRE(value_and_index.second == std::numeric_limits<size_t>::max());

        }
        SECTION("Distances") {
            ZipCode zip1;
            zip1.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            ZipCode zip7;
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));

            REQUIRE(ZipCode::minimum_distance_between(zip1, make_pos_t(n1->id(), false, 0),
                                                      zip2, make_pos_t(n2->id(), false, 0),
                                                     distance_index)
                    == 3);
            REQUIRE(ZipCode::is_farther_than(zip1, zip6, 3));
            REQUIRE(!ZipCode::is_farther_than(zip1, zip6, 5));
            REQUIRE(ZipCode::is_farther_than(zip1, zip7, 8));
            REQUIRE(!ZipCode::is_farther_than(zip1, zip7, 10));
            REQUIRE(!ZipCode::is_farther_than(zip2, zip7, 10));
            REQUIRE(ZipCode::is_farther_than(zip2, zip7, 8));
        }
        SECTION("n1 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n7 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.byte_count() <= 15) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("serialization") {
            ZipCodeCollection zipcodes;
            for (size_t i = 1 ; i <= 7 ; i++) {
                ZipCode zip;
                zip.fill_in_zipcode(distance_index, make_pos_t(i, 0, false));
                zipcodes.emplace_back(zip);
            }
            ofstream out ("zipcodes");
            zipcodes.serialize(out);
            out.close();

            ifstream in("zipcodes");
            ZipCodeCollection new_zipcodes;
            new_zipcodes.deserialize(in);
            in.close();

            REQUIRE(zipcodes.size() == new_zipcodes.size());
            for (size_t i = 0 ; i < zipcodes.size() ; i++) {
                REQUIRE(zipcodes.at(i).zipcode == new_zipcodes.at(i).zipcode);
            }
            
        }
    }
}
}
