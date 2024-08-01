#include "catch.hpp"
#include <stdio.h>
#include <iostream>
#include "../zip_code.hpp"
#include "../integrated_snarl_finder.hpp"

namespace vg{
namespace unittest{
using namespace std;

    TEST_CASE("One node zipcode", "[zipcode][bug]") {
        VG graph;
 
        Node* n1 = graph.create_node("GCAAACAGATT");
 
        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);

        SECTION ("zip code") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));

            //1st value is 1 to indicate that it's a chain
            REQUIRE(zipcode.zipcode.at(0) == 1);

            //Second value is the rank of the node (chain) in the root-snarl
            REQUIRE(zipcode.zipcode.at(1) == 0);

            //Third value is the length of the node
            REQUIRE(zipcode.zipcode.at(2) == 11+1);

            //Connectivity
            REQUIRE(zipcode.zipcode.at(3) == 0);

            //That's it
            REQUIRE(zipcode.zipcode.size() == 4);


        }
        SECTION("decoder") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 1);
            REQUIRE(zipcode.decoder.front().first == 1);
            REQUIRE(zipcode.decoder.front().second == 0);
        }
        SECTION("decoded code") {
            cerr << "New code" << endl;
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            zipcode.fill_in_full_decoder();

            net_handle_t chain1 = distance_index.get_parent(distance_index.get_node_net_handle(n1->id()));

            REQUIRE(zipcode.get_length(0) == distance_index.minimum_length(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_NODE);
        }
        SECTION("n1 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            }
        }
        SECTION("Distances within one node") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            zipcode.fill_in_full_decoder();
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
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 2);

            //1st value is 1 to indicate that it's a chain
            REQUIRE(zipcode.zipcode.at(0) == 1);
            REQUIRE(zipcode.decoder[0] == std::make_pair(true, (size_t)0));

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == 0);
            
            //Component count of the chain
            REQUIRE(zipcode.zipcode.at(2) == 0);

            //Connectivity of the chain
            REQUIRE(zipcode.zipcode.at(3) == 0);

            //Next is the node code
            //Third value is the prefix sum of the node

            REQUIRE(zipcode.decoder[1] == std::make_pair(true,(size_t)4));
            REQUIRE(zipcode.zipcode.at(4) == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n1->id()))+1);

            //Fourth is the node length
            REQUIRE(zipcode.zipcode.at(5) == 3+1);

            //Fifth is if the node is reversed
            REQUIRE(zipcode.zipcode.at(6) == distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id())));

            //The component
            REQUIRE(zipcode.zipcode.at(7) == 0);

            //That's it
            REQUIRE(zipcode.zipcode.size() == 7);

        }
        SECTION ("decoded zip code for node on top-level chain") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            zipcode.fill_in_full_decoder();
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
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 3);

            //1st value is 1 to indicate that it's a chain
            REQUIRE(zipcode.zipcode.at(0) == 1);
            REQUIRE(zipcode.decoder.at(0) == std::make_pair(true, (size_t)0));

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == 0);

            //Chain component count
            REQUIRE(zipcode.zipcode.at(2) == 0);

            //Connectivity of the chain
            REQUIRE(zipcode.zipcode.at(3) == 0);

            //Next is the snarl code

            //1 for a regular snarl
            REQUIRE(zipcode.decoder.at(1) == std::make_pair(false, (size_t)4));
            REQUIRE(zipcode.zipcode.at(4) == 1);

            //prefix sum of the snarl
            REQUIRE(zipcode.zipcode.at(5) == (chain_is_reversed ? 5 : 6)+1);

            //length of the snarl
            REQUIRE(zipcode.zipcode.at(6) == 1+1);

            //Child count
            REQUIRE(zipcode.zipcode.at(7) == 2);

            //Chain component
            REQUIRE(zipcode.zipcode.at(8) == 0);

            //node is reversed in the snarl
            net_handle_t chain4 = distance_index.get_parent(distance_index.get_node_net_handle(n4->id()));
            net_handle_t snarl = distance_index.get_parent(chain4);
            bool is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(chain4)) != 0;
            REQUIRE(zipcode.zipcode.at(9) == is_rev);

            //Next is the chain code
            //rank of the chain in the snarl
            REQUIRE(zipcode.decoder[2] == std::make_pair(true, (size_t)10));
            REQUIRE(zipcode.zipcode.at(10) == distance_index.get_rank_in_parent(distance_index.get_parent(
                                                distance_index.get_node_net_handle(n4->id()))));

            //node length
            REQUIRE(zipcode.zipcode.at(11) == 2+1);

            //chain component count
            REQUIRE(zipcode.zipcode.at(12) == 0);

            //That's it
            REQUIRE(zipcode.zipcode.size() == 12);


        }
        SECTION ("decoded zip code for node in simple snarl") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            zipcode.fill_in_full_decoder();


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
            zip1.fill_in_full_decoder();
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            zip2.fill_in_full_decoder();
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            zip3.fill_in_full_decoder();
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            zip4.fill_in_full_decoder();
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            zip5.fill_in_full_decoder();
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            zip6.fill_in_full_decoder();

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
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                 ZipCode decoded;
                 decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            }
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
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
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 2); 

            REQUIRE(zipcode.decoder[0] == std::make_pair(true, (size_t)0));
            //1st value is 1 to indicate that it's a chain
            REQUIRE(zipcode.zipcode.at(0) == 1);

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == 0);

            //Third value is the chain component count
            REQUIRE(zipcode.zipcode.at(2) == 0);

            //Connectivity of the chain
            REQUIRE(zipcode.zipcode.at(3) == 0);

            //Next is the node code
            //Third value is the prefix sum of the node

            REQUIRE(zipcode.decoder[1] == std::make_pair(true, (size_t) 4));

            REQUIRE(zipcode.zipcode.at(4) == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n1->id()))+1);

            //Fourth is the node length
            REQUIRE(zipcode.zipcode.at(5) == 3+1);

            //Fifth is if the node is reversed
            REQUIRE(zipcode.zipcode.at(6) == distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id())));

            //component
            REQUIRE(zipcode.zipcode.at(7) == distance_index.get_chain_component(
                                                distance_index.get_node_net_handle(n1->id())));

            //That's it
            REQUIRE(zipcode.zipcode.size() == 7);


        }
        SECTION ("decode zip code for node on top-level chain") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            zipcode.fill_in_full_decoder();

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
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 4); 

            REQUIRE(zipcode.decoder[0] == std::make_pair(true, (size_t)0));
            //1st value is 1 to indicate that it's a chain
            REQUIRE(zipcode.zipcode.at(0) == 1);

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == 0);

            //Third value is the chain component count of the chain
            REQUIRE(zipcode.zipcode.at(2) == 0);

            //Connectivity of the chain
            REQUIRE(zipcode.zipcode.at(3) == 0);

            //Next is the regular snarl code
            REQUIRE(zipcode.decoder[1] == std::make_pair(false, (size_t) 4));

            //1 for regular snarl tag
            REQUIRE(zipcode.zipcode.at(4) == 1);

            //Prefix sum of the snarl
            REQUIRE(zipcode.zipcode.at(5) == (chain_is_reversed ? 4 : 3)+1);

            //snarl length
            REQUIRE(zipcode.zipcode.at(6) == 0+1);

            //Snarl child count
            REQUIRE(zipcode.zipcode.at(7) == 1);

            //chain component
            REQUIRE(zipcode.zipcode.at(8) == distance_index.get_chain_component(distance_index.get_node_net_handle(n2->id())));

            //Is the chain is reversed in the snarl 
            net_handle_t chain2 = distance_index.get_parent(distance_index.get_node_net_handle(n2->id()));
            net_handle_t snarl = distance_index.get_parent(chain2);
            bool is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain2))) != 0;
            REQUIRE(zipcode.zipcode.at(9) == is_rev);

            //Next is the chain code
            REQUIRE(zipcode.decoder[2] == std::make_pair(true, (size_t) 10));

            //rank in snarl
            REQUIRE(zipcode.zipcode.at(10) == distance_index.get_rank_in_parent(
                                                distance_index.get_parent(distance_index.get_node_net_handle(n2->id()))));
            
            //chain length
            REQUIRE(zipcode.zipcode.at(11) == 3+1);

            //chain component count
            REQUIRE(zipcode.zipcode.at(12) == 0);

            //Next is the node code
            REQUIRE(zipcode.decoder[3] == std::make_pair(true, (size_t) 13));
            //Offset of the node in the chain
            REQUIRE(zipcode.zipcode.at(13) == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n2->id()))+1);

            //length of the node
            REQUIRE(zipcode.zipcode.at(14) == 1+1);

            //is the node reversed in the parent
            REQUIRE(zipcode.zipcode.at(15) == distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n2->id())));

            //chain component
            REQUIRE(zipcode.zipcode.at(16) == distance_index.get_chain_component(distance_index.get_node_net_handle(n2->id())));

            //That's it
            REQUIRE(zipcode.zipcode.size() == 16);


        }
        SECTION ("decode zip code for node on in nested chain") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            zipcode.fill_in_full_decoder();

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
            zipcode.fill_in_full_decoder();
            REQUIRE(zipcode.decoder_length() == 7); 

            REQUIRE(zipcode.decoder[0] == std::make_pair(true, (size_t)0));

            //1st value is 1 to indicate that it's a chain
            REQUIRE(zipcode.zipcode.at(0) == 1);

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == 0);

            //Second value is the chain component count of the chain
            REQUIRE(zipcode.zipcode.at(2) == 0);

            //Connectivity of the chain
            REQUIRE(zipcode.zipcode.at(3) == 0);

            //Next is the regular snarl code for snarl 1-8
            REQUIRE(zipcode.decoder[1] == std::make_pair(false, (size_t) 3));

            //1 for regular snarl tag
            REQUIRE(zipcode.zipcode.at(4) == 1);

            //Prefix sum of the snarl
            REQUIRE(zipcode.zipcode.at(5) == (chain_is_reversed ? 4 : 3)+1);

            //snarl length
            REQUIRE(zipcode.zipcode.at(6) == 0+1);

            //snarl child count
            REQUIRE(zipcode.zipcode.at(7) == 1);

            //Chain component
            REQUIRE(zipcode.zipcode.at(8) == distance_index.get_chain_component(distance_index.get_node_net_handle(n2->id())));

            //Is the chain is reversed in the snarl 
            net_handle_t chain2 = distance_index.get_parent(distance_index.get_node_net_handle(n2->id()));
            net_handle_t snarl = distance_index.get_parent(chain2);
            bool is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain2))) != 0;
            REQUIRE(zipcode.zipcode.at(9) == is_rev);
            //Next is the chain code for chain 2-7
            REQUIRE(zipcode.decoder[2] == std::make_pair(true, (size_t) 10));
            //rank in snarl
            REQUIRE(zipcode.zipcode.at(10) == distance_index.get_rank_in_parent(
                                                distance_index.get_parent(distance_index.get_node_net_handle(n2->id()))));
            
            //chain length
            REQUIRE(zipcode.zipcode.at(11) == 3+1);

            //chain component_count
            REQUIRE(zipcode.zipcode.at(12) == 0);

            //Next is the regular snarl code for snarl 2-7
            REQUIRE(zipcode.decoder[3] == std::make_pair(false, (size_t) 13));
            //1 as tag for regular snarl
            REQUIRE(zipcode.zipcode.at(13) == 1);

            //offset in chain
            REQUIRE(zipcode.zipcode.at(14) == 1+1);

            //length
            REQUIRE(zipcode.zipcode.at(15) == 1+1);

            //child count
            REQUIRE(zipcode.zipcode.at(16) == 2);

            //is_reversed
            net_handle_t chain3 = distance_index.get_parent(distance_index.get_node_net_handle(n3->id()));
            snarl = distance_index.get_parent(chain3);
            is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain3))) != 0;
            REQUIRE(zipcode.zipcode.at(17) == is_rev);

            REQUIRE(zipcode.zipcode.at(18) == distance_index.get_chain_component(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true))));

            //Chain code for chain 3-5
            REQUIRE(zipcode.decoder[4] == std::make_pair(true, (size_t) 19));
            //Rank in parent
            REQUIRE(zipcode.zipcode.at(19) == distance_index.get_rank_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n3->id()))) );

            //length
            REQUIRE(zipcode.zipcode.at(20) == distance_index.minimum_length(distance_index.get_parent(distance_index.get_node_net_handle(n3->id()))) +1);

            //component_count
            REQUIRE(zipcode.zipcode.at(21) == 0);

            //REgular snarl code for snarl 3-5
            REQUIRE(zipcode.decoder[5] == std::make_pair(false, (size_t) 22));
            REQUIRE(zipcode.zipcode.at(22) == 1);

            //offset in chain
            REQUIRE(zipcode.zipcode.at(23) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n3->id())) ? 3 : 1)+1);

            //length
            REQUIRE(zipcode.zipcode.at(24) == 0+1);

            //child count
            REQUIRE(zipcode.zipcode.at(25) == 1);

            net_handle_t chain4 = distance_index.get_parent(distance_index.get_node_net_handle(n4->id()));
            snarl = distance_index.get_parent(chain4);

            REQUIRE(zipcode.zipcode.at(26) == distance_index.get_chain_component(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true))));

            //is_reversed
            is_rev = distance_index.distance_in_parent(snarl, distance_index.get_bound(snarl, false, true),
                                                                   distance_index.flip(distance_index.canonical(chain4))) != 0;
            REQUIRE(zipcode.zipcode.at(27) == is_rev);

            //Chain code for node 4
            REQUIRE(zipcode.decoder[6] == std::make_pair(true, (size_t) 28));
            //rank in snarl
            REQUIRE(zipcode.zipcode.at(28) == distance_index.get_rank_in_parent(distance_index.get_node_net_handle(n4->id()))) ;

            //length
            REQUIRE(zipcode.zipcode.at(29) == 4+1) ;

            //Chain component
            REQUIRE(zipcode.zipcode.at(30) == 0) ;

            //That's it
            REQUIRE(zipcode.zipcode.size() == 30);


        }

        SECTION ("decoded zip code for more deeply nested node") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            zipcode.fill_in_full_decoder();

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
            zip1.fill_in_full_decoder();
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            zip2.fill_in_full_decoder();
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            zip3.fill_in_full_decoder();
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            zip4.fill_in_full_decoder();
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            zip5.fill_in_full_decoder();
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            zip6.fill_in_full_decoder();
            ZipCode zip7;
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            zip7.fill_in_full_decoder();
            ZipCode zip8;
            zip8.fill_in_zipcode(distance_index, make_pos_t(n8->id(), 0, false));
            zip8.fill_in_full_decoder();

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
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n7 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n8 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n8->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
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

        SECTION ("zip code for node in irregular snarl") { 
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 3); 

            REQUIRE(zipcode.decoder[0] == std::make_pair(true, (size_t)0));

            //1st value is 1 to indicate that it's a chain
            REQUIRE(zipcode.zipcode.at(0) == 1);

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == 0);

            //Third is the chain component count
            REQUIRE(zipcode.zipcode.at(2) == 0);

            //Connectivity of the chain
            REQUIRE(zipcode.zipcode.at(3) == 0);

            //Irregular snarl code for snarl 1-4
            REQUIRE(zipcode.decoder[1] == std::make_pair(false, (size_t) 4));
            //0 as tag for irregular snarl
            REQUIRE(zipcode.zipcode.at(4) == 2);

            net_handle_t irregular_snarl = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(n2->id())));

            //Snarl prefix sum
            net_handle_t bound = distance_index.get_node_from_sentinel(distance_index.get_bound(irregular_snarl, false, true));
            REQUIRE(zipcode.zipcode.at(5) == SnarlDistanceIndex::sum(distance_index.get_prefix_sum_value(bound),
                                                                     distance_index.minimum_length(bound))+1);

            //Snarl length
            REQUIRE(zipcode.zipcode.at(6) == distance_index.minimum_length(irregular_snarl)+1);

            size_t child_count = 0 ; 
            distance_index.for_each_child(irregular_snarl, [&] (const net_handle_t& child) { child_count++; });
            //Snarl child count
            REQUIRE(zipcode.zipcode.at(7) == child_count);

            //component
            REQUIRE(zipcode.zipcode.at(8) == distance_index.get_chain_component(distance_index.get_node_from_sentinel(distance_index.get_bound(irregular_snarl, false, false))));

            //Snarl record offset
            REQUIRE(zipcode.zipcode.at(9) == distance_index.get_record_offset(irregular_snarl));

            //Distance from left side of child to snarl start
            //REQUIRE(zipcode.zipcode.at(10) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 0 : 1));

            //Distance from right side of child to snarl start
            //REQUIRE(zipcode.zipcode.at(11) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 0 : 1));

            //Distance from left side of child to snarl end
            //REQUIRE(zipcode.zipcode.at(12) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 1 : 0));

            //Distance from right side of child to snarl end
            //REQUIRE(zipcode.zipcode.at(13) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 1 : 0));

            //Node 3 as a chain
            REQUIRE(zipcode.decoder[2] == std::make_pair(true, (size_t) 14));
            //Rank in snarl
            REQUIRE(zipcode.zipcode.at(14) == distance_index.get_rank_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n3->id()))));

            //Length
            REQUIRE(zipcode.zipcode.at(15) == 1+1);

            //Component count
            REQUIRE(zipcode.zipcode.at(16) == 0);

            //That's it
            REQUIRE(zipcode.zipcode.size() == 16);
        }
        SECTION ("decode zip code for node in irregular snarl") { 
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            zipcode.fill_in_full_decoder();

            net_handle_t chain3 = distance_index.get_parent(distance_index.get_node_net_handle(n3->id()));
            net_handle_t snarl1 = distance_index.get_parent(chain3);
            net_handle_t chain1 = distance_index.get_parent(snarl1);


            REQUIRE(distance_index.canonical(zipcode.get_net_handle(0, &distance_index)) == 
                    distance_index.canonical(chain1));
            REQUIRE(zipcode.get_code_type(0) == ZipCode::ROOT_CHAIN);

            //Snarl1 at depth 1
            REQUIRE(zipcode.get_offset_in_chain(1, &distance_index) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id())) ? 6 : 3));
            REQUIRE(zipcode.get_code_type(1) == ZipCode::CYCLIC_SNARL);

            //chain3 at depth 3
            REQUIRE(zipcode.get_length(2) == 1);
            REQUIRE(zipcode.get_rank_in_snarl(2) == distance_index.get_rank_in_parent(chain3));
            REQUIRE(zipcode.get_code_type(2) == ZipCode::CHAIN);
            bool snarl_is_rev = distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n1->id()));
            bool chain_is_rev = distance_index.is_reversed_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n3->id())));
            //node1 to left side of node 3
            REQUIRE(zipcode.get_distance_to_snarl_bound(2, !snarl_is_rev, true) ==  1);
            //Node 1 to right side of node 3
            REQUIRE(zipcode.get_distance_to_snarl_bound(2, !snarl_is_rev, false) == 2);
            //node4 to left side of node 3
            REQUIRE(zipcode.get_distance_to_snarl_bound(2, snarl_is_rev, true) ==  std::numeric_limits<size_t>::max());
            //Node 4 to right side of node 3
            REQUIRE(zipcode.get_distance_to_snarl_bound(2, snarl_is_rev, false) == 0);
        }
        SECTION("Distances") {
            ZipCode zip1;
            zip1.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            zip1.fill_in_full_decoder();
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            zip2.fill_in_full_decoder();
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            zip3.fill_in_full_decoder();
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            zip4.fill_in_full_decoder();
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            zip5.fill_in_full_decoder();
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            zip6.fill_in_full_decoder();
            ZipCode zip7;
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            zip7.fill_in_full_decoder();


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
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n7 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
    }

    TEST_CASE("Top-level snarl zipcode", "[zipcode]") {
 
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
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 2); 

            REQUIRE(zipcode.decoder[0] == std::make_pair(false, (size_t)0));

            //0 to indicate that it's a top-level snarl
            REQUIRE(zipcode.zipcode.at(0) == 0);

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == distance_index.get_connected_component_number(distance_index.get_node_net_handle(n1->id())));

            //Next is node 1 as a chain
            REQUIRE(zipcode.decoder[1] == std::make_pair(true, (size_t) 2));
            //rank in snarl
            REQUIRE(zipcode.zipcode.at(2) == distance_index.get_rank_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n1->id()))));
            //length
            REQUIRE(zipcode.zipcode.at(3) == 3+1);
        }
        SECTION ("decoded zip code for node in top-level snarl") { 
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n1->id(), 0, false));
            zipcode.fill_in_full_decoder();


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
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 3); 

            REQUIRE(zipcode.decoder[0] == std::make_pair(false, (size_t)0));

            //0 to indicate that it's a top-level snarl
            REQUIRE(zipcode.zipcode.at(0) == 0);

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == distance_index.get_connected_component_number(distance_index.get_node_net_handle(n1->id())));

            //Next is chain 2-3
            REQUIRE(zipcode.decoder[1] == std::make_pair(true, (size_t) 2));
            //rank in snarl
            REQUIRE(zipcode.zipcode.at(2) == distance_index.get_rank_in_parent(distance_index.get_parent(distance_index.get_node_net_handle(n3->id()))));
            //length
            REQUIRE(zipcode.zipcode.at(3) == 2+1);
            //component count
            REQUIRE(zipcode.zipcode.at(4) == 0);

            //Node 3
            REQUIRE(zipcode.decoder[2] == std::make_pair(true, (size_t) 5));
            //rank in snarl
            REQUIRE(zipcode.zipcode.at(5) == (distance_index.is_reversed_in_parent(distance_index.get_node_net_handle(n3->id())) ? 0 : 1)+1);
            //length
            REQUIRE(zipcode.zipcode.at(6) == 1+1);
        }
        SECTION ("decode zip code for node in chain in top-level snarl") { 
            net_handle_t node3 = distance_index.get_node_net_handle(n3->id());
            net_handle_t chain2 = distance_index.get_parent(node3);
            net_handle_t root_snarl = distance_index.get_parent(chain2);

            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            zipcode.fill_in_full_decoder();


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
            zip1.fill_in_full_decoder();
            ZipCode zip2;
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            zip2.fill_in_full_decoder();
            ZipCode zip3;
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            zip3.fill_in_full_decoder();
            ZipCode zip4;
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            zip4.fill_in_full_decoder();
            ZipCode zip5;
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            zip5.fill_in_full_decoder();
            ZipCode zip6;
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            zip6.fill_in_full_decoder();
            ZipCode zip7;
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            zip7.fill_in_full_decoder();

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
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n7 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
    }
    TEST_CASE("Top-level chain zipcode", "[zipcode]") {
 
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
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.decoder_length() == 2);

            //1st value is 1 to indicate that it's a chain
            REQUIRE(zipcode.zipcode.at(0) == 1);
            REQUIRE(zipcode.decoder[0] == std::make_pair(true, (size_t)0));

            //Second value is the connected component number of the chain
            REQUIRE(zipcode.zipcode.at(1) == 0);

            //Third value is the chain component count
            REQUIRE(zipcode.zipcode.at(2) == 0);

            //Connectivity of the chain
            REQUIRE(zipcode.zipcode.at(3) == 0);

            //Next is the node code
            //Third value is the prefix sum of the node

            REQUIRE(zipcode.decoder[1] == std::make_pair(true, (size_t) 4));
            REQUIRE(zipcode.zipcode.at(4) == distance_index.get_prefix_sum_value(distance_index.get_node_net_handle(n1->id()))+1);

            //Fourth is the node length
            REQUIRE(zipcode.zipcode.at(5) == 3+1);

            //Fifth is if the node is reversed
            REQUIRE(zipcode.zipcode.at(6) == distance_index.is_reversed_in_parent(
                                                distance_index.get_node_net_handle(n1->id())));

            //Chain component
            REQUIRE(zipcode.zipcode.at(7) == distance_index.get_chain_component(
                                                distance_index.get_node_net_handle(n1->id())));

            //That's it
            REQUIRE(zipcode.zipcode.size() == 7);

        }
        SECTION("Distances") {
            ZipCode zip1;
            zip1.fill_in_zipcode(distance_index, make_pos_t(n1->id(), false, 0));
            zip1.fill_in_full_decoder();
            ZipCode zip2;                                                     
            zip2.fill_in_zipcode(distance_index, make_pos_t(n2->id(), false, 0));
            zip2.fill_in_full_decoder();
            ZipCode zip3;                                                     
            zip3.fill_in_zipcode(distance_index, make_pos_t(n3->id(), false, 0));
            ZipCode zip4;                                                     
            zip4.fill_in_zipcode(distance_index, make_pos_t(n4->id(), false, 0));
            ZipCode zip5;                                                     
            zip5.fill_in_zipcode(distance_index, make_pos_t(n5->id(), false, 0));
            ZipCode zip6;                                                     
            zip6.fill_in_zipcode(distance_index, make_pos_t(n6->id(), false, 0));
            ZipCode zip7;                                                     
            zip7.fill_in_zipcode(distance_index, make_pos_t(n7->id(), false, 0));

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
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n2 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n3 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n3->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n4 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n4->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n5 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n6 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n6->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
                ZipCode decoded;
                decoded.fill_in_zipcode_from_payload(payload);
                REQUIRE(zipcode == decoded);
            };
        }
        SECTION("n7 as payload") {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n7->id(), 0, false));
            gbwtgraph::Payload payload = zipcode.get_payload_from_zip();
            if (zipcode.bit_count() <= 112) {
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
    TEST_CASE( "Looping chain zipcode", "[zipcode]" ) {
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
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, &graph, &snarl_finder);
        
        //graph.to_dot(cerr);

        SECTION( "node2" ) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n2->id(), 0, false));
            zipcode.fill_in_full_decoder();
            net_handle_t node2 = distance_index.get_node_net_handle(n2->id());
            net_handle_t parent = distance_index.get_parent(node2);
            net_handle_t bound = distance_index.get_bound(parent, true, false);

            REQUIRE(zipcode.decoder_length() == 2);

            REQUIRE(distance_index.minimum_length(node2) == zipcode.get_length(1));
            REQUIRE(zipcode.get_chain_component(1) == distance_index.get_chain_component(node2));
            REQUIRE(zipcode.get_last_chain_component(0, true) == distance_index.get_chain_component(bound, true));
            REQUIRE(zipcode.get_last_chain_component(0, false) == distance_index.get_chain_component(bound, false));
            REQUIRE(zipcode.get_is_looping_chain(0));
        }

        SECTION( "node5" ) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(distance_index, make_pos_t(n5->id(), 0, false));
            zipcode.fill_in_full_decoder();
            net_handle_t node = distance_index.get_node_net_handle(n5->id());
            net_handle_t parent = distance_index.get_parent(node);
            net_handle_t bound = distance_index.get_bound(parent, true, false);


            REQUIRE(distance_index.minimum_length(node) == zipcode.get_length(zipcode.max_depth()));
        }
    }
    TEST_CASE( "Chain with external connectivity zipcode","[zipcode]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n1, n1, true, false);

        ofstream out ("testGraph.hg");
        graph.serialize(out);

        IntegratedSnarlFinder snarl_finder(graph);
        SnarlDistanceIndex dist_index;
        fill_in_distance_index(&dist_index, &graph, &snarl_finder);


        SECTION( "Check connectivity" ) {
            ZipCode zipcode;
            zipcode.fill_in_zipcode(dist_index, make_pos_t(n2->id(), false, 0));
            zipcode.fill_in_full_decoder();

            REQUIRE(zipcode.get_length(1) == 1);

            if (dist_index.is_reversed_in_parent(dist_index.get_node_net_handle(n1->id()))) {
                REQUIRE(zipcode.is_externally_end_end_connected(0));
            } else {
                REQUIRE(zipcode.is_externally_start_start_connected(0));
            }

        }
    }
}
}
