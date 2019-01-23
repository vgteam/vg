/// \file packed_graph.cpp
///  
/// Unit tests for the PackedGraph class.
///

#include <iostream>

#include "../json2pb.h"
#include "../packed_graph.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

    TEST_CASE("PackedGraph can perform handle graph operations", "[packed][handle]") {
        
        PackedGraph graph;
        
        REQUIRE(graph.node_size() == 0);
        
        handle_t h = graph.create_handle("ATG", 2);
        
        SECTION("PackedGraph has correct structure after creating a node") {
            REQUIRE(graph.get_sequence(h) == "ATG");
            REQUIRE(graph.get_sequence(graph.flip(h)) == "CAT");
            REQUIRE(graph.get_length(h) == 3);
            REQUIRE(graph.has_node(graph.get_id(h)));
            REQUIRE(!graph.has_node(graph.get_id(h) + 1));
            
            REQUIRE(graph.get_handle(graph.get_id(h)) == h);
            REQUIRE(!graph.get_is_reverse(h));
            REQUIRE(graph.get_is_reverse(graph.flip(h)));
            
            REQUIRE(graph.node_size() == 1);
            REQUIRE(graph.min_node_id() == graph.get_id(h));
            REQUIRE(graph.max_node_id() == graph.get_id(h));
            
            graph.follow_edges(h, true, [](const handle_t& prev) {
                REQUIRE(false);
                return true;
            });
            graph.follow_edges(h, false, [](const handle_t& next) {
                REQUIRE(false);
                return true;
            });
        }
        
        handle_t h2 = graph.create_handle("CT", 1);
        
        SECTION("PackedGraph has correct structure after creating a node at the beginning of ID space") {
            
            REQUIRE(graph.get_sequence(h2) == "CT");
            REQUIRE(graph.get_sequence(graph.flip(h2)) == "AG");
            REQUIRE(graph.get_length(h2) == 2);
            REQUIRE(graph.has_node(graph.get_id(h2)));
            REQUIRE(!graph.has_node(max(graph.get_id(h), graph.get_id(h2)) + 1));
            
            REQUIRE(graph.get_handle(graph.get_id(h2)) == h2);
            
            REQUIRE(graph.node_size() == 2);
            REQUIRE(graph.min_node_id() == graph.get_id(h2));
            REQUIRE(graph.max_node_id() == graph.get_id(h));
            
            graph.follow_edges(h2, true, [](const handle_t& prev) {
                REQUIRE(false);
                return true;
            });
            graph.follow_edges(h2, false, [](const handle_t& next) {
                REQUIRE(false);
                return true;
            });
        }
        
        // creating and accessing a node at the end of ID space
        
        handle_t h3 = graph.create_handle("GAC", 4);
        
        SECTION("PackedGraph has correct structure after creating a node at the end of ID space") {
            REQUIRE(graph.get_sequence(h3) == "GAC");
            REQUIRE(graph.get_sequence(graph.flip(h3)) == "GTC");
            REQUIRE(graph.get_length(h3) == 3);
            
            REQUIRE(graph.get_handle(graph.get_id(h3)) == h3);
            
            REQUIRE(graph.node_size() == 3);
            REQUIRE(graph.min_node_id() == graph.get_id(h2));
            REQUIRE(graph.max_node_id() == graph.get_id(h3));
            
            graph.follow_edges(h3, true, [](const handle_t& prev) {
                REQUIRE(false);
                return true;
            });
            graph.follow_edges(h3, false, [](const handle_t& next) {
                REQUIRE(false);
                return true;
            });
        }
        
        
        // creating and accessing in the middle of ID space
        
        handle_t h4 = graph.create_handle("T", 3);
        
        SECTION("PackedGraph has correct structure after creating a node in the middle of ID space") {
            REQUIRE(graph.get_sequence(h4) == "T");
            REQUIRE(graph.get_sequence(graph.flip(h4)) == "A");
            REQUIRE(graph.get_length(h4) == 1);
            
            REQUIRE(graph.get_handle(graph.get_id(h4)) == h4);
            
            REQUIRE(graph.node_size() == 4);
            REQUIRE(graph.min_node_id() == graph.get_id(h2));
            REQUIRE(graph.max_node_id() == graph.get_id(h3));
            
            graph.follow_edges(h4, true, [](const handle_t& prev) {
                REQUIRE(false);
                return true;
            });
            graph.follow_edges(h4, false, [](const handle_t& next) {
                REQUIRE(false);
                return true;
            });
        }
        
        graph.create_edge(h, h2);
        
        bool found1 = false, found2 = false, found3 = false, found4 = false;
        int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
        
        SECTION("PackedGraph has correct structure after creating an edge") {
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found4 = true;
                }
                count4++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            
            count1 = count2 = count3 = count4 = 0;
            found1 = found2 = found3 = found4 = false;
        }
        
        graph.create_edge(h, graph.flip(h3));
        
        bool found5 = false, found6 = false, found7 = false, found8 = false;
        int count5 = 0, count6 = 0;
        
        SECTION("PackedGraph has correct structure after creating an edge with a traversal") {
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
        }
        
        graph.create_edge(h4, graph.flip(h4));
        
        SECTION("PackedGraph has correct structure after creating a reversing self-loop") {
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        graph.create_edge(h, graph.flip(h4));
        graph.create_edge(graph.flip(h3), h4);
        
        graph.destroy_edge(h, graph.flip(h4));
        graph.destroy_edge(graph.flip(h3), h4);
        
        SECTION("PackedGraph has correct structure after creating and deleting edges") {
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        handle_t h5 = graph.create_handle("GGACC");
        
        // make some edges to ensure that deleting is difficult
        graph.create_edge(h, h5);
        graph.create_edge(h5, h);
        graph.create_edge(graph.flip(h5), h2);
        graph.create_edge(h3, graph.flip(h5));
        graph.create_edge(h3, h5);
        graph.create_edge(h5, h4);
        
        graph.destroy_handle(h5);
        
        SECTION("PackedGraph has correct structure after creating and deleting a node") {
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        graph.swap_handles(h, h2);
        graph.swap_handles(h2, h3);
        
        SECTION("PackedGraph has correct structure after swapping nodes") {
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        SECTION("PackedGraph visits all nodes with for_each_handle") {
            graph.for_each_handle([&](const handle_t& handle) {
                if (handle == h) {
                    found1 = true;
                }
                else if (handle == h2) {
                    found2 = true;
                }
                else if (handle == h3) {
                    found3 = true;
                }
                else if (handle == h4) {
                    found4 = true;
                }
                else {
                    REQUIRE(false);
                }
                return true;
            });
            
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            
            found1 = found2 = found3 = found4 = false;
        }
        
        // to make sure the sequence reverse complemented correctly
        auto check_rev_comp = [](const std::string& seq1, const std::string& seq2) {
            REQUIRE(seq1.size() == seq2.size());
            auto it = seq1.begin();
            auto rit = seq2.rbegin();
            for (; it != seq1.end(); it++) {
                if (*it == 'A') {
                    REQUIRE(*rit == 'T');
                }
                else if (*it == 'C') {
                    REQUIRE(*rit == 'G');
                }
                else if (*it == 'G') {
                    REQUIRE(*rit == 'C');
                }
                else if (*it == 'T') {
                    REQUIRE(*rit == 'A');
                }
                else if (*it == 'N') {
                    REQUIRE(*rit == 'N');
                }
                else {
                    REQUIRE(false);
                }
                
                rit++;
            }
        };
        
        
        int count7 = 0, count8 = 0;
        
        SECTION("PackedGraph correctly reverses a node") {
            
            string seq1 = graph.get_sequence(h);
            h = graph.apply_orientation(graph.flip(h));
            
            // check the sequence
            string rev_seq1 = graph.get_sequence(h);
            check_rev_comp(seq1, rev_seq1);
            
            // check that the edges are what we expect
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                count1++;
                return true;
            });
            graph.follow_edges(h, true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found1 = true;
                }
                else if (prev == h3) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& next) {
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h), false, [&](const handle_t& prev) {
                if (prev == h2) {
                    found3 = true;
                }
                else if (prev == graph.flip(h3)) {
                    found4 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == graph.flip(h)) {
                    found5 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == h) {
                    found6 = true;
                }
                count6++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h)) {
                    found7 = true;
                }
                count7++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == h) {
                    found8 = true;
                }
                count8++;
                return true;
            });
            REQUIRE(count1 == 0);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 0);
            REQUIRE(count4 == 2);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(count7 == 1);
            REQUIRE(count8 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = count7 = count8 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            
            // and now switch it back to the same orientation and repeat the topology checks
            
            h = graph.apply_orientation(graph.flip(h));
            
            graph.swap_handles(h, h2);
            graph.swap_handles(h2, h3);
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        vector<size_t> offsets{1, 2};
        vector<handle_t> parts = graph.divide_handle(h, offsets);
        
        int count9 = 0, count10 = 0, count11 = 0, count12 = 0;
        bool found9 = false, found10 = false, found11 = false, found12 = false, found13 = false, found14 = false;
        
        SECTION("PackedGraph can correctly divide a node") {
            
            REQUIRE(parts.size() == 3);
            
            REQUIRE(graph.get_sequence(parts[0]) == "A");
            REQUIRE(graph.get_length(parts[0]) == 1);
            REQUIRE(graph.get_sequence(parts[1]) == "T");
            REQUIRE(graph.get_length(parts[1]) == 1);
            REQUIRE(graph.get_sequence(parts[2]) == "G");
            REQUIRE(graph.get_length(parts[2]) == 1);
            
            
            graph.follow_edges(parts[0], false, [&](const handle_t& next) {
                if (next == parts[1]) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(parts[0], true, [&](const handle_t& prev) {
                count2++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[0]), true, [&](const handle_t& prev) {
                if (prev == graph.flip(parts[1])) {
                    found2 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[0]), false, [&](const handle_t& next) {
                count4++;
                return true;
            });
            
            graph.follow_edges(parts[1], false, [&](const handle_t& next) {
                if (next == parts[2]) {
                    found3 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(parts[1], true, [&](const handle_t& prev) {
                if (prev == parts[0]) {
                    found4 = true;
                }
                count6++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[1]), true, [&](const handle_t& prev) {
                if (prev == graph.flip(parts[2])) {
                    found5 = true;
                }
                count7++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[1]), false, [&](const handle_t& next) {
                if (next == graph.flip(parts[0])) {
                    found6 = true;
                }
                count8++;
                return true;
            });
            
            graph.follow_edges(parts[2], false, [&](const handle_t& next) {
                if (next == h2) {
                    found7 = true;
                }
                else if (next == graph.flip(h3)) {
                    found8 = true;
                }
                count9++;
                return true;
            });
            graph.follow_edges(parts[2], true, [&](const handle_t& prev) {
                if (prev == parts[1]) {
                    found9 = true;
                }
                count10++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[2]), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found10 = true;
                }
                else if (prev == h3) {
                    found11 = true;
                }
                count11++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[2]), false, [&](const handle_t& next) {
                if (next == graph.flip(parts[1])) {
                    found12 = true;
                }
                count12++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == parts[2]) {
                    found13 = true;
                }
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == parts[2]) {
                    found14 = true;
                }
                return true;
            });
            
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 0);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 0);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(count7 == 1);
            REQUIRE(count8 == 1);
            REQUIRE(count9 == 2);
            REQUIRE(count10 == 1);
            REQUIRE(count11 == 2);
            REQUIRE(count12 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            REQUIRE(found9);
            REQUIRE(found10);
            REQUIRE(found11);
            REQUIRE(found12);
            REQUIRE(found13);
            REQUIRE(found14);
            
            count1 = count2 = count3 = count4 = count5 = count6 = count7 = count8 = count9 = count10 = count11 = count12 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = found9 = found10 = found11 = found12 = false;
        }
        
        std::vector<handle_t> rev_parts = graph.divide_handle(graph.flip(h3), {1});
        
        SECTION("PackedGraph can correctly divide a node on the reverse strand") {
            
            REQUIRE(graph.get_sequence(rev_parts[0]) == "G");
            REQUIRE(graph.get_length(rev_parts[0]) == 1);
            REQUIRE(graph.get_is_reverse(rev_parts[0]));
            REQUIRE(graph.get_sequence(rev_parts[1]) == "TC");
            REQUIRE(graph.get_length(rev_parts[1]) == 2);
            REQUIRE(graph.get_is_reverse(rev_parts[1]));
            
            graph.follow_edges(rev_parts[0], false, [&](const handle_t& next) {
                if (next == rev_parts[1]) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(rev_parts[1], true, [&](const handle_t& prev) {
                if (prev == rev_parts[0]) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(graph.flip(rev_parts[1]), false, [&](const handle_t& next) {
                if (next == graph.flip(rev_parts[0])) {
                    found3 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(rev_parts[0]), true, [&](const handle_t& prev) {
                if (prev == graph.flip(rev_parts[1])) {
                    found4 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(rev_parts[0], true, [&](const handle_t& prev) {
                if (prev == parts[2]) {
                    found5 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(rev_parts[1], false, [&](const handle_t& next) {
                count6++;
                return true;
            });
            
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 0);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
        }
    }
    
    TEST_CASE("PackedGraph can perform path handle graph operations", "[packed][handle]") {
        
        PackedGraph graph;
        
        auto check_path = [&](const path_handle_t& p, const vector<handle_t>& occs) {
            
            occurrence_handle_t occ;
            for (int i = 0; i < occs.size(); i++){
                if (i == 0) {
                    occ = graph.get_first_occurrence(p);
                }
                
                REQUIRE(graph.get_path_handle_of_occurrence(occ) == p);
                REQUIRE(graph.get_occurrence(occ) == occs[i]);
                REQUIRE(graph.has_previous_occurrence(occ) == (i > 0));
                REQUIRE(graph.has_next_occurrence(occ) == (i < occs.size() - 1));
                
                if (i != occs.size() - 1) {
                    occ = graph.get_next_occurrence(occ);
                }
            }
            
            for (int i = occs.size() - 1; i >= 0; i--){
                if (i == occs.size() - 1) {
                    occ = graph.get_last_occurrence(p);
                }
                
                REQUIRE(graph.get_path_handle_of_occurrence(occ) == p);
                REQUIRE(graph.get_occurrence(occ) == occs[i]);
                REQUIRE(graph.has_previous_occurrence(occ) == (i > 0));
                REQUIRE(graph.has_next_occurrence(occ) == (i < occs.size() - 1));
                
                if (i != 0) {
                    occ = graph.get_previous_occurrence(occ);
                }
            }
        };
        
        auto check_flips = [&](const path_handle_t& p, const vector<handle_t>& occs) {
            auto flipped = occs;
            for (size_t i = 0; i < occs.size(); i++) {
                
                graph.apply_orientation(graph.flip(graph.forward(flipped[i])));
                flipped[i] = graph.flip(flipped[i]);
                check_path(p, flipped);
                
                graph.apply_orientation(graph.flip(graph.forward(flipped[i])));
                flipped[i] = graph.flip(flipped[i]);
                check_path(p, flipped);
            }
        };
        
        
        handle_t h1 = graph.create_handle("AC");
        handle_t h2 = graph.create_handle("CAGTGA");
        handle_t h3 = graph.create_handle("GT");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h1, graph.flip(h2));
        graph.create_edge(graph.flip(h2), h3);
        
        REQUIRE(!graph.has_path("1"));
        REQUIRE(graph.get_path_count() == 0);
        
        path_handle_t p1 = graph.create_path_handle("1");
        
        REQUIRE(graph.has_path("1"));
        REQUIRE(graph.get_path_count() == 1);
        REQUIRE(graph.get_path_handle("1") == p1);
        REQUIRE(graph.get_path_name(p1) == "1");
        REQUIRE(graph.get_occurrence_count(p1) == 0);
        REQUIRE(graph.is_empty(p1));
        
        graph.append_occurrence(p1, h1);
        
        REQUIRE(graph.get_occurrence_count(p1) == 1);
        REQUIRE(!graph.is_empty(p1));
        
        graph.append_occurrence(p1, h2);
        graph.append_occurrence(p1, h3);
        
        REQUIRE(graph.get_occurrence_count(p1) == 3);
        
        SECTION("PackedGraph can traverse a path") {
            check_path(p1, {h1, h2, h3});
        }
        
        SECTION("PackedGraph preserves paths when reversing nodes") {
            check_flips(p1, {h1, h2, h3});
        }
        
        path_handle_t p2 = graph.create_path_handle("2");
        REQUIRE(graph.get_path_count() == 2);
        
        graph.append_occurrence(p2, h1);
        graph.append_occurrence(p2, graph.flip(h2));
        graph.append_occurrence(p2, h3);
        
        check_path(p2, {h1, graph.flip(h2), h3});
        
        SECTION("PackedGraph can query occurrences of a node on paths") {
            
            bool found1 = false, found2 = false;
            vector<occurrence_handle_t> occs = graph.occurrences_of_handle(h1);
            for (auto& occ : occs) {
                if (graph.get_path_handle_of_occurrence(occ) == p1 &&
                    graph.get_occurrence(occ) == h1) {
                    found1 = true;
                }
                else if (graph.get_path_handle_of_occurrence(occ) == p2 &&
                         graph.get_occurrence(occ) == h1) {
                    found2 = true;
                }
                else {
                    REQUIRE(false);
                }
            }
            REQUIRE(found1);
            REQUIRE(found2);
            found1 = found2 = false;
            
            occs = graph.occurrences_of_handle(h1, true);
            for (auto& occ : occs) {
                if (graph.get_path_handle_of_occurrence(occ) == p1 &&
                    graph.get_occurrence(occ) == h1) {
                    found1 = true;
                }
                else if (graph.get_path_handle_of_occurrence(occ) == p2 &&
                         graph.get_occurrence(occ) == h1) {
                    found2 = true;
                }
                else {
                    REQUIRE(false);
                }
            }
            REQUIRE(found1);
            REQUIRE(found2);
            found1 = found2 = false;
            
            occs = graph.occurrences_of_handle(graph.flip(h1), true);
            for (auto& occ : occs) {
                REQUIRE(false);
            }
            
            occs = graph.occurrences_of_handle(h2, true);
            for (auto& occ : occs) {
                if (graph.get_path_handle_of_occurrence(occ) == p1 &&
                    graph.get_occurrence(occ) == h2) {
                    found1 = true;
                }
                else {
                    REQUIRE(false);
                }
            }
            occs = graph.occurrences_of_handle(graph.flip(h2), true);
            for (auto& occ : occs) {
                if (graph.get_path_handle_of_occurrence(occ) == p2 &&
                    graph.get_occurrence(occ) == graph.flip(h2)) {
                    found2 = true;
                }
                else {
                    REQUIRE(false);
                }
            }
            REQUIRE(found1);
            REQUIRE(found2);
            found1 = found2 = false;
            
        }
        
        vector<handle_t> segments = graph.divide_handle(h2, {size_t(2), size_t(4)});
        
        SECTION("PackedGraph preserves paths when dividing nodes") {
            check_path(p1, {h1, segments[0], segments[1], segments[2], h3});
            check_path(p2, {h1, graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0]), h3});
        }
        
        path_handle_t p3 = graph.create_path_handle("3");
        graph.append_occurrence(p3, h1);
        graph.append_occurrence(p3, segments[0]);
        
        REQUIRE(graph.has_path("3"));
        REQUIRE(graph.get_path_count() == 3);
        
        SECTION("PackedGraph can destroy paths") {
            graph.destroy_path(p3);
            
            REQUIRE(!graph.has_path("3"));
            REQUIRE(graph.get_path_count() == 2);
            
            bool found1 = false, found2 = false, found3 = false;
            
            graph.for_each_path_handle([&](const path_handle_t& p) {
                if (graph.get_path_name(p) == "1") {
                    found1 = true;
                }
                else if (graph.get_path_name(p) == "2") {
                    found2 = true;
                }
                else if (graph.get_path_name(p) == "3") {
                    found3 = true;
                }
                else {
                    REQUIRE(false);
                }
            });
            
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(!found3);
            
            
            // check flips to see if membership records are still functional
            check_flips(p1, {h1, segments[0], segments[1], segments[2], h3});
            check_flips(p2, {h1, graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0]), h3});
            
            graph.destroy_path(p1);
            
            REQUIRE(!graph.has_path("1"));
            REQUIRE(graph.get_path_count() == 1);
            
            found1 = found2 = found3 = false;
            
            graph.for_each_path_handle([&](const path_handle_t& p) {
                if (graph.get_path_name(p) == "1") {
                    found1 = true;
                }
                else if (graph.get_path_name(p) == "2") {
                    found2 = true;
                }
                else if (graph.get_path_name(p) == "3") {
                    found3 = true;
                }
                else {
                    REQUIRE(false);
                }
            });
            
            REQUIRE(!found1);
            REQUIRE(found2);
            REQUIRE(!found3);
            
            // check flips to see if membership records are still functional
            check_flips(p2, {h1, graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0]), h3});
        }
    }
    
    TEST_CASE("PackedGraph's defragmentation preserves topology", "[packed][handle]") {
        
        PackedGraph graph;
        
        auto check_path = [&](const path_handle_t& p, const vector<handle_t>& occs) {
            
            occurrence_handle_t occ;
            for (int i = 0; i < occs.size(); i++){
                if (i == 0) {
                    occ = graph.get_first_occurrence(p);
                }
                
                REQUIRE(graph.get_path_handle_of_occurrence(occ) == p);
                REQUIRE(graph.get_occurrence(occ) == occs[i]);
                REQUIRE(graph.has_previous_occurrence(occ) == (i > 0));
                REQUIRE(graph.has_next_occurrence(occ) == (i < occs.size() - 1));
                
                if (i != occs.size() - 1) {
                    occ = graph.get_next_occurrence(occ);
                }
            }
            
            for (int i = occs.size() - 1; i >= 0; i--){
                if (i == occs.size() - 1) {
                    occ = graph.get_last_occurrence(p);
                }
                
                REQUIRE(graph.get_path_handle_of_occurrence(occ) == p);
                REQUIRE(graph.get_occurrence(occ) == occs[i]);
                REQUIRE(graph.has_previous_occurrence(occ) == (i > 0));
                REQUIRE(graph.has_next_occurrence(occ) == (i < occs.size() - 1));
                
                if (i != 0) {
                    occ = graph.get_previous_occurrence(occ);
                }
            }
        };
        
        auto check_flips = [&](const path_handle_t& p, const vector<handle_t>& occs) {
            auto flipped = occs;
            for (size_t i = 0; i < occs.size(); i++) {
                
                graph.apply_orientation(graph.flip(graph.forward(flipped[i])));
                flipped[i] = graph.flip(flipped[i]);
                check_path(p, flipped);
                
                graph.apply_orientation(graph.flip(graph.forward(flipped[i])));
                flipped[i] = graph.flip(flipped[i]);
                check_path(p, flipped);
            }
        };
        
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("C");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("G");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h1, h3);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h5);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        
        path_handle_t p0 = graph.create_path_handle("0");
        path_handle_t p1 = graph.create_path_handle("1");
        path_handle_t p2 = graph.create_path_handle("2");
        
        
        graph.append_occurrence(p0, h3);
        graph.append_occurrence(p0, h4);
        graph.append_occurrence(p0, h5);
        
        graph.append_occurrence(p1, h1);
        graph.append_occurrence(p1, h3);
        graph.append_occurrence(p1, h5);
        
        graph.append_occurrence(p2, h1);
        graph.append_occurrence(p2, h2);
        graph.append_occurrence(p2, h3);
        graph.append_occurrence(p2, h4);
        graph.append_occurrence(p2, h5);
        
        
        // delete enough nodes/edges/path memberships to trigger defrag
        
        graph.destroy_path(p0);
        graph.destroy_path(p2);
        graph.destroy_handle(h2);
        graph.destroy_handle(h4);
        
        REQUIRE(graph.get_sequence(h1) == "A");
        REQUIRE(graph.get_sequence(h3) == "C");
        REQUIRE(graph.get_sequence(h5) == "G");
        
        bool found = false;
        graph.follow_edges(h1, false, [&](const handle_t& next) {
            if (next == h3) {
                found = true;
            }
            else {
                REQUIRE(false);
            }
            return true;
        });
        REQUIRE(found);
        
        found = false;
        graph.follow_edges(h3, false, [&](const handle_t& next) {
            if (next == h5) {
                found = true;
            }
            else {
                REQUIRE(false);
            }
            return true;
        });
        REQUIRE(found);
        
        check_flips(p1, {h1, h3, h5});
    }
}
}
        
