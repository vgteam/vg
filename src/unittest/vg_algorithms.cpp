//
//  vg_algorithms.cpp
//  
//
//  Created by Jordan Eizenga on 2/7/17.
//
//

#include <stdio.h>
#include <set>
#include <random>
#include "catch.hpp"
#include "algorithms/extract_connecting_graph.hpp"
#include "algorithms/extract_containing_graph.hpp"
#include "algorithms/extract_extending_graph.hpp"
#include "algorithms/topological_sort.hpp"
#include "algorithms/weakly_connected_components.hpp"
#include "algorithms/is_acyclic.hpp"
#include "algorithms/split_strands.hpp"
#include "algorithms/is_single_stranded.hpp"
#include "algorithms/distance_to_head.hpp"
#include "algorithms/distance_to_tail.hpp"
#include "algorithms/apply_bulk_modifications.hpp"
#include "algorithms/count_walks.hpp"
#include "algorithms/strongly_connected_components.hpp"
#include "algorithms/a_star.hpp"
#include "unittest/random_graph.hpp"
#include "vg.hpp"
#include "json2pb.h"


using namespace google::protobuf;

namespace vg {
    
    namespace unittest {
        
        Node* find_node(Graph& g, int64_t id) {
            for (int i = 0; i < g.node_size(); i++) {
                if (g.node(i).id() == id) {
                    return g.mutable_node(i);
                }
            }
            return nullptr;
        }
        
        TEST_CASE( "Connecting graph extraction produces expected results in an acyclic graph",
                  "[algorithms]" ) {
            
            VG vg;
            
            Node* n0 = vg.create_node("CGA");
            Node* n1 = vg.create_node("TTGG");
            Node* n2 = vg.create_node("CCGT");
            Node* n3 = vg.create_node("C");
            Node* n4 = vg.create_node("GT");
            Node* n5 = vg.create_node("GATAA");
            Node* n6 = vg.create_node("CGG");
            Node* n7 = vg.create_node("ACA");
            Node* n8 = vg.create_node("GCCG");
            Node* n9 = vg.create_node("ATATAAC");
            
            vg.create_edge(n1, n0, true, true); // a doubly reversing edge to keep it interesting
            vg.create_edge(n1, n2);
            vg.create_edge(n2, n3);
            vg.create_edge(n2, n4);
            vg.create_edge(n3, n5);
            vg.create_edge(n4, n5);
            vg.create_edge(n5, n6);
            vg.create_edge(n5, n8);
            vg.create_edge(n6, n7);
            vg.create_edge(n6, n8);
            vg.create_edge(n7, n9);
            vg.create_edge(n8, n9);
            
            SECTION( "Graph extraction can extract a section of one node" ) {
                
                pos_t pos_1 = make_pos_t(n1->id(), false, 1);
                pos_t pos_2 = make_pos_t(n1->id(), false, 3);
                
                int64_t max_len = 5;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 1 );
                REQUIRE( g.edge_size() == 0 );
                REQUIRE( trans[g.node(0).id()] == n1->id() );
                REQUIRE( g.node(0).sequence() == "TG" );
            }
            
            SECTION( "Graph extraction can extract a section between two nodes" ) {
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 2);
                pos_t pos_2 = make_pos_t(n1->id(), false, 2);
                
                int64_t max_len = 5;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n0->id(), n1->id()};
                
                REQUIRE( g.node_size() == 2 );
                REQUIRE( g.edge_size() == 1 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n0->id()) {
                        REQUIRE( n.sequence() == "A" );
                    }
                    else {
                        REQUIRE( n.sequence() == "TT" );
                    }
                }
                
                const Edge& e = g.edge(0);
                REQUIRE(expected_node_ids.count(trans[e.from()]));
                REQUIRE(expected_node_ids.count(trans[e.to()]));
                REQUIRE(trans[e.from()] != trans[e.to()]);
                REQUIRE(e.from_start() == e.to_end());
            }
            
            SECTION( "Graph extraction can extract a bubble" ) {
                
                pos_t pos_1 = make_pos_t(n2->id(), false, 2);
                pos_t pos_2 = make_pos_t(n5->id(), false, 3);
                
                int64_t max_len = 10;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n2->id(), n3->id(), n4->id(), n5->id()};
                
                REQUIRE( g.node_size() == 4 );
                REQUIRE( g.edge_size() == 4 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n2->id()) {
                        REQUIRE( n.sequence() == "GT" );
                    }
                    else if (trans[n.id()] == n3->id()) {
                        REQUIRE( n.sequence() == "C" );
                    }
                    else if (trans[n.id()] == n4->id()) {
                        REQUIRE( n.sequence() == "GT" );
                    }
                    else if (trans[n.id()] == n5->id()) {
                        REQUIRE( n.sequence() == "GAT" );
                    }
                }
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    
                    const Edge& e = g.edge(i);
                    REQUIRE(expected_node_ids.count(trans[e.from()]));
                    REQUIRE(expected_node_ids.count(trans[e.to()]));
                    REQUIRE(e.from_start() == e.to_end());
                    
                    if ((trans[e.from()] == n2->id() && trans[e.to()] == n3->id())
                        || (trans[e.from()] == n3->id() && trans[e.to()] == n2->id()) ) {
                        found_edge_0 = true;
                    }
                    else if ((trans[e.from()] == n2->id() && trans[e.to()] == n4->id())
                             || (trans[e.from()] == n4->id() && trans[e.to()] == n2->id()) ) {
                        found_edge_1 = true;
                    }
                    else if ((trans[e.from()] == n5->id() && trans[e.to()] == n3->id())
                             || (trans[e.from()] == n3->id() && trans[e.to()] == n5->id()) ) {
                        found_edge_2 = true;
                    }
                    else if ((trans[e.from()] == n5->id() && trans[e.to()] == n4->id())
                             || (trans[e.from()] == n4->id() && trans[e.to()] == n5->id()) ) {
                        found_edge_3 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
            }
            
            SECTION( "Graph extraction can trim nodes that are not a path between the two positions" ) {
                
                pos_t pos_1 = make_pos_t(n2->id(), false, 1);
                pos_t pos_2 = make_pos_t(n4->id(), false, 1);
                
                int64_t max_len = 10;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n2->id(), n4->id()};
                
                REQUIRE( g.node_size() == 2 );
                REQUIRE( g.edge_size() == 1 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n2->id()) {
                        REQUIRE( n.sequence() == "CGT" );
                    }
                    else {
                        REQUIRE( n.sequence() == "G" );
                    }
                }
                
                const Edge& e = g.edge(0);
                REQUIRE(expected_node_ids.count(trans[e.from()]));
                REQUIRE(expected_node_ids.count(trans[e.to()]));
                REQUIRE(trans[e.from()] != trans[e.to()]);
                REQUIRE(e.from_start() == e.to_end());
            }
            
            SECTION( "Graph extraction returns an empty graph if there are no paths under the max length and the positions are on the same node" ) {
                
                pos_t pos_1 = make_pos_t(n5->id(), true, 1);
                pos_t pos_2 = make_pos_t(n5->id(), true, 4);
                
                int64_t max_len = 2;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 0 );
                REQUIRE( g.edge_size() == 0 );
            }
            
            SECTION( "Graph extraction returns an empty graph if there are no paths under the max length and the positions are on different nodes" ) {
                
                pos_t pos_1 = make_pos_t(n5->id(), true, 1);
                pos_t pos_2 = make_pos_t(n0->id(), true, 2);
                
                int64_t max_len = 7;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 0 );
                REQUIRE( g.edge_size() == 0 );
            }
            
            SECTION( "Graph extraction can detect a one-node path at the maximum distance" ) {
                
                pos_t pos_1 = make_pos_t(n5->id(), true, 1);
                pos_t pos_2 = make_pos_t(n5->id(), true, 4);
                
                int64_t max_len = 3;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 1 );
                REQUIRE( g.edge_size() == 0 );
                REQUIRE( trans[g.node(0).id()] == n5->id() );
                REQUIRE( g.node(0).sequence() == "ATA" );
            }
            
            SECTION( "Graph extraction can detect a two-node path at the maximum distance" ) {
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 1);
                pos_t pos_2 = make_pos_t(n1->id(), false, 3);
                
                int64_t max_len = 5;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n0->id(), n1->id()};
                
                REQUIRE( g.node_size() == 2 );
                REQUIRE( g.edge_size() == 1 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n0->id()) {
                        REQUIRE( n.sequence() == "GA" );
                    }
                    else {
                        REQUIRE( n.sequence() == "TTG" );
                    }
                }
                
                const Edge& e = g.edge(0);
                REQUIRE(expected_node_ids.count(trans[e.from()]));
                REQUIRE(expected_node_ids.count(trans[e.to()]));
                REQUIRE(trans[e.from()] != trans[e.to()]);
                REQUIRE(e.from_start() == e.to_end());
            }
            
            SECTION( "Graph extraction can detect a multi-node path at the maximum distance" ) {
                
                pos_t pos_1 = make_pos_t(n2->id(), false, 1);
                pos_t pos_2 = make_pos_t(n5->id(), false, 4);
                
                int64_t max_len = 8;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n2->id(), n3->id(), n5->id()};
                
                REQUIRE( g.node_size() == 3 );
                REQUIRE( g.edge_size() == 2 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n2->id()) {
                        REQUIRE( n.sequence() == "CGT" );
                    }
                    else if (trans[n.id()] == n3->id()) {
                        REQUIRE( n.sequence() == "C" );
                    }
                    else {
                        REQUIRE( n.sequence() == "GATA" );
                    }
                }
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    
                    const Edge& e = g.edge(i);
                    REQUIRE(expected_node_ids.count(trans[e.from()]));
                    REQUIRE(expected_node_ids.count(trans[e.to()]));
                    REQUIRE(e.from_start() == e.to_end());
                    
                    if ((trans[e.from()] == n2->id() && trans[e.to()] == n3->id())
                        || (trans[e.from()] == n3->id() && trans[e.to()] == n2->id()) ) {
                        found_edge_0 = true;
                    }
                    else if ((trans[e.from()] == n3->id() && trans[e.to()] == n5->id())
                             || (trans[e.from()] == n5->id() && trans[e.to()] == n3->id()) ) {
                        found_edge_1 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
            }
            
            SECTION( "Graph extraction can remove paths that are above maximum length" ) {
                
                pos_t pos_1 = make_pos_t(n9->id(), true, 7);
                pos_t pos_2 = make_pos_t(n5->id(), true, 0);
                
                int64_t max_len = 4;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n5->id(), n8->id(), n9->id()};
                
                REQUIRE( g.node_size() == 3 );
                REQUIRE( g.edge_size() == 2 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n5->id()) {
                        REQUIRE( n.sequence() == "" );
                    }
                    else if (trans[n.id()] == n8->id()) {
                        REQUIRE( n.sequence() == "GCCG" );
                    }
                    else {
                        REQUIRE( n.sequence() == "" );
                    }
                }
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    
                    const Edge& e = g.edge(i);
                    REQUIRE(expected_node_ids.count(trans[e.from()]));
                    REQUIRE(expected_node_ids.count(trans[e.to()]));
                    REQUIRE(e.from_start() == e.to_end());
                    
                    if ((trans[e.from()] == n5->id() && trans[e.to()] == n8->id())
                        || (trans[e.from()] == n8->id() && trans[e.to()] == n5->id()) ) {
                        found_edge_0 = true;
                    }
                    else if ((trans[e.from()] == n9->id() && trans[e.to()] == n8->id())
                             || (trans[e.from()] == n8->id() && trans[e.to()] == n9->id()) ) {
                        found_edge_1 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
            }
        }
        
        TEST_CASE( "Connecting graph extraction produces expected results in cyclic graphs",
                  "[algorithms]" ) {
            
            SECTION( "Graph extraction can detect a loop onto the same node" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("ATATAAC");
                
                vg.create_edge(n0, n0);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 5);
                pos_t pos_2 = make_pos_t(n0->id(), false, 3);
                
                int64_t max_len = 6;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n0->id()};
                
                bool found_node_copy_0 = false;
                bool found_node_copy_1 = false;
                
                REQUIRE( g.node_size() == 2 );
                REQUIRE( g.edge_size() == 1 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n0->id()) {
                        if (n.sequence() == "AC") {
                            found_node_copy_0 = true;
                        }
                        else if (n.sequence() == "ATA") {
                            found_node_copy_1 = true;
                        }
                    }
                }
                
                REQUIRE(found_node_copy_0);
                REQUIRE(found_node_copy_1);
                
                const Edge& e = g.edge(0);
                REQUIRE(expected_node_ids.count(trans[e.from()]));
                REQUIRE(expected_node_ids.count(trans[e.to()]));
                REQUIRE(e.from_start() == e.to_end());
            }
            
            SECTION( "Graph extraction can detect a self-loop onto the same node in the opposite orientation" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("ATATAAC");
                
                vg.create_edge(n0, n0, false, true);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 5);
                pos_t pos_2 = make_pos_t(n0->id(), true, 3);
                
                int64_t max_len = 10;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, false, true, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n0->id()};
                
                bool found_node_copy_0 = false;
                bool found_node_copy_1 = false;
                
                REQUIRE( g.node_size() == 2 );
                REQUIRE( g.edge_size() == 1 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n0->id()) {
                        if (n.sequence() == "AC") {
                            found_node_copy_0 = true;
                        }
                        else if (n.sequence() == "AAC") {
                            found_node_copy_1 = true;
                        }
                    }
                }
                
                REQUIRE(found_node_copy_0);
                REQUIRE(found_node_copy_1);
                
                const Edge& e = g.edge(0);
                REQUIRE(expected_node_ids.count(trans[e.from()]));
                REQUIRE(expected_node_ids.count(trans[e.to()]));
                REQUIRE(e.from_start() != e.to_end());
            }
            
            SECTION( "Graph extraction can detect a self-loop to the start node" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("TGAT");
                Node* n1 = vg.create_node("AC");
                
                vg.create_edge(n0, n0);
                vg.create_edge(n0, n1);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 3);
                pos_t pos_2 = make_pos_t(n1->id(), false, 1);
                
                int64_t max_len = 15;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, true);
                
                Graph& g = extractor.graph;
                
                set<int64_t> expected_node_ids{n0->id(), n1->id()};
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                
                REQUIRE( g.node_size() == 3 );
                REQUIRE( g.edge_size() == 4 );
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    REQUIRE( expected_node_ids.count(trans[n.id()]) );
                    if (trans[n.id()] == n0->id()) {
                        if (n.sequence() == "T") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "TGAT") {
                            found_node_1 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[n.id()] == n1->id()) {
                        if (n.sequence() == "A") {
                            found_node_2 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    REQUIRE(expected_node_ids.count(trans[e.from()]));
                    REQUIRE(expected_node_ids.count(trans[e.to()]));
                    
                    if (trans[e.from()] == n0->id() && trans[e.to()] == n0->id() ) {
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        if ((from_seq == "T" && to_seq == "TGAT" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "T" && from_seq == "TGAT" && e.from_start() && e.to_end())) {
                            found_edge_0 = true;
                        }
                        else if ((from_seq == "TGAT" && to_seq == "TGAT" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TGAT" && from_seq == "TGAT" && e.from_start() && e.to_end())) {
                            found_edge_1 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if ((trans[e.from()] == n0->id() && trans[e.to()] == n1->id())
                             || (trans[e.from()] == n1->id() && trans[e.to()] == n0->id()) ) {
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        if ((from_seq == "T" && to_seq == "A" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "T" && from_seq == "A" && e.from_start() && e.to_end())) {
                            found_edge_2 = true;
                        }
                        else if ((from_seq == "TGAT" && to_seq == "A" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TGAT" && from_seq == "A" && e.from_start() && e.to_end())) {
                            found_edge_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                
            }
            
            SECTION( "Graph extraction can detect simultaneous self-loops on start and end nodes" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("TGAT");
                Node* n1 = vg.create_node("AC");
                
                vg.create_edge(n0, n0);
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n1);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 3);
                pos_t pos_2 = make_pos_t(n1->id(), false, 1);
                
                int64_t max_len = 20;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 4 );
                REQUIRE( g.edge_size() == 8 );
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (trans[n.id()] == n0->id()) {
                        if (n.sequence() == "T") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "TGAT") {
                            found_node_1 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[n.id()] == n1->id()) {
                        if (n.sequence() == "A") {
                            found_node_2 = true;
                        }
                        else if (n.sequence() == "AC") {
                            found_node_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                bool found_edge_4 = false;
                bool found_edge_5 = false;
                bool found_edge_6 = false;
                bool found_edge_7 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (trans[e.from()] == n0->id() && trans[e.to()] == n0->id() ) {
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        if ((from_seq == "T" && to_seq == "TGAT" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "T" && from_seq == "TGAT" && e.from_start() && e.to_end())) {
                            found_edge_0 = true;
                        }
                        else if ((from_seq == "TGAT" && to_seq == "TGAT" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TGAT" && from_seq == "TGAT" && e.from_start() && e.to_end())) {
                            found_edge_1 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[e.from()] == n1->id() && trans[e.to()] == n1->id() ) {
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        if ((from_seq == "AC" && to_seq == "A" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "AC" && from_seq == "A" && e.from_start() && e.to_end())) {
                            found_edge_2 = true;
                        }
                        else if ((from_seq == "AC" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "AC" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if ((trans[e.from()] == n0->id() && trans[e.to()] == n1->id())
                             || (trans[e.from()] == n1->id() && trans[e.to()] == n0->id()) ) {
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        if ((from_seq == "T" && to_seq == "A" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "T" && from_seq == "A" && e.from_start() && e.to_end())) {
                            found_edge_4 = true;
                        }
                        else if ((from_seq == "TGAT" && to_seq == "A" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TGAT" && from_seq == "A" && e.from_start() && e.to_end())) {
                            found_edge_5 = true;
                        }
                        else if ((from_seq == "T" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "T" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_6 = true;
                        }
                        else if ((from_seq == "TGAT" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TGAT" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_7 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                REQUIRE(found_edge_4);
                REQUIRE(found_edge_5);
                REQUIRE(found_edge_6);
                REQUIRE(found_edge_7);
                
            }
            
            SECTION( "Graph extraction can detect simultaneous non-self loops on start and end nodes" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("TGAT");
                Node* n1 = vg.create_node("AC");
                Node* n2 = vg.create_node("GT");
                Node* n3 = vg.create_node("CG");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n0, n2);
                vg.create_edge(n1, n3);
                vg.create_edge(n2, n0);
                vg.create_edge(n3, n1);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 3);
                pos_t pos_2 = make_pos_t(n1->id(), false, 1);
                
                int64_t max_len = 30;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 6 );
                REQUIRE( g.edge_size() == 10 );
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                bool found_node_4 = false;
                bool found_node_5 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (trans[n.id()] == n0->id()) {
                        if (n.sequence() == "T") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "TGAT") {
                            found_node_1 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[n.id()] == n1->id()) {
                        if (n.sequence() == "A") {
                            found_node_2 = true;
                        }
                        else if (n.sequence() == "AC") {
                            found_node_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[n.id()] == n2->id()) {
                        found_node_4 = true;
                    }
                    else if (trans[n.id()] == n3->id()) {
                        found_node_5 = true;
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                REQUIRE(found_node_4);
                REQUIRE(found_node_5);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                bool found_edge_4 = false;
                bool found_edge_5 = false;
                bool found_edge_6 = false;
                bool found_edge_7 = false;
                bool found_edge_8 = false;
                bool found_edge_9 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if ((trans[e.from()] == n0->id() && trans[e.to()] == n1->id())
                        || (trans[e.from()] == n1->id() && trans[e.to()] == n0->id())) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "T" && to_seq == "A" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "T" && from_seq == "A" && e.from_start() && e.to_end())) {
                            found_edge_0 = true;
                        }
                        else if ((from_seq == "TGAT" && to_seq == "A" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TGAT" && from_seq == "A" && e.from_start() && e.to_end())) {
                            found_edge_1 = true;
                        }
                        else if ((from_seq == "T" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "T" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_2 = true;
                        }
                        else if ((from_seq == "TGAT" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TGAT" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                        
                    }
                    else if ((trans[e.from()] == n0->id() && trans[e.to()] == n2->id())
                             || (trans[e.from()] == n2->id() && trans[e.to()] == n0->id()) ) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "T" && to_seq == "GT" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "T" && from_seq == "GT" && e.from_start() && e.to_end())) {
                            found_edge_4 = true;
                        }
                        else if ((from_seq == "TGAT" && to_seq == "GT" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TGAT" && from_seq == "GT" && e.from_start() && e.to_end())) {
                            found_edge_5 = true;
                        }
                        else if ((from_seq == "GT" && to_seq == "TGAT" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "GT" && from_seq == "TGAT" && e.from_start() && e.to_end())) {
                            found_edge_6 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if ((trans[e.from()] == n1->id() && trans[e.to()] == n3->id())
                             || (trans[e.from()] == n3->id() && trans[e.to()] == n1->id()) ) {
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "CG" && to_seq == "A" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "CG" && from_seq == "A" && e.from_start() && e.to_end())) {
                            found_edge_7 = true;
                        }
                        else if ((from_seq == "CG" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "CG" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_8 = true;
                        }
                        else if ((from_seq == "AC" && to_seq == "CG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "AC" && from_seq == "CG" && e.from_start() && e.to_end())) {
                            found_edge_9    = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }

                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                REQUIRE(found_edge_4);
                REQUIRE(found_edge_5);
                REQUIRE(found_edge_6);
                REQUIRE(found_edge_7);
                REQUIRE(found_edge_8);
                REQUIRE(found_edge_9);
                
            }
            
            SECTION( "Graph extraction can detect loops on positions on the same node that can reach each other" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("TCGAT");
                Node* n1 = vg.create_node("AC");
                
                vg.create_edge(n0, n0);
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n0);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 2);
                pos_t pos_2 = make_pos_t(n0->id(), false, 3);
                
                int64_t max_len = 40;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 5 );
                REQUIRE( g.edge_size() == 8 );
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                bool found_node_4 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (trans[n.id()] == n0->id()) {
                        if (n.sequence() == "G") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "TCGAT") {
                            found_node_1 = true;
                        }
                        else if (n.sequence() == "GAT") {
                            found_node_2 = true;
                        }
                        else if (n.sequence() == "TCG") {
                            found_node_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[n.id()] == n1->id()) {
                        if (n.sequence() == "AC") {
                            found_node_4 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                REQUIRE(found_node_4);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                bool found_edge_4 = false;
                bool found_edge_5 = false;
                bool found_edge_6 = false;
                bool found_edge_7 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if ((trans[e.from()] == n0->id() && trans[e.to()] == n0->id())
                        || (trans[e.from()] == n0->id() && trans[e.to()] == n0->id())) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "GAT" && to_seq == "TCG" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "GAT" && from_seq == "TCG" && e.from_start() && e.to_end())) {
                            found_edge_0 = true;
                        }
                        else if ((from_seq == "GAT" && to_seq == "TCGAT" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "GAT" && from_seq == "TCGAT" && e.from_start() && e.to_end())) {
                            found_edge_1 = true;
                        }
                        else if ((from_seq == "TCGAT" && to_seq == "TCG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TCGAT" && from_seq == "TCG" && e.from_start() && e.to_end())) {
                            found_edge_2 = true;
                        }
                        else if ((from_seq == "TCGAT" && to_seq == "TCGAT" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TCGAT" && from_seq == "TCGAT" && e.from_start() && e.to_end())) {
                            found_edge_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                        
                    }
                    else if ((trans[e.from()] == n0->id() && trans[e.to()] == n1->id())
                             || (trans[e.from()] == n1->id() && trans[e.to()] == n0->id()) ) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "GAT" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "GAT" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_4 = true;
                        }
                        else if ((from_seq == "TCGAT" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TCGAT" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_5 = true;
                        }
                        else if ((from_seq == "AC" && to_seq == "TCGAT" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "AC" && from_seq == "TCGAT" && e.from_start() && e.to_end())) {
                            found_edge_6 = true;
                        }
                        else if ((from_seq == "AC" && to_seq == "TCG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "AC" && from_seq == "TCG" && e.from_start() && e.to_end())) {
                            found_edge_7 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                REQUIRE(found_edge_4);
                REQUIRE(found_edge_5);
                REQUIRE(found_edge_6);
                REQUIRE(found_edge_7);
                
            }
            
            SECTION( "Graph extraction can detect loops on positions on the same node that cannot reach each other" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("TCGAG");
                Node* n1 = vg.create_node("AC");
                
                vg.create_edge(n0, n0);
                vg.create_edge(n0, n1, false, true);
                vg.create_edge(n1, n0, true, false);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 4);
                pos_t pos_2 = make_pos_t(n0->id(), false, 1);
                
                int64_t max_len = 40;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 4 );
                REQUIRE( g.edge_size() == 8 );
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (trans[n.id()] == n0->id()) {
                        if (n.sequence() == "G") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "T") {
                            found_node_1 = true;
                        }
                        else if (n.sequence() == "TCGAG") {
                            found_node_2 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[n.id()] == n1->id()) {
                        if (n.sequence() == "AC") {
                            found_node_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                bool found_edge_4 = false;
                bool found_edge_5 = false;
                bool found_edge_6 = false;
                bool found_edge_7 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if ((trans[e.from()] == n0->id() && trans[e.to()] == n0->id())
                        || (trans[e.from()] == n0->id() && trans[e.to()] == n0->id())) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "G" && to_seq == "T" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "G" && from_seq == "T" && e.from_start() && e.to_end())) {
                            found_edge_0 = true;
                        }
                        else if ((from_seq == "G" && to_seq == "TCGAG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "G" && from_seq == "TCGAG" && e.from_start() && e.to_end())) {
                            found_edge_1 = true;
                        }
                        else if ((from_seq == "TCGAG" && to_seq == "T" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TCGAG" && from_seq == "T" && e.from_start() && e.to_end())) {
                            found_edge_2 = true;
                        }
                        else if ((from_seq == "TCGAG" && to_seq == "TCGAG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TCGAG" && from_seq == "TCGAG" && e.from_start() && e.to_end())) {
                            found_edge_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                        
                    }
                    else if ((trans[e.from()] == n0->id() && trans[e.to()] == n1->id())
                             || (trans[e.from()] == n1->id() && trans[e.to()] == n0->id()) ) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "G" && to_seq == "AC" && !e.from_start() && e.to_end()) ||
                            (to_seq == "G" && from_seq == "AC" && !e.from_start() && e.to_end())) {
                            found_edge_4 = true;
                        }
                        else if ((from_seq == "TCGAG" && to_seq == "AC" && !e.from_start() && e.to_end()) ||
                                 (to_seq == "TCGAG" && from_seq == "AC" && !e.from_start() && e.to_end())) {
                            found_edge_5 = true;
                        }
                        else if ((from_seq == "AC" && to_seq == "TCGAG" && e.from_start() && !e.to_end()) ||
                                 (to_seq == "AC" && from_seq == "TCGAG" && e.from_start() && !e.to_end())) {
                            found_edge_6 = true;
                        }
                        else if ((from_seq == "AC" && to_seq == "T" && e.from_start() && !e.to_end()) ||
                                 (to_seq == "AC" && from_seq == "T" && e.from_start() && !e.to_end())) {
                            found_edge_7 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                REQUIRE(found_edge_4);
                REQUIRE(found_edge_5);
                REQUIRE(found_edge_6);
                REQUIRE(found_edge_7);
                
            }
            
            SECTION( "Graph extraction can detect loops on positions on the same node on the opposite strand" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("TCGAG");
                Node* n1 = vg.create_node("AC");
                Node* n2 = vg.create_node("GTA");
                
                vg.create_edge(n0, n0);
                vg.create_edge(n0, n0, true, false);
                vg.create_edge(n0, n0, false, true);
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n0);
                vg.create_edge(n0, n2);
                vg.create_edge(n2, n0, false, true);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 4);
                pos_t pos_2 = make_pos_t(n0->id(), true, 3);
                
                int64_t max_len = 40;
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, max_len, pos_1, pos_2, true, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE( g.node_size() == 5 );
                REQUIRE( g.edge_size() == 18 );
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                bool found_node_4 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (trans[n.id()] == n0->id()) {
                        if (n.sequence() == "G") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "GAG") {
                            found_node_1 = true;
                        }
                        else if (n.sequence() == "TCGAG") {
                            found_node_2 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[n.id()] == n1->id()) {
                        if (n.sequence() == "AC") {
                            found_node_3 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if (trans[n.id()] == n2->id()) {
                        if (n.sequence() == "GTA") {
                            found_node_4 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                REQUIRE(found_node_4);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                bool found_edge_4 = false;
                bool found_edge_5 = false;
                bool found_edge_6 = false;
                bool found_edge_7 = false;
                bool found_edge_8 = false;
                bool found_edge_9 = false;
                bool found_edge_10 = false;
                bool found_edge_11 = false;
                bool found_edge_12 = false;
                bool found_edge_13 = false;
                bool found_edge_14 = false;
                bool found_edge_15 = false;
                bool found_edge_16 = false;
                bool found_edge_17 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if ((trans[e.from()] == n0->id() && trans[e.to()] == n0->id())
                        || (trans[e.from()] == n0->id() && trans[e.to()] == n0->id())) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "G" && to_seq == "GAG" && !e.from_start() && e.to_end()) ||
                            (to_seq == "G" && from_seq == "GAG" && !e.from_start() && e.to_end())) {
                            found_edge_0 = true;
                        }
                        else if ((from_seq == "G" && to_seq == "TCGAG" && !e.from_start() && e.to_end()) ||
                                 (to_seq == "G" && from_seq == "TCGAG" && !e.from_start() && e.to_end())) {
                            found_edge_1 = true;
                        }
                        else if ((from_seq == "G" && to_seq == "TCGAG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "G" && from_seq == "TCGAG" && e.from_start() && e.to_end())) {
                            found_edge_2 = true;
                        }
                        else if ((from_seq == "GAG" && to_seq == "TCGAG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "GAG" && from_seq == "TCGAG" && e.from_start() && e.to_end())) {
                            found_edge_3 = true;
                        }
                        else if ((from_seq == "GAG" && to_seq == "TCGAG" && !e.from_start() && e.to_end()) ||
                                 (to_seq == "GAG" && from_seq == "TCGAG" && !e.from_start() && e.to_end())) {
                            found_edge_4 = true;
                        }
                        else if ((from_seq == "TCGAG" && to_seq == "TCGAG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TCGAG" && from_seq == "TCGAG" && e.from_start() && e.to_end())) {
                            found_edge_5 = true;
                        }
                        else if ((from_seq == "TCGAG" && to_seq == "TCGAG" && !e.from_start() && e.to_end()) ||
                                 (to_seq == "TCGAG" && from_seq == "TCGAG" && !e.from_start() && e.to_end())) {
                            found_edge_6 = true;
                        }
                        else if ((from_seq == "TCGAG" && to_seq == "TCGAG" && e.from_start() && !e.to_end()) ||
                                 (to_seq == "TCGAG" && from_seq == "TCGAG" && e.from_start() && !e.to_end())) {
                            found_edge_7 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                        
                    }
                    else if ((trans[e.from()] == n0->id() && trans[e.to()] == n1->id())
                             || (trans[e.from()] == n1->id() && trans[e.to()] == n0->id()) ) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "G" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "G" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_8 = true;
                        }
                        else if ((from_seq == "GAG" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "GAG" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_9 = true;
                        }
                        else if ((from_seq == "TCGAG" && to_seq == "AC" && !e.from_start() && !e.to_end()) ||
                             (to_seq == "TCGAG" && from_seq == "AC" && e.from_start() && e.to_end())) {
                            found_edge_10 = true;
                        }
                        else if ((from_seq == "AC" && to_seq == "TCGAG" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "AC" && from_seq == "TCGAG" && e.from_start() && e.to_end())) {
                            found_edge_11 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else if ((trans[e.from()] == n0->id() && trans[e.to()] == n2->id())
                             || (trans[e.from()] == n2->id() && trans[e.to()] == n0->id()) ) {
                        
                        string from_seq = find_node(g, e.from())->sequence();
                        string to_seq = find_node(g, e.to())->sequence();
                        
                        if ((from_seq == "G" && to_seq == "GTA" && !e.from_start() && !e.to_end()) ||
                            (to_seq == "G" && from_seq == "GTA" && e.from_start() && e.to_end())) {
                            found_edge_12 = true;
                        }
                        else if ((from_seq == "GAG" && to_seq == "GTA" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "GAG" && from_seq == "GTA" && e.from_start() && e.to_end())) {
                            found_edge_13 = true;
                        }
                        else if ((from_seq == "TCGAG" && to_seq == "GTA" && !e.from_start() && !e.to_end()) ||
                                 (to_seq == "TCGAG" && from_seq == "GTA" && e.from_start() && e.to_end())) {
                            found_edge_14 = true;
                        }
                        else if ((from_seq == "GTA" && to_seq == "G" && !e.from_start() && e.to_end()) ||
                                 (to_seq == "GTA" && from_seq == "G" && !e.from_start() && e.to_end())) {
                            found_edge_15 = true;
                        }
                        else if ((from_seq == "GTA" && to_seq == "GAG" && !e.from_start() && e.to_end()) ||
                                 (to_seq == "GTA" && from_seq == "GAG" && !e.from_start() && e.to_end())) {
                            found_edge_16 = true;
                        }
                        else if ((from_seq == "GTA" && to_seq == "TCGAG" && !e.from_start() && e.to_end()) ||
                                 (to_seq == "GTA" && from_seq == "TCGAG" && !e.from_start() && e.to_end())) {
                            found_edge_17 = true;
                        }
                        else {
                            REQUIRE(false);
                        }
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                REQUIRE(found_edge_4);
                REQUIRE(found_edge_5);
                REQUIRE(found_edge_6);
                REQUIRE(found_edge_7);
                REQUIRE(found_edge_8);
                REQUIRE(found_edge_9);
                REQUIRE(found_edge_10);
                REQUIRE(found_edge_11);
                REQUIRE(found_edge_12);
                REQUIRE(found_edge_13);
                REQUIRE(found_edge_14);
                REQUIRE(found_edge_15);
                REQUIRE(found_edge_16);
                REQUIRE(found_edge_17);
                
            }
        }
        
        TEST_CASE( "Connecting graph extraction pruning options perform as expected", "[algorithms][broken]" ) {
            VG vg;
            
            Node* n0 = vg.create_node("CGA");
            Node* n1 = vg.create_node("TTGG");
            Node* n2 = vg.create_node("GT"); // a tip
            Node* n3 = vg.create_node("ATG"); // a self-contained cycle
            Node* n4 = vg.create_node("TGAG"); // a longer path
            Node* n5 = vg.create_node("CA"); // a shorter path
            Node* n6 = vg.create_node("T");
            
            vg.create_edge(n1, n0, true, true);
            vg.create_edge(n1, n2);
            vg.create_edge(n1, n3);
            vg.create_edge(n1, n4, false, true);
            vg.create_edge(n1, n5);
            vg.create_edge(n3, n3, true, true);
            vg.create_edge(n4, n6, true, false);
            vg.create_edge(n5, n6);
            
            SECTION("Connecting graph extraction can extract a graph without pruning") {
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, 10,
                                                                  make_pos_t(n1->id(), false, 3),
                                                                  make_pos_t(n6->id(), false, 1));
                
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 6);
                REQUIRE(g.edge_size() == 7);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                bool found_node_4 = false;
                bool found_node_5 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (n.sequence() == "G") {
                        found_node_0 = true;
                    }
                    else if (n.sequence() == "T") {
                        found_node_1 = true;
                    }
                    else if (n.sequence() == "GT") {
                        found_node_2 = true;
                    }
                    else if (n.sequence() == "ATG") {
                        found_node_3 = true;
                    }
                    else if (n.sequence() == "TGAG") {
                        found_node_4 = true;
                    }
                    else if (n.sequence() == "CA") {
                        found_node_5 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                REQUIRE(found_node_4);
                REQUIRE(found_node_5);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                bool found_edge_4 = false;
                bool found_edge_5 = false;
                bool found_edge_6 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((trans[e.from()] == n1->id() && !e.from_start() && trans[e.to()] == n2->id() && !e.to_end())
                         ||(trans[e.from()] == n2->id() && e.from_start() && trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((trans[e.from()] == n1->id() && !e.from_start() && trans[e.to()] == n3->id() && !e.to_end())
                              ||(trans[e.from()] == n3->id() && e.from_start() && trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_1 = true;
                    }
                    else if (((trans[e.from()] == n1->id() && !e.from_start() && trans[e.to()] == n4->id() && e.to_end())
                              ||(trans[e.from()] == n4->id() && !e.from_start() && trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_2 = true;
                    }
                    else if (((trans[e.from()] == n1->id() && !e.from_start() && trans[e.to()] == n5->id() && !e.to_end())
                              ||(trans[e.from()] == n5->id() && e.from_start() && trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_3 = true;
                    }
                    else if (((trans[e.from()] == n3->id() && !e.from_start() && trans[e.to()] == n3->id() && !e.to_end())
                              ||(trans[e.from()] == n3->id() && e.from_start() && trans[e.to()] == n3->id() && e.to_end()))) {
                        found_edge_4 = true;
                    }
                    else if (((trans[e.from()] == n4->id() && e.from_start() && trans[e.to()] == n6->id() && !e.to_end())
                              ||(trans[e.from()] == n6->id() && e.from_start() && trans[e.to()] == n4->id() && !e.to_end()))) {
                        found_edge_5 = true;
                    }
                    else if (((trans[e.from()] == n5->id() && !e.from_start() && trans[e.to()] == n6->id() && !e.to_end())
                              ||(trans[e.from()] == n6->id() && e.from_start() && trans[e.to()] == n5->id() && e.to_end()))) {
                        found_edge_6 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                REQUIRE(found_edge_4);
                REQUIRE(found_edge_5);
                REQUIRE(found_edge_6);
            }
            
            SECTION("Connecting graph extraction can prune a graph to only nodes on walks between the two positions") {
                
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, 10,
                                                                  make_pos_t(n1->id(), false, 3),
                                                                  make_pos_t(n6->id(), false, 1),
                                                                  false, true);
                
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 4);
                REQUIRE(g.edge_size() == 4);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (n.sequence() == "G") {
                        found_node_0 = true;
                    }
                    else if (n.sequence() == "T") {
                        found_node_1 = true;
                    }
                    else if (n.sequence() == "TGAG") {
                        found_node_2 = true;
                    }
                    else if (n.sequence() == "CA") {
                        found_node_3 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((trans[e.from()] == n1->id() && !e.from_start() && trans[e.to()] == n4->id() && e.to_end())
                         ||(trans[e.from()] == n4->id() && !e.from_start() && trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((trans[e.from()] == n1->id() && !e.from_start() && trans[e.to()] == n5->id() && !e.to_end())
                              ||(trans[e.from()] == n5->id() && e.from_start() && trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_1 = true;
                    }
                    else if (((trans[e.from()] == n4->id() && e.from_start() && trans[e.to()] == n6->id() && !e.to_end())
                              ||(trans[e.from()] == n6->id() && e.from_start() && trans[e.to()] == n4->id() && !e.to_end()))) {
                        found_edge_2 = true;
                    }
                    else if (((trans[e.from()] == n5->id() && !e.from_start() && trans[e.to()] == n6->id() && !e.to_end())
                              ||(trans[e.from()] == n6->id() && e.from_start() && trans[e.to()] == n5->id() && e.to_end()))) {
                        found_edge_3 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
            }
            
            SECTION("Connecting graph extraction can prune a graph to only nodes on sufficiently short walks between the two positions") {
                VG extractor;
                auto trans = algorithms::extract_connecting_graph(&vg, &extractor, 4,
                                                                  make_pos_t(n1->id(), false, 3),
                                                                  make_pos_t(n6->id(), false, 1),
                                                                  false, true, true);
                Graph g = extractor.graph;
                
                REQUIRE(g.node_size() == 3);
                REQUIRE(g.edge_size() == 2);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (n.sequence() == "G") {
                        found_node_0 = true;
                    }
                    else if (n.sequence() == "T") {
                        found_node_1 = true;
                    }
                    else if (n.sequence() == "CA") {
                        found_node_2 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((trans[e.from()] == n1->id() && !e.from_start() && trans[e.to()] == n5->id() && !e.to_end())
                         ||(trans[e.from()] == n5->id() && e.from_start() && trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((trans[e.from()] == n5->id() && !e.from_start() && trans[e.to()] == n6->id() && !e.to_end())
                              ||(trans[e.from()] == n6->id() && e.from_start() && trans[e.to()] == n5->id() && e.to_end()))) {
                        found_edge_1 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
            }
        }
        
        TEST_CASE( "Connecting graph extraction works on a particular case without leaving dangling edges",
                  "[algorithms]" ) {
                  
            string graph_json = R"(
{
  "node": [
    {
      "sequence": "GTTATAGCCTC",
      "id": 1393985
    },
    {
      "sequence": "TTTATAGTTGCTTATGTACTTGACATTGATTT",
      "id": 1393984
    },
    {
      "sequence": "A",
      "id": 1393982
    },
    {
      "sequence": "G",
      "id": 1393983
    },
    {
      "sequence": "AGTGTGTGTGTGTATATCTCATCTATCTATCT",
      "id": 1393981
    },
    {
      "sequence": "A",
      "id": 1393980
    },
    {
      "sequence": "G",
      "id": 1393979
    },
    {
      "sequence": "TCTATCTAT",
      "id": 1393978
    },
    {
      "sequence": "AT",
      "id": 1393977
    },
    {
      "sequence": "C",
      "id": 1393976
    },
    {
      "sequence": "TATC",
      "id": 1393975
    },
    {
      "sequence": "TATC",
      "id": 1393974
    },
    {
      "sequence": "TATC",
      "id": 1393973
    },
    {
      "sequence": "TATC",
      "id": 1393966
    },
    {
      "sequence": "TA",
      "id": 1393972
    },
    {
      "sequence": "CCTATCTA",
      "id": 1393971
    },
    {
      "sequence": "T",
      "id": 1393970
    },
    {
      "sequence": "TC",
      "id": 1393962
    },
    {
      "sequence": "C",
      "id": 1393969
    },
    {
      "sequence": "TATCTATC",
      "id": 1393967
    },
    {
      "sequence": "CTATCTATCTATCTAT",
      "id": 1393968
    },
    {
      "sequence": "TATC",
      "id": 1393964
    },
    {
      "sequence": "TATCTATC",
      "id": 1393965
    },
    {
      "sequence": "A",
      "id": 1393963
    },
    {
      "sequence": "TTGATCTACCTATGA",
      "id": 1393961
    },
    {
      "sequence": "T",
      "id": 1393960
    },
    {
      "sequence": "C",
      "id": 1393959
    },
    {
      "sequence": "ATTCTCATC",
      "id": 1393958
    },
    {
      "sequence": "TTCACTCTTAAATAGAGAAATTGAAGCTGTTG",
      "id": 1393957
    },
    {
      "sequence": "T",
      "id": 1393955
    },
    {
      "sequence": "G",
      "id": 1393956
    },
    {
      "sequence": "GGAGTTTGA",
      "id": 1393954
    },
    {
      "sequence": "C",
      "id": 1393953
    },
    {
      "sequence": "T",
      "id": 1393952
    },
    {
      "sequence": "TCAAAATGGTTGATCTCCAATCATAT",
      "id": 1393951
    },
    {
      "sequence": "CACAATTCTTCTCATAATATTGACATATTGAC",
      "id": 1393950
    },
    {
      "sequence": "T",
      "id": 1393948
    },
    {
      "sequence": "C",
      "id": 1393949
    },
    {
      "sequence": "TCCTCATGT",
      "id": 1393947
    },
    {
      "sequence": "AGAAGCTACACATTTCAAAAAATCTGAGTAAA",
      "id": 1393946
    }
  ],
  "edge": [
    {
      "to_end": true,
      "to": 1393947,
      "from": 1393946,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393985,
      "from": 1393984,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393950,
      "from": 1393949,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393983,
      "from": 1393981,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393950,
      "from": 1393948,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393982,
      "from": 1393981,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393976,
      "from": 1393968,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393981,
      "from": 1393979,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393980,
      "from": 1393978,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393974,
      "from": 1393967,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393967,
      "from": 1393965,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393969,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393956,
      "from": 1393954,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393961,
      "from": 1393959,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393978,
      "from": 1393976,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393949,
      "from": 1393947,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393978,
      "from": 1393962,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393968,
      "from": 1393965,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393964,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393977,
      "from": 1393969,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393957,
      "from": 1393956,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393965,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393970,
      "from": 1393969,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393963,
      "from": 1393961,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393955,
      "from": 1393954,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393958,
      "from": 1393957,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393974,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393975,
      "from": 1393974,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393954,
      "from": 1393952,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393976,
      "from": 1393975,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393961,
      "from": 1393960,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393951,
      "from": 1393950,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393953,
      "from": 1393951,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393960,
      "from": 1393958,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393984,
      "from": 1393983,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393952,
      "from": 1393951,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393977,
      "from": 1393976,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393972,
      "from": 1393971,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393966,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393967,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393973,
      "from": 1393972,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393979,
      "from": 1393978,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393969,
      "from": 1393965,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393957,
      "from": 1393955,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393959,
      "from": 1393958,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393962,
      "from": 1393961,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393968,
      "from": 1393964,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393978,
      "from": 1393977,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393973,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393975,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393970,
      "from": 1393962,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393968,
      "from": 1393963,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393954,
      "from": 1393953,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393977,
      "from": 1393962,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393984,
      "from": 1393982,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393978,
      "from": 1393969,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393969,
      "from": 1393964,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393948,
      "from": 1393947,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393971,
      "from": 1393970,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393972,
      "from": 1393970,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393974,
      "from": 1393973,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393973,
      "from": 1393966,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393981,
      "from": 1393980,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393966,
      "from": 1393964,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393966,
      "from": 1393965,
      "from_start": true
    },
    {
      "to_end": true,
      "to": 1393967,
      "from": 1393964,
      "from_start": true
    }
  ]
}
        
            )";
            
            Graph source;
            json2pb(source, graph_json.c_str(), graph_json.size());
            
            VG vg;
            vg.extend(source);
            
            VG extractor;
            
            auto max_dist = 38;
            pos_t src_pos = make_pos_t(1393981, false, 31);
            pos_t dest_pos = make_pos_t(1393958, false, 1);
            
            unordered_map<id_t, id_t> connect_trans = algorithms::extract_connecting_graph(&vg,              // DAG with split strands
                                                                                           &extractor,       // graph to extract into
                                                                                           max_dist,         // longest distance necessary
                                                                                           src_pos,          // end of earlier match
                                                                                           dest_pos,         // beginning of later match
                                                                                           false,            // do not bother finding all cycles (it's a DAG)
                                                                                           true,             // only include nodes on connecting paths
                                                                                           true);            // enforce max distance strictly
            
            Graph& dest = extractor.graph;
            
            
            // Make sure there are no dangling edges
            unordered_set<id_t> dest_nodes;
            for (const auto& node : dest.node()) {
                // Find all the nodes
                dest_nodes.insert(node.id());
            }
            for (const auto& edge : dest.edge()) {
                // Complain if an edge goes somewhere else
                REQUIRE(dest_nodes.count(edge.from()));
                REQUIRE(dest_nodes.count(edge.to()));
            }
        
        
        }
        
        TEST_CASE( "Containing graph extraction algorithm produces expected results", "[algorithms]" ) {
            
            VG vg;
            
            Node* n0 = vg.create_node("TGAG");
            Node* n1 = vg.create_node("CAGATCCACCACA");
            Node* n2 = vg.create_node("GAT");
            Node* n3 = vg.create_node("A");
            Node* n4 = vg.create_node("TGAG");
            Node* n5 = vg.create_node("AAT");
            Node* n6 = vg.create_node("TAAAC");
            Node* n7 = vg.create_node("CCGTA");
            
            vg.create_edge(n7, n0);
            vg.create_edge(n0, n1, false, true);
            vg.create_edge(n2, n0, true, true);
            vg.create_edge(n1, n2, true, false);
            vg.create_edge(n2, n3);
            vg.create_edge(n3, n4);
            vg.create_edge(n4, n5, false, true);
            vg.create_edge(n6, n4, true, true);
            
            SECTION( "Containing graph extraction works with a single maximum distance" ) {
                
                
                size_t max_len = 3;
                vector<pos_t> positions{make_pos_t(n0->id(), false, 2), make_pos_t(n5->id(), true, 1)};
                
                VG extractor;
                algorithms::extract_containing_graph(&vg, &extractor, positions, max_len);
                
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 6);
                REQUIRE(g.edge_size() == 4);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                bool found_node_4 = false;
                bool found_node_5 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    auto id = g.node(i).id();
                    if (id == n0->id()) {
                        found_node_0 = true;
                    }
                    else if (id == n1->id()) {
                        found_node_1 = true;
                    }
                    else if (id == n2->id()) {
                        found_node_2 = true;
                    }
                    else if (id == n4->id()) {
                        found_node_3 = true;
                    }
                    else if (id == n5->id()) {
                        found_node_4 = true;
                    }
                    else if (id == n7->id()) {
                        found_node_5 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                REQUIRE(found_node_4);
                REQUIRE(found_node_5);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    Edge e = g.edge(i);
                    if (((e.from() == n0->id() && e.from_start() && e.to() == n7->id() && e.to_end())
                         ||(e.from() == n7->id() && !e.from_start() && e.to() == n0->id() && !e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((e.from() == n0->id() && !e.from_start() && e.to() == n1->id() && e.to_end())
                              ||(e.from() == n1->id() && !e.from_start() && e.to() == n0->id() && e.to_end()))) {
                        found_edge_1 = true;
                    }
                    else if (((e.from() == n0->id() && !e.from_start() && e.to() == n2->id() && !e.to_end())
                              ||(e.from() == n2->id() && e.from_start() && e.to() == n0->id() && e.to_end()))) {
                        found_edge_2 = true;
                    }
                    else if (((e.from() == n4->id() && !e.from_start() && e.to() == n5->id() && e.to_end())
                              ||(e.from() == n5->id() && !e.from_start() && e.to() == n4->id() && e.to_end()))) {
                        found_edge_3 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
            }
            
            SECTION( "Containing graph extraction works with position specific maximum distances" ) {
                
                
                vector<size_t> max_lens{2, 3};
                vector<pos_t> positions{make_pos_t(n0->id(), false, 2), make_pos_t(n5->id(), true, 1)};
                
                VG extractor;
                algorithms::extract_containing_graph(&vg, &extractor, positions, max_lens);
                
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 3);
                REQUIRE(g.edge_size() == 1);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    auto id = g.node(i).id();
                    if (id == n0->id()) {
                        found_node_0 = true;
                    }
                    else if (id == n4->id()) {
                        found_node_1 = true;
                    }
                    else if (id == n5->id()) {
                        found_node_2 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                
                bool found_edge_0 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    Edge e = g.edge(i);
                    if (((e.from() == n4->id() && !e.from_start() && e.to() == n5->id() && e.to_end())
                              ||(e.from() == n5->id() && !e.from_start() && e.to() == n4->id() && e.to_end()))) {
                        found_edge_0 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
            }
            
            
            
            SECTION( "Containing graph extraction works with position and orientation specific maximum distances" ) {
                
                
                vector<size_t> forward_max_lens{3, 3};
                vector<size_t> backward_max_lens{2, 3};
                vector<pos_t> positions{make_pos_t(n0->id(), true, 2), make_pos_t(n5->id(), false, 1)};
                
                VG extractor;
                algorithms::extract_containing_graph(&vg, &extractor, positions, forward_max_lens, backward_max_lens);
                
                Graph& g = extractor.graph;
                                
                REQUIRE(g.node_size() == 4);
                REQUIRE(g.edge_size() == 2);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                
                for (int i = 0; i < g.node_size(); i++) {
                    auto id = g.node(i).id();
                    if (id == n0->id()) {
                        found_node_0 = true;
                    }
                    else if (id == n4->id()) {
                        found_node_1 = true;
                    }
                    else if (id == n5->id()) {
                        found_node_2 = true;
                    }
                    else if (id == n7->id()) {
                        found_node_3 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                
                for (int i = 0; i < g.edge_size(); i++) {
                    Edge e = g.edge(i);
                    if (((e.from() == n0->id() && e.from_start() && e.to() == n7->id() && e.to_end())
                         ||(e.from() == n7->id() && !e.from_start() && e.to() == n0->id() && !e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((e.from() == n4->id() && !e.from_start() && e.to() == n5->id() && e.to_end())
                         ||(e.from() == n5->id() && !e.from_start() && e.to() == n4->id() && e.to_end()))) {
                        found_edge_1 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
            }
        }
        
        TEST_CASE( "Extending graph extraction algorithm produces expected results", "[algorithms]" ) {
            
            VG vg;
            
            Node* n0 = vg.create_node("TGAG");
            Node* n1 = vg.create_node("CAGATCCACCACA");
            Node* n2 = vg.create_node("GAT");
            Node* n3 = vg.create_node("A");
            Node* n4 = vg.create_node("TGAG");
            Node* n5 = vg.create_node("AAT");
            Node* n6 = vg.create_node("TAAAC");
            Node* n7 = vg.create_node("CCGTA");
            
            vg.create_edge(n7, n0);
            vg.create_edge(n0, n1, false, true);
            vg.create_edge(n2, n0, true, true);
            vg.create_edge(n1, n2, true, false);
            vg.create_edge(n2, n3);
            vg.create_edge(n3, n4);
            vg.create_edge(n4, n5, false, true);
            vg.create_edge(n6, n4, true, true);
            vg.create_edge(n6, n6, false, true);
            vg.create_edge(n7, n7);
            
            SECTION( "Extending graph extraction algorithm can stay within a node and returns all of it after the cut point" ) {
                
                pos_t pos = make_pos_t(n1->id(), false, 5);
                int64_t max_dist = 1;
                bool search_backward = false;
                bool preserve_cycles = false;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 1);
                REQUIRE(g.edge_size() == 0);
                
                bool found_node_0 = false;
                
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n1->id() && n.sequence() == "CCACCACA") {
                        found_node_0 = true;
                    }
                }
                
                REQUIRE(found_node_0);
            }
            
            SECTION( "Extending graph extraction algorithm only finds nodes within the maximum distance" ) {
                
                pos_t pos = make_pos_t(n3->id(), false, 0);
                int64_t max_dist = 5;
                bool search_backward = false;
                bool preserve_cycles = false;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 2);
                REQUIRE(g.edge_size() == 1);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n3->id() && n.sequence() == "A") {
                        found_node_0 = true;
                    }
                    else if (id_trans[n.id()] == n4->id() && n.sequence() == "TGAG") {
                        found_node_1 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                
                bool found_edge_0 = false;
                
                for (size_t i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((id_trans[e.from()] == n3->id() && !e.from_start()) || (id_trans[e.to()] == n3->id() && e.to_end())) &&
                        ((id_trans[e.from()] == n4->id() && e.from_start()) || (id_trans[e.to()] == n4->id() && !e.to_end()))) {
                        found_edge_0 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
            }
            
            SECTION( "Extending graph extraction algorithm produces same output when searching in opposite direction on opposite strand" ) {
                
                pos_t pos = make_pos_t(n3->id(), true, 1);
                int64_t max_dist = 5;
                bool search_backward = true;
                bool preserve_cycles = false;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 2);
                REQUIRE(g.edge_size() == 1);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n3->id() && n.sequence() == "A") {
                        found_node_0 = true;
                    }
                    else if (id_trans[n.id()] == n4->id() && n.sequence() == "TGAG") {
                        found_node_1 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                
                bool found_edge_0 = false;
                
                for (size_t i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((id_trans[e.from()] == n3->id() && !e.from_start()) || (id_trans[e.to()] == n3->id() && e.to_end())) &&
                        ((id_trans[e.from()] == n4->id() && e.from_start()) || (id_trans[e.to()] == n4->id() && !e.to_end()))) {
                        found_edge_0 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
            }
            
            SECTION( "Extending graph extraction algorithm produces same output when searching in opposite direction on opposite strand" ) {
                
                pos_t pos = make_pos_t(n3->id(), true, 1);
                int64_t max_dist = 5;
                bool search_backward = true;
                bool preserve_cycles = false;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 2);
                REQUIRE(g.edge_size() == 1);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n3->id() && n.sequence() == "A") {
                        found_node_0 = true;
                    }
                    else if (id_trans[n.id()] == n4->id() && n.sequence() == "TGAG") {
                        found_node_1 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                
                bool found_edge_0 = false;
                
                for (size_t i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((id_trans[e.from()] == n3->id() && !e.from_start()) || (id_trans[e.to()] == n3->id() && e.to_end())) &&
                        ((id_trans[e.from()] == n4->id() && e.from_start()) || (id_trans[e.to()] == n4->id() && !e.to_end()))) {
                        found_edge_0 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
            }
            
            SECTION( "Extending graph extraction algorithm produces expected output when searching backwards" ) {
                
                pos_t pos = make_pos_t(n4->id(), false, 1);
                int64_t max_dist = 3;
                bool search_backward = true;
                bool preserve_cycles = false;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 3);
                REQUIRE(g.edge_size() == 2);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n3->id() && n.sequence() == "A") {
                        found_node_0 = true;
                    }
                    else if (id_trans[n.id()] == n2->id() && n.sequence() == "GAT") {
                        found_node_1 = true;
                    }
                    else if (id_trans[n.id()] == n4->id() && n.sequence() == "T") {
                        found_node_2 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                
                for (size_t i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((id_trans[e.from()] == n3->id() && !e.from_start()) || (id_trans[e.to()] == n3->id() && e.to_end())) &&
                        ((id_trans[e.from()] == n4->id() && e.from_start()) || (id_trans[e.to()] == n4->id() && !e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((id_trans[e.from()] == n2->id() && !e.from_start()) || (id_trans[e.to()] == n2->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n3->id() && e.from_start()) || (id_trans[e.to()] == n3->id() && !e.to_end()))) {
                        found_edge_1 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
            }
            
            SECTION( "Extending graph extraction algorithm produces same output when searching on reverse strand as searching backwards" ) {
                
                pos_t pos = make_pos_t(n4->id(), true, 3);
                int64_t max_dist = 3;
                bool search_backward = false;
                bool preserve_cycles = false;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 3);
                REQUIRE(g.edge_size() == 2);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n3->id() && n.sequence() == "A") {
                        found_node_0 = true;
                    }
                    else if (id_trans[n.id()] == n2->id() && n.sequence() == "GAT") {
                        found_node_1 = true;
                    }
                    else if (id_trans[n.id()] == n4->id() && n.sequence() == "T") {
                        found_node_2 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                
                for (size_t i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((id_trans[e.from()] == n3->id() && !e.from_start()) || (id_trans[e.to()] == n3->id() && e.to_end())) &&
                        ((id_trans[e.from()] == n4->id() && e.from_start()) || (id_trans[e.to()] == n4->id() && !e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((id_trans[e.from()] == n2->id() && !e.from_start()) || (id_trans[e.to()] == n2->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n3->id() && e.from_start()) || (id_trans[e.to()] == n3->id() && !e.to_end()))) {
                        found_edge_1 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
            }
            
            SECTION( "Extending graph extraction algorithm can traverse every type of edge" ) {
                
                pos_t pos = make_pos_t(n5->id(), false, 0);
                int64_t max_dist = 100;
                bool search_backward = false;
                bool preserve_cycles = false;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 7);
                REQUIRE(g.edge_size() == 8);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                bool found_node_4 = false;
                bool found_node_5 = false;
                bool found_node_6 = false;

                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n0->id() && n.sequence() == "TGAG") {
                        found_node_0 = true;
                    }
                    else if (id_trans[n.id()] == n1->id() && n.sequence() == "CAGATCCACCACA") {
                        found_node_1 = true;
                    }
                    else if (id_trans[n.id()] == n2->id() && n.sequence() == "GAT") {
                        found_node_2 = true;
                    }
                    else if (id_trans[n.id()] == n3->id() && n.sequence() == "A") {
                        found_node_3 = true;
                    }
                    else if (id_trans[n.id()] == n4->id() && n.sequence() == "TGAG") {
                        found_node_4 = true;
                    }
                    else if (id_trans[n.id()] == n5->id() && n.sequence() == "AAT") {
                        found_node_5 = true;
                    }
                    else if (id_trans[n.id()] == n7->id() && n.sequence() == "CCGTA") {
                        found_node_6 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                REQUIRE(found_node_4);
                REQUIRE(found_node_5);
                REQUIRE(found_node_6);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                bool found_edge_4 = false;
                bool found_edge_5 = false;
                bool found_edge_6 = false;
                bool found_edge_7 = false;
                
                for (size_t i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((id_trans[e.from()] == n7->id() && !e.from_start()) || (id_trans[e.to()] == n7->id() && e.to_end())) &&
                        ((id_trans[e.from()] == n7->id() && e.from_start()) || (id_trans[e.to()] == n7->id() && !e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((id_trans[e.from()] == n7->id() && !e.from_start()) || (id_trans[e.to()] == n7->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n0->id() && e.from_start()) || (id_trans[e.to()] == n0->id() && !e.to_end()))) {
                        found_edge_1 = true;
                    }
                    else if (((id_trans[e.from()] == n0->id() && !e.from_start()) || (id_trans[e.to()] == n0->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n1->id() && !e.from_start()) || (id_trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_2 = true;
                    }
                    else if (((id_trans[e.from()] == n0->id() && !e.from_start()) || (id_trans[e.to()] == n0->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n2->id() && e.from_start()) || (id_trans[e.to()] == n2->id() && !e.to_end()))) {
                        found_edge_3 = true;
                    }
                    else if (((id_trans[e.from()] == n1->id() && e.from_start()) || (id_trans[e.to()] == n1->id() && !e.to_end())) &&
                             ((id_trans[e.from()] == n2->id() && e.from_start()) || (id_trans[e.to()] == n2->id() && !e.to_end()))) {
                        found_edge_4 = true;
                    }
                    else if (((id_trans[e.from()] == n2->id() && !e.from_start()) || (id_trans[e.to()] == n2->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n3->id() && e.from_start()) || (id_trans[e.to()] == n3->id() && !e.to_end()))) {
                        found_edge_5 = true;
                    }
                    else if (((id_trans[e.from()] == n3->id() && !e.from_start()) || (id_trans[e.to()] == n3->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n4->id() && e.from_start()) || (id_trans[e.to()] == n4->id() && !e.to_end()))) {
                        found_edge_6 = true;
                    }
                    else if (((id_trans[e.from()] == n4->id() && !e.from_start()) || (id_trans[e.to()] == n4->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n5->id() && !e.from_start()) || (id_trans[e.to()] == n5->id() && e.to_end()))) {
                        found_edge_7 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                REQUIRE(found_edge_4);
                REQUIRE(found_edge_5);
                REQUIRE(found_edge_6);
                REQUIRE(found_edge_7);
            }
            
            SECTION( "Extending graph extraction algorithm correctly duplicates a node in a reversing cycle" ) {
                
                pos_t pos = make_pos_t(n6->id(), false, 2);
                int64_t max_dist = 10;
                bool search_backward = false;
                bool preserve_cycles = true;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                
                REQUIRE(g.node_size() == 3);
                REQUIRE(g.edge_size() == 4);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                
                id_t original_src;
                id_t duplicate_src;
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n6->id() && n.sequence() == "TAAAC") {
                        found_node_0 = true;
                        duplicate_src = n.id();
                    }
                    else if (id_trans[n.id()] == n6->id() && n.sequence() == "AAC") {
                        found_node_1 = true;
                        original_src = n.id();
                    }
                    else if (id_trans[n.id()] == n4->id() && n.sequence() == "TGAG") {
                        found_node_2 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                
                for (size_t i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((e.from() == duplicate_src && !e.from_start()) || (e.to() == duplicate_src && e.to_end())) &&
                        ((e.from() == original_src && !e.from_start()) || (e.to() == original_src && e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((e.from() == original_src && !e.from_start()) || (e.to() == original_src && e.to_end())) &&
                             ((e.from() == original_src && !e.from_start()) || (e.to() == original_src && e.to_end()))) {
                        found_edge_1 = true;
                    }
                    else if (((e.from() == duplicate_src && !e.from_start()) || (e.to() == duplicate_src && e.to_end())) &&
                             ((e.from() == duplicate_src && !e.from_start()) || (e.to() == duplicate_src && e.to_end()))) {
                        found_edge_2 = true;
                    }
                    else if (((e.from() == duplicate_src && e.from_start()) || (e.to() == duplicate_src && !e.to_end())) &&
                             ((id_trans[e.from()] == n4->id() && !e.from_start()) || (id_trans[e.to()] == n4->id() && e.to_end()))) {
                        found_edge_3 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
            }
            
            
            
            SECTION( "Extending graph extraction algorithm correctly duplicates a node in a non-reversing cycle" ) {
                
                pos_t pos = make_pos_t(n7->id(), false, 3);
                int64_t max_dist = 7;
                bool search_backward = false;
                bool preserve_cycles = true;
                
                VG extractor;
                auto id_trans = algorithms::extract_extending_graph(&vg, &extractor, max_dist, pos, search_backward, preserve_cycles);
                Graph& g = extractor.graph;
                                
                REQUIRE(g.node_size() == 5);
                REQUIRE(g.edge_size() == 6);
                
                bool found_node_0 = false;
                bool found_node_1 = false;
                bool found_node_2 = false;
                bool found_node_3 = false;
                bool found_node_4 = false;
                
                id_t original_src;
                id_t duplicate_src;
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    if (id_trans[n.id()] == n7->id() && n.sequence() == "CCGTA") {
                        found_node_0 = true;
                        duplicate_src = n.id();
                    }
                    else if (id_trans[n.id()] == n7->id() && n.sequence() == "TA") {
                        found_node_1 = true;
                        original_src = n.id();
                    }
                    else if (id_trans[n.id()] == n0->id() && n.sequence() == "TGAG") {
                        found_node_2 = true;
                    }
                    else if (id_trans[n.id()] == n1->id() && n.sequence() == "CAGATCCACCACA") {
                        found_node_3 = true;
                    }
                    else if (id_trans[n.id()] == n2->id() && n.sequence() == "GAT") {
                        found_node_4 = true;
                    }
                }
                
                REQUIRE(found_node_0);
                REQUIRE(found_node_1);
                REQUIRE(found_node_2);
                REQUIRE(found_node_3);
                REQUIRE(found_node_4);
                
                bool found_edge_0 = false;
                bool found_edge_1 = false;
                bool found_edge_2 = false;
                bool found_edge_3 = false;
                bool found_edge_4 = false;
                bool found_edge_5 = false;
                
                for (size_t i = 0; i < g.edge_size(); i++) {
                    const Edge& e = g.edge(i);
                    if (((e.from() == duplicate_src && e.from_start()) || (e.to() == duplicate_src && !e.to_end())) &&
                        ((e.from() == original_src && !e.from_start()) || (e.to() == original_src && e.to_end()))) {
                        found_edge_0 = true;
                    }
                    else if (((e.from() == duplicate_src && e.from_start()) || (e.to() == duplicate_src && !e.to_end())) &&
                             ((e.from() == duplicate_src && !e.from_start()) || (e.to() == duplicate_src && e.to_end()))) {
                        found_edge_1 = true;
                    }
                    else if (((e.from() == original_src && !e.from_start()) || (e.to() == original_src && e.to_end())) &&
                             ((id_trans[e.from()] == n0->id() && e.from_start()) || (id_trans[e.to()] == n0->id() && !e.to_end()))) {
                        found_edge_2 = true;
                    }
                    else if (((e.from() == duplicate_src && !e.from_start()) || (e.to() == duplicate_src && e.to_end())) &&
                             ((id_trans[e.from()] == n0->id() && e.from_start()) || (id_trans[e.to()] == n0->id() && !e.to_end()))) {
                        found_edge_3 = true;
                    }
                    else if (((id_trans[e.from()] == n0->id() && !e.from_start()) || (id_trans[e.to()] == n0->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n1->id() && !e.from_start()) || (id_trans[e.to()] == n1->id() && e.to_end()))) {
                        found_edge_4 = true;
                    }
                    else if (((id_trans[e.from()] == n0->id() && !e.from_start()) || (id_trans[e.to()] == n0->id() && e.to_end())) &&
                             ((id_trans[e.from()] == n2->id() && e.from_start()) || (id_trans[e.to()] == n2->id() && !e.to_end()))) {
                        found_edge_5 = true;
                    }
                }
                
                REQUIRE(found_edge_0);
                REQUIRE(found_edge_1);
                REQUIRE(found_edge_2);
                REQUIRE(found_edge_3);
                REQUIRE(found_edge_4);
                REQUIRE(found_edge_5);
            }
        }
        
        TEST_CASE( "Topological sort works on a small graph",
                  "[algorithms][topologicalsort]" ) {
            
            VG vg;
            
            Node* n0 = vg.create_node("CGA");
            Node* n1 = vg.create_node("TTGG");
            Node* n2 = vg.create_node("CCGT");
            Node* n3 = vg.create_node("C");
            Node* n4 = vg.create_node("GT");
            Node* n5 = vg.create_node("GATAA");
            Node* n6 = vg.create_node("CGG");
            Node* n7 = vg.create_node("ACA");
            Node* n8 = vg.create_node("GCCG");
            Node* n9 = vg.create_node("ATATAAC");
            
            vg.create_edge(n1, n0, true, true); // a doubly reversing edge to keep it interesting
            vg.create_edge(n1, n2);
            vg.create_edge(n2, n3);
            vg.create_edge(n2, n4);
            vg.create_edge(n3, n5);
            vg.create_edge(n4, n5);
            vg.create_edge(n5, n6);
            vg.create_edge(n5, n8);
            vg.create_edge(n6, n7);
            vg.create_edge(n6, n8);
            vg.create_edge(n7, n9);
            vg.create_edge(n8, n9);
            
            SECTION( "algorithms::topological_order produces a consistent total ordering and orientation" ) {
                auto handle_sort = algorithms::topological_order(&vg);

                SECTION( "Ordering and orientation is consistent" ) {
                
                    unordered_map<id_t, handle_t> oriented;
                    
                    for (auto& handle : handle_sort) {
                        // For each oriented node
                        
                        // What was before it
                        vector<handle_t> prev_handles;
                        vg.follow_edges(handle, true, [&](const handle_t& prev) {
                            prev_handles.push_back(prev);
                        });
                        
                        for (auto& prev : prev_handles) {
                            // We should have visited this already
                            REQUIRE(oriented.count(vg.get_id(prev)) != 0);
                            // And it should have been in the correct orientation
                            REQUIRE(oriented.at(vg.get_id(prev)) == prev);
                        }
                        
                        // Remember we were here in this orientation
                        oriented.insert(make_pair(vg.get_id(handle), handle));
                    }
                    
                }
                
                SECTION( "All nodes are placed" ) {
                    REQUIRE(handle_sort.size() == vg.node_size());
                
                    unordered_set<id_t> found;
                    
                    for (auto& handle : handle_sort) {
                        found.insert(vg.get_id(handle));
                    }
                    
                    REQUIRE(found.size() == vg.node_size());
                    
                }
               
            }
        }
        
        TEST_CASE( "Topological sort works on a more complex graph",
                  "[algorithms][topologicalsort]" ) {


            string graph_json = R"(
            {"node": [{"id": 1, "sequence": "GTATTTTTAGTA"}, {"id": 2, "sequence": "G"}, {"id": 3, "sequence": "GAGACGGGGTTTCACCATGTT"}, {"id": 4, "sequence": "T"}, {"id": 5, "sequence": "CTAATTTTT"}, {"id": 6, "sequence": "CA"}, {"id": 7, "sequence": "GG"}, {"id": 8, "sequence": "ACGCCC"}, {"id": 9, "sequence": "C"}, {"id": 10, "sequence": "T"}, {"id": 11, "sequence": "C"}, {"id": 12, "sequence": "GCCA"}, {"id": 13, "sequence": "A"}, {"id": 14, "sequence": "GGGATTACAGGCGCACACC"}, {"id": 15, "sequence": "CCACACC"}, {"id": 16, "sequence": "AT"}, {"id": 17, "sequence": "CC"}, {"id": 18, "sequence": "GGTCAGGCTGGTCTCGACTCC"}, {"id": 19, "sequence": "TGACCTCCTGATCTGCCCCCC"}, {"id": 20, "sequence": "A"}, {"id": 21, "sequence": "G"}, {"id": 22, "sequence": "TATTTTTAGTA"}, {"id": 23, "sequence": "A"}, {"id": 24, "sequence": "G"}, {"id": 25, "sequence": "GA"}], "edge": [{"from": 4, "to": 1}, {"from": 5, "to": 1}, {"from": 1, "to": 2}, {"from": 1, "to": 3}, {"from": 22, "to": 2}, {"from": 2, "to": 20}, {"from": 2, "to": 21}, {"from": 3, "to": 18}, {"from": 5, "to": 4}, {"from": 6, "to": 5}, {"from": 7, "to": 5}, {"from": 8, "to": 6}, {"from": 8, "to": 7}, {"from": 9, "to": 8}, {"from": 10, "to": 8}, {"from": 11, "to": 9}, {"from": 11, "to": 10}, {"from": 12, "to": 11}, {"from": 13, "to": 11}, {"from": 16, "to": 12}, {"from": 17, "to": 12}, {"from": 12, "to": 15}, {"from": 14, "to": 13}, {"from": 18, "to": 19}, {"from": 20, "to": 25}, {"from": 21, "to": 25}, {"from": 23, "to": 22}, {"from": 24, "to": 22}]}
            )";
            
            // Load the JSON
            Graph proto_graph;
            json2pb(proto_graph, graph_json.c_str(), graph_json.size());
            
            // Make it into a VG
            VG vg;
            vg.extend(proto_graph);
            
            SECTION( "algorithms::topological_order produces a consistent total ordering and orientation" ) {
                auto handle_sort = algorithms::topological_order(&vg);

                SECTION( "Ordering and orientation is consistent" ) {
                
                    unordered_map<id_t, handle_t> oriented;
                    
                    for (auto& handle : handle_sort) {
                        // For each oriented node
                        
                        // What was before it
                        vector<handle_t> prev_handles;
                        vg.follow_edges(handle, true, [&](const handle_t& prev) {
                            prev_handles.push_back(prev);
                        });
                        
                        for (auto& prev : prev_handles) {
                            // We should have visited this already
                            REQUIRE(oriented.count(vg.get_id(prev)) != 0);
                            // And it should have been in the correct orientation
                            REQUIRE(oriented.at(vg.get_id(prev)) == prev);
                        }
                        
                        // Remember we were here in this orientation
                        oriented.insert(make_pair(vg.get_id(handle), handle));
                    }
                    
                }
                
                SECTION( "All nodes are placed" ) {
                    REQUIRE(handle_sort.size() == vg.node_size());
                
                    unordered_set<id_t> found;
                    
                    for (auto& handle : handle_sort) {
                        found.insert(vg.get_id(handle));
                    }
                    
                    REQUIRE(found.size() == vg.node_size());
                    
                }
               
            }
                  
        }
        
        TEST_CASE( "Weakly connected components works",
                  "[algorithms]" ) {
            
            VG vg;
            
            Node* n0 = vg.create_node("CGA");
            Node* n1 = vg.create_node("TTGG");
            Node* n2 = vg.create_node("CCGT");
            Node* n3 = vg.create_node("C");
            Node* n4 = vg.create_node("GT");
            Node* n5 = vg.create_node("GATAA");
            Node* n6 = vg.create_node("CGG");
            Node* n7 = vg.create_node("ACA");
            Node* n8 = vg.create_node("GCCG");
            Node* n9 = vg.create_node("ATATAAC");
            
            vg.create_edge(n1, n0, true, true); // a doubly reversing edge to keep it interesting
            vg.create_edge(n1, n2);
            vg.create_edge(n2, n3);
            vg.create_edge(n2, n4);
            
            // Skip edges from 3 and 4 to 5 so we have 2 components
            
            vg.create_edge(n5, n6);
            vg.create_edge(n5, n8);
            vg.create_edge(n6, n7);
            vg.create_edge(n6, n8);
            vg.create_edge(n7, n9);
            vg.create_edge(n8, n9);
            
            SECTION( "algorithms::weakly_connected_components finds two components" ) {
                auto components = algorithms::weakly_connected_components(&vg);

                REQUIRE(components.size() == 2);
                
                SECTION( "The components are the first and second parts of the graph" ){
                    set<size_t> sizes;
                    for (auto& component : components) {
                        sizes.insert(component.size());
                    }
                    
                    // Actually these are both the same size
                    REQUIRE(sizes == set<size_t>{5});
                }
            }
            
            vg.create_edge(n3, n5);
            vg.create_edge(n4, n5);
            
            SECTION( "Components can be joined" ) {
            
                auto components = algorithms::weakly_connected_components(&vg);

                REQUIRE(components.size() == 1);
                REQUIRE(components.front().size() == 10);
            
            }
        }
        TEST_CASE("distance_to_head() using HandleGraph produces expected results", "[vg]") {
            VG vg;
            Node* n0 = vg.create_node("AA");
            Node* n1 = vg.create_node("ACTGA");
            Node* n2 = vg.create_node("AG");
            Node* n3 = vg.create_node("ATC");
            unordered_set<handle_t> trav;
            handle_t n;

            SECTION("distance_to_head() using HandleGraph works when node is at head") {
                // Set handle to head node
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_head(n, limit, 0, trav, &vg);
                
                REQUIRE(check == 0);
            }
            SECTION("distance_to_head() using HandleGraph works when node is not at head") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n3->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_head(n, limit, 0, trav, &vg);
                
                REQUIRE(check == 9);
            }
            SECTION("distance_to_head() using HandleGraph works when limit is less than the distance") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n3->id(),false);
                int32_t limit = 4;

                int32_t check = algorithms::distance_to_head(n, limit, 0, trav, &vg);
                
                REQUIRE(check == -1);
            }
            SECTION("distance_to_head() using HandleGraph works when limit is equal to the distance") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n3->id(),false);
                int32_t limit = 6;

                int32_t check = algorithms::distance_to_head(n, limit, 0, trav, &vg);
                
                REQUIRE(check == -1);
            }
              SECTION("distance_to_head() using HandleGraph works when there is no head") {
                vg.create_edge(n0, n1, true, false);
                // Set handle to the node you are currently on
                n = vg.get_handle(n1->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_head(n, limit, 0, trav, &vg);
                
                REQUIRE(check == -1);
            }
              SECTION("distance_to_head() using HandleGraph works when there are 2 previous nodes") {
                Node* n4 = vg.create_node("AG");
                Node* n5 = vg.create_node("AT");
                Node* n6 = vg.create_node("GATTACA");
                vg.create_edge(n0, n4);
                vg.create_edge(n4, n5);
                vg.create_edge(n0, n6);
                vg.create_edge(n6, n1);
                vg.create_edge(n5, n1);
                // Set handle to the node you are currently on
                n = vg.get_handle(n1->id(),false);
                int32_t limit = 100;
                // if there are 2 previous nodes, gets the distance of previous node that was made first
                int32_t check = algorithms::distance_to_head(n, limit, 0, trav, &vg);
                
                REQUIRE(check == 9);
            }
        }

        TEST_CASE("Simplified distance_to_head() using HandleGraph produces expected results", "[vg]") {
            VG vg;
            Node* n0 = vg.create_node("AA");
            Node* n1 = vg.create_node("ACTGA");
            Node* n2 = vg.create_node("AG");
            Node* n3 = vg.create_node("ATC");
            handle_t n;

            SECTION("Simplified distance_to_head() using HandleGraph works when node is at head") {
                // Set handle to head node
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_head(n, limit, &vg);
                
                REQUIRE(check == 0);
            }
            SECTION("Simplified distance_to_head() using HandleGraph works when node is not at head") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n3->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_head(n, limit, &vg);
                
                REQUIRE(check == 9);
            }
            SECTION("Simplified distance_to_head() using HandleGraph works when limit is less than the distance") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n3->id(),false);
                int32_t limit = 4;

                int32_t check = algorithms::distance_to_head(n, limit, &vg);
                
                REQUIRE(check == -1);
            }
            SECTION("Simplified distance_to_head() using HandleGraph works when limit is equal to the distance") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n3->id(),false);
                int32_t limit = 6;

                int32_t check = algorithms::distance_to_head(n, limit, &vg);
                
                REQUIRE(check == -1);
            }
              SECTION("Simplified distance_to_head() using HandleGraph works when there is no head") {
                vg.create_edge(n0, n1, true, false);
                // Set handle to the node you are currently on
                n = vg.get_handle(n1->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_head(n, limit, &vg);
                
                REQUIRE(check == -1);
            }
              SECTION("Simplified distance_to_head() using HandleGraph works when there are 2 previous nodes") {
                Node* n4 = vg.create_node("AG");
                Node* n5 = vg.create_node("AT");
                Node* n6 = vg.create_node("GATTACA");
                vg.create_edge(n0, n4);
                vg.create_edge(n4, n5);
                vg.create_edge(n0, n6);
                vg.create_edge(n5, n1);
                vg.create_edge(n6, n1);
        
                // Set handle to the node you are currently on
                n = vg.get_handle(n1->id(),false);
                int32_t limit = 100;
                // if there are 2 previous nodes, gets the distance of previous node that was made first
                int32_t check = algorithms::distance_to_head(n, limit, &vg);
                
                REQUIRE(check == 6);
            }
        }

        TEST_CASE("distance_to_tail() using HandleGraph produces expected results", "[vg]") {
            VG vg;
            Node* n0 = vg.create_node("AA");
            Node* n1 = vg.create_node("ACTGA");
            Node* n2 = vg.create_node("AG");
            Node* n3 = vg.create_node("ATC");
            unordered_set<handle_t> trav;
            handle_t n;

            SECTION("distance_to_tail() using HandleGraph works when node is at tail") {
                // Set handle to head node
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_tail(n, limit, 0, trav, &vg);
                
                REQUIRE(check == 0);
            }
            SECTION("distance_to_tail() using HandleGraph works when node is not at tail") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_tail(n, limit, 0, trav, &vg);
                
                REQUIRE(check == 10);
            }
            SECTION("distance_to_tail() using HandleGraph works when limit is less than the distance") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 4;

                int32_t check = algorithms::distance_to_tail(n, limit, 0, trav, &vg);
                
                REQUIRE(check == -1);
            }
            SECTION("distance_to_tail() using HandleGraph works when limit is equal to the distance") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n3->id(),false);
                int32_t limit = 6;

                int32_t check = algorithms::distance_to_tail(n, limit, 0, trav, &vg);
                
                REQUIRE(check == 0);
            }
              SECTION("distance_to_tail() using HandleGraph works when there is no tail") {
                vg.create_edge(n0, n1, false, true);
                // Set handle to the node you are currently on
                n = vg.get_handle(n1->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_tail(n, limit, 0, trav, &vg);
                
                REQUIRE(check == -1);
            }
              SECTION("distance_to_tail() using HandleGraph works when there are 3 next nodes") {
                Node* n4 = vg.create_node("AG");
                Node* n5 = vg.create_node("AT");
                Node* n6 = vg.create_node("GATTACA");
                Node* n7 = vg.create_node("ATG");
                vg.create_edge(n0, n6); 
                vg.create_edge(n6, n1); 
                vg.create_edge(n0, n2); 
                vg.create_edge(n2, n7);
                vg.create_edge(n7, n1); 
                vg.create_edge(n0, n4);
                vg.create_edge(n4, n5);
                vg.create_edge(n5, n1); 
                
                // Set handle to the node you are currently on
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 100;
                // if there are 3 next nodes, gets the distance of next node that was made first
                int32_t check = algorithms::distance_to_tail(n, limit, 0, trav, &vg);
                
                REQUIRE(check == 12);
            }
        }
        TEST_CASE("Simplified distance_to_tail() using HandleGraph produces expected results", "[vg]") {
            VG vg;
            Node* n0 = vg.create_node("AA");
            Node* n1 = vg.create_node("ACTGA");
            Node* n2 = vg.create_node("AG");
            Node* n3 = vg.create_node("ATC");
            unordered_set<handle_t> trav;
            handle_t n;

            SECTION("Simplified distance_to_tail() using HandleGraph works when node is at tail") {
                // Set handle to head node
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_tail(n, limit, &vg);
                
                REQUIRE(check == 0);
            }
            SECTION("Simplified distance_to_tail() using HandleGraph works when node is not at tail") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_tail(n, limit, &vg);
                
                REQUIRE(check == 10);
            }
            SECTION("Simplified distance_to_tail() using HandleGraph works when limit is less than the distance") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 4;

                int32_t check = algorithms::distance_to_tail(n, limit, &vg);
                
                REQUIRE(check == -1);
            }
            SECTION("Simplified distance_to_tail() using HandleGraph works when limit is equal to the distance") {
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n3);
                // Set handle to the node you are currently on
                n = vg.get_handle(n3->id(),false);
                int32_t limit = 6;

                int32_t check = algorithms::distance_to_tail(n, limit, &vg);
                
                REQUIRE(check == 0);
            }
              SECTION("Simplified distance_to_tail() using HandleGraph works when there is no tail") {
                vg.create_edge(n0, n1, false, true);
                // Set handle to the node you are currently on
                n = vg.get_handle(n1->id(),false);
                int32_t limit = 100;

                int32_t check = algorithms::distance_to_tail(n, limit, &vg);
                
                REQUIRE(check == -1);
            }
              SECTION("Simplified distance_to_tail() using HandleGraph works when there are 3 next nodes") {
                Node* n4 = vg.create_node("AG");
                Node* n5 = vg.create_node("AT");
                Node* n6 = vg.create_node("GATTACA");
                Node* n7 = vg.create_node("ATG");
                vg.create_edge(n0, n6); 
                vg.create_edge(n6, n1); 
                vg.create_edge(n0, n2); 
                vg.create_edge(n2, n7);
                vg.create_edge(n7, n1); 
                vg.create_edge(n0, n4);
                vg.create_edge(n4, n5);
                vg.create_edge(n5, n1); 
                
                // Set handle to the node you are currently on
                n = vg.get_handle(n0->id(),false);
                int32_t limit = 100;
                // if there are 3 next nodes, gets the distance of next node that was made first
                int32_t check = algorithms::distance_to_tail(n, limit, &vg);
                
                REQUIRE(check == 12);
            }
        }

        TEST_CASE("is_directed_acyclic() should return whether the graph is directed acyclic", "[algorithms]") {
            SECTION("is_directed_acyclic() works on a single node") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                
                xg::XG xg1(vg.graph);
                
                // the graph has no edges
                REQUIRE(algorithms::is_directed_acyclic(&vg));
                REQUIRE(algorithms::is_directed_acyclic(&xg1));
                
                vg.create_edge(n0, n0, false, true);
                
                xg::XG xg2(vg.graph);
                
                // the graph has a reversing cycle, but no directed cycles
                REQUIRE(algorithms::is_directed_acyclic(&vg));
                REQUIRE(algorithms::is_directed_acyclic(&xg2));
                
                vg.create_edge(n0, n0, true, false);
                
                xg::XG xg3(vg.graph);
                
                // the graph now has a directed cycle
                REQUIRE(!algorithms::is_directed_acyclic(&vg));
                REQUIRE(!algorithms::is_directed_acyclic(&xg3));
            }
            
            SECTION("is_directed_acyclic() works on DAG with only simple edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                Node* n6 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n0, n2);
                vg.create_edge(n2, n3);
                vg.create_edge(n1, n3);
                vg.create_edge(n3, n4);
                vg.create_edge(n0, n4);
                vg.create_edge(n4, n5);
                vg.create_edge(n0, n5);
                vg.create_edge(n4, n6);
                vg.create_edge(n3, n6);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::is_directed_acyclic(&vg));
                REQUIRE(algorithms::is_directed_acyclic(&xg1));
            }
            
            SECTION("is_directed_acyclic() works on DAG with some doubly reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                Node* n6 = vg.create_node("A");
                
                vg.create_edge(n1, n0, true, true);
                vg.create_edge(n0, n2);
                vg.create_edge(n2, n3);
                vg.create_edge(n3, n1, true, true);
                vg.create_edge(n4, n3, true, true);
                vg.create_edge(n0, n4);
                vg.create_edge(n4, n5);
                vg.create_edge(n0, n5);
                vg.create_edge(n6, n4, true, true);
                vg.create_edge(n3, n6);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::is_directed_acyclic(&vg));
                REQUIRE(algorithms::is_directed_acyclic(&xg1));
            }
            
            SECTION("is_directed_acyclic() works on DAG with doubly reversing and singly reversing eges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                Node* n6 = vg.create_node("A");
                
                vg.create_edge(n1, n0, true, true);
                vg.create_edge(n0, n2, false, true);
                vg.create_edge(n2, n3, true, false);
                vg.create_edge(n3, n1, true, true);
                vg.create_edge(n4, n3, true, true);
                vg.create_edge(n0, n4);
                vg.create_edge(n4, n5, false, true);
                vg.create_edge(n0, n5, false, true);
                vg.create_edge(n6, n4, true, true);
                vg.create_edge(n3, n6);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::is_directed_acyclic(&vg));
                REQUIRE(algorithms::is_directed_acyclic(&xg1));
            }
            
            SECTION("is_directed_acyclic() works on a non trivial graph a reversing cycle but no directed cycles") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n1, n3);
                vg.create_edge(n2, n3, false, true);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::is_directed_acyclic(&vg));
                REQUIRE(algorithms::is_directed_acyclic(&xg1));
            }
            
            SECTION("is_directed_acyclic() works on a simple directed cycle") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n0);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(!algorithms::is_directed_acyclic(&vg));
                REQUIRE(!algorithms::is_directed_acyclic(&xg1));
            }
            
            SECTION("is_directed_acyclic() works on a non trivial graph with a directed cycle") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n0);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n2, n5);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(!algorithms::is_directed_acyclic(&vg));
                REQUIRE(!algorithms::is_directed_acyclic(&xg1));
            }
        }
        
        
        TEST_CASE("is_single_stranded() correctly identifies graphs with reversing edges", "[algorithms]") {
            
            SECTION("is_single_stranded() works a trivial graph with no edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::is_single_stranded(&vg));
                REQUIRE(algorithms::is_single_stranded(&xg1));
            }
            
            SECTION("is_single_stranded() works a non trivial graph with no reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n0);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n2, n5);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::is_single_stranded(&vg));
                REQUIRE(algorithms::is_single_stranded(&xg1));
            }
            
            SECTION("is_single_stranded() works a non trivial graph with a directed cycle but no reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n0);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n2, n5);
                vg.create_edge(n4, n0);
                vg.create_edge(n5, n0);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::is_single_stranded(&vg));
                REQUIRE(algorithms::is_single_stranded(&xg1));
            }
            
            SECTION("is_single_stranded() works a non trivial graph with no reversing edges, but with doubly reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n0, n2, true, true);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n5, n2, true, true);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::is_single_stranded(&vg));
                REQUIRE(algorithms::is_single_stranded(&xg1));
            }
            
            SECTION("is_single_stranded() works a non trivial graph with reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n0, n2, true, false);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n5, n2, true, true);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(!algorithms::is_single_stranded(&vg));
                REQUIRE(!algorithms::is_single_stranded(&xg1));
            }
            
            SECTION("is_single_stranded() works a non trivial graph with reversing edges in the opposite orientation") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n0, n2, false, true);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n5, n2, true, true);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(!algorithms::is_single_stranded(&vg));
                REQUIRE(!algorithms::is_single_stranded(&xg1));
            }
        }
        
        TEST_CASE("single_stranded_orientation() can identify orientations of nodes that produce a single stranded graph", "[algorithms]") {
            
            auto validate_single_stranded_orientation = [](const HandleGraph* g, const vector<handle_t>& orientation) {
                unordered_map<id_t, bool> orientation_by_id;
                for (handle_t handle : orientation) {
                    orientation_by_id[g->get_id(handle)] = g->get_is_reverse(handle);
                }
                
                bool pass = true;
                function<bool(const handle_t&)> lam = [&](const handle_t& next) {
                    pass = pass && (g->get_is_reverse(next) == orientation_by_id[g->get_id(next)]);
                    return pass;
                };
                
                for (handle_t handle : orientation) {
                    g->follow_edges(handle, true, lam);
                    g->follow_edges(handle, false, lam);
                }
                
                return pass;
            };
            
            SECTION("single_stranded_orientation() works a trivial graph with no edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(validate_single_stranded_orientation(&vg, algorithms::single_stranded_orientation(&vg)));
                REQUIRE(validate_single_stranded_orientation(&xg1, algorithms::single_stranded_orientation(&xg1)));
            }
            
            SECTION("single_stranded_orientation() works a non trivial graph with no reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n0);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n2, n5);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(validate_single_stranded_orientation(&vg, algorithms::single_stranded_orientation(&vg)));
                REQUIRE(validate_single_stranded_orientation(&xg1, algorithms::single_stranded_orientation(&xg1)));
            }
            
            SECTION("single_stranded_orientation() works a non trivial graph with a directed cycle but no reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n2, n0);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n2, n5);
                vg.create_edge(n4, n0);
                vg.create_edge(n5, n0);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(validate_single_stranded_orientation(&vg, algorithms::single_stranded_orientation(&vg)));
                REQUIRE(validate_single_stranded_orientation(&xg1, algorithms::single_stranded_orientation(&xg1)));
            }
            
            SECTION("single_stranded_orientation() works a non trivial graph with no reversing edges, but with doubly reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n0, n2, true, true);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n5, n2, true, true);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(validate_single_stranded_orientation(&vg, algorithms::single_stranded_orientation(&vg)));
                REQUIRE(validate_single_stranded_orientation(&xg1, algorithms::single_stranded_orientation(&xg1)));
            }
            
            SECTION("single_stranded_orientation() works a non trivial, single-strand-able graph with reversing edges") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2, false, true);
                vg.create_edge(n2, n0, true, false);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n2, n5, true, false);
                vg.create_edge(n4, n0);
                vg.create_edge(n5, n0);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(validate_single_stranded_orientation(&vg, algorithms::single_stranded_orientation(&vg)));
                REQUIRE(validate_single_stranded_orientation(&xg1, algorithms::single_stranded_orientation(&xg1)));
            }
            
            SECTION("single_stranded_orientation() correctly identifies a graph with no single stranded orientation") {
                
                VG vg;
                
                Node* n0 = vg.create_node("A");
                Node* n1 = vg.create_node("A");
                Node* n2 = vg.create_node("A");
                Node* n3 = vg.create_node("A");
                Node* n4 = vg.create_node("A");
                Node* n5 = vg.create_node("A");
                
                vg.create_edge(n0, n1);
                vg.create_edge(n1, n2);
                vg.create_edge(n0, n2, false, true);
                vg.create_edge(n0, n3);
                vg.create_edge(n1, n4);
                vg.create_edge(n5, n2, true, true);
                
                xg::XG xg1(vg.graph);
                
                REQUIRE(algorithms::single_stranded_orientation(&vg).empty());
                REQUIRE(algorithms::single_stranded_orientation(&xg1).empty());
            }
        }

        TEST_CASE("lazy_topological_sort() and lazier_topological_sort() should put a DAG in topological order", "[algorithms][sort]") {
            
            auto is_in_topological_order = [](const Graph& graph) {
                
                unordered_map<id_t, bool> node_orientation;
                
                unordered_map<id_t, size_t> id_to_idx;
                for (size_t i = 0; i < graph.node_size(); i++) {
                    id_to_idx[graph.node(i).id()] = i;
                }
                
                bool return_val = true;
                for (size_t i = 0; i < graph.edge_size() && return_val; i++) {
                    const Edge& e = graph.edge(i);
                    
                    if (id_to_idx[e.from()] < id_to_idx[e.to()]) {
                        if (node_orientation.count(e.from())) {
                            return_val = return_val && (e.from_start() == node_orientation[e.from()]);
                        }
                        else {
                            node_orientation[e.from()] = e.from_start();
                        }
                        
                        if (node_orientation.count(e.to())) {
                            return_val = return_val && (e.to_end() == node_orientation[e.to()]);
                        }
                        else {
                            node_orientation[e.to()] = e.to_end();
                        }
                    }
                    else {
                        if (node_orientation.count(e.from())) {
                            return_val = return_val && (e.from_start() != node_orientation[e.from()]);
                        }
                        else {
                            node_orientation[e.from()] = !e.from_start();
                        }
                        
                        if (node_orientation.count(e.to())) {
                            return_val = return_val && (e.to_end() != node_orientation[e.to()]);
                        }
                        else {
                            node_orientation[e.to()] = !e.to_end();
                        }
                    }
                    
                    if (e.from() == e.to()) {
                        return_val = false;
                    }
                    else if (id_to_idx[e.from()] < id_to_idx[e.to()]) {
                        return_val = return_val && (e.from_start() == node_orientation[e.from()] && e.to_end() == node_orientation[e.to()]);
                    }
                    else {
                        return_val = return_val && (e.from_start() != node_orientation[e.from()] && e.to_end() != node_orientation[e.to()]);
                    }
                }
                return return_val;
            };
            
            SECTION("laz[y/ier]_topological_sort() works on a simple graph that's already in topological order") {
                
                VG vg1;
                
                Node* n0 = vg1.create_node("A");
                Node* n1 = vg1.create_node("A");
                
                vg1.create_edge(n0, n1);
                
                VG vg2 = vg1;
                // make the second graph have some locally stored nodes in the reverse orientation
                vg2.apply_orientation(vg2.get_handle(n1->id(), true));
                
                algorithms::lazier_topological_sort(&vg1);
                algorithms::lazy_topological_sort(&vg2);
                
                REQUIRE(is_in_topological_order(vg1.graph));
                REQUIRE(is_in_topological_order(vg2.graph));
            }
            
            SECTION("laz[y/ier]_topological_sort() works on a simple graph that's not already in topological order") {
                
                VG vg1;
                
                Node* n0 = vg1.create_node("A");
                Node* n1 = vg1.create_node("A");
                
                vg1.create_edge(n1, n0);
                
                VG vg2 = vg1;
                // make the second graph have some locally stored nodes in the reverse orientation
                vg2.apply_orientation(vg2.get_handle(n1->id(), true));
                
                algorithms::lazier_topological_sort(&vg1);
                algorithms::lazy_topological_sort(&vg2);
                
                REQUIRE(is_in_topological_order(vg1.graph));
                REQUIRE(is_in_topological_order(vg2.graph));

            }
            
            SECTION("laz[y/ier]_topological_sort() works on a more complex graph that's not already in topological order") {
                
                VG vg1;
                
                Node* n0 = vg1.create_node("A");
                Node* n9 = vg1.create_node("A");
                Node* n2 = vg1.create_node("A");
                Node* n1 = vg1.create_node("A");
                Node* n4 = vg1.create_node("A");
                Node* n8 = vg1.create_node("A");
                Node* n5 = vg1.create_node("A");
                Node* n3 = vg1.create_node("A");
                Node* n7 = vg1.create_node("A");
                Node* n6 = vg1.create_node("A");
                
                vg1.create_edge(n0, n1);
                vg1.create_edge(n0, n3);
                vg1.create_edge(n0, n8);
                vg1.create_edge(n1, n2);
                vg1.create_edge(n1, n9);
                vg1.create_edge(n2, n5);
                vg1.create_edge(n3, n4);
                vg1.create_edge(n3, n7);
                vg1.create_edge(n5, n6);
                vg1.create_edge(n5, n8);
                vg1.create_edge(n6, n8);
                vg1.create_edge(n8, n9);
                
                VG vg2 = vg1;
                // make the second graph have some locally stored nodes in the reverse orientation
                vg2.apply_orientation(vg2.get_handle(n1->id(), true));
                vg2.apply_orientation(vg2.get_handle(n3->id(), true));
                vg2.apply_orientation(vg2.get_handle(n8->id(), true));
                vg2.apply_orientation(vg2.get_handle(n6->id(), true));
                
                algorithms::lazier_topological_sort(&vg1);
                algorithms::lazy_topological_sort(&vg2);
                
                REQUIRE(is_in_topological_order(vg1.graph));
                REQUIRE(is_in_topological_order(vg2.graph));
            }
        }
        
        TEST_CASE("apply_ordering and apply_orientations work as expected", "[algorithms]") {
            
            VG vg;
            
            Node* n0 = vg.create_node("A");
            Node* n1 = vg.create_node("A");
            Node* n2 = vg.create_node("A");
            Node* n3 = vg.create_node("A");
            Node* n4 = vg.create_node("A");
            Node* n5 = vg.create_node("A");
            
            random_device rd;
            default_random_engine prng(rd());
            
            SECTION("apply_ordering works with random permutations") {
                
                for (int i = 0; i < 20; i++) {
                    
                    vector<handle_t> order;
                    vg.for_each_handle([&](const handle_t& handle) {
                        order.push_back(handle);
                    });
                    shuffle(order.begin(), order.end(), prng);
                    
                    vector<id_t> id_order;
                    for (handle_t handle : order){
                        id_order.push_back(vg.get_id(handle));
                    }
                    
                    algorithms::apply_ordering(&vg, order);
                    
                    int idx = 0;
                    vg.for_each_handle([&](const handle_t& handle) {
                        REQUIRE(vg.get_id(handle) == id_order[idx]);
                        idx++;
                    });
                }
            }
            
            SECTION("apply_orientations works with random orientations") {
                
                for (int i = 0; i < 20; i++) {
                    
                    vector<handle_t> orientation;
                    unordered_map<id_t, string> sequence_by_id;
                    vg.for_each_handle([&](const handle_t& handle) {
                        orientation.push_back(bernoulli_distribution()(prng) ? handle : vg.flip(handle));
                        sequence_by_id[vg.get_id(handle)] = vg.get_sequence(orientation.back());
                    });
                    
                    algorithms::apply_orientations(&vg, orientation);
                    
                    int idx = 0;
                    vg.for_each_handle([&](const handle_t& handle) {
                        REQUIRE(vg.get_sequence(handle) == sequence_by_id[vg.get_id(handle)]);
                        idx++;
                    });
                }
            }
        }
        
        TEST_CASE("is_acyclic can detect cyclic graphs", "[algorithms][cycles]") {
            
            SECTION("is_acyclic works on a graph with one node") {
                
                VG vg1;
                
                // empty graph is acyclic
                REQUIRE(algorithms::is_acyclic(&vg1));
                
                handle_t n1 = vg1.create_handle("GATTACA");
                
                // no edges, still acyclic
                REQUIRE(algorithms::is_acyclic(&vg1));
                
                VG vg2 = vg1;
                VG vg3 = vg1;
                
                handle_t n2 = vg2.get_handle(vg1.get_id(n1), false);
                handle_t n3 = vg3.get_handle(vg1.get_id(n1), false);
                
                // add each type of self edge
                vg1.create_edge(n1, n1);
                vg2.create_edge(n2, vg2.flip(n2));
                vg3.create_edge(vg3.flip(n3), n3);
                
                // all of these are now cyclic
                REQUIRE(!algorithms::is_acyclic(&vg1));
                REQUIRE(!algorithms::is_acyclic(&vg2));
                REQUIRE(!algorithms::is_acyclic(&vg3));
                
            }
            
            SECTION("is_acyclic works on a graph with multiple nodes") {
                
                VG vg;
                
                handle_t n1 = vg.create_handle("GATTACA");
                handle_t n2 = vg.create_handle("GATTACA");
                handle_t n3 = vg.create_handle("GATTACA");
                handle_t n4 = vg.create_handle("GATTACA");
                handle_t n5 = vg.create_handle("GATTACA");
                
                vg.create_edge(n1, n2);
                vg.create_edge(n1, vg.flip(n3));
                vg.create_edge(vg.flip(n3), n4);
                vg.create_edge(n2, n4);
                vg.create_edge(n2, n5);
                
                // the base is a DAG
                REQUIRE(algorithms::is_acyclic(&vg));
                
                // add a non-reversing cycle
                {
                    VG cyclic = vg;
                    cyclic.create_edge(cyclic.get_handle(vg.get_id(n5), false), cyclic.get_handle(vg.get_id(n2), false));
                    REQUIRE(!algorithms::is_acyclic(&cyclic));
                }
                
                // add a reversing cycle
                {
                    VG cyclic = vg;
                    cyclic.create_edge(cyclic.get_handle(vg.get_id(n5), false), cyclic.get_handle(vg.get_id(n3), false));
                    REQUIRE(!algorithms::is_acyclic(&cyclic));
                }
            }
        }
        
        TEST_CASE("split_strands() should properly split the forward and reverse strands", "[vg][split]") {
            
            VG graph;
            handle_t n1 = graph.create_handle("ATA", 1);
            handle_t n2 = graph.create_handle("CT", 2);
            handle_t n3 = graph.create_handle("TGA", 3);
            
            graph.create_edge(n1, n2);
            graph.create_edge(graph.flip(n3), graph.flip(n2));
            graph.create_edge(n1, graph.flip(n2));
            graph.create_edge(graph.flip(n2), n3);
            
            VG split;
            unordered_map<id_t, pair<id_t, bool> > node_translation = algorithms::split_strands(&graph, &split);
                    
            Graph& g = split.graph;
            
            REQUIRE(g.node_size() == 6);
            REQUIRE(g.edge_size() == 8);
            
            int64_t node_1 = 0;
            int64_t node_2 = 0;
            int64_t node_3 = 0;
            int64_t node_4 = 0;
            int64_t node_5 = 0;
            int64_t node_6 = 0;
            
            for (int i = 0; i < g.node_size(); i++) {
                const Node& n = g.node(i);
                int64_t orig_id = node_translation[n.id()].first;
                bool flipped =  node_translation[n.id()].second;
                if (orig_id == 1 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                    node_1 = n.id();
                }
                else if (orig_id == 1 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                    node_2 = n.id();
                }
                else if (orig_id == 2 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                    node_3 = n.id();
                }
                else if (orig_id == 2 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                    node_4 = n.id();
                }
                else if (orig_id == 3 && !flipped && n.sequence() == graph.get_node(orig_id)->sequence()) {
                    node_5 = n.id();
                }
                else if (orig_id == 3 && flipped && n.sequence() == reverse_complement(graph.get_node(orig_id)->sequence())) {
                    node_6 = n.id();
                }
            }
            
            REQUIRE(node_1 != 0);
            REQUIRE(node_2 != 0);
            REQUIRE(node_3 != 0);
            REQUIRE(node_4 != 0);
            REQUIRE(node_5 != 0);
            REQUIRE(node_6 != 0);
            
            bool found_edge_1 = false;
            bool found_edge_2 = false;
            bool found_edge_3 = false;
            bool found_edge_4 = false;
            bool found_edge_5 = false;
            bool found_edge_6 = false;
            bool found_edge_7 = false;
            bool found_edge_8 = false;
            
            for (int i = 0; i < g.edge_size(); i++) {
                const Edge& e = g.edge(i);
                if ((e.from() == node_1 && e.to() == node_3 && !e.from_start() && !e.to_end()) ||
                    (e.from() == node_3 && e.to() == node_1 && e.from_start() && e.to_end())) {
                    found_edge_1 = true;
                }
                else if ((e.from() == node_1 && e.to() == node_4 && !e.from_start() && !e.to_end()) ||
                         (e.from() == node_4 && e.to() == node_1 && e.from_start() && e.to_end())) {
                    found_edge_2 = true;
                }
                else if ((e.from() == node_6 && e.to() == node_3 && !e.from_start() && !e.to_end()) ||
                         (e.from() == node_3 && e.to() == node_6 && e.from_start() && e.to_end())) {
                    found_edge_3 = true;
                }
                else if ((e.from() == node_6 && e.to() == node_4 && !e.from_start() && !e.to_end()) ||
                         (e.from() == node_4 && e.to() == node_6 && e.from_start() && e.to_end())) {
                    found_edge_4 = true;
                }
                else if ((e.from() == node_3 && e.to() == node_5 && !e.from_start() && !e.to_end()) ||
                         (e.from() == node_5 && e.to() == node_3 && e.from_start() && e.to_end())) {
                    found_edge_5 = true;
                }
                else if ((e.from() == node_3 && e.to() == node_2 && !e.from_start() && !e.to_end()) ||
                         (e.from() == node_2 && e.to() == node_3 && e.from_start() && e.to_end())) {
                    found_edge_6 = true;
                }
                else if ((e.from() == node_4 && e.to() == node_5 && !e.from_start() && !e.to_end()) ||
                         (e.from() == node_5 && e.to() == node_4 && e.from_start() && e.to_end())) {
                    found_edge_7 = true;
                }
                else if ((e.from() == node_4 && e.to() == node_2 && !e.from_start() && !e.to_end()) ||
                         (e.from() == node_2 && e.to() == node_4 && e.from_start() && e.to_end())) {
                    found_edge_8 = true;
                }
            }
            
            REQUIRE(found_edge_1);
            REQUIRE(found_edge_2);
            REQUIRE(found_edge_3);
            REQUIRE(found_edge_4);
            REQUIRE(found_edge_5);
            REQUIRE(found_edge_6);
            REQUIRE(found_edge_7);
            REQUIRE(found_edge_8);
        }
        
        TEST_CASE("count_walks() can count the source-to-sink walks in a DAG", "[algorithms][walks]") {
            
            VG vg;
            
            handle_t n0 = vg.create_handle("GATTACA");
            handle_t n1 = vg.create_handle("GATTACA");
            handle_t n2 = vg.create_handle("GATTACA");
            handle_t n3 = vg.create_handle("GATTACA");
            handle_t n4 = vg.create_handle("GATTACA");
            handle_t n5 = vg.create_handle("GATTACA");
            handle_t n6 = vg.create_handle("GATTACA");
            handle_t n7 = vg.create_handle("GATTACA");
            
            vg.create_edge(n0, n5);
            vg.create_edge(n1, n2);
            vg.create_edge(n1, n3);
            vg.create_edge(n3, n4);
            vg.create_edge(n2, n4);
            vg.create_edge(n2, n5);
            vg.create_edge(n4, n6);
            vg.create_edge(n4, n7);
            
            REQUIRE(algorithms::count_walks(&vg) == 6);
        }
        
        TEST_CASE("strongly_connected_components() works in a connected graph with reversing edges and no tips", "[algorithms][components]") {
        
            VG graph;
                    
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            Node* n4 = graph.create_node("CTGA");
            Node* n5 = graph.create_node("GCA");
            
            Edge* e1 = graph.create_edge(n1, n2);
            Edge* e2 = graph.create_edge(n1, n2, true, false);
            Edge* e3 = graph.create_edge(n2, n3);
            Edge* e4 = graph.create_edge(n3, n4);
            Edge* e5 = graph.create_edge(n3, n5);
            Edge* e6 = graph.create_edge(n4, n5);
            Edge* e7 = graph.create_edge(n5, n3, false, true);
            
            auto components = algorithms::strongly_connected_components(&graph);
            
            REQUIRE(components.size() == 1);
            REQUIRE(components.at(0).size() == 5);
                
        
        }
        
        // a heuristic that can be manually provided using a map
        struct TestDistHeuristic {
        public:
            TestDistHeuristic(const HandleGraph* graph, const unordered_map<handle_t, int64_t>& heuristic_values) :
            heuristic_values(heuristic_values), graph(graph) { }
            
            unordered_map<handle_t, int64_t> heuristic_values;
            const HandleGraph* graph;
            
            int64_t operator()(const pos_t& pos_1, const pos_t& pos_2) const {
                return heuristic_values.at(graph->get_handle(id(pos_1), is_rev(pos_1))) - offset(pos_1) + offset(pos_2);
            }
        };
        
        TEST_CASE("A* search can detect minimum length paths", "[algorithms][a-star]") {
            
            SECTION("A* runs on a simple linear graph") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                
                graph.create_edge(n1, n2);
                
                pos_t pos_1 = make_pos_t(graph.get_id(n1), false, 0);
                pos_t pos_2 = make_pos_t(graph.get_id(n2), false, 4);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 3},
                    {n2, 0}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic);
                
                REQUIRE(path.size() == 2);
                REQUIRE(path[0] == n1);
                REQUIRE(path[1] == n2);
            }
            
            SECTION("A* finds a minimum distance path in a DAG") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                
                pos_t pos_1 = make_pos_t(graph.get_id(n1), false, 2);
                pos_t pos_2 = make_pos_t(graph.get_id(n6), false, 1);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 5},
                    {n2, 6},
                    {n3, 2},
                    {n4, 2},
                    {n5, 1},
                    {n6, 0}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic);
                
                REQUIRE(path.size() == 4);
                REQUIRE(path[0] == n1);
                REQUIRE(path[1] == n3);
                REQUIRE(path[2] == n5);
                REQUIRE(path[3] == n6);
            }
            
            SECTION("A* finds a minimum distance path in a DAG with an uninformative heuristic") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                
                pos_t pos_1 = make_pos_t(graph.get_id(n1), false, 2);
                pos_t pos_2 = make_pos_t(graph.get_id(n6), false, 1);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 0},
                    {n2, 0},
                    {n3, 0},
                    {n4, 0},
                    {n5, 0},
                    {n6, 0}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic);
                
                REQUIRE(path.size() == 4);
                REQUIRE(path[0] == n1);
                REQUIRE(path[1] == n3);
                REQUIRE(path[2] == n5);
                REQUIRE(path[3] == n6);
            }
            
            SECTION("A* finds a minimum distance path in a cyclic graph") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, graph.flip(n1));
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                graph.create_edge(n5, graph.flip(n6));
                
                pos_t pos_1 = make_pos_t(graph.get_id(n1), false, 2);
                pos_t pos_2 = make_pos_t(graph.get_id(n1), true, 1);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 3},
                    {n2, 5},
                    {n3, 0},
                    {n4, 1000},
                    {n5, 6},
                    {n6, 3},
                    {graph.flip(n1), 0},
                    {graph.flip(n2), 1000},
                    {graph.flip(n3), 1000},
                    {graph.flip(n4), 2},
                    {graph.flip(n5), 2},
                    {graph.flip(n6), 4}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic);
                
                REQUIRE(path.size() == 3);
                REQUIRE(path[0] == n1);
                REQUIRE(path[1] == n3);
                REQUIRE(path[2] == graph.flip(n1));
            }
            
            SECTION("A* finds a minimum distance path on the same node") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, graph.flip(n1));
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                graph.create_edge(n5, graph.flip(n6));
                
                pos_t pos_1 = make_pos_t(graph.get_id(n6), true, 1);
                pos_t pos_2 = make_pos_t(graph.get_id(n6), true, 2);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 0},
                    {n2, 0},
                    {n3, 0},
                    {n4, 0},
                    {n5, 0},
                    {n6, 0},
                    {graph.flip(n1), 0},
                    {graph.flip(n2), 0},
                    {graph.flip(n3), 0},
                    {graph.flip(n4), 0},
                    {graph.flip(n5), 0},
                    {graph.flip(n6), 0}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic);
                
                REQUIRE(path.size() == 1);
                REQUIRE(path[0] == graph.flip(n6));
            }
            
            SECTION("A* will not find a minimum distance above the pruning distance") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, graph.flip(n1));
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                graph.create_edge(n5, graph.flip(n6));
                
                pos_t pos_1 = make_pos_t(graph.get_id(n1), false, 2);
                pos_t pos_2 = make_pos_t(graph.get_id(n3), false, 1);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 3},
                    {n2, 5},
                    {n3, 0},
                    {n4, 1000},
                    {n5, 6},
                    {n6, 3},
                    {graph.flip(n1), 0},
                    {graph.flip(n2), 1000},
                    {graph.flip(n3), 1000},
                    {graph.flip(n4), 2},
                    {graph.flip(n5), 2},
                    {graph.flip(n6), 4}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic, true, 2);
                
                REQUIRE(path.empty());
            }
            
            SECTION("A* will not find a minimum distance above the pruning distance on the same node") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, graph.flip(n1));
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                graph.create_edge(n5, graph.flip(n6));
                
                pos_t pos_1 = make_pos_t(graph.get_id(n4), false, 0);
                pos_t pos_2 = make_pos_t(graph.get_id(n4), false, 3);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 3},
                    {n2, 5},
                    {n3, 0},
                    {n4, 1000},
                    {n5, 6},
                    {n6, 3},
                    {graph.flip(n1), 0},
                    {graph.flip(n2), 1000},
                    {graph.flip(n3), 1000},
                    {graph.flip(n4), 2},
                    {graph.flip(n5), 2},
                    {graph.flip(n6), 4}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic, true, 2);
                
                REQUIRE(path.empty());
            }
            
            SECTION("A* finds a minimum distance path in a cyclic graph and an uninformative heuristic") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, graph.flip(n1));
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                graph.create_edge(n5, graph.flip(n6));
                
                pos_t pos_1 = make_pos_t(graph.get_id(n1), false, 2);
                pos_t pos_2 = make_pos_t(graph.get_id(n6), true, 1);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 0},
                    {n2, 0},
                    {n3, 0},
                    {n4, 0},
                    {n5, 0},
                    {n6, 0},
                    {graph.flip(n1), 0},
                    {graph.flip(n2), 0},
                    {graph.flip(n3), 0},
                    {graph.flip(n4), 0},
                    {graph.flip(n5), 0},
                    {graph.flip(n6), 0}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic);
                
                REQUIRE(path.size() == 4);
                REQUIRE(path[0] == n1);
                REQUIRE(path[1] == n3);
                REQUIRE(path[2] == n5);
                REQUIRE(path[3] == graph.flip(n6));
            }
        }
        
        TEST_CASE("A* search can find maximum length paths","[a-star][algorithms]") {
            
            SECTION("A* finds a maximum length path in a DAG") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                
                pos_t pos_1 = make_pos_t(graph.get_id(n6), true, 2);
                pos_t pos_2 = make_pos_t(graph.get_id(n1), true, 1);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 0},
                    {n2, 0},
                    {n3, 0},
                    {n4, 0},
                    {n5, 0},
                    {n6, 0},
                    {graph.flip(n1), 0},
                    {graph.flip(n2), 5},
                    {graph.flip(n3), 3},
                    {graph.flip(n4), 8},
                    {graph.flip(n5), 4},
                    {graph.flip(n6), 11}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic, false, numeric_limits<int64_t>::lowest());
                
                REQUIRE(path.size() == 4);
                REQUIRE(path[0] == graph.flip(n6));
                REQUIRE(path[1] == graph.flip(n4));
                REQUIRE(path[2] == graph.flip(n2));
                REQUIRE(path[3] == graph.flip(n1));
            }
            
            SECTION("A* finds a maximum length path in a DAG with an uninformative heuristic") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("AGGTA");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TGAC");
                handle_t n6 = graph.create_handle("CCG");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                
                pos_t pos_1 = make_pos_t(graph.get_id(n6), true, 2);
                pos_t pos_2 = make_pos_t(graph.get_id(n1), true, 1);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 100},
                    {n2, 100},
                    {n3, 100},
                    {n4, 100},
                    {n5, 100},
                    {n6, 100},
                    {graph.flip(n1), 0},
                    {graph.flip(n2), 100},
                    {graph.flip(n3), 100},
                    {graph.flip(n4), 100},
                    {graph.flip(n5), 100},
                    {graph.flip(n6), 100}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic, false, numeric_limits<int64_t>::lowest());
                
                REQUIRE(path.size() == 4);
                REQUIRE(path[0] == graph.flip(n6));
                REQUIRE(path[1] == graph.flip(n5));
                REQUIRE(path[2] == graph.flip(n3));
                REQUIRE(path[3] == graph.flip(n1));
            }
        }
        
        TEST_CASE("A* search works with non-monotonic heuristics","[a-star][algorithms]") {
            
            SECTION("A* finds shortest path with a non-monotonic heuristic") {
                
                VG graph;
                
                handle_t n1 = graph.create_handle("GGGA");
                handle_t n2 = graph.create_handle("TACC");
                handle_t n3 = graph.create_handle("A");
                handle_t n4 = graph.create_handle("ATG");
                handle_t n5 = graph.create_handle("TG");
                handle_t n6 = graph.create_handle("CCG");
                handle_t n7 = graph.create_handle("C");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n7);
                graph.create_edge(n7, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n6);
                graph.create_edge(n5, n6);
                
                pos_t pos_1 = make_pos_t(graph.get_id(n1), false, 2);
                pos_t pos_2 = make_pos_t(graph.get_id(n6), false, 1);
                
                unordered_map<handle_t, int64_t> heuristic_values{
                    {n1, 4},
                    {n2, 0},
                    {n3, 1},
                    {n4, 1},
                    {n5, 2},
                    {n6, 0},
                    {n7, 5}
                };
                TestDistHeuristic heuristic(&graph, heuristic_values);
                
                vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic, true, numeric_limits<int64_t>::max(), false);
                
                REQUIRE(path.size() == 5);
                REQUIRE(path[0] == n1);
                REQUIRE(path[1] == n7);
                REQUIRE(path[2] == n3);
                REQUIRE(path[3] == n5);
                REQUIRE(path[4] == n6);
            }
            
        }
        
        TEST_CASE("A* search works on random graphs","[a-star][algorithms]") {
            
            SECTION("A* finds shortest path on a random graph with a perfect heuristic") {
                
                int64_t seq_size = 100;
                int64_t avg_struct_var_len = 6;
                int64_t var_count = 10;
                
                size_t num_graphs = 10;
                size_t num_trials_per_graph = 10;
                
                for (size_t graph_iter = 0; graph_iter < num_graphs; graph_iter++) {
                    
                    VG graph = randomGraph(seq_size, avg_struct_var_len, var_count);
                    
                    size_t total_seq_len = 0;
                    vector<handle_t> all_handles;
                    graph.for_each_handle([&](const handle_t& handle) {
                        all_handles.push_back(handle);
                        total_seq_len += graph.get_length(handle) * 2;
                    });
                    random_device rd;
                    default_random_engine gen(rd());
                    
                    function<pos_t(void)> random_pos = [&](void) {
                        handle_t h = all_handles[uniform_int_distribution<int>(0, all_handles.size() - 1)(gen)];
                        return make_pos_t(graph.get_id(h),
                                          uniform_int_distribution<int>(0, 1)(gen),
                                          uniform_int_distribution<size_t>(0, graph.get_length(h) - 1)(gen));
                    };
                    
                    for (size_t pos_iter = 0; pos_iter < num_trials_per_graph; pos_iter++) {
                        
                        pos_t pos_1 = random_pos();
                        pos_t pos_2 = random_pos();
                        // shortest path calculations get complicated in the same node unreachable case, so
                        // we just forbit it for simplicity here
                        while (id(pos_2) == id(pos_1) && is_rev(pos_2) == is_rev(pos_1) && offset(pos_2) < offset(pos_1)) {
                            pos_1 = random_pos();
                            pos_2 = random_pos();
                        }
                        
                        handle_t h1 = graph.get_handle(id(pos_1), is_rev(pos_1));
                        handle_t h2 = graph.get_handle(id(pos_2), is_rev(pos_2));
                        
                        unordered_map<handle_t, size_t> shortest_paths = algorithms::find_shortest_paths(&graph, h2, true);
                        
                        unordered_map<handle_t, int64_t> heuristic_values;
                        for (const auto& handle : all_handles) {
                            auto flipped = graph.flip(handle);
                            if (shortest_paths.count(handle)) {
                                heuristic_values[handle] = shortest_paths[handle] + graph.get_length(handle);
                            }
                            else {
                                heuristic_values[handle] = total_seq_len + 1;
                            }
                            if (shortest_paths.count(flipped)) {
                                heuristic_values[flipped] = shortest_paths[flipped] + graph.get_length(flipped);
                            }
                            else {
                                heuristic_values[flipped] = total_seq_len + 1;
                            }
                        }
                        TestDistHeuristic heuristic(&graph, heuristic_values);
                        
                        vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic);

                        size_t path_length = 0;
                        // find_shortest_path does not include the end nodes
                        for (size_t i = 1; i + 1 < path.size(); i++) {
                            path_length += graph.get_length(path[i]);
                        }
                        
                        if (shortest_paths.count(h1)) {
                            if (path_length != shortest_paths[h1]) {
                                cerr << "A* failed to find shortest path from " << pos_1 << " to " << pos_2 << " in the following graph:" << endl;
                                cerr << pb2json(graph.graph) << endl;
                            }
                            
                            REQUIRE(path_length == shortest_paths[h1]);
                        }
                        else {
                            if (!path.empty()) {
                                cerr << "A* failed to find identify unreachable path " << pos_1 << " to " << pos_2 << " in the following graph:" << endl;
                                cerr << pb2json(graph.graph) << endl;
                            }
                            
                            REQUIRE(path.empty());
                        }
                    }
                }
            }
            
            SECTION("A* finds shortest path on a random graph with an uninformative heuristic") {
                
                int64_t seq_size = 100;
                int64_t avg_struct_var_len = 6;
                int64_t var_count = 10;
                
                size_t num_graphs = 10;
                size_t num_trials_per_graph = 10;
                
                for (size_t graph_iter = 0; graph_iter < num_graphs; graph_iter++) {
                    VG graph = randomGraph(seq_size, avg_struct_var_len, var_count);
                    
                    size_t total_seq_len = 0;
                    vector<handle_t> all_handles;
                    graph.for_each_handle([&](const handle_t& handle) {
                        all_handles.push_back(handle);
                        total_seq_len += graph.get_length(handle) * 2;
                    });
                    random_device rd;
                    default_random_engine gen(rd());
                    
                    function<pos_t(void)> random_pos = [&](void) {
                        handle_t h = all_handles[uniform_int_distribution<int>(0, all_handles.size() - 1)(gen)];
                        return make_pos_t(graph.get_id(h),
                                          uniform_int_distribution<int>(0, 1)(gen),
                                          uniform_int_distribution<size_t>(0, graph.get_length(h) - 1)(gen));
                    };
                    
                    for (size_t pos_iter = 0; pos_iter < num_trials_per_graph; pos_iter++) {
                        
                        pos_t pos_1 = random_pos();
                        pos_t pos_2 = random_pos();
                        // shortest path calculations get complicated in the same node unreachable case, so
                        // we just forbit it for simplicity here
                        while (id(pos_2) == id(pos_1) && is_rev(pos_2) == is_rev(pos_1) && offset(pos_2) < offset(pos_1)) {
                            pos_1 = random_pos();
                            pos_2 = random_pos();
                        }
                        
                        handle_t h1 = graph.get_handle(id(pos_1), is_rev(pos_1));
                        handle_t h2 = graph.get_handle(id(pos_2), is_rev(pos_2));
                        
                        unordered_map<handle_t, size_t> shortest_paths = algorithms::find_shortest_paths(&graph, h2, true);
                        
                        unordered_map<handle_t, int64_t> heuristic_values;
                        for (const auto& handle : all_handles) {
                            auto flipped = graph.flip(handle);
                            heuristic_values[handle] = 0;
                            heuristic_values[flipped] = 0;
                        }
                        TestDistHeuristic heuristic(&graph, heuristic_values);
                        
                        vector<handle_t> path = algorithms::a_star(&graph, pos_1, pos_2, heuristic);
                        
                        size_t path_length = 0;
                        // find_shortest_path does not include the end nodes
                        for (size_t i = 1; i + 1 < path.size(); i++) {
                            path_length += graph.get_length(path[i]);
                        }
                        
                        if (shortest_paths.count(h1)) {
                            if (path_length != shortest_paths[h1]) {
                                cerr << "A* failed to find shortest path from " << pos_1 << " to " << pos_2 << " in the following graph:" << endl;
                                cerr << pb2json(graph.graph) << endl;
                            }
                            
                            REQUIRE(path_length == shortest_paths[h1]);
                        }
                        else {
                            if (!path.empty()) {
                                cerr << "A* failed to find identify unreachable path " << pos_1 << " to " << pos_2 << " in the following graph:" << endl;
                                cerr << pb2json(graph.graph) << endl;
                            }
                            
                            REQUIRE(path.empty());
                        }
                    }
                }
            }
        }
    }
}
