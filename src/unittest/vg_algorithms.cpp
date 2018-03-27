//
//  vg_algorithms.cpp
//  
//
//  Created by Jordan Eizenga on 2/7/17.
//
//

#include <stdio.h>
#include <set>
#include "catch.hpp"
#include "algorithms/extract_connecting_graph.hpp"
#include "algorithms/extract_containing_graph.hpp"
#include "algorithms/extract_extending_graph.hpp"
#include "algorithms/topological_sort.hpp"
#include "algorithms/weakly_connected_components.hpp"
#include "algorithms/distance_to_head.hpp"
#include "algorithms/distance_to_tail.hpp"
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
                
                pos_t pos_1 = make_pos_t(n1->id(), false, 0);
                pos_t pos_2 = make_pos_t(n1->id(), false, 3);
                
                int64_t max_len = 5;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
                REQUIRE( g.node_size() == 1 );
                REQUIRE( g.edge_size() == 0 );
                REQUIRE( g.node(0).id() == n1->id() );
                REQUIRE( g.node(0).sequence() == "TG" );
            }
            
            SECTION( "Graph extraction can extract a section between two nodes" ) {
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 1);
                pos_t pos_2 = make_pos_t(n1->id(), false, 2);
                
                int64_t max_len = 5;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
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
                
                pos_t pos_1 = make_pos_t(n2->id(), false, 1);
                pos_t pos_2 = make_pos_t(n5->id(), false, 3);
                
                int64_t max_len = 10;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
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
                
                pos_t pos_1 = make_pos_t(n2->id(), false, 0);
                pos_t pos_2 = make_pos_t(n4->id(), false, 1);
                
                int64_t max_len = 10;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
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
                
                pos_t pos_1 = make_pos_t(n5->id(), true, 0);
                pos_t pos_2 = make_pos_t(n5->id(), true, 4);
                
                int64_t max_len = 2;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
                REQUIRE( g.node_size() == 0 );
                REQUIRE( g.edge_size() == 0 );
            }
            
            SECTION( "Graph extraction returns an empty graph if there are no paths under the max length and the positions are on different nodes" ) {
                
                pos_t pos_1 = make_pos_t(n5->id(), true, 0);
                pos_t pos_2 = make_pos_t(n0->id(), true, 2);
                
                int64_t max_len = 5;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
                REQUIRE( g.node_size() == 0 );
                REQUIRE( g.edge_size() == 0 );
            }
            
            SECTION( "Graph extraction can detect a one-node path at the maximum distance" ) {
                
                pos_t pos_1 = make_pos_t(n5->id(), true, 0);
                pos_t pos_2 = make_pos_t(n5->id(), true, 4);
                
                int64_t max_len = 4;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
                REQUIRE( g.node_size() == 1 );
                REQUIRE( g.edge_size() == 0 );
                REQUIRE( g.node(0).id() == n5->id() );
                REQUIRE( g.node(0).sequence() == "ATA" );
            }
            
            SECTION( "Graph extraction can detect a two-node path at the maximum distance" ) {
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 0);
                pos_t pos_2 = make_pos_t(n1->id(), false, 3);
                
                int64_t max_len = 6;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
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
                
                pos_t pos_1 = make_pos_t(n2->id(), false, 0);
                pos_t pos_2 = make_pos_t(n5->id(), false, 4);
                
                int64_t max_len = 9;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
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
                
                pos_t pos_1 = make_pos_t(n9->id(), true, 6);
                pos_t pos_2 = make_pos_t(n5->id(), true, 0);
                
                int64_t max_len = 5;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
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
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 4);
                pos_t pos_2 = make_pos_t(n0->id(), false, 3);
                
                int64_t max_len = 10;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
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
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 4);
                pos_t pos_2 = make_pos_t(n0->id(), true, 3);
                
                int64_t max_len = 10;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, false, true, true, true);
                
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
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 2);
                pos_t pos_2 = make_pos_t(n1->id(), false, 1);
                
                int64_t max_len = 15;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, true);
                
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
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 2);
                pos_t pos_2 = make_pos_t(n1->id(), false, 1);
                
                int64_t max_len = 20;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, true);
                
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
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 2);
                pos_t pos_2 = make_pos_t(n1->id(), false, 1);
                
                int64_t max_len = 30;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, true);
                
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
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 1);
                pos_t pos_2 = make_pos_t(n0->id(), false, 3);
                
                int64_t max_len = 40;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, true);
                
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
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 3);
                pos_t pos_2 = make_pos_t(n0->id(), false, 1);
                
                int64_t max_len = 40;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, true);
                
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
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 3);
                pos_t pos_2 = make_pos_t(n0->id(), true, 3);
                
                int64_t max_len = 40;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(&vg, g, max_len, pos_1, pos_2, false, true);
                
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
        
        TEST_CASE( "Connecting graph extraction options perform as expected", "[algorithms]" ) {
            VG vg;
            
            Node* n0 = vg.create_node("CGA");
            Node* n1 = vg.create_node("TTGG");
            Node* n2 = vg.create_node("GT"); // a tip
            Node* n3 = vg.create_node("ATG"); // a self-contained cycle
            Node* n4 = vg.create_node("TGAG"); // a longer path
            Node* n5 = vg.create_node("CA"); // a shorter path
            Node* n6 = vg.create_node("T");
            
            vg.create_edge(n1, n0, true, true); // a doubly reversing edge to keep it interesting
            vg.create_edge(n1, n2);
            vg.create_edge(n1, n3);
            vg.create_edge(n1, n4, false, true);
            vg.create_edge(n1, n5);
            vg.create_edge(n3, n3, true, true);
            vg.create_edge(n4, n6, true, false);
            vg.create_edge(n5, n6);
            
            SECTION("Cutting behind end positions works on a single node path") {
                
                int64_t max_len = 40;
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, max_len,
                                                                      make_pos_t(n0->id(), false, 0),
                                                                      make_pos_t(n0->id(), false, 2),
                                                                      true);
                    
                    REQUIRE(g.node_size() == 1);
                    REQUIRE(g.edge_size() == 0);
                    REQUIRE(g.node(0).sequence() == "CGA");
                }

                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, max_len,
                                                                      make_pos_t(n0->id(), true, 0),
                                                                      make_pos_t(n0->id(), true, 2),
                                                                      true);
                    
                    REQUIRE(g.node_size() == 1);
                    REQUIRE(g.edge_size() == 0);
                    REQUIRE(g.node(0).sequence() == "CGA");
                }
            }
            
            SECTION("Cutting behind end positions works on a two node path") {
                
                int64_t max_len = 40;
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, max_len,
                                                                      make_pos_t(n0->id(), false, 0),
                                                                      make_pos_t(n1->id(), false, 2),
                                                                      true);
                    
                    REQUIRE(g.node_size() == 2);
                    REQUIRE(g.edge_size() == 1);
                    
                    bool found_node_0 = false;
                    bool found_node_1 = false;
                    
                    for (int i = 0; i < g.node_size(); i++) {
                        const Node& n = g.node(i);
                        if (n.sequence() == "CGA") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "TTG") {
                            found_node_1 = true;
                        }
                    }
                    
                    REQUIRE(found_node_0);
                    REQUIRE(found_node_1);
                    
                    bool found_edge_0 = false;
                    const Edge& e = g.edge(0);
                    if ((e.from() == n0->id() && !e.from_start() && e.to() == n1->id() && !e.to_end())
                        ||(e.from() == n1->id() && e.from_start() && e.to() == n0->id() && e.to_end())){
                        found_edge_0 = true;
                    }
                    REQUIRE(found_edge_0);
                }
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, max_len,
                                                                      make_pos_t(n1->id(), true, 1),
                                                                      make_pos_t(n0->id(), true, 2),
                                                                      true);
                    
                    REQUIRE(g.node_size() == 2);
                    REQUIRE(g.edge_size() == 1);
                    
                    bool found_node_0 = false;
                    bool found_node_1 = false;
                    
                    for (int i = 0; i < g.node_size(); i++) {
                        const Node& n = g.node(i);
                        if (n.sequence() == "CGA") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "TTG") {
                            found_node_1 = true;
                        }
                    }
                    
                    REQUIRE(found_node_0);
                    REQUIRE(found_node_1);
                    
                    bool found_edge_0 = false;
                    const Edge& e = g.edge(0);
                    if ((e.from() == n0->id() && !e.from_start() && e.to() == n1->id() && !e.to_end())
                        ||(e.from() == n1->id() && e.from_start() && e.to() == n0->id() && e.to_end())){
                        found_edge_0 = true;
                    }
                    REQUIRE(found_edge_0);
                }
            }
            
            SECTION("Distance filtering works when cutting behind end positions") {
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, 0,
                                                                      make_pos_t(n0->id(), true, 0),
                                                                      make_pos_t(n0->id(), true, 0),
                                                                      true);
                    
                    REQUIRE(g.node_size() == 1);
                    REQUIRE(g.edge_size() == 0);
                    REQUIRE(g.node(0).sequence() == "A");
                }
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, 5,
                                                                      make_pos_t(n0->id(), false, 0),
                                                                      make_pos_t(n1->id(), false, 2),
                                                                      true);
                    
                    REQUIRE(g.node_size() == 2);
                    REQUIRE(g.edge_size() == 1);
                    
                    bool found_node_0 = false;
                    bool found_node_1 = false;
                    
                    for (int i = 0; i < g.node_size(); i++) {
                        const Node& n = g.node(i);
                        if (n.sequence() == "CGA") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "TTG") {
                            found_node_1 = true;
                        }
                    }
                    
                    REQUIRE(found_node_0);
                    REQUIRE(found_node_1);
                    
                    bool found_edge_0 = false;
                    const Edge& e = g.edge(0);
                    if ((e.from() == n0->id() && !e.from_start() && e.to() == n1->id() && !e.to_end())
                        ||(e.from() == n1->id() && e.from_start() && e.to() == n0->id() && e.to_end())){
                        found_edge_0 = true;
                    }
                    REQUIRE(found_edge_0);
                }
                
                
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, 4,
                                                                      make_pos_t(n0->id(), false, 0),
                                                                      make_pos_t(n1->id(), false, 2),
                                                                      true);
                    
                    REQUIRE(g.node_size() == 0);
                    REQUIRE(g.edge_size() == 0);
                }
            }
            
            SECTION("Pruning options produce expected output") {
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, 10,
                                                                      make_pos_t(n1->id(), false, 3),
                                                                      make_pos_t(n6->id(), false, 0),
                                                                      true, false);
                    
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
                        if (((e.from() == n1->id() && !e.from_start() && e.to() == n2->id() && !e.to_end())
                             ||(e.from() == n2->id() && e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_0 = true;
                        }
                        else if (((e.from() == n1->id() && !e.from_start() && e.to() == n3->id() && !e.to_end())
                                  ||(e.from() == n3->id() && e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_1 = true;
                        }
                        else if (((e.from() == n1->id() && !e.from_start() && e.to() == n4->id() && e.to_end())
                                  ||(e.from() == n4->id() && !e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_2 = true;
                        }
                        else if (((e.from() == n1->id() && !e.from_start() && e.to() == n5->id() && !e.to_end())
                                  ||(e.from() == n5->id() && e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_3 = true;
                        }
                        else if (((e.from() == n3->id() && !e.from_start() && e.to() == n3->id() && !e.to_end())
                                  ||(e.from() == n3->id() && e.from_start() && e.to() == n3->id() && e.to_end()))) {
                            found_edge_4 = true;
                        }
                        else if (((e.from() == n4->id() && e.from_start() && e.to() == n6->id() && !e.to_end())
                                  ||(e.from() == n6->id() && e.from_start() && e.to() == n4->id() && !e.to_end()))) {
                            found_edge_5 = true;
                        }
                        else if (((e.from() == n5->id() && !e.from_start() && e.to() == n6->id() && !e.to_end())
                                  ||(e.from() == n6->id() && e.from_start() && e.to() == n5->id() && e.to_end()))) {
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
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, 10,
                                                                      make_pos_t(n1->id(), false, 3),
                                                                      make_pos_t(n6->id(), false, 0),
                                                                      true, false, true);
                    
                    REQUIRE(g.node_size() == 5);
                    REQUIRE(g.edge_size() == 6);
                    
                    bool found_node_0 = false;
                    bool found_node_1 = false;
                    bool found_node_2 = false;
                    bool found_node_3 = false;
                    bool found_node_4 = false;
                    
                    for (int i = 0; i < g.node_size(); i++) {
                        const Node& n = g.node(i);
                        if (n.sequence() == "G") {
                            found_node_0 = true;
                        }
                        else if (n.sequence() == "T") {
                            found_node_1 = true;
                        }
                        else if (n.sequence() == "ATG") {
                            found_node_2 = true;
                        }
                        else if (n.sequence() == "TGAG") {
                            found_node_3 = true;
                        }
                        else if (n.sequence() == "CA") {
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
                    
                    for (int i = 0; i < g.edge_size(); i++) {
                        const Edge& e = g.edge(i);
                         if (((e.from() == n1->id() && !e.from_start() && e.to() == n3->id() && !e.to_end())
                                  ||(e.from() == n3->id() && e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_0 = true;
                        }
                        else if (((e.from() == n1->id() && !e.from_start() && e.to() == n4->id() && e.to_end())
                                  ||(e.from() == n4->id() && !e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_1 = true;
                        }
                        else if (((e.from() == n1->id() && !e.from_start() && e.to() == n5->id() && !e.to_end())
                                  ||(e.from() == n5->id() && e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_2 = true;
                        }
                        else if (((e.from() == n3->id() && !e.from_start() && e.to() == n3->id() && !e.to_end())
                                  ||(e.from() == n3->id() && e.from_start() && e.to() == n3->id() && e.to_end()))) {
                            found_edge_3 = true;
                        }
                        else if (((e.from() == n4->id() && e.from_start() && e.to() == n6->id() && !e.to_end())
                                  ||(e.from() == n6->id() && e.from_start() && e.to() == n4->id() && !e.to_end()))) {
                            found_edge_4 = true;
                        }
                        else if (((e.from() == n5->id() && !e.from_start() && e.to() == n6->id() && !e.to_end())
                                  ||(e.from() == n6->id() && e.from_start() && e.to() == n5->id() && e.to_end()))) {
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
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, 10,
                                                                      make_pos_t(n1->id(), false, 3),
                                                                      make_pos_t(n6->id(), false, 0),
                                                                      true, false, true, true);
                    
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
                        if (((e.from() == n1->id() && !e.from_start() && e.to() == n4->id() && e.to_end())
                                  ||(e.from() == n4->id() && !e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_0 = true;
                        }
                        else if (((e.from() == n1->id() && !e.from_start() && e.to() == n5->id() && !e.to_end())
                                  ||(e.from() == n5->id() && e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_1 = true;
                        }
                        else if (((e.from() == n4->id() && e.from_start() && e.to() == n6->id() && !e.to_end())
                                  ||(e.from() == n6->id() && e.from_start() && e.to() == n4->id() && !e.to_end()))) {
                            found_edge_2 = true;
                        }
                        else if (((e.from() == n5->id() && !e.from_start() && e.to() == n6->id() && !e.to_end())
                                  ||(e.from() == n6->id() && e.from_start() && e.to() == n5->id() && e.to_end()))) {
                            found_edge_3 = true;
                        }
                    }
                    
                    REQUIRE(found_edge_0);
                    REQUIRE(found_edge_1);
                    REQUIRE(found_edge_2);
                    REQUIRE(found_edge_3);
                }
                
                {
                    Graph g;
                    
                    auto trans = algorithms::extract_connecting_graph(&vg, g, 4,
                                                                      make_pos_t(n1->id(), false, 3),
                                                                      make_pos_t(n6->id(), false, 0),
                                                                      true, false, true, true, true);
                    
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
                        if (((e.from() == n1->id() && !e.from_start() && e.to() == n5->id() && !e.to_end())
                                  ||(e.from() == n5->id() && e.from_start() && e.to() == n1->id() && e.to_end()))) {
                            found_edge_0 = true;
                        }
                        else if (((e.from() == n5->id() && !e.from_start() && e.to() == n6->id() && !e.to_end())
                                  ||(e.from() == n6->id() && e.from_start() && e.to() == n5->id() && e.to_end()))) {
                            found_edge_1 = true;
                        }
                    }
                    
                    REQUIRE(found_edge_0);
                    REQUIRE(found_edge_1);
                }
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
            
            Graph dest;
            
            auto max_dist = 38;
            pos_t src_pos = make_pos_t(1393981, false, 31);
            pos_t dest_pos = make_pos_t(1393958, false, 0);
            
            unordered_map<id_t, id_t> connect_trans = algorithms::extract_connecting_graph(&vg,              // DAG with split strands
                                                                                           dest,             // graph to extract into
                                                                                           max_dist,         // longest distance necessary
                                                                                           src_pos,          // end of earlier match
                                                                                           dest_pos,         // beginning of later match
                                                                                           false,            // do not extract the end positions in the matches
                                                                                           false,            // do not bother finding all cycles (it's a DAG)
                                                                                           true,             // remove tips
                                                                                           true,             // only include nodes on connecting paths
                                                                                           true);            // enforce max distance strictly
            
            
            
            
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
                
                Graph g;
                
                size_t max_len = 3;
                vector<pos_t> positions{make_pos_t(n0->id(), false, 2), make_pos_t(n5->id(), true, 1)};
                
                algorithms::extract_containing_graph(&vg, g, positions, max_len);
                
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
                
                Graph g;
                
                vector<size_t> max_lens{2, 3};
                vector<pos_t> positions{make_pos_t(n0->id(), false, 2), make_pos_t(n5->id(), true, 1)};
                
                algorithms::extract_containing_graph(&vg, g, positions, max_lens);
                
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
                
                Graph g;
                
                vector<size_t> forward_max_lens{3, 3};
                vector<size_t> backward_max_lens{2, 3};
                vector<pos_t> positions{make_pos_t(n0->id(), true, 2), make_pos_t(n5->id(), false, 1)};
                
                algorithms::extract_containing_graph(&vg, g, positions, forward_max_lens, backward_max_lens);
                                
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
        
        TEST_CASE( "Extending graph extraction algorithm produces expected results",
                  "[algorithms]" ) {
            
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                
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
                
                Graph g;
                
                auto id_trans = algorithms::extract_extending_graph(&vg, g, max_dist, pos, search_backward, preserve_cycles);
                                
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
            
            SECTION( "algorithms::topological_sort produces a consistent total ordering and orientation" ) {
                auto handle_sort = algorithms::topological_sort(&vg);

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
            
            SECTION( "algorithms::topological_sort produces a consistent total ordering and orientation" ) {
                auto handle_sort = algorithms::topological_sort(&vg);

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

        
}
    
    

}
