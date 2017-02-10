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
#include "algorithms/vg_algorithms.hpp"
#include "vg.hpp"
#include "json2pb.h"


using namespace google::protobuf;

namespace vg {
    
    namespace unittest {
        
        TEST_CASE( "Graph extraction produces expected results in an acyclic graph",
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
                REQUIRE( g.node_size() == 0 );
                REQUIRE( g.edge_size() == 0 );
            }
            
            SECTION( "Graph extraction returns an empty graph if there are no paths under the max length and the positions are on different nodes" ) {
                
                pos_t pos_1 = make_pos_t(n5->id(), true, 0);
                pos_t pos_2 = make_pos_t(n0->id(), true, 2);
                
                int64_t max_len = 5;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
                REQUIRE( g.node_size() == 0 );
                REQUIRE( g.edge_size() == 0 );
            }
            
            SECTION( "Graph extraction can detect a one-node path at the maximum distance" ) {
                
                pos_t pos_1 = make_pos_t(n5->id(), true, 0);
                pos_t pos_2 = make_pos_t(n5->id(), true, 4);
                
                int64_t max_len = 4;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
                cerr << pb2json(g) << endl;
                
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
        TEST_CASE( "Graph extraction produces expected results in cyclic graphs",
                  "[algorithms]" ) {
            
            SECTION( "Graph extraction can detect a loop onto the same node" ) {
                
                VG vg;
                
                Node* n0 = vg.create_node("ATATAAC");
                
                vg.create_edge(n0, n0);
                
                pos_t pos_1 = make_pos_t(n0->id(), false, 4);
                pos_t pos_2 = make_pos_t(n0->id(), false, 3);
                
                int64_t max_len = 10;
                
                Graph g;
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
                
                auto trans = algorithms::extract_connecting_graph(vg, g, max_len, pos_1, pos_2);
                
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
        }
    }
}
