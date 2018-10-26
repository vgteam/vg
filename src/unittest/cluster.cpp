#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include "json2pb.h"
#include "vg.hpp"
#include "catch.hpp"
#include "snarls.hpp"
#include "cluster.hpp"
#include "distance.hpp"
#include "genotypekit.hpp"
#include "random_graph.hpp"
#include <fstream>
#include <random>
#include <time.h> 

//#define print

namespace vg {
namespace unittest {


    TEST_CASE( "TVS for simple nested snarl",
                   "[tvs][bug]" ) {
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

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
        const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
        const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

        DistanceIndex di (&graph, &snarl_manager, 20);
        TargetValueSearch tvs(graph, new TipAnchoredMaxDistance(di), 
                               new SnarlMinDistance(di)); 
    
        
        SECTION( "Test tvs" ) {
            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos1r = make_pos_t(1, true, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos6r = make_pos_t(6, true, 0);
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);

            REQUIRE(tvs.tv_path(pos1, pos2, 3, 5).size() == 2);
            REQUIRE(tvs.tv_path(pos1, pos2, 4, 5).size() == 2);
            REQUIRE(tvs.tv_path(pos1, pos2, 9, 4).size() == 0);
            REQUIRE(tvs.tv_path(pos1, pos5, 5, 4).size() == 4);
            REQUIRE(tvs.tv_path(pos1, pos5, 6, 4).size() == 4);
            REQUIRE(tvs.tv_path(pos1, pos5, 7, 1).size() == 0);

            REQUIRE(tvs.tv_path(pos1, pos5, 8, 4).size() == 5);
            REQUIRE(tvs.tv_path(pos1, pos5, 9, 4).size() == 5);
            REQUIRE(tvs.tv_path(pos1, pos8, 4, 4).size() == 2);
            REQUIRE(tvs.tv_path(pos1, pos8, 3, 4).size() == 2);
            REQUIRE(tvs.tv_path(pos1, pos8, 7, 4).size() == 5);
            REQUIRE(tvs.tv_path(pos1, pos8, 6, 4).size() == 5);
            REQUIRE(tvs.tv_path(pos1, pos8, 9, 4).size() == 6);
            REQUIRE(tvs.tv_path(pos1, pos8, 13, 4).size() == 7);
            REQUIRE(tvs.tv_path(pos1, pos8, 12, 4).size() == 7);
            REQUIRE(tvs.tv_path(pos6, pos1, 22, 4).size() == 0);
            REQUIRE(tvs.tv_path(pos6r, pos1r, 2, 4).size() == 3);
           
        }
    }//End test case
    TEST_CASE( "TVS for loopy snarls",
                   "[tvs][bug]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCAAAAAAAAAA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("CA");
        Node* n10 = graph.create_node("C");
        Node* n11 = graph.create_node("C");
        Node* n12 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n11);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n6, n12);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n9, n2);
        Edge* e15 = graph.create_edge(n11, n9);
        Edge* e16 = graph.create_edge(n12, n8);
        Edge* e17 = graph.create_edge(n9, n9, false, true);
        Edge* e18 = graph.create_edge(n2, n2, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
        const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
        const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

        DistanceIndex di (&graph, &snarl_manager, 50);
        TargetValueSearch tvs(graph, new TipAnchoredMaxDistance(di), 
                               new SnarlMinDistance(di)); 
    
        
        SECTION( "Test tvs" ) {
            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos1r = make_pos_t(1, true, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos6r = make_pos_t(6, true, 0);
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);
            pos_t pos9 = make_pos_t(9, false, 0);
            pos_t pos10 = make_pos_t(10, false, 0);
            pos_t pos11 = make_pos_t(11, false, 0);
            pos_t pos12 = make_pos_t(12, false, 0);

            REQUIRE(tvs.tv_path(pos1, pos2, 3, 5).size() == 2);
            REQUIRE(tvs.tv_path(pos1, pos11, 18, 5).size() == 10);
//TODO: This depends on the max distance            REQUIRE(tvs.tv_path(pos1, pos11, 26, 5).size() == 10);
            REQUIRE(tvs.tv_path(pos1, pos11, 8, 5).size() == 6);
            REQUIRE(tvs.tv_path(pos1, pos11, 7, 5).size() == 6);
            REQUIRE(tvs.tv_path(pos2, pos11, 8, 5).size() == 8);
            REQUIRE(tvs.tv_path(pos11, pos6r, 8, 5).size() == 6);
            REQUIRE(tvs.tv_path(pos11, pos6r, 9, 1).size() == 6);
            REQUIRE(tvs.tv_path(pos11, pos6r, 10, 0).size() == 6);
            REQUIRE(tvs.tv_path(pos11, pos6r, 9, 0).size() == 0);
            REQUIRE(tvs.tv_path(pos1, pos12, 3, 3).size() == 0);
            REQUIRE(tvs.tv_path(pos1, pos9, 9, 5).size() == 7);
            REQUIRE(tvs.tv_path(pos1, pos9, 60, 5).size() == 0);
           
        }
    }//End test case
    TEST_CASE( "TVS for unary snarl",
                   "[tvs][bug]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n7);
        Edge* e4 = graph.create_edge(n3, n5);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n1, n1, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        DistanceIndex di (&graph, &snarl_manager, 50);
        TargetValueSearch tvs(graph, new TipAnchoredMaxDistance(di), 
                               new SnarlMinDistance(di)); 
    
        
        SECTION( "Test tvs" ) {
            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos1r = make_pos_t(1, true, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos6r = make_pos_t(6, true, 0);
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);
            pos_t pos9 = make_pos_t(9, false, 0);
            pos_t pos10 = make_pos_t(10, false, 0);
            pos_t pos11 = make_pos_t(11, false, 0);
            pos_t pos12 = make_pos_t(12, false, 0);

            REQUIRE(tvs.tv_path(pos1, pos2, 3, 5).size() == 2);
            REQUIRE(tvs.tv_path(pos1, pos2, 5, 5).size() == 2);
            REQUIRE(tvs.tv_path(pos1r, pos2, 6, 5).size() == 3);
            REQUIRE(tvs.tv_path(pos6r, pos3, 11, 5).size() == 6);
        }
    }//End test case

}

}
