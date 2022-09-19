#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include "vg/io/json2pb.h"
#include "../vg.hpp"
#include "../xg.hpp"
#include "catch.hpp"
#include "../snarls.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../cluster.hpp"
#include "../snarl_distance_index.hpp"
#include "../genotypekit.hpp"
#include "random_graph.hpp"
#include <fstream>
#include <random>
#include <time.h> 

//#define print

namespace vg {
namespace unittest {


    TEST_CASE( "TVS for simple nested snarl",
                   "[tvs]" ) {
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

        SnarlDistanceIndex di;
        IntegratedSnarlFinder snarl_finder(graph);
        fill_in_distance_index (&di, &graph, &snarl_finder);
        TargetValueSearch tvs(graph, new TipAnchoredMaxDistance(di), 
                               new SnarlMinDistance(di)); 
    
        
        SECTION( "Test tvs" ) {
            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos1_2 = make_pos_t(1, false, 2);
            pos_t pos1r = make_pos_t(1, true, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos6r = make_pos_t(6, true, 0);
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);

            REQUIRE(tvs.tv_path(pos1, pos1_2, 2, 5).size() == 1);
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
                   "[tvs]" ) {
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

        SnarlDistanceIndex di;
        IntegratedSnarlFinder snarl_finder(graph);
        fill_in_distance_index (&di, &graph, &snarl_finder);
        TargetValueSearch tvs(graph, new TipAnchoredMaxDistance(di), 
                               new SnarlMinDistance(di)); 
    
        
        SECTION( "Test tvs" ) {

            REQUIRE(tvs.tv_path(make_pos_t(1, false, 0), 
                                make_pos_t(2, false, 0), 3, 5).size() == 2);
            REQUIRE(tvs.tv_path(make_pos_t(1, false, 0), 
                                make_pos_t(11, false, 0), 18, 5).size() == 10);
            REQUIRE(tvs.tv_path(make_pos_t(1, false, 0), 
                                make_pos_t(11, false, 0), 26, 5).size() != 0);
            REQUIRE(tvs.tv_path(make_pos_t(1, false, 0), 
                                make_pos_t(11, false, 0), 8, 5).size() == 6);
            REQUIRE(tvs.tv_path(make_pos_t(1, false, 0), 
                                make_pos_t(11, false, 0), 7, 5).size() == 6);
            REQUIRE(tvs.tv_path(make_pos_t(2, false, 0), 
                                make_pos_t(11, false, 0), 8, 5).size() == 8);
            REQUIRE(tvs.tv_path(make_pos_t(11, false, 0), 
                                make_pos_t(6, true, 0), 8, 5).size() == 6);
            REQUIRE(tvs.tv_path(make_pos_t(11, false, 0), 
                                make_pos_t(6, true, 0), 9, 1).size() == 6);
            REQUIRE(tvs.tv_path(make_pos_t(11, false, 0), 
                                make_pos_t(6, true, 0), 10, 0).size() == 6);
            REQUIRE(tvs.tv_path(make_pos_t(11, false, 0), 
                                make_pos_t(6, true, 0), 9, 0).size() == 0);
            REQUIRE(tvs.tv_path(make_pos_t(1, false, 0), 
                                make_pos_t(12, false, 0), 3, 3).size() == 0);
            REQUIRE(tvs.tv_path(make_pos_t(1, false, 0), 
                                make_pos_t(9, false, 0), 9, 5).size() == 7);
        }
    }//End test case
    TEST_CASE( "TVS for unary snarl",
                   "[tvs]" ) {
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

        SnarlDistanceIndex di;
        IntegratedSnarlFinder snarl_finder(graph);
        fill_in_distance_index (&di, &graph, &snarl_finder);
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

    TEST_CASE( "TVS for random graph",
                   "[tvs]" ) {
        for (int i = 0; i  < 0; i++) {
            VG graph;
            random_graph(1000, 20, 100, &graph);
            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls();
         
            SnarlDistanceIndex di;
            IntegratedSnarlFinder snarl_finder(graph);
            fill_in_distance_index (&di, &graph, &snarl_finder);
            TargetValueSearch tvs(graph, new TipAnchoredMaxDistance(di), 
                               new SnarlMinDistance(di));    

            vector<const Snarl*> allSnarls;
            auto addSnarl = [&] (const Snarl* s) {
                allSnarls.push_back(s);
            };
            snarl_manager.for_each_snarl_preorder(addSnarl);

            uniform_int_distribution<int> randSnarlIndex(0, allSnarls.size()-1);
            default_random_engine generator(time(NULL));
            for (int j = 0; j < 100; j++) {
                const Snarl* snarl1 = allSnarls[randSnarlIndex(generator)];
                const Snarl* snarl2 = allSnarls[randSnarlIndex(generator)];
                 
                pair<unordered_set<id_t>, unordered_set<edge_t>> contents1 = 
                           snarl_manager.shallow_contents(snarl1, graph, true);
                pair<unordered_set<id_t>, unordered_set<edge_t>> contents2 = 
                           snarl_manager.shallow_contents(snarl2, graph, true);
 
                vector<id_t> nodes1 (contents1.first.begin(), contents1.first.end());
                vector<id_t> nodes2 (contents2.first.begin(), contents2.first.end());

                
                uniform_int_distribution<int> randNodeIndex2(0,nodes2.size()-1);
                uniform_int_distribution<int> randNodeIndex1(0,nodes1.size()-1);

                id_t nodeID1 = nodes1[randNodeIndex1(generator)];
                id_t nodeID2 = nodes2[randNodeIndex2(generator)];
                handle_t node1 = graph.get_handle(nodeID1);
                handle_t node2 = graph.get_handle(nodeID2);
 
                off_t offset1 = uniform_int_distribution<int>(0,graph.get_length(node1) - 1)(generator);
                off_t offset2 = uniform_int_distribution<int>(0,graph.get_length(node2) - 1)(generator);

                pos_t pos1 = make_pos_t(nodeID1, 
                  uniform_int_distribution<int>(0,1)(generator) == 0,offset1 );
                pos_t pos2 = make_pos_t(nodeID2, 
                  uniform_int_distribution<int>(0,1)(generator) == 0, offset2 );
 
                int64_t minDist = minimum_distance(di, pos1, pos2);
                int64_t maxDist = maximum_distance(di, pos1, pos2);
                if (minDist != -1 && maxDist != 20 && minDist <= maxDist) {

                    REQUIRE(tvs.tv_path(pos1, pos2, minDist-10, 11).size() !=0);
                    REQUIRE(tvs.tv_path(pos1, pos2, minDist+10, 11).size() !=0);
                    REQUIRE(tvs.tv_path(pos1, pos2, minDist, 1).size() !=0);
                }

            }
        }
    }//End test case
    
    TEST_CASE("Path oriented distance functions perform as expected", "[cluster][mapping]") {
        
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
        Node* n9 = vg.create_node("A");
        Node* n10 = vg.create_node("C");
        Node* n11 = vg.create_node("G");
        Node* n12 = vg.create_node("T");
        Node* n13 = vg.create_node("A");
        Node* n14 = vg.create_node("C");
        Node* n15 = vg.create_node("C");
        
        vg.create_edge(n0, n1);
        vg.create_edge(n2, n0, true, true);
        vg.create_edge(n1, n3);
        vg.create_edge(n2, n3);
        vg.create_edge(n3, n4, false, true);
        vg.create_edge(n4, n5, true, false);
        vg.create_edge(n5, n6);
        vg.create_edge(n8, n6, false, true);
        vg.create_edge(n6, n7, false, true);
        vg.create_edge(n7, n9, true, true);
        vg.create_edge(n9, n10, true, false);
        vg.create_edge(n10, n11, false, false);
        vg.create_edge(n12, n11, false, true);
        vg.create_edge(n13, n12, false, false);
        vg.create_edge(n14, n13, true, false);
        vg.create_edge(n15, n14, true, true);
        
        Graph graph = vg.graph;
        
        Path* path = graph.add_path();
        path->set_name("path");
        Mapping* mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(n0->id());
        mapping->set_rank(1);
        mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(n2->id());
        mapping->set_rank(2);
        mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(n3->id());
        mapping->set_rank(3);
        mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(n4->id());
        mapping->mutable_position()->set_is_reverse(true);
        mapping->set_rank(4);
        mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(n5->id());
        mapping->set_rank(5);
        mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(n6->id());
        mapping->set_rank(6);
        mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(n8->id());
        mapping->mutable_position()->set_is_reverse(true);
        mapping->set_rank(7);
        
        xg::XG xg_index;
        xg_index.from_path_handle_graph(VG(graph));
        PathOrientedDistanceMeasurer measurer(&xg_index);
        
        SECTION("Distance approxmation produces exactly correct path distances when positions are on path") {
            
            measurer.max_walk = 10;
            int64_t dist = measurer.oriented_distance(make_pos_t(n2->id(), false, 0),
                                                      make_pos_t(n5->id(), false, 0));
            REQUIRE(dist == 7);
        }
        
        SECTION("Distance approxmation produces negative distance when order is reversed when positions are on path") {
            
            measurer.max_walk = 10;
            int64_t dist = measurer.oriented_distance(make_pos_t(n5->id(), false, 0),
                                                      make_pos_t(n2->id(), false, 0));
            REQUIRE(dist == -7);
            
        }
        
        SECTION("Distance approxmation produces correctly signed distances on reverse strand when positions are on path") {
            
            measurer.max_walk = 10;
            int64_t dist = measurer.oriented_distance(make_pos_t(n5->id(), true, n5->sequence().size()),
                                                      make_pos_t(n2->id(), true, n2->sequence().size()));

            REQUIRE(dist == 7);
            dist = measurer.oriented_distance(make_pos_t(n2->id(), true, n2->sequence().size()),
                                              make_pos_t(n5->id(), true, n5->sequence().size()));
            REQUIRE(dist == -7);
        }
        
        SECTION("Distance approxmation produces expected distances when positions are not on path") {
            
            measurer.max_walk = 10;
            int64_t dist = measurer.oriented_distance(make_pos_t(n1->id(), false, 3),
                                                      make_pos_t(n7->id(), true, 3));
            REQUIRE(dist == 15);
            
            dist = measurer.oriented_distance(make_pos_t(n7->id(), true, 3),
                                              make_pos_t(n1->id(), false, 3));
            REQUIRE(dist == -15);
            
        }
        
        SECTION("Distance approxmation produces expected output when search to path must traverse every type of edge") {
            
            measurer.max_walk = 20;
            int64_t dist = measurer.oriented_distance(make_pos_t(n1->id(), false, 3),
                                                      make_pos_t(n15->id(), false, 1));
            REQUIRE(dist == 22);
            
        }
        
        SECTION("Distance approxmation produces expected output when positions are on opposite strands") {
            
            measurer.max_walk = 10;
            int64_t dist = measurer.oriented_distance(make_pos_t(n1->id(), true, 3),
                                                      make_pos_t(n7->id(), true, 3));
            REQUIRE(dist == std::numeric_limits<int64_t>::max());
        }
        
        SECTION("Distance approxmation produces expected output when cannot reach each other within the maximum distance") {
            
            measurer.max_walk = 0;
            int64_t dist = measurer.oriented_distance(make_pos_t(n1->id(), false, 3),
                                                      make_pos_t(n7->id(), true, 3));
            REQUIRE(dist == std::numeric_limits<int64_t>::max());
        }
    }
}

}
