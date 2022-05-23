#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include "vg/io/json2pb.h"
#include "../vg.hpp"
#include "catch.hpp"
#include "../snarls.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../position.hpp"
#include "../min_distance.hpp"
#include "../genotypekit.hpp"
#include "random_graph.hpp"
#include "randomness.hpp"
#include <fstream>
#include <random>
#include <time.h> 

//#define print

namespace vg {
namespace unittest {

static pair<unordered_set<Node*>, unordered_set<Edge*> > pb_contents(
    VG& graph, const pair<unordered_set<id_t>, unordered_set<edge_t> >& contents) {
    pair<unordered_set<Node*>, unordered_set<Edge*> > ret;
    for (id_t node_id : contents.first) {
        ret.first.insert(graph.get_node(node_id));
    }
    for (const edge_t& edge_handle : contents.second) {
        Edge* edge = graph.get_edge(NodeTraversal(graph.get_node(graph.get_id(edge_handle.first)),
                                                  graph.get_is_reverse(edge_handle.first)),
                                    NodeTraversal(graph.get_node(graph.get_id(edge_handle.second)),
                                                  graph.get_is_reverse(edge_handle.second)));
        ret.second.insert(edge);
    }
    return ret;
}

int64_t min_distance(VG* graph, pos_t pos1, pos_t pos2){
    //Distance using djikstras algorithm

    auto cmp = [] (pair<pair<id_t, bool> , int64_t> x, 
                   pair<pair<id_t, bool>, int64_t> y ) {
        return (x.second > y.second);
    };
 
    int64_t shortestDistance = -1;
    if (get_id(pos1) == get_id(pos2) && is_rev(pos1) == is_rev(pos2)) { //if positions are on the same node

        int64_t nodeSize = graph->get_node(get_id(pos1))->sequence().size();
        int64_t offset1 = get_offset(pos1);
        int64_t offset2 = get_offset(pos2);

        if (offset1 <= offset2) {
            shortestDistance = offset2-offset1+1; //+1 to be consistent
        }

    }


    priority_queue< pair<pair<id_t, bool> , int64_t>, 
                    vector<pair<pair<id_t, bool>, int64_t>>,
                          decltype(cmp)> reachable(cmp); 
    handle_t currHandle = graph->get_handle(get_id(pos1), is_rev(pos1));

    int64_t dist = graph->get_length(currHandle) - get_offset(pos1);

    auto addFirst = [&](const handle_t& h) -> bool {
        pair<id_t, bool> node = make_pair(graph->get_id(h), 
                                          graph->get_is_reverse(h));
        reachable.push(make_pair(node, dist));
        return true;
    };
  
    graph->follow_edges(currHandle, false, addFirst);
    unordered_set<pair<id_t, bool>> seen;
    seen.insert(make_pair(get_id(pos1), is_rev(pos1)));
    while (reachable.size() > 0) {
    
        pair<pair<id_t, bool>, int64_t> next = reachable.top();
        reachable.pop();
        pair<id_t, bool> currID = next.first;
        dist = next.second;
        if (seen.count(currID) == 0) {
 
            seen.insert(currID);
            currHandle = graph->get_handle(currID.first, currID.second);
            int64_t currDist = graph->get_length(currHandle); 
        
            auto addNext = [&](const handle_t& h) -> bool {
                pair<id_t, bool> node = make_pair(graph->get_id(h), 
                                                graph->get_is_reverse(h));
                reachable.push(make_pair(node, currDist + dist));
                return true;
            };
            graph->follow_edges(currHandle, false, addNext);
        
        }

        if (currID.first == get_id(pos2) && currID.second == is_rev(pos2)){
        //Dist is distance to beginning or end of node containing pos2
        
            if (is_rev(pos2) == currID.second) { 
                dist = dist + get_offset(pos2) + 1;
            } else {
                dist = dist +  graph->get_node(get_id(pos2))->sequence().size() - 
                            get_offset(pos2);
            }
            if (shortestDistance == -1) {shortestDistance = dist;}
            else {shortestDistance = min(dist, shortestDistance);}
        }

    }
    return shortestDistance == -1 ? -1 : shortestDistance-1;
};

    TEST_CASE( "One node",
                   "[min_dist]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCATGGAGCGTTGAGTGCGGGTG");

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 


        
        SECTION( "Create distance index" ) {

            MinimumDistanceIndex di (&graph, &snarl_manager);
#ifdef print
            di.print_self();
#endif

            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(1, false, 3) ) == 3);

        }

    }
    TEST_CASE( "Create min distance index for simple nested snarl",
                   "[min_dist]" ) {
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

        
        SECTION( "Create distance index" ) {

            MinimumDistanceIndex di (&graph, &snarl_manager);
#ifdef print
            di.print_self();
#endif

            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(2, false, 0) ) == 3);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0) ) == 3);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0) ) == 4);
            REQUIRE(di.min_distance(make_pos_t(2, false, 0),
                                   make_pos_t(6, false, 0) ) == 1);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0) ) == 5);

        } 
        SECTION( "Create max distance index" ) {

            MinimumDistanceIndex di (&graph, &snarl_manager, 20);
#ifdef print
            di.print_self();
#endif

            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(2, false, 0) ) == 3);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0) ) == 3);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0) ) == 4);
            REQUIRE(di.min_distance(make_pos_t(2, false, 0),
                                   make_pos_t(6, false, 0) ) == 1);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0) ) == 5);

            REQUIRE(di.max_distance(make_pos_t(1, false, 0),
                                   make_pos_t(2, false, 0) ) >= 3);
            REQUIRE(di.max_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0) ) >= 13);
            REQUIRE(di.max_distance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0) ) >= 5);
            REQUIRE(di.max_distance(make_pos_t(2, false, 0),
                                   make_pos_t(6, false, 0) ) >= 1);
            REQUIRE(di.max_distance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0) ) >= 12);

        }

    }//End test case

    TEST_CASE( "Create min distance index disconnected graph",
                   "[min_dist]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
//Disconnected
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("CTGA");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);

        Edge* e11 = graph.create_edge(n9, n10);
        Edge* e12 = graph.create_edge(n9, n11);
        Edge* e13 = graph.create_edge(n10, n11);
        Edge* e14 = graph.create_edge(n11, n12);
        Edge* e15 = graph.create_edge(n11, n13);
        Edge* e16 = graph.create_edge(n12, n13);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        MinimumDistanceIndex di (&graph, &snarl_manager);
        
        SECTION( "Create distance index" ) {
            //REQUIRE(di.chainDistances.size() == 0);
#ifdef print
            di.print_self();
#endif
        }


        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);
            const Snarl* snarl9 = snarl_manager.into_which_snarl(9, false);
            const Snarl* snarl11 = snarl_manager.into_which_snarl(11, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos9 = make_pos_t(9, false, 0);
            pos_t pos11 = make_pos_t(11, false, 0);
            pos_t pos12 = make_pos_t(12, false, 0);

            REQUIRE(di.min_distance( pos1, pos4) == 4);
            REQUIRE(di.min_distance( pos1, pos6) == 7);
            REQUIRE(di.min_distance( pos6, pos7) == -1);

            REQUIRE(di.min_distance( pos9, pos12) == 5);
            REQUIRE(di.min_distance(pos9, pos11) == 1);
            REQUIRE(di.min_distance( pos4, pos9) == -1);
            REQUIRE(di.min_distance( pos4, pos12) == -1);
            REQUIRE(di.min_distance( pos1, pos11) == -1);


        }
    }//End test case


    TEST_CASE("Simple chain min distance", "[min_dist]") {
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
        Edge* e4 = graph.create_edge(n5, n6);
        Edge* e5 = graph.create_edge(n2, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        SECTION("Min distance") {

            MinimumDistanceIndex di (&graph, &snarl_manager);
#ifdef print
            di.print_self();
#endif

            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) == 5);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0)) == 8);
            REQUIRE(di.min_distance(make_pos_t(6, true, 0),
                                   make_pos_t(1, true, 0)) == 6);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(6, false, 0)) == 4);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) == 4);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(8, false, 0)) == 5);

        }

        SECTION("Max distance") {

            MinimumDistanceIndex di (&graph, &snarl_manager, 20);
#ifdef print
            di.print_self();
#endif

            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) == 5);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0)) == 8);
            REQUIRE(di.min_distance(make_pos_t(6, true, 0),
                                   make_pos_t(1, true, 0)) == 6);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(6, false, 0)) == 4);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) == 4);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(8, false, 0)) == 5);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) == 3);

            REQUIRE(di.max_distance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) >= 8);
            REQUIRE(di.max_distance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0)) >= 11);
            REQUIRE(di.max_distance(make_pos_t(6, true, 0),
                                   make_pos_t(1, true, 0)) >= 9);
            REQUIRE(di.max_distance(make_pos_t(3, false, 0),
                                   make_pos_t(6, false, 0)) >= 4);
            REQUIRE(di.max_distance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) >= 5);
            REQUIRE(di.max_distance(make_pos_t(3, false, 0),
                                   make_pos_t(8, false, 0)) >= 6);
            REQUIRE(di.max_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) >= 13);

        }
        
    }//end test case

    TEST_CASE("Chain with reversing edge min distance", "[min_dist]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("TTTTTTTTTTTT"); //12 T's
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n8);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n5, n6);
        Edge* e5 = graph.create_edge(n2, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n5, n5, false, true);
        Edge* e12 = graph.create_edge(n6, n6, false, true);
        Edge* e13 = graph.create_edge(n7, n7, false, true);


        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        SECTION("Min distance") {

            MinimumDistanceIndex di (&graph, &snarl_manager);
            #ifdef print
            di.print_self();
            #endif
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) == 3);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) == 4);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(4, true, 0)) == 7);
            REQUIRE(di.min_distance(make_pos_t(5, false, 0),
                                   make_pos_t(6, true, 0)) == 5);



        }
        SECTION("Max distance") {
            MinimumDistanceIndex di (&graph, &snarl_manager, 30);
            #ifdef print
            di.print_self();
            #endif
            REQUIRE(di.max_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) >= 13);
            REQUIRE(di.max_distance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) >= 5);
            REQUIRE(di.max_distance(make_pos_t(3, false, 0),
                                   make_pos_t(4, true, 0)) >= 27);
            REQUIRE(di.max_distance(make_pos_t(5, false, 0),
                                   make_pos_t(6, true, 0)) >= 15);



        }
    }//end test case

    TEST_CASE("Top level chain min", "[min_dist]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n8);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n5, n6);
        Edge* e5 = graph.create_edge(n2, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n8, n9);
        Edge* e12 = graph.create_edge(n9, n10);
        Edge* e13 = graph.create_edge(n8, n10);
        Edge* e14 = graph.create_edge(n5, n5, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 


        SECTION("Create distance index") {
            MinimumDistanceIndex di (&graph, &snarl_manager);
            #ifdef print
            di.print_self();
            #endif
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0)) == 8);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(10, false, 0)) == 7);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) == 5);
            REQUIRE(di.min_distance(make_pos_t(2, false, 0),
                                   make_pos_t(9, false, 0)) == 10);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(9, false, 0)) == 9);

        }
    }//end test case
    TEST_CASE("Interior chain min", "[min_dist]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n4, n7);
        Edge* e9 = graph.create_edge(n5, n7);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n7, n8, false, true);
        Edge* e13 = graph.create_edge(n8, n9);
        Edge* e14 = graph.create_edge(n8, n9, true, false);
        Edge* e15 = graph.create_edge(n9, n10);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();


        SECTION("Create distance index") { 
            MinimumDistanceIndex di (&graph, &snarl_manager);
            #ifdef print
            di.print_self();
            #endif
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0)) == 8);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) == 9);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(8, true, 0)) == 9);
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) == 8);
            REQUIRE(di.min_distance(make_pos_t(8, false, 0),
                                   make_pos_t(5, true, 0)) == 5);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(6, false, 0)) == 5);
            REQUIRE(di.min_distance(make_pos_t(3, false, 0),
                                   make_pos_t(10, false, 0)) == 11);
            REQUIRE(di.min_distance(make_pos_t(7, false, 0),
                                   make_pos_t(1, true, 0)) == 11);

        }
    }//end test case
    TEST_CASE("Top level loop creates unary snarl min", "[min_dist]") {
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
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n1, n1, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 


        // We end up with a big unary snarl of 7 rev -> 7 rev
        // Inside that we have a chain of two normal snarls 2 rev -> 3 fwd, and 3 fwd -> 6 fwd
        // And inside 2 rev -> 3 fwd, we get 1 rev -> 1 rev as another unary snarl.
        
        // We name the snarls for the distance index by their start nodes.

        SECTION("Minimum distance") {
            MinimumDistanceIndex di (&graph, &snarl_manager);
            #ifdef print
            di.print_self();
            #endif
            REQUIRE(di.min_distance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0)) == 4);
            REQUIRE(di.min_distance(make_pos_t(2, true, 0),
                                   make_pos_t(3, false, 0)) == 7);
            REQUIRE(di.min_distance(make_pos_t(4, true, 0),
                                   make_pos_t(2, false, 0)) == 11);
        }
    }//end test case
    TEST_CASE( "Shortest path exits common ancestor min","[min_dist]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("A");
        Node* n10 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n9);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n6, n8);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n2, n2, true, false);
        Edge* e15 = graph.create_edge(n5, n5);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        MinimumDistanceIndex di (&graph, &snarl_manager);
        
        SECTION( "Min distance" ) {
            REQUIRE(di.min_distance(make_pos_t(5, false, 1), make_pos_t(5, false, 0)) == 11); 
            REQUIRE(di.min_distance(make_pos_t(5, false, 0), make_pos_t(5, false, 0)) == 0); 
        }

    }//End test case

    TEST_CASE( "More loops min",
                   "[min_dist]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n7);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n3, n3, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        MinimumDistanceIndex di (&graph, &snarl_manager);
        


        SECTION ("Min distance") {

            REQUIRE(di.min_distance(make_pos_t(1, false, 0), make_pos_t(4, false, 0)) == 4);
            REQUIRE(di.min_distance( make_pos_t(5, false, 0), make_pos_t(6, false, 0)) == -1);
            REQUIRE(di.min_distance( make_pos_t(5, false, 0), make_pos_t(6, true, 0)) == -1);
            REQUIRE(di.min_distance( make_pos_t(5, true, 2), make_pos_t(6, false, 0)) == 11);
            REQUIRE(di.min_distance(make_pos_t(2, false, 0), make_pos_t(7, false, 0)) == 6);

            REQUIRE(min_distance(&graph, make_pos_t(1, false, 0), make_pos_t(4, false, 0)) == 4);
            REQUIRE(min_distance(&graph, make_pos_t(5, true, 2), make_pos_t(6, false, 0)) == 11);
            REQUIRE(min_distance(&graph, make_pos_t(2, false, 0), make_pos_t(7, false, 0)) == 6);

      
        }
    }//End test case

    TEST_CASE( "Exit common ancestor min","[min_dist]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("AA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("G");
        Node* n12 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n11);
        Edge* e5 = graph.create_edge(n11, n9);
        Edge* e6 = graph.create_edge(n3, n4);
        Edge* e7 = graph.create_edge(n3, n5);
        Edge* e8 = graph.create_edge(n4, n6);
        Edge* e9 = graph.create_edge(n5, n6);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n6, n12);
        Edge* e12 = graph.create_edge(n12, n8);
        Edge* e13 = graph.create_edge(n7, n8);
        Edge* e14 = graph.create_edge(n8, n9);
        Edge* e15 = graph.create_edge(n9, n10);
        Edge* e16 = graph.create_edge(n2, n2, true, false);
        Edge* e17 = graph.create_edge(n9, n9, false, true);
        Edge* e18 = graph.create_edge(n2, n9, true, true);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        


        SECTION ("Min distance") {
            MinimumDistanceIndex di (&graph, &snarl_manager);
#ifdef print
            di.print_self();
#endif

            REQUIRE(di.min_distance( make_pos_t(2, true, 0), 
                                    make_pos_t(9, true, 1)) == 2);
            REQUIRE(di.min_distance( make_pos_t(2, false, 0), 
                                    make_pos_t(9, true, 1)) == 5);
            REQUIRE(di.min_distance( make_pos_t(3, true, 0),
                                    make_pos_t(9, false, 0)) == 4);
            REQUIRE(di.min_distance( make_pos_t(3, true, 0), 
                                    make_pos_t(9, true, 1)) == 3);
            REQUIRE(di.min_distance( make_pos_t(3, false, 0), 
                                    make_pos_t(9, false, 0)) == 11);
            REQUIRE(di.min_distance( make_pos_t(4, false, 0), 
                                    make_pos_t(5, false, 0)) == 14);
            REQUIRE(di.min_distance( make_pos_t(4, true, 0), 
                                    make_pos_t(5, false, 0)) == 8);
            REQUIRE(di.min_distance(make_pos_t(7, false, 0), 
                                   make_pos_t(12, false, 0)) == 14);
            REQUIRE(di.min_distance(make_pos_t(7, false, 0), 
                                   make_pos_t(12, true, 0)) == 13);
            REQUIRE(di.min_distance(make_pos_t(7, true, 0),
                                   make_pos_t(12, false, 0)) == 15);
            REQUIRE(di.min_distance(make_pos_t(8, false, 0), 
                                   make_pos_t(11, false, 0)) == 7);

        }
    }//End test case

    TEST_CASE( "Exit common ancestor chain min","[min_dist]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");
        Node* n9 = graph.create_node("AA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("G");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("GA");
        Node* n14 = graph.create_node("G");
        Node* n15 = graph.create_node("G");
        Node* n16 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n13);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n16);
        Edge* e27 = graph.create_edge(n16, n9);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n6, n8);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n9, n11);
        Edge* e15 = graph.create_edge(n10, n11);
        Edge* e16 = graph.create_edge(n11, n12);
        Edge* e17 = graph.create_edge(n11, n2);
        Edge* e18 = graph.create_edge(n12, n1);
        Edge* e19 = graph.create_edge(n13, n14);
        Edge* e20 = graph.create_edge(n13, n15);
        Edge* e21 = graph.create_edge(n14, n15);
        Edge* e22 = graph.create_edge(n15, n12);
        Edge* e23 = graph.create_edge(n2, n2, true, false);
        Edge* e24 = graph.create_edge(n11, n11, false, true);
        Edge* e25 = graph.create_edge(n1, n1, true, false);
        Edge* e26 = graph.create_edge(n12, n12, false, true);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        


        SECTION ("Min distance") {
        MinimumDistanceIndex di (&graph, &snarl_manager);

#ifdef print
            di.print_self();
#endif

            REQUIRE(di.min_distance(make_pos_t(2, true, 0), 
                                   make_pos_t(10, true, 0)) == 2);
            REQUIRE(di.min_distance(make_pos_t(2, false, 0), 
                                   make_pos_t(10, false, 0)) == 4);
            REQUIRE(di.min_distance(make_pos_t(4, true, 3), 
                                   make_pos_t(5, false, 0)) == 5);
            REQUIRE(di.min_distance(make_pos_t(4, false, 0),
                                   make_pos_t(5, false, 0)) == 11);
            REQUIRE(di.min_distance(make_pos_t(4, true, 3),
                                   make_pos_t(5, true, 0)) == 8);
            REQUIRE(di.min_distance(make_pos_t(14, false, 0),
                                   make_pos_t(10, false, 0)) == 10);
            REQUIRE(di.min_distance(make_pos_t(14, false, 0), 
                                   make_pos_t(10, true, 0)) == 5);
            REQUIRE(di.min_distance(make_pos_t(14, false, 0), 
                                   make_pos_t(3, false, 0)) == 7);
            REQUIRE(di.min_distance(make_pos_t(16, false, 0),
                                   make_pos_t(3, false, 0)) == 5);
        }
    }//End test case

    TEST_CASE("Top level loops min", "[min_dist]") {
        VG graph;

        Node* n1 = graph.create_node("G");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGAAAAAAAAAAAA"); //15
        Node* n5 = graph.create_node("GCAA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("A");
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n9, n1);
        Edge* e2 = graph.create_edge(n9, n10);
        Edge* e3 = graph.create_edge(n1, n2);
        Edge* e4 = graph.create_edge(n1, n8);
        Edge* e5 = graph.create_edge(n2, n3);
        Edge* e6 = graph.create_edge(n2, n4);
        Edge* e7 = graph.create_edge(n3, n5);
        Edge* e8 = graph.create_edge(n4, n5);
        Edge* e9 = graph.create_edge(n5, n6);
        Edge* e10 = graph.create_edge(n5, n7);
        Edge* e11 = graph.create_edge(n6, n7);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n8, n10);
        Edge* e14 = graph.create_edge(n10, n10, false, true);
        Edge* e15 = graph.create_edge(n9, n9, true, false);
        Edge* e16 = graph.create_edge(n10, n9);
        Edge* e17 = graph.create_edge(n2, n2, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        SECTION ("Min distance") {

            MinimumDistanceIndex di (&graph, &snarl_manager);
            #ifdef print
                di.print_self();
            #endif

            REQUIRE(di.min_distance(make_pos_t(4, false, 0), 
                                   make_pos_t(6, false, 0)) == 19);
            REQUIRE(di.min_distance(make_pos_t(4, false, 0), 
                                   make_pos_t(6, true, 0)) == 25);
            REQUIRE(di.min_distance(make_pos_t(2, true, 0), 
                                   make_pos_t(7, false, 0)) == 7);
            REQUIRE(di.min_distance(make_pos_t(2, false, 0), 
                                   make_pos_t(7, false, 0)) == 6);
            REQUIRE(di.min_distance(make_pos_t(2, false, 0), 
                                   make_pos_t(7, true, 0)) == 11);

            REQUIRE(min_distance(&graph, make_pos_t(4, false, 0),  
                                         make_pos_t(6, true, 0)) == 25);
            REQUIRE(min_distance(&graph, make_pos_t(2, false, 0),
                                         make_pos_t(7, true, 0)) == 11);
        }
        SECTION ("Max distance") {

            MinimumDistanceIndex di (&graph, &snarl_manager, 50);
            #ifdef print
                di.print_self();
            #endif

            REQUIRE(di.max_distance(make_pos_t(4, false, 0), 
                                   make_pos_t(6, false, 0)) >= 40);
            REQUIRE(di.max_distance(make_pos_t(4, false, 0), 
                                   make_pos_t(6, true, 0)) >= 25);
            REQUIRE(di.max_distance(make_pos_t(2, false, 0), 
                                   make_pos_t(7, false, 0)) >= 22);

        }
    }//end test case


    TEST_CASE( "Snarl that loops on itself min",
                  "[min_dist]" ) {
        
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
            
            // Define the snarls for the top level
           
            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls();


       SECTION ("Distance functions") {
            MinimumDistanceIndex di (&graph, &snarl_manager);
            #ifdef print
                di.print_self();
            #endif
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, true);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);

            REQUIRE(di.min_distance(pos1, pos2) == min_distance(&graph, pos1, pos2));



        }
 
    }

    TEST_CASE( "Get connected component and length of root chains",
                  "[min_dist]" ) {
        
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
//Disconnected
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("CTGA");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("CTGA");
//Disconnected
        Node* n14 = graph.create_node("T");
        Node* n15 = graph.create_node("G");
        Node* n16 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);

        Edge* e11 = graph.create_edge(n9, n10);
        Edge* e12 = graph.create_edge(n9, n11);
        Edge* e13 = graph.create_edge(n10, n11);
        Edge* e14 = graph.create_edge(n11, n12);
        Edge* e15 = graph.create_edge(n11, n13);
        Edge* e16 = graph.create_edge(n12, n13);
 
        Edge* e17 = graph.create_edge(n14, n15);
        Edge* e18 = graph.create_edge(n14, n16);
        Edge* e19 = graph.create_edge(n15, n16);           
        // Define the snarls for the top level
       
        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();


       SECTION ("Simple disconnected graph") {
            MinimumDistanceIndex di (&graph, &snarl_manager);
            #ifdef print
                di.print_self();
            #endif

            //Connected components should be have unique identifiers

            REQUIRE (std::get<0>(di.get_minimizer_distances(make_pos_t(1 , false, 0))) );
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(2 , false, 0))));
            REQUIRE (std::get<1>(di.get_minimizer_distances(make_pos_t(3 , false, 0))) == 
                     std::get<1>(di.get_minimizer_distances(make_pos_t(1, false, 0))) );
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(4 , false, 0))));
            REQUIRE (std::get<1>(di.get_minimizer_distances(make_pos_t(5 , false, 0))) == 
                     std::get<1>(di.get_minimizer_distances(make_pos_t(1, false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(6 , false, 0))));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(7 , false, 0))));
            REQUIRE (std::get<1>(di.get_minimizer_distances(make_pos_t(8 , false, 0))) == 
                     std::get<1>(di.get_minimizer_distances(make_pos_t(1, false, 0)) ));
            REQUIRE (std::get<0>(di.get_minimizer_distances(make_pos_t(9 , false, 0))));
            REQUIRE (std::get<1>(di.get_minimizer_distances(make_pos_t(9 , false, 0))) != 
                     std::get<1>(di.get_minimizer_distances(make_pos_t(1, false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(10, false, 0))));
            REQUIRE (std::get<1>(di.get_minimizer_distances(make_pos_t(11, false, 0))) == 
                     std::get<1>(di.get_minimizer_distances(make_pos_t(9, false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(12, false, 0))));
            REQUIRE (std::get<1>(di.get_minimizer_distances(make_pos_t(13, false, 0))) == 
                     std::get<1>(di.get_minimizer_distances(make_pos_t(9, false, 0)) ));

            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(14, false, 0))));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(15, false, 0))));
            //Offsets should be correct
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(1 , false, 0))) == 1);
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(3 , false, 0))) == 4);
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(5 , false, 0))) == 5);
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(8 , false, 0))) == 9);
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(9 , false, 0))) == 1);
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(11, false, 0))) == 2);
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(13, false, 0))) == 6);


            //Top-level bubble assignments should be correct
            REQUIRE (std::get<3>(di.get_minimizer_distances(make_pos_t(2 , false, 0))) );
            REQUIRE (std::get<3>(di.get_minimizer_distances(make_pos_t(4 , false, 0))) );
            REQUIRE (std::get<3>(di.get_minimizer_distances(make_pos_t(6 , false, 0))) );
            REQUIRE (std::get<3>(di.get_minimizer_distances(make_pos_t(7 , false, 0))) );

            //ranks of top-level bubbles
            REQUIRE (std::get<4>(di.get_minimizer_distances(make_pos_t(2 , false, 0)))== 0 );
            REQUIRE (std::get<4>(di.get_minimizer_distances(make_pos_t(4 , false, 0)))== 1 );
            REQUIRE (std::get<4>(di.get_minimizer_distances(make_pos_t(6 , false, 0)))== 2 );
            REQUIRE (std::get<4>(di.get_minimizer_distances(make_pos_t(7 , false, 0)))== 2 );

            //start node length
            REQUIRE ((std::get<5>(di.get_minimizer_distances(make_pos_t(2 , false, 0)))== 3 || std::get<5>(di.get_minimizer_distances(make_pos_t(2 , false, 0)))== 1));
            REQUIRE ((std::get<5>(di.get_minimizer_distances(make_pos_t(4 , false, 0)))== 1 ||std::get<5>(di.get_minimizer_distances(make_pos_t(4 , false, 0)))== 3) );
            REQUIRE ((std::get<5>(di.get_minimizer_distances(make_pos_t(6 , false, 0)))== 3 || std::get<5>(di.get_minimizer_distances(make_pos_t(6 , false, 0)))== 4));
            REQUIRE ((std::get<5>(di.get_minimizer_distances(make_pos_t(7 , false, 0)))== 3 || std::get<5>(di.get_minimizer_distances(make_pos_t(7 , false, 0)))== 4));
            //end node length
            REQUIRE ((std::get<6>(di.get_minimizer_distances(make_pos_t(2 , false, 0)))== 3 || std::get<6>(di.get_minimizer_distances(make_pos_t(2 , false, 0)))== 1));
            REQUIRE ((std::get<6>(di.get_minimizer_distances(make_pos_t(4 , false, 0)))== 1 ||std::get<6>(di.get_minimizer_distances(make_pos_t(4 , false, 0)))== 3) );
            REQUIRE ((std::get<6>(di.get_minimizer_distances(make_pos_t(6 , false, 0)))== 3 || std::get<6>(di.get_minimizer_distances(make_pos_t(6 , false, 0)))== 4));
            REQUIRE ((std::get<6>(di.get_minimizer_distances(make_pos_t(7 , false, 0)))== 3 || std::get<6>(di.get_minimizer_distances(make_pos_t(7 , false, 0)))== 4));

            //Node lengths
            REQUIRE (std::get<7>(di.get_minimizer_distances(make_pos_t(2 , false, 0)))== 1 );
            REQUIRE (std::get<7>(di.get_minimizer_distances(make_pos_t(4 , false, 0)))== 4 );
            REQUIRE (std::get<7>(di.get_minimizer_distances(make_pos_t(6 , false, 0)))== 1 );
            REQUIRE (std::get<7>(di.get_minimizer_distances(make_pos_t(7 , false, 0)))== 1 );
        }
 
    }
    TEST_CASE( "Get connected component and length of root chains in nested graph",
                  "[min_dist]" ) {
        
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");

//Disconnected
        Node* n11 = graph.create_node("CTGA");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("CTGA");
        Node* n14 = graph.create_node("T");
        Node* n15 = graph.create_node("G");
        Node* n16 = graph.create_node("G");
        Node* n17 = graph.create_node("CTGA");
        Node* n18 = graph.create_node("T");
        Node* n19 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n5);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n10);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n7, n9);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);

        Edge* e14 = graph.create_edge(n11, n12);
        Edge* e15 = graph.create_edge(n11, n13);
        Edge* e16 = graph.create_edge(n12, n13);
        Edge* e17 = graph.create_edge(n13, n14);
        Edge* e18 = graph.create_edge(n13, n19);
        Edge* e19 = graph.create_edge(n14, n15);
        Edge* e20 = graph.create_edge(n14, n16);
        Edge* e21 = graph.create_edge(n15, n16);
        Edge* e22 = graph.create_edge(n16, n17);
        Edge* e23 = graph.create_edge(n16, n18);
        Edge* e24 = graph.create_edge(n17, n18);
        Edge* e25 = graph.create_edge(n18, n19);
        // Define the snarls for the top level
       
        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();


       SECTION ("Nested disconnected graph") {
            MinimumDistanceIndex di (&graph, &snarl_manager);
            #ifdef print
                di.print_self();
            #endif

            //Connected components should be have unique identifiers

            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(1 , false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(2 , false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(3 , false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(4 , true, 3))  ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(5 , false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(6 , false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(7 , false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(8 , false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(9 , false, 0)) ));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(10 , false, 0))));
            REQUIRE (std::get<0>(di.get_minimizer_distances(make_pos_t(11 , false, 0))));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(12 , false, 0))));
            REQUIRE (std::get<1>(di.get_minimizer_distances(make_pos_t(13, false, 0)) ) == 
                     std::get<1>(di.get_minimizer_distances(make_pos_t(11 , false, 0))));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(14 , false, 0))));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(15 , false, 0))));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(16 , false, 0))));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(17 , false, 0))));
            REQUIRE (!std::get<0>(di.get_minimizer_distances(make_pos_t(18 , false, 0))));
            REQUIRE (std::get<1>(di.get_minimizer_distances(make_pos_t(19, false, 0)) ) == 
                     std::get<1>(di.get_minimizer_distances(make_pos_t(11 , false, 0))));

            //Offsets should be correct
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(11, true, 3)) ) == 1);
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(13, false, 0))) == 5);
            REQUIRE (std::get<2>(di.get_minimizer_distances(make_pos_t(19, false, 0))) == 9);

            //Check simple snarl assignments
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(1 , false, 0)) ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(2 , false, 0)) ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(3 , false, 0)) ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(4 , true, 3))  ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(5 , false, 0)) ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(6 , false, 0)) ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(7 , false, 0)) ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(8 , false, 0)) ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(9 , false, 0)) ));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(10 , false, 0))));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(11 , false, 0))));
            REQUIRE (std::get<3>(di.get_minimizer_distances(make_pos_t(12 , false, 0))));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(13 , false, 0))));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(14 , false, 0))));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(15 , false, 0))));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(16 , false, 0))));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(17 , false, 0))));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(18 , false, 0))));
            REQUIRE (!std::get<3>(di.get_minimizer_distances(make_pos_t(19 , false, 0))));

            REQUIRE (std::get<4>(di.get_minimizer_distances(make_pos_t(12 , false, 0))) == 0);
            REQUIRE (std::get<5>(di.get_minimizer_distances(make_pos_t(12 , false, 0))) == 4);
            REQUIRE (std::get<6>(di.get_minimizer_distances(make_pos_t(12 , false, 0))) == 4);
            REQUIRE (std::get<7>(di.get_minimizer_distances(make_pos_t(12 , false, 0))) == 1);


        }
 
    }

    TEST_CASE("Random test min", "[min_dist][rand]") {

/*
        ifstream vg_stream("testGraph");
        VG vg(vg_stream);
        vg_stream.close();
        CactusSnarlFinder bubble_finder(vg);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        MinimumDistanceIndex di (&vg, &snarl_manager);
        di.print_self();
        pos_t pos1 = make_pos_t(208, true, 4);
        pos_t pos2 = make_pos_t(256, false, 2); 
        REQUIRE(di.min_distance(pos1, pos2) == min_distance(&vg, pos1, pos2));
*/

        for (int i = 0; i < 0; i++) {
            //1000 different graphs
            VG graph;
            random_graph(1000, 20, 100, &graph);

            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls(); 

            uint64_t cap = 50;
            MinimumDistanceIndex di (&graph, &snarl_manager, cap);
            #ifdef print        
                di.print_self();

            #endif



            vector<const Snarl*> allSnarls;
            auto addSnarl = [&] (const Snarl* s) {
                allSnarls.push_back(s);
            };
            snarl_manager.for_each_snarl_preorder(addSnarl);

            uniform_int_distribution<int> randSnarlIndex(0, allSnarls.size()-1);
            default_random_engine generator(time(NULL));
            for (int j = 0; j < 100; j++) {
                //Check distances for random pairs of positions 
                const Snarl* snarl1 = allSnarls[randSnarlIndex(generator)];
                const Snarl* snarl2 = allSnarls[randSnarlIndex(generator)];

                id_t start1_id = snarl1->start().node_id();
                bool start1_rev = snarl1->start().backward();
                id_t end1_id = snarl1->end().node_id();
                bool end1_rev = !snarl1->end().backward();

                REQUIRE((di.into_which_snarl(start1_id, start1_rev) == make_tuple(start1_id, start1_rev, snarl_manager.is_trivial(snarl1, graph)) ||
                         di.into_which_snarl(start1_id, start1_rev) == make_tuple(end1_id, end1_rev,     snarl_manager.is_trivial(snarl1, graph))));

                 
                pair<unordered_set<Node*>, unordered_set<Edge*>> contents1 = 
                    pb_contents(graph, snarl_manager.shallow_contents(snarl1, graph, true));
                pair<unordered_set<Node*>, unordered_set<Edge*>> contents2 = 
                    pb_contents(graph, snarl_manager.shallow_contents(snarl2, graph, true));
 
                vector<Node*> nodes1 (contents1.first.begin(), contents1.first.end());
                vector<Node*> nodes2 (contents2.first.begin(), contents2.first.end());

                
                uniform_int_distribution<int> randNodeIndex2(0,nodes2.size()-1);
                uniform_int_distribution<int> randNodeIndex1(0,nodes1.size()-1);
 
                Node* node1 = nodes1[randNodeIndex1(generator)];
                Node* node2 = nodes2[randNodeIndex2(generator)];
                id_t nodeID1 = node1->id();
                id_t nodeID2 = node2->id();
 
                off_t offset1 = uniform_int_distribution<int>(0,node1->sequence().size() - 1)(generator);
                off_t offset2 = uniform_int_distribution<int>(0,node2->sequence().size() - 1)(generator);

                pos_t pos1 = make_pos_t(nodeID1, uniform_int_distribution<int>(0,1)(generator) == 0, offset1 );
                pos_t pos2 = make_pos_t(nodeID2, uniform_int_distribution<int>(0,1)(generator) == 0, offset2 );
                //auto chain_info = di.get_minimizer_distances(pos1);
                //auto encoded_chain_info = MIPayload::decode(MIPayload::encode(di.get_minimizer_distances(pos1)));
                //assert (chain_info == encoded_chain_info);

                const Snarl* pos1_snarl = snarl_manager.into_which_snarl(nodeID1, false);

                if (pos1_snarl == nullptr) {
                    REQUIRE(di.into_which_snarl(nodeID1, false) == make_tuple((id_t)0, false, false));
                } else {
                    auto result = di.into_which_snarl(nodeID1, false);
                    //cerr << "node: " << nodeID1 << ": guess: " << std::get<0>(result) << " " << std::get<1>(result) << " " << std::get<2>(result) 
                    //     << " actual snarl start: "<< pos1_snarl->start().node_id() << " " <<  pos1_snarl->start().backward() 
                    //     << ", actual snarl end: " << pos1_snarl->end().node_id() << " " <<  pos1_snarl->end().backward() 
                    //     << " trivial? " <<  snarl_manager.is_trivial(pos1_snarl, graph) << endl; 
                    REQUIRE((di.into_which_snarl(nodeID1, false) == make_tuple((id_t)pos1_snarl->start().node_id(), pos1_snarl->start().backward(), snarl_manager.is_trivial(pos1_snarl, graph)) ||
                             di.into_which_snarl(nodeID1, false) == make_tuple((id_t)pos1_snarl->end().node_id(), !pos1_snarl->end().backward(), snarl_manager.is_trivial(pos1_snarl, graph))));
                }

 
                if (!(nodeID1 != snarl1->start().node_id() && 
                    (snarl_manager.into_which_snarl(nodeID1, false) != NULL ||
                      snarl_manager.into_which_snarl(nodeID1, true) != NULL)) &&
                    ! (nodeID2 != snarl2->start().node_id() &&
                        (snarl_manager.into_which_snarl(nodeID2, false) != NULL
                       || snarl_manager.into_which_snarl(nodeID2, true) != NULL
                  ))){

                    //If the nodes aren't child snarls

                    int64_t myDist = di.min_distance(pos1, pos2);
                    int64_t maxDist = di.max_distance(pos1, pos2);
                    int64_t actDist = min_distance(&graph, pos1, pos2);


                    int64_t dist2 = di.min_distance(pos1, make_pos_t(nodeID2, !is_rev(pos2), node2->sequence().size() - offset2 - 1) ); 
                    dist2 = myDist == -1 ? dist2 : (dist2 == -1 ? myDist : min(myDist, dist2));
                    int64_t dist3 = di.min_distance(make_pos_t(nodeID1, !is_rev(pos1), node1->sequence().size() - offset1 - 1), pos2 ); 
                    dist3 = dist2 == -1 ? dist3 : (dist3 == -1 ? dist2 : min(dist2, dist3));
                    int64_t dist4 = di.min_distance(make_pos_t(nodeID1, !is_rev(pos1), node1->sequence().size() - offset1 - 1), 
                                                   make_pos_t(nodeID2, !is_rev(pos2), node2->sequence().size() - offset2 - 1) ); 
                    dist4 = dist3 == -1 ? dist4 : (dist4 == -1 ? dist3 : min(dist3, dist4));

                    auto root_offset1 = di.get_minimizer_distances (pos1);
                    auto root_offset2 = di.get_minimizer_distances (pos2);
                    if (std::get<0>(root_offset1) && std::get<0>(root_offset2) &&
                        std::get<1>(root_offset1) == std::get<1>(root_offset2) &&
                        dist4 != -1) {
                        graph.serialize_to_file("testGraph");
                        //If these positions are both on the same root chain, then the minimum (unoriented) distance between them
                        //should be the difference between the offsets we store
                        if (std::get<2>(root_offset1) < std::get<2>(root_offset2)) {
                            REQUIRE(std::get<2>(root_offset2) - std::get<2>(root_offset1) == dist4);
                        } else {
                            REQUIRE(std::get<2>(root_offset1) - std::get<2>(root_offset2) == dist4);
                        }

                        
                    }
         

 
                    bool found = false;
                    pair<id_t, bool> next;
                    auto addFirst = [&](const handle_t& h) -> bool {
                        next = make_pair(graph.get_id(h), 
                                          graph.get_is_reverse(h));
                        found = true;
                        return false;
                    };
                    handle_t currHandle = graph.get_handle(nodeID1, false); 
                    graph.follow_edges(currHandle, false, addFirst);
 
                    bool passed = (myDist == actDist) &&
                                  (myDist > cap || maxDist > actDist);


                    if (!passed) { 
                        graph.serialize_to_file("testGraph");
                        di.print_self();
                        cerr << "Failed on random test: " << endl;
                        
                        cerr << "Position 1 on snarl " << 
                               snarl1->start().node_id() << " " << " Node " <<
                               nodeID1 << " is rev? " << is_rev(pos1) << 
                               " offset: " << offset1 << endl;
                        cerr << "Position 2 on snarl " << 
                              snarl2->start().node_id() << " Node " << nodeID2 
                              << " is rev? " << is_rev(pos2) << " offset: " 
                              << offset2 << endl;

                        cerr << "Actual min distance: " << actDist << "    " <<
                                "Guessed min distance: " << myDist << "    " <<
                                "Guessed max dist " << maxDist << endl;
                        cerr << endl;
                    }
                    REQUIRE(passed);
                }
            }
             
        }
    } //end test case

/*
    TEST_CASE("From serialized index", "[dist]"){

       ifstream vg_stream("primary-BRCA1.vg");
       VG vg(vg_stream);
       vg_stream.close();

       ifstream xg_stream("primary-BRCA1.xg");
       xg::XG xg(xg_stream);
       xg_stream.close();
 
       CactusSnarlFinder bubble_finder(vg);
       //SnarlManager snarl_manager1 = bubble_finder.find_snarls();
       SnarlManager* snarl_manager = nullptr; 
       ifstream snarl_stream("primary-snarls.pb");
       snarl_manager = new SnarlManager(snarl_stream);
       snarl_stream.close();





       ifstream dist_stream("primary-dist");
       
//       DistanceIndex di(&vg, snarl_manager);
       DistanceIndex di(&vg, snarl_manager, dist_stream);
       dist_stream.close();

       default_random_engine generator(test_seed_source());
       for (size_t i = 0; i < 1000 ; i++) {
    
           size_t maxPos = xg.seq_length;
           size_t offset1 = uniform_int_distribution<int>(1, maxPos)(generator);
           size_t offset2 = uniform_int_distribution<int>(1, maxPos)(generator); 
           id_t nodeID1 = xg.node_at_seq_pos(offset1);
           id_t nodeID2 = xg.node_at_seq_pos(offset2);

           pos_t pos1 = make_pos_t(nodeID1, false, 0);
           pos_t pos2 = make_pos_t(nodeID2, false, 0);
           REQUIRE(di.min_distance(pos1, pos2) == distance(&vg, pos1, pos2));
        }

    }
        
*/

    TEST_CASE("Serialize distance index", "[min_dist][serial]") {
        for (int i = 0; i < 0; i++) {

            VG graph;
            random_graph(1000, 20, 100, &graph);

            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls(); 

            MinimumDistanceIndex di (&graph, &snarl_manager);
  
            filebuf buf;
            ofstream out("distanceIndex");
           
            di.serialize(out);
            out.close();

            //buf.close();

            buf.open("distanceIndex", ios::in);
            istream in(&buf);
            MinimumDistanceIndex sdi (in);
            buf.close();
            #ifdef print        
                di.print_self();
                sdi.print_self();

            #endif
 
            vector<const Snarl*> allSnarls;
            auto addSnarl = [&] (const Snarl* s) {
                allSnarls.push_back(s);
            };
            snarl_manager.for_each_snarl_preorder(addSnarl);

            uniform_int_distribution<int> randSnarlIndex(0, allSnarls.size()-1);
            default_random_engine generator(time(NULL));
            for (int j = 0; j < 100; j++) {
                //Check distances for random pairs of positions 
                const Snarl* snarl1 = allSnarls[randSnarlIndex(generator)];
                const Snarl* snarl2 = allSnarls[randSnarlIndex(generator)];
                 
                pair<unordered_set<Node*>, unordered_set<Edge*>> contents1 = 
                    pb_contents(graph, snarl_manager.shallow_contents(snarl1, graph, true));
                pair<unordered_set<Node*>, unordered_set<Edge*>> contents2 = 
                    pb_contents(graph, snarl_manager.shallow_contents(snarl2, graph, true));
 
                vector<Node*> nodes1 (contents1.first.begin(), contents1.first.end());
                vector<Node*> nodes2 (contents2.first.begin(), contents2.first.end());

                
                uniform_int_distribution<int> randNodeIndex2(0,nodes2.size()-1);
                uniform_int_distribution<int> randNodeIndex1(0,nodes1.size()-1);
 
                Node* node1 = nodes1[randNodeIndex1(generator)];
                Node* node2 = nodes2[randNodeIndex2(generator)];
                id_t nodeID1 = node1->id();
                id_t nodeID2 = node2->id();
 
                off_t offset1 = uniform_int_distribution<int>(0,node1->sequence().size() - 1)(generator);
                off_t offset2 = uniform_int_distribution<int>(0,node2->sequence().size() - 1)(generator);

                pos_t pos1 = make_pos_t(nodeID1, 
                  uniform_int_distribution<int>(0,1)(generator) == 0,offset1 );
                pos_t pos2 = make_pos_t(nodeID2, 
                  uniform_int_distribution<int>(0,1)(generator) == 0, offset2 );
 
                if (!(nodeID1 != snarl1->start().node_id() && 
                    (snarl_manager.into_which_snarl(nodeID1, false) != NULL ||
                      snarl_manager.into_which_snarl(nodeID1, true) != NULL)) &&
                    ! (nodeID2 != snarl2->start().node_id() &&
                        (snarl_manager.into_which_snarl(nodeID2, false) != NULL
                       || snarl_manager.into_which_snarl(nodeID2, true) != NULL
                  ))){

                    //If the nodes aren't child snarls

                    int64_t myDist = di.min_distance(pos1, pos2);
                    int64_t serialDist = sdi.min_distance(pos1, pos2);

                    bool passed = myDist == serialDist;

                    if (!passed) { 
                        graph.serialize_to_file("testGraph");
                        di.print_self();
                        cerr << "Failed on random test: " << endl;
                        
                        cerr << "Position 1 on snarl " << 
                               snarl1->start().node_id() << " " << " Node " <<
                               nodeID1 << " is rev? " << is_rev(pos1) << 
                               " offset: " << offset1 << endl;
                        cerr << "Position 2 on snarl " << 
                              snarl2->start().node_id() << " Node " << nodeID2 
                              << " is rev? " << is_rev(pos2) << " offset: " 
                              << offset2 << endl;

                        cerr << "Serial distance: " << serialDist << "    " <<
                                "Guessed distance: " << myDist << endl;
                    }
                    REQUIRE(passed);
                }
            }
             
        }
    } //end test case
    /*
    TEST_CASE("Serialize only minimum distance index", "[dist][serial]") {
        for (int i = 0; i < 100; i++) {

            VG graph;
            random_graph(1000, 20, 100, &graph);

            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls(); 

            MinimumDistanceIndex di (&graph, &snarl_manager, 50, false);
  
            filebuf buf;
            ofstream out("distanceIndex");
           
            di.serialize(out);
            out.close();

            //buf.close();

            buf.open("distanceIndex", ios::in);
            istream in(&buf);
            MinimumDistanceIndex sdi (&graph, &snarl_manager, in);
            buf.close();
            #ifdef print        
                di.print_self();
                sdi.print_self();

            #endif
            REQUIRE (sdi.maxIndex.min_distances.size() == 0);
 
            vector<const Snarl*> allSnarls;
            auto addSnarl = [&] (const Snarl* s) {
                allSnarls.push_back(s);
            };
            snarl_manager.for_each_snarl_preorder(addSnarl);

            uniform_int_distribution<int> randSnarlIndex(0, allSnarls.size()-1);
            default_random_engine generator(time(NULL));
            for (int j = 0; j < 100; j++) {
                //Check distances for random pairs of positions 
                const Snarl* snarl1 = allSnarls[randSnarlIndex(generator)];
                const Snarl* snarl2 = allSnarls[randSnarlIndex(generator)];
                 
                pair<unordered_set<Node*>, unordered_set<Edge*>> contents1 = 
                           pb_contents(graph, snarl_manager.shallow_contents(snarl1, graph, true));
                pair<unordered_set<Node*>, unordered_set<Edge*>> contents2 = 
                           snarl_manager.shallow_contents(snarl2, graph, true);
 
                vector<Node*> nodes1 (contents1.first.begin(), contents1.first.end());
                vector<Node*> nodes2 (contents2.first.begin(), contents2.first.end());

                
                uniform_int_distribution<int> randNodeIndex2(0,nodes2.size()-1);
                uniform_int_distribution<int> randNodeIndex1(0,nodes1.size()-1);
 
                Node* node1 = nodes1[randNodeIndex1(generator)];
                Node* node2 = nodes2[randNodeIndex2(generator)];
                id_t nodeID1 = node1->id();
                id_t nodeID2 = node2->id();
 
                off_t offset1 = uniform_int_distribution<int>(0,node1->sequence().size() - 1)(generator);
                off_t offset2 = uniform_int_distribution<int>(0,node2->sequence().size() - 1)(generator);

                pos_t pos1 = make_pos_t(nodeID1, 
                  uniform_int_distribution<int>(0,1)(generator) == 0,offset1 );
                pos_t pos2 = make_pos_t(nodeID2, 
                  uniform_int_distribution<int>(0,1)(generator) == 0, offset2 );
 
                if (!(nodeID1 != snarl1->start().node_id() && 
                    (snarl_manager.into_which_snarl(nodeID1, false) != NULL ||
                      snarl_manager.into_which_snarl(nodeID1, true) != NULL)) &&
                    ! (nodeID2 != snarl2->start().node_id() &&
                        (snarl_manager.into_which_snarl(nodeID2, false) != NULL
                       || snarl_manager.into_which_snarl(nodeID2, true) != NULL
                  ))){

                    //If the nodes aren't child snarls

                    int64_t myDist = di.min_distance(pos1, pos2);
                    int64_t serialDist = sdi.min_distance(snarl1, snarl2,pos1, pos2);

                    bool passed = myDist == serialDist;

                    if (!passed) { 
                        graph.serialize_to_file("testGraph");
                        di.print_self();
                        cerr << "Failed on random test: " << endl;
                        
                        cerr << "Position 1 on snarl " << 
                               snarl1->start().node_id() << " " << " Node " <<
                               nodeID1 << " is rev? " << is_rev(pos1) << 
                               " offset: " << offset1 << endl;
                        cerr << "Position 2 on snarl " << 
                              snarl2->start().node_id() << " Node " << nodeID2 
                              << " is rev? " << is_rev(pos2) << " offset: " 
                              << offset2 << endl;

                        cerr << "Serial distance: " << serialDist << "    " <<
                                "Guessed distance: " << myDist << endl;
                    }
                    REQUIRE(passed);
                }
            }
             
        }
    } //end test case
*/
    TEST_CASE( "Simple snarl subgraph",
                   "[min_dist][min_subgraph]" ) {
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

        SECTION("Subgraph extraction") {

            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(2, false);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 5, 7, sub_graph, true); 

            REQUIRE(!sub_graph.count(3));
            REQUIRE(sub_graph.count(4));
            REQUIRE(sub_graph.count(5));
            REQUIRE(!sub_graph.count(6));
            REQUIRE(!sub_graph.count(7));
            REQUIRE(sub_graph.count(8));
        } 
        SECTION("Subgraph extraction same node") {

            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(3, false);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 4, 7, sub_graph, true); 

            REQUIRE(!sub_graph.count(3));
            REQUIRE(sub_graph.count(4));
            REQUIRE(sub_graph.count(5));
            REQUIRE(!sub_graph.count(6));
            REQUIRE(sub_graph.count(7));
            REQUIRE(sub_graph.count(8));
        } 
    } //end test case
    TEST_CASE("Chain subgraph", "[min_dist][min_subgraph]") {
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
        Edge* e4 = graph.create_edge(n5, n6);
        Edge* e5 = graph.create_edge(n2, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        SECTION("Subgraph extraction") {

            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(2, false);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 4, 7, sub_graph, true); 

            REQUIRE(!sub_graph.count(3));
            REQUIRE(sub_graph.count(4));
            REQUIRE(sub_graph.count(5));
            REQUIRE(sub_graph.count(6));
            REQUIRE(sub_graph.count(7));
            REQUIRE(sub_graph.count(8));
        }
        SECTION ("Another subgraph") {
            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(2, false);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 10, 10, sub_graph, true); 

            REQUIRE(!sub_graph.count(3));
            REQUIRE(!sub_graph.count(4));
            REQUIRE(!sub_graph.count(5));
            REQUIRE(!sub_graph.count(6));
            REQUIRE(!sub_graph.count(7));
            REQUIRE(sub_graph.count(8));
            REQUIRE(!sub_graph.count(3));
        }
        SECTION ("Skip snarl") {
            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(2, false);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 6, 6, sub_graph, true); 

            REQUIRE(!sub_graph.count(3));
            REQUIRE(!sub_graph.count(4));
            REQUIRE(!sub_graph.count(5));
            REQUIRE(sub_graph.count(6));
            REQUIRE(sub_graph.count(7));
            REQUIRE(!sub_graph.count(8));
            REQUIRE(!sub_graph.count(3));
        }

    }//end test case
    TEST_CASE("Top level chain subgraph", "[min_dist][min_subgraph]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("G");
        Node* n12 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n8);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e5 = graph.create_edge(n2, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e4 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n8, n9);
        Edge* e12 = graph.create_edge(n9, n10);
        Edge* e13 = graph.create_edge(n8, n10);
        Edge* e14 = graph.create_edge(n5, n5, true, false);
        Edge* e15 = graph.create_edge(n10, n11);
        Edge* e16 = graph.create_edge(n10, n12);
        Edge* e17 = graph.create_edge(n11, n12);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 


        SECTION("Skip right in chain") {

            MinimumDistanceIndex di (&graph, &snarl_manager);
            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(1, false);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 9, 9, sub_graph, true); 

            REQUIRE(!sub_graph.count(3));
            REQUIRE(!sub_graph.count(8));
            REQUIRE(!sub_graph.count(9));
            REQUIRE(!sub_graph.count(10));
            REQUIRE(sub_graph.count(11));
            REQUIRE(sub_graph.count(12));

        }

        SECTION("Skip right in root chain") {

            MinimumDistanceIndex di (&graph, &snarl_manager);
            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(2, false);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 12, 14, sub_graph, true); 

            REQUIRE(!sub_graph.count(3));
            REQUIRE(!sub_graph.count(8));
            REQUIRE(!sub_graph.count(9));
            REQUIRE(!sub_graph.count(10));
            REQUIRE(sub_graph.count(11));
            REQUIRE(sub_graph.count(12));

        }
        SECTION("Skip left in root chain") {

            MinimumDistanceIndex di (&graph, &snarl_manager);
            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(11, true);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 8, 11, sub_graph, true); 

            REQUIRE(sub_graph.count(1));
            REQUIRE(!sub_graph.count(2));
            REQUIRE(sub_graph.count(3));
            REQUIRE(sub_graph.count(4));
            REQUIRE(sub_graph.count(5));
            REQUIRE(sub_graph.count(6));
            REQUIRE(!sub_graph.count(7));

        }
        SECTION("Take loop") {

            MinimumDistanceIndex di (&graph, &snarl_manager);
            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(8, true);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 12, 13, sub_graph, true); 

            REQUIRE(!sub_graph.count(1));
            REQUIRE(!sub_graph.count(2));
            REQUIRE(!sub_graph.count(3));
            REQUIRE(sub_graph.count(4));
            REQUIRE(!sub_graph.count(5));
            REQUIRE(sub_graph.count(6));
            REQUIRE(sub_graph.count(7));
            REQUIRE(sub_graph.count(8));
            REQUIRE(!sub_graph.count(9));

        }
        SECTION("Take loop again") {

            MinimumDistanceIndex di (&graph, &snarl_manager);
            std::unordered_set<id_t> sub_graph;
            handle_t handle = graph.get_handle(8, true);
            path_handle_t path_handle = graph.create_path_handle("path");
            graph.append_step(path_handle, handle);
            Path path = path_from_path_handle(graph, path_handle);

            MinimumDistanceIndex dist_index (&graph, &snarl_manager);
            dist_index.subgraph_in_range(path, &graph, 16, 17, sub_graph, true); 

            REQUIRE(!sub_graph.count(1));
            REQUIRE(!sub_graph.count(2));
            REQUIRE(!sub_graph.count(3));
            REQUIRE(!sub_graph.count(4));
            REQUIRE(!sub_graph.count(5));
            REQUIRE(!sub_graph.count(6));
            REQUIRE(!sub_graph.count(7));
            REQUIRE(sub_graph.count(8));
            REQUIRE(sub_graph.count(9));
            REQUIRE(sub_graph.count(10));
            REQUIRE(!sub_graph.count(11));
            REQUIRE(!sub_graph.count(12));

        }
    }//end test case

    TEST_CASE("Random test subgraph", "[min_dist][min_subgraph][rand]") {

        int64_t min = 20; int64_t max = 50;

//        ifstream vg_stream("testGraph");
//        VG vg(vg_stream);
//        vg_stream.close();
//        CactusSnarlFinder bubble_finder(vg);
//        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
//
//        MinimumDistanceIndex di (&vg, &snarl_manager);
//        di.print_self();
//
//        handle_t handle = vg.get_handle(30, false);
//        path_handle_t path_handle = vg.create_path_handle("test_path");
//        vg.append_step(path_handle, handle);
//        Path path = path_from_path_handle(vg, path_handle);
//
//        std::unordered_set<id_t> sub_graph;
//        di.subgraph_in_range(path, &vg, min, max, sub_graph, false); 
//
//        REQUIRE(sub_graph.count(27));

        for (int i = 0; i < 0; i++) {
            //1000 different graphs
            VG graph;
            random_graph(1000, 10, 15, &graph);

            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls(); 
            MinimumDistanceIndex dist_index (&graph, &snarl_manager);

            vector<const Snarl*> allSnarls;
            auto addSnarl = [&] (const Snarl* s) {
                allSnarls.push_back(s);
            };
            snarl_manager.for_each_snarl_preorder(addSnarl);

            uniform_int_distribution<int> randSnarlIndex(0, allSnarls.size()-1);
            default_random_engine generator(time(NULL));
            for (int j = 0; j < 100; j++) {
                //Check distances for random pairs of positions 
                const Snarl* snarl1 = allSnarls[randSnarlIndex(generator)];
                 
                pair<unordered_set<Node*>, unordered_set<Edge*>> contents1 = 
                    pb_contents(graph, snarl_manager.shallow_contents(snarl1, graph, true));
 
                vector<Node*> nodes1 (contents1.first.begin(), contents1.first.end());

                uniform_int_distribution<int> randNodeIndex1(0,nodes1.size()-1);
 
                Node* node1 = nodes1[randNodeIndex1(generator)];
                id_t nodeID1 = node1->id();
                handle_t handle = graph.get_handle(nodeID1, false);
                path_handle_t path_handle = graph.create_path_handle("test_path");
                graph.prepend_step(path_handle, handle);
                Path path = path_from_path_handle(graph, path_handle);
                pos_t pos1 = make_pos_t(nodeID1, false, 0 );

                std::unordered_set<id_t> sub_graph;
                size_t node_len = graph.get_length(handle);
                dist_index.subgraph_in_range(path, &graph, min+1, max+1, sub_graph, true); 

                graph.for_each_handle([&] (const handle_t h ) { 
                    id_t node_id = graph.get_id(h);
                    int64_t len = graph.get_length(h);
                    int64_t dist_start_fd = dist_index.min_distance(pos1, make_pos_t(node_id, false, 0));
                    int64_t dist_end_fd = dist_start_fd == -1 ? -1 : dist_start_fd + len - 1;

                    bool start_forward = dist_start_fd != -1 && (dist_start_fd >= min && dist_start_fd <= max);
                    bool end_forward = dist_end_fd != -1 && (dist_end_fd >= min && dist_end_fd <= max);
                    bool in_forward = dist_start_fd != -1 && dist_end_fd == -1 || (dist_start_fd <= min && dist_end_fd >= max);
                
                    int64_t dist_start_bk = dist_index.min_distance(pos1, make_pos_t(node_id, true, 0));
                    int64_t dist_end_bk = dist_start_bk == -1 ? -1 : dist_start_bk + len - 1;

                    bool start_backward = dist_start_bk != -1 && (dist_start_bk >= min && dist_start_bk <= max);
                    bool end_backward = dist_end_bk != -1 && (dist_end_bk >= min && dist_end_bk <= max);
                    bool in_backward = dist_start_bk != -1 && dist_end_bk == -1 || (dist_start_bk <= min && dist_end_bk >= max);
                    if (sub_graph.count(node_id)) {
                        //If this node is in the subgraph, then the node must be within the range

                        if (!(start_forward || end_forward || in_forward || start_backward || end_backward || in_backward)) {
                            cerr << "Node " << node_id << " from pos " << pos1 << " with distances " 
                                 << dist_index.min_distance(pos1, make_pos_t(node_id, false, 0)) << " and " 
                                 << dist_index.min_distance(pos1, make_pos_t(node_id, true, 0)) 
                                 << " (" << dist_start_fd << " " << dist_end_fd << " " << dist_start_bk << " " << dist_end_bk << ") "
                                 << " is in the subgraph but shouldn't be " << endl;
                            graph.serialize_to_file("testGraph");
                        }
                        REQUIRE((start_forward || end_forward || in_forward || start_backward || end_backward || in_backward));
                    } else {
                        if ((start_forward || end_forward || in_forward || start_backward || end_backward || in_backward) &&
                             node_id != get_id(pos1)) {
                            cerr << "Node " << node_id << " from pos " << pos1 <<" with distances " 
                                 << dist_index.min_distance(pos1, make_pos_t(node_id, false, 0)) << " and " 
                                 << dist_index.min_distance(pos1, make_pos_t(node_id, true, 0)) 
                                 << " (" << dist_start_fd << " " << dist_end_fd << " " << dist_start_bk << " " << dist_end_bk << ") "
                                 << " is not in the subgraph but should be " << endl;
                            graph.serialize_to_file("testGraph");
                            REQUIRE(!(start_forward || end_forward || in_forward || start_backward || end_backward || in_backward));
                        }
                    }
                });
            }
        }
    }//End test case
}
}
