#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include "json2pb.h"
#include "vg.hpp"
#include "catch.hpp"
#include "snarls.hpp"
#include "min_distance.hpp"
#include "genotypekit.hpp"
#include "random_graph.hpp"
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
            di.printSelf();
#endif

            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(2, false, 0) ) == 3);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0) ) == 3);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0) ) == 4);
            REQUIRE(di.minDistance(make_pos_t(2, false, 0),
                                   make_pos_t(6, false, 0) ) == 1);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0) ) == 5);

        } 
        SECTION( "Create max distance index" ) {

            MinimumDistanceIndex di (&graph, &snarl_manager, 20);
#ifdef print
            di.printSelf();
#endif

            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(2, false, 0) ) == 3);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0) ) == 3);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0) ) == 4);
            REQUIRE(di.minDistance(make_pos_t(2, false, 0),
                                   make_pos_t(6, false, 0) ) == 1);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0) ) == 5);

            REQUIRE(di.maxDistance(make_pos_t(1, false, 0),
                                   make_pos_t(2, false, 0) ) >= 3);
            REQUIRE(di.maxDistance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0) ) >= 13);
            REQUIRE(di.maxDistance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0) ) >= 5);
            REQUIRE(di.maxDistance(make_pos_t(2, false, 0),
                                   make_pos_t(6, false, 0) ) >= 1);
            REQUIRE(di.maxDistance(make_pos_t(1, false, 0),
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
            di.printSelf();
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

            REQUIRE(di.minDistance( pos1, pos4) == 4);
            REQUIRE(di.minDistance( pos1, pos6) == 7);
            REQUIRE(di.minDistance( pos6, pos7) == -1);

            REQUIRE(di.minDistance( pos9, pos12) == 5);
            REQUIRE(di.minDistance(pos9, pos11) == 1);
            REQUIRE(di.minDistance( pos4, pos9) == -1);
            REQUIRE(di.minDistance( pos4, pos12) == -1);
            REQUIRE(di.minDistance( pos1, pos11) == -1);


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
            di.printSelf();
#endif

            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) == 5);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0)) == 8);
            REQUIRE(di.minDistance(make_pos_t(6, true, 0),
                                   make_pos_t(1, true, 0)) == 6);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(6, false, 0)) == 4);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) == 4);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(8, false, 0)) == 5);

        }

        SECTION("Max distance") {

            MinimumDistanceIndex di (&graph, &snarl_manager, 20);
#ifdef print
            di.printSelf();
#endif

            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) == 5);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0)) == 8);
            REQUIRE(di.minDistance(make_pos_t(6, true, 0),
                                   make_pos_t(1, true, 0)) == 6);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(6, false, 0)) == 4);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) == 4);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(8, false, 0)) == 5);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) == 3);

            REQUIRE(di.maxDistance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) >= 8);
            REQUIRE(di.maxDistance(make_pos_t(1, false, 0),
                                   make_pos_t(6, false, 0)) >= 11);
            REQUIRE(di.maxDistance(make_pos_t(6, true, 0),
                                   make_pos_t(1, true, 0)) >= 9);
            REQUIRE(di.maxDistance(make_pos_t(3, false, 0),
                                   make_pos_t(6, false, 0)) >= 4);
            REQUIRE(di.maxDistance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) >= 5);
            REQUIRE(di.maxDistance(make_pos_t(3, false, 0),
                                   make_pos_t(8, false, 0)) >= 6);
            REQUIRE(di.maxDistance(make_pos_t(1, false, 0),
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
            di.printSelf();
            #endif
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) == 3);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) == 4);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(4, true, 0)) == 7);
            REQUIRE(di.minDistance(make_pos_t(5, false, 0),
                                   make_pos_t(6, true, 0)) == 5);



        }
        SECTION("Max distance") {
            MinimumDistanceIndex di (&graph, &snarl_manager, 30);
            #ifdef print
            di.printSelf();
            #endif
            REQUIRE(di.maxDistance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) >= 13);
            REQUIRE(di.maxDistance(make_pos_t(3, false, 0),
                                   make_pos_t(7, false, 0)) >= 5);
            REQUIRE(di.maxDistance(make_pos_t(3, false, 0),
                                   make_pos_t(4, true, 0)) >= 27);
            REQUIRE(di.maxDistance(make_pos_t(5, false, 0),
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
            di.printSelf();
            #endif
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0)) == 8);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(10, false, 0)) == 7);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) == 5);
            REQUIRE(di.minDistance(make_pos_t(2, false, 0),
                                   make_pos_t(9, false, 0)) == 10);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
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
            di.printSelf();
            #endif
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0)) == 8);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(8, false, 0)) == 9);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(8, true, 0)) == 9);
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(5, false, 0)) == 8);
            REQUIRE(di.minDistance(make_pos_t(8, false, 0),
                                   make_pos_t(5, true, 0)) == 5);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(6, false, 0)) == 5);
            REQUIRE(di.minDistance(make_pos_t(3, false, 0),
                                   make_pos_t(10, false, 0)) == 11);
            REQUIRE(di.minDistance(make_pos_t(7, false, 0),
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
            di.printSelf();
            #endif
            REQUIRE(di.minDistance(make_pos_t(1, false, 0),
                                   make_pos_t(7, false, 0)) == 4);
            REQUIRE(di.minDistance(make_pos_t(2, true, 0),
                                   make_pos_t(3, false, 0)) == 7);
            REQUIRE(di.minDistance(make_pos_t(4, true, 0),
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
            REQUIRE(di.minDistance(make_pos_t(5, false, 1), make_pos_t(5, false, 0)) == 11); 
            REQUIRE(di.minDistance(make_pos_t(5, false, 0), make_pos_t(5, false, 0)) == 0); 
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

            REQUIRE(di.minDistance(make_pos_t(1, false, 0), make_pos_t(4, false, 0)) == 4);
            REQUIRE(di.minDistance( make_pos_t(5, false, 0), make_pos_t(6, false, 0)) == -1);
            REQUIRE(di.minDistance( make_pos_t(5, false, 0), make_pos_t(6, true, 0)) == -1);
            REQUIRE(di.minDistance( make_pos_t(5, true, 2), make_pos_t(6, false, 0)) == 11);
            REQUIRE(di.minDistance(make_pos_t(2, false, 0), make_pos_t(7, false, 0)) == 6);

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
            di.printSelf();
#endif

            REQUIRE(di.minDistance( make_pos_t(2, true, 0), 
                                    make_pos_t(9, true, 1)) == 2);
            REQUIRE(di.minDistance( make_pos_t(2, false, 0), 
                                    make_pos_t(9, true, 1)) == 5);
            REQUIRE(di.minDistance( make_pos_t(3, true, 0),
                                    make_pos_t(9, false, 0)) == 4);
            REQUIRE(di.minDistance( make_pos_t(3, true, 0), 
                                    make_pos_t(9, true, 1)) == 3);
            REQUIRE(di.minDistance( make_pos_t(3, false, 0), 
                                    make_pos_t(9, false, 0)) == 11);
            REQUIRE(di.minDistance( make_pos_t(4, false, 0), 
                                    make_pos_t(5, false, 0)) == 14);
            REQUIRE(di.minDistance( make_pos_t(4, true, 0), 
                                    make_pos_t(5, false, 0)) == 8);
            REQUIRE(di.minDistance(make_pos_t(7, false, 0), 
                                   make_pos_t(12, false, 0)) == 14);
            REQUIRE(di.minDistance(make_pos_t(7, false, 0), 
                                   make_pos_t(12, true, 0)) == 13);
            REQUIRE(di.minDistance(make_pos_t(7, true, 0),
                                   make_pos_t(12, false, 0)) == 15);
            REQUIRE(di.minDistance(make_pos_t(8, false, 0), 
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
            di.printSelf();
#endif

            REQUIRE(di.minDistance(make_pos_t(2, true, 0), 
                                   make_pos_t(10, true, 0)) == 2);
            REQUIRE(di.minDistance(make_pos_t(2, false, 0), 
                                   make_pos_t(10, false, 0)) == 4);
            REQUIRE(di.minDistance(make_pos_t(4, true, 3), 
                                   make_pos_t(5, false, 0)) == 5);
            REQUIRE(di.minDistance(make_pos_t(4, false, 0),
                                   make_pos_t(5, false, 0)) == 11);
            REQUIRE(di.minDistance(make_pos_t(4, true, 3),
                                   make_pos_t(5, true, 0)) == 8);
            REQUIRE(di.minDistance(make_pos_t(14, false, 0),
                                   make_pos_t(10, false, 0)) == 10);
            REQUIRE(di.minDistance(make_pos_t(14, false, 0), 
                                   make_pos_t(10, true, 0)) == 5);
            REQUIRE(di.minDistance(make_pos_t(14, false, 0), 
                                   make_pos_t(3, false, 0)) == 7);
            REQUIRE(di.minDistance(make_pos_t(16, false, 0),
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
                di.printSelf();
            #endif

            REQUIRE(di.minDistance(make_pos_t(4, false, 0), 
                                   make_pos_t(6, false, 0)) == 19);
            REQUIRE(di.minDistance(make_pos_t(4, false, 0), 
                                   make_pos_t(6, true, 0)) == 25);
            REQUIRE(di.minDistance(make_pos_t(2, true, 0), 
                                   make_pos_t(7, false, 0)) == 7);
            REQUIRE(di.minDistance(make_pos_t(2, false, 0), 
                                   make_pos_t(7, false, 0)) == 6);
            REQUIRE(di.minDistance(make_pos_t(2, false, 0), 
                                   make_pos_t(7, true, 0)) == 11);

            REQUIRE(min_distance(&graph, make_pos_t(4, false, 0),  
                                         make_pos_t(6, true, 0)) == 25);
            REQUIRE(min_distance(&graph, make_pos_t(2, false, 0),
                                         make_pos_t(7, true, 0)) == 11);
        }
        SECTION ("Max distance") {

            MinimumDistanceIndex di (&graph, &snarl_manager, 50);
            #ifdef print
                di.printSelf();
            #endif

            REQUIRE(di.maxDistance(make_pos_t(4, false, 0), 
                                   make_pos_t(6, false, 0)) >= 40);
            REQUIRE(di.maxDistance(make_pos_t(4, false, 0), 
                                   make_pos_t(6, true, 0)) >= 25);
            REQUIRE(di.maxDistance(make_pos_t(2, false, 0), 
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
                di.printSelf();
            #endif
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, true);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);

            REQUIRE(di.minDistance(pos1, pos2) == min_distance(&graph, pos1, pos2));



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
        di.printSelf();
        pos_t pos1 = make_pos_t(208, true, 4);
        pos_t pos2 = make_pos_t(256, false, 2); 
        REQUIRE(di.minDistance(pos1, pos2) == min_distance(&vg, pos1, pos2));
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
                di.printSelf();

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

                    int64_t myDist = di.minDistance(pos1, pos2);
                    int64_t maxDist = di.maxDistance(pos1, pos2);
                    int64_t actDist = min_distance(&graph, pos1, pos2);


         

 
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
                        di.printSelf();
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

       random_device seed_source;
       default_random_engine generator(seed_source());
       for (size_t i = 0; i < 1000 ; i++) {
    
           size_t maxPos = xg.seq_length;
           size_t offset1 = uniform_int_distribution<int>(1, maxPos)(generator);
           size_t offset2 = uniform_int_distribution<int>(1, maxPos)(generator); 
           id_t nodeID1 = xg.node_at_seq_pos(offset1);
           id_t nodeID2 = xg.node_at_seq_pos(offset2);

           pos_t pos1 = make_pos_t(nodeID1, false, 0);
           pos_t pos2 = make_pos_t(nodeID2, false, 0);
           REQUIRE(di.minDistance(pos1, pos2) == distance(&vg, pos1, pos2));
        }

    }
        
*/

    TEST_CASE("Serialize distance index", "[dist][serial]") {
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
                di.printSelf();
                sdi.printSelf();

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

                    int64_t myDist = di.minDistance(pos1, pos2);
                    int64_t serialDist = sdi.minDistance(pos1, pos2);

                    bool passed = myDist == serialDist;

                    if (!passed) { 
                        graph.serialize_to_file("testGraph");
                        di.printSelf();
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
                di.printSelf();
                sdi.printSelf();

            #endif
            REQUIRE (sdi.maxIndex.minDistances.size() == 0);
 
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

                    int64_t myDist = di.minDistance(pos1, pos2);
                    int64_t serialDist = sdi.minDistance(snarl1, snarl2,pos1, pos2);

                    bool passed = myDist == serialDist;

                    if (!passed) { 
                        graph.serialize_to_file("testGraph");
                        di.printSelf();
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
}

}
