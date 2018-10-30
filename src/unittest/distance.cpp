#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include "json2pb.h"
#include "vg.hpp"
#include "catch.hpp"
#include "snarls.hpp"
#include "distance.hpp"
#include "genotypekit.hpp"
#include "random_graph.hpp"
#include <fstream>
#include <random>
#include <time.h> 

//#define print

namespace vg {
namespace unittest {

int64_t loopDistance(VG* graph, const Snarl* snarl1, const Snarl* snarl2, pair<id_t, bool> node1, pair<id_t, bool> node2){
    //Loop distance using djikstras algorithm

    auto cmp = [] (pair<pair<id_t, bool> , int64_t> x, 
                   pair<pair<id_t, bool>, int64_t> y ) {
        return (x.second > y.second);
    };
 
    int64_t shortestDistance = -1;


    priority_queue< pair<pair<id_t, bool> , int64_t>, 
                    vector<pair<pair<id_t, bool>, int64_t>>,
                          decltype(cmp)> reachable(cmp); 
    handle_t currHandle = graph->get_handle(node2.first, node2.second);

    int64_t dist;
    dist = graph->get_length(currHandle);

    auto addFirst = [&](const handle_t& h) -> bool {
        pair<id_t, bool> node = make_pair(graph->get_id(h), 
                                          graph->get_is_reverse(h));
        reachable.push(make_pair(node, dist));
        return true;
    };
  
    graph->follow_edges(currHandle, false, addFirst);
    unordered_set<pair<id_t, bool>> seen;
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

        if (currID == node1){
        //Dist is distance to beginning or end of node containing pos2
            if (node1 != node2){        
                dist = dist +  graph->get_node(node1.first)->sequence().size();
            }
            if (shortestDistance == -1) {shortestDistance = dist;}
            else {shortestDistance = min(dist, shortestDistance);}
        }

    }
    return shortestDistance;
};


int64_t distance(VG* graph, pos_t pos1, pos_t pos2){
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
class TestDistanceIndex : public DistanceIndex {

    public:
        using DistanceIndex::DistanceIndex;
        using DistanceIndex::SnarlIndex;
        using DistanceIndex::MaxDistanceIndex;
        using DistanceIndex::maxIndex;
        using DistanceIndex::snarlDistances;
        using DistanceIndex::ChainIndex;
        using DistanceIndex::loopDistance;
        using DistanceIndex::chainDistances;
        using DistanceIndex::distToCommonAncestor;
        using DistanceIndex::checkChainDist;
        using DistanceIndex::checkChainLoopFd;
        using DistanceIndex::checkChainLoopRev;
        using DistanceIndex::printSelf;
};



    TEST_CASE( "Create distance index for simple nested snarl",
                   "[dist]" ) {
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

        TestDistanceIndex di (&graph, &snarl_manager, 20);
        
        SECTION( "Create distance index" ) {

#ifdef print
            di.printSelf();
#endif
            const vector<const Snarl*> topSnarls = snarl_manager.top_level_snarls();
            pair<id_t, bool> start1 = make_pair(snarl1->start().node_id(), 
                                                    snarl1->start().backward());
            pair<id_t, bool> end1 = make_pair(snarl1->end().node_id(), 
                                                    snarl1->end().backward());
            pair<id_t, bool> end1r = make_pair(snarl1->end().node_id(), 
                                                    !snarl1->end().backward());
            TestDistanceIndex::SnarlIndex& sd1 = di.snarlDistances.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));

            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlIndex& sd2 = di.snarlDistances.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));

            pair<id_t, bool> start3 = make_pair(snarl3->start().node_id(), 
                                                    snarl3->start().backward());
            pair<id_t, bool> start3r = make_pair(snarl3->start().node_id(), 
                                                   !snarl3->start().backward());
            pair<id_t, bool> end3 = make_pair(snarl3->end().node_id(), 
                                                    snarl3->end().backward());
            pair<id_t, bool> end3r = make_pair(snarl3->end().node_id(), 
                                                    !snarl3->end().backward());
            TestDistanceIndex::SnarlIndex& sd3 = di.snarlDistances.at(make_pair(snarl3->start().node_id(),
                                                 snarl3->start().backward()));

            NetGraph ng1 = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);
            NetGraph ng2 = NetGraph(snarl2->start(), snarl2->end(), snarl_manager.chains_of(snarl2), &graph);
            NetGraph ng3 = NetGraph(snarl3->start(), snarl3->end(), snarl_manager.chains_of(snarl3), &graph);

            REQUIRE(sd1.snarlDistance(&graph, &ng1,start1, make_pair(2, false)) == 3); 
            REQUIRE(sd1.snarlDistance(&graph, &ng1,start1, make_pair(2, true)) == -1); 
            REQUIRE(sd1.snarlDistanceShort(start1, make_pair(2, false)) == 0); 
            REQUIRE(sd1.snarlDistance(&graph , &ng1, start2r, make_pair(1,true)) == 3); 
            REQUIRE(sd1.snarlDistanceShort(make_pair(8, true), make_pair(1, true)) == 0); 
            REQUIRE(sd1.snarlDistanceShort(make_pair(1, false), make_pair(8, false)) == 0); 
 

             REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(3, false), make_pair(7, false)) == 4);          
             REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(7, true), make_pair(2, true)) == 2);

             REQUIRE(((sd2.snarlDistance(&graph, &ng2, end2r, start3r) == 1 && //3 is start
                      sd2.snarlDistance(&graph, &ng2, start2, start3) == 1) || 
                     (sd2.snarlDistance(&graph, &ng2, end2r, start3) == 1 &&   //5 is start
                      sd2.snarlDistance(&graph, &ng2, start2, start3r) == 1)
              ));

             REQUIRE(((sd1.snarlDistance(&graph, &ng1, end1r, start2r) == 4 && //3 is start
                      sd1.snarlDistance(&graph, &ng1, start1, start2) == 3) || 
                     (sd1.snarlDistance(&graph, &ng1, end1r, start2) == 4 &&   //5 is start
                      sd1.snarlDistance(&graph, &ng1, start1, start2r) == 3)
              ));

             REQUIRE(sd3.snarlDistance(&graph, &ng3, make_pair(5, true), make_pair(4, true))
                                                                          == 3);
             REQUIRE((sd2.snarlDistance(&graph, &ng2, make_pair(2, false), start3) == 1 ||
                       sd2.snarlDistance(&graph, &ng2, make_pair(2, false), start3r) == 1 ) );
             REQUIRE(( sd2.snarlDistance(&graph, &ng2, start2, start3) == 1 || 
                       sd2.snarlDistance(&graph, &ng2, start2, start3r) == 1)   );

             REQUIRE(( sd2.snarlDistance(&graph, &ng2, end2r, start3) == 1 || 
                       sd2.snarlDistance(&graph, &ng2, end2r, start3r) == 1)   );
            REQUIRE(di.chainDistances.size() == 0);
        }


        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

            pos_t pos1 = make_pos_t(4, false, 1);
            pos_t pos2 = make_pos_t(5, false, 2);
            pos_t pos3 = make_pos_t(5, false, 1);
            pos_t pos4 = make_pos_t(4, false, 2);
            pos_t pos5 = make_pos_t(2, false, 0);
            pos_t pos6 = make_pos_t(2, true, 0);
            pos_t pos7 = make_pos_t(1, true, 1);
            pos_t pos8r = make_pos_t(8, false, 2);
            pos_t pos8 = make_pos_t(8, true, 1);
            pos_t pos9 = make_pos_t(6, true, 0);

            REQUIRE(di.minDistance(pos3, pos2) == 1);

            REQUIRE(di.minDistance(pos1,pos2) == 5);

            REQUIRE(di.minDistance(pos5,pos2) == 4);
            REQUIRE(di.minDistance(pos6,pos2) == -1);
            REQUIRE(di.minDistance(pos2, pos5) == -1);
            REQUIRE(di.minDistance( pos7, pos8r) == -1);
            REQUIRE(di.minDistance( pos8, pos7) == 4);
            REQUIRE(di.minDistance( pos7, pos8) == -1);
            REQUIRE(di.minDistance(pos4, pos9) == -1);


            REQUIRE(distance(&graph, pos1,pos2) == 5);

            REQUIRE(distance(&graph, pos5,pos2) == 4);
            REQUIRE(distance(&graph, pos6,pos2) == -1);
            REQUIRE(distance(&graph, pos2, pos5) == -1);
            REQUIRE(distance( &graph, pos7, pos8r) == -1);
            REQUIRE(distance( &graph, pos8, pos7) == 4);
            REQUIRE(distance( &graph, pos7, pos8) == -1);
            REQUIRE(distance(&graph, pos4, pos9) == -1);
 
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(loopDistance(&graph, snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(loopDistance(&graph, snarl2, snarl2, make_pair(2, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(3, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(loopDistance(&graph, snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(4, false),  
                                        make_pair(4, false)) == -1);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(4, false),  
                                        make_pair(5, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(8, false),  
                                        make_pair(8, false)) == -1);

/*
            for (auto x : di.maxIndex.nodeToComponent) {
                REQUIRE(x == 1);
            }
*/
        }
    }//End test case

    TEST_CASE( "Create distance index disconnected graph",
                   "[dist]" ) {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);
        
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


    TEST_CASE("Simple chain", "[dist]") {
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

        const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
        const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
        const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);



        TestDistanceIndex di (&graph, &snarl_manager, 20);


        SECTION("Create distance index") {
            const Chain* chain = snarl_manager.chain_of(snarl2);

            pair<id_t, bool> start1 = make_pair(snarl1->start().node_id(), 
                                                    snarl1->start().backward());
            pair<id_t, bool> start1r = make_pair(snarl1->start().node_id(), 
                                                   !snarl1->start().backward());
            pair<id_t, bool> end1 = make_pair(snarl1->end().node_id(), 
                                                    snarl1->end().backward());
            pair<id_t, bool> end1r = make_pair(snarl1->end().node_id(), 
                                                    !snarl1->end().backward());
            TestDistanceIndex::SnarlIndex& sd1 = di.snarlDistances.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlIndex& sd2 = di.snarlDistances.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
 
            pair<id_t, bool> start5 = make_pair(snarl5->start().node_id(), 
                                                    snarl5->start().backward());
            pair<id_t, bool> start5r = make_pair(snarl5->start().node_id(), 
                                                   !snarl5->start().backward());
            pair<id_t, bool> end5 = make_pair(snarl5->end().node_id(), 
                                                    snarl5->end().backward());
            pair<id_t, bool> end5r = make_pair(snarl5->end().node_id(), 
                                                   !snarl5->end().backward());
            TestDistanceIndex::SnarlIndex& sd5 = di.snarlDistances.at(make_pair(snarl5->start().node_id(),
                                     snarl5->start().backward()));

            #ifdef print
            di.printSelf();
            #endif
            NetGraph ng = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 0) == 0);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 1) == 2);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 2) == 5);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 0) == -1);

            REQUIRE(sd1.snarlDistance(&graph, &ng, make_pair(1, false), make_pair(8, false))
                                                                         == 3);
            REQUIRE(sd1.snarlDistance(&graph, &ng, make_pair(8, true), make_pair(1, true))
                                                                          == 4);
            REQUIRE(di.chainDistances.size() == 1);
    
            TestDistanceIndex::ChainIndex& cd = di.chainDistances.at(get_start_of(*chain).node_id());
            REQUIRE(cd.chainDistance(start2, start2, snarl2, snarl2) == 0);
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(5, false), snarl2, snarl5) == 2);
            REQUIRE(cd.chainDistance(make_pair(5, true), make_pair(2, true), snarl5, snarl2) == 4);

            }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);

            pos_t pos1 = make_pos_t(4, false, 1);
            pos_t pos2 = make_pos_t(6, false, 0);
            pos_t pos3 = make_pos_t(2, false, 0);
            pos_t pos4 = make_pos_t(7, false, 0);
            pos_t pos5 = make_pos_t(1, false, 0);
            pos_t pos6 = make_pos_t(8, true, 0);

            REQUIRE(di.minDistance(pos1, pos2) == 6);
            REQUIRE(di.minDistance(pos2, pos1) == -1);
            REQUIRE(di.minDistance(pos3, pos1) == 2);
            REQUIRE(di.minDistance( pos1, pos3) == -1);
            REQUIRE(di.minDistance(pos3, pos2) == 5);
            REQUIRE(di.minDistance(pos3, pos4) == 5);
            REQUIRE(di.minDistance(pos5, pos2) == 8);
            REQUIRE(di.minDistance(pos5, pos6) == -1);

            REQUIRE(distance(&graph, pos1, pos2) == 6);
            REQUIRE(distance(&graph, pos2, pos1) == -1);
            REQUIRE(distance(&graph, pos3, pos1) == 2);
            REQUIRE(distance(&graph, pos1, pos3) == -1);
            REQUIRE(distance(&graph, pos3, pos2) == 5);
            REQUIRE(distance(&graph, pos3, pos4) == 5);
            REQUIRE(distance(&graph, pos5, pos2) == 8);
            REQUIRE(distance(&graph, pos5, pos6) == -1);


            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(3, true),  
                                        make_pair(3, true)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(3, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(loopDistance(&graph, snarl5, snarl5, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl5, snarl5, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl5, snarl5, make_pair(5, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(5, false),  
                                        make_pair(5, false)) == -1);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(8, false),  
                                        make_pair(8, false)) == -1);
/*
            for (auto x : di.maxIndex.nodeToComponent) {
                REQUIRE(x == 1);
            }
*/
        }
    }//end test case

    TEST_CASE("Chain with reversing edge", "[dist]") {
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

        const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
        const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
        const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);



        TestDistanceIndex di (&graph, &snarl_manager, 20);


        SECTION("Create distance index") {
            const Chain* chain = snarl_manager.chain_of(snarl2);

            pair<id_t, bool> start1 = make_pair(snarl1->start().node_id(), 
                                                    snarl1->start().backward());
            pair<id_t, bool> start1r = make_pair(snarl1->start().node_id(), 
                                                   !snarl1->start().backward());
            pair<id_t, bool> end1 = make_pair(snarl1->end().node_id(), 
                                                    snarl1->end().backward());
            pair<id_t, bool> end1r = make_pair(snarl1->end().node_id(), 
                                                    !snarl1->end().backward());
            TestDistanceIndex::SnarlIndex& sd1 = di.snarlDistances.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlIndex& sd2 = di.snarlDistances.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
 
            pair<id_t, bool> start5 = make_pair(snarl5->start().node_id(), 
                                                    snarl5->start().backward());
            pair<id_t, bool> start5r = make_pair(snarl5->start().node_id(), 
                                                   !snarl5->start().backward());
            pair<id_t, bool> end5 = make_pair(snarl5->end().node_id(), 
                                                    snarl5->end().backward());
            pair<id_t, bool> end5r = make_pair(snarl5->end().node_id(), 
                                                   !snarl5->end().backward());
            TestDistanceIndex::SnarlIndex& sd5 = di.snarlDistances.at(make_pair(snarl5->start().node_id(),
                                     snarl5->start().backward()));

            NetGraph ng1 = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);
            NetGraph ng2 = NetGraph(snarl2->start(), snarl2->end(), snarl_manager.chains_of(snarl2), &graph);
            NetGraph ng5 = NetGraph(snarl5->start(), snarl5->end(), snarl_manager.chains_of(snarl5), &graph);
    
            TestDistanceIndex::ChainIndex& cd = di.chainDistances.at(get_start_of(*chain).node_id());
            REQUIRE(cd.chainDistance(start2, start2, snarl2, snarl2) == 0);
            #ifdef print
            di.printSelf();
            #endif
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 0) == 0);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 1) == 2);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 2) == 5);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 0) == 9);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 1) == 3);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 2) == 1);
            REQUIRE(di.checkChainLoopRev(get_start_of(*chain).node_id(), 0) == -1);
            REQUIRE(di.checkChainLoopRev(get_start_of(*chain).node_id(), 1) == -1);

            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(1, false), make_pair(8, false))
                                                                         == 3);
            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(8, true), make_pair(1, true))
                                                                          == 4);
            REQUIRE(sd5.snarlDistance(&graph, &ng5, make_pair(5, false), make_pair(5, true))
                                                                          == 3);
            REQUIRE(sd5.snarlDistance(&graph, &ng5, make_pair(6, false), make_pair(6, true))
                                                                        == 12);
            REQUIRE(sd5.snarlDistance(&graph, &ng5, make_pair(5, false), make_pair(6, true))
                                                                        == 5);
            REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(3, false), make_pair(4, true))
                                                                         == 7);
            REQUIRE(di.chainDistances.size() == 1);
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(5, false), snarl2, snarl5) == 2);
            REQUIRE(cd.chainDistance(make_pair(5, true), make_pair(2, true), snarl5, snarl2) == 4);

            }
        SECTION ("Distance functions") {

            pos_t pos1 = make_pos_t(4, false, 1);
            pos_t pos2 = make_pos_t(6, false, 0);
            pos_t pos3 = make_pos_t(2, false, 0);
            pos_t pos4 = make_pos_t(7, false, 0);
            pos_t pos5 = make_pos_t(1, false, 0);
            pos_t pos6 = make_pos_t(8, true, 0);
            pos_t pos7 = make_pos_t(6, false, 10);
            pos_t pos9 = make_pos_t(6, true, 0);
            pos_t pos8 = make_pos_t(3, false, 0);

            REQUIRE(di.minDistance(pos1, pos2) == 6);
            REQUIRE(di.minDistance(pos2, pos1) == -1);
            REQUIRE(di.minDistance(pos3, pos1) == 2);
            REQUIRE(di.minDistance(pos1, pos3) == -1);
            REQUIRE(di.minDistance(pos3, pos2) == 5);
            REQUIRE(di.minDistance(pos3, pos4) == 5);
            REQUIRE(di.minDistance(pos5, pos2) == 8);
            REQUIRE(di.minDistance(pos2, pos5) == -1);
            REQUIRE(di.minDistance(pos5, pos6) == -1);
            REQUIRE(di.minDistance(pos1, pos7) == 16);
            REQUIRE(di.minDistance(pos1, pos9) == 8);
            REQUIRE(di.minDistance(pos1, pos8) == -1);

            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(3, false)) == -1);
            REQUIRE(di.loopDistance(snarl2,  snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(di.loopDistance(snarl5, snarl5, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl5, snarl5, make_pair(6, false),  
                                        make_pair(7, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(5, false),  
                                        make_pair(5, false)) == -1);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(8, false),  
                                        make_pair(8, false)) == -1);

/*
            for (auto x : di.maxIndex.nodeToComponent) {
                REQUIRE(x == 1);
            }
*/
        }
    }//end test case

    TEST_CASE("Top level chain", "[dist]") {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);


        SECTION("Create distance index") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);
            const Snarl* snarl8 = snarl_manager.into_which_snarl(8, false);
            const Chain* chain = snarl_manager.chain_of(snarl2);
            const Chain* topChain = snarl_manager.chain_of(snarl8);

            pair<id_t, bool> start1 = make_pair(snarl1->start().node_id(), 
                                                    snarl1->start().backward());
            pair<id_t, bool> start1r = make_pair(snarl1->start().node_id(), 
                                                   !snarl1->start().backward());
            pair<id_t, bool> end1 = make_pair(snarl1->end().node_id(), 
                                                    snarl1->end().backward());
            pair<id_t, bool> end1r = make_pair(snarl1->end().node_id(), 
                                                    !snarl1->end().backward());
            TestDistanceIndex::SnarlIndex& sd1 = 
                     di.snarlDistances.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlIndex& sd2 = 
                   di.snarlDistances.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
 
            pair<id_t, bool> start5 = make_pair(snarl5->start().node_id(), 
                                                    snarl5->start().backward());
            pair<id_t, bool> start5r = make_pair(snarl5->start().node_id(), 
                                                   !snarl5->start().backward());
            pair<id_t, bool> end5 = make_pair(snarl5->end().node_id(), 
                                                    snarl5->end().backward());
            pair<id_t, bool> end5r = make_pair(snarl5->end().node_id(), 
                                                   !snarl5->end().backward());
            TestDistanceIndex::SnarlIndex& sd5 = 
                di.snarlDistances.at(make_pair(snarl5->start().node_id(),
                                     snarl5->start().backward()));

            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 0) == 0);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 1) == 2);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 2) == 5);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 0) == -1);

            NetGraph ng1 = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);
            NetGraph ng5 = NetGraph(snarl5->start(), snarl5->end(), snarl_manager.chains_of(snarl5), &graph);
            REQUIRE(di.checkChainLoopRev(get_start_of(*chain).node_id(),1)== 3);
            REQUIRE(di.checkChainLoopRev(get_start_of(*chain).node_id(),0)==-1);

            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(1, false), make_pair(8, false))
                                                                         == 3);
            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(8, true), make_pair(1, true))
                                                                          == 4);
            REQUIRE(sd5.snarlDistance(&graph, &ng5, make_pair(6, true), make_pair(6, false))
                                                                          == 7);
            REQUIRE(di.chainDistances.size() == 2);
    
            TestDistanceIndex::ChainIndex& cd = di.chainDistances.at(get_start_of(*chain).node_id());
            REQUIRE(cd.chainDistance(start2, start2, snarl2, snarl2) == 0);
            #ifdef print
            di.printSelf();
            #endif
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(5, false), snarl2, snarl5) == 2);
            REQUIRE(cd.chainDistance(make_pair(5, true), make_pair(2, true), snarl5, snarl2) == 4);

            }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);
            const Snarl* snarl8 = snarl_manager.into_which_snarl(8, false);

            pos_t pos1 = make_pos_t(4, false, 1);
            pos_t pos2 = make_pos_t(6, false, 0);
            pos_t pos3 = make_pos_t(2, false, 0);
            pos_t pos4 = make_pos_t(7, false, 0);
            pos_t pos5 = make_pos_t(1, false, 0);
            pos_t pos6 = make_pos_t(8, true, 0);
            pos_t pos7 = make_pos_t(9, true, 0);

            REQUIRE(di.minDistance(pos1, pos2) == 6);
            REQUIRE(di.minDistance(pos2, pos1) == -1);
            REQUIRE(di.minDistance(pos3, pos1) == 2);
            REQUIRE(di.minDistance(pos1, pos3) == -1);
            REQUIRE(di.minDistance(pos3, pos2) == 5);
            REQUIRE(di.minDistance(pos3, pos4) == 5);
            REQUIRE(di.minDistance(pos5, pos2) == 8);
            REQUIRE(di.minDistance(pos5, pos6) == -1);
            REQUIRE(di.minDistance(pos3, pos7) == -1);


            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(di.loopDistance(snarl5, snarl5, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(5, false),  
                                        make_pair(5, false)) == -1);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(8, false),  
                                        make_pair(8, false)) == -1);
            REQUIRE(di.loopDistance(snarl8, snarl8, make_pair(9, false),  
                                        make_pair(9, false)) == -1);

/*
            for (auto x : di.maxIndex.nodeToComponent) {
                REQUIRE(x == 1);
            }
*/
        }
    }//end test case
    TEST_CASE("Interior chain", "[dist]") {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);


        SECTION("Create distance index") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl4 = snarl_manager.into_which_snarl(4, false);
            const Snarl* snarl7 = snarl_manager.into_which_snarl(7, false);

            const Chain* chain = snarl_manager.chain_of(snarl2);

            pair<id_t, bool> start1 = make_pair(snarl1->start().node_id(), 
                                                    snarl1->start().backward());
            pair<id_t, bool> start1r = make_pair(snarl1->start().node_id(), 
                                                   !snarl1->start().backward());
            pair<id_t, bool> end1 = make_pair(snarl1->end().node_id(), 
                                                    snarl1->end().backward());
            pair<id_t, bool> end1r = make_pair(snarl1->end().node_id(), 
                                                    !snarl1->end().backward());
            TestDistanceIndex::SnarlIndex& sd1 = 
                     di.snarlDistances.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlIndex& sd2 = 
                   di.snarlDistances.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
 
            pair<id_t, bool> start4 = make_pair(snarl4->start().node_id(), 
                                                    snarl4->start().backward());
            pair<id_t, bool> start4r = make_pair(snarl4->start().node_id(), 
                                                   !snarl4->start().backward());
            pair<id_t, bool> end4 = make_pair(snarl4->end().node_id(), 
                                                    snarl4->end().backward());
            pair<id_t, bool> end4r = make_pair(snarl4->end().node_id(), 
                                                   !snarl4->end().backward());
            TestDistanceIndex::SnarlIndex& sd4 = 
                di.snarlDistances.at(make_pair(snarl4->start().node_id(),
                                     snarl4->start().backward()));
            TestDistanceIndex::SnarlIndex& sd7 = 
                di.snarlDistances.at(make_pair(snarl7->start().node_id(),
                                     snarl7->start().backward()));
    
            TestDistanceIndex::ChainIndex& cd = 
                               di.chainDistances.at(get_start_of(*chain).node_id());

            NetGraph ng1 = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);
            #ifdef print
            di.printSelf();
            #endif

            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 0) == 0);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 1) == 1);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 2) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 3) == 10);

            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 0) == 15);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 1) == 10);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 2) == 5);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 3) == -1);

            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(1, false), make_pair(10, false))
                                                                         == 3);
            REQUIRE(di.chainDistances.size() == 1);

            REQUIRE(cd.chainDistance(start2, start2, snarl2, snarl2) == 0);
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(7, false), snarl2, snarl7) == 5);

            }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl4 = snarl_manager.into_which_snarl(4, false);
            const Snarl* snarl7 = snarl_manager.into_which_snarl(7, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos4 = make_pos_t(4, false, 1);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5r = make_pos_t(5, true, 2);
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);
            pos_t pos8r = make_pos_t(8, true, 0);
            pos_t pos9 = make_pos_t(9, false, 0);

            REQUIRE(di.minDistance(pos4, pos8) == 4);
            REQUIRE(di.minDistance(pos4, pos8) == 4);
            REQUIRE(di.minDistance(pos8, pos4) == -1);
            REQUIRE(di.minDistance(pos1, pos4) == 5);
            REQUIRE(di.minDistance(pos1, pos4) == 5);
            REQUIRE(di.minDistance(pos8, pos8r) == -1);
            REQUIRE(di.minDistance(pos7, pos8r) == 1);
            REQUIRE(di.minDistance(pos7, pos8) == 1);
            REQUIRE(di.minDistance(pos8r, pos7) == -1);
            REQUIRE(di.minDistance(pos2, pos8) == 6);
            REQUIRE(di.minDistance(pos6, pos5) == -1);
            REQUIRE(di.minDistance(pos6, pos5r) == 9);

            REQUIRE(distance(&graph, pos4, pos8) == 4);
            REQUIRE(distance(&graph, pos4, pos8) == 4);
            REQUIRE(distance(&graph, pos8, pos4) == -1);
            REQUIRE(distance(&graph, pos1, pos4) == 5);
            REQUIRE(distance(&graph, pos1, pos4) == 5);
            REQUIRE(distance(&graph, pos8, pos8r) == -1);
            REQUIRE(distance(&graph, pos7, pos8r) == 1);
            REQUIRE(distance(&graph, pos7, pos8) == 1);
            REQUIRE(distance(&graph, pos8r, pos7) == -1);
            REQUIRE(distance(&graph, pos2, pos8) == 6);
            REQUIRE(distance(&graph, pos6, pos5r) == 9);

            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(di.loopDistance(snarl4, snarl4, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl4, snarl4, make_pair(5, false),  
                                        make_pair(5, false)) == -1);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(10, false),  
                                        make_pair(10, false)) == -1);
            REQUIRE(di.loopDistance(snarl7, snarl7, make_pair(9, false),  
                                        make_pair(9, false)) == -1);

/*
            for (auto x : di.maxIndex.nodeToComponent) {
                REQUIRE(x == 1);
            }
*/
        }
    }//end test case
    TEST_CASE("Top level loop creates unary snarl", "[dist]") {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);


        // We end up with a big unary snarl of 7 rev -> 7 rev
        // Inside that we have a chain of two normal snarls 2 rev -> 3 fwd, and 3 fwd -> 6 fwd
        // And inside 2 rev -> 3 fwd, we get 1 rev -> 1 rev as another unary snarl.
        
        // We name the snarls for the distance index by their start nodes.

        SECTION("Create distance index") {
            #ifdef print
            di.printSelf();
            #endif
            auto snarl1 = snarl_manager.into_which_snarl(1, true);
            auto snarl2 = snarl_manager.into_which_snarl(2, true);
            auto snarl3 = snarl_manager.into_which_snarl(3, false);
            auto snarl7 = snarl_manager.into_which_snarl(7, true);

            const Chain* chain = snarl_manager.chain_of(snarl2);
            TestDistanceIndex::SnarlIndex& sd1 = 
                     di.snarlDistances.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            TestDistanceIndex::SnarlIndex& sd2 = 
                     di.snarlDistances.at(make_pair(snarl2->start().node_id(),
                                              snarl2->start().backward()));
            TestDistanceIndex::SnarlIndex& sd7 = 
                     di.snarlDistances.at(make_pair(snarl7->start().node_id(),

                                              snarl7->start().backward()));  


            NetGraph ng1 = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);
            NetGraph ng2 = NetGraph(snarl2->start(), snarl2->end(), snarl_manager.chains_of(snarl2), &graph);
            NetGraph ng7 = NetGraph(snarl7->start(), snarl7->end(), snarl_manager.chains_of(snarl7), &graph);
            TestDistanceIndex::ChainIndex& cd = di.chainDistances.at(get_start_of(*chain).node_id());

            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(1, true), make_pair(1, false))
                                                                         == 3);
            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(1, false), make_pair(1, true))
                                                                         == -1);
            REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(2, true), make_pair(1, true))
                                                                         == 1);
            REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(2, true), make_pair(1, false))
                                                                         == -1);
            REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(1, false), make_pair(3, false))
                                                                         == 6);
            REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(3, true), make_pair(1, true))
                                                                         == 1);
            REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(3, true), make_pair(1, false))
                                                                         == -1);
            REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(3, true), make_pair(3, false))
                                                                         == 7);
             
            REQUIRE(cd.chainDistance(make_pair(6, true), make_pair(3, true), snarl3, snarl3) 
                                                                    == 4);
            REQUIRE(cd.chainDistance(make_pair(3, false), make_pair(6, false), snarl3, snarl3) 
                                                                    == 4);
            REQUIRE(cd.chainDistance(make_pair(6, true), make_pair(3, false), snarl3, snarl3) 
                                                                    == 11);
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(2, true), snarl2, snarl2) 
                                                                    == -1);
            REQUIRE(cd.chainDistance(make_pair(2, true), make_pair(2, false), snarl2, snarl2) 
                                                                    == 7);
            REQUIRE(cd.chainDistance(make_pair(6, true), make_pair(2, true), snarl3, snarl2) 
                                                                    == -1);
            REQUIRE(cd.chainDistance(make_pair(6, true), make_pair(2, false), snarl3, snarl2) 
                                                                    == 11);
            REQUIRE(di.chainDistances.size() == 1);

            REQUIRE(sd7.snarlDistance(&graph, &ng7, make_pair(7, true), make_pair(2, false))
                                                                         == 1);
            REQUIRE(sd7.snarlDistance(&graph, &ng7, make_pair(7, true), make_pair(2, true))
                                                                         == 1);
            REQUIRE(sd7.snarlDistance(&graph, &ng7, make_pair(7, true), make_pair(7, false))
                                                                       == 9);
            REQUIRE(sd7.snarlDistance(&graph, &ng7, make_pair(2, true), make_pair(7, false))
                                                                        == 12);
            REQUIRE(sd7.snarlDistance(&graph, &ng7, make_pair(2, false), make_pair(7, true))
                                                                        == -1);

            }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, true);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, true);
            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, true);
            const Snarl* snarl7 = snarl_manager.into_which_snarl(7, true);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos1r = make_pos_t(1, true, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5r = make_pos_t(5, true, 2);
            pos_t pos7 = make_pos_t(7, false, 0);

            REQUIRE(di.minDistance(pos1, pos7) == 4);
            REQUIRE(di.minDistance(pos2, pos3) == -1);
            REQUIRE(di.minDistance(pos4, pos7) == 5);
            REQUIRE(di.minDistance(pos1, pos5) == 4);
            REQUIRE(di.minDistance(pos5, pos1) == -1);
            REQUIRE(di.minDistance(pos5r, pos1) == 5);
            REQUIRE(di.minDistance(pos5r, pos1r) == 2);

            REQUIRE(distance(&graph, pos1, pos7) == 4);
            REQUIRE(distance(&graph, pos2, pos3) == -1);
            REQUIRE(distance(&graph, pos4, pos7) == 5);
            REQUIRE(distance(&graph, pos1, pos5) == 4);

            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(di.loopDistance(snarl6, snarl6, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl6, snarl6, make_pair(5, false),  
                                        make_pair(5, false)) == -1);
            REQUIRE(di.loopDistance(snarl7, snarl7, make_pair(7, false),  
                                        make_pair(7, false)) == -1);

/*
            for (auto x : di.maxIndex.nodeToComponent) {
                REQUIRE(x == 1);
            }
*/

        }
    }//end test case
    TEST_CASE( "Shortest path exits common ancestor","[dist]" ) {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);
        
        SECTION( "Create distance index" ) {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlIndex& sd1 = di.snarlDistances.at(make_pair(
                  snarl1->start().node_id(),snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlIndex& sd2 = di.snarlDistances.at(make_pair(
                       snarl2->start().node_id(),snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            TestDistanceIndex::SnarlIndex& sd3 = di.snarlDistances.at(make_pair(
                       snarl3->start().node_id(),snarl3->start().backward()));

            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);
            TestDistanceIndex::SnarlIndex& sd6 = di.snarlDistances.at(make_pair(
                       snarl6->start().node_id(),snarl6->start().backward()));

#ifdef print
            di.printSelf();
#endif
            NetGraph ng = NetGraph(snarl3->start(), snarl3->end(), snarl_manager.chains_of(snarl3), &graph);

            REQUIRE(sd3.snarlDistance(&graph, &ng, make_pair(5, false), make_pair(5, false)) == 12); 
        }


        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos4r = make_pos_t(4, true, 3);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5e = make_pos_t(5, false, 11);
            pos_t pos5r = make_pos_t(5, true, 11);
            pos_t pos7 = make_pos_t(7, false, 0);

            REQUIRE(di.minDistance(pos5e, pos5) == 1);
            REQUIRE(di.minDistance(pos5, pos5e) == 11);
            REQUIRE(distance(&graph, pos5, pos5e) == 11);
            REQUIRE(di.minDistance(pos4, pos5) == -1);
            REQUIRE(di.minDistance(pos4r, pos5) == 5);
            REQUIRE(di.minDistance(pos5, pos5) == 0);
            REQUIRE(di.minDistance(pos5r, pos5) == 5);
            REQUIRE(di.minDistance(pos5, pos5r) == -1);
            REQUIRE(di.minDistance(pos5, pos7) == 13);
            REQUIRE(di.minDistance(pos5r, pos7) == 10);


            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(di.loopDistance(snarl6, snarl6, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(3, false),  
                                        make_pair(3, false)) == -1);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(5, false),  
                                        make_pair(5, false)) == 12);
            REQUIRE(di.loopDistance(snarl6, snarl6, make_pair(7, false),  
                                        make_pair(7, false)) == -1);

//            REQUIRE(di.maxIndex.nodeToComponent[4] == 1);
        }
    }//End test case
    TEST_CASE( "Simple nested snarl with loop",
                   "[dist]" ) {
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
        Edge* e11 = graph.create_edge(n7, n7, false, true);
        Edge* e12 = graph.create_edge(n3, n3, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        TestDistanceIndex di (&graph, &snarl_manager, 20);
        
        SECTION( "Create distance index" ) {

#ifdef print
            di.printSelf();
#endif
            const vector<const Snarl*> topSnarls = snarl_manager.top_level_snarls();
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlIndex& sd1 = di.snarlDistances.at(
                                    make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlIndex& sd2 = di.snarlDistances.at(
                            make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            TestDistanceIndex::SnarlIndex& sd3 = di.snarlDistances.at(
                                       make_pair(snarl3->start().node_id(),
                                                 snarl3->start().backward()));

            NetGraph ng = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);

            REQUIRE(sd1.snarlDistance(&graph, &ng, make_pair(1, false), make_pair(2, false)) == 3); 
            REQUIRE(sd1.snarlDistance(&graph, &ng, make_pair(1, false), make_pair(2, true)) == 6); 
            REQUIRE(sd1.snarlDistanceShort(make_pair(1, false), make_pair(2, false)) == 0); 
            REQUIRE(sd1.snarlDistance(&graph, &ng, make_pair(2, true), make_pair(1,true)) == 3); 

        }


        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

            pos_t pos1 = make_pos_t(4, false, 1);
            pos_t pos2 = make_pos_t(5, false, 2);
            pos_t pos3 = make_pos_t(5, false, 1);
            pos_t pos4 = make_pos_t(4, false, 2);
            pos_t pos5 = make_pos_t(2, false, 0);
            pos_t pos6 = make_pos_t(2, true, 0);
            pos_t pos6r = make_pos_t(6, true, 0);
            pos_t pos7 = make_pos_t(1, true, 1);
            pos_t pos8 = make_pos_t(8, true, 1);
            pos_t pos9 = make_pos_t(6, false, 0);
            pos_t pos10 = make_pos_t(4, false, 0);
            pos_t pos4r = make_pos_t(4, true, 2);

            REQUIRE(di.minDistance(pos3, pos2) == 1);
            REQUIRE(di.minDistance(pos1,pos2) == 5);

            REQUIRE(di.minDistance(pos5,pos2) == 4);
            REQUIRE(di.minDistance(pos2, pos5) == -1);
            REQUIRE(di.minDistance(pos6,pos2) == -1);
            REQUIRE(di.minDistance(pos8, pos7) == 4);
            REQUIRE(di.minDistance(pos4, pos6r) == 7);
            REQUIRE(di.minDistance(pos4r, pos9) == -1);
            REQUIRE(di.minDistance(pos9, pos4r) == 8);


            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(6, false),  
                                        make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(3, false),  
                                        make_pair(3, false)) == 10);
            REQUIRE(loopDistance(&graph, snarl3, snarl3, make_pair(3, false),  
                                        make_pair(3, false)) == 10);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(5, false),  
                                        make_pair(5, false)) == 10);
            REQUIRE(loopDistance(&graph, snarl3, snarl3, make_pair(5, false),  
                                        make_pair(5, false)) == 10);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(4, false),  
                                        make_pair(4, false)) == 14);
            REQUIRE(loopDistance(&graph, snarl3, snarl3, make_pair(4, false),  
                                        make_pair(4, false)) == 14);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(7, false),  
                                        make_pair(7, false)) == 10);

//            REQUIRE(di.maxIndex.nodeToComponent[2] == 1);
//            REQUIRE(di.maxIndex.nodeToComponent[3] == 1);
//            REQUIRE(di.maxIndex.nodeToComponent[4] == 1);
//            REQUIRE(di.maxIndex.nodeToComponent[6] == 1);
        }
    }//End test case

    TEST_CASE( "Three loops",
                   "[dist]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");
        Node* n9 = graph.create_node("G");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("G");
        Node* n12 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n12);
        Edge* e3 = graph.create_edge(n1, n10);
        Edge* e4 = graph.create_edge(n2, n3);
        Edge* e5 = graph.create_edge(n2, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n5, n2);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n6, n8);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n7, n7);
        Edge* e14 = graph.create_edge(n8, n9);
        Edge* e15 = graph.create_edge(n9, n8);
        Edge* e16 = graph.create_edge(n9, n11);
        Edge* e17 = graph.create_edge(n11, n10);
        Edge* e18 = graph.create_edge(n12, n2);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        TestDistanceIndex di (&graph, &snarl_manager, 20);
        
        SECTION( "Create distance index" ) {
#ifdef print
            di.printSelf();
#endif

            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            REQUIRE(di.loopDistance(make_pair(7, false),  
                                                    make_pair(8, false)) == -1);
            
            REQUIRE(di.loopDistance(make_pair(9, false),  
                                                   make_pair(11, false)) == -1);
            REQUIRE(di.loopDistance(make_pair(11, false),  
                                                   make_pair(10, false)) == -1);
            REQUIRE(loopDistance(&graph, snarl1, snarl1, make_pair(11, false),  
                                                   make_pair(10, false)) == -1);
            REQUIRE(di.loopDistance(make_pair(12, false),  
                                                   make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(make_pair(1, false),  
                                                   make_pair(12, false)) == -1);

            int64_t c1 = di.maxIndex.nodeToComponent[1];
            int64_t c2 = di.maxIndex.nodeToComponent[6];
            int64_t c3 = di.maxIndex.nodeToComponent[7];
            int64_t c4 = di.maxIndex.nodeToComponent[0];
            int64_t c5 = di.maxIndex.nodeToComponent[5];
            REQUIRE(di.maxIndex.nodeToComponent[0] == c4);
            REQUIRE(di.maxIndex.nodeToComponent[1] == c1);
            REQUIRE(di.maxIndex.nodeToComponent[2] == c1);
            REQUIRE(di.maxIndex.nodeToComponent[3] == c1);
            REQUIRE(di.maxIndex.nodeToComponent[4] == c1);
            REQUIRE(di.maxIndex.nodeToComponent[5] == c5);
            REQUIRE(di.maxIndex.nodeToComponent[6] == c2);
            REQUIRE(di.maxIndex.nodeToComponent[7] == c3);
            REQUIRE(di.maxIndex.nodeToComponent[8] == c3);
            REQUIRE(di.maxIndex.nodeToComponent[9] == c4);
            REQUIRE(di.maxIndex.nodeToComponent[10] == c4);
            REQUIRE(di.maxIndex.nodeToComponent[11] == c4);
            REQUIRE(c4 != c5);

        }
    }//End test case

    TEST_CASE( "More loops",
                   "[dist]" ) {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);
        
        SECTION( "Create distance index" ) {

#ifdef print
            di.printSelf();
#endif
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlIndex& sd1 = di.snarlDistances.at(
                                    make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlIndex& sd2 = di.snarlDistances.at(
                            make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
            NetGraph ng1 = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);

            NetGraph ng2 = NetGraph(snarl2->start(), snarl2->end(), snarl_manager.chains_of(snarl2), &graph);

            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(1, false), make_pair(7, false)) == 3);
            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(1, false), make_pair(7, false)) == 3);
            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(
                       snarl2->start().node_id(), !snarl2->start().backward()), 
                      make_pair(snarl2->start().node_id(), 
                                           snarl2->start().backward())) == -1);
            REQUIRE(sd2.snarlDistance(&graph, &ng2, make_pair(4, true), make_pair(4, false)) == 6);
            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(5, true), make_pair(5, false)) == 13);
            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(5, true), make_pair(7, false)) == 14);
            REQUIRE(sd1.snarlDistance(&graph, &ng1, make_pair(7, true), make_pair(7, false)) == 13);

        }


        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5r = make_pos_t(5, true, 2);
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos6r = make_pos_t(6, true, 0);
            pos_t pos7 = make_pos_t(7, false, 0);

            REQUIRE(di.minDistance(pos1, pos4) == 4);
            REQUIRE(di.minDistance( pos5, pos6) == -1);
            REQUIRE(di.minDistance( pos5, pos6r) == -1);
            REQUIRE(di.minDistance( pos5r, pos6) == 11);
            REQUIRE(di.minDistance(pos2, pos7) == 6);

            REQUIRE(distance(&graph, pos1, pos4) == 4);
            REQUIRE(distance(&graph, pos5r, pos6) == 11);
            REQUIRE(distance(&graph, pos2, pos7) == 6);

            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                                    make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                                    make_pair(2, true)) == -1);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(6, false),  
                                                    make_pair(6, false)) == -1);
            REQUIRE(di.loopDistance( make_pair(3, false),  
                                                    make_pair(3, false)) == -1);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(5, false),  
                                                    make_pair(5, false)) == -1);

      
//            for (auto x : di.maxIndex.nodeToComponent) {
//                REQUIRE(x == 1);
//            }
        }
    }//End test case

    TEST_CASE( "Exit common ancestor","[dist]" ) {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);
        
        SECTION( "Create distance index" ) {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlIndex& sd1 = di.snarlDistances.at(make_pair(
                  snarl1->start().node_id(),snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlIndex& sd2 = di.snarlDistances.at(make_pair(
                       snarl2->start().node_id(),snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            TestDistanceIndex::SnarlIndex& sd3 = di.snarlDistances.at(make_pair(
                       snarl3->start().node_id(),snarl3->start().backward()));

            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);
            TestDistanceIndex::SnarlIndex& sd6 = di.snarlDistances.at(make_pair(
                       snarl6->start().node_id(),snarl6->start().backward()));

#ifdef print
            di.printSelf();
#endif
        }


        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2r = make_pos_t(2, true, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos3r = make_pos_t(3, true, 0);
            pos_t pos4r = make_pos_t(4, true, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5r = make_pos_t(5, false, 11);
            pos_t pos7r = make_pos_t(7, true, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);
            pos_t pos9 = make_pos_t(9, false, 0);
            pos_t pos9r = make_pos_t(9, true, 1);
            pos_t pos11 = make_pos_t(11, false, 0);
            pos_t pos12 = make_pos_t(12, false, 0);
            pos_t pos12r = make_pos_t(12, true, 0);

            REQUIRE(di.minDistance( pos2r, pos9r) == 2);
            REQUIRE(di.minDistance( pos2, pos9r) == 5);
            REQUIRE(di.minDistance( pos3r, pos9) == 4);
            REQUIRE(di.minDistance( pos3r, pos9r) == 3);
            REQUIRE(di.minDistance( pos3, pos9) == 11);
            REQUIRE(di.minDistance( pos4, pos5) == 14);
            REQUIRE(di.minDistance( pos4r, pos5) == 8);
            REQUIRE(di.minDistance(pos7, pos12) == 14);
            REQUIRE(di.minDistance(pos7, pos12r) == 13);
            REQUIRE(di.minDistance(pos7r, pos12) == 15);
            REQUIRE(di.minDistance(pos8, pos11) == 7);



            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(1, false),  
                                        make_pair(1, false)) == -1);
            REQUIRE(di.loopDistance(snarl1, snarl2, make_pair(1, false),  
                                        make_pair(2, false)) == -1);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == 4);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == 4);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(11, false)) == 4);
            REQUIRE(di.loopDistance(snarl2, snarl3, make_pair(2, false),  
                                        make_pair(3, false)) == 14);
            REQUIRE(di.loopDistance(snarl3, snarl2, make_pair(3, true),  
                                        make_pair(2, true)) == 14);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(5, false),  
                                        make_pair(5, false)) == 22);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(3, false),  
                                        make_pair(3, false)) == 14);
            REQUIRE(di.loopDistance(snarl6, snarl2, make_pair(8, false),  
                                        make_pair(9, false)) == 14);
            REQUIRE(di.loopDistance(snarl2, snarl6, make_pair(9, true),  
                                        make_pair(8, true)) == 14);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(11, true),  
                                        make_pair(2, true)) == 4);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(11, false)) == 4);
            REQUIRE(di.loopDistance(snarl6, snarl6, make_pair(7, true),  
                                        make_pair(6, true)) == 14);

        }
    }//End test case

    TEST_CASE( "Exit common ancestor chain","[dist]" ) {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);
        
        SECTION( "Create distance index" ) {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlIndex& sd1 = di.snarlDistances.at(make_pair(
                  snarl1->start().node_id(),snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlIndex& sd2 = di.snarlDistances.at(make_pair(
                       snarl2->start().node_id(),snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            TestDistanceIndex::SnarlIndex& sd3 = di.snarlDistances.at(make_pair(
                       snarl3->start().node_id(),snarl3->start().backward()));

            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);
            TestDistanceIndex::SnarlIndex& sd6 = di.snarlDistances.at(make_pair(
                       snarl6->start().node_id(),snarl6->start().backward()));

#ifdef print
            di.printSelf();
#endif
        }


        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);
            const Snarl* snarl9 = snarl_manager.into_which_snarl(9, false);
            const Snarl* snarl13 = snarl_manager.into_which_snarl(13, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2r = make_pos_t(2, true, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos4r = make_pos_t(4, true, 3);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5r = make_pos_t(5, true, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);
            pos_t pos9 = make_pos_t(9, false, 0);
            pos_t pos9r = make_pos_t(9, false, 1);
            pos_t pos10 = make_pos_t(10, false, 0);
            pos_t pos10r = make_pos_t(10, true, 0);
            pos_t pos11 = make_pos_t(11, false, 0);
            pos_t pos12 = make_pos_t(12, false, 0);
            pos_t pos14 = make_pos_t(14, false, 0);
            pos_t pos16 = make_pos_t(16, false, 0);

            REQUIRE(di.minDistance(pos2r, pos10r) == 2);
            REQUIRE(di.minDistance(pos2, pos10) == 4);
            REQUIRE(di.minDistance(pos4r, pos5) == 5);
            REQUIRE(di.minDistance(pos4, pos5) == 11);
            REQUIRE(di.minDistance(pos4r, pos5r) == 8);
            REQUIRE(di.minDistance(pos14, pos10) == 10);
            REQUIRE(di.minDistance(pos14, pos10r) == 5);
            REQUIRE(di.minDistance(pos14, pos3) == 7);
            REQUIRE(di.minDistance(pos16, pos3) == 5);

            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == 5);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == 5);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(5, false),  
                                        make_pair(5, false)) == 19);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(3, false),  
                                        make_pair(3, false)) == 11);
            REQUIRE(di.loopDistance(snarl3, snarl6, make_pair(5, false),  
                                        make_pair(6, false)) == 19);
            REQUIRE(di.loopDistance(snarl1, snarl2, make_pair(1, false),  
                                        make_pair(2, false)) == 9);
            REQUIRE(di.loopDistance(make_pair(1, false),  
                                        make_pair(2, false)) == 9);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(16, false)) == 5);
            REQUIRE(di.loopDistance(snarl2, snarl3, make_pair(2, false),  
                                        make_pair(3, false)) == 11);
            REQUIRE(di.loopDistance(snarl1, snarl2, make_pair(1, false),  
                                        make_pair(2, false)) == 9);
            REQUIRE(di.loopDistance(snarl9, snarl1, make_pair(11, false),  
                                        make_pair(12, false)) == 9);
            REQUIRE(di.loopDistance(snarl1, snarl9, make_pair(12, true),  
                                        make_pair(11, true)) == 9);
            REQUIRE(di.loopDistance(snarl13, snarl13, make_pair(13, false),  
                                        make_pair(15, false)) == 7);
            REQUIRE(di.loopDistance(snarl6, snarl6, make_pair(6, false),  
                                        make_pair(8, false)) == 11);
            REQUIRE(di.loopDistance(snarl6, snarl6, make_pair(6, false),  
                                        make_pair(7, false)) == 12);
            REQUIRE(di.loopDistance(snarl9, snarl2, make_pair(11, false),  
                                        make_pair(2, false)) == 5);

//            for (auto x : di.maxIndex.nodeToComponent) {
//                REQUIRE(x == 1);
//            }
        }
    }//End test case

    TEST_CASE("Top level loops", "[dist]") {
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
        TestDistanceIndex di (&graph, &snarl_manager, 20);


        SECTION("Create distance index") {

            #ifdef print
                di.printSelf();
            #endif
 
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Chain* chain = snarl_manager.chain_of(snarl1);

            const Snarl* snarl8 = snarl_manager.into_which_snarl(8, false);
            TestDistanceIndex::ChainIndex& cd = di.chainDistances.at(get_start_of(*chain).node_id());
            
                // Chain should be 1 fwd -> 8 fwd, 8 fwd -> 1 fwd
                // or 8 fwd -> 1 fwd, 1 fwd -> 8 fwd
                // after snarls have been flipped into their in-chain orientations.
                // Node side flags are normal.
                
                // Distance from start of node 1 to start of node 8 should be 1
                REQUIRE(cd.chainDistance(make_pair(1, false), make_pair(8, false), snarl1, snarl8) == 1);
                // Distance from start of node 8 to start of node 1 should be 3
                REQUIRE(cd.chainDistance(make_pair(8, false), make_pair(1, false), snarl8, snarl1) == 3);
                // Distance from end of node 1 to end of node 8, reading backward, should be 3
                REQUIRE(cd.chainDistance(make_pair(1, true), make_pair(8, true), snarl1, snarl8) == 3);
                // Distance from end of node 8 to end of node 1, reading backward, should be 1
                REQUIRE(cd.chainDistance(make_pair(8, true), make_pair(1, true),
                                                       snarl8, snarl1) == 1);
                
                // Distance from side opposite start of node 8 to start of node 1 should be 2
                REQUIRE(cd.chainDistanceShort(&graph, make_pair(8, false), make_pair(1, false), snarl8, snarl1) == 2);
                
          
        }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);
            const Snarl* snarl9 = snarl_manager.into_which_snarl(9, false);

            pos_t pos1 = make_pos_t(1, false, 1);
            pos_t pos2r = make_pos_t(2, true, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos6r = make_pos_t(6, true, 0);
            pos_t pos7r = make_pos_t(7, true, 0);
            pos_t pos7 = make_pos_t(7, false, 0);

            REQUIRE(di.minDistance(pos4, pos6) == 19);
            REQUIRE(di.minDistance(pos4, pos6r) == 25);
            REQUIRE(di.minDistance(pos2r, pos7) == 7);
            REQUIRE(di.minDistance(pos2, pos7) == 6);
            REQUIRE(di.minDistance(pos2, pos7r) == 11);

            REQUIRE(distance(&graph, pos4, pos6r) == 25);
            REQUIRE(distance(&graph, pos2, pos7r) == 11);

            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == 11);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == 11);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(5, false),  
                                        make_pair(5, false)) == 11);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(3, false),  
                                        make_pair(3, false)) == 11);
            REQUIRE(di.loopDistance(snarl5, snarl5, make_pair(6, false),
                                        make_pair(6, false)) == 12);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(8, false),
                                        make_pair(8, false)) == 4);

            REQUIRE(di.loopDistance(snarl5, snarl1, make_pair(7, false),  
                                        make_pair(8, false)) == 11);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(3, false)) == 11);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(1, false),  
                                        make_pair(8, false)) == 4);

//            for (auto x : di.maxIndex.nodeToComponent) {
//                REQUIRE(x == 1);
//            }
        }
    }//end test case


    TEST_CASE( "Snarl that loops on itself",
                  "[dist]" ) {
        
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
            TestDistanceIndex di (&graph, &snarl_manager, 20);
            #ifdef print
                di.printSelf();
            #endif
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, true);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);

            REQUIRE(di.minDistance(pos1, pos2) == distance(&graph, pos1, pos2));



            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, false),  
                                        make_pair(2, false)) == 10);
            REQUIRE(di.loopDistance(snarl2, snarl2, make_pair(2, true),  
                                        make_pair(2, true)) == 10);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(4, false),  
                                        make_pair(4, false)) == 14);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(5, false),  
                                        make_pair(5, false)) == 10);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(3, false),  
                                        make_pair(3, false)) == 10);
            REQUIRE(di.loopDistance(snarl1, snarl1, make_pair(1, false),
                                        make_pair(1, false)) == 10);
            REQUIRE(di.loopDistance(snarl1, snarl2, make_pair(1, false),
                                        make_pair(2, false)) == 10);
            REQUIRE(di.loopDistance(snarl2, snarl3, make_pair(2, false),
                                        make_pair(3, false)) == 10);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(3, false),
                                        make_pair(5, true)) == 10);
            REQUIRE(di.loopDistance(snarl3, snarl3, make_pair(3, false),
                                        make_pair(4, false)) == 14);

/*
            for (auto x : di.maxIndex.nodeToComponent) {
                REQUIRE(x == 1);
            }
*/
        }
 
    }

    TEST_CASE("Random test", "[dist][rand]") {

/*
        ifstream vg_stream("testGraph");
        VG vg(vg_stream);
        vg_stream.close();
        CactusSnarlFinder bubble_finder(vg);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 

        TestDistanceIndex di (&vg, &snarl_manager, 50);
        pos_t pos1 = make_pos_t(98, false, 4);
        pos_t pos2 = make_pos_t(98, true, 6); 
        REQUIRE(di.maxDistance(pos1, pos2) >= 41);

            for (size_t i = 0 ; i < vg.max_node_id(); i++) {
                if (vg.has_node(i+1)) {
                    REQUIRE(di.maxIndex.nodeToComponent[i] > 0);
                    REQUIRE(di.maxIndex.maxDistances[i] >= 
                            di.maxIndex.minDistances[i]);
                    REQUIRE(di.maxIndex.minDistances[i] >= 0);
                    REQUIRE(di.maxIndex.maxDistances[i] >= 0);
                    if (di.maxIndex.nodeToComponent[i] > 
                                di.maxIndex.numCycles) {
                        REQUIRE(di.maxIndex.maxDistances[i] > 0);
                        REQUIRE(di.maxIndex.minDistances[i] > 0);
                    }
                }
            }

*/

        for (int i = 0; i < 0; i++) {
            //1000 different graphs
            VG graph = randomGraph(1000, 20, 100); 


            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls(); 

            uint64_t cap = 50;
            TestDistanceIndex di (&graph, &snarl_manager, cap);
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
                           snarl_manager.shallow_contents(snarl1, graph, true);
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
                    int64_t actDist = distance(&graph, pos1, pos2);

                    int64_t maxDist = di.maxDistance(pos1, pos2);

         
                    int64_t myLoop = di.loopDistance(snarl1, snarl1, make_pair(nodeID1, false), make_pair(nodeID1, false));
                    int64_t actLoop = loopDistance(&graph, snarl1, snarl1, make_pair(nodeID1, false), make_pair(nodeID1, false));

 
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
                    int64_t myEdge;
                    int64_t actEdge; 
                    if (found) {
                        myEdge = di.loopDistance(make_pair(nodeID1, false),
                                                     next); 
                    
                        actEdge = loopDistance(&graph, snarl1, di.snarlOf(nodeID1), make_pair(nodeID1, false), next);
           
                    } else {myEdge = 0; actEdge = 0;}
 
                    if (snarl_manager.into_which_snarl(nodeID1, false) == NULL ||
                        snarl_manager.into_which_snarl(nodeID1, true) == NULL) {
                        REQUIRE( di.snarlOf(nodeID1) == snarl1); 
                    }
                    if (snarl_manager.into_which_snarl(nodeID2, false) == NULL ||
                        snarl_manager.into_which_snarl(nodeID2, true) == NULL) {
                        REQUIRE( di.snarlOf(nodeID2) == snarl2);
                    }
                    bool passed = (myDist == actDist) && (myLoop == actLoop) &&
                                   (myEdge == actEdge);

                    if (myDist < cap) {
                        passed = passed && maxDist > 0;
                        passed = passed &&  (myDist <= maxDist);
                    }	    


                    if (!passed) { 
                        graph.serialize_to_file("testGraph");
                        cerr << "Failed on random test: " << endl;
                        
                        cerr << "Position 1 on snarl " << 
                               snarl1->start().node_id() << " " << " Node " <<
                               nodeID1 << " is rev? " << is_rev(pos1) << 
                               " offset: " << offset1 << endl;
                        cerr << "Position 2 on snarl " << 
                              snarl2->start().node_id() << " Node " << nodeID2 
                              << " is rev? " << is_rev(pos2) << " offset: " 
                              << offset2 << endl;

                        cerr << "Actual distance: " << actDist << "    " <<
                                "Guessed distance: " << myDist << endl;
                        cerr << "Guessed max distance: " << maxDist << endl;
                        cerr << endl;

                        cerr << "Actual loop distance: "  << actLoop << "    " 
                             << "Guessed loop distance: " << myLoop << endl;
                        cerr << "Actual edge loop distance: "  << actEdge
                             << "    " << "Guessed edge loop distance: " << 
                             myEdge << " From " << nodeID1 << " TO " 
                             << next.first << endl;
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

/*
    TEST_CASE("Serialize distance index", "[dist][serial]") {
        for (int i = 0; i < 100; i++) {

            VG graph = randomGraph(1000, 20, 100); 

            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls(); 

            TestDistanceIndex di (&graph, &snarl_manager, 50);
  
            filebuf buf;
            ofstream out("distanceIndex");
           
            di.serialize(out);
            out.close();

            //buf.close();

            buf.open("distanceIndex", ios::in);
            istream in(&buf);
            TestDistanceIndex sdi (&graph, &snarl_manager, in);
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
                           snarl_manager.shallow_contents(snarl1, graph, true);
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

                    int64_t maxDist = di.maxDistance(pos1, pos2);
                    int64_t serialMaxDist = sdi.maxDistance(pos1, pos2);
                    REQUIRE(maxDist == serialMaxDist);
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
