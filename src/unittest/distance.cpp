#include <stdio.h>
#include <stdlib.h>
#include <iostream>
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

int64_t distanceHelp(VG* graph, pos_t pos1, pos_t pos2,
                       bool rev){
    //Distance using djikstras algorithm

    auto cmp = [] (pair<pair<id_t, bool> , int64_t> x, 
                   pair<pair<id_t, bool>, int64_t> y ) {
        return (x.second > y.second);
    };
 
    int64_t shortestDistance = -1;
    if (get_id(pos1) == get_id(pos2)) { //if positions are on the same node
        int64_t nodeSize = graph->get_node(get_id(pos1))->sequence().size();
        int64_t offset1;
        if (is_rev(pos1)) {
            offset1 = nodeSize -get_offset(pos1) - 1;//Len of node - offset 
        } else {
            offset1 = get_offset(pos1);
        }

        int64_t offset2;
        if (is_rev(pos2)) {
            offset2 = nodeSize - get_offset(pos2) - 1;
        } else {
            offset2 = get_offset(pos2);
        }

        if (graph->has_edge(node_start(get_id(pos1)), node_end(get_id(pos2)))){
            //If there is an edge from start to end of node

            shortestDistance = min(   abs(offset1-offset2)+1,
                          nodeSize - abs(offset1-offset2) + 1  );

        } else {

            shortestDistance = abs(offset1-offset2)+1; //+1 to be consistent

        }
    }


    priority_queue< pair<pair<id_t, bool> , int64_t>, 
                    vector<pair<pair<id_t, bool>, int64_t>>,
                          decltype(cmp)> reachable(cmp); 
    handle_t currHandle = graph->get_handle(get_id(pos1), rev);

    int64_t dist;
    if (is_rev(pos1) != rev) { 
        dist = get_offset(pos1) + 1;
    } else {
        dist = graph->get_length(currHandle) - get_offset(pos1);
    }

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

        if (currID.first == get_id(pos2)){
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
    return shortestDistance;
};

int64_t distance(VG* graph, pos_t pos1, pos_t pos2){
    //Find the distance between two positions
    int64_t d1 = distanceHelp(graph, pos1, pos2, true);
    int64_t d2 = distanceHelp(graph, pos1, pos2, false);
   
    if (d1 == -1) {return d2;}
    else if (d2 == -1) {return d1;}
    else {return min(d1, d2);}
};
class TestDistanceIndex : public DistanceIndex {

    public:
        using DistanceIndex::DistanceIndex;
        using DistanceIndex::SnarlDistances;
        using DistanceIndex::snarlIndex;
        using DistanceIndex::ChainDistances;
        using DistanceIndex::chainIndex;
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
        TestDistanceIndex di (&graph, &snarl_manager);
        
        SECTION( "Create distance index" ) {
            const vector<const Snarl*> topSnarls = snarl_manager.top_level_snarls();
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            pair<id_t, bool> start1 = make_pair(snarl1->start().node_id(), 
                                                    snarl1->start().backward());
            pair<id_t, bool> end1 = make_pair(snarl1->end().node_id(), 
                                                    snarl1->end().backward());
            pair<id_t, bool> end1r = make_pair(snarl1->end().node_id(), 
                                                    !snarl1->end().backward());
            TestDistanceIndex::SnarlDistances& sd1 = di.snarlIndex.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlDistances& sd2 = di.snarlIndex.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            pair<id_t, bool> start3 = make_pair(snarl3->start().node_id(), 
                                                    snarl3->start().backward());
            pair<id_t, bool> start3r = make_pair(snarl3->start().node_id(), 
                                                   !snarl3->start().backward());
            pair<id_t, bool> end3 = make_pair(snarl3->end().node_id(), 
                                                    snarl3->end().backward());
            pair<id_t, bool> end3r = make_pair(snarl3->end().node_id(), 
                                                    !snarl3->end().backward());
            TestDistanceIndex::SnarlDistances& sd3 = di.snarlIndex.at(make_pair(snarl3->start().node_id(),
                                                 snarl3->start().backward()));

#ifdef print
            di.printSelf();
#endif
            NetGraph ng = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);

            REQUIRE(sd1.snarlDistance(start1, make_pair(2, false)) == 3); 
            REQUIRE(sd1.snarlDistance(start1, make_pair(2, true)) == -1); 
            REQUIRE(sd1.snarlDistanceShort(&graph , &ng, start1, make_pair(2, false)) == 0); 
            REQUIRE(sd1.snarlDistance(start2r, make_pair(1,true)) == 3); 
            REQUIRE(sd1.snarlDistanceShort(&graph , &ng, make_pair(8, true), make_pair(1, true)) == 0); 
            REQUIRE(sd1.snarlDistanceShort(&graph , &ng, make_pair(1, false), make_pair(8, false)) == 0); 
 

             //REQUIRE(sd2.snarlDistance(make_pair(3, false), make_pair(7, false)) == 4);          
             REQUIRE(sd2.snarlDistance(make_pair(7, true), make_pair(2, true)) == 2);

             REQUIRE(((sd2.snarlDistance(end2r, start3r) == 1 && //3 is start
                      sd2.snarlDistance(start2, start3) == 1) || 
                     (sd2.snarlDistance(end2r, start3) == 1 &&   //5 is start
                      sd2.snarlDistance(start2, start3r) == 1)
              ));

             REQUIRE(((sd1.snarlDistance(end1r, start2r) == 4 && //3 is start
                      sd1.snarlDistance(start1, start2) == 3) || 
                     (sd1.snarlDistance(end1r, start2) == 4 &&   //5 is start
                      sd1.snarlDistance(start1, start2r) == 3)
              ));

             REQUIRE(sd3.snarlDistance(make_pair(5, true), make_pair(4, true))
                                                                          == 3);
             REQUIRE((sd2.snarlDistance(make_pair(2, false), start3) == 1 ||
                       sd2.snarlDistance(make_pair(2, false), start3r) == 1 ) );
             REQUIRE(( sd2.snarlDistance(start2, start3) == 1 || 
                       sd2.snarlDistance(start2, start3r) == 1)   );

             REQUIRE(( sd2.snarlDistance(end2r, start3) == 1 || 
                       sd2.snarlDistance(end2r, start3r) == 1)   );
            REQUIRE(di.chainIndex.size() == 0);
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

            pair<int64_t, int64_t> distances = di.distToCommonAncestor(snarl3, snarl1, pos4).first;
            REQUIRE(((distances.first == 5)||(distances.first == 6)));
            REQUIRE(((distances.first == 5)||(distances.second == 5)));

            pair<int64_t, int64_t> distances1= di.distToCommonAncestor(snarl3, snarl2, pos2).first;
            REQUIRE(((distances1.first == 4)||(distances1.second == 4)));
            REQUIRE(((distances1.first == 4)||(distances1.first == 1)));
            REQUIRE(((distances1.second == 4)||(distances1.second == 1)));

            REQUIRE(di.distance(snarl3, snarl3,pos3, pos2) == 2);

            REQUIRE(di.distance(snarl3, snarl3, pos1,pos2) == 6);

            REQUIRE(di.distance(snarl2, snarl3, pos5,pos2) == 5);
            REQUIRE(di.distance(snarl2, snarl3, pos6,pos2) == 5);
            REQUIRE(di.distance(snarl3, snarl2, pos2, pos5) == 5);
            REQUIRE(di.distance(snarl1, snarl1, pos7, pos8r) == 5);
            REQUIRE(di.distance(snarl1, snarl1, pos7, pos8) == 5);
            REQUIRE(di.distance(snarl3, snarl2, pos4, pos9) == -1);

            REQUIRE(distance(&graph, pos3, pos2) == 2);
            REQUIRE(distance(&graph, pos1,pos2) == 6);
            REQUIRE(distance(&graph, pos5,pos2) == 5);
            REQUIRE(distance(&graph, pos6,pos2) == 5);
            REQUIRE(distance(&graph, pos2, pos5) == 5);
            REQUIRE(distance(&graph, pos7, pos8r) == 5);
            REQUIRE(distance(&graph, pos7, pos8) == 5);
            REQUIRE(distance(&graph, pos4, pos9) == -1);


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
        TestDistanceIndex di (&graph, &snarl_manager);


        SECTION("Create distance index") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);
            const Chain* chain = snarl_manager.chain_of(snarl2);

            pair<id_t, bool> start1 = make_pair(snarl1->start().node_id(), 
                                                    snarl1->start().backward());
            pair<id_t, bool> start1r = make_pair(snarl1->start().node_id(), 
                                                   !snarl1->start().backward());
            pair<id_t, bool> end1 = make_pair(snarl1->end().node_id(), 
                                                    snarl1->end().backward());
            pair<id_t, bool> end1r = make_pair(snarl1->end().node_id(), 
                                                    !snarl1->end().backward());
            TestDistanceIndex::SnarlDistances& sd1 = di.snarlIndex.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlDistances& sd2 = di.snarlIndex.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
 
            pair<id_t, bool> start5 = make_pair(snarl5->start().node_id(), 
                                                    snarl5->start().backward());
            pair<id_t, bool> start5r = make_pair(snarl5->start().node_id(), 
                                                   !snarl5->start().backward());
            pair<id_t, bool> end5 = make_pair(snarl5->end().node_id(), 
                                                    snarl5->end().backward());
            pair<id_t, bool> end5r = make_pair(snarl5->end().node_id(), 
                                                   !snarl5->end().backward());
            TestDistanceIndex::SnarlDistances& sd5 = di.snarlIndex.at(make_pair(snarl5->start().node_id(),
                                     snarl5->start().backward()));

            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 0) == 0);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 1) == 1);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 2) == 2);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 3) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 4) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 5) == 6);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 0) == -1);

            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(8, false))
                                                                         == 3);
            REQUIRE(sd1.snarlDistance(make_pair(8, true), make_pair(1, true))
                                                                          == 4);
            REQUIRE(di.chainIndex.size() == 1);
    
            TestDistanceIndex::ChainDistances& cd = di.chainIndex.at(get_start_of(*chain).node_id());
            REQUIRE(cd.chainDistance(start2, start2) == 0);
            #ifdef print
            di.printSelf();
            #endif
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(5, false)) == 2);
            REQUIRE(cd.chainDistance(make_pair(5, true), make_pair(2, true)) == 4);

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

            REQUIRE(di.distance(snarl2, snarl5, pos1, pos2) == 7);
            REQUIRE(di.distance(snarl5, snarl2, pos2, pos1) == 7);
            REQUIRE(di.distance(snarl2, snarl2, pos3, pos1) == 3);
            REQUIRE(di.distance(snarl2, snarl2, pos1, pos3) == 3);
            REQUIRE(di.distance(snarl2, snarl5, pos3, pos2) == 6);
            REQUIRE(di.distance(snarl2, snarl5, pos3, pos4) == 6);
            REQUIRE(di.distance(snarl1, snarl5, pos5, pos2) == 9);
            REQUIRE(di.distance(snarl1, snarl1, pos5, pos6) == 7);

            REQUIRE(distance(&graph, pos1, pos2) == 7);
            REQUIRE(distance(&graph, pos2, pos1) == 7);
            REQUIRE(distance(&graph, pos3, pos1) == 3);
            REQUIRE(distance(&graph, pos1, pos3) == 3);
            REQUIRE(distance(&graph, pos3, pos2) == 6);
            REQUIRE(distance(&graph, pos3, pos4) == 6);
            REQUIRE(distance(&graph, pos5, pos2) == 9);
            REQUIRE(distance(&graph, pos5, pos6) == 7);

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
        TestDistanceIndex di(&graph, &snarl_manager);


        SECTION("Create distance index") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);
            const Chain* chain = snarl_manager.chain_of(snarl2);

            pair<id_t, bool> start1 = make_pair(snarl1->start().node_id(), 
                                                    snarl1->start().backward());
            pair<id_t, bool> start1r = make_pair(snarl1->start().node_id(), 
                                                   !snarl1->start().backward());
            pair<id_t, bool> end1 = make_pair(snarl1->end().node_id(), 
                                                    snarl1->end().backward());
            pair<id_t, bool> end1r = make_pair(snarl1->end().node_id(), 
                                                    !snarl1->end().backward());
            TestDistanceIndex::SnarlDistances& sd1 = di.snarlIndex.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlDistances& sd2 = di.snarlIndex.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
 
            pair<id_t, bool> start5 = make_pair(snarl5->start().node_id(), 
                                                    snarl5->start().backward());
            pair<id_t, bool> start5r = make_pair(snarl5->start().node_id(), 
                                                   !snarl5->start().backward());
            pair<id_t, bool> end5 = make_pair(snarl5->end().node_id(), 
                                                    snarl5->end().backward());
            pair<id_t, bool> end5r = make_pair(snarl5->end().node_id(), 
                                                   !snarl5->end().backward());
            TestDistanceIndex::SnarlDistances& sd5 = di.snarlIndex.at(make_pair(snarl5->start().node_id(),
                                     snarl5->start().backward()));

    
            TestDistanceIndex::ChainDistances& cd = di.chainIndex.at(get_start_of(*chain).node_id());
            REQUIRE(cd.chainDistance(start2, start2) == 0);
            #ifdef print
            di.printSelf();
            #endif
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 0) == 0);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 1) == 1);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 2) == 2);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 3) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 4) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 5) == 6);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 0) == 9);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 1) == 3);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 2) == 1);
            REQUIRE(di.checkChainLoopRev(get_start_of(*chain).node_id(), 0) == -1);
            REQUIRE(di.checkChainLoopRev(get_start_of(*chain).node_id(), 1) == -1);

            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(8, false))
                                                                         == 3);
            REQUIRE(sd1.snarlDistance(make_pair(8, true), make_pair(1, true))
                                                                          == 4);
            REQUIRE(sd5.snarlDistance(make_pair(5, false), make_pair(5, true))
                                                                          == 3);
            REQUIRE(sd5.snarlDistance(make_pair(6, false), make_pair(6, true))
                                                                        == 12);
            REQUIRE(sd5.snarlDistance(make_pair(5, false), make_pair(6, true))
                                                                        == 5);
            REQUIRE(sd2.snarlDistance(make_pair(3, false), make_pair(4, true))
                                                                         == 7);
            REQUIRE(di.chainIndex.size() == 1);
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(5, false)) == 2);
            REQUIRE(cd.chainDistance(make_pair(5, true), make_pair(2, true)) == 4);

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
            pos_t pos7 = make_pos_t(6, false, 10);
            pos_t pos8 = make_pos_t(3, false, 0);

            REQUIRE(di.distance(snarl2, snarl5, pos1, pos2) == 7);
            REQUIRE(di.distance(snarl5, snarl2, pos2, pos1) == 7);
            REQUIRE(di.distance(snarl2, snarl2, pos3, pos1) == 3);
            REQUIRE(di.distance(snarl2, snarl2, pos1, pos3) == 3);
            REQUIRE(di.distance(snarl2, snarl5, pos3, pos2) == 6);
            REQUIRE(di.distance(snarl2, snarl5, pos3, pos4) == 6);
            REQUIRE(di.distance(snarl1, snarl5, pos5, pos2) == 9);
            REQUIRE(di.distance(snarl5, snarl1, pos2, pos5) == 9);
            REQUIRE(di.distance(snarl1, snarl1, pos5, pos6) == 7);
            REQUIRE(di.distance(snarl2, snarl5, pos1, pos7) == 10);
            REQUIRE(di.distance(snarl2, snarl2, pos1, pos8) == 10);

            REQUIRE(distance(&graph, pos1, pos2) == 7);
            REQUIRE(distance(&graph, pos2, pos1) == 7);
            REQUIRE(distance(&graph, pos3, pos1) == 3);
            REQUIRE(distance(&graph, pos1, pos3) == 3);
            REQUIRE(distance(&graph, pos3, pos2) == 6);
            REQUIRE(distance(&graph, pos3, pos4) == 6);
            REQUIRE(distance(&graph, pos5, pos2) == 9);
            REQUIRE(distance(&graph, pos2, pos5) == 9);
            REQUIRE(distance(&graph, pos5, pos6) == 7);
            REQUIRE(distance(&graph, pos1, pos7) == 10);
            REQUIRE(distance(&graph, pos1, pos8) == 10);
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
        TestDistanceIndex di (&graph, &snarl_manager);


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
            TestDistanceIndex::SnarlDistances& sd1 = 
                     di.snarlIndex.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlDistances& sd2 = 
                   di.snarlIndex.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
 
            pair<id_t, bool> start5 = make_pair(snarl5->start().node_id(), 
                                                    snarl5->start().backward());
            pair<id_t, bool> start5r = make_pair(snarl5->start().node_id(), 
                                                   !snarl5->start().backward());
            pair<id_t, bool> end5 = make_pair(snarl5->end().node_id(), 
                                                    snarl5->end().backward());
            pair<id_t, bool> end5r = make_pair(snarl5->end().node_id(), 
                                                   !snarl5->end().backward());
            TestDistanceIndex::SnarlDistances& sd5 = 
                di.snarlIndex.at(make_pair(snarl5->start().node_id(),
                                     snarl5->start().backward()));

            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 0) == 0);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 1) == 1);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 2) == 2);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 3) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 4) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 5) == 6);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 0) == -1);

            REQUIRE(di.checkChainLoopRev(get_start_of(*chain).node_id(),1)== 3);
            REQUIRE(di.checkChainLoopRev(get_start_of(*chain).node_id(),0)==-1);

            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(8, false))
                                                                         == 3);
            REQUIRE(sd1.snarlDistance(make_pair(8, true), make_pair(1, true))
                                                                          == 4);
            REQUIRE(sd5.snarlDistance(make_pair(6, true), make_pair(6, false))
                                                                          == 7);
            REQUIRE(di.chainIndex.size() == 2);
    
            TestDistanceIndex::ChainDistances& cd = di.chainIndex.at(get_start_of(*chain).node_id());
            REQUIRE(cd.chainDistance(start2, start2) == 0);
            #ifdef print
            di.printSelf();
            #endif
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(5, false)) == 2);
            REQUIRE(cd.chainDistance(make_pair(5, true), make_pair(2, true)) == 4);

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

            REQUIRE(di.distance(snarl2, snarl5, pos1, pos2) == 7);
            REQUIRE(di.distance(snarl5, snarl2, pos2, pos1) == 7);
            REQUIRE(di.distance(snarl2, snarl2, pos3, pos1) == 3);
            REQUIRE(di.distance(snarl2, snarl2, pos1, pos3) == 3);
            REQUIRE(di.distance(snarl2, snarl5, pos3, pos2) == 6);
            REQUIRE(di.distance(snarl2, snarl5, pos3, pos4) == 6);
            REQUIRE(di.distance(snarl1, snarl5, pos5, pos2) == 9);
            REQUIRE(di.distance(snarl1, snarl1, pos5, pos6) == 7);
            REQUIRE(di.distance(snarl2, snarl8, pos3, pos7) == 11);

            REQUIRE(distance(&graph, pos1, pos2) == 7);
            REQUIRE(distance(&graph, pos2, pos1) == 7);
            REQUIRE(distance(&graph, pos3, pos1) == 3);
            REQUIRE(distance(&graph, pos1, pos3) == 3);
            REQUIRE(distance(&graph, pos3, pos2) == 6);
            REQUIRE(distance(&graph, pos3, pos4) == 6);
            REQUIRE(distance(&graph, pos5, pos2) == 9);
            REQUIRE(distance(&graph, pos5, pos6) == 7);
            REQUIRE(distance(&graph, pos3, pos7) == 11);
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
        TestDistanceIndex di (&graph, &snarl_manager);


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
            TestDistanceIndex::SnarlDistances& sd1 = 
                     di.snarlIndex.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            pair<id_t, bool> start2 = make_pair(snarl2->start().node_id(), 
                                                    snarl2->start().backward());
            pair<id_t, bool> start2r = make_pair(snarl2->start().node_id(), 
                                                   !snarl2->start().backward());
            pair<id_t, bool> end2 = make_pair(snarl2->end().node_id(), 
                                                    snarl2->end().backward());
            pair<id_t, bool> end2r = make_pair(snarl2->end().node_id(), 
                                                   !snarl2->end().backward());
            TestDistanceIndex::SnarlDistances& sd2 = 
                   di.snarlIndex.at(make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));
 
            pair<id_t, bool> start4 = make_pair(snarl4->start().node_id(), 
                                                    snarl4->start().backward());
            pair<id_t, bool> start4r = make_pair(snarl4->start().node_id(), 
                                                   !snarl4->start().backward());
            pair<id_t, bool> end4 = make_pair(snarl4->end().node_id(), 
                                                    snarl4->end().backward());
            pair<id_t, bool> end4r = make_pair(snarl4->end().node_id(), 
                                                   !snarl4->end().backward());
            TestDistanceIndex::SnarlDistances& sd4 = 
                di.snarlIndex.at(make_pair(snarl4->start().node_id(),
                                     snarl4->start().backward()));
            TestDistanceIndex::SnarlDistances& sd7 = 
                di.snarlIndex.at(make_pair(snarl7->start().node_id(),
                                     snarl7->start().backward()));
    
            TestDistanceIndex::ChainDistances& cd = 
                               di.chainIndex.at(get_start_of(*chain).node_id());

            #ifdef print
            di.printSelf();
            #endif

            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 0) == 0);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 1) == 1);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 2) == 1);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 3) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 4) == 5);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 5) == 6);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 6) == 10);
            REQUIRE(di.checkChainDist(get_start_of(*chain).node_id(), 7) == 11);

            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 0) == 15);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 1) == 10);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 2) == 5);
            REQUIRE(di.checkChainLoopFd(get_start_of(*chain).node_id(), 3) == -1);

            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(10, false))
                                                                         == 3);
            REQUIRE(di.chainIndex.size() == 1);

            REQUIRE(cd.chainDistance(start2, start2) == 0);
            REQUIRE(cd.chainDistance(make_pair(2, false), make_pair(7, false)) == 5);

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
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);
            pos_t pos8r = make_pos_t(8, true, 0);
            pos_t pos9 = make_pos_t(9, false, 0);

            REQUIRE(di.distance(snarl2, snarl7, pos4, pos8) == 5);
            REQUIRE(di.distance(snarl4, snarl7, pos4, pos8) == 5);
            REQUIRE(di.distance(snarl7, snarl4, pos8, pos4) == 5);
            REQUIRE(di.distance(snarl1, snarl4, pos1, pos4) == 6);
            REQUIRE(di.distance(snarl1, snarl4, pos1, pos4) == 6);
            REQUIRE(di.distance(snarl7, snarl7, pos8, pos8r) == 4);
            REQUIRE(di.distance(snarl7, snarl7, pos7, pos8r) == 2);
            REQUIRE(di.distance(snarl7, snarl7, pos7, pos8) == 2);
            REQUIRE(di.distance(snarl7, snarl7, pos8r, pos7) == 2);
            REQUIRE(di.distance(snarl2, snarl7, pos2, pos8) == 7);
            REQUIRE(di.distance(snarl4, snarl4, pos6, pos5) == 10);

            REQUIRE(distance(&graph, pos4, pos8) == 5);
            REQUIRE(distance(&graph, pos4, pos8) == 5);
            REQUIRE(distance(&graph, pos8, pos4) == 5);
            REQUIRE(distance(&graph, pos1, pos4) == 6);
            REQUIRE(distance(&graph, pos1, pos4) == 6);
            REQUIRE(distance(&graph, pos8, pos8r) == 4);
            REQUIRE(distance(&graph, pos7, pos8r) == 2);
            REQUIRE(distance(&graph, pos7, pos8) == 2);
            REQUIRE(distance(&graph, pos8r, pos7) == 2);
            REQUIRE(distance(&graph, pos2, pos8) == 7);
            REQUIRE(distance(&graph, pos6, pos5) == 10);
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
        TestDistanceIndex di (&graph, &snarl_manager);


        SECTION("Create distance index") {
            #ifdef print
            di.printSelf();
            #endif
            auto snarl1 = snarl_manager.into_which_snarl(1, true);
            auto snarl2 = snarl_manager.into_which_snarl(2, true);
            auto snarl6 = snarl_manager.into_which_snarl(6, true);
            auto snarl7 = snarl_manager.into_which_snarl(7, true);

            const Chain* chain = snarl_manager.chain_of(snarl2);
            TestDistanceIndex::SnarlDistances& sd1 = 
                     di.snarlIndex.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            TestDistanceIndex::SnarlDistances& sd2 = 
                     di.snarlIndex.at(make_pair(snarl2->start().node_id(),
                                              snarl2->start().backward()));
            TestDistanceIndex::SnarlDistances& sd6 = 
                     di.snarlIndex.at(make_pair(snarl6->start().node_id(),
                                              snarl6->start().backward()));
            TestDistanceIndex::SnarlDistances& sd7 = 
                     di.snarlIndex.at(make_pair(snarl7->start().node_id(),
                                              snarl7->start().backward()));  


            TestDistanceIndex::ChainDistances& cd = di.chainIndex.at(get_start_of(*chain).node_id());

            REQUIRE(sd1.snarlDistance(make_pair(1, true), make_pair(1, false))
                                                                         == 3);
            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(1, true))
                                                                         == -1);
            REQUIRE(sd2.snarlDistance(make_pair(2, true), make_pair(1, true))
                                                                         == 1);
            REQUIRE(sd2.snarlDistance(make_pair(2, true), make_pair(1, false))
                                                                         == -1);
            REQUIRE(sd2.snarlDistance(make_pair(1, false), make_pair(3, false))
                                                                         == 6);
            REQUIRE(sd2.snarlDistance(make_pair(3, true), make_pair(1, true))
                                                                         == 1);
            REQUIRE(sd2.snarlDistance(make_pair(3, true), make_pair(1, false))
                                                                         == -1);
            REQUIRE(sd2.snarlDistance(make_pair(3, true), make_pair(3, false))
                                                                         == 7);
            REQUIRE(sd6.snarlDistance(make_pair(3, false), make_pair(4, false))
                                                                         == 1);
            REQUIRE(sd6.snarlDistance(make_pair(3, false), make_pair(6, false))
                                                                         == 4);
            REQUIRE(sd6.snarlDistance(make_pair(6, true), make_pair(3, true))
                                                                         == 4);
            REQUIRE(sd6.snarlDistance(make_pair(5, true), make_pair(3, false))
                                                                         == -1);
             
            REQUIRE(cd.chainDistance(make_pair(6, false), make_pair(3, false)) 
                                                                    == 4);
            REQUIRE(cd.chainDistance(make_pair(3, true), make_pair(6, true)) 
                                                                    == 4);
            REQUIRE(cd.chainDistance(make_pair(6, false), make_pair(3, true)) 
                                                                    == 11);
            REQUIRE(cd.chainDistance(make_pair(2, true), make_pair(2, false)) 
                                                                    == 7);
            REQUIRE(cd.chainDistance(make_pair(6, false), make_pair(2, false)) 
                                                                    == 11);
            REQUIRE(cd.chainDistance(make_pair(6, false), make_pair(2, true)) 
                                                                    == -1);
            REQUIRE(di.chainIndex.size() == 1);

            REQUIRE(sd7.snarlDistance(make_pair(7, true), make_pair(6, false))
                                                                         == 1);
            REQUIRE(sd7.snarlDistance(make_pair(7, true), make_pair(6, true))
                                                                         == 1);
            REQUIRE(sd7.snarlDistance(make_pair(7, true), make_pair(7, false))
                                                                         == 9);
            REQUIRE(sd7.snarlDistance(make_pair(6, true), make_pair(7, false))
                                                                        == 12);
            REQUIRE(sd7.snarlDistance(make_pair(6, false), make_pair(7, true))
                                                                        == -1);

            }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, true);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, true);
            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, true);
            const Snarl* snarl7 = snarl_manager.into_which_snarl(7, true);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);

            REQUIRE(di.distance(snarl1, snarl7, pos1, pos7) == 5);
            REQUIRE(di.distance(snarl2, snarl6, pos2, pos3) == 8);
            REQUIRE(di.distance(snarl2, snarl2, pos2, pos3) == 8);
            REQUIRE(di.distance(snarl6, snarl7, pos4, pos7) == 6);
            REQUIRE(di.distance(snarl1, snarl6, pos1, pos5) == 5);

            REQUIRE(distance(&graph, pos1, pos7) == 5);
            REQUIRE(distance(&graph, pos2, pos3) == 8);
            REQUIRE(distance(&graph, pos2, pos3) == 8);
            REQUIRE(distance(&graph, pos4, pos7) == 6);
            REQUIRE(distance(&graph, pos1, pos5) == 5);
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
        TestDistanceIndex di (&graph, &snarl_manager);
        
        SECTION( "Create distance index" ) {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlDistances& sd1 = di.snarlIndex.at(make_pair(
                  snarl1->start().node_id(),snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlDistances& sd2 = di.snarlIndex.at(make_pair(
                       snarl2->start().node_id(),snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            TestDistanceIndex::SnarlDistances& sd3 = di.snarlIndex.at(make_pair(
                       snarl3->start().node_id(),snarl3->start().backward()));

            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);
            TestDistanceIndex::SnarlDistances& sd6 = di.snarlIndex.at(make_pair(
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
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5r = make_pos_t(5, false, 11);
            pos_t pos7 = make_pos_t(7, false, 0);

            REQUIRE(di.distance(snarl3, snarl3, pos4, pos5) == 6);
            REQUIRE(di.distance(snarl3, snarl3, pos5, pos5r) == 2);
            REQUIRE(di.distance(snarl3, snarl6, pos5, pos7) == 11);

            REQUIRE(distance(&graph, pos4, pos5) == 6);
            REQUIRE(distance(&graph, pos5, pos5r) == 2);
            REQUIRE(distance(&graph, pos5, pos7) == 11);
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
        TestDistanceIndex di (&graph, &snarl_manager);
        
        SECTION( "Create distance index" ) {

#ifdef print
            di.printSelf();
#endif
            const vector<const Snarl*> topSnarls = snarl_manager.top_level_snarls();
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlDistances& sd1 = di.snarlIndex.at(
                                    make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlDistances& sd2 = di.snarlIndex.at(
                            make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            TestDistanceIndex::SnarlDistances& sd3 = di.snarlIndex.at(
                                       make_pair(snarl3->start().node_id(),
                                                 snarl3->start().backward()));

            NetGraph ng = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);

            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(2, false)) == 3); 
            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(2, true)) == 6); 
            REQUIRE(sd1.snarlDistanceShort(&graph , &ng, make_pair(1, false), make_pair(2, false)) == 0); 
            REQUIRE(sd1.snarlDistance(make_pair(2, true), make_pair(1,true)) == 3); 

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
            pos_t pos8 = make_pos_t(8, true, 1);
            pos_t pos9 = make_pos_t(6, false, 0);
            pos_t pos10 = make_pos_t(4, false, 0);

            pair<int64_t, int64_t> distances = di.distToCommonAncestor(snarl3, snarl1, pos4).first;
            REQUIRE(((distances.first == 5)||(distances.first == 6)));
            REQUIRE(((distances.first == 5)||(distances.second == 5)));

            pair<int64_t, int64_t> distances1= di.distToCommonAncestor(snarl3, snarl2, pos2).first;
            REQUIRE(((distances1.first == 4)||(distances1.second == 4)));
            REQUIRE(((distances1.first == 4)||(distances1.first == 1)));
            REQUIRE(((distances1.second == 4)||(distances1.second == 1)));

            REQUIRE(di.distance(snarl3, snarl3,pos3, pos2) == 2);

            REQUIRE(di.distance(snarl3, snarl3, pos1,pos2) == 6);

            REQUIRE(di.distance(snarl2, snarl3, pos5,pos2) == 5);
            REQUIRE(di.distance(snarl2, snarl3, pos6,pos2) == 5);
            REQUIRE(di.distance(snarl3, snarl2, pos2, pos5) == 5);
            REQUIRE(di.distance(snarl1, snarl1, pos7, pos8) == 5);
            REQUIRE(di.distance(snarl3, snarl2, pos4, pos9) == 8);
            REQUIRE(di.distance(snarl3, snarl2, pos10, pos9) == 9);
            REQUIRE(di.distance(snarl2, snarl3, pos9, pos10) == 9);

            REQUIRE(distance(&graph,pos3, pos2) == 2);
            REQUIRE(distance(&graph, pos1,pos2) == 6);
            REQUIRE(distance(&graph, pos5,pos2) == 5);
            REQUIRE(distance(&graph, pos6,pos2) == 5);
            REQUIRE(distance(&graph, pos2, pos5) == 5);
            REQUIRE(distance(&graph, pos7, pos8) == 5);
            REQUIRE(distance(&graph, pos4, pos9) == 8);
            REQUIRE(distance(&graph, pos10, pos9) == 9);
            REQUIRE(distance(&graph, pos9, pos10) == 9);
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
        TestDistanceIndex di (&graph, &snarl_manager);
        
        SECTION( "Create distance index" ) {

#ifdef print
            di.printSelf();
#endif
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlDistances& sd1 = di.snarlIndex.at(
                                    make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlDistances& sd2 = di.snarlIndex.at(
                            make_pair(snarl2->start().node_id(),
                                     snarl2->start().backward()));


            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(7, false)) == 3);
            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(7, false)) == 3);
            REQUIRE(sd1.snarlDistance(make_pair(
                       snarl2->start().node_id(), !snarl2->start().backward()), 
                      make_pair(snarl2->start().node_id(), 
                                           snarl2->start().backward())) == -1);
            REQUIRE(sd2.snarlDistance(make_pair(4, true), make_pair(4, false)) == 6);
            REQUIRE(sd1.snarlDistance(make_pair(5, true), make_pair(5, false)) == 13);
            REQUIRE(sd1.snarlDistance(make_pair(5, true), make_pair(7, false)) == 14);
            REQUIRE(sd1.snarlDistance(make_pair(7, true), make_pair(7, false)) == 13);

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
            pos_t pos6 = make_pos_t(6, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);

            REQUIRE(di.distance(snarl1, snarl2, pos1, pos4) == 5);
            REQUIRE(di.distance(snarl1, snarl1, pos5, pos6) == 12);
            REQUIRE(di.distance(snarl2, snarl1, pos2, pos7) == 7);

            REQUIRE(distance(&graph, pos1, pos4) == 5);
            REQUIRE(distance(&graph, pos5, pos6) == 12);
            REQUIRE(distance(&graph, pos2, pos7) == 7);
        }
    }//End test case
    TEST_CASE("Nested unary snarls", "[dist]") {
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
        Edge* e2 = graph.create_edge(n1, n4);
        Edge* e3 = graph.create_edge(n1, n8);
        Edge* e4 = graph.create_edge(n2, n3);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n7);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n7, n7, false, true);
        Edge* e10 = graph.create_edge(n8, n4, true, false);
        Edge* e11 = graph.create_edge(n8, n2, true, false);



        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        TestDistanceIndex di (&graph, &snarl_manager);


        SECTION("Create distance index") {
            #ifdef print
            di.printSelf();
            #endif
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, true);
/*            const Snarl* snarl4 = snarl_manager.into_which_snarl(4, false);

            TestDistanceIndex::SnarlDistances& sd1 = 
                     di.snarlIndex.at(make_pair(snarl1->start().node_id(),
                                              snarl1->start().backward()));
 
            TestDistanceIndex::SnarlDistances& sd2 = 
                     di.snarlIndex.at(make_pair(snarl2->start().node_id(),
                                              snarl2->start().backward()));
            TestDistanceIndex::SnarlDistances& sd4 = 
                     di.snarlIndex.at(make_pair(snarl4->start().node_id(),
                                              snarl4->start().backward()));
           REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(1, true))
                                                                       == 15);
           REQUIRE(sd1.snarlDistance(make_pair(1, true), make_pair(1, false))
                                                                       == -1);
*/
            }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, true);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, true);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, true);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);

            REQUIRE(di.distance(snarl1, snarl1, pos1, pos1) == 1);
            REQUIRE(di.distance(snarl1, snarl2, pos1, pos2) == 4);
            REQUIRE(di.distance(snarl1, snarl2, pos1, pos8) == 4);
            REQUIRE(di.distance(snarl1, snarl5, pos1, pos5) == 8);
            REQUIRE(di.distance(snarl3, snarl5, pos3, pos5) == -1);
            REQUIRE(di.distance(snarl3, snarl2, pos3, pos8) == 3);

            REQUIRE(distance(&graph, pos1, pos1) == 1);
            REQUIRE(distance(&graph, pos1, pos2) == 4);
            REQUIRE(distance(&graph, pos1, pos8) == 4);
            REQUIRE(distance(&graph, pos1, pos5) == 8);
            REQUIRE(distance(&graph, pos3, pos5) == -1);
            REQUIRE(distance(&graph, pos3, pos8) == 3);
 
        }
    }//end test case

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
        TestDistanceIndex di (&graph, &snarl_manager);
        
        SECTION( "Create distance index" ) {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlDistances& sd1 = di.snarlIndex.at(make_pair(
                  snarl1->start().node_id(),snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlDistances& sd2 = di.snarlIndex.at(make_pair(
                       snarl2->start().node_id(),snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            TestDistanceIndex::SnarlDistances& sd3 = di.snarlIndex.at(make_pair(
                       snarl3->start().node_id(),snarl3->start().backward()));

            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);
            TestDistanceIndex::SnarlDistances& sd6 = di.snarlIndex.at(make_pair(
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
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5r = make_pos_t(5, false, 11);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);
            pos_t pos9 = make_pos_t(9, false, 0);
            pos_t pos9r = make_pos_t(9, false, 1);
            pos_t pos11 = make_pos_t(11, false, 0);
            pos_t pos12 = make_pos_t(12, false, 0);

            REQUIRE(di.distance(snarl2, snarl2, pos2, pos9r) == 2);
            REQUIRE(di.distance(snarl3, snarl2, pos3, pos9) == 4);
            REQUIRE(di.distance(snarl3, snarl3, pos4, pos5) == 6);
            REQUIRE(di.distance(snarl3, snarl3, pos4, pos5) == 6);
            REQUIRE(di.distance(snarl6, snarl6, pos7, pos12) == 14);
            REQUIRE(di.distance(snarl6, snarl2, pos8, pos11) == 8);


            REQUIRE(distance(&graph, pos2, pos9r) == 2);
            REQUIRE(distance(&graph, pos3, pos9) == 4);
            REQUIRE(distance(&graph, pos4, pos5) == 6);
            REQUIRE(distance(&graph, pos4, pos5) == 6);
            REQUIRE(distance(&graph, pos7, pos12) == 14);
            REQUIRE(distance(&graph, pos8, pos11) == 8);

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
        TestDistanceIndex di (&graph, &snarl_manager);
        
        SECTION( "Create distance index" ) {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            TestDistanceIndex::SnarlDistances& sd1 = di.snarlIndex.at(make_pair(
                  snarl1->start().node_id(),snarl1->start().backward()));

            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            TestDistanceIndex::SnarlDistances& sd2 = di.snarlIndex.at(make_pair(
                       snarl2->start().node_id(),snarl2->start().backward()));

            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);
            TestDistanceIndex::SnarlDistances& sd3 = di.snarlIndex.at(make_pair(
                       snarl3->start().node_id(),snarl3->start().backward()));

            const Snarl* snarl6 = snarl_manager.into_which_snarl(6, false);
            TestDistanceIndex::SnarlDistances& sd6 = di.snarlIndex.at(make_pair(
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
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos4r = make_pos_t(4, false, 3);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos5r = make_pos_t(5, false, 11);
            pos_t pos7 = make_pos_t(7, false, 0);
            pos_t pos8 = make_pos_t(8, false, 0);
            pos_t pos9 = make_pos_t(9, false, 0);
            pos_t pos9r = make_pos_t(9, false, 1);
            pos_t pos10 = make_pos_t(10, false, 0);
            pos_t pos11 = make_pos_t(11, false, 0);
            pos_t pos12 = make_pos_t(12, false, 0);
            pos_t pos14 = make_pos_t(14, false, 0);
            pos_t pos16 = make_pos_t(16, false, 0);

            REQUIRE(di.distance(snarl2, snarl9, pos2, pos10) == 3);
            REQUIRE(di.distance(snarl3, snarl3, pos4, pos5) == 6);
            REQUIRE(di.distance(snarl3, snarl3, pos4r, pos5r) == 12);
            REQUIRE(di.distance(snarl13, snarl9, pos14, pos10) == 6);
            REQUIRE(di.distance(snarl13, snarl3, pos14, pos3) == 8);
            REQUIRE(di.distance(snarl2, snarl3, pos16, pos3) == 4);

            REQUIRE(distance(&graph, pos2, pos10) == 3);
            REQUIRE(distance(&graph, pos4, pos5) == 6);
            REQUIRE(distance(&graph, pos4r, pos5r) == 12);
            REQUIRE(distance(&graph, pos14, pos10) == 6);
            REQUIRE(distance(&graph, pos14, pos3) == 8);
            REQUIRE(distance(&graph, pos16, pos3) == 4);
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
        TestDistanceIndex di (&graph, &snarl_manager);


        SECTION("Create distance index") {

            #ifdef print
                di.printSelf();
            #endif
 
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Chain* chain = snarl_manager.chain_of(snarl1);
            TestDistanceIndex::ChainDistances& cd = di.chainIndex.at(get_start_of(*chain).node_id());
            REQUIRE(cd.chainDistance(make_pair(1, true), make_pair(8, true)) == 3);
            REQUIRE(cd.chainDistance(make_pair(8, false), make_pair(1, false)) == 3);
            REQUIRE(cd.chainDistanceShort(&graph, make_pair(8, false), make_pair(1, false)) == 2);
          
        }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
            const Snarl* snarl5 = snarl_manager.into_which_snarl(5, false);
            const Snarl* snarl9 = snarl_manager.into_which_snarl(9, false);

            pos_t pos1 = make_pos_t(1, false, 1);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos6 = make_pos_t(6, true, 0);
            pos_t pos7 = make_pos_t(7, true, 0);

            REQUIRE(di.distance(snarl2, snarl5, pos4, pos6) == 8);
            REQUIRE(di.distance(snarl2, snarl5, pos2, pos7) == 6);

            REQUIRE(distance(&graph, pos4, pos6) == 8);
            REQUIRE(distance(&graph, pos2, pos7) == 6);
        }
    }//end test case

    TEST_CASE("Two tip start", "[dist]") {
        VG graph;

        Node* n1 = graph.create_node("G");
        Node* n2 = graph.create_node("TA");
        Node* n3 = graph.create_node("GGG");

        Edge* e1 = graph.create_edge(n1, n3);
        Edge* e2 = graph.create_edge(n2, n3);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        TestDistanceIndex di (&graph, &snarl_manager);


       SECTION ("Distance functions") {
            #ifdef print
                di.printSelf();
            #endif
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, true);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(2, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);

            REQUIRE(di.distance(snarl1, snarl3, pos1, pos3) == 2);

        }
    }//end test case
    TEST_CASE("Random test", "[dist]") {
        for (int i = 0; i < 100; i++) {
            //1000 different graphs
            VG graph = randomGraph(1000, 20, 100); 

            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls(); 
            TestDistanceIndex di (&graph, &snarl_manager);
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
 
                off_t offset1 = rand() % node1->sequence().size();
                off_t offset2 = rand() % node2->sequence().size();

                pos_t pos1 = make_pos_t(nodeID1, rand()%2 == 0, offset1 );
                pos_t pos2 = make_pos_t(nodeID2, rand()%2 == 0, offset2 );
 
                if (!(nodeID1 != snarl1->start().node_id() && 
                    (snarl_manager.into_which_snarl(nodeID1, false) != NULL ||
                      snarl_manager.into_which_snarl(nodeID1, true) != NULL)) &&
                    ! (nodeID2 != snarl2->start().node_id() &&
                        (snarl_manager.into_which_snarl(nodeID2, false) != NULL
                       || snarl_manager.into_which_snarl(nodeID2, true) != NULL
                  ))){

                    //If the nodes aren't child snarls

                    int64_t myDist = di.distance(snarl1, snarl2,pos1, pos2);
                    int64_t actDist = distance(&graph, pos1, pos2);
                    bool passed = myDist == actDist;

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

                        cerr << "Actual distance: " << actDist << "    " <<
                                "Guessed distance: " << myDist << endl;
                    }
                    REQUIRE(passed);
                }
            }
             
        }
    } //end test case
}

}
