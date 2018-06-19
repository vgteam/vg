#include <stdio.h>
#include <iostream>
#include <set>
#include "json2pb.h"
#include "vg.hpp"
#include "catch.hpp"
#include "snarls.hpp"
#include "distance.hpp"
#include "testTraversal.cpp"
#include "genotypekit.hpp"
#include <fstream>
 
//#define print
namespace vg {
namespace unittest {

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
            sd1.printSelf();
            sd2.printSelf();
            sd3.printSelf();
#endif
            NetGraph ng = NetGraph(snarl1->start(), snarl1->end(), snarl_manager.chains_of(snarl1), &graph);

            REQUIRE(sd1.snarlDistance(start1, make_pair(2, false)) == 3); 
            REQUIRE(sd1.snarlDistance(start1, make_pair(2, true)) == -1); 
            REQUIRE(sd1.snarlDistanceShort(&ng, &snarl_manager, start1, make_pair(2, false)) == 0); 
            REQUIRE(sd1.snarlDistance(start2r, make_pair(1,true)) == 1);         
 

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
            REQUIRE(di.distance(snarl1, snarl1, pos7, pos8) == 4);
            REQUIRE(di.distance(snarl3, snarl2, pos4, pos9) == -1);

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
            cd.printSelf();
            sd1.printSelf();
            sd2.printSelf();
            sd5.printSelf();
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
            REQUIRE(di.distance(snarl1, snarl1, pos5, pos6) == 4);
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
            cd.printSelf();
            sd1.printSelf();
            sd2.printSelf();
            sd5.printSelf();
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
            REQUIRE(di.distance(snarl1, snarl1, pos5, pos6) == 4);
            REQUIRE(di.distance(snarl2, snarl5, pos1, pos7) == 10);
            REQUIRE(di.distance(snarl2, snarl2, pos1, pos8) == 10);
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
            REQUIRE(di.chainIndex.size() == 2);
    
            TestDistanceIndex::ChainDistances& cd = di.chainIndex.at(get_start_of(*chain).node_id());
            REQUIRE(cd.chainDistance(start2, start2) == 0);
            #ifdef print
            cd.printSelf();
            sd1.printSelf();
            sd2.printSelf();
            sd5.printSelf();
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
            REQUIRE(di.distance(snarl1, snarl1, pos5, pos6) == 4);
            REQUIRE(di.distance(snarl2, snarl8, pos3, pos7) == 11);
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
            cd.printSelf();
            sd1.printSelf();
            sd2.printSelf();
            sd4.printSelf();
            sd7.printSelf();
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
            pos_t pos4 = make_pos_t(4, false, 1);
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
            REQUIRE(di.distance(snarl7, snarl7, pos7, pos9) == 6);
 
        }
        }//end test case
    TEST_CASE("Top level loop breaks", "[dist]") {
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
        //Edge* e9 = graph.create_edge(n7, n7, false, true);
        Edge* e9 = graph.create_edge(n1, n1, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls(); 
        TestDistanceIndex di (&graph, &snarl_manager);


        SECTION("Create distance index") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);


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
 
            pair<id_t, bool> start3 = make_pair(snarl3->start().node_id(), 
                                                    snarl3->start().backward());
            pair<id_t, bool> start3r = make_pair(snarl3->start().node_id(), 
                                                   !snarl3->start().backward());
            pair<id_t, bool> end3 = make_pair(snarl3->end().node_id(), 
                                                    snarl3->end().backward());
            pair<id_t, bool> end3r = make_pair(snarl3->end().node_id(), 
                                                   !snarl3->end().backward());
            TestDistanceIndex::SnarlDistances& sd3 = 
                   di.snarlIndex.at(make_pair(snarl3->start().node_id(),
                                     snarl3->start().backward()));

            #ifdef print
            sd1.printSelf();
            sd3.printSelf();
            #endif

            REQUIRE(sd1.snarlDistance(make_pair(1, false), make_pair(7, false))
                                                                         == 4);
            REQUIRE(di.chainIndex.size() == 0);


            }
        SECTION ("Distance functions") {
            const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
            const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

            pos_t pos1 = make_pos_t(1, false, 0);
            pos_t pos2 = make_pos_t(2, false, 0);
            pos_t pos3 = make_pos_t(3, false, 0);
            pos_t pos4 = make_pos_t(4, false, 0);
            pos_t pos5 = make_pos_t(5, false, 0);
            pos_t pos7 = make_pos_t(7, false, 0);

            REQUIRE(di.distance(snarl1, snarl1, pos1, pos7) == 5);
            //REQUIRE(di.distance(snarl1, snarl3, pos2, pos3) == 5);
 
        }
        }//end test case
    TEST_CASE( "Shortest path exits chain or snarl","[dist]" ) {
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
            sd1.printSelf();
            sd2.printSelf();
            sd3.printSelf();
            sd6.printSelf();
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
            pos_t pos7 = make_pos_t(7, false, 0);

            //REQUIRE(di.distance(snarl3, snarl3, pos4, pos5) == 4);
            //REQUIRE(di.distance(snarl3, snarl6, pos5, pos7) == 10);

        }
        }//End test case



    }

}
