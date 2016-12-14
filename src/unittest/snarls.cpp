//
//  snarls.cpp
//
//  Unit tests for SnarlManager and related functions
//

#include <stdio.h>
#include <iostream>
#include <set>
#include "json2pb.h"
#include "vg.pb.h"
#include "catch.hpp"
#include "snarls.hpp"
#include "genotypekit.hpp"

namespace vg {
    namespace unittest {
        TEST_CASE( "SnarlManager functions return expected answers",
                  "[sites][bubbles][snarls]" ) {
            
            SECTION( "SnarlManager can be constructed with cactus ultrabubbles") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                CactusUltrabubbleFinder bubble_finder(graph, "");
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
            }
            
            SECTION( "SnarlManager can correctly navigate tree relationships") {
                
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
                
                CactusUltrabubbleFinder bubble_finder(graph, "");
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                const vector<const Snarl*>& top_level_snarls = snarl_manager.top_level_snarls();
                
                REQUIRE(top_level_snarls.size() == 1);
                
                bool top_level_correct = ((top_level_snarls[0]->start().node_id() == n1->id() &&
                                           top_level_snarls[0]->end().node_id() == n8->id() &&
                                           !top_level_snarls[0]->start().backward() &&
                                           !top_level_snarls[0]->end().backward()) ||
                                          (top_level_snarls[0]->start().node_id() == n8->id() &&
                                           top_level_snarls[0]->end().node_id() == n1->id() &&
                                           top_level_snarls[0]->start().backward() &&
                                           top_level_snarls[0]->end().backward()));
                REQUIRE(top_level_correct);
                
                const vector<const Snarl*>& middle_level_snarls = snarl_manager.children_of(*top_level_snarls[0]);
                REQUIRE(middle_level_snarls.size() == 1);
                
                bool middle_level_correct = ((middle_level_snarls[0]->start().node_id() == n2->id() &&
                                              middle_level_snarls[0]->end().node_id() == n7->id() &&
                                              !middle_level_snarls[0]->start().backward() &&
                                              !middle_level_snarls[0]->end().backward()) ||
                                             (middle_level_snarls[0]->start().node_id() == n7->id() &&
                                              middle_level_snarls[0]->end().node_id() == n2->id() &&
                                              middle_level_snarls[0]->start().backward() &&
                                              middle_level_snarls[0]->end().backward()));
                REQUIRE(middle_level_correct);
                
                const vector<const Snarl*>& bottom_level_snarls = snarl_manager.children_of(*middle_level_snarls[0]);
                REQUIRE(bottom_level_snarls.size() == 1);
                
                bool bottom_level_correct = ((bottom_level_snarls[0]->start().node_id() == n3->id() &&
                                              bottom_level_snarls[0]->end().node_id() == n5->id() &&
                                              !bottom_level_snarls[0]->start().backward() &&
                                              !bottom_level_snarls[0]->end().backward()) ||
                                             (bottom_level_snarls[0]->start().node_id() == n5->id() &&
                                              bottom_level_snarls[0]->end().node_id() == n3->id() &&
                                              bottom_level_snarls[0]->start().backward() &&
                                              bottom_level_snarls[0]->end().backward()));
                REQUIRE(bottom_level_correct);
                
            }
            
            SECTION( "SnarlManager can correctly extract the shallow contents of a snarl") {
                
                VG graph;
                
                Node* n0 = graph.create_node("CA");
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("C");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("GCA");
                Node* n6 = graph.create_node("T");
                Node* n7 = graph.create_node("G");
                Node* n8 = graph.create_node("AGTA");
                Node* n9 = graph.create_node("CCC");
                
                Edge* e0 = graph.create_edge(n0, n1);
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
                Edge* e11 = graph.create_edge(n8, n9);
                
                CactusUltrabubbleFinder bubble_finder(graph, "");
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                const vector<const Snarl*>& top_level_snarls = snarl_manager.top_level_snarls();
                
                pair<vector<Node*>, vector<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n1->id() && snarl->end().node_id() == n8->id()) ||
                        (snarl->start().node_id() == n8->id() && snarl->end().node_id() == n1->id())) {
                        contents = snarl_manager.shallow_contents(*snarl, graph);
                    }
                }
                
                set<Node*> nodes;
                set<Edge*> edges;
                
                for (Node* node : contents.first) {
                    nodes.insert(node);
                }
                
                for (Edge* edge : contents.second) {
                    edges.insert(edge);
                }
                
                REQUIRE(nodes.size() == 2);
                
                REQUIRE(nodes.count(n2));
                REQUIRE(nodes.count(n7));
                
                REQUIRE(edges.size() == 3);
                
                REQUIRE(edges.count(e1));
                REQUIRE(edges.count(e2));
                REQUIRE(edges.count(e10));
            }
            
            SECTION( "SnarlManager can correctly extract the full contents of a snarl") {
                
                VG graph;
                
                Node* n0 = graph.create_node("CA");
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("C");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("GCA");
                Node* n6 = graph.create_node("T");
                Node* n7 = graph.create_node("G");
                Node* n8 = graph.create_node("AGTA");
                Node* n9 = graph.create_node("CCC");
                
                Edge* e0 = graph.create_edge(n0, n1);
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
                Edge* e11 = graph.create_edge(n8, n9);
                
                CactusUltrabubbleFinder bubble_finder(graph, "");
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                const vector<const Snarl*>& top_level_snarls = snarl_manager.top_level_snarls();
                
                pair<vector<Node*>, vector<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n1->id() && snarl->end().node_id() == n8->id()) ||
                        (snarl->start().node_id() == n8->id() && snarl->end().node_id() == n1->id())) {
                        contents = snarl_manager.deep_contents(*snarl, graph);
                    }
                }
                
                set<Node*> nodes;
                set<Edge*> edges;
                
                for (Node* node : contents.first) {
                    nodes.insert(node);
                }
                
                for (Edge* edge : contents.second) {
                    edges.insert(edge);
                }
                
                REQUIRE(nodes.size() == 6);
                
                REQUIRE(nodes.count(n2));
                REQUIRE(nodes.count(n3));
                REQUIRE(nodes.count(n4));
                REQUIRE(nodes.count(n5));
                REQUIRE(nodes.count(n6));
                REQUIRE(nodes.count(n7));
                
                REQUIRE(edges.size() == 10);
                
                REQUIRE(edges.count(e1));
                REQUIRE(edges.count(e2));
                REQUIRE(edges.count(e3));
                REQUIRE(edges.count(e4));
                REQUIRE(edges.count(e5));
                REQUIRE(edges.count(e6));
                REQUIRE(edges.count(e7));
                REQUIRE(edges.count(e8));
                REQUIRE(edges.count(e9));
                REQUIRE(edges.count(e10));
            }
            
            SECTION( "SnarlManager can correctly extract the shallow contents of a snarl containing cycles and tips") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("C");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("GCA");
                Node* n6 = graph.create_node("T");
                Node* n7 = graph.create_node("G");
                Node* n8 = graph.create_node("AGTA");
                
                Edge* e1 = graph.create_edge(n1, n2);
                Edge* e2 = graph.create_edge(n2, n7);
                Edge* e3 = graph.create_edge(n2, n4);
                Edge* e4 = graph.create_edge(n3, n4);
                Edge* e5 = graph.create_edge(n4, n5);
                Edge* e6 = graph.create_edge(n4, n6);
                Edge* e7 = graph.create_edge(n5, n6);
                Edge* e8 = graph.create_edge(n6, n2, false, true);
                Edge* e9 = graph.create_edge(n6, n7);
                Edge* e10 = graph.create_edge(n7, n8);
                
                Snarl snarl1;
                snarl1.mutable_start()->set_node_id(n2->id());
                snarl1.mutable_end()->set_node_id(n7->id());
                snarl1.set_type(UNCLASSIFIED);
                
                Snarl snarl2;
                snarl2.mutable_start()->set_node_id(n4->id());
                snarl2.mutable_end()->set_node_id(n6->id());
                snarl2.set_type(ULTRABUBBLE);
                *snarl2.mutable_parent() = snarl1;
                
                list<Snarl> snarls;
                snarls.push_back(snarl1);
                snarls.push_back(snarl2);
                
                SnarlManager snarl_manager(snarls.begin(), snarls.end());
                
                const vector<const Snarl*>& top_level_snarls = snarl_manager.top_level_snarls();
                
                pair<vector<Node*>, vector<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n2->id() && snarl->end().node_id() == n7->id()) ||
                        (snarl->start().node_id() == n7->id() && snarl->end().node_id() == n2->id())) {
                        contents = snarl_manager.shallow_contents(*snarl, graph);
                    }
                }
                
                set<Node*> nodes;
                set<Edge*> edges;
                
                for (Node* node : contents.first) {
                    nodes.insert(node);
                }
                
                for (Edge* edge : contents.second) {
                    edges.insert(edge);
                }
                
                REQUIRE(nodes.size() == 3);
                
                REQUIRE(nodes.count(n3));
                REQUIRE(nodes.count(n4));
                REQUIRE(nodes.count(n6));
                
                REQUIRE(edges.size() == 5);
                
                REQUIRE(edges.count(e2));
                REQUIRE(edges.count(e3));
                REQUIRE(edges.count(e4));
                REQUIRE(edges.count(e8));
                REQUIRE(edges.count(e9));
            }
            
            
            SECTION( "SnarlManager can correctly extract the full contents of a snarl containing cycles and tips") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("C");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("GCA");
                Node* n6 = graph.create_node("T");
                Node* n7 = graph.create_node("G");
                Node* n8 = graph.create_node("AGTA");
                
                Edge* e1 = graph.create_edge(n1, n2);
                Edge* e2 = graph.create_edge(n2, n7);
                Edge* e3 = graph.create_edge(n2, n4);
                Edge* e4 = graph.create_edge(n3, n4);
                Edge* e5 = graph.create_edge(n4, n5);
                Edge* e6 = graph.create_edge(n4, n6);
                Edge* e7 = graph.create_edge(n5, n6);
                Edge* e8 = graph.create_edge(n6, n2, false, true);
                Edge* e9 = graph.create_edge(n6, n7);
                Edge* e10 = graph.create_edge(n7, n8);
                
                Snarl snarl1;
                snarl1.mutable_start()->set_node_id(n2->id());
                snarl1.mutable_end()->set_node_id(n7->id());
                snarl1.set_type(UNCLASSIFIED);
                
                Snarl snarl2;
                snarl2.mutable_start()->set_node_id(n4->id());
                snarl2.mutable_end()->set_node_id(n6->id());
                snarl2.set_type(ULTRABUBBLE);
                *snarl2.mutable_parent() = snarl1;
                
                list<Snarl> snarls;
                snarls.push_back(snarl1);
                snarls.push_back(snarl2);
                
                SnarlManager snarl_manager(snarls.begin(), snarls.end());
                
                const vector<const Snarl*>& top_level_snarls = snarl_manager.top_level_snarls();
                
                pair<vector<Node*>, vector<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n2->id() && snarl->end().node_id() == n7->id()) ||
                        (snarl->start().node_id() == n7->id() && snarl->end().node_id() == n2->id())) {
                        contents = snarl_manager.deep_contents(*snarl, graph);
                    }
                }
                
                set<Node*> nodes;
                set<Edge*> edges;
                
                for (Node* node : contents.first) {
                    nodes.insert(node);
                }
                
                for (Edge* edge : contents.second) {
                    edges.insert(edge);
                }
                
                REQUIRE(nodes.size() == 4);
                
                REQUIRE(nodes.count(n3));
                REQUIRE(nodes.count(n4));
                REQUIRE(nodes.count(n5));
                REQUIRE(nodes.count(n6));
                
                REQUIRE(edges.size() == 8);
                
                REQUIRE(edges.count(e2));
                REQUIRE(edges.count(e3));
                REQUIRE(edges.count(e4));
                REQUIRE(edges.count(e5));
                REQUIRE(edges.count(e6));
                REQUIRE(edges.count(e7));
                REQUIRE(edges.count(e8));
                REQUIRE(edges.count(e9));
            }
        }
    }
}