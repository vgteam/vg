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
                
                const vector<const Snarl*>& middle_level_snarls = snarl_manager.children_of(top_level_snarls[0]);
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
                
                const vector<const Snarl*>& bottom_level_snarls = snarl_manager.children_of(middle_level_snarls[0]);
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
                
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n1->id() && snarl->end().node_id() == n8->id()) ||
                        (snarl->start().node_id() == n8->id() && snarl->end().node_id() == n1->id())) {
                        contents = snarl_manager.shallow_contents(snarl, graph, false);
                        break;
                    }
                }
                
                unordered_set<Node*>& nodes = contents.first;
                unordered_set<Edge*>& edges = contents.second;
                
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
                
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n1->id() && snarl->end().node_id() == n8->id()) ||
                        (snarl->start().node_id() == n8->id() && snarl->end().node_id() == n1->id())) {
                        contents = snarl_manager.deep_contents(snarl, graph, false);
                        break;
                    }
                }
                
                unordered_set<Node*>& nodes = contents.first;
                unordered_set<Edge*>& edges = contents.second;
                
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
                
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n2->id() && snarl->end().node_id() == n7->id()) ||
                        (snarl->start().node_id() == n7->id() && snarl->end().node_id() == n2->id())) {
                        contents = snarl_manager.shallow_contents(snarl, graph, false);
                        break;
                    }
                }
                
                unordered_set<Node*>& nodes = contents.first;
                unordered_set<Edge*>& edges = contents.second;
                
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
                
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n2->id() && snarl->end().node_id() == n7->id()) ||
                        (snarl->start().node_id() == n7->id() && snarl->end().node_id() == n2->id())) {
                        contents = snarl_manager.deep_contents(snarl, graph, false);
                        break;
                    }
                }
                
                unordered_set<Node*>& nodes = contents.first;
                unordered_set<Edge*>& edges = contents.second;
                
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
            
            SECTION( "SnarlManager can flip the internal representation of a snarl") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("C");
                Node* n3 = graph.create_node("A");
                Node* n4 = graph.create_node("CTGA");
                
                Edge* e1 = graph.create_edge(n1, n2);
                Edge* e2 = graph.create_edge(n1, n3);
                Edge* e3 = graph.create_edge(n2, n4);
                Edge* e4 = graph.create_edge(n3, n4);
                
                Snarl snarl;
                snarl.mutable_start()->set_node_id(n1->id());
                snarl.mutable_end()->set_node_id(n4->id());
                snarl.set_type(ULTRABUBBLE);
                
                list<Snarl> snarls;
                snarls.push_back(snarl);
                
                SnarlManager snarl_manager(snarls.begin(), snarls.end());
                
                const Snarl* snarl_ref = snarl_manager.top_level_snarls()[0];
                
                auto id_start = snarl_ref->start().node_id();
                auto id_end = snarl_ref->end().node_id();
                
                snarl_manager.flip(snarl_ref);
                
                REQUIRE(snarl_ref->start().node_id() == id_end);
                REQUIRE(snarl_ref->end().node_id() == id_start);
            }
            
            SECTION( "SnarlManager can create boundary indices") {
                
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
                
                const Snarl* top_snarl = snarl_manager.top_level_snarls()[0];
                const Snarl* child_snarl = snarl_manager.children_of(top_snarl)[0];
                
                auto boundary_index = snarl_manager.child_boundary_index(top_snarl, graph);
                auto start_index = snarl_manager.child_start_index(top_snarl, graph);
                auto end_index = snarl_manager.child_end_index(top_snarl, graph);
                
                NodeTraversal child_start = to_node_traversal(child_snarl->start(), graph);
                NodeTraversal child_end = to_rev_node_traversal(child_snarl->end(), graph);
                
                REQUIRE(boundary_index.count(child_start));
                REQUIRE(boundary_index.count(child_end));
                REQUIRE(start_index.count(child_start));
                REQUIRE(end_index.count(child_end));
                REQUIRE(boundary_index.size() == 2);
                REQUIRE(start_index.size() == 1);
                REQUIRE(end_index.size() == 1);
            }
            
            
            
            SECTION( "SnarlManager can correctly extract the full contents of a reversing-edge snarl") {
                
                string graph_json = R"(
                {
                    "node": [
                             {
                             "sequence": "TTTTTG",
                             "id": 6462830
                             },
                             {
                             "sequence": "AAAAAAAAAAAAAA",
                             "id": 8480141
                             },
                             {
                             "sequence": "A",
                             "id": 6462831
                             },
                             {
                             "sequence": "T",
                             "id": 6462832
                             },
                             {
                             "sequence": "G",
                             "id": 8480142
                             },
                             {
                             "sequence": "A",
                             "id": 8480143
                             }
                             ],
                    "edge": [
                             {
                             "to": 8480141,
                             "from": 6462830,
                             "from_start": true
                             },
                             {
                             "to": 6462831,
                             "from": 6462830
                             },
                             {
                             "to": 6462832,
                             "from": 6462830
                             },
                             {
                             "to": 8480142,
                             "from": 8480141
                             },
                             {
                             "to": 8480143,
                             "from": 8480141
                             }
                             ]
                }
                )";
                
                VG graph;
                
                // Load up the graph
                Graph g;
                json2pb(g, graph_json.c_str(), graph_json.size());
                graph.extend(g);
                
                // Define the one snarl
                Snarl snarl1;
                snarl1.mutable_start()->set_node_id(6462830);
                snarl1.mutable_start()->set_backward(true);
                snarl1.mutable_end()->set_node_id(8480141);
                snarl1.set_type(ULTRABUBBLE);
                
                list<Snarl> snarls;
                snarls.push_back(snarl1);
                
                SnarlManager snarl_manager(snarls.begin(), snarls.end());
                
                // Find the snarl again
                const Snarl* snarl = snarl_manager.top_level_snarls()[0];
                
                // Get its contents
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents = snarl_manager.deep_contents(snarl, graph, true);
                
                // We need the right snarl
                REQUIRE(snarl->start().node_id() == 6462830);
                REQUIRE(snarl->start().backward());
                REQUIRE(snarl->end().node_id() == 8480141);
                REQUIRE(!snarl->end().backward());
                
                // And it needs to contain just those two nodes and the edges connecting them.
                REQUIRE(contents.first.size() == 2);
                REQUIRE(contents.second.size() == 1);
                
            }  
        }
    }
}
