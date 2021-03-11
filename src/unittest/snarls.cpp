//
//  snarls.cpp
//
//  Unit tests for SnarlManager and related functions
//

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <set>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "catch.hpp"
#include "random_graph.hpp"
#include "randomness.hpp"
#include "../snarls.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../genotypekit.hpp"
#include "../traversal_finder.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>

//#define debug

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
    
        TEST_CASE( "NetGraph can allow traversal of a simple net graph",
                  "[snarls][netgraph]" ) {
        
        
            // This graph will have a snarl from 1 to 8, a snarl from 2 to 7,
            // and a snarl from 3 to 5, all nested in each other.
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
            
            // Define the snarls for the top level
            Snarl top_snarl;
            top_snarl.mutable_start()->set_node_id(n1->id());
            top_snarl.mutable_end()->set_node_id(n8->id());
            
            // We have a chain with a snarl in it
            vector<vector<pair<Snarl, bool>>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back().first;
            
            // And that snarl has these characteristics
            nested_snarl1.mutable_start()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_node_id(n7->id());
            nested_snarl1.set_type(ULTRABUBBLE);
            *nested_snarl1.mutable_parent() = top_snarl;
            nested_snarl1.set_directed_acyclic_net_graph(true);
            nested_snarl1.set_start_self_reachable(false);
            nested_snarl1.set_end_self_reachable(false);
            nested_snarl1.set_start_end_reachable(true);
            
            // We have an empty vector of top-level unary snarls.
            // TODO: should we have any?
            vector<Snarl> top_unary_snarls;
            
            // Make a net graph
            NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, false);
            
            SECTION( "The top-level NetGraph has 3 nodes" ) {
                unordered_set<handle_t> nodes;
                size_t node_count = 0;
                net_graph.for_each_handle([&](const handle_t& handle) {
                    nodes.insert(handle);
                    node_count++;
                });
                REQUIRE(nodes.size() == node_count);
                
                REQUIRE(nodes.size() == 3);
                
                SECTION( "The nodes are numbered 1, 2, and 8" ) {
                    REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                    REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                    REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                }
            }
            
            SECTION( "The top-level NetGraph has 3 edges" ) {
                unordered_set<pair<handle_t, handle_t>> edges;
                for (auto& id : {1, 2, 8}) {
                    // Go through the nodes we should have manually.
                    handle_t handle = net_graph.get_handle(id, false);
                
                    // Save all the edges off of each node
                    net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                        edges.insert(net_graph.edge_handle(handle, other));
                    });
                    net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                        edges.insert(net_graph.edge_handle(other, handle));
                    });
                }
                
                REQUIRE(edges.size() == 3);
            }
        
        }
        
        TEST_CASE( "NetGraph can handle start-start connectivity",
                  "[snarls][netgraph]" ) {
        
        
            // This graph will have a snarl from 1 to 8, a snarl from 2 to 7,
            // and a snarl from 3 to 5, all nested in each other.
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
            // Add an extra reversing edge so you can leave out the start
            Edge* eRev = graph.create_edge(n3, n3, false, true);
            Edge* e6 = graph.create_edge(n3, n5);
            Edge* e7 = graph.create_edge(n4, n5);
            Edge* e8 = graph.create_edge(n5, n7);
            Edge* e9 = graph.create_edge(n6, n7);
            Edge* e10 = graph.create_edge(n7, n8);
            
            // Define the snarls for the top level
            Snarl top_snarl;
            top_snarl.mutable_start()->set_node_id(n1->id());
            top_snarl.mutable_end()->set_node_id(n8->id());
            
            // We have a chain with a snarl in it
            vector<vector<pair<Snarl, bool>>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back().first;
            
            // And that snarl has these characteristics
            nested_snarl1.mutable_start()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_node_id(n7->id());
            nested_snarl1.set_type(UNCLASSIFIED);
            *nested_snarl1.mutable_parent() = top_snarl;
            nested_snarl1.set_directed_acyclic_net_graph(true);
            // Because we can reverse inside a child, we are start-start reachable
            nested_snarl1.set_start_self_reachable(true);
            nested_snarl1.set_end_self_reachable(false);
            nested_snarl1.set_start_end_reachable(true);
            
            // We have an empty vector of top-level unary snarls.
            // TODO: should we have any?
            vector<Snarl> top_unary_snarls;
            
            SECTION( "A connectivity-ignoring net graph ignores connectivity" ) {
            
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, false);
                
                SECTION( "The top-level NetGraph has 3 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 3);
                    
                    SECTION( "The nodes are numbered 1, 2, and 8" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 3 edges" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 2, 8}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    
                    REQUIRE(edges.size() == 3);
                }
            }
            
            SECTION( "A connectivity-respecting net graph allows more traversals" ) {
            
                // Make a connectivity-respecting graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, true);
                
                SECTION( "The top-level NetGraph has 3 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 3);
                    
                    SECTION( "The nodes are numbered 1, 2, and 8" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 4 edges" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 2, 8}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    
                    REQUIRE(edges.size() == 4);
                    
#ifdef debug
                    for (auto& edge : edges) {
                        cerr << net_graph.get_id(edge.first) << " " << net_graph.get_is_reverse(edge.first) << " -> "
                            << net_graph.get_id(edge.second) << " " << net_graph.get_is_reverse(edge.second) << endl;
                    }
#endif
                    
                    SECTION( "One edge is from node 1 to node 2 in reverse" ) {
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(1, false),
                            net_graph.get_handle(2, true))) == 1);
                    }
                }
            
            }
        
        }
        
        TEST_CASE( "NetGraph can handle end-end connectivity",
                  "[snarls][netgraph]" ) {
        
        
            // This graph will have a snarl from 1 to 8, a snarl from 2 to 7,
            // and a snarl from 3 to 5, all nested in each other.
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
            // Add an extra reversing edge so you can leave out the end
            Edge* eRev = graph.create_edge(n4, n4, true, false);
            Edge* e6 = graph.create_edge(n3, n5);
            Edge* e7 = graph.create_edge(n4, n5);
            Edge* e8 = graph.create_edge(n5, n7);
            Edge* e9 = graph.create_edge(n6, n7);
            Edge* e10 = graph.create_edge(n7, n8);
            
            // Define the snarls for the top level
            Snarl top_snarl;
            top_snarl.mutable_start()->set_node_id(n1->id());
            top_snarl.mutable_end()->set_node_id(n8->id());
            
            // We have a chain with a snarl in it
            vector<vector<pair<Snarl, bool>>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back().first;
            
            // And that snarl has these characteristics
            nested_snarl1.mutable_start()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_node_id(n7->id());
            nested_snarl1.set_type(UNCLASSIFIED);
            *nested_snarl1.mutable_parent() = top_snarl;
            nested_snarl1.set_directed_acyclic_net_graph(true);
            // Because we can reverse inside a child, we are end-end reachable
            nested_snarl1.set_start_self_reachable(false);
            nested_snarl1.set_end_self_reachable(true);
            nested_snarl1.set_start_end_reachable(true);
            
            // We have an empty vector of top-level unary snarls.
            // TODO: should we have any?
            vector<Snarl> top_unary_snarls;
            
            SECTION( "A connectivity-ignoring net graph ignores connectivity" ) {
            
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, false);
                
                SECTION( "The top-level NetGraph has 3 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 3);
                    
                    SECTION( "The nodes are numbered 1, 2, and 8" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 3 edges" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 2, 8}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    
                    REQUIRE(edges.size() == 3);
                }
            }
            
            SECTION( "A connectivity-respecting net graph allows more traversals" ) {
            
                // Make a connectivity-respecting graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, true);
                
                SECTION( "The top-level NetGraph has 3 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 3);
                    
                    SECTION( "The nodes are numbered 1, 2, and 8" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 4 edges" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 2, 8}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    
                    REQUIRE(edges.size() == 4);
             
#ifdef debug
                    for (auto& edge : edges) {
                        cerr << net_graph.get_id(edge.first) << " " << net_graph.get_is_reverse(edge.first) << " -> "
                            << net_graph.get_id(edge.second) << " " << net_graph.get_is_reverse(edge.second) << endl;
                    }
#endif
                    
                    SECTION( "One edge is from node 2 in reverse to node 8" ) {
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(2, true),
                            net_graph.get_handle(8, false))) == 1);
                    }
                }
            
            }
        
        }
        
        TEST_CASE( "NetGraph can handle disconnected chain bounds",
                  "[snarls][netgraph]" ) {
        
        
            // We will have a snarl 1 to 2, and within it a chain of 3 to 4 and
            // 4 to 5 (both trivial)
            VG graph;
                
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            Node* n4 = graph.create_node("CTGA");
            Node* n5 = graph.create_node("GCA");
            
            Edge* e1 = graph.create_edge(n1, n2);
            Edge* e2 = graph.create_edge(n3, n4);
            Edge* e3 = graph.create_edge(n4, n5);
            
            for (bool do_start : {false, true}) {
            
                // Attach only start or end of chain to parent snarl
                Edge* e4 = do_start ? graph.create_edge(n1, n3) : graph.create_edge(n5, n2);
                
                // Define the snarls for the top level
                Snarl top_snarl;
                top_snarl.mutable_start()->set_node_id(n1->id());
                top_snarl.mutable_end()->set_node_id(n2->id());
                
                // We have a chain with two snarls in it
                vector<vector<pair<Snarl, bool>>> top_chains;
                top_chains.emplace_back();
                auto& top_chain1 = top_chains.back();
                top_chain1.emplace_back();
                top_chain1.emplace_back();
                auto& nested_snarl1 = top_chain1[0].first;
                auto& nested_snarl2 = top_chain1[1].first;
                
                // Which have these characteristics
                nested_snarl1.mutable_start()->set_node_id(n3->id());
                nested_snarl1.mutable_end()->set_node_id(n4->id());
                nested_snarl1.set_type(ULTRABUBBLE);
                *nested_snarl1.mutable_parent() = top_snarl;
                nested_snarl1.set_directed_acyclic_net_graph(true);
                nested_snarl1.set_start_self_reachable(false);
                nested_snarl1.set_end_self_reachable(false);
                nested_snarl1.set_start_end_reachable(true);
                
                nested_snarl2.mutable_start()->set_node_id(n4->id());
                nested_snarl2.mutable_end()->set_node_id(n5->id());
                nested_snarl2.set_type(ULTRABUBBLE);
                *nested_snarl2.mutable_parent() = top_snarl;
                nested_snarl2.set_directed_acyclic_net_graph(true);
                nested_snarl2.set_start_self_reachable(false);
                nested_snarl2.set_end_self_reachable(false);
                nested_snarl2.set_start_end_reachable(true);
                
                // We have an empty vector of top-level unary snarls.
                // TODO: These kind of won't exist anymore.
                vector<Snarl> top_unary_snarls;
                
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, false);
                
                SECTION( "The top-level NetGraph has 3 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 3);
                    
                    SECTION( "The nodes are numbered 1, 2, and 3" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(3, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 2 edges" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 2, 3}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    
                    REQUIRE(edges.size() == 2);
                }
            }
        }
        
        TEST_CASE( "NetGraph can handle no connectivity",
                  "[snarls][netgraph]" ) {
        
        
            // This graph will have a snarl from 1 to 8, a snarl from 2 to 7,
            // and a snarl from 3 to 5, all nested in each other.
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
            // Make node 7's start connect to what node 2's end connects to, so the 2 to 7 snarl has two unary children.
            Edge* e8 = graph.create_edge(n7, n3, true, false);
            Edge* e9 = graph.create_edge(n7, n6, true, false);
            Edge* e10 = graph.create_edge(n7, n8);
            
            // Define the snarls for the top level
            Snarl top_snarl;
            top_snarl.mutable_start()->set_node_id(n1->id());
            top_snarl.mutable_end()->set_node_id(n8->id());
            
            // We have a chain with a snarl in it
            vector<vector<pair<Snarl, bool>>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back().first;
            
            // And that snarl has these characteristics
            nested_snarl1.mutable_start()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_node_id(n7->id());
            nested_snarl1.set_type(UNCLASSIFIED);
            *nested_snarl1.mutable_parent() = top_snarl;
            nested_snarl1.set_directed_acyclic_net_graph(true);
            // We aren't actually traversable at all
            nested_snarl1.set_start_self_reachable(false);
            nested_snarl1.set_end_self_reachable(false);
            nested_snarl1.set_start_end_reachable(false);
            
            // We have an empty vector of top-level unary snarls.
            // TODO: should we have any?
            vector<Snarl> top_unary_snarls;
            
            SECTION( "A connectivity-ignoring net graph ignores connectivity" ) {
            
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, false);
                
                SECTION( "The top-level NetGraph has 3 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 3);
                    
                    SECTION( "The nodes are numbered 1, 2, and 8" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 3 edges" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 2, 8}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    
                    REQUIRE(edges.size() == 3);
                }
            }
            
            SECTION( "A connectivity-respecting net graph doesn't allow through-traversal" ) {
            
                // Make a connectivity-respecting graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, true);
                
                SECTION( "The top-level NetGraph has 3 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 3);
                    
                    SECTION( "The nodes are numbered 1, 2, and 8" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 3 edges looking around node 2" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 8}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    
                    REQUIRE(edges.size() == 3);
             
#ifdef debug
                    for (auto& edge : edges) {
                        cerr << net_graph.get_id(edge.first) << " " << net_graph.get_is_reverse(edge.first) << " -> "
                            << net_graph.get_id(edge.second) << " " << net_graph.get_is_reverse(edge.second) << endl;
                    }
#endif
                }
                
                SECTION( "The top-level NetGraph has 0 edges looking at node 2" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {2}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    
                    REQUIRE(edges.size() == 0);
             
#ifdef debug
                    for (auto& edge : edges) {
                        cerr << net_graph.get_id(edge.first) << " " << net_graph.get_is_reverse(edge.first) << " -> "
                            << net_graph.get_id(edge.second) << " " << net_graph.get_is_reverse(edge.second) << endl;
                    }
#endif
                }
            
            }
        
        }
        
        TEST_CASE( "NetGraph finds all edges correctly when traversing in all directions",
                  "[snarls][netgraph]" ) {
        
        
            // This graph will have a snarl from 1 to 8, a snarl from 2 to 7,
            // and a snarl from 3 to 5, all nested in each other.
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
            // Add an extra reversing edge so you can leave out the start
            Edge* eRev = graph.create_edge(n3, n3, false, true);
            // Add an extra reversing edge so you can leave out the end
            Edge* eRev2 = graph.create_edge(n5, n5, false, true);
            Edge* e6 = graph.create_edge(n3, n5);
            Edge* e7 = graph.create_edge(n4, n5);
            Edge* e8 = graph.create_edge(n5, n7);
            Edge* e9 = graph.create_edge(n6, n7);
            Edge* e10 = graph.create_edge(n7, n8);
            
            // Define the snarls for the top level
            Snarl top_snarl;
            top_snarl.mutable_start()->set_node_id(n1->id());
            top_snarl.mutable_end()->set_node_id(n8->id());
            
            // We have a chain with a snarl in it
            vector<vector<pair<Snarl, bool>>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back().first;
            
            // And that snarl has these characteristics
            nested_snarl1.mutable_start()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_node_id(n7->id());
            nested_snarl1.set_type(UNCLASSIFIED);
            *nested_snarl1.mutable_parent() = top_snarl;
            nested_snarl1.set_directed_acyclic_net_graph(true);
            // Because we can reverse inside a child, we are start-start and end-end reachable
            nested_snarl1.set_start_self_reachable(true);
            nested_snarl1.set_end_self_reachable(true);
            nested_snarl1.set_start_end_reachable(true);
            
            // We have an empty vector of top-level unary snarls.
            // TODO: should we have any?
            vector<Snarl> top_unary_snarls;
            
            // Make a connectivity-respecting graph
            NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, true);
        
            // Make a handle to the nested child snarl
            handle_t nested_fwd = net_graph.get_handle(n2->id(), false);
            
            // ANd reverse it
            auto nested_rev = net_graph.flip(nested_fwd);
            
            SECTION( "Forward handle sees 1 fwd and 8 rev on its left" ) {
                unordered_set<handle_t> seen;
                net_graph.follow_edges(nested_fwd, true, [&](const handle_t& other) {
                    seen.insert(other);
                });
                
                REQUIRE(seen.size() == 2);
                REQUIRE(seen.count(net_graph.get_handle(n1->id(), false)) == 1);
                REQUIRE(seen.count(net_graph.get_handle(n8->id(), true)) == 1);
            }
            
            SECTION( "Forward handle sees 1 rev and 8 fwd on its right" ) {
                unordered_set<handle_t> seen;
                net_graph.follow_edges(nested_fwd, false, [&](const handle_t& other) {
                    seen.insert(other);
                });
                
                REQUIRE(seen.size() == 2);
                REQUIRE(seen.count(net_graph.get_handle(n1->id(), true)) == 1);
                REQUIRE(seen.count(net_graph.get_handle(n8->id(), false)) == 1);
            }
            
            SECTION( "Reverse handle sees 1 fwd and 8 rev on its left" ) {
                unordered_set<handle_t> seen;
                net_graph.follow_edges(nested_rev, true, [&](const handle_t& other) {
                    seen.insert(other);
                });
                
                REQUIRE(seen.size() == 2);
                REQUIRE(seen.count(net_graph.get_handle(n1->id(), false)) == 1);
                REQUIRE(seen.count(net_graph.get_handle(n8->id(), true)) == 1);
            }
            
            SECTION( "Reverse handle sees 1 rev and 8 fwd on its right" ) {
                unordered_set<handle_t> seen;
                net_graph.follow_edges(nested_rev, false, [&](const handle_t& other) {
                    seen.insert(other);
                });
                
                REQUIRE(seen.size() == 2);
                REQUIRE(seen.count(net_graph.get_handle(n1->id(), true)) == 1);
                REQUIRE(seen.count(net_graph.get_handle(n8->id(), false)) == 1);
            }
        
        }
        
        TEST_CASE( "NetGraph can handle nested unary snarls",
                  "[snarls][netgraph]" ) {
        
        
            // This graph will have a snarl from 1 to 8, a nested unary snarl on 2, and a nested unary snarl on 4.
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
            
            // Define the snarls for the top level
            Snarl top_snarl;
            top_snarl.mutable_start()->set_node_id(n1->id());
            top_snarl.mutable_end()->set_node_id(n8->id());
            
            // We have no chains
            vector<vector<pair<Snarl, bool>>> top_chains;
            
            // We havetwo child unary snarls.
            vector<Snarl> top_unary_snarls;

            top_unary_snarls.emplace_back();
            auto& nested_snarl1 = top_unary_snarls.back(); 
            
            // The first one is a tip of nodes 2 and 3, based at node 2
            nested_snarl1.mutable_start()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_backward(true);
            nested_snarl1.set_type(UNARY);
            *nested_snarl1.mutable_parent() = top_snarl;
            nested_snarl1.set_directed_acyclic_net_graph(true);
            // We aren't actually traversable at all
            nested_snarl1.set_start_self_reachable(false);
            nested_snarl1.set_end_self_reachable(false);
            nested_snarl1.set_start_end_reachable(false);
            
            top_unary_snarls.emplace_back();
            auto& nested_snarl2 = top_unary_snarls.back(); 
            
            // The second one is a big assemblage rooted at 4 with a child snarl and a reversing edge in it.
            nested_snarl2.mutable_start()->set_node_id(n4->id());
            nested_snarl2.mutable_end()->set_node_id(n4->id());
            nested_snarl2.mutable_end()->set_backward(true);
            nested_snarl2.set_type(UNARY);
            *nested_snarl2.mutable_parent() = top_snarl;
            // The reversal is actually in a child unary snarl.
            nested_snarl2.set_directed_acyclic_net_graph(true);
            // We are traversable in all directions
            nested_snarl2.set_start_self_reachable(true);
            nested_snarl2.set_end_self_reachable(true);
            nested_snarl2.set_start_end_reachable(true);
            
            SECTION( "A connectivity-ignoring net graph ignores connectivity" ) {
            
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, false);
                
                SECTION( "The top-level NetGraph has 4 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 4);
                    
                    SECTION( "The nodes are numbered 1, 2, 4, and 8" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(4, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 5 edges" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 2, 4, 8}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
#ifdef debug
                    for (auto& edge : edges) {
                        cerr << net_graph.get_id(edge.first) << " " << net_graph.get_is_reverse(edge.first) << " -> "
                            << net_graph.get_id(edge.second) << " " << net_graph.get_is_reverse(edge.second) << endl;
                    }
#endif
                    
                    REQUIRE(edges.size() == 5);
                    
                    // The edges are 1 to 8, 1 to 2, 1 to 4, 2 rev to 8, and 4 rev to 8.
                    SECTION( "All of the expected edges exist" ) {
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(1, false),
                            net_graph.get_handle(8, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(1, false),
                            net_graph.get_handle(2, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(1, false),
                            net_graph.get_handle(4, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(2, true),
                            net_graph.get_handle(8, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(4, true),
                            net_graph.get_handle(8, false))) == 1);
                    }
                }
            }
            
             SECTION( "A connectivity-respecting net graph respects connectivity" ) {
            
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, true);
                
                SECTION( "The top-level NetGraph has 4 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 4);
                    
                    SECTION( "The nodes are numbered 1, 2, 4, and 8" ) {
                        REQUIRE(nodes.count(net_graph.get_handle(1, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(2, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(4, false)) == 1);
                        REQUIRE(nodes.count(net_graph.get_handle(8, false)) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph has 5 edges when looking at its ends" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {1, 8}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
#ifdef debug
                    for (auto& edge : edges) {
                        cerr << net_graph.get_id(edge.first) << " " << net_graph.get_is_reverse(edge.first) << " -> "
                            << net_graph.get_id(edge.second) << " " << net_graph.get_is_reverse(edge.second) << endl;
                    }
#endif
                    
                    REQUIRE(edges.size() == 5);
                    
                    // The edges are 1 to 8, 1 to 2, 1 to 4, 2 rev to 8, and 4 rev to 8.
                    SECTION( "All of the expected edges exist" ) {
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(1, false),
                            net_graph.get_handle(8, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(1, false),
                            net_graph.get_handle(2, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(1, false),
                            net_graph.get_handle(4, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(2, true),
                            net_graph.get_handle(8, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(4, true),
                            net_graph.get_handle(8, false))) == 1);
                    }
                }
                
                SECTION( "The top-level NetGraph only has edges for the connected child when looking at its child snarls" ) {
                    unordered_set<pair<handle_t, handle_t>> edges;
                    for (auto& id : {2, 4}) {
                        // Go through the nodes we should have manually.
                        handle_t handle = net_graph.get_handle(id, false);
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
#ifdef debug
                    for (auto& edge : edges) {
                        cerr << net_graph.get_id(edge.first) << " " << net_graph.get_is_reverse(edge.first) << " -> "
                            << net_graph.get_id(edge.second) << " " << net_graph.get_is_reverse(edge.second) << endl;
                    }
#endif
                    
                    // We ignore the edges we would get if we could do a real
                    // predecessor lookup. In this mode we only allow real
                    // traversals, so looking left from entering the unary snarl
                    // gives you nothing, not the edges you could have come from.
                    REQUIRE(edges.size() == 2);
                    
                    // The edges are the fake traversal-allowing ones: 4 forward to 8 forward, and 4 forward to 1 reverse.
                    SECTION( "All of the expected edges exist" ) {
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(4, false),
                            net_graph.get_handle(8, false))) == 1);
                        REQUIRE(edges.count(net_graph.edge_handle(net_graph.get_handle(4, false),
                            net_graph.get_handle(1, true))) == 1);
                    }
                }
            }
            
        }
        
        TEST_CASE( "NetGraph can find hard-to-reach snarl contents",
                  "[snarls][netgraph]" ) {
        
        
            // This graph will have a snarl from 1 to 8, a nested unary snarl on 2, and a nested unary snarl on 4.
            // Node 3 will be attached only through these non-traversable unary snarls.
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
            Edge* e4 = graph.create_edge(n3, n2);
            Edge* e5 = graph.create_edge(n4, n5);
            Edge* e6 = graph.create_edge(n4, n6);
            Edge* e7 = graph.create_edge(n5, n7);
            Edge* e8 = graph.create_edge(n6, n7);
            Edge* e9 = graph.create_edge(n8, n4, true, false);
            Edge* e10 = graph.create_edge(n3, n4);
            
            // Define the snarls for the top level
            Snarl top_snarl;
            top_snarl.mutable_start()->set_node_id(n1->id());
            top_snarl.mutable_end()->set_node_id(n8->id());
            
            // We have no chains
            vector<vector<pair<Snarl, bool>>> top_chains;
            
            // We havetwo child unary snarls.
            vector<Snarl> top_unary_snarls;

            top_unary_snarls.emplace_back();
            auto& nested_snarl1 = top_unary_snarls.back(); 
            
            // The first one is a tip of just node 2, based at node 2
            nested_snarl1.mutable_start()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_node_id(n2->id());
            nested_snarl1.mutable_end()->set_backward(true);
            nested_snarl1.set_type(UNARY);
            *nested_snarl1.mutable_parent() = top_snarl;
            nested_snarl1.set_directed_acyclic_net_graph(true);
            // We aren't actually traversable at all
            nested_snarl1.set_start_self_reachable(false);
            nested_snarl1.set_end_self_reachable(false);
            nested_snarl1.set_start_end_reachable(false);
            
            top_unary_snarls.emplace_back();
            auto& nested_snarl2 = top_unary_snarls.back(); 
            
            // The second one is a big assemblage rooted at 4 that also isn't traversable
            nested_snarl2.mutable_start()->set_node_id(n4->id());
            nested_snarl2.mutable_end()->set_node_id(n4->id());
            nested_snarl2.mutable_end()->set_backward(true);
            nested_snarl2.set_type(UNARY);
            *nested_snarl2.mutable_parent() = top_snarl;
            nested_snarl2.set_directed_acyclic_net_graph(true);
            // We aren't actually traversable at all
            nested_snarl2.set_start_self_reachable(false);
            nested_snarl2.set_end_self_reachable(false);
            nested_snarl2.set_start_end_reachable(false);
            
            SECTION( "A connectivity-ignoring net graph lets you find all the nodes" ) {
            
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, false);
                
                SECTION( "The top-level NetGraph has 5 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 5);
                    
                }
                
            }
            
            SECTION( "A connectivity-respecting net graph still lets you find all the nodes" ) {
            
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, true);
                
                SECTION( "The top-level NetGraph has 5 nodes" ) {
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    REQUIRE(nodes.size() == node_count);
                    
                    REQUIRE(nodes.size() == 5);
                    
                }
                
            }
            
        }
    
        TEST_CASE( "SnarlManager functions return expected answers",
                  "[sites][snarls]" ) {
            
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
                
                CactusSnarlFinder bubble_finder(graph);
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
                
                CactusSnarlFinder bubble_finder(graph);
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
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                const vector<const Snarl*>& top_level_snarls = snarl_manager.top_level_snarls();
                
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n1->id() && snarl->end().node_id() == n8->id()) ||
                        (snarl->start().node_id() == n8->id() && snarl->end().node_id() == n1->id())) {
                        contents = pb_contents(graph, snarl_manager.shallow_contents(snarl, graph, false));
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
                
                CactusSnarlFinder bubble_finder(graph);
                SnarlManager snarl_manager = bubble_finder.find_snarls();
                
                const vector<const Snarl*>& top_level_snarls = snarl_manager.top_level_snarls();
                
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents;
                for (const Snarl* snarl : top_level_snarls) {
                    if ((snarl->start().node_id() == n1->id() && snarl->end().node_id() == n8->id()) ||
                        (snarl->start().node_id() == n8->id() && snarl->end().node_id() == n1->id())) {
                        contents = pb_contents(graph, snarl_manager.deep_contents(snarl, graph, false));
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
                        contents = pb_contents(graph, snarl_manager.shallow_contents(snarl, graph, false));
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
                        contents = pb_contents(graph, snarl_manager.deep_contents(snarl, graph, false));
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
                
                for (Snarl snarl : snarls) {
                    const Snarl* from_start = snarl_manager.into_which_snarl(snarl.start().node_id(),
                                                                             snarl.start().backward());
                    const Snarl* from_end = snarl_manager.into_which_snarl(snarl.end().node_id(),
                                                                           !snarl.end().backward());
                    REQUIRE(from_start != nullptr);
                    REQUIRE(from_end != nullptr);
                    
                    REQUIRE(from_start->start().node_id() == snarl.start().node_id());
                    REQUIRE(from_start->start().backward() == snarl.start().backward());
                    REQUIRE(from_start->end().node_id() == snarl.end().node_id());
                    REQUIRE(from_start->end().backward() == snarl.end().backward());
                    
                    REQUIRE(from_end->start().node_id() == snarl.start().node_id());
                    REQUIRE(from_end->start().backward() == snarl.start().backward());
                    REQUIRE(from_end->end().node_id() == snarl.end().node_id());
                    REQUIRE(from_end->end().backward() == snarl.end().backward());
                }
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
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents = pb_contents(graph, snarl_manager.deep_contents(snarl, graph, true));
                
                // We need the right snarl
                REQUIRE(snarl->start().node_id() == 6462830);
                REQUIRE(snarl->start().backward());
                REQUIRE(snarl->end().node_id() == 8480141);
                REQUIRE(!snarl->end().backward());

                // And it needs to contain just those two nodes and the edges connecting them.
                REQUIRE(contents.first.size() == 2);
                REQUIRE(contents.second.size() == 1);
                
            }
            
            SECTION( "SnarlManager does not include child snarls' edges in parent snarls") {
                
                // This graph is 3 nodes in a row, with two anchoring nodes on
                // the end, and an edge deleting the three in the middle and
                // just linking the anchoring nodes.
                string graph_json = R"(
                {
                  "node": [
                    {
                      "sequence": "A",
                      "id": 178895
                    },
                    {
                      "sequence": "G",
                      "id": 178896
                    },
                    {
                      "sequence": "A",
                      "id": 187209
                    },
                    {
                      "sequence": "TCTCAAAAAAAAAAAAAAAAAAAAAAAAAA",
                      "id": 178894
                    },
                    {
                      "sequence": "AATGTGTCTTCCTGGGT",
                      "id": 187208
                    }
                  ],
                  "edge": [
                    {
                      "from": 187209,
                      "to": 178895
                    },
                    {
                      "from": 178895,
                      "to": 178896
                    },
                    {
                      "from": 178896,
                      "to": 187208
                    },
                    {
                      "from": 178894,
                      "to": 187209
                    },
                    {
                      "from": 178894,
                      "to": 187208
                    }
                  ],
                  "path": [
                    {
                      "name": "5",
                      "mapping": [
                        {
                          "position": {
                            "node_id": 178894
                          },
                          "rank": 98372
                        },
                        {
                          "position": {
                            "node_id": 187209
                          },
                          "rank": 98373
                        },
                        {
                          "position": {
                            "node_id": 178895
                          },
                          "rank": 98374
                        },
                        {
                          "position": {
                            "node_id": 178896
                          },
                          "rank": 98375
                        },
                        {
                          "position": {
                            "node_id": 187208
                          },
                          "rank": 98376
                        }
                      ]
                    }
                  ]
                }
                )";
                
                // We have one parent snarl for the deletion, with two back-to-back trivial child snarls.
                string snarl1_json = R"({"type": 1, "end": {"node_id": 187208}, "start": {"node_id": 178894}})";
                string snarl2_json = R"({"type": 1, "end": {"node_id": 187209, "backward": true}, "start": {"node_id": 178895, "backward": true}, "parent": {"end": {"node_id": 187208}, "start": {"node_id": 178894}}})";
                string snarl3_json = R"({"type": 1, "end": {"node_id": 178896}, "start": {"node_id": 178895}, "parent": {"end": {"node_id": 187208}, "start": {"node_id": 178894}}})";
                
                VG graph;
                
                // Load up the graph
                Graph g;
                json2pb(g, graph_json.c_str(), graph_json.size());
                graph.extend(g);
                
                // Load the snarls
                Snarl snarl1, snarl2, snarl3;
                json2pb(snarl1, snarl1_json.c_str(), snarl1_json.size());
                json2pb(snarl2, snarl2_json.c_str(), snarl2_json.size());
                json2pb(snarl3, snarl3_json.c_str(), snarl3_json.size());
                
                // Put them in a list
                list<Snarl> snarls;
                snarls.push_back(snarl1);
                snarls.push_back(snarl2);
                snarls.push_back(snarl3);
                
                SnarlManager snarl_manager(snarls.begin(), snarls.end());
                
                // Find the root snarl again
                const Snarl* snarl = snarl_manager.manage(snarl1);
                
                // Get its contents0
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents = pb_contents(graph, snarl_manager.shallow_contents(snarl, graph, true));
                
                // We need the right snarl
                REQUIRE(snarl->start().node_id() == 178894);
                REQUIRE(!snarl->start().backward());
                REQUIRE(snarl->end().node_id() == 187208);
                REQUIRE(!snarl->end().backward());
                
                SECTION("The top-level snarl contains all 5 nodes") {
                    REQUIRE(contents.first.size() == 5);
                }
                
                SECTION("The top-level snarl only contains the three edges not in any child snarl") {
                    REQUIRE(contents.second.size() == 3);
                }
                
            }  
        }
        
        TEST_CASE("snarls can be found", "[snarls]") {
    
            // Build a toy graph
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "G"},
                    {"id": 2, "sequence": "A"},
                    {"id": 3, "sequence": "T"},
                    {"id": 4, "sequence": "GGG"},
                    {"id": 5, "sequence": "T"},
                    {"id": 6, "sequence": "A"},
                    {"id": 7, "sequence": "C"},
                    {"id": 8, "sequence": "A"},
                    {"id": 9, "sequence": "A"}
                ],
                "edge": [
                    {"from": 1, "to": 2},
                    {"from": 1, "to": 6},
                    {"from": 2, "to": 3},
                    {"from": 2, "to": 4},
                    {"from": 3, "to": 5},
                    {"from": 4, "to": 5},
                    {"from": 5, "to": 6},
                    {"from": 6, "to": 7},
                    {"from": 6, "to": 8},
                    {"from": 7, "to": 9},
                    {"from": 8, "to": 9}
                    
                ],
                "path": [
                    {"name": "hint", "mapping": [
                        {"position": {"node_id": 1}, "rank" : 1 },
                        {"position": {"node_id": 6}, "rank" : 2 },
                        {"position": {"node_id": 8}, "rank" : 3 },
                        {"position": {"node_id": 9}, "rank" : 4 }
                    ]}
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            // We need to see the path.
            REQUIRE(graph.paths.size() == 1);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();
            
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
                
            SECTION("There are 2 top level snarls") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 2);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                const Snarl* child2 = snarl_manager.top_level_snarls()[1];
                
                if (child1->start().node_id() > child1->end().node_id()) {
                    snarl_manager.flip(child1);
                }
                
                if (child2->start().node_id() > child2->end().node_id()) {
                    snarl_manager.flip(child2);
                }
                
                SECTION("First child is from 1 end to 6 start") {
                    
                    {
                        bool found_in_forward_orientation = (child1->start().node_id() == 1 &&
                                                             child1->start().backward() == false &&
                                                             child1->end().node_id() == 6 &&
                                                             child1->end().backward() == false);
                        bool found_in_reverse_orientation = (child1->start().node_id() == 6 &&
                                                             child1->start().backward() == true &&
                                                             child1->end().node_id() == 1 &&
                                                             child1->end().backward() == true);
                        bool found_snarl = found_in_forward_orientation || found_in_reverse_orientation;
                        REQUIRE(found_snarl);
                    }
                    
                    SECTION("First child has a child from 2 end to 5 start") {
                        REQUIRE(snarl_manager.children_of(child1).size() == 1);
                        
                        const Snarl* subchild = snarl_manager.children_of(child1)[0];
                        
                        {
                            bool found_in_forward_orientation = (subchild->start().node_id() == 2 &&
                                                                 subchild->start().backward() == false &&
                                                                 subchild->end().node_id() == 5 &&
                                                                 subchild->end().backward() == false);
                            bool found_in_reverse_orientation = (subchild->start().node_id() == 5 &&
                                                                 subchild->start().backward() == true &&
                                                                 subchild->end().node_id() == 2 &&
                                                                 subchild->end().backward() == true);
                            bool found_snarl = found_in_forward_orientation || found_in_reverse_orientation;
                            REQUIRE(found_snarl);
                        }
                        
                        SECTION("Subchild has no children") {
                            REQUIRE(snarl_manager.children_of(subchild).size() == 0);
                        }
                            
                    }
                    
                }
                
                SECTION("Second child is from 6 end to 9 start") {
                    {
                        bool found_in_forward_orientation = (child2->start().node_id() == 6 &&
                                                             child2->start().backward() == false &&
                                                             child2->end().node_id() == 9 &&
                                                             child2->end().backward() == false);
                        bool found_in_reverse_orientation = (child2->start().node_id() == 9 &&
                                                             child2->start().backward() == true &&
                                                             child2->end().node_id() == 6 &&
                                                             child2->end().backward() == true);
                        bool found_snarl = found_in_forward_orientation || found_in_reverse_orientation;
                        REQUIRE(found_snarl);
                    }
                    
                    SECTION("Second child has no children") {
                        REQUIRE(snarl_manager.children_of(child2).size() == 0);
                    }
                }
                
            }
                
        }

        TEST_CASE("bubbles can be found in graphs with only heads", "[bubbles]") {
            
            // Build a toy graph
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "G"},
                    {"id": 2, "sequence": "A"}
                ],
                "edge": [
                    {"from": 1, "to": 2, "to_end": true}
                    
                ],
                "path": [
                    {"name": "hint", "mapping": [
                        {"position": {"node_id": 1}, "rank" : 1 },
                        {"position": {"node_id": 2, "is_reverse": true}, "rank" : 2 }
                    ]}
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();

#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
        
            SECTION("Root node has 1 child bubble") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 1);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                
                if (child1->start().node_id() > child1->end().node_id()) {
                    snarl_manager.flip(child1);
                }
                
                SECTION("First child is from 1 end to 2 end") {
                    REQUIRE(child1->start().node_id() == 1);
                    REQUIRE(!child1->start().backward() == true);
                    REQUIRE(child1->end().node_id() == 2);
                    REQUIRE(child1->end().backward() == true);
                    
                    SECTION("First child has no children") {
                        REQUIRE(snarl_manager.children_of(child1).size() == 0);
                    }
                    
                }
            }
            
        }


        TEST_CASE("bubbles can be found in bigger graphs with only heads", "[bubbles][broken]") {
            
            // Build a toy graph
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "G"},
                    {"id": 2, "sequence": "A"},
                    {"id": 3, "sequence": "T"},
                    {"id": 4, "sequence": "GGG"},
                    {"id": 5, "sequence": "T"},
                    {"id": 6, "sequence": "A"},
                    {"id": 7, "sequence": "C"},
                    {"id": 8, "sequence": "A"},
                    {"id": 9, "sequence": "T"}
                ],
                "edge": [
                    {"from": 1, "to": 2},
                    {"from": 1, "to": 6},
                    {"from": 2, "to": 3},
                    {"from": 2, "to": 4},
                    {"from": 3, "to": 5},
                    {"from": 4, "to": 5},
                    {"from": 5, "to": 6},
                    {"from": 6, "to": 7},
                    {"from": 6, "to": 8},
                    {"from": 7, "to": 9, "to_end": true},
                    {"from": 8, "to": 9, "to_end": true}
                    
                ],
                "path": [
                    {"name": "hint", "mapping": [
                        {"position": {"node_id": 1}, "rank" : 1 },
                        {"position": {"node_id": 6}, "rank" : 2 },
                        {"position": {"node_id": 8}, "rank" : 3 },
                        {"position": {"node_id": 9, "is_reverse": true}, "rank" : 4 }
                    ]}
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();
            
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
            SECTION("Root node has 2 child bubbles") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 2);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                const Snarl* child2 = snarl_manager.top_level_snarls()[1];
                
                if (child1->start().node_id() > child1->end().node_id()) {
                    snarl_manager.flip(child1);
                }
                
                if (child2->start().node_id() > child2->end().node_id()) {
                    snarl_manager.flip(child2);
                }
                
                SECTION("First child is from 1 end to 6 start") {
                    REQUIRE(child1->start().node_id() == 1);
                    REQUIRE(!child1->start().backward() == true);
                    REQUIRE(child1->end().node_id() == 6);
                    REQUIRE(child1->end().backward() == false);
                    
                    SECTION("First child has a child from 2 end to 5 start") {
                        REQUIRE(snarl_manager.children_of(child1).size() == 1);
                        
                        const Snarl* subchild = snarl_manager.children_of(child1)[0];
                        
                        if (subchild->start().node_id() > subchild->end().node_id()) {
                            snarl_manager.flip(subchild);
                        }
                        
                        REQUIRE(subchild->start().node_id() == 2);
                        REQUIRE(subchild->start().backward() == false);
                        REQUIRE(subchild->end().node_id() == 5);
                        REQUIRE(subchild->end().backward() == false);
                        
                        SECTION("Subchild has no children") {
                            REQUIRE(snarl_manager.children_of(subchild).size() == 0);
                        }
                            
                    }
                    
                }
                
                SECTION("Second child is from 6 end to 9 end") {
                    REQUIRE(child2->start().node_id() == 6);
                    REQUIRE(!child2->start().backward() == true);
                    REQUIRE(child2->end().node_id() == 9);
                    REQUIRE(child2->end().backward() == true);
                    
                    SECTION("Second child has no children") {
                        REQUIRE(snarl_manager.children_of(child2).size() == 0);
                    }
                }
                
            }
                
        }

        TEST_CASE("bubbles can be found in graphs with only tails", "[bubbles]") {
            
            // Build a toy graph
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "C"},
                    {"id": 2, "sequence": "A"},
                    {"id": 3, "sequence": "T"},
                    {"id": 4, "sequence": "GGG"},
                    {"id": 5, "sequence": "T"},
                    {"id": 6, "sequence": "A"},
                    {"id": 7, "sequence": "C"},
                    {"id": 8, "sequence": "A"},
                    {"id": 9, "sequence": "A"}
                ],
                "edge": [
                    {"from": 1, "to": 2, "from_start": true},
                    {"from": 1, "to": 6, "from_start": true},
                    {"from": 2, "to": 3},
                    {"from": 2, "to": 4},
                    {"from": 3, "to": 5},
                    {"from": 4, "to": 5},
                    {"from": 5, "to": 6},
                    {"from": 6, "to": 7},
                    {"from": 6, "to": 8},
                    {"from": 7, "to": 9},
                    {"from": 8, "to": 9}
                    
                ],
                "path": [
                    {"name": "hint", "mapping": [
                        {"position": {"node_id": 1, "is_reverse": true}, "rank" : 1 },
                        {"position": {"node_id": 6}, "rank" : 2 },
                        {"position": {"node_id": 8}, "rank" : 3 },
                        {"position": {"node_id": 9}, "rank" : 4 }
                    ]}
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();
            
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
                
            SECTION("Root node has 2 child bubbles") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 2);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                const Snarl* child2 = snarl_manager.top_level_snarls()[1];
                
                if (child1->start().node_id() > child1->end().node_id()) {
                    snarl_manager.flip(child1);
                }
                
                if (child2->start().node_id() > child2->end().node_id()) {
                    snarl_manager.flip(child2);
                }
                
                SECTION("First child is from 1 start to 6 start") {
                    REQUIRE(child1->start().node_id() == 1);
                    REQUIRE(!child1->start().backward() == false);
                    REQUIRE(child1->end().node_id() == 6);
                    REQUIRE(child1->end().backward() == false);
                    
                    SECTION("First child should have all the contained nodes in its contents, including contents of its children") {
                        REQUIRE(snarl_manager.deep_contents(child1, graph, true).first.size() == 6);
                    };
                    
                    SECTION("First child has a child from 2 end to 5 start") {
                        REQUIRE(snarl_manager.children_of(child1).size() == 1);
                        
                        const Snarl* subchild = snarl_manager.children_of(child1)[0];
                        
                        if (subchild->start().node_id() > subchild->end().node_id()) {
                            snarl_manager.flip(subchild);
                        }
                        
                        REQUIRE(subchild->start().node_id() == 2);
                        REQUIRE(!subchild->start().backward() == true);
                        REQUIRE(subchild->end().node_id() == 5);
                        REQUIRE(subchild->end().backward() == false);
                        
                        SECTION("Subchild has no children") {
                            REQUIRE(snarl_manager.children_of(subchild).size() == 0);
                        }
                            
                    }
                    
                }
                
                SECTION("Second child is from 6 end to 9 start") {
                    REQUIRE(child2->start().node_id() == 6);
                    REQUIRE(!child2->start().backward() == true);
                    REQUIRE(child2->end().node_id() == 9);
                    REQUIRE(child2->end().backward() == false);
                    
                    SECTION("Second child has no children") {
                        REQUIRE(snarl_manager.children_of(child2).size() == 0);
                    }
                }
                
            }
        }

        TEST_CASE("bubbles can be found when heads cannot reach tails", "[bubbles]") {
            
            // Build a toy graph
            // Looks like:
            //
            //    1\
            //  2---3
            //  v\4 v
            //
            // The head is 1, the tail is 4, 2 and 3 have self loops, and the graph is connected.
            // But the head cannot reach the tail.
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "A"},
                    {"id": 2, "sequence": "A"},
                    {"id": 3, "sequence": "A"},
                    {"id": 4, "sequence": "A"}
                ],
                "edge": [
                    {"from": 1, "to": 3},
                    {"from": 2, "to": 4},
                    {"from": 2, "to": 3},
                    {"from": 2, "to": 2},
                    {"from": 3, "to": 3}            
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();
            
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
            SECTION("Root should have one child actual bubble") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 1);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                
                SECTION("The child should contain all the nodes as contents") {
                    REQUIRE(snarl_manager.deep_contents(child1, graph, true).first.size() == 4);
                }
                
                // TODO: When unary snarls are exposed, make sure we have unary children.
            }
            
        }

        TEST_CASE("bubbles can be found in a graph with no heads or tails", "[bubbles][snarls]") {
            
            // Build a toy graph
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "G"},
                    {"id": 2, "sequence": "A"}
                ],
                "edge": [
                    {"from": 1, "to": 2},
                    {"from": 2, "to": 1}
                    
                ],
                "path": [
                    {"name": "hint", "mapping": [
                        {"position": {"node_id": 1}, "rank" : 1 },
                        {"position": {"node_id": 2}, "rank" : 2 }
                    ]}
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();
            
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
            SECTION("Root node has 2 child bubbles") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 2);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                const Snarl* child2 = snarl_manager.top_level_snarls()[1];
                
                if (child1->start().node_id() > child1->end().node_id()) {
                    snarl_manager.flip(child1);
                }
                
                if (child2->start().node_id() > child2->end().node_id()) {
                    snarl_manager.flip(child2);
                }
                
                SECTION("First child is from 1 start/end to 2 end/start") {
                    REQUIRE(child1->start().node_id() == 1);
                    REQUIRE(child1->end().node_id() == 2);
                    REQUIRE(!child1->start().backward() != child1->end().backward());
                    
                    SECTION("First child has no children") {
                        REQUIRE(snarl_manager.children_of(child1).size() == 0);
                    }
                    
                }
                
                SECTION("Second child is from 1 start/end to 2 end/start") {
                    REQUIRE(child2->start().node_id() == 1);
                    REQUIRE(child2->end().node_id() == 2);
                    REQUIRE(!child1->start().backward() != child1->end().backward());
                    
                    SECTION("Second child has no children") {
                        REQUIRE(snarl_manager.children_of(child2).size() == 0);
                    }
                    
                }
                
            }
            
        }

        TEST_CASE("bubbles can be found in a graph with no heads or tails or paths", "[bubbles]") {
            
            // Build a toy graph
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "G"},
                    {"id": 2, "sequence": "A"}
                ],
                "edge": [
                    {"from": 1, "to": 2},
                    {"from": 2, "to": 1}
                    
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();
            
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
            SECTION("Root node has 2 child bubbles") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 2);
                
                // We should have 2 bubbles, one looking around the cycle in each direction.
                // TODO: can't really say much about its contents.
                
            }
            
        }

        TEST_CASE("bubbles are created based on most distant connected tips", "[bubbles]") {
            
            // Build a toy graph
            // Looks like:
            //
            //       1\
            //  5--2---3--6
            //      \4  
            //
            // The head is 1, the tail is 4, 2 and 3 have self loops, and the graph is connected.
            // But the head cannot reach the tail.
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "A"},
                    {"id": 2, "sequence": "A"},
                    {"id": 3, "sequence": "A"},
                    {"id": 4, "sequence": "A"},
                    {"id": 5, "sequence": "A"},
                    {"id": 6, "sequence": "A"}
                ],
                "edge": [
                    {"from": 1, "to": 3},
                    {"from": 2, "to": 4},
                    {"from": 2, "to": 3},
                    {"from": 5, "to": 2},
                    {"from": 3, "to": 6}
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();
            
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
                
            SECTION("Root should have three child actual bubbles") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 3);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                const Snarl* child2 = snarl_manager.top_level_snarls()[1];
                const Snarl* child3 = snarl_manager.top_level_snarls()[2];
                
                if (child1->start().node_id() > child1->end().node_id()) {
                    snarl_manager.flip(child1);
                }
                
                if (child2->start().node_id() > child2->end().node_id()) {
                    snarl_manager.flip(child2);
                }
                
                if (child3->start().node_id() > child3->end().node_id()) {
                    snarl_manager.flip(child3);
                }
                
                SECTION("First child is trivial snarl from 2 start to 5 end") {
                    // Not 5 end to 2 start because we seem to be sorting by node ID.
                    REQUIRE(child1->start().node_id() == 2);
                    REQUIRE(!child1->start().backward() == false);
                    REQUIRE(child1->end().node_id() == 5);
                    REQUIRE(child1->end().backward() == true);
                    
                }
                
                SECTION("Second child is tip-containing snarl from 2 end to 3 start") {
                    REQUIRE(child2->start().node_id() == 2);
                    REQUIRE(!child2->start().backward() == true);
                    REQUIRE(child2->end().node_id() == 3);
                    REQUIRE(child2->end().backward() == false);
                    
                }
                
                SECTION("Third child is trivial snarl from 3 end to 6 start") {
                    REQUIRE(child3->start().node_id() == 3);
                    REQUIRE(!child3->start().backward() == true);
                    REQUIRE(child3->end().node_id() == 6);
                    REQUIRE(child3->end().backward() == false);
                    
                }
                
            }
            
        }
        
        TEST_CASE("SnarlManager accepts chain input", "[snarls]") {
            // Make a little graph where snarl1 and snarl2 are a top-level
            // chain, and snarl3 and snarl4 are trivial chains inside snarl1
      
            Snarl snarl1;
            snarl1.mutable_start()->set_node_id(1);
            snarl1.mutable_end()->set_node_id(6);
            snarl1.set_start_end_reachable(true);
            
            Snarl snarl2;
            snarl2.mutable_start()->set_node_id(6);
            snarl2.mutable_end()->set_node_id(7);
            snarl2.set_start_end_reachable(true);
            
            Snarl snarl3;
            snarl3.mutable_start()->set_node_id(2);
            snarl3.mutable_end()->set_node_id(3);
            snarl3.set_start_end_reachable(true);
            transfer_boundary_info(snarl1, *snarl3.mutable_parent());
            
            Snarl snarl4;
            snarl4.mutable_start()->set_node_id(4);
            snarl4.mutable_end()->set_node_id(5);
            snarl4.set_start_end_reachable(true);
            transfer_boundary_info(snarl1, *snarl4.mutable_parent());
            
            // Load all this into the SnarlManager
            SnarlManager snarl_manager;
            
            auto ptr3 = snarl_manager.add_snarl(snarl3);
            auto ptr4 = snarl_manager.add_snarl(snarl4);
            
            auto ptr1 = snarl_manager.add_snarl(snarl1);
            
            auto ptr2 = snarl_manager.add_snarl(snarl2);
            
            snarl_manager.finish();
 
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << endl;
            });
#endif

            SECTION("There should be two top-level snarls") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 2);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                const Snarl* child2 = snarl_manager.top_level_snarls()[1];
                
                if (child1->start().node_id() > child1->end().node_id()) {
                    snarl_manager.flip(child1);
                }
                
                if (child2->start().node_id() > child2->end().node_id()) {
                    snarl_manager.flip(child2);
                }
                
                SECTION("First child is from 1 to 6") {
                    REQUIRE(child1->start().node_id() == 1);
                    REQUIRE(child1->start().backward() == false);
                    REQUIRE(child1->end().node_id() == 6);
                    REQUIRE(child1->end().backward() == false);
                    
                    SECTION("First child has two children") {
                        REQUIRE(snarl_manager.children_of(child1).size() == 2);
                        
                        const Snarl* subchild1 = snarl_manager.children_of(child1)[0];
                        const Snarl* subchild2 = snarl_manager.children_of(child1)[1];
                        
                        SECTION("First child is from 2 to 3") {
                            REQUIRE(subchild1->start().node_id() == 2);
                            REQUIRE(subchild1->start().backward() == false);
                            REQUIRE(subchild1->end().node_id() == 3);
                            REQUIRE(subchild1->end().backward() == false);
                        }
                        
                        SECTION("Second child is from 4 to 5") {
                            REQUIRE(subchild2->start().node_id() == 4);
                            REQUIRE(subchild2->start().backward() == false);
                            REQUIRE(subchild2->end().node_id() == 5);
                            REQUIRE(subchild2->end().backward() == false);
                        }
                        
                    }
                    
                }
                
                SECTION("Second child from 6 to 7") {
                    REQUIRE(child2->start().node_id() == 6);
                    REQUIRE(child2->start().backward() == false);
                    REQUIRE(child2->end().node_id() == 7);
                    REQUIRE(child2->end().backward() == false);
                    
                }
                
            }
            
            
        }
        
        TEST_CASE("SnarlManager exposes chains correctly", "[snarls]") {
            
            // We need a graph a chain in it.
            // Looks like:
            //
            //    2   5   8
            //  1   4   7   10 
            //    3   6   9
            //
            const string graph_json = R"(
            
            {
                "node": [
                    {"id": 1, "sequence": "A"},
                    {"id": 2, "sequence": "A"},
                    {"id": 3, "sequence": "A"},
                    {"id": 4, "sequence": "A"},
                    {"id": 5, "sequence": "A"},
                    {"id": 6, "sequence": "A"},
                    {"id": 7, "sequence": "A"},
                    {"id": 8, "sequence": "A"},
                    {"id": 9, "sequence": "A"},
                    {"id": 10, "sequence": "A"}
                ],
                "edge": [
                    {"from": 1, "to": 2},
                    {"from": 1, "to": 3},
                    {"from": 2, "to": 4},
                    {"from": 3, "to": 4},
                    {"from": 4, "to": 5},
                    {"from": 4, "to": 6},
                    {"from": 5, "to": 7},
                    {"from": 6, "to": 7},
                    {"from": 7, "to": 8},
                    {"from": 7, "to": 9},
                    {"from": 8, "to": 10},
                    {"from": 9, "to": 10}
                ]
            }
            
            )";
            
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            
            SnarlManager snarl_manager = CactusSnarlFinder(graph).find_snarls();
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                    << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
                
            SECTION("There should be three top-level snarls") {
                REQUIRE(snarl_manager.top_level_snarls().size() == 3);
                
                const Snarl* child1 = snarl_manager.top_level_snarls()[0];
                const Snarl* child2 = snarl_manager.top_level_snarls()[1];
                const Snarl* child3 = snarl_manager.top_level_snarls()[2];
                
                if (child1->start().node_id() > child1->end().node_id()) {
                    snarl_manager.flip(child1);
                }
                
                if (child2->start().node_id() > child2->end().node_id()) {
                    snarl_manager.flip(child2);
                }
                
                if (child3->start().node_id() > child3->end().node_id()) {
                    snarl_manager.flip(child3);
                }
                
                SECTION("They should be in a chain") {
                    REQUIRE(snarl_manager.in_nontrivial_chain(child1));
                    REQUIRE(snarl_manager.in_nontrivial_chain(child2));
                    REQUIRE(snarl_manager.in_nontrivial_chain(child3));
                }
                
                SECTION("They should be in the same chain") {
                    REQUIRE(snarl_manager.chain_of(child1) == snarl_manager.chain_of(child2));
                    REQUIRE(snarl_manager.chain_of(child2) == snarl_manager.chain_of(child3));
                }
                
                SECTION("We can traverse the chain with iterators") {
                    auto& chains = snarl_manager.chains_of(nullptr);
                    
                    REQUIRE(chains.size() == 1);
                    
                    auto& chain = chains.front();
                    
                    // Make sure we are talking about the same exact chain
                    REQUIRE(&chain == snarl_manager.chain_of(child1));
                    REQUIRE(&chain == snarl_manager.chain_of(child2));
                    REQUIRE(&chain == snarl_manager.chain_of(child3));
                    
                    auto begin = chain_begin(chain);
                    auto rbegin = chain_rbegin(chain);
                    auto rcbegin = chain_rcbegin(chain);
                    auto end = chain_end(chain);
                    auto rend = chain_rend(chain);
                    auto rcend = chain_rcend(chain);
                    
                    SECTION("Iterator equality works") {
                    
                        REQUIRE(begin == begin);
                        REQUIRE(begin != end);
                        REQUIRE(begin != rbegin);
                        REQUIRE(begin != rend);
                        
                        REQUIRE(rbegin == rbegin);
                        REQUIRE(rbegin != end);
                        REQUIRE(rbegin != begin);
                        REQUIRE(rbegin != rend);
                        
                        REQUIRE(end == end);
                        REQUIRE(end != begin);
                        REQUIRE(end != rbegin);
                        REQUIRE(end != rend);
                        
                        REQUIRE(rend == rend);
                        REQUIRE(rend != end);
                        REQUIRE(rend != rbegin);
                        REQUIRE(rend != begin);
                        
                    }
                    
                    SECTION("Iterators traverse the chain left to right and right to left") {
                        auto it = begin;
                        auto rit = rbegin;
                        
                        REQUIRE(*it == make_pair(child1, false));
                        REQUIRE(*rit == make_pair(child3, false));
                        
                        REQUIRE(it->first == child1);
                        REQUIRE(it->second == false);
                        REQUIRE(rit->first == child3);
                        REQUIRE(rit->second == false);
                        
                        ++it;
                        ++rit;
                        
                        REQUIRE(*it == make_pair(child2, false));
                        REQUIRE(*rit == make_pair(child2, false));
                        
                        REQUIRE(it->first == child2);
                        REQUIRE(it->second == false);
                        REQUIRE(rit->first == child2);
                        REQUIRE(rit->second == false);
                        
                        ++it;
                        ++rit;
                        
                        REQUIRE(*it == make_pair(child3, false));
                        REQUIRE(*rit == make_pair(child1, false));
                        
                        REQUIRE(it->first == child3);
                        REQUIRE(it->second == false);
                        REQUIRE(rit->first == child1);
                        REQUIRE(rit->second == false);
                        
                        ++it;
                        ++rit;
                        
                        REQUIRE(it == end);
                        REQUIRE(rit == rend);
                        
                    }
                    
                    SECTION("Reverse complement iterators traverse the chain correctly") {
                        auto it = rcbegin;
                        
                        REQUIRE(*it == make_pair(child3, true));
                        
                        REQUIRE(it->first == child3);
                        REQUIRE(it->second == true);
                        
                        ++it;
                        
                        REQUIRE(*it == make_pair(child2, true));
                        
                        REQUIRE(it->first == child2);
                        REQUIRE(it->second == true);
                        
                        ++it;
                        
                        REQUIRE(*it == make_pair(child1, true));
                        
                        REQUIRE(it->first == child1);
                        REQUIRE(it->second == true);
                        
                        ++it;
                        
                        REQUIRE(it == rcend);
                        
                    }
                    
                    SECTION("We can view the chain from each end") {
#ifdef debug
                        for (auto& pair : chain) {
                            cerr << "In chain " << &chain << ": " << pair.first->start().node_id() << " - " << pair.first->end().node_id()
                                << " orientation " << pair.second << endl;
                        }
#endif
                    
                        REQUIRE(!start_backward(chain));
                        REQUIRE(chain_begin_from(chain, child1, false) == begin);
                        REQUIRE(chain_end_from(chain, child1, false) == end);
                        
                        REQUIRE(!end_backward(chain));
                        REQUIRE(chain_begin_from(chain, child3, true) == rcbegin);
                        REQUIRE(chain_end_from(chain, child3, true) == rcend);
                        
                        snarl_manager.flip(child1);
                        snarl_manager.flip(child3);
                        
#ifdef debug
                        for (auto& pair : chain) {
                            cerr << "In chain " << &chain << ": " << pair.first->start().node_id() << " - " << pair.first->end().node_id()
                                << " orientation " << pair.second << endl;
                        }
#endif
                        
                        REQUIRE(start_backward(chain));
                        REQUIRE(chain_begin_from(chain, child1, true) == chain_begin(chain));
                        REQUIRE(chain_end_from(chain, child1, true) == chain_end(chain));
                        
                        REQUIRE(end_backward(chain));
                        REQUIRE(chain_begin_from(chain, child3, false) == chain_rcbegin(chain));
                        REQUIRE(chain_end_from(chain, child3, false) == chain_rcend(chain));
                        
                    }
                    
                    SECTION("Empty chains have proper iterators") {
                        Chain empty;
                        
                        REQUIRE(chain_begin(empty) == chain_end(empty));
                        REQUIRE(chain_rbegin(empty) == chain_rend(empty));
                    }
                    
                }
                
                SECTION("We can still see the chain if we flip the snarls around") {
                    snarl_manager.flip(child1);
                    snarl_manager.flip(child2);
                    
                    auto& chains = snarl_manager.chains_of(nullptr);
                    
                    REQUIRE(chains.size() == 1);
                    
                    auto& chain = chains.front();
                    
                    REQUIRE(chain.size() == 3);
                    
                    auto it = chain_begin(chain);
                    
                    // Chain should be in the same order but with some orientations flipped.
                    REQUIRE(*it == make_pair(child1, true));
                    ++it;
                    REQUIRE(*it == make_pair(child2, true));
                    ++it;
                    REQUIRE(*it == make_pair(child3, false));
                    ++it;
                    REQUIRE(it == chain_end(chain));
                }
                
                SECTION("We can look around from a snarl") {
                    const Chain* chain = snarl_manager.chain_of(child1);
                    
                    ChainIterator here = chain_begin_from(*chain, child1, false);
                    ChainIterator end = chain_end_from(*chain, child1, false);
                    
                    SECTION("Looking right into the chain gives us the Snarl to the right") {
                        ChainIterator right = here;
                        ++right;
                        
                        REQUIRE(right != end);
                        REQUIRE(right->first == child2);
                        // Must not be backward
                        REQUIRE(right->second == false);
                    }
                    
                    // Now look from the other end
                    here = chain_begin_from(*chain, child3, true);
                    end = chain_end_from(*chain, child3, true);
                    
                    SECTION("Looking left into the chain gives us a Visit to the right Snarl") {
                        ChainIterator left = here;
                        ++left;
                        
                        REQUIRE(left != end);
                        REQUIRE(left->first == child2);
                        // Must be backward
                        REQUIRE(left->second == true);
                    }
                    
                    
                }
                
            }
            
        }
        
        TEST_CASE("Chain start and end functions work on difficult chains", "[snarls]") {
            // This graph will have a snarl from 1 to 8, a snarl from 2 to 4, and a
            // snarl from 4 to 7, with a chain in the top snarl. The snarl from 4 to 7
            // will have a child snarl from 5 to 6 (an insertion of node 9)
            VG graph;
                
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("G");
            Node* n4 = graph.create_node("CTGA");
            Node* n5 = graph.create_node("GCA");
            Node* n6 = graph.create_node("T");
            Node* n7 = graph.create_node("G");
            Node* n8 = graph.create_node("CTGA");
            Node* n9 = graph.create_node("GCA");
            
            Edge* e1 = graph.create_edge(n1, n2);
            Edge* e2 = graph.create_edge(n1, n8);
            Edge* e3 = graph.create_edge(n2, n3);
            Edge* e4 = graph.create_edge(n2, n4);
            Edge* e5 = graph.create_edge(n3, n4);
            Edge* e6 = graph.create_edge(n4, n5);
            Edge* e7 = graph.create_edge(n4, n7);
            Edge* e8 = graph.create_edge(n5, n6);
            Edge* e9 = graph.create_edge(n5, n9);
            Edge* e10 = graph.create_edge(n9, n6);
            Edge* e11 = graph.create_edge(n6, n7);
            Edge* e12 = graph.create_edge(n7, n8);
            
            // Work out its snarls
            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls();
            
#ifdef debug
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                cerr << "Found snarl " << snarl->start().node_id() << " " << snarl->start().backward()
                << " to " << snarl->end().node_id() << " " << snarl->end().backward() << " containing ";
                for (auto& node : snarl_manager.shallow_contents(snarl, graph, false).first) {
                    cerr << node->id() << " ";
                }
                cerr << endl;
            });
#endif
            
            // Get the top snarl
            const Snarl* top_snarl = snarl_manager.top_level_snarls().at(0);
            
            if (top_snarl->end().node_id() < top_snarl->start().node_id()) {
                // Put it a consistent way around
                snarl_manager.flip(top_snarl);
            }
            
            // Make sure it's what we expect.
            REQUIRE(top_snarl->start().node_id() == 1);
            REQUIRE(top_snarl->end().node_id() == 8);

            // Make sure we have one chain    
            auto& chains = snarl_manager.chains_of(top_snarl);
            REQUIRE(chains.size() == 1);
            
            // Get the chain
            auto& chain = chains.at(0);
            REQUIRE(chain.size() == 2);
            
            // And the snarls in the chain
            const Snarl* left_child = chain.at(0).first;
            const Snarl* right_child = chain.at(1).first;
            
            bool found_chain_orientation_1 = false, found_chain_orientation_2 = false;
            
            {
                bool found_left_snarl_orientation_1 = (left_child->start().node_id() == 2 &&
                                                       left_child->start().backward() == false &&
                                                       left_child->end().node_id() == 4 &&
                                                       left_child->end().backward() == false);
                bool found_left_snarl_orientation_2 = (left_child->start().node_id() == 4 &&
                                                       left_child->start().backward() == true &&
                                                       left_child->end().node_id() == 2 &&
                                                       left_child->end().backward() == true);
                bool found_left = found_left_snarl_orientation_1 | found_left_snarl_orientation_2;
                
                bool found_right_snarl_orientation_1 = (right_child->start().node_id() == 4 &&
                                                        right_child->start().backward() == false &&
                                                        right_child->end().node_id() == 7 &&
                                                        right_child->end().backward() == false);
                bool found_right_snarl_orientation_2 = (right_child->start().node_id() == 7 &&
                                                        right_child->start().backward() == true &&
                                                        right_child->end().node_id() == 4 &&
                                                        right_child->end().backward() == true);
                bool found_right = found_right_snarl_orientation_1 | found_right_snarl_orientation_2;
                
                bool found_sub_child = false;
                const auto& children = snarl_manager.children_of(right_child);
                if (children.size() == 1){
                    const Snarl* sub_child = children[0];
                    
                    bool found_sub_child_orientation_1 = (sub_child->start().node_id() == 5 &&
                                                          sub_child->start().backward() == false &&
                                                          sub_child->end().node_id() == 6 &&
                                                          sub_child->end().backward() == false);
                    bool found_sub_child_orientation_2 = (sub_child->start().node_id() == 6 &&
                                                          sub_child->start().backward() == true &&
                                                          sub_child->end().node_id() == 5 &&
                                                          sub_child->end().backward() == true);
                    found_sub_child = found_sub_child_orientation_1 || found_sub_child_orientation_2;
                }
                
                found_chain_orientation_1 = (found_left && found_right && found_sub_child);
            }
            
            {
                bool found_left_snarl_orientation_1 = (right_child->start().node_id() == 2 &&
                                                       right_child->start().backward() == false &&
                                                       right_child->end().node_id() == 4 &&
                                                       right_child->end().backward() == false);
                bool found_left_snarl_orientation_2 = (right_child->start().node_id() == 4 &&
                                                       right_child->start().backward() == true &&
                                                       right_child->end().node_id() == 2 &&
                                                       right_child->end().backward() == true);
                bool found_left = found_left_snarl_orientation_1 | found_left_snarl_orientation_2;
                
                bool found_right_snarl_orientation_1 = (left_child->start().node_id() == 4 &&
                                                        left_child->start().backward() == false &&
                                                        left_child->end().node_id() == 7 &&
                                                        left_child->end().backward() == false);
                bool found_right_snarl_orientation_2 = (left_child->start().node_id() == 7 &&
                                                        left_child->start().backward() == true &&
                                                        left_child->end().node_id() == 4 &&
                                                        left_child->end().backward() == true);
                bool found_right = found_right_snarl_orientation_1 | found_right_snarl_orientation_2;
                
                bool found_sub_child = false;
                const auto& children = snarl_manager.children_of(left_child);
                if (children.size() == 1){
                    const Snarl* sub_child = children[0];
                    
                    bool found_sub_child_orientation_1 = (sub_child->start().node_id() == 5 &&
                                                          sub_child->start().backward() == false &&
                                                          sub_child->end().node_id() == 6 &&
                                                          sub_child->end().backward() == false);
                    bool found_sub_child_orientation_2 = (sub_child->start().node_id() == 6 &&
                                                          sub_child->start().backward() == true &&
                                                          sub_child->end().node_id() == 5 &&
                                                          sub_child->end().backward() == true);
                    found_sub_child = found_sub_child_orientation_1 || found_sub_child_orientation_2;
                }
                
                found_chain_orientation_2 = (found_left && found_right && found_sub_child);
            }
            
            bool chain_correct = found_chain_orientation_1 || found_chain_orientation_2;
            REQUIRE(chain_correct);
            
            // Make sure the right child is BACKWARD
            snarl_manager.flip(right_child);
            
            SECTION("A chain can be found from a backward member snarl") {
                const Chain* chain = snarl_manager.chain_of(right_child);
                
                REQUIRE(chain != nullptr);
                
                SECTION("The chain has the two snarls in it") {
                    REQUIRE(chain->size() == 2);
                    REQUIRE(chain->at(0).first == left_child);
                    REQUIRE(chain->at(1).first == right_child);
                }
                
                SECTION("The chain end orientations are correct") {
                    REQUIRE(start_backward(*chain) == false);
                    REQUIRE(end_backward(*chain) == true);
                    REQUIRE(snarl_manager.chain_orientation_of(left_child) == false);
                    REQUIRE(snarl_manager.chain_orientation_of(right_child) == true);
                }
                
                SECTION("The chain ranks are correct") {
                    REQUIRE(snarl_manager.chain_rank_of(left_child) < snarl_manager.chain_rank_of(right_child));
                }
                    
                SECTION("The chain ends are correct") {
                    REQUIRE(get_start_of(*chain) == left_child->start());
                    REQUIRE(get_end_of(*chain) == reverse(right_child->start()));
                }
                
            }
            
        }

        TEST_CASE( "NetGraph can traverse unary snarls",
                  "[snarls][netgraph]" ) {
        
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
            
            // Define the snarls for the top level
           
            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls();

            const Snarl* snarl = snarl_manager.into_which_snarl(3, true); 
            
            
            SECTION( "Traverse unary snarl in netgraph with internal connectivity" ) {

                // Make a net graph
                NetGraph net_graph(snarl->start(), snarl->end(), snarl_manager.chains_of(snarl), &graph, true);
                handle_t node_in = net_graph.get_handle(1, true);
                handle_t node_out = net_graph.get_handle(1, false);
            
                unordered_set<pair<id_t, bool>> seen_in;
                net_graph.follow_edges(node_in, false, [&](const handle_t& other) {
                    seen_in.insert(make_pair(net_graph.get_id(other), 
                                             net_graph.get_is_reverse(other)));
                });
                
                REQUIRE(seen_in.size() == 2);
                REQUIRE(seen_in.find(make_pair(2, false)) != seen_in.end());
                REQUIRE(seen_in.find(make_pair(3, false)) != seen_in.end());

                unordered_set<pair<id_t, bool>> seen_out;
                net_graph.follow_edges(node_out, false, [&](const handle_t& other) {
                    seen_out.insert(make_pair(net_graph.get_id(other), 
                                             net_graph.get_is_reverse(other)));
                });

                unordered_set<pair<id_t, bool>> seen_out_rev;
                net_graph.follow_edges(node_out, true, [&](const handle_t& other) {
                    seen_out_rev.insert(make_pair(net_graph.get_id(other), 
                                             net_graph.get_is_reverse(other)));
});
                //Predecessors
                REQUIRE(seen_out_rev.size() == 2);
                REQUIRE(seen_out_rev.find(make_pair(2, true)) != seen_out_rev.end());
                REQUIRE(seen_out_rev.find(make_pair(3, true)) != seen_out_rev.end());
                
                unordered_set<pair<id_t, bool>> seen_in_rev;
                net_graph.follow_edges(node_in, true, [&](const handle_t& other) {
                    seen_in_rev.insert(make_pair(net_graph.get_id(other), 
                                             net_graph.get_is_reverse(other)));
                });
                REQUIRE(seen_in_rev.size() == 0);
            }

            
            
             SECTION( "Traverse unary snarl in netgraph without internal connectivity" ) {
                
                NetGraph net_graph(snarl->start(), snarl->end(), snarl_manager.chains_of(snarl), &graph);
                handle_t node_in = net_graph.get_handle(1, true);
                handle_t node_out = net_graph.get_handle(1, false);
                unordered_set<pair<id_t, bool>> seen_in;
                net_graph.follow_edges(node_in, false, [&](const handle_t& other) {
                    seen_in.insert(make_pair(net_graph.get_id(other), 
                                             net_graph.get_is_reverse(other)));
                });
                
                REQUIRE(seen_in.size() == 0);

                unordered_set<pair<id_t, bool>> seen_out;
                net_graph.follow_edges(node_out, false, [&](const handle_t& other) {
                    seen_out.insert(make_pair(net_graph.get_id(other), 
                                             net_graph.get_is_reverse(other)));
                });
                REQUIRE(seen_out.size() == 2);
                REQUIRE(seen_out.find(make_pair(2, false)) != seen_out.end());
                REQUIRE(seen_out.find(make_pair(3, false)) != seen_out.end());

                //Predecessors
                unordered_set<pair<id_t, bool>> seen_out_rev;
                net_graph.follow_edges(node_out, true, [&](const handle_t& other) {
                    seen_out_rev.insert(make_pair(net_graph.get_id(other), 
                                             net_graph.get_is_reverse(other)));
                });
                
                REQUIRE(seen_out_rev.size() == 0);
                
                unordered_set<pair<id_t, bool>> seen_in_rev;
                net_graph.follow_edges(node_in, true, [&](const handle_t& other) {
                    seen_in_rev.insert(make_pair(net_graph.get_id(other), 
                                             net_graph.get_is_reverse(other)));
                });
                REQUIRE(seen_in_rev.size() == 2);
                REQUIRE(seen_in_rev.find(make_pair(2, true)) != seen_in_rev.end());
                REQUIRE(seen_in_rev.find(make_pair(3, true)) != seen_in_rev.end());
             }
        }

        TEST_CASE( "Snarls can be found for a graph with no ordinary cycles", "[snarls]" ) {
            VG graph;
                
            // We have this dumbell-shaped graph, where you have to break open
            // a cycle but just saying you go from a node to itself isn't a
            // valid snarl.
            
            Node* n1 = graph.create_node("A");
            Node* n2 = graph.create_node("G");
            
            Edge* e1 = graph.create_edge(n1, n2);
            Edge* e2 = graph.create_edge(n1, n1, true, false);
            Edge* e3 = graph.create_edge(n2, n2, false, true);
            
            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls();
            
            // There must be something in the top level snarls
            REQUIRE(!snarl_manager.top_level_snarls().empty());
            
            // The decomposition must cover all the nodes
            unordered_set<id_t> seen_nodes;
            
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                // Get the contents of each snarl
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents = pb_contents(graph, snarl_manager.shallow_contents(snarl, graph, true));
            
                for (auto& node_ptr : contents.first) {
                    // And record all the nodes
                    seen_nodes.insert(node_ptr->id());
                }
            });
            
            // Make sure both nodes appear.
            REQUIRE(seen_nodes.size() == 2);
            
            
        }
        
        TEST_CASE( "Snarls can be found for a bigger graph with no ordinary cycles", "[snarls]" ) {
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
            
            CactusSnarlFinder bubble_finder(graph);
            SnarlManager snarl_manager = bubble_finder.find_snarls();
            
            // There must be something in the top level snarls
            REQUIRE(!snarl_manager.top_level_snarls().empty());
            
            // The decomposition must cover all the nodes
            unordered_set<id_t> seen_nodes;
            
            snarl_manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                // Get the contents of each snarl
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents = pb_contents(graph, snarl_manager.shallow_contents(snarl, graph, true));
            
                for (auto& node_ptr : contents.first) {
                    // And record all the nodes
                    seen_nodes.insert(node_ptr->id());
                }
            });
            
            // Make sure all nodes appear.
            REQUIRE(seen_nodes.size() == 5);
            
            
        }

        TEST_CASE( "NetGraph can traverse looping snarls",
                  "[snarls][netgraph]" ) {
        
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
            
            // This test depends on the snarl decomposition being rooted at
            // node 1. So we force that to happen by specifying the snarls
            // manually.
            
            vector<Snarl> to_manage;
            
            // We only need the top snarl and its direct child for this test.
            to_manage.emplace_back();
            to_manage.back().mutable_start()->set_node_id(1);
            to_manage.back().mutable_start()->set_backward(true);
            to_manage.back().mutable_end()->set_node_id(1);
            to_manage.back().mutable_end()->set_backward(true);
            to_manage.back().set_type(SnarlType::UNCLASSIFIED);
            to_manage.back().set_start_self_reachable(true);
            to_manage.back().set_end_self_reachable(true);
            to_manage.back().set_start_end_reachable(true);
            
            
            to_manage.emplace_back();
            to_manage.back().mutable_start()->set_node_id(2);
            to_manage.back().mutable_start()->set_backward(false);
            to_manage.back().mutable_end()->set_node_id(2);
            to_manage.back().mutable_end()->set_backward(true);
            to_manage.back().set_type(SnarlType::UNARY);
            to_manage.back().set_start_self_reachable(true);
            to_manage.back().set_end_self_reachable(true);
            to_manage.back().set_start_end_reachable(true);
            to_manage.back().mutable_parent()->mutable_start()->set_node_id(1);
            to_manage.back().mutable_parent()->mutable_start()->set_backward(true);
            to_manage.back().mutable_parent()->mutable_end()->set_node_id(1);
            to_manage.back().mutable_parent()->mutable_end()->set_backward(true);
            
           
            SnarlManager snarl_manager(to_manage.begin(), to_manage.end());

            const vector<const Snarl*>& snarls = snarl_manager.top_level_snarls();
            
            const Snarl* topSnarl = snarls[0];
            //top level snarl starting and ending at node 1 in reverse - 
            //Not a unary snarl but starts and ends at opposite sides of the same node  
            
            
            SECTION( "Traverse unary snarl in netgraph with internal connectivity" ) {

                /* Make a net graph - contains two nodes - node1 and the unary 
                   snarl at node2. node1 has two edges to the beginning of node2
                   Top snarl starts at node1 pointing left and ends at node1 
                   pointing left
                */ 

                NetGraph ng(topSnarl->start(), topSnarl->end(), snarl_manager.chains_of(topSnarl), &graph, true);
  
                handle_t startHandle = ng.get_handle(topSnarl->start().node_id()
                          , topSnarl->start().backward());
                handle_t startHandleR = ng.get_handle(
                    topSnarl->start().node_id(), !topSnarl->start().backward());


                const vector<const Snarl*>& children= snarl_manager.children_of(
                                                                   topSnarl);
                //One child snarl starting at node 2 forward
                const Snarl* childSnarl = children.at(0);
                pair<id_t, bool> childNode (childSnarl->start().node_id(),
                                            childSnarl->start().backward());

 
                pair<id_t, bool> nextFd;
                auto getNextFd = [&](const handle_t& h) -> bool {
                    nextFd.first =  ng.get_id(h); 
                    nextFd.second = ng.get_is_reverse(h);
                    return true;
                };          

                ng.follow_edges(startHandle, false, getNextFd);
              
                //Following edges going forward from start will find child node
                REQUIRE(nextFd == childNode);

                pair<id_t, bool> nextRev;
                auto getNextRev = [&](const handle_t& h) -> bool {
                    nextRev.first =  ng.get_id(h);
                    nextRev.second = ng.get_is_reverse(h);
                    return true;
                };          

                //Following edges going backward from start will find child node
                ng.follow_edges(startHandleR, false, getNextRev);
                REQUIRE(nextRev == childNode);

 
                unordered_set<id_t> allHandles;
                auto addHandle = [&](const handle_t& h) -> bool {
                    allHandles.insert(ng.get_id(h));
                    return true;
                };          
                ng.for_each_handle(addHandle);                   


                //for_each_handle finds both nodes in net graph
                REQUIRE(allHandles.size() == 2);
                REQUIRE(allHandles.count(childSnarl->start().node_id()) == 1);


                //Check following edges from child snarl
                handle_t childFd = ng.get_handle(childNode.first, childNode.second);
                handle_t childRev = ng.get_handle(childNode.first, !childNode.second);

                //Following edges going fd from child node will find start
                unordered_set<pair<id_t, bool>> seenFd;

                auto childNextFd = [&](const handle_t& h) -> bool {
                    seenFd.insert(make_pair( ng.get_id(h), 
                                             ng.get_is_reverse(h)));
                    return true;
                };          
                ng.follow_edges(childFd, false, childNextFd);
                REQUIRE(seenFd.size() == 2);
                REQUIRE(seenFd.count(make_pair(topSnarl->start().node_id(),
                                           topSnarl->start().backward())) == 1);
                REQUIRE(seenFd.count(make_pair(topSnarl->start().node_id(),
                                          !topSnarl->start().backward())) == 1);
                
                //Following edges going back from child node will find nothing
                unordered_set<pair<id_t, bool>> seenRev;

                auto childNextRev = [&](const handle_t& h) -> bool {
                    seenRev.insert(make_pair( ng.get_id(h), 
                                             ng.get_is_reverse(h)));
                    return true;
                };          
                ng.follow_edges(childRev, false, childNextRev);
                REQUIRE(seenRev.size() == 0);
                
    
                //Following edges fd from child going left (predecessors) finds nothing 
                unordered_set<pair<id_t, bool>> seenFdPred;
                ng.follow_edges(childFd, true, [&](const handle_t& other) {
                    seenFdPred.insert(make_pair(ng.get_id(other), 
                                             ng.get_is_reverse(other)));
                });
                REQUIRE(seenFdPred.size() == 0);

                //Following edges rev from child going left (predecessors) finds nothing 
                unordered_set<pair<id_t, bool>> seenRevPred;
                ng.follow_edges(childRev, true, [&](const handle_t& other) {
                    seenRevPred.insert(make_pair(ng.get_id(other), 
                                             ng.get_is_reverse(other)));
                });
                REQUIRE(seenRevPred.size() == 2);
                REQUIRE(seenRevPred.count(make_pair(topSnarl->start().node_id(),
                                           topSnarl->start().backward())) == 1);
                REQUIRE(seenRevPred.count(make_pair(topSnarl->start().node_id(),
                                          !topSnarl->start().backward())) == 1);

            }
            SECTION( "Traverse unary snarl in netgraph without internal connectivity" ) {

                /* Make a net graph - contains two nodes - node1 and the unary 
                   snarl at node2. node1 has two edges to the beginning of node2
                   Top snarl starts at node1 pointing left and ends at node1 
                   pointing left
                */ 

                NetGraph ng(topSnarl->start(), topSnarl->end(), snarl_manager.chains_of(topSnarl), &graph);
  
                handle_t startHandle = ng.get_handle(topSnarl->start().node_id()
                          , topSnarl->start().backward());
                handle_t startHandleR = ng.get_handle(
                    topSnarl->start().node_id(), !topSnarl->start().backward());


                const vector<const Snarl*>& children= snarl_manager.children_of(
                                                                   topSnarl);
                //One child snarl starting at node 2 forward
                const Snarl* childSnarl = children.at(0);
                pair<id_t, bool> childNode (childSnarl->start().node_id(),
                                            childSnarl->start().backward());

 
                pair<id_t, bool> nextFd;
                auto getNextFd = [&](const handle_t& h) -> bool {
                    nextFd.first =  ng.get_id(h); 
                    nextFd.second = ng.get_is_reverse(h);
                    return true;
                };          

                ng.follow_edges(startHandle, false, getNextFd);
              
                //Following edges going forward from start will find child node
                REQUIRE(nextFd == childNode);

                pair<id_t, bool> nextRev;
                auto getNextRev = [&](const handle_t& h) -> bool {
                    nextRev.first =  ng.get_id(h);
                    nextRev.second = ng.get_is_reverse(h);
                    return true;
                };          

                //Following edges going backward from start will find child node
                ng.follow_edges(startHandleR, false, getNextRev);
                REQUIRE(nextRev == childNode);

 
                unordered_set<id_t> allHandles;
                auto addHandle = [&](const handle_t& h) -> bool {
                    allHandles.insert(ng.get_id(h));
                    return true;
                };          
                ng.for_each_handle(addHandle);                   


                //for_each_handle finds both nodes in net graph
                REQUIRE(allHandles.size() == 2);
                REQUIRE(allHandles.count(childSnarl->start().node_id()) == 1);


                //Check following edges from child snarl
                handle_t childFd = ng.get_handle(childNode.first, childNode.second);
                handle_t childRev = ng.get_handle(childNode.first, !childNode.second);


                //Following edges going fd from child node will find start
                unordered_set<pair<id_t, bool>> seenFd;

                auto childNextFd = [&](const handle_t& h) -> bool {
                    seenFd.insert(make_pair( ng.get_id(h), 
                                             ng.get_is_reverse(h)));
                    return true;
                };          
                ng.follow_edges(childFd, false, childNextFd);
                REQUIRE(seenFd.size() == 0);
                
                //Following edges going back from child node will find start
                unordered_set<pair<id_t, bool>> seenRev;

                auto childNextRev = [&](const handle_t& h) -> bool {
                    seenRev.insert(make_pair( ng.get_id(h), 
                                             ng.get_is_reverse(h)));
                    return true;
                };          
                ng.follow_edges(childRev, false, childNextRev);
                REQUIRE(seenRev.size() == 2);
                REQUIRE(seenRev.count(make_pair(topSnarl->start().node_id(),
                                           topSnarl->start().backward())) == 1);
                REQUIRE(seenRev.count(make_pair(topSnarl->start().node_id(),
                                          !topSnarl->start().backward())) == 1);
    
                //Following edges fd from child going left (predecessors) finds start 
                unordered_set<pair<id_t, bool>> seenFdPred;
                ng.follow_edges(childFd, true, [&](const handle_t& other) {
                    seenFdPred.insert(make_pair(ng.get_id(other), 
                                             ng.get_is_reverse(other)));
                });
                REQUIRE(seenFdPred.size() == 2);
                REQUIRE(seenFdPred.count(make_pair(topSnarl->start().node_id(),
                                           topSnarl->start().backward())) == 1);
                REQUIRE(seenFdPred.count(make_pair(topSnarl->start().node_id(),
                                          !topSnarl->start().backward())) == 1);



                //Following edges rev from child going left (predecessors) finds nothing 
                unordered_set<pair<id_t, bool>> seenRevPred;
                ng.follow_edges(childRev, true, [&](const handle_t& other) {
                    seenRevPred.insert(make_pair(ng.get_id(other), 
                                             ng.get_is_reverse(other)));
                });
                REQUIRE(seenRevPred.size() == 0);
            }

            
            
             
        }
        
        TEST_CASE( "NetGraph can traverse all the snarls in random graphs",
                  "[snarls][netgraph][integrated-snarl-finder]" ) {
        
            // Each actual graph takes a fairly long time to do so we randomize sizes...
            
            default_random_engine generator(test_seed_source());
            
            for (size_t repeat = 0; repeat < 100; repeat++) {
            
                uniform_int_distribution<size_t> bases_dist(100, 10000);
                size_t bases = bases_dist(generator);
                uniform_int_distribution<size_t> variant_bases_dist(1, bases/2);
                size_t variant_bases = variant_bases_dist(generator);
                uniform_int_distribution<size_t> variant_count_dist(1, bases/2);
                size_t variant_count = variant_count_dist(generator);
                        
#ifdef debug
                cerr << repeat << ": Do graph of " << bases << " bp with ~" << variant_bases << " bp large variant length and " << variant_count << " events" << endl;
#endif
            
                VG graph;
                random_graph(bases, variant_bases, variant_count, &graph);
                IntegratedSnarlFinder finder(graph); 
                auto manager = finder.find_snarls_parallel();
                
                size_t snarls_seen = 0;
                
                manager.for_each_snarl_preorder([&](const Snarl* snarl) {
                    snarls_seen++;
                    // Make sure we don't follow internal-snarl edges
                    NetGraph net_graph = manager.net_graph_of(snarl, &graph, false);
                    
                    // Make sure we like the nodes in this snarl
                    unordered_set<handle_t> nodes;
                    size_t node_count = 0;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                        node_count++;
                    });
                    if (nodes.size() != node_count) {
                        cerr << "Problem with graph: " << pb2json(graph.graph) << endl;
                    }
                    REQUIRE(nodes.size() == node_count);
                    
                    // Track edges as followed from, followed to
                    unordered_set<pair<handle_t, handle_t>> edges;
                    // And in canonical order
                    unordered_set<pair<handle_t, handle_t>> canonical_edges;
                    for (auto& handle : nodes) {
                    
                        // Save all the edges off of each node
                        net_graph.follow_edges(handle, false, [&](const handle_t& other) {
                            edges.emplace(handle, other);
                            canonical_edges.insert(net_graph.edge_handle(handle, other));
                        });
                        net_graph.follow_edges(handle, true, [&](const handle_t& other) {
                            edges.emplace(net_graph.flip(handle), net_graph.flip(other));
                            canonical_edges.insert(net_graph.edge_handle(other, handle));
                        });
                    }
                    if (canonical_edges.size() != net_graph.get_edge_count()) {
                        cerr << "Problem with graph: " << pb2json(graph.graph) << endl;
                        cerr << "In snarl: " << snarl->start().node_id() << " -> " << snarl->end().node_id() << endl;
                        cerr << "Observed edges: " << endl;
                        for (auto& edge : canonical_edges) {
                            cerr << "\t" << net_graph.get_id(edge.first) << (net_graph.get_is_reverse(edge.first) ? "-" : "+")
                                << " -> " << net_graph.get_id(edge.second) << (net_graph.get_is_reverse(edge.second) ? "-" : "+") << endl;
                        }
                    }
                    REQUIRE(canonical_edges.size() == net_graph.get_edge_count());
                    for (auto& edge : edges) {
                        // Each edge should be seen from both directions
                        REQUIRE(edges.count(make_pair(net_graph.flip(edge.second), net_graph.flip(edge.first))));
                    }
                });
                    
            }
        
            
        }
        
        TEST_CASE( "SnarlManager IO works correctly",
                  "[sites][snarls]" ) {
            
            SECTION( "SnarlManager can be loaded from an empty snarls file") {
                
                stringstream buff;
                
                {
                    // Make an emitter
                    vg::io::ProtobufEmitter<Snarl> emitter(buff);
                    // Emit nothing
                }
                
                // Now load it back
                unique_ptr<SnarlManager> manager = vg::io::VPKG::try_load_one<SnarlManager>(buff);
                
                REQUIRE(manager.get() != nullptr);
                
            }
            
        }

        TEST_CASE( "PathTraversalFinder works correctly",
               "[sites][snarls]" ) {


            // Build a toy graph
            const string graph_json = R"(
            
                {
                    "node": [
                        {"id": 1, "sequence": "G"},
                        {"id": 2, "sequence": "A"},
                        {"id": 3, "sequence": "G"},
                        {"id": 4, "sequence": "A"},
                        {"id": 5, "sequence": "G"},
                        {"id": 6, "sequence": "A"},
                        {"id": 7, "sequence": "G"},
                        {"id": 8, "sequence": "A"},
                        {"id": 9, "sequence": "A"},
                        {"id": 10, "sequence": "G"},
                        {"id": 11, "sequence": "A"}
                    ],
                    "edge": [
                        {"from": 1, "to": 2},
                        {"from": 2, "to": 3},
                        {"from": 2, "to": 4},
                        {"from": 3, "to": 5},
                        {"from": 3, "to": 6},
                        {"from": 4, "to": 5},
                        {"from": 4, "to": 6},
                        {"from": 5, "to": 7},
                        {"from": 6, "to": 7},
                        {"from": 7, "to": 8},
                        {"from": 8, "to": 9},
                        {"from": 8, "to": 10, "to_end": "true"},
                        {"from": 9, "to": 10},
                        {"from": 9, "to": 10, "from_start": "true"},
                        {"from": 9, "to": 11, "to_end": "true"},
                        {"from": 9, "to": 11, "from_start": "true", "to_end": "true"},
                        {"from": 10, "to": 11, "to_end": "true"}
                    ],
                    "path": [
                        {"name": "ref", "mapping": [
                            {"position": {"node_id": 1}, "rank" : 1 },
                            {"position": {"node_id": 2}, "rank" : 2 },
                            {"position": {"node_id": 4}, "rank" : 3 },
                            {"position": {"node_id": 6}, "rank" : 4 },
                            {"position": {"node_id": 7}, "rank" : 5 },
                            {"position": {"node_id": 8}, "rank" : 6 },
                            {"position": {"node_id": 9}, "rank" : 7 },
                            {"position": {"node_id": 10}, "rank" : 8 },
                            {"position": {"node_id": 11, "is_reverse" : "true"}, "rank" : 9 }
                        ]},
                        {"name": "alt1", "mapping": [
                            {"position": {"node_id": 1}, "rank" : 1 },
                            {"position": {"node_id": 2}, "rank" : 2 },
                            {"position": {"node_id": 3}, "rank" : 3 },
                            {"position": {"node_id": 5}, "rank" : 4 },
                            {"position": {"node_id": 7}, "rank" : 5 },
                            {"position": {"node_id": 8}, "rank" : 6 },
                            {"position": {"node_id": 10, "is_reverse" : "true"}, "rank" : 7 },
                            {"position": {"node_id": 9}, "rank" : 8 },
                            {"position": {"node_id": 11, "is_reverse" : "true"}, "rank" : 9 }
                        ]},
                        {"name": "alt1a", "mapping": [
                            {"position": {"node_id": 2}, "rank" : 2 },
                            {"position": {"node_id": 3}, "rank" : 3 },
                            {"position": {"node_id": 5}, "rank" : 4 },
                            {"position": {"node_id": 7}, "rank" : 5 }
                        ]},
                        {"name": "alt2", "mapping": [
                            {"position": {"node_id": 8, "is_reverse" : "true"}, "rank" : 1 },
                            {"position": {"node_id": 7, "is_reverse" : "true"}, "rank" : 2 },
                            {"position": {"node_id": 6, "is_reverse" : "true"}, "rank" : 3 },
                            {"position": {"node_id": 3, "is_reverse" : "true"}, "rank" : 4 },
                            {"position": {"node_id": 2, "is_reverse" : "true"}, "rank" : 5 },
                            {"position": {"node_id": 1, "is_reverse" : "true"}, "rank" : 6 }
                        ]},
                        {"name": "shorty", "mapping": [
                            {"position": {"node_id": 1}, "rank" : 1 },
                            {"position": {"node_id": 2}, "rank" : 2 },
                            {"position": {"node_id": 3}, "rank" : 3 },
                            {"position": {"node_id": 6}, "rank" : 4 }
                        ]},
                        {"name": "alt3", "mapping": [
                            {"position": {"node_id": 11}, "rank" : 1 },
                            {"position": {"node_id": 9}, "rank" : 2 },
                            {"position": {"node_id": 10}, "rank" : 3 },
                            {"position": {"node_id": 8, "is_reverse" : "true"}, "rank" : 4 },
                            {"position": {"node_id": 7, "is_reverse" : "true"}, "rank" : 5 }
                        ]},
                        {"name": "alt4", "mapping": [
                            {"position": {"node_id": 11}, "rank" : 1 },
                            {"position": {"node_id": 10, "is_reverse" : "true"}, "rank" : 2 },
                            {"position": {"node_id": 9, "is_reverse" : "true"}, "rank" : 3 },
                            {"position": {"node_id": 8, "is_reverse" : "true"}, "rank" : 4 },
                            {"position": {"node_id": 7, "is_reverse" : "true"}, "rank" : 5 }
                        ]}
                    ]
                }            
                )";
                
            // Make an actual graph
            VG graph;
            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            graph.extend(chunk);
            assert(graph.is_valid());
            
            SECTION( "PathTraversalFinder can find simple forward traversals") {

                CactusSnarlFinder snarl_finder(graph);
                SnarlManager snarl_manager = snarl_finder.find_snarls();
                PathTraversalFinder trav_finder(graph, snarl_manager);

                Snarl snarl;
                snarl.mutable_start()->set_node_id(2);
                snarl.mutable_end()->set_node_id(7);

                auto trav_results = trav_finder.find_path_traversals(snarl);

                // get a path for ref, atl1, alt1a and alt2
                REQUIRE(trav_results.first.size() == 4);

                set<string> correct_names = {"ref", "alt1", "alt1a", "alt2"};
                for (auto step_pair : trav_results.second) {
                    string name1 = graph.get_path_name(graph.get_path_handle_of_step(step_pair.first));
                    string name2 = graph.get_path_name(graph.get_path_handle_of_step(step_pair.second));
                    REQUIRE(name1 == name2);
                    REQUIRE(correct_names.count(name1));
                }
                
                map<string, string> true_trav_strings = {
                    {"ref", R"({"visit":[{"node_id":"2"},{"node_id":"4"},{"node_id":"6"},{"node_id":"7"}]})"},
                    {"alt1", R"({"visit":[{"node_id":"2"},{"node_id":"3"},{"node_id":"5"},{"node_id":"7"}]})"},
                    {"alt1a", R"({"visit":[{"node_id":"2"},{"node_id":"3"},{"node_id":"5"},{"node_id":"7"}]})"},
                    {"alt2", R"({"visit":[{"node_id":"2"},{"node_id":"3"},{"node_id":"6"},{"node_id":"7"}]})"}
                };
                for (int i = 0; i < trav_results.first.size(); ++i) {
                    string name = graph.get_path_name(graph.get_path_handle_of_step(trav_results.second[i].first));
                    REQUIRE(true_trav_strings.count(name));
                    SnarlTraversal true_trav;
                    json2pb(true_trav, true_trav_strings[name]);
                    bool trav_is_correct = trav_results.first[i] == true_trav;
                    REQUIRE(trav_is_correct);
                }
            }

            SECTION( "PathTraversalFinder can find simple traversals when snarl is backward") {

                CactusSnarlFinder snarl_finder(graph);
                SnarlManager snarl_manager = snarl_finder.find_snarls();
                PathTraversalFinder trav_finder(graph, snarl_manager);

                Snarl snarl;
                snarl.mutable_start()->set_node_id(7);
                snarl.mutable_start()->set_backward(true);
                snarl.mutable_end()->set_node_id(2);
                snarl.mutable_end()->set_backward(true);
                
                auto trav_results = trav_finder.find_path_traversals(snarl);

                // get a path for ref, atl1, alt1a and alt2
                REQUIRE(trav_results.first.size() == 4);

                set<string> correct_names = {"ref", "alt1", "alt1a", "alt2"};
                for (auto step_pair : trav_results.second) {
                    string name1 = graph.get_path_name(graph.get_path_handle_of_step(step_pair.first));
                    string name2 = graph.get_path_name(graph.get_path_handle_of_step(step_pair.second));
                    REQUIRE(name1 == name2);
                    REQUIRE(correct_names.count(name1));
                }

                map<string, string> true_trav_strings = {
                    {"ref", R"({"visit":[{"node_id":"7","backward":true},{"node_id":"6","backward":true},{"node_id":"4","backward":true},{"node_id":"2","backward":true}]})"},
                    {"alt1", R"({"visit":[{"node_id":"7","backward":true},{"node_id":"5","backward":true},{"node_id":"3","backward":true},{"node_id":"2","backward":true}]})"},
                    {"alt1a", R"({"visit":[{"node_id":"7","backward":true},{"node_id":"5","backward":true},{"node_id":"3","backward":true},{"node_id":"2","backward":true}]})"},
                    {"alt2", R"({"visit":[{"node_id":"7","backward":true},{"node_id":"6","backward":true},{"node_id":"3","backward":true},{"node_id":"2","backward":true}]})"}
                };
                for (int i = 0; i < trav_results.first.size(); ++i) {
                    string name = graph.get_path_name(graph.get_path_handle_of_step(trav_results.second[i].first));                
                    REQUIRE(true_trav_strings.count(name));
                    SnarlTraversal true_trav;
                    json2pb(true_trav, true_trav_strings[name]);
                    bool trav_is_correct = trav_results.first[i] == true_trav;
                    REQUIRE(trav_is_correct);
                }
            }

            SECTION( "PathTraversalFinder can find forward traversals in snarl with inversion") {

                CactusSnarlFinder snarl_finder(graph);
                SnarlManager snarl_manager = snarl_finder.find_snarls();
                PathTraversalFinder trav_finder(graph, snarl_manager);

                Snarl snarl;
                snarl.mutable_start()->set_node_id(8);
                snarl.mutable_end()->set_node_id(11);
                snarl.mutable_end()->set_backward(true);
                
                auto trav_results = trav_finder.find_path_traversals(snarl);

                // get a path for ref, atl1, alt1a and alt2
                REQUIRE(trav_results.first.size() == 4);

                set<string> correct_names = {"ref", "alt1", "alt3", "alt4"};
                for (auto step_pair : trav_results.second) {
                    string name1 = graph.get_path_name(graph.get_path_handle_of_step(step_pair.first));
                    string name2 = graph.get_path_name(graph.get_path_handle_of_step(step_pair.second));
                    REQUIRE(name1 == name2);
                    REQUIRE(correct_names.count(name1));
                }

                map<string, string> true_trav_strings = {
                    {"ref", R"({"visit":[{"node_id":"8"},{"node_id":"9"},{"node_id":"10"},{"node_id":"11","backward":true}]})"},
                    {"alt1", R"({"visit":[{"node_id":"8"},{"node_id":"10","backward":true},{"node_id":"9"},{"node_id":"11","backward":true}]})"},
                    {"alt3", R"({"visit":[{"node_id":"8"},{"node_id":"10","backward":true},{"node_id":"9","backward":true},{"node_id":"11","backward":true}]})"},
                    {"alt4", R"({"visit":[{"node_id":"8"},{"node_id":"9"},{"node_id":"10"},{"node_id":"11","backward":true}]})"}
                };
                for (int i = 0; i < trav_results.first.size(); ++i) {
                    string name = graph.get_path_name(graph.get_path_handle_of_step(trav_results.second[i].first));                
                    REQUIRE(true_trav_strings.count(name));
                    SnarlTraversal true_trav;
                    json2pb(true_trav, true_trav_strings[name]);
                    bool trav_is_correct = trav_results.first[i] == true_trav;
                    REQUIRE(trav_is_correct);
                }
            }

            SECTION( "PathTraversalFinder can find traversals in backward snarl with inversion") {

                CactusSnarlFinder snarl_finder(graph);
                SnarlManager snarl_manager = snarl_finder.find_snarls();
                PathTraversalFinder trav_finder(graph, snarl_manager);

                Snarl snarl;
                snarl.mutable_start()->set_node_id(11);
                snarl.mutable_end()->set_node_id(8);
                snarl.mutable_end()->set_backward(true);
                
                auto trav_results = trav_finder.find_path_traversals(snarl);

                // get a path for ref, atl1, alt1a and alt2
                REQUIRE(trav_results.first.size() == 4);

                set<string> correct_names = {"ref", "alt1", "alt3", "alt4"};
                for (auto step_pair : trav_results.second) {
                    string name1 = graph.get_path_name(graph.get_path_handle_of_step(step_pair.first));
                    string name2 = graph.get_path_name(graph.get_path_handle_of_step(step_pair.second));
                    REQUIRE(name1 == name2);
                    REQUIRE(correct_names.count(name1));
                }
                    
                map<string, string> true_trav_strings = {
                    {"ref", R"({"visit":[{"node_id":"11"},{"node_id":"10","backward":true},{"node_id":"9","backward":true},{"node_id":"8","backward":true}]})"},
                    {"alt1", R"({"visit":[{"node_id":"11"},{"node_id":"9","backward":true},{"node_id":"10"},{"node_id":"8","backward":true}]})"},
                    {"alt3", R"({"visit":[{"node_id":"11"},{"node_id":"9"},{"node_id":"10"},{"node_id":"8","backward":true}]})"},
                    {"alt4", R"({"visit":[{"node_id":"11"},{"node_id":"10","backward":true},{"node_id":"9","backward":true},{"node_id":"8","backward":true}]})"}
                };
                for (int i = 0; i < trav_results.first.size(); ++i) {
                    string name = graph.get_path_name(graph.get_path_handle_of_step(trav_results.second[i].first));                                
                    REQUIRE(true_trav_strings.count(name));
                    SnarlTraversal true_trav;
                    json2pb(true_trav, true_trav_strings[name]);
                    bool trav_is_correct = trav_results.first[i] == true_trav;
                    REQUIRE(trav_is_correct);
                }
            }
                
        }
        
        TEST_CASE("snarls and chains can be found in a graph with lots of root snarl connectivity", "[snarls]") {
    
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
           
            IntegratedSnarlFinder snarl_finder(graph);
            
            SECTION("Endpoints should only be seen once") {
                unordered_set<pair<id_t, bool>> seen_chain_sides;
                unordered_set<pair<id_t, bool>> seen_snarl_sides;
                snarl_finder.traverse_decomposition(
                [&](handle_t chain_start_handle){
#ifdef debug
                    cerr << "Start new chain at " << graph.get_id(chain_start_handle) << (graph.get_is_reverse(chain_start_handle) ? "rev" : "fd") << endl;
#endif
                    REQUIRE(seen_chain_sides.count(make_pair(graph.get_id(chain_start_handle), graph.get_is_reverse(chain_start_handle)))==0);
                    seen_chain_sides.emplace(graph.get_id(chain_start_handle), graph.get_is_reverse(chain_start_handle));
                },
                [&](handle_t chain_end_handle) {
#ifdef debug
                    cerr << "End new chain at " << graph.get_id(chain_end_handle) << (graph.get_is_reverse(chain_end_handle) ? "rev" : "fd") << endl;
#endif
                    REQUIRE(seen_chain_sides.count(make_pair(graph.get_id(chain_end_handle), !graph.get_is_reverse(chain_end_handle)))==0);
                    seen_chain_sides.emplace(graph.get_id(chain_end_handle), !graph.get_is_reverse(chain_end_handle));
                },
                [&](handle_t snarl_start_handle) {
#ifdef debug
                    cerr << "Start new snarl at " << graph.get_id(snarl_start_handle) << (graph.get_is_reverse(snarl_start_handle) ? "rev" : "fd") << endl;
#endif
                    REQUIRE(seen_snarl_sides.count(make_pair(graph.get_id(snarl_start_handle), graph.get_is_reverse(snarl_start_handle)))==0);
                    seen_snarl_sides.emplace(graph.get_id(snarl_start_handle), graph.get_is_reverse(snarl_start_handle));
                },
                [&](handle_t snarl_end_handle) {
#ifdef debug
                    cerr << "End new snarl at " << graph.get_id(snarl_end_handle) << (graph.get_is_reverse(snarl_end_handle) ? "rev" : "fd") << endl;
#endif
                    REQUIRE(seen_snarl_sides.count(make_pair(graph.get_id(snarl_end_handle), !graph.get_is_reverse(snarl_end_handle)))==0);
                    seen_snarl_sides.emplace(graph.get_id(snarl_end_handle), !graph.get_is_reverse(snarl_end_handle));
                });
            }
        }
    }
}
