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
            vector<vector<Snarl>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back();
            
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
                net_graph.for_each_handle([&](const handle_t& handle) {
                    nodes.insert(handle);
                });
                
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
            vector<vector<Snarl>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back();
            
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
            vector<vector<Snarl>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back();
            
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
            vector<vector<Snarl>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back();
            
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
            vector<vector<Snarl>> top_chains;
            top_chains.emplace_back();
            auto& top_chain1 = top_chains.back();
            top_chain1.emplace_back();
            auto& nested_snarl1 = top_chain1.back();
            
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
            vector<vector<Snarl>> top_chains;
            
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
            vector<vector<Snarl>> top_chains;
            
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
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
                    REQUIRE(nodes.size() == 5);
                    
                }
                
            }
            
            SECTION( "A connectivity-respecting net graph still lets you find all the nodes" ) {
            
                // Make a net graph
                NetGraph net_graph(top_snarl.start(), top_snarl.end(), top_chains, top_unary_snarls, &graph, true);
                
                SECTION( "The top-level NetGraph has 5 nodes" ) {
                    unordered_set<handle_t> nodes;
                    net_graph.for_each_handle([&](const handle_t& handle) {
                        nodes.insert(handle);
                    });
                    
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
                
                CactusSnarlFinder bubble_finder(graph);
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
                
                // Get its contents
                pair<unordered_set<Node*>, unordered_set<Edge*> > contents = snarl_manager.shallow_contents(snarl, graph, true);
                
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
    }
}
