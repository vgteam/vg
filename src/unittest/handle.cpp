/**
 * \file 
 * unittest/handle.cpp: test cases for the implementations of the HandleGraph class.
 */

#include "catch.hpp"

#include "../handle.hpp"
#include "../vg.hpp"
#include "xg.hpp"
#include "vg/io/json2pb.h"


#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"

#include <handlegraph/util.hpp>
#include <handlegraph/path_metadata.hpp>

#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <unordered_set>

namespace vg {
namespace unittest {

using namespace std;
using namespace handlegraph;

TEST_CASE( "Handle utility functions work", "[handle]" ) {

    SECTION("Handles work like ints") {
        
        SECTION("Handles are int-sized") {
            REQUIRE(sizeof(handle_t) == sizeof(int64_t));
        }
        
        SECTION("Handles can hold a range of positive integer values") {
            for (int64_t i = 0; i < 100; i++) {
                REQUIRE(as_integer(as_handle(i)) == i);
            }      
            
            REQUIRE(as_integer(as_handle(numeric_limits<int64_t>::max())) == numeric_limits<int64_t>::max());
            
        }
        
    }
    
    SECTION("Handle equality works") {
        vector<handle_t> handles;
        
        for (size_t i = 0; i < 100; i++) {
            handles.push_back(as_handle(i));
        }
        
        for (size_t i = 0; i < handles.size(); i++) {
            for (size_t j = 0; j < handles.size(); j++) {
                if (i == j) {
                    REQUIRE(handles[i] == handles[j]);
                    REQUIRE(!(handles[i] != handles[j]));
                } else {
                    REQUIRE(handles[i] != handles[j]);
                    REQUIRE(!(handles[i] == handles[j]));
                }
            }
        }
    }
    
    
    SECTION("Path handle equality works") {
        vector<path_handle_t> handles;
        
        for (size_t i = 0; i < 100; i++) {
            handles.push_back(as_path_handle(i));
        }
        
        for (size_t i = 0; i < handles.size(); i++) {
            for (size_t j = 0; j < handles.size(); j++) {
                if (i == j) {
                    REQUIRE(handles[i] == handles[j]);
                    REQUIRE(!(handles[i] != handles[j]));
                } else {
                    REQUIRE(handles[i] != handles[j]);
                    REQUIRE(!(handles[i] == handles[j]));
                }
            }
        }
    }
    
    SECTION("Step handle equality works") {
        vector<step_handle_t> handles;
        
        for (size_t i = 0; i < 10; i++) {
            for (size_t j = 0; j < 10; j++) {
                step_handle_t handle;
                as_integers(handle)[0] = i;
                as_integers(handle)[1] = j;
                handles.push_back(handle);
            }
        }
        
        for (size_t i = 0; i < handles.size(); i++) {
            for (size_t j = 0; j < handles.size(); j++) {
                if (i == j) {
                    REQUIRE(handles[i] == handles[j]);
                    REQUIRE(!(handles[i] != handles[j]));
                } else {
                    REQUIRE(handles[i] != handles[j]);
                    REQUIRE(!(handles[i] == handles[j]));
                }
            }
        }
    }
    
}

TEST_CASE("VG and XG handle implementations are correct", "[handle][vg][xg]") {
    
    // Make a vg graph
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
    Node* n9 = vg.create_node("ATATAAC");
    
    vg.create_edge(n1, n0, true, true); // a doubly reversing edge to keep it interesting
    vg.create_edge(n1, n2);
    vg.create_edge(n2, n3);
    vg.create_edge(n2, n4);
    vg.create_edge(n3, n5);
    vg.create_edge(n4, n5);
    vg.create_edge(n5, n6);
    vg.create_edge(n5, n8);
    vg.create_edge(n6, n7);
    vg.create_edge(n6, n8);
    vg.create_edge(n7, n9);
    vg.create_edge(n8, n9);
    
    // Make an xg out of it
    xg::XG xg_index;
    xg_index.from_path_handle_graph(vg);
    
    SECTION("Each graph exposes the right nodes") {
        
        for (const HandleGraph* g : {(HandleGraph*) &vg, (HandleGraph*) &xg_index}) {
            for (Node* node : {n0, n1, n2, n3, n4, n5, n6, n7, n8, n9}) {
                
                handle_t node_handle = g->get_handle(node->id(), false);
                
                SECTION("We see each node correctly forward") {
                    REQUIRE(g->get_id(node_handle) == node->id());
                    REQUIRE(g->get_is_reverse(node_handle) == false);
                    REQUIRE(g->get_sequence(node_handle) == node->sequence());
                    REQUIRE(g->get_length(node_handle) == node->sequence().size());
                }
                
                
                handle_t rev1 = g->flip(node_handle);
                handle_t rev2 = g->get_handle(node->id(), true);
                
                SECTION("We see each node correctly reverse") {
                    REQUIRE(rev1 == rev2);
                    
                    REQUIRE(g->get_id(rev1) == node->id());
                    REQUIRE(g->get_is_reverse(rev1) == true);
                    REQUIRE(g->get_sequence(rev1) == reverse_complement(node->sequence()));
                    REQUIRE(g->get_length(rev1) == node->sequence().size());
                    
                    // Check it again for good measure!
                    REQUIRE(g->get_id(rev2) == node->id());
                    REQUIRE(g->get_is_reverse(rev2) == true);
                    REQUIRE(g->get_sequence(rev2) == reverse_complement(node->sequence()));
                    REQUIRE(g->get_length(rev2) == node->sequence().size());
                    
                }
                
                
            }
        }
    
    }
    
    SECTION("Each graph exposes the right edges") {
        for (const HandleGraph* g : {(HandleGraph*) &vg, (HandleGraph*) &xg_index}) {
            // For each graph type
            for (Node* node : {n0, n1, n2, n3, n4, n5, n6, n7, n8, n9}) {
                // For each node
                for (bool orientation : {false, true}) {
                    // In each orientation
            
                    handle_t node_handle = g->get_handle(node->id(), orientation);
                    
                    vector<handle_t> next_handles;
                    vector<handle_t> prev_handles;
                    
                    // Load handles from the handle graph
                    g->follow_edges(node_handle, false, [&](const handle_t& next) {
                        next_handles.push_back(next);
                        // Exercise both returning and non-returning syntaxes
                        return true;
                    });
                    
                    g->follow_edges(node_handle, true, [&](const handle_t& next) {
                        prev_handles.push_back(next);
                        // Exercise both returning and non-returning syntaxes
                    });
                    
                    // Make sure all the entries are unique
                    REQUIRE(unordered_set<handle_t>(next_handles.begin(), next_handles.end()).size() == next_handles.size());
                    REQUIRE(unordered_set<handle_t>(prev_handles.begin(), prev_handles.end()).size() == prev_handles.size());
                    
                    // Go look up the true prev/next neighbors as NodeTraversals
                    NodeTraversal trav(node, orientation);
                    vector<NodeTraversal> true_next = vg.nodes_next(trav);
                    vector<NodeTraversal> true_prev = vg.nodes_prev(trav);
                    
                    REQUIRE(next_handles.size() == true_next.size());
                    REQUIRE(prev_handles.size() == true_prev.size());
                    
                    for (auto& handle : next_handles) {
                        // Each next handle becomes a NodeTraversal.
                        NodeTraversal handle_trav(vg.get_node(g->get_id(handle)), g->get_is_reverse(handle));
                        // And we insist on finding that traversal in the truth set.
                        REQUIRE(find(true_next.begin(), true_next.end(), handle_trav) != true_next.end());
                    }
                    
                    for (auto& handle : prev_handles) {
                        // Each next handle becomes a NodeTraversal.
                        NodeTraversal handle_trav(vg.get_node(g->get_id(handle)), g->get_is_reverse(handle));
                        // And we insist on finding that traversal in the truth set.
                        REQUIRE(find(true_prev.begin(), true_prev.end(), handle_trav) != true_prev.end());
                    }
                }
            }
        }
    }
    
    SECTION("Edge iteratees can stop early") {
        for (const HandleGraph* g : {(HandleGraph*) &vg, (HandleGraph*) &xg_index}) {
            for (Node* node : {n0, n1, n2, n3, n4, n5, n6, n7, n8, n9}) {
            
                // How many edges are we given?
                size_t loop_count = 0;
                
                handle_t node_handle = g->get_handle(node->id(), false);
                
                g->follow_edges(node_handle, false, [&](const handle_t& next) {
                    loop_count++;
                    // Never ask for more edges
                    return false;
                });
                
                // We have 1 or fewer edges on the right viewed.
                REQUIRE(loop_count <= 1);
                
                loop_count = 0;
                
                g->follow_edges(node_handle, true, [&](const handle_t& next) {
                    loop_count++;
                    // Never ask for more edges
                    return false;
                });
                
                // We have 1 or fewer edges on the left viewed.
                REQUIRE(loop_count <= 1);
            }
        }
    }
    
    SECTION("Converting handles to the forward strand works") {
        for (const HandleGraph* g : {(HandleGraph*) &vg, (HandleGraph*) &xg_index}) {
            // For each graph type
            for (Node* node : {n0, n1, n2, n3, n4, n5, n6, n7, n8, n9}) {
                // For each node
                for (bool orientation : {false, true}) {
                    // In each orientation
            
                    handle_t node_handle = g->get_handle(node->id(), orientation);
                    
                    REQUIRE(g->get_id(g->forward(node_handle)) == node->id());
                    REQUIRE(g->get_is_reverse(g->forward(node_handle)) == false);
                    
                    if (orientation) {
                        // We're reverse, so forward is our opposite
                        REQUIRE(g->forward(node_handle) == g->flip(node_handle));
                    } else {
                        // Already forward
                        REQUIRE(g->forward(node_handle) == node_handle);
                    }
                    
                }
            }
        }
    
    }
    
    SECTION("Handle pair edge cannonicalization works") {
        for (const HandleGraph* g : {(HandleGraph*) &vg, (HandleGraph*) &xg_index}) {

            SECTION("Two versions of the same edge are recognized as equal") {
                // Make the edge as it was added            
                handle_t h1 = g->get_handle(n0->id(), true);
                handle_t h2 = g->get_handle(n1->id(), true);
                pair<handle_t, handle_t> edge_as_added = g->edge_handle(h1, h2);

                // Make the edge in its simpler form
                handle_t h3 = g->get_handle(n1->id(), false);
                handle_t h4 = g->get_handle(n0->id(), false);
                pair<handle_t, handle_t> easier_edge = g->edge_handle(h3, h4);
                
                // Looking at the edge both ways must return the same result
                REQUIRE(edge_as_added == easier_edge);
                // And that result must be one of the ways of looking at the edge
                bool is_first = (edge_as_added.first == h1 && edge_as_added.second == h2);
                bool is_second = (easier_edge.first == h3 && easier_edge.second == h4);
                REQUIRE((is_first || is_second) == true);
            }
            
            SECTION("Single-sided self loops work") {
                handle_t h1 = g->get_handle(n5->id(), true);
                handle_t h2 = g->flip(h1);
                
                // Flipping this edge the other way produces the same edge.
                pair<handle_t, handle_t> only_version = make_pair(h1, h2);
                REQUIRE(g->edge_handle(only_version.first, only_version.second) == only_version);

                // We also need to handle the other end's loop                
                pair<handle_t, handle_t> other_end_loop = make_pair(h2, h1);
                REQUIRE(g->edge_handle(other_end_loop.first, other_end_loop.second) == other_end_loop);
                
                
            }
            
        }
    
    }
    
    SECTION("Node iteration works") {
        for (const HandleGraph* g : {(HandleGraph*) &vg, (HandleGraph*) &xg_index}) {
            vector<handle_t> found;
            g->for_each_handle([&](const handle_t& handle) {
                // Everything should be in its local forward orientation.
                REQUIRE(g->get_is_reverse(handle) == false);
                
                found.push_back(handle);
            });
            
            // We should have all the nodes and they should all be unique
            REQUIRE(found.size() == 10);
            REQUIRE(unordered_set<handle_t>(found.begin(), found.end()).size() == found.size());
            // They should be in the order we added them
            REQUIRE(g->get_id(found[0]) == n0->id());
            REQUIRE(g->get_id(found[1]) == n1->id());
            REQUIRE(g->get_id(found[2]) == n2->id());
            REQUIRE(g->get_id(found[3]) == n3->id());
            REQUIRE(g->get_id(found[4]) == n4->id());
            REQUIRE(g->get_id(found[5]) == n5->id());
            REQUIRE(g->get_id(found[6]) == n6->id());
            REQUIRE(g->get_id(found[7]) == n7->id());
            REQUIRE(g->get_id(found[8]) == n8->id());
            REQUIRE(g->get_id(found[9]) == n9->id());
        }
    }

    SECTION("Handle graph subsequence accessors for VG and XG work") {
        
        for (const HandleGraph* g : {(HandleGraph*) &vg, (HandleGraph*) &xg_index}) {
            REQUIRE(g->get_sequence(g->get_handle(n0->id(), false)) == "CGA");
            REQUIRE(g->get_sequence(g->get_handle(n1->id(), false)) == "TTGG");
            REQUIRE(g->get_sequence(g->get_handle(n0->id(), true)) == "TCG");
            REQUIRE(g->get_sequence(g->get_handle(n1->id(), true)) == "CCAA");
            REQUIRE(g->get_base(g->get_handle(n0->id(), false), 2) == 'A');
            REQUIRE(g->get_base(g->get_handle(n1->id(), false), 1) == 'T');
            REQUIRE(g->get_base(g->get_handle(n0->id(), true), 2) == 'G');
            REQUIRE(g->get_base(g->get_handle(n1->id(), true), 0) == 'C');
            REQUIRE(g->get_subsequence(g->get_handle(n0->id(), false), 1, 2) == "GA");
            REQUIRE(g->get_subsequence(g->get_handle(n1->id(), false), 2, 3) == "GG");
            REQUIRE(g->get_subsequence(g->get_handle(n0->id(), true), 0, 1) == "T");
            REQUIRE(g->get_subsequence(g->get_handle(n1->id(), true), 0, 3) == "CCA");
        }
    }
}

TEST_CASE("Deletable handle graphs work", "[handle][vg]") {
    
    vector<DeletableHandleGraph*> implementations;
    
    // Test the VG implementation
    VG vg;
    implementations.push_back(&vg);
    
    for(auto* g : implementations) {
    
        SECTION("No nodes exist by default") {
            size_t node_count = 0;
            g->for_each_handle([&](const handle_t& ignored) {
                node_count++;
            });
            REQUIRE(node_count == 0);
        }
    
        SECTION("A node can be added") {
            
            handle_t handle = g->create_handle("GATTACA");
            
            REQUIRE(g->get_is_reverse(handle) == false);
            REQUIRE(g->get_sequence(handle) == "GATTACA");
            REQUIRE(g->get_handle(g->get_id(handle)) == handle);
            
            SECTION("Its orientation can be changed") {
                handle_t modified = g->apply_orientation(g->flip(handle));
                
                REQUIRE(g->get_is_reverse(modified) == false);
                REQUIRE(g->get_sequence(modified) == reverse_complement("GATTACA"));
                // We don't check the ID. It's possible the ID can change.
                
                size_t node_count = 0;
                g->for_each_handle([&](const handle_t& ignored) {
                    node_count++;
                });
                REQUIRE(node_count == 1);
            }
            
            SECTION("Another node can be added") {
                handle_t handle2 = g->create_handle("CATTAG");
                
                REQUIRE(g->get_is_reverse(handle2) == false);
                REQUIRE(g->get_sequence(handle2) == "CATTAG");
                REQUIRE(g->get_handle(g->get_id(handle2)) == handle2);
                
                SECTION("The graph finds the right number of nodes") {
                    size_t node_count = 0;
                    g->for_each_handle([&](const handle_t& ignored) {
                        node_count++;
                    });
                    REQUIRE(node_count == 2);
                }
                
                SECTION("No edges exist by default") {
                    
                    // Grab all the edges            
                    vector<pair<handle_t, handle_t>> edges;
                    g->follow_edges(handle, false, [&](const handle_t& other) {
                        edges.push_back(g->edge_handle(handle, other));
                    });
                    g->follow_edges(handle, true, [&](const handle_t& other) {
                        edges.push_back(g->edge_handle(other, handle));
                    });
                    
                    REQUIRE(edges.size() == 0);    
                    
                }
                
                SECTION("Edges can be added") {

                    // Test deduplication
                    g->create_edge(handle, handle2);
                    g->create_edge(g->flip(handle2), g->flip(handle));
                    g->create_edge(handle, handle2);

                    // Grab all the edges            
                    vector<pair<handle_t, handle_t>> edges;
                    g->follow_edges(handle, false, [&](const handle_t& other) {
                        edges.push_back(g->edge_handle(handle, other));
                    });
                    g->follow_edges(handle, true, [&](const handle_t& other) {
                        edges.push_back(g->edge_handle(other, handle));
                    });
                    
                    REQUIRE(edges.size() == 1);
                    REQUIRE(edges.front() == g->edge_handle(handle, handle2));
                    
                    SECTION("Reorienting nodes modifies edges") {
                        handle_t modified = g->apply_orientation(g->flip(handle));
                        
                        // Grab all the edges            
                        vector<pair<handle_t, handle_t>> edges;
                        g->follow_edges(modified, false, [&](const handle_t& other) {
                            edges.push_back(g->edge_handle(modified, other));
                        });
                        g->follow_edges(modified, true, [&](const handle_t& other) {
                            edges.push_back(g->edge_handle(other, modified));
                        });
                        
                        REQUIRE(edges.size() == 1);
                        REQUIRE(edges.front() == g->edge_handle(g->flip(handle2), modified));
                        
                        
                    }
                
                }
                
            }
            
            SECTION("A node can be split") {
                // Should get GATT and ACA, but in reverse (TGT, AATC)
                auto parts = g->divide_handle(g->flip(handle), 3);
                
                REQUIRE(g->get_sequence(parts.first) == "TGT");
                REQUIRE(g->get_is_reverse(parts.first) == true);
                REQUIRE(g->get_sequence(parts.second) == "AATC");
                REQUIRE(g->get_is_reverse(parts.second) == true);
                
                SECTION("The original node is gone") {
                    size_t node_count = 0;
                    g->for_each_handle([&](const handle_t& ignored) {
                        node_count++;
                    });
                    REQUIRE(node_count == 2);
                }
                
                SECTION("Splitting creates the appropriate edge") {
                    vector<handle_t> found;
                    g->follow_edges(parts.first, false, [&](const handle_t& other) {
                        found.push_back(other);
                    });
                    
                    REQUIRE(found.size() == 1);
                    REQUIRE(found.front() == parts.second);
                }
                
                SECTION("An edge can be removed") {
                    g->destroy_edge(parts.first, parts.second);
                    
                    vector<handle_t> found;
                    g->follow_edges(parts.first, false, [&](const handle_t& other) {
                        found.push_back(other);
                    });
                    
                    REQUIRE(found.size() == 0);
                }
                
                SECTION("A node can be removed") {
                    g->destroy_handle(parts.second);
                    
                    vector<handle_t> found;
                    g->follow_edges(parts.first, false, [&](const handle_t& other) {
                        found.push_back(other);
                    });
                    
                    REQUIRE(found.size() == 0);
                    
                    size_t node_count = 0;
                    g->for_each_handle([&](const handle_t& ignored) {
                        node_count++;
                    });
                    REQUIRE(node_count == 1);
                }
            }
        }
    }
}

TEST_CASE("divide_handle works on both strands", "[handle][vg][packed][hashgraph]") {
    vector<pair<MutablePathMutableHandleGraph*, MutablePathMutableHandleGraph*>> implementations;
    
    // Add implementations
    
    VG vg, vg2;
    implementations.push_back(make_pair(&vg, &vg2));

    bdsg::HashGraph hg, hg2;
    implementations.push_back(make_pair(&hg, &hg2));

    bdsg::PackedGraph pg, pg2;
    implementations.push_back(make_pair(&pg, &pg2));
        
    // And test them
    for (int imp = 0; imp < implementations.size(); ++imp) {

        for (bool backwards : {false, true}) {

            MutablePathMutableHandleGraph* g = backwards ? implementations[imp].first : implementations[imp].second;
        
            REQUIRE(g->get_node_count() == 0);

            handle_t handle = g->create_handle("TTATATTCCAACTCTCTG");
            if (backwards) {
                handle = g->flip(handle);
            }
            path_handle_t path_handle = g->create_path_handle("Path");
            g->append_step(path_handle, handle);
            string seq = g->get_sequence(handle);
            vector<string> true_parts = { seq.substr(0, 1), seq.substr(1, 4), seq.substr(5, 5), seq.substr(10) };

            // Should get (C,AGAG,AGTTG,GAATATAA)  (forward)
            // Should get (T,TATA,TTCCA,ACTCTCTG)  (reverse)
            auto parts = g->divide_handle(handle, {1, 5, 10});
            
            REQUIRE(parts.size() == true_parts.size());
            for (int i = 0; i < parts.size(); ++i) {
                REQUIRE(g->get_sequence(parts[i]) == true_parts[i]);
                REQUIRE(g->get_is_reverse(parts[i]) == backwards);
            }
        
            vector<handle_t> steps;
            g->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                    steps.push_back(g->get_handle_of_step(step_handle));
                });
            REQUIRE(steps.size() == true_parts.size());
            for (int i = 0; i < parts.size(); ++i) {
                REQUIRE(g->get_sequence(steps[i]) == true_parts[i]);
                REQUIRE(g->get_is_reverse(steps[i]) == backwards);
            }
        }
    }
}

TEST_CASE("divide_handle and destroy_handle work together on both strands", "[handle][vg][packed][hashgraph]") {
    vector<pair<MutablePathMutableHandleGraph*, MutablePathMutableHandleGraph*>> implementations;
    
    // Add implementations
    
    VG vg, vg2;
    implementations.push_back(make_pair(&vg, &vg2));

    bdsg::HashGraph hg, hg2;
    implementations.push_back(make_pair(&hg, &hg2));

    bdsg::PackedGraph pg, pg2;
    implementations.push_back(make_pair(&pg, &pg2));
        
    // And test them
    for (int imp = 0; imp < implementations.size(); ++imp) {

        for (bool backwards : {false, true}) {

            MutablePathMutableHandleGraph* g = backwards ? implementations[imp].first : implementations[imp].second;
        
            REQUIRE(g->get_node_count() == 0);

            handle_t handle1 = g->create_handle("CAAATAAGGCTTGGAAATTTTCTGGAGTTCTA");
            handle_t handle2 = g->create_handle("TTATATTCCAACTCTCTG");
            path_handle_t path_handle = g->create_path_handle("x");
            g->create_edge(handle1, handle2);
            
            if (backwards) {
                handle1 = g->flip(handle1);
                handle2 = g->flip(handle2);
                g->append_step(path_handle, handle2);
                g->append_step(path_handle, handle1);                
            } else {
                g->append_step(path_handle, handle1);
                g->append_step(path_handle, handle2);                
            }

            g->divide_handle(handle2, vector<size_t>({1, 5, 10}));
            g->divide_handle(handle1, vector<size_t>({2, 7, 22, 31}));

            vector<handle_t> steps;
            g->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                    steps.push_back(g->get_handle_of_step(step_handle));
                });
            REQUIRE(steps.size() == 9);
            int i = 0;
            vector<handle_t> to_delete;
            g->append_step(g->create_path_handle(to_string(i)), steps[i]);
            i++;
            to_delete.push_back(steps[i++]);
            g->append_step(g->create_path_handle(to_string(i)), steps[i]);
            i++;
            to_delete.push_back(steps[i++]);
            to_delete.push_back(steps[i++]);
            to_delete.push_back(steps[i++]);
            g->append_step(g->create_path_handle(to_string(i)), steps[i]);
            i++;
            to_delete.push_back(steps[i++]);
            g->append_step(g->create_path_handle(to_string(i)), steps[i]);
            i++;
                                
            g->destroy_path(path_handle);

            for (auto handle : to_delete) {
                dynamic_cast<DeletableHandleGraph*>(g)->destroy_handle(handle);
            }

            REQUIRE(g->get_node_count() == 4);
            REQUIRE(g->get_path_count() == 4);
        }
    }
}

TEST_CASE("DeletableHandleGraphs are correct", "[handle][vg][packed][hashgraph]") {
    
    vector<DeletableHandleGraph*> implementations;
    
    // Add implementations
    
    VG vg;
    implementations.push_back(&vg);
    
    bdsg::PackedGraph pg;
    implementations.push_back(&pg);
    
    bdsg::HashGraph hg;
    implementations.push_back(&hg);
    
    // And test them
    
    for (DeletableHandleGraph* implementation : implementations) {
        
        DeletableHandleGraph& graph = *implementation;
        
        REQUIRE(graph.get_node_count() == 0);
        
        handle_t h = graph.create_handle("ATG", 2);
        
        // DeletableHandleGraph has correct structure after creating a node
        {
            REQUIRE(graph.get_sequence(h) == "ATG");
            REQUIRE(graph.get_sequence(graph.flip(h)) == "CAT");
            REQUIRE(graph.get_base(h, 1) == 'T');
            REQUIRE(graph.get_base(graph.flip(h), 2) == 'T');
            REQUIRE(graph.get_subsequence(h, 1, 3) == "TG");
            REQUIRE(graph.get_subsequence(graph.flip(h), 0, 2) == "CA");
            REQUIRE(graph.get_length(h) == 3);
            REQUIRE(graph.has_node(graph.get_id(h)));
            REQUIRE(!graph.has_node(graph.get_id(h) + 1));
            
            REQUIRE(graph.get_handle(graph.get_id(h)) == h);
            REQUIRE(!graph.get_is_reverse(h));
            REQUIRE(graph.get_is_reverse(graph.flip(h)));
            
            REQUIRE(graph.get_node_count() == 1);
            REQUIRE(graph.min_node_id() == graph.get_id(h));
            REQUIRE(graph.max_node_id() == graph.get_id(h));
            
            graph.follow_edges(h, true, [](const handle_t& prev) {
                REQUIRE(false);
                return true;
            });
            graph.follow_edges(h, false, [](const handle_t& next) {
                REQUIRE(false);
                return true;
            });
        }
        
        handle_t h2 = graph.create_handle("CT", 1);
        
        // DeletableHandleGraph has correct structure after creating a node at the beginning of ID space
        {
            
            REQUIRE(graph.get_sequence(h2) == "CT");
            REQUIRE(graph.get_sequence(graph.flip(h2)) == "AG");
            REQUIRE(graph.get_base(h2, 1) == 'T');
            REQUIRE(graph.get_base(graph.flip(h2), 0) == 'A');
            REQUIRE(graph.get_subsequence(h2, 1, 10) == "T");
            REQUIRE(graph.get_subsequence(graph.flip(h2), 0, 2) == "AG");
            REQUIRE(graph.get_length(h2) == 2);
            REQUIRE(graph.has_node(graph.get_id(h2)));
            REQUIRE(!graph.has_node(max(graph.get_id(h), graph.get_id(h2)) + 1));
            
            REQUIRE(graph.get_handle(graph.get_id(h2)) == h2);
            
            REQUIRE(graph.get_node_count() == 2);
            REQUIRE(graph.min_node_id() == graph.get_id(h2));
            REQUIRE(graph.max_node_id() == graph.get_id(h));
            
            graph.follow_edges(h2, true, [](const handle_t& prev) {
                REQUIRE(false);
                return true;
            });
            graph.follow_edges(h2, false, [](const handle_t& next) {
                REQUIRE(false);
                return true;
            });
        }
        
        // creating and accessing a node at the end of ID space
        
        handle_t h3 = graph.create_handle("GAC", 4);
        
        // DeletableHandleGraph has correct structure after creating a node at the end of ID space
        {
            REQUIRE(graph.get_sequence(h3) == "GAC");
            REQUIRE(graph.get_sequence(graph.flip(h3)) == "GTC");
            REQUIRE(graph.get_base(h3, 1) == 'A');
            REQUIRE(graph.get_base(graph.flip(h3), 0) == 'G');
            REQUIRE(graph.get_subsequence(h3, 1, 1) == "A");
            REQUIRE(graph.get_subsequence(graph.flip(h3), 0, 5) == "GTC");
            REQUIRE(graph.get_length(h3) == 3);
            
            REQUIRE(graph.get_handle(graph.get_id(h3)) == h3);
            
            REQUIRE(graph.get_node_count() == 3);
            REQUIRE(graph.min_node_id() == graph.get_id(h2));
            REQUIRE(graph.max_node_id() == graph.get_id(h3));
            
            graph.follow_edges(h3, true, [](const handle_t& prev) {
                REQUIRE(false);
                return true;
            });
            graph.follow_edges(h3, false, [](const handle_t& next) {
                REQUIRE(false);
                return true;
            });
        }
        
        
        // creating and accessing in the middle of ID space
        
        handle_t h4 = graph.create_handle("T", 3);
        
        // DeletableHandleGraph has correct structure after creating a node in the middle of ID space
        {
            REQUIRE(graph.get_sequence(h4) == "T");
            REQUIRE(graph.get_sequence(graph.flip(h4)) == "A");
            REQUIRE(graph.get_length(h4) == 1);
            
            REQUIRE(graph.get_handle(graph.get_id(h4)) == h4);
            
            REQUIRE(graph.get_node_count() == 4);
            REQUIRE(graph.min_node_id() == graph.get_id(h2));
            REQUIRE(graph.max_node_id() == graph.get_id(h3));
            
            graph.follow_edges(h4, true, [](const handle_t& prev) {
                REQUIRE(false);
                return true;
            });
            graph.follow_edges(h4, false, [](const handle_t& next) {
                REQUIRE(false);
                return true;
            });
        }
        
        graph.create_edge(h, h2);
        
        bool found1 = false, found2 = false, found3 = false, found4 = false;
        int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
        
        // DeletableHandleGraph has correct structure after creating an edge
        {
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found4 = true;
                }
                count4++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            
            count1 = count2 = count3 = count4 = 0;
            found1 = found2 = found3 = found4 = false;
        }
        
        graph.create_edge(h, graph.flip(h3));
        
        bool found5 = false, found6 = false, found7 = false, found8 = false;
        int count5 = 0, count6 = 0;
        
        // DeletableHandleGraph has correct structure after creating an edge with a traversal
        {
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
        }
        
        graph.create_edge(h4, graph.flip(h4));
        
        // DeletableHandleGraph has correct structure after creating a reversing self-loop
        {
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        graph.create_edge(h, graph.flip(h4));
        graph.create_edge(graph.flip(h3), h4);
        
        graph.destroy_edge(h, graph.flip(h4));
        graph.destroy_edge(graph.flip(h3), h4);
        
        // DeletableHandleGraph has correct structure after creating and deleting edges
        {
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        handle_t h5 = graph.create_handle("GGACC");
        
        // make some edges to ensure that deleting is difficult
        graph.create_edge(h, h5);
        graph.create_edge(h5, h);
        graph.create_edge(graph.flip(h5), h2);
        graph.create_edge(h3, graph.flip(h5));
        graph.create_edge(h3, h5);
        graph.create_edge(h5, h4);
        
        graph.destroy_handle(h5);
        
        // DeletableHandleGraph has correct structure after creating and deleting a node
        {
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        // DeletableHandleGraph has correct structure after swapping nodes
        {
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        // DeletableHandleGraph visits all nodes with for_each_handle
        {
            graph.for_each_handle([&](const handle_t& handle) {
                if (handle == h) {
                    found1 = true;
                }
                else if (handle == h2) {
                    found2 = true;
                }
                else if (handle == h3) {
                    found3 = true;
                }
                else if (handle == h4) {
                    found4 = true;
                }
                else {
                    REQUIRE(false);
                }
                return true;
            });
            
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            
            found1 = found2 = found3 = found4 = false;
        }
        
        // to make sure the sequence reverse complemented correctly
        int i = 0;
        auto check_rev_comp = [&](const std::string& seq1, const std::string& seq2) {
            i++;
            REQUIRE(seq1.size() == seq2.size());
            auto it = seq1.begin();
            auto rit = seq2.rbegin();
            for (; it != seq1.end(); it++) {
                if (*it == 'A') {
                    REQUIRE(*rit == 'T');
                }
                else if (*it == 'C') {
                    REQUIRE(*rit == 'G');
                }
                else if (*it == 'G') {
                    REQUIRE(*rit == 'C');
                }
                else if (*it == 'T') {
                    REQUIRE(*rit == 'A');
                }
                else if (*it == 'N') {
                    REQUIRE(*rit == 'N');
                }
                else {
                    REQUIRE(false);
                }
                
                rit++;
            }
        };
        
        
        int count7 = 0, count8 = 0;
        
        // DeletableHandleGraph correctly reverses a node
        {
            
            string seq1 = graph.get_sequence(h);
            h = graph.apply_orientation(graph.flip(h));
            
            // check the sequence
            string rev_seq1 = graph.get_sequence(h);
            check_rev_comp(seq1, rev_seq1);
            
            // check that the edges are what we expect
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                count1++;
                return true;
            });
            graph.follow_edges(h, true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found1 = true;
                }
                else if (prev == h3) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& next) {
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h), false, [&](const handle_t& prev) {
                if (prev == h2) {
                    found3 = true;
                }
                else if (prev == graph.flip(h3)) {
                    found4 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == graph.flip(h)) {
                    found5 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == h) {
                    found6 = true;
                }
                count6++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h)) {
                    found7 = true;
                }
                count7++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == h) {
                    found8 = true;
                }
                count8++;
                return true;
            });
            REQUIRE(count1 == 0);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 0);
            REQUIRE(count4 == 2);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(count7 == 1);
            REQUIRE(count8 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = count7 = count8 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            
            // and now switch it back to the same orientation and repeat the topology checks
            
            h = graph.apply_orientation(graph.flip(h));
            
            graph.follow_edges(h, false, [&](const handle_t& next) {
                if (next == h2) {
                    found1 = true;
                }
                else if (next == graph.flip(h3)) {
                    found2 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found3 = true;
                }
                else if (prev == h3) {
                    found4 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == h) {
                    found5 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found6 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == h) {
                    found7 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(h3, false, [&](const handle_t& next) {
                if (next == graph.flip(h)) {
                    found8 = true;
                }
                count6++;
                return true;
            });
            REQUIRE(count1 == 2);
            REQUIRE(count2 == 2);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            
            count1 = count2 = count3 = count4 = count5 = count6 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
            
            graph.follow_edges(h4, false, [&](const handle_t& next) {
                if (next == graph.flip(h4)) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
                if (prev == h4) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            
            count1 = count2 = 0;
            found1 = found2 = false;
        }
        
        vector<handle_t> parts = graph.divide_handle(h, vector<size_t>{1, 2});
        
        int count9 = 0, count10 = 0, count11 = 0, count12 = 0;
        bool found9 = false, found10 = false, found11 = false, found12 = false, found13 = false, found14 = false;
        
        // DeletableHandleGraph can correctly divide a node
        {
            
            REQUIRE(parts.size() == 3);
            
            REQUIRE(graph.get_sequence(parts[0]) == "A");
            REQUIRE(graph.get_length(parts[0]) == 1);
            REQUIRE(graph.get_sequence(parts[1]) == "T");
            REQUIRE(graph.get_length(parts[1]) == 1);
            REQUIRE(graph.get_sequence(parts[2]) == "G");
            REQUIRE(graph.get_length(parts[2]) == 1);
            
            
            graph.follow_edges(parts[0], false, [&](const handle_t& next) {
                if (next == parts[1]) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(parts[0], true, [&](const handle_t& prev) {
                count2++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[0]), true, [&](const handle_t& prev) {
                if (prev == graph.flip(parts[1])) {
                    found2 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[0]), false, [&](const handle_t& next) {
                count4++;
                return true;
            });
            
            graph.follow_edges(parts[1], false, [&](const handle_t& next) {
                if (next == parts[2]) {
                    found3 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(parts[1], true, [&](const handle_t& prev) {
                if (prev == parts[0]) {
                    found4 = true;
                }
                count6++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[1]), true, [&](const handle_t& prev) {
                if (prev == graph.flip(parts[2])) {
                    found5 = true;
                }
                count7++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[1]), false, [&](const handle_t& next) {
                if (next == graph.flip(parts[0])) {
                    found6 = true;
                }
                count8++;
                return true;
            });
            
            graph.follow_edges(parts[2], false, [&](const handle_t& next) {
                if (next == h2) {
                    found7 = true;
                }
                else if (next == graph.flip(h3)) {
                    found8 = true;
                }
                count9++;
                return true;
            });
            graph.follow_edges(parts[2], true, [&](const handle_t& prev) {
                if (prev == parts[1]) {
                    found9 = true;
                }
                count10++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[2]), true, [&](const handle_t& prev) {
                if (prev == graph.flip(h2)) {
                    found10 = true;
                }
                else if (prev == h3) {
                    found11 = true;
                }
                count11++;
                return true;
            });
            graph.follow_edges(graph.flip(parts[2]), false, [&](const handle_t& next) {
                if (next == graph.flip(parts[1])) {
                    found12 = true;
                }
                count12++;
                return true;
            });
            graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
                if (prev == parts[2]) {
                    found13 = true;
                }
                return true;
            });
            graph.follow_edges(h2, true, [&](const handle_t& prev) {
                if (prev == parts[2]) {
                    found14 = true;
                }
                return true;
            });
            
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 0);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 0);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 1);
            REQUIRE(count7 == 1);
            REQUIRE(count8 == 1);
            REQUIRE(count9 == 2);
            REQUIRE(count10 == 1);
            REQUIRE(count11 == 2);
            REQUIRE(count12 == 1);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
            REQUIRE(found6);
            REQUIRE(found7);
            REQUIRE(found8);
            REQUIRE(found9);
            REQUIRE(found10);
            REQUIRE(found11);
            REQUIRE(found12);
            REQUIRE(found13);
            REQUIRE(found14);
            
            count1 = count2 = count3 = count4 = count5 = count6 = count7 = count8 = count9 = count10 = count11 = count12 = 0;
            found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = found9 = found10 = found11 = found12 = false;
        }
        
        vector<handle_t> rev_parts = graph.divide_handle(graph.flip(h3), vector<size_t>{1});
        
        // DeletableHandleGraph can correctly divide a node on the reverse strand
        {
            
            REQUIRE(graph.get_sequence(rev_parts[0]) == "G");
            REQUIRE(graph.get_length(rev_parts[0]) == 1);
            REQUIRE(graph.get_is_reverse(rev_parts[0]));
            REQUIRE(graph.get_sequence(rev_parts[1]) == "TC");
            REQUIRE(graph.get_length(rev_parts[1]) == 2);
            REQUIRE(graph.get_is_reverse(rev_parts[1]));
            
            graph.follow_edges(rev_parts[0], false, [&](const handle_t& next) {
                if (next == rev_parts[1]) {
                    found1 = true;
                }
                count1++;
                return true;
            });
            graph.follow_edges(rev_parts[1], true, [&](const handle_t& prev) {
                if (prev == rev_parts[0]) {
                    found2 = true;
                }
                count2++;
                return true;
            });
            graph.follow_edges(graph.flip(rev_parts[1]), false, [&](const handle_t& next) {
                if (next == graph.flip(rev_parts[0])) {
                    found3 = true;
                }
                count3++;
                return true;
            });
            graph.follow_edges(graph.flip(rev_parts[0]), true, [&](const handle_t& prev) {
                if (prev == graph.flip(rev_parts[1])) {
                    found4 = true;
                }
                count4++;
                return true;
            });
            graph.follow_edges(rev_parts[0], true, [&](const handle_t& prev) {
                if (prev == parts[2]) {
                    found5 = true;
                }
                count5++;
                return true;
            });
            graph.follow_edges(rev_parts[1], false, [&](const handle_t& next) {
                count6++;
                return true;
            });
            
            REQUIRE(count1 == 1);
            REQUIRE(count2 == 1);
            REQUIRE(count3 == 1);
            REQUIRE(count4 == 1);
            REQUIRE(count5 == 1);
            REQUIRE(count6 == 0);
            REQUIRE(found1);
            REQUIRE(found2);
            REQUIRE(found3);
            REQUIRE(found4);
            REQUIRE(found5);
        }
    }
}
    
TEST_CASE("VG and XG path handle implementations are correct", "[handle][vg][xg]") {
    
    // Make a vg graph
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
    Node* n9 = vg.create_node("ATATAAC");
    
    vg.create_edge(n1, n0, true, true); // a doubly reversing edge to keep it interesting
    vg.create_edge(n1, n2);
    vg.create_edge(n2, n3);
    vg.create_edge(n2, n4);
    vg.create_edge(n3, n5);
    vg.create_edge(n4, n5);
    vg.create_edge(n5, n6);
    vg.create_edge(n5, n8);
    vg.create_edge(n6, n7);
    vg.create_edge(n6, n8);
    vg.create_edge(n7, n9);
    vg.create_edge(n8, n9);
    
    Path path1 ,path2, path3;
    path1.set_name("1");
    path2.set_name("2");
    path3.set_name("3");
    
    Mapping* m10 = path1.add_mapping();
    m10->mutable_position()->set_node_id(n0->id());
    m10->set_rank(1);
    Mapping* m11 = path1.add_mapping();
    m11->mutable_position()->set_node_id(n1->id());
    m11->set_rank(2);
    Mapping* m12 = path1.add_mapping();
    m12->mutable_position()->set_node_id(n2->id());
    m12->set_rank(3);
    Mapping* m13 = path1.add_mapping();
    m13->mutable_position()->set_node_id(n4->id());
    m13->set_rank(4);
    Mapping* m14 = path1.add_mapping();
    m14->mutable_position()->set_node_id(n5->id());
    m14->set_rank(5);
    
    Mapping* m20 = path2.add_mapping();
    m20->mutable_position()->set_node_id(n3->id());
    m20->set_rank(1);
    Mapping* m21 = path2.add_mapping();
    m21->mutable_position()->set_node_id(n5->id());
    m21->set_rank(2);
    Mapping* m22 = path2.add_mapping();
    m22->mutable_position()->set_node_id(n6->id());
    m22->set_rank(3);
    Mapping* m23 = path2.add_mapping();
    m23->mutable_position()->set_node_id(n7->id());
    m23->set_rank(4);
    Mapping* m24 = path2.add_mapping();
    m24->mutable_position()->set_node_id(n9->id());
    m24->set_rank(5);
    path2.set_is_circular(true);
    
    Mapping* m30 = path3.add_mapping();
    m30->mutable_position()->set_node_id(n8->id());
    m30->mutable_position()->set_is_reverse(true);
    m30->set_rank(1);
    Mapping* m31 = path3.add_mapping();
    m31->mutable_position()->set_node_id(n5->id());
    m31->mutable_position()->set_is_reverse(true);
    m31->set_rank(2);
    Mapping* m32 = path3.add_mapping();
    m32->mutable_position()->set_node_id(n3->id());
    m32->mutable_position()->set_is_reverse(true);
    m32->set_rank(3);
    
    vg.paths.extend(path1);
    vg.paths.extend(path2);
    vg.paths.extend(path3);
    
    // also add the paths to the Protobuf graph so that they're XG'able
    vg.paths.to_graph(vg.graph);
    
    xg::XG xg_index;
    xg_index.from_path_handle_graph(vg);
    
    SECTION("Handles can find all paths") {
        
        int vg_path_count = 0, xg_path_count = 0;
        
        function<void(const path_handle_t&)> count_vg_paths = [&](const path_handle_t& ph) {
            vg_path_count++;
        };
        
        function<void(const path_handle_t&)> count_xg_paths = [&](const path_handle_t& ph) {
            xg_path_count++;
        };
        
        vg.for_each_path_handle(count_vg_paths);
        xg_index.for_each_path_handle(count_xg_paths);
        
        REQUIRE(vg_path_count == 3);
        REQUIRE(xg_path_count == 3);
        REQUIRE(vg.get_path_count() == 3);
        REQUIRE(xg_index.get_path_count() == 3);
    }
    
    // check that a path is correctly accessible and traversible
    auto check_path_traversal = [&](const PathHandleGraph& graph, const Path& path) {
        
        path_handle_t path_handle = graph.get_path_handle(path.name());
        
        REQUIRE(graph.get_step_count(path_handle) == path.mapping_size());
        
        // check that steps is pointing to the same index along the path
        auto check_step = [&](const step_handle_t& step_handle,
                              int mapping_idx) {
            REQUIRE(graph.get_path_handle_of_step(step_handle) == path_handle);
            
            const Mapping& mapping = path.mapping(mapping_idx);
            
            handle_t handle = graph.get_handle_of_step(step_handle);
            
            REQUIRE(graph.get_id(handle) == mapping.position().node_id());
            REQUIRE(graph.get_is_reverse(handle) == mapping.position().is_reverse());
            
            bool should_have_next = path.is_circular() || mapping_idx + 1 < path.mapping_size();
            bool should_have_prev = path.is_circular() || mapping_idx > 0;
            
            REQUIRE(graph.has_next_step(step_handle) == should_have_next);
            REQUIRE(graph.has_previous_step(step_handle) == should_have_prev);
        };
        
        step_handle_t step_handle = graph.path_begin(path_handle);
        
        // iterate front to back
        for (int i = 0; i < path.mapping_size(); i++) {
            if (i > 0) {
                step_handle_t previous = graph.get_previous_step(step_handle);
                check_step(previous, i - 1);
            }
            else if (path.is_circular()) {
                step_handle_t previous = graph.get_previous_step(step_handle);
                check_step(previous, path.mapping_size() - 1);
            }
            
            check_step(step_handle, i);
            step_handle = graph.get_next_step(step_handle);
        }
        
        if (path.is_circular() && path.mapping_size() > 0) {
            REQUIRE(step_handle == graph.path_begin(path_handle));
        }
        else {
            REQUIRE(step_handle == graph.path_end(path_handle));
        }
        
        // iterate front to back with the iteration function
        {
            int i = 0;
            graph.for_each_step_in_path(path_handle, [&i, &check_step](const step_handle_t& step_handle) {
                check_step(step_handle, i);
                i++;
            });
            REQUIRE(i == path.mapping_size());
        }
    };
    
    SECTION("Handles can traverse paths") {
        
        check_path_traversal(vg, path1);
        check_path_traversal(vg, path2);
        check_path_traversal(vg, path3);
        check_path_traversal(xg_index, path1);
        check_path_traversal(xg_index, path2);
        check_path_traversal(xg_index, path3);
    }
    
    SECTION("VG can rewrite a segment of a path") {
        
        // get the segments we want to replace
        
        path_handle_t path_handle_1 = vg.get_path_handle("1");
        step_handle_t segment_begin_1 = vg.get_next_step(vg.get_next_step(vg.path_begin(path_handle_1)));
        step_handle_t segment_end_1 = vg.get_next_step(vg.get_next_step(segment_begin_1));
        
        path_handle_t path_handle_2 = vg.get_path_handle("2");
        step_handle_t segment_begin_2 = vg.path_begin(path_handle_2);
        step_handle_t segment_end_2 = vg.get_next_step(segment_begin_2);
        
        path_handle_t path_handle_3 = vg.get_path_handle("3");
        step_handle_t segment_end_3 = vg.path_end(path_handle_3);
        step_handle_t segment_begin_3 = vg.get_previous_step(vg.get_previous_step(segment_end_3));
        
        // replace them
        auto new_segment_1 = vg.rewrite_segment(segment_begin_1, segment_end_1, {vg.get_handle(n2->id()), vg.get_handle(n3->id())});
        auto new_segment_2 = vg.rewrite_segment(segment_begin_2, segment_end_2, {vg.get_handle(n4->id())});
        auto new_segment_3 = vg.rewrite_segment(segment_begin_3, segment_end_3, {vg.get_handle(n6->id(), true), vg.get_handle(n5->id(), true), vg.get_handle(n4->id(), true)});
        
        // update the Path objects to match our edits
        path1.mutable_mapping(3)->mutable_position()->set_node_id(n3->id());
        
        path2.mutable_mapping(0)->mutable_position()->set_node_id(n4->id());
        
        path3.mutable_mapping(1)->mutable_position()->set_node_id(n6->id());
        path3.mutable_mapping(2)->mutable_position()->set_node_id(n5->id());
        Mapping* mapping = path3.add_mapping();
        Position* position = mapping->mutable_position();
        position->set_node_id(n4->id());
        position->set_is_reverse(true);
        
        // check to make sure the handle backings are correct
        check_path_traversal(vg, path1);
        check_path_traversal(vg, path2);
        check_path_traversal(vg, path3);
        
        
        // maintains validity after serialization (had a bug with this after rewrites at
        // one point)
        stringstream strm;
        vg.serialize_to_ostream(strm);
        strm.seekg(0);
        VG copy(strm);
        
        REQUIRE(handlealgs::are_equivalent_with_paths(&vg, &copy));
        
        // TODO: check the replaced segments' handles
    }
}
    
TEST_CASE("Deletable handle graphs behave correctly when a graph has multiple edges between the same pair of nodes", "[handle][vg][packed][hashgraph]") {
    
    vector<DeletableHandleGraph*> implementations;
    
    bdsg::PackedGraph pg;
    implementations.push_back(&pg);
    
    bdsg::HashGraph hg;
    implementations.push_back(&hg);
    
    VG vg;
    implementations.push_back(&vg);
    
    for(DeletableHandleGraph* implementation : implementations) {
        
        DeletableHandleGraph& graph = *implementation;
        
        // initialize the graph
        
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("C");
        
        graph.create_edge(h1, h2);
        graph.create_edge(graph.flip(h1), h2);
        
        // test for the right initial topology
        bool found1 = false, found2 = false, found3 = false, found4 = false, found5 = false, found6 = false;
        int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
        
        graph.follow_edges(h1, false, [&](const handle_t& other) {
            if (other == h2) {
                found1 = true;
            }
            count1++;
        });
        graph.follow_edges(h1, true, [&](const handle_t& other) {
            if (other == graph.flip(h2)) {
                found2 = true;
            }
            count2++;
        });
        graph.follow_edges(h2, false, [&](const handle_t& other) {
            count3++;
        });
        graph.follow_edges(h2, true, [&](const handle_t& other) {
            if (other == h1) {
                found3 = true;
            }
            else if (other == graph.flip(h1)) {
                found4 = true;
            }
            count4++;
        });
        REQUIRE(found1);
        REQUIRE(found2);
        REQUIRE(found3);
        REQUIRE(found4);
        REQUIRE(count1 == 1);
        REQUIRE(count2 == 1);
        REQUIRE(count3 == 0);
        REQUIRE(count4 == 2);
        found1 = found2 = found3 = found4 = found5 = found6 = false;
        count1 = count2 = count3 = count4 = 0;
        
        // flip a node and check if the orientation is correct
        h1 = graph.apply_orientation(graph.flip(h1));
        
        graph.follow_edges(h1, false, [&](const handle_t& other) {
            if (other == h2) {
                found1 = true;
            }
            count1++;
        });
        graph.follow_edges(h1, true, [&](const handle_t& other) {
            if (other == graph.flip(h2)) {
                found2 = true;
            }
            count2++;
        });
        graph.follow_edges(h2, false, [&](const handle_t& other) {
            count3++;
        });
        graph.follow_edges(h2, true, [&](const handle_t& other) {
            if (other == h1) {
                found3 = true;
            }
            else if (other == graph.flip(h1)) {
                found4 = true;
            }
            count4++;
        });
        REQUIRE(found1);
        REQUIRE(found2);
        REQUIRE(found3);
        REQUIRE(found4);
        REQUIRE(count1 == 1);
        REQUIRE(count2 == 1);
        REQUIRE(count3 == 0);
        REQUIRE(count4 == 2);
        found1 = found2 = found3 = found4 = found5 = found6 = false;
        count1 = count2 = count3 = count4 = 0;
        
        // create a new edge
        
        graph.create_edge(h1, graph.flip(h2));
        
        // check the topology
        
        graph.follow_edges(h1, false, [&](const handle_t& other) {
            if (other == h2) {
                found1 = true;
            }
            else if (other == graph.flip(h2)) {
                found2 = true;
            }
            count1++;
        });
        graph.follow_edges(h1, true, [&](const handle_t& other) {
            if (other == graph.flip(h2)) {
                found3 = true;
            }
            count2++;
        });
        graph.follow_edges(h2, false, [&](const handle_t& other) {
             if (other == graph.flip(h1)) {
                found4 = true;
            }
            count3++;
        });
        graph.follow_edges(h2, true, [&](const handle_t& other) {
            if (other == h1) {
                found5 = true;
            }
            else if (other == graph.flip(h1)) {
                found6 = true;
            }
            count4++;
        });
        REQUIRE(found1);
        REQUIRE(found2);
        REQUIRE(found3);
        REQUIRE(found4);
        REQUIRE(found5);
        REQUIRE(found6);
        REQUIRE(count1 == 2);
        REQUIRE(count2 == 1);
        REQUIRE(count3 == 1);
        REQUIRE(count4 == 2);
        found1 = found2 = found3 = found4 = found5 = found6 = false;
        count1 = count2 = count3 = count4 = 0;
        
        // now another node and check to make sure that the edges are updated appropriately

        h2 = graph.apply_orientation(graph.flip(h2));
        
        graph.follow_edges(h1, false, [&](const handle_t& other) {
            if (other == h2) {
                found1 = true;
            }
            else if (other == graph.flip(h2)) {
                found2 = true;
            }
            count1++;
        });
        graph.follow_edges(h1, true, [&](const handle_t& other) {
            if (other == h2) {
                found3 = true;
            }
            count2++;
        });
        graph.follow_edges(h2, false, [&](const handle_t& other) {
            if (other == h1) {
                found4 = true;
            }
            else if (other == graph.flip(h1)) {
                found5 = true;
            }
            count3++;
        });
        graph.follow_edges(h2, true, [&](const handle_t& other) {
            if (other == h1) {
                found6 = true;
            }
            count4++;
        });
        REQUIRE(found1);
        REQUIRE(found2);
        REQUIRE(found3);
        REQUIRE(found4);
        REQUIRE(found5);
        REQUIRE(found6);
        REQUIRE(count1 == 2);
        REQUIRE(count2 == 1);
        REQUIRE(count3 == 2);
        REQUIRE(count4 == 1);
    }
}
    
TEST_CASE("Mutable handle graphs with mutable paths work", "[handle][packed][hashgraph]") {
    
    vector<MutablePathDeletableHandleGraph*> implementations;
    
    bdsg::PackedGraph pg;
    implementations.push_back(&pg);
    
    bdsg::HashGraph hg;
    implementations.push_back(&hg);
    
    // These tests include assertions that embedded paths are maintained to be
    // consistent during graph edits, which VG doesn't currently do
    //VG vg;
    //implementations.push_back(&vg);
    
    for(MutablePathDeletableHandleGraph* implementation : implementations) {
        
        MutablePathDeletableHandleGraph& graph = *implementation;
        
        auto check_path = [&](const path_handle_t& p, const vector<handle_t>& steps) {
            
            REQUIRE(graph.get_step_count(p) == steps.size());
            
            step_handle_t step = graph.path_begin(p);
            for (int i = 0; i < steps.size(); i++) {
                
                REQUIRE(graph.get_path_handle_of_step(step) == p);
                REQUIRE(graph.get_handle_of_step(step) == steps[i]);
                
                if (graph.get_is_circular(p)) {
                    REQUIRE(graph.has_next_step(step));
                    REQUIRE(graph.has_previous_step(step));
                }
                else {
                    REQUIRE(graph.has_next_step(step) == i + 1 < steps.size());
                    REQUIRE(graph.has_previous_step(step) == i > 0);
                }
                
                step = graph.get_next_step(step);
            }
            
            if (graph.get_is_circular(p) && !graph.is_empty(p)) {
                REQUIRE(step == graph.path_begin(p));
            }
            else {
                REQUIRE(step == graph.path_end(p));
            }
            
            step = graph.path_back(p);
            
            for (int i = steps.size() - 1; i >= 0; i--) {
                
                REQUIRE(graph.get_path_handle_of_step(step) == p);
                REQUIRE(graph.get_handle_of_step(step) == steps[i]);
                
                if (graph.get_is_circular(p)) {
                    REQUIRE(graph.has_next_step(step));
                    REQUIRE(graph.has_previous_step(step));
                }
                else {
                    REQUIRE(graph.has_next_step(step) == i + 1 < steps.size());
                    REQUIRE(graph.has_previous_step(step) == i > 0);
                }
                
                step = graph.get_previous_step(step);
            }
            
            if (graph.get_is_circular(p) && !graph.is_empty(p)) {
                REQUIRE(step == graph.path_back(p));
            }
            else {
                REQUIRE(step == graph.path_front_end(p));
            }
        };
        
        auto check_flips = [&](const path_handle_t& p, const vector<handle_t>& steps) {
            auto flipped = steps;
            for (size_t i = 0; i < steps.size(); i++) {
                graph.apply_orientation(graph.flip(graph.forward(flipped[i])));
                flipped[i] = graph.flip(flipped[i]);
                check_path(p, flipped);
                
                graph.apply_orientation(graph.flip(graph.forward(flipped[i])));
                flipped[i] = graph.flip(flipped[i]);
                check_path(p, flipped);
            }
        };
        
        handle_t h1 = graph.create_handle("AC");
        handle_t h2 = graph.create_handle("CAGTGA");
        handle_t h3 = graph.create_handle("GT");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h1, graph.flip(h2));
        graph.create_edge(graph.flip(h2), h3);
        
        REQUIRE(!graph.has_path("1"));
        REQUIRE(graph.get_path_count() == 0);
        
        path_handle_t p1 = graph.create_path_handle("1");
        
        REQUIRE(graph.has_path("1"));
        REQUIRE(graph.get_path_count() == 1);
        REQUIRE(graph.get_path_handle("1") == p1);
        REQUIRE(graph.get_path_name(p1) == "1");
        REQUIRE(graph.get_step_count(p1) == 0);
        REQUIRE(graph.is_empty(p1));
        
        graph.append_step(p1, h1);
        
        REQUIRE(graph.get_step_count(p1) == 1);
        REQUIRE(!graph.is_empty(p1));
        
        graph.append_step(p1, h2);
        graph.append_step(p1, h3);
        
        REQUIRE(graph.get_step_count(p1) == 3);
        
        // graph can traverse a path
        check_path(p1, {h1, h2, h3});
        
        // graph preserves paths when reversing nodes
        check_flips(p1, {h1, h2, h3});
        
        // make a circular path
        path_handle_t p2 = graph.create_path_handle("2", true);
        REQUIRE(graph.get_path_count() == 2);
        
        graph.append_step(p2, h1);
        graph.append_step(p2, graph.flip(h2));
        graph.append_step(p2, h3);
        
        check_path(p2, {h1, graph.flip(h2), h3});
        
        // graph can query steps of a node on paths
        
        bool found1 = false, found2 = false;
        vector<step_handle_t> steps = graph.steps_of_handle(h1);
        for (auto& step : steps) {
            if (graph.get_path_handle_of_step(step) == p1 &&
                graph.get_handle_of_step(step) == h1) {
                found1 = true;
            }
            else if (graph.get_path_handle_of_step(step) == p2 &&
                     graph.get_handle_of_step(step) == h1) {
                found2 = true;
            }
            else {
                REQUIRE(false);
            }
        }
        REQUIRE(found1);
        REQUIRE(found2);
        found1 = found2 = false;
        
        steps = graph.steps_of_handle(h1, true);
        for (auto& step : steps) {
            if (graph.get_path_handle_of_step(step) == p1 &&
                graph.get_handle_of_step(step) == h1) {
                found1 = true;
            }
            else if (graph.get_path_handle_of_step(step) == p2 &&
                     graph.get_handle_of_step(step) == h1) {
                found2 = true;
            }
            else {
                REQUIRE(false);
            }
        }
        REQUIRE(found1);
        REQUIRE(found2);
        found1 = found2 = false;
        
        steps = graph.steps_of_handle(graph.flip(h1), true);
        for (auto& step : steps) {
            REQUIRE(false);
        }
        
        steps = graph.steps_of_handle(h2, true);
        for (auto& step : steps) {
            if (graph.get_path_handle_of_step(step) == p1 &&
                graph.get_handle_of_step(step) == h2) {
                found1 = true;
            }
            else {
                REQUIRE(false);
            }
        }
        steps = graph.steps_of_handle(graph.flip(h2), true);
        for (auto& step : steps) {
            if (graph.get_path_handle_of_step(step) == p2 &&
                graph.get_handle_of_step(step) == graph.flip(h2)) {
                found2 = true;
            }
            else {
                REQUIRE(false);
            }
        }
        REQUIRE(found1);
        REQUIRE(found2);
        found1 = found2 = false;
        
        vector<handle_t> segments = graph.divide_handle(h2, {size_t(2), size_t(4)});
        
        // graph preserves paths when dividing nodes
        
        check_path(p1, {h1, segments[0], segments[1], segments[2], h3});
        check_path(p2, {h1, graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0]), h3});
        
        path_handle_t p3 = graph.create_path_handle("3");
        graph.append_step(p3, h1);
        graph.append_step(p3, segments[0]);
        
        REQUIRE(graph.has_path("3"));
        REQUIRE(graph.get_path_count() == 3);
        
        // graph can toggle circularity
        
        graph.for_each_path_handle([&](const path_handle_t& p) {
            
            vector<handle_t> steps;
            
            for (handle_t h : graph.scan_path(p)) {
                steps.push_back(h);
            }
            
            bool starting_circularity = graph.get_is_circular(p);
            
            // make every transition occur
            for (bool circularity : {true, true, false, false, true}) {
                graph.set_circularity(p, circularity);
                REQUIRE(graph.get_is_circular(p) == circularity);
                check_path(p, steps);
            }
            
            graph.set_circularity(p, starting_circularity);
        });
        
        // graph can destroy paths
        
        graph.destroy_path(p3);
        
        REQUIRE(!graph.has_path("3"));
        REQUIRE(graph.get_path_count() == 2);
        
        bool found3 = false;
        
        graph.for_each_path_handle([&](const path_handle_t& p) {
            if (graph.get_path_name(p) == "1") {
                found1 = true;
            }
            else if (graph.get_path_name(p) == "2") {
                found2 = true;
            }
            else if (graph.get_path_name(p) == "3") {
                found3 = true;
            }
            else {
                REQUIRE(false);
            }
        });
        
        REQUIRE(found1);
        REQUIRE(found2);
        REQUIRE(!found3);
        
        // check flips to see if membership records are still functional
        check_flips(p1, {h1, segments[0], segments[1], segments[2], h3});
        check_flips(p2, {h1, graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0]), h3});
        
        graph.destroy_path(p1);
        
        REQUIRE(!graph.has_path("1"));
        REQUIRE(graph.get_path_count() == 1);
        
        found1 = found2 = found3 = false;
        
        graph.for_each_path_handle([&](const path_handle_t& p) {
            if (graph.get_path_name(p) == "1") {
                found1 = true;
            }
            else if (graph.get_path_name(p) == "2") {
                found2 = true;
            }
            else if (graph.get_path_name(p) == "3") {
                found3 = true;
            }
            else {
                REQUIRE(false);
            }
        });
        
        REQUIRE(!found1);
        REQUIRE(found2);
        REQUIRE(!found3);
        
        // check flips to see if membership records are still functional
        check_flips(p2, {h1, graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0]), h3});
        
        // make a path to rewrite
        path_handle_t p4 = graph.create_path_handle("4");
        graph.prepend_step(p4, h3);
        graph.prepend_step(p4, segments[2]);
        graph.prepend_step(p4, segments[1]);
        graph.prepend_step(p4, segments[0]);
        graph.prepend_step(p4, h1);
        
        check_flips(p4, {h1, segments[0], segments[1], segments[2], h3});
        
        auto check_rewritten_segment = [&](const pair<step_handle_t, step_handle_t>& new_segment,
                                           const vector<handle_t>& steps) {
            int i = 0;
            for (auto step = new_segment.first; step != new_segment.second; step = graph.get_next_step(step)) {
                REQUIRE(graph.get_handle_of_step(step) == steps[i]);
                i++;
            }
            REQUIRE(i == steps.size());
        };
        
        // rewrite the middle portion of a path
        
        step_handle_t s1 = graph.get_next_step(graph.path_begin(p4));
        step_handle_t s2 = graph.get_next_step(graph.get_next_step(graph.get_next_step(s1)));
        
        auto new_segment = graph.rewrite_segment(s1, s2, {graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0])});
        
        check_flips(p4, {h1, graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0]), h3});
        check_rewritten_segment(new_segment, {graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0])});
        
        // rewrite around the end of a circular path to delete
        
        graph.create_edge(h3, h1);
        graph.create_edge(segments[2], segments[0]);
        graph.set_circularity(p4, true);
        
        s1 = graph.get_previous_step(graph.path_begin(p4));
        s2 = graph.get_next_step(graph.path_begin(p4));
        
        new_segment = graph.rewrite_segment(s1, s2, vector<handle_t>());
        
        check_flips(p4, {graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0])});
        check_rewritten_segment(new_segment, vector<handle_t>());
        
        // add into an empty slot
        
        new_segment = graph.rewrite_segment(new_segment.first, new_segment.second, {graph.flip(h1), graph.flip(h3)});
        
        check_flips(p4, {graph.flip(h1), graph.flip(h3), graph.flip(segments[2]), graph.flip(segments[1]), graph.flip(segments[0])});
        check_rewritten_segment(new_segment, {graph.flip(h1), graph.flip(h3)});
        
    }
    
}

TEST_CASE("handlegraph PathMetadata name format preserves ranges on generic paths", "[handle]") {

    std::string path_name = PathMetadata::create_path_name(
        PathSense::GENERIC,
        PathMetadata::NO_SAMPLE_NAME,
        "randompath",
        PathMetadata::NO_HAPLOTYPE,
        PathMetadata::NO_PHASE_BLOCK,
        {10, PathMetadata::NO_END_POSITION}
    );
    
    REQUIRE(PathMetadata::parse_sense(path_name) == PathSense::GENERIC);
    REQUIRE(PathMetadata::parse_sample_name(path_name) == PathMetadata::NO_SAMPLE_NAME);
    REQUIRE(PathMetadata::parse_locus_name(path_name) == "randompath");
    REQUIRE(PathMetadata::parse_haplotype(path_name) == PathMetadata::NO_HAPLOTYPE);
    REQUIRE(PathMetadata::parse_phase_block(path_name) == PathMetadata::NO_PHASE_BLOCK);
    auto subrange = PathMetadata::parse_subrange(path_name);
    REQUIRE(subrange.first == 10);
    REQUIRE(subrange.second == PathMetadata::NO_END_POSITION);
}

TEST_CASE("handlegraph PathMetadata name format can parse two-part names", "[handle]") {
    std::string path_name = "GRCh38#chr1";

    PathSense sense;
    string sample;
    string locus;
    size_t haplotype;
    size_t phase_block;
    subrange_t subrange;
    PathMetadata::parse_path_name(path_name,
                                  sense,
                                  sample,
                                  locus,
                                  haplotype,
                                  phase_block,
                                  subrange);

    REQUIRE(PathMetadata::parse_sample_name(path_name) == sample);
    REQUIRE(PathMetadata::parse_locus_name(path_name) == locus);

    REQUIRE(sense == PathSense::REFERENCE); 
    REQUIRE(sample == "GRCh38");
    REQUIRE(locus == "chr1");
    REQUIRE(phase_block == PathMetadata::NO_PHASE_BLOCK);
    REQUIRE(subrange == PathMetadata::NO_SUBRANGE);
}

}
}
