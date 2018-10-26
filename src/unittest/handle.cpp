/**
 * \file 
 * unittest/handle.cpp: test cases for the implementations of the HandleGraph class.
 */

#include "catch.hpp"

#include "../handle.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../json2pb.h"

#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <unordered_set>

namespace vg {
namespace unittest {

using namespace std;

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
    
    SECTION("Occurrence handle equality works") {
        vector<occurrence_handle_t> handles;
        
        for (size_t i = 0; i < 10; i++) {
            for (size_t j = 0; j < 10; j++) {
                occurrence_handle_t handle;
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
    xg::XG xg_index(vg.graph);
    
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

}

TEST_CASE("Mutable handle graphs work", "[handle][vg]") {
    
    vector<MutableHandleGraph*> implementations;
    
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
                
                SECTION("Nodes can be swapped") {
                    // Get all the nodes
                    vector<handle_t> order;
                    g->for_each_handle([&](const handle_t& found) {
                        order.push_back(found);
                    });
                    REQUIRE(order.size() == 2);
                    
                    // Swap the two
                    g->swap_handles(order.front(), order.back());
                    
                    // Get all the nodes again
                    vector<handle_t> swapped;
                    g->for_each_handle([&](const handle_t& found) {
                        swapped.push_back(found);
                    });
                    REQUIRE(swapped.size() == 2);
                    
                    // Make sure they are in the opposite order when iterated
                    // after being swapped.
                    REQUIRE(swapped.front() == order.back());
                    REQUIRE(swapped.back() == order.front());
                    
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
    
    xg::XG xg_index(vg.graph);
    
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
    
    SECTION("Handles can traverse paths") {
        
        // check that a path is correctly accessible and traversible
        auto check_path_traversal = [&](const PathHandleGraph& graph, const Path& path) {
            
            path_handle_t path_handle = graph.get_path_handle(path.name());
            
            REQUIRE(graph.get_occurrence_count(path_handle) == path.mapping_size());
            
            // check that occurrence is pointing to the same index along the path
            auto check_occurrence = [&](const occurrence_handle_t& occurrence_handle,
                                        int mapping_idx) {
                
                REQUIRE(graph.get_path_handle_of_occurrence(occurrence_handle) == path_handle);
                
                const Mapping& mapping = path.mapping(mapping_idx);
                
                handle_t handle = graph.get_occurrence(occurrence_handle);
                
                REQUIRE(graph.get_id(handle) == mapping.position().node_id());
                REQUIRE(graph.get_is_reverse(handle) == mapping.position().is_reverse());
                REQUIRE(graph.get_ordinal_rank_of_occurrence(occurrence_handle) == mapping_idx);
            };
            
            occurrence_handle_t occurrence_handle;
            
            // iterate front to back
            occurrence_handle = graph.get_first_occurrence(path_handle);
            for (int i = 0; i < path.mapping_size(); i++) {
                if (i + 1 < path.mapping_size()) {
                    REQUIRE(graph.has_next_occurrence(occurrence_handle));
                }
                else {
                    REQUIRE(!graph.has_next_occurrence(occurrence_handle));
                }
                
                if (i > 0) {
                    REQUIRE(graph.has_previous_occurrence(occurrence_handle));
                }
                else {
                    REQUIRE(!graph.has_previous_occurrence(occurrence_handle));
                }
                
                check_occurrence(occurrence_handle, i);
                occurrence_handle = graph.get_next_occurrence(occurrence_handle);
            }

            // iterate front to back with a while
            {
                occurrence_handle = graph.get_first_occurrence(path_handle);
                int i = 0;
                check_occurrence(occurrence_handle, i);
                i++;
                while(graph.has_next_occurrence(occurrence_handle)) {
                    occurrence_handle = graph.get_next_occurrence(occurrence_handle);
                    check_occurrence(occurrence_handle, i);
                    i++;
                }
                REQUIRE(i == path.mapping_size());
            }

            // iterate front to back with the iteration function
            {
                int i = 0;
                graph.for_each_occurrence_in_path(path_handle, [&i, &check_occurrence](const occurrence_handle_t& occurrence_handle) {
                    check_occurrence(occurrence_handle, i);
                    i++;
                });
                REQUIRE(i == path.mapping_size());
            }

            // iterate back to front
            occurrence_handle = graph.get_last_occurrence(path_handle);
            for (int i = path.mapping_size() - 1; i >= 0; i--) {
                if (i + 1 < path.mapping_size()) {
                    REQUIRE(graph.has_next_occurrence(occurrence_handle));
                }
                else {
                    REQUIRE(!graph.has_next_occurrence(occurrence_handle));
                }
                
                if (i > 0) {
                    REQUIRE(graph.has_previous_occurrence(occurrence_handle));
                }
                else {
                    REQUIRE(!graph.has_previous_occurrence(occurrence_handle));
                }
                
                check_occurrence(occurrence_handle, i);
                occurrence_handle = graph.get_previous_occurrence(occurrence_handle);
            }

            // iterate back to front with a while
            {
                occurrence_handle = graph.get_last_occurrence(path_handle);
                int i = path.mapping_size() - 1;
                check_occurrence(occurrence_handle, i);
                i--;
                while(graph.has_previous_occurrence(occurrence_handle)) {
                    occurrence_handle = graph.get_previous_occurrence(occurrence_handle);
                    check_occurrence(occurrence_handle, i);
                    i--;
                }
                REQUIRE(i == -1);
            }
        };
        
        check_path_traversal(vg, path1);
        check_path_traversal(vg, path2);
        check_path_traversal(vg, path3);
        check_path_traversal(xg_index, path1);
        check_path_traversal(xg_index, path2);
        check_path_traversal(xg_index, path3);
    }
}

}
}
