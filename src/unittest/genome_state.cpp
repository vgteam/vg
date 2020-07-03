//
//  genome_state.cpp
//
//  Unit tests for SnarlState, GenomeState and related functions
//

#include <stdio.h>
#include <iostream>
#include <set>
#include "catch.hpp"
#include "vg/io/json2pb.h"
#include "../genotypekit.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../genome_state.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("SnarlState can hold and manipulate haplotypes", "[snarlstate][genomestate]") {
    

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
    
    // Work out its snarls
    CactusSnarlFinder bubble_finder(graph);
    SnarlManager snarl_manager = bubble_finder.find_snarls();
    
    // Get the top snarl
    const Snarl* top_snarl = snarl_manager.top_level_snarls().at(0);
    
    // Make sure it's what we expect.
    REQUIRE(top_snarl->start().node_id() == 1);
    REQUIRE(top_snarl->end().node_id() == 8);
    
    // And get its net graph
    NetGraph net_graph = snarl_manager.net_graph_of(top_snarl, &graph, true);
    
    // Make a SnarlState for it
    SnarlState state(&net_graph);
    
    SECTION("The state starts empty") {
        REQUIRE(state.size() == 0);
    }
    
    SECTION("A haplotype with lane numbers can be added in lane 0 when the SnarlState is empty") {
        // Make a haplotype
        vector<pair<handle_t, size_t>> annotated_haplotype;
        
        // Say we go 1, 2 (which is a child chain), 8
        annotated_haplotype.emplace_back(net_graph.get_handle(1, false), 0);
        annotated_haplotype.emplace_back(net_graph.get_handle(2, false), 0);
        annotated_haplotype.emplace_back(net_graph.get_handle(8, false), 0);
        
        // Put it in the state
        state.insert(annotated_haplotype);
        
        SECTION("The state now has 1 haplotype") {
            REQUIRE(state.size() == 1);
        }
        
        SECTION("The haplotype can be traced back again") {
            vector<pair<handle_t, size_t>> recovered;
            
            state.trace(0, false, [&](const handle_t& visit, size_t local_lane) {
                recovered.emplace_back(visit, local_lane);
            });
            
            REQUIRE(recovered == annotated_haplotype);
        }
        
        SECTION("The haplotype can be traced in reverse") {
            vector<pair<handle_t, size_t>> recovered;
            
            state.trace(0, true, [&](const handle_t& visit, size_t local_lane) {
                recovered.emplace_back(visit, local_lane);
            });
            
            REQUIRE(recovered.size() == 3);
            REQUIRE(recovered[0].first == net_graph.get_handle(8, true));
            REQUIRE(recovered[0].second == 0);
            // We use the snarl's leading node in reverse for the snarl in
            // reverse, even though we entered the end of the snarl and then
            // skipped tot he front.
            REQUIRE(recovered[1].first == net_graph.get_handle(2, true));
            REQUIRE(recovered[1].second == 0);
            REQUIRE(recovered[2].first == net_graph.get_handle(1, true));
            REQUIRE(recovered[2].second == 0);
        }
        
        SECTION("A haplotype with lane numbers can be inserted before an existing haplotype") {
            
            vector<pair<handle_t, size_t>> hap2;
            // Say we go 1, 8 directly in lane 0
            hap2.emplace_back(net_graph.get_handle(1, false), 0);
            hap2.emplace_back(net_graph.get_handle(8, false), 0);
            
            // Put it in the state
            state.insert(hap2);
            
            SECTION("The state now has 2 haplotypes") {
                REQUIRE(state.size() == 2);
            }
            
            SECTION("The new haplotype can be traced back again") {
                vector<pair<handle_t, size_t>> recovered;
                
                state.trace(0, false, [&](const handle_t& visit, size_t local_lane) {
                    recovered.emplace_back(visit, local_lane);
                });
                
                REQUIRE(recovered == hap2);
            }
            
            SECTION("The old haplotype can be traced back again") {
                vector<pair<handle_t, size_t>> recovered;
                
                state.trace(1, false, [&](const handle_t& visit, size_t local_lane) {
                    recovered.emplace_back(visit, local_lane);
                });
                
                REQUIRE(recovered.size() == 3);
                REQUIRE(recovered[0].first == annotated_haplotype[0].first);
                REQUIRE(recovered[0].second == annotated_haplotype[0].second + 1);
                // The second mapping should not get bumped up.
                REQUIRE(recovered[1].first == annotated_haplotype[1].first);
                REQUIRE(recovered[1].second == annotated_haplotype[1].second);
                REQUIRE(recovered[2].first == annotated_haplotype[2].first);
                REQUIRE(recovered[2].second == annotated_haplotype[2].second + 1);
                
            }
            
            SECTION("A haplotype without lane numbers can be appended") {
                // Make a haplotype
                vector<handle_t> hap3;
                hap3.emplace_back(net_graph.get_handle(1, false));
                hap3.emplace_back(net_graph.get_handle(2, false));
                hap3.emplace_back(net_graph.get_handle(8, false));
                
                // Put it in the state
                auto added = state.append(hap3);
                
                SECTION("The state now has 3 haplotypes") {
                    REQUIRE(state.size() == 3);
                }
                
                SECTION("The returned annotated haplotype is correct") {
                    REQUIRE(added.size() == 3);
                    REQUIRE(added[0].first == hap3[0]);
                    REQUIRE(added[0].second == 2);
                    REQUIRE(added[1].first == hap3[1]);
                    REQUIRE(added[1].second == 1);
                    REQUIRE(added[2].first == hap3[2]);
                    REQUIRE(added[2].second == 2);
                }
                
                SECTION("The new haplotype can be traced back again") {
                    vector<pair<handle_t, size_t>> recovered;
                    
                    state.trace(2, false, [&](const handle_t& visit, size_t local_lane) {
                        recovered.emplace_back(visit, local_lane);
                    });
                    
                    REQUIRE(recovered == added);
                }
                
                SECTION("It can be deleted again") {
                    state.erase(2);
                    
                    REQUIRE(state.size() == 2);
                        
                    SECTION("The second haplotype can be traced back again") {
                        vector<pair<handle_t, size_t>> recovered;
                        
                        state.trace(0, false, [&](const handle_t& visit, size_t local_lane) {
                            recovered.emplace_back(visit, local_lane);
                        });
                        
                        REQUIRE(recovered == hap2);
                    }
                    
                    SECTION("The first haplotype can be traced back again") {
                        vector<pair<handle_t, size_t>> recovered;
                        
                        state.trace(1, false, [&](const handle_t& visit, size_t local_lane) {
                            recovered.emplace_back(visit, local_lane);
                        });
                        
                        REQUIRE(recovered.size() == 3);
                        REQUIRE(recovered[0].first == annotated_haplotype[0].first);
                        REQUIRE(recovered[0].second == annotated_haplotype[0].second + 1);
                        // The second mapping should not get bumped up.
                        REQUIRE(recovered[1].first == annotated_haplotype[1].first);
                        REQUIRE(recovered[1].second == annotated_haplotype[1].second);
                        REQUIRE(recovered[2].first == annotated_haplotype[2].first);
                        REQUIRE(recovered[2].second == annotated_haplotype[2].second + 1);
                        
                    }
                    
                }
                
                SECTION("It can be swapped with another haplotype") {
                    state.swap(0, 2);
                    
                    SECTION("Swapping haplotypes dows not change the overall size") {
                        REQUIRE(state.size() == 3);
                    }
                    
                    SECTION("The second haplotype added is now in overall lane 2") {
                        vector<pair<handle_t, size_t>> recovered;
                        
                        
                        state.trace(2, false, [&](const handle_t& visit, size_t local_lane) {
                            recovered.emplace_back(visit, local_lane);
                        });
                        
                        
                        
                        REQUIRE(recovered.size() == 2);
                        REQUIRE(recovered[0].first == hap2[0].first);
                        REQUIRE(recovered[0].second == 2);
                        REQUIRE(recovered[1].first == hap2[1].first);
                        REQUIRE(recovered[1].second == 2);
                    }
                    
                    SECTION("The third haplotype added is now in overall lane 0") {
                        vector<pair<handle_t, size_t>> recovered;
                        
                        state.trace(0, false, [&](const handle_t& visit, size_t local_lane) {
                            recovered.emplace_back(visit, local_lane);
                        });
                        
                        REQUIRE(recovered.size() == 3);
                        REQUIRE(recovered[0].first == hap3[0]);
                        REQUIRE(recovered[0].second == 0);
                        SECTION("Swapping haplotypes does not change interior child chain lane assignments") {
                            // The second mapping should stay in place in its assigned lane
                            REQUIRE(recovered[1].first == hap3[1]);
                            REQUIRE(recovered[1].second == 1);
                        }
                        REQUIRE(recovered[2].first == hap3[2]);
                        REQUIRE(recovered[2].second == 0);
                    }
                    
                    SECTION("The first haplotype added is unaffected in overall lane 1") {
                        vector<pair<handle_t, size_t>> recovered;
                        
                        state.trace(1, false, [&](const handle_t& visit, size_t local_lane) {
                            recovered.emplace_back(visit, local_lane);
                        });
                        
                        REQUIRE(recovered.size() == 3);
                        REQUIRE(recovered[0].first == annotated_haplotype[0].first);
                        REQUIRE(recovered[0].second == annotated_haplotype[0].second + 1);
                        // The second mapping should not get bumped up.
                        REQUIRE(recovered[1].first == annotated_haplotype[1].first);
                        REQUIRE(recovered[1].second == annotated_haplotype[1].second);
                        REQUIRE(recovered[2].first == annotated_haplotype[2].first);
                        REQUIRE(recovered[2].second == annotated_haplotype[2].second + 1);
                        
                    }
                    
                    
                }
                
            }
            
            SECTION("A haplotype without lane numbers can be inserted") {
                // Make a haplotype
                vector<handle_t> hap3;
                hap3.emplace_back(net_graph.get_handle(1, false));
                hap3.emplace_back(net_graph.get_handle(2, false));
                hap3.emplace_back(net_graph.get_handle(8, false));
                
                // Put it in the state, at lane 1
                auto added = state.insert(1, hap3);
                
                SECTION("The returned annotated haplotype is correct") {
                    REQUIRE(added.size() == 3);
                    REQUIRE(added[0].first == hap3[0]);
                    REQUIRE(added[0].second == 1);
                    REQUIRE(added[1].first == hap3[1]);
                    // We don't actually care about our middle node lane assignment.
                    // It could be before or after other haplotypes; they're allowed to cross over each other arbitrarily.
                    REQUIRE(added[1].second >= 0);
                    REQUIRE(added[1].second <= 1);
                    REQUIRE(added[2].first == hap3[2]);
                    REQUIRE(added[2].second == 1);
                }
                
                SECTION("The new haplotype can be traced back again") {
                    vector<pair<handle_t, size_t>> recovered;
                    
                    state.trace(1, false, [&](const handle_t& visit, size_t local_lane) {
                        recovered.emplace_back(visit, local_lane);
                    });
                    
                    REQUIRE(recovered == added);
                }
                
                SECTION("The bumped-up haplotype can be traced back again") {
                    vector<pair<handle_t, size_t>> recovered;
                    
                    state.trace(2, false, [&](const handle_t& visit, size_t local_lane) {
                        recovered.emplace_back(visit, local_lane);
                    });
                    
                    REQUIRE(recovered.size() == 3);
                    REQUIRE(net_graph.get_id(recovered[0].first) == 1);
                    REQUIRE(recovered[0].second == 2);
                    REQUIRE(net_graph.get_id(recovered[1].first) == 2);
                    // Lane assignment at the middle visit may or may not have been pushed up.
                    REQUIRE(recovered[1].second >= 0);
                    REQUIRE(recovered[1].second <= 1);
                    REQUIRE(recovered[1].second != added[1].second);
                    REQUIRE(net_graph.get_id(recovered[2].first) == 8);
                    REQUIRE(recovered[2].second == 2);
                }
                
            }
        
        }
        
    }
    
    SECTION("Haplotypes cannot be added in reverse") {
        // Make a haplotype
        vector<pair<handle_t, size_t>> annotated_haplotype;
        
        // Say we go 8 rev, 2 rev (which is a child chain), 1 rev
        annotated_haplotype.emplace_back(net_graph.get_handle(8, true), 0);
        annotated_haplotype.emplace_back(net_graph.get_handle(2, true), 0);
        annotated_haplotype.emplace_back(net_graph.get_handle(1, true), 0);
        
        // Try and fail to put it in the state
        REQUIRE_THROWS(state.insert(annotated_haplotype));
        
    }
}


TEST_CASE("SnarlState works on snarls with nontrivial child chains", "[snarlstate][genomestate]") {
    

    // This graph will have a snarl from 1 to 8, a snarl from 2 to 4,
    // and a snarl from 4 to 7, with a chain in the top snarl.
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
    Edge* e4 = graph.create_edge(n2, n4);
    Edge* e5 = graph.create_edge(n3, n4);
    Edge* e6 = graph.create_edge(n4, n5);
    Edge* e7 = graph.create_edge(n4, n6);
    Edge* e8 = graph.create_edge(n5, n7);
    Edge* e9 = graph.create_edge(n6, n7);
    Edge* e10 = graph.create_edge(n7, n8);
    
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
    
    if (left_child->end().node_id() < left_child->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(left_child);
    }
    
    if (right_child->end().node_id() < right_child->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(right_child);
    }
    
    if (left_child->start().node_id() > right_child->start().node_id()) {
        // Put them in a consistent orientation
        std::swap(left_child, right_child);
    }
    
    REQUIRE(left_child->start().node_id() == 2);
    REQUIRE(left_child->end().node_id() == 4);
    REQUIRE(right_child->start().node_id() == 4);
    REQUIRE(right_child->end().node_id() == 7);
    
    // And get its net graph
    NetGraph net_graph = snarl_manager.net_graph_of(top_snarl, &graph, true);
    
    // Make a SnarlState for the top snarl
    SnarlState state(&net_graph);
    
    SECTION("The state starts empty") {
        REQUIRE(state.size() == 0);
    }
    
    SECTION("We can add a traversal through the chain") {
        // Make a haplotype
        vector<pair<handle_t, size_t>> annotated_haplotype;
        
        // Say we go 1, 2 (which is a child chain), 8
        annotated_haplotype.emplace_back(net_graph.get_handle(1, false), 0);
        annotated_haplotype.emplace_back(net_graph.get_handle(2, false), 0);
        annotated_haplotype.emplace_back(net_graph.get_handle(8, false), 0);
        
        // Put it in the state
        state.insert(annotated_haplotype);
        
        SECTION("The state now has 1 haplotype") {
            REQUIRE(state.size() == 1);
        }
        
        SECTION("The haplotype can be traced back again") {
            vector<pair<handle_t, size_t>> recovered;
            
            state.trace(0, false, [&](const handle_t& visit, size_t local_lane) {
                recovered.emplace_back(visit, local_lane);
            });
            
            REQUIRE(recovered == annotated_haplotype);
        }
        
    }
}

TEST_CASE("GenomeState can hold and manipulate haplotypes", "[genomestate]") {
    

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
    
    // Work out its snarls
    CactusSnarlFinder bubble_finder(graph);
    SnarlManager snarl_manager = bubble_finder.find_snarls();

    // Get the top snarl
    const Snarl* top_snarl = snarl_manager.top_level_snarls().at(0);
    
    if (top_snarl->end().node_id() < top_snarl->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(top_snarl);
    }
    
    // Make sure it's what we expect.
    REQUIRE(top_snarl->start().node_id() == 1);
    REQUIRE(top_snarl->end().node_id() == 8);
    
    // Get the middle snarl
    const Snarl* middle_snarl = snarl_manager.children_of(top_snarl).at(0);
    
    if (middle_snarl->end().node_id() < middle_snarl->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(middle_snarl);
    }
    
    REQUIRE(middle_snarl->start().node_id() == 2);
    REQUIRE(middle_snarl->end().node_id() == 7);
    
    // And the bottom snarl
    const Snarl* bottom_snarl = snarl_manager.children_of(middle_snarl).at(0);
    
    if (bottom_snarl->end().node_id() < bottom_snarl->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(bottom_snarl);
    }
    
    REQUIRE(bottom_snarl->start().node_id() == 3);
    REQUIRE(bottom_snarl->end().node_id() == 5);
    
    // Define the chromosome by telomere snarls (first and last)
    auto chromosome = make_pair(top_snarl, top_snarl);
    
    // Make a genome state for this genome, with the only top snarl being the
    // telomere snarls for the main chromosome
    GenomeState state(snarl_manager, &graph, {chromosome});
    
    // Make sure to get all the net graphs
    const NetGraph* top_graph = state.get_net_graph(top_snarl);
    const NetGraph* middle_graph = state.get_net_graph(middle_snarl);
    const NetGraph* bottom_graph = state.get_net_graph(bottom_snarl);
    
    SECTION("NetGraphs for snarls can be obtained") {
        REQUIRE(top_graph != nullptr);
        REQUIRE(middle_graph != nullptr);
        REQUIRE(bottom_graph != nullptr);
        
#ifdef debug
        cerr << "Top graph " << top_snarl->start().node_id() << " - " << top_snarl->end().node_id() << " contents: " << endl;
        top_graph->for_each_handle([&](const handle_t h) {
            cerr << "\t" << top_graph->get_id(h) << endl;
        });
        
        cerr << "Middle graph " << middle_snarl->start().node_id() << " - " << middle_snarl->end().node_id() << " contents: " << endl;
        middle_graph->for_each_handle([&](const handle_t h) {
            cerr << "\t" << middle_graph->get_id(h) << endl;
        });
        
        cerr << "Bottom graph " << bottom_snarl->start().node_id() << " - " << bottom_snarl->end().node_id() << " contents: " << endl;
        bottom_graph->for_each_handle([&](const handle_t h) {
            cerr << "\t" << bottom_graph->get_id(h) << endl;
        });
#endif
    }
    
    SECTION("GenomeState starts empty") {
        REQUIRE(state.count_haplotypes(chromosome) == 0);
    }
    
    SECTION("A haplotype can be added") {
        // Define a haplotype across the entire graph, in three levels
        InsertHaplotypeCommand insert;
        
        // We fill these with handles from the appropriate net graphs.
        // TODO: really we could just use the backing graph. Would that be better???
        
        // For the top snarl we go 1, 2 (child), and 8
        insert.insertions.emplace(top_snarl, vector<vector<pair<handle_t, size_t>>>{{
            {top_graph->get_handle(1, false), 0},
            {top_graph->get_handle_from_inward_backing_handle(graph.get_handle(2, false)), 0},
            {top_graph->get_handle(8, false), 0}
        }});
        
        // For the middle snarl we go 2, 3 (child), 7
        insert.insertions.emplace(middle_snarl, vector<vector<pair<handle_t, size_t>>>{{
            {middle_graph->get_handle(2, false), 0},
            {middle_graph->get_handle_from_inward_backing_handle(graph.get_handle(3, false)), 0},
            {middle_graph->get_handle(7, false), 0}
        }});
        
        // For the bottom snarl we go 3, 5 (skipping over 4)
        insert.insertions.emplace(bottom_snarl, vector<vector<pair<handle_t, size_t>>>{{
            {bottom_graph->get_handle(3, false), 0},
            {bottom_graph->get_handle(5, false), 0}
        }});
        
        // Execute the command and get the undo command
        GenomeStateCommand* undo = state.execute(&insert);
        
        SECTION("The added haplotype is counted") {
            REQUIRE(state.count_haplotypes(chromosome) == 1);
        }
        
        SECTION("The haplotype can be traced") {
            // We trace out all the handles in the backing graph
            vector<handle_t> traced;
            
            state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                // Put every handle in the vector
                traced.push_back(visit);
            });
         
            // We should visit nodes 1, 2, 3, 5, 7, 8 in that order
            REQUIRE(traced.size() == 6);
            REQUIRE(traced[0] == graph.get_handle(1, false));
            REQUIRE(traced[1] == graph.get_handle(2, false));
            REQUIRE(traced[2] == graph.get_handle(3, false));
            REQUIRE(traced[3] == graph.get_handle(5, false));
            REQUIRE(traced[4] == graph.get_handle(7, false));
            REQUIRE(traced[5] == graph.get_handle(8, false));
        }
        
        SECTION("A second haplotype can be added") {
            // We're just going to steal lane 0 in everything so we know things match up.
        
            InsertHaplotypeCommand insert2;
        
            // For the top snarl we go 1, 2 (child), and 8
            insert2.insertions.emplace(top_snarl, vector<vector<pair<handle_t, size_t>>>{{
                {top_graph->get_handle(1, false), 0},
                {top_graph->get_handle_from_inward_backing_handle(graph.get_handle(2, false)), 0},
                {top_graph->get_handle(8, false), 0}
            }});
            
            // For the middle snarl we go 2, 3 (child), 7
            insert2.insertions.emplace(middle_snarl, vector<vector<pair<handle_t, size_t>>>{{
                {middle_graph->get_handle(2, false), 0},
                {middle_graph->get_handle_from_inward_backing_handle(graph.get_handle(3, false)), 0},
                {middle_graph->get_handle(7, false), 0}
            }});
            
            // For the bottom snarl we go 3, 4, 5
            insert2.insertions.emplace(bottom_snarl, vector<vector<pair<handle_t, size_t>>>{{
                {bottom_graph->get_handle(3, false), 0},
                {bottom_graph->get_handle(4, false), 0},
                {bottom_graph->get_handle(5, false), 0}
            }});
            
            // Execute the command and get the undo command
            GenomeStateCommand* undo2 = state.execute(&insert2);
            
            SECTION("The added haplotype is counted") {
                REQUIRE(state.count_haplotypes(chromosome) == 2);
            }
            
            SECTION("The new haplotype can be traced") {
                // We trace out all the handles in the backing graph
                vector<handle_t> traced;
                
                state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                    // Put every handle in the vector
                    traced.push_back(visit);
                });
             
                // We should visit nodes 1, 2, 3, 4, 5, 7, 8 in that order
                REQUIRE(traced.size() == 7);
                REQUIRE(traced[0] == graph.get_handle(1, false));
                REQUIRE(traced[1] == graph.get_handle(2, false));
                REQUIRE(traced[2] == graph.get_handle(3, false));
                REQUIRE(traced[3] == graph.get_handle(4, false));
                REQUIRE(traced[4] == graph.get_handle(5, false));
                REQUIRE(traced[5] == graph.get_handle(7, false));
                REQUIRE(traced[6] == graph.get_handle(8, false));
            }
            
            SECTION("The old haplotype can be traced") {
                // We trace out all the handles in the backing graph
                vector<handle_t> traced;
                
                state.trace_haplotype(chromosome, 1, [&](const handle_t& visit) {
                    // Put every handle in the vector
                    traced.push_back(visit);
                });
             
                // We should visit nodes 1, 2, 3, 5, 7, 8 in that order
                REQUIRE(traced.size() == 6);
                REQUIRE(traced[0] == graph.get_handle(1, false));
                REQUIRE(traced[1] == graph.get_handle(2, false));
                REQUIRE(traced[2] == graph.get_handle(3, false));
                REQUIRE(traced[3] == graph.get_handle(5, false));
                REQUIRE(traced[4] == graph.get_handle(7, false));
                REQUIRE(traced[5] == graph.get_handle(8, false));
            }
            
            SECTION("The two haplotypes can be swapped") {
                SwapHaplotypesCommand swapper;
                swapper.telomere_pair = chromosome;
                swapper.to_swap = make_pair(0, 1);
                
                GenomeStateCommand* unswapper = state.execute(&swapper);
                
                SECTION("There are still two haplotypes") {
                    REQUIRE(state.count_haplotypes(chromosome) == 2);
                }
                
                SECTION("The new haplotype can be traced in its new position") {
                    // We trace out all the handles in the backing graph
                    vector<handle_t> traced;
                    
                    state.trace_haplotype(chromosome, 1, [&](const handle_t& visit) {
                        // Put every handle in the vector
                        traced.push_back(visit);
                    });
                 
                    // We should visit nodes 1, 2, 3, 4, 5, 7, 8 in that order
                    REQUIRE(traced.size() == 7);
                    REQUIRE(traced[0] == graph.get_handle(1, false));
                    REQUIRE(traced[1] == graph.get_handle(2, false));
                    REQUIRE(traced[2] == graph.get_handle(3, false));
                    REQUIRE(traced[3] == graph.get_handle(4, false));
                    REQUIRE(traced[4] == graph.get_handle(5, false));
                    REQUIRE(traced[5] == graph.get_handle(7, false));
                    REQUIRE(traced[6] == graph.get_handle(8, false));
                }
                
                SECTION("The old haplotype can be traced in its new position") {
                    // We trace out all the handles in the backing graph
                    vector<handle_t> traced;
                    
                    state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                        // Put every handle in the vector
                        traced.push_back(visit);
                    });
                 
                    // We should visit nodes 1, 2, 3, 5, 7, 8 in that order
                    REQUIRE(traced.size() == 6);
                    REQUIRE(traced[0] == graph.get_handle(1, false));
                    REQUIRE(traced[1] == graph.get_handle(2, false));
                    REQUIRE(traced[2] == graph.get_handle(3, false));
                    REQUIRE(traced[3] == graph.get_handle(5, false));
                    REQUIRE(traced[4] == graph.get_handle(7, false));
                    REQUIRE(traced[5] == graph.get_handle(8, false));
                }
                
                SECTION("The swap can be undone") {
                    
                    delete state.execute(unswapper);
                    
                    SECTION("There are still two haplotypes") {
                        REQUIRE(state.count_haplotypes(chromosome) == 2);
                    }
                    
                    SECTION("The new haplotype can be traced in its original position") {
                        // We trace out all the handles in the backing graph
                        vector<handle_t> traced;
                        
                        state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                            // Put every handle in the vector
                            traced.push_back(visit);
                        });
                     
                        // We should visit nodes 1, 2, 3, 4, 5, 7, 8 in that order
                        REQUIRE(traced.size() == 7);
                        REQUIRE(traced[0] == graph.get_handle(1, false));
                        REQUIRE(traced[1] == graph.get_handle(2, false));
                        REQUIRE(traced[2] == graph.get_handle(3, false));
                        REQUIRE(traced[3] == graph.get_handle(4, false));
                        REQUIRE(traced[4] == graph.get_handle(5, false));
                        REQUIRE(traced[5] == graph.get_handle(7, false));
                        REQUIRE(traced[6] == graph.get_handle(8, false));
                    }
                    
                    SECTION("The old haplotype can be traced in its original position") {
                        // We trace out all the handles in the backing graph
                        vector<handle_t> traced;
                        
                        state.trace_haplotype(chromosome, 1, [&](const handle_t& visit) {
                            // Put every handle in the vector
                            traced.push_back(visit);
                        });
                     
                        // We should visit nodes 1, 2, 3, 5, 7, 8 in that order
                        REQUIRE(traced.size() == 6);
                        REQUIRE(traced[0] == graph.get_handle(1, false));
                        REQUIRE(traced[1] == graph.get_handle(2, false));
                        REQUIRE(traced[2] == graph.get_handle(3, false));
                        REQUIRE(traced[3] == graph.get_handle(5, false));
                        REQUIRE(traced[4] == graph.get_handle(7, false));
                        REQUIRE(traced[5] == graph.get_handle(8, false));
                    }               
                    
                }
                
                delete unswapper;
                
            }
            
            // We can't delete the old haplotype with its original uninsert
            // command because now it has been moved. We have to do the
            // deletions in reverse insertion order.
            
            SECTION("The new haplotype can be deleted") {
                delete state.execute(undo2);
                
                SECTION("Only one haplotype is left") {
                    REQUIRE(state.count_haplotypes(chromosome) == 1);
                }
                
                SECTION("The remaining, original haplotype can be traced") {
                    // We trace out all the handles in the backing graph
                    vector<handle_t> traced;
                    
                    state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                        // Put every handle in the vector
                        traced.push_back(visit);
                    });
                 
                    // We should visit nodes 1, 2, 3, 5, 7, 8 in that order
                    REQUIRE(traced.size() == 6);
                    REQUIRE(traced[0] == graph.get_handle(1, false));
                    REQUIRE(traced[1] == graph.get_handle(2, false));
                    REQUIRE(traced[2] == graph.get_handle(3, false));
                    REQUIRE(traced[3] == graph.get_handle(5, false));
                    REQUIRE(traced[4] == graph.get_handle(7, false));
                    REQUIRE(traced[5] == graph.get_handle(8, false));
                }
                
                SECTION("The old haplotype can also be deleted") {
                    delete state.execute(undo);
                    
                    SECTION("No haplotypes are left") {
                        REQUIRE(state.count_haplotypes(chromosome) == 0);
                    }
                    
                }
                
            }
            
            delete undo2;
            
        }
        
        SECTION("The added haplotype can be deleted again") {
            GenomeStateCommand* undelete = state.execute(undo);
            
            REQUIRE(state.count_haplotypes(chromosome) == 0);
        
            SECTION("The undelete command matches the insert command") {
                REQUIRE(((InsertHaplotypeCommand*)undelete)->insertions == insert.insertions);
            }
        
            SECTION("And then it can be un-deleted") {
                
                GenomeStateCommand* redelete = state.execute(undelete);
                
                REQUIRE(state.count_haplotypes(chromosome) == 1);
                
                delete redelete;
            
            }
        
            delete undelete;
        
        }
        
        // Free the undo command
        delete undo;
    }
    
}

TEST_CASE("GenomeSate works on snarls with nontrivial child chains", "[genomestate][broken]") {
    

    // This graph will have a snarl from 1 to 8, a snarl from 2 to 4,
    // and a snarl from 4 to 7, with a chain in the top snarl.
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
    Edge* e4 = graph.create_edge(n2, n4);
    Edge* e5 = graph.create_edge(n3, n4);
    Edge* e6 = graph.create_edge(n4, n5);
    Edge* e7 = graph.create_edge(n4, n6);
    Edge* e8 = graph.create_edge(n5, n7);
    Edge* e9 = graph.create_edge(n6, n7);
    Edge* e10 = graph.create_edge(n7, n8);
    
    // Work out its snarls
    CactusSnarlFinder bubble_finder(graph);
    SnarlManager snarl_manager = bubble_finder.find_snarls();
    
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
    
    if (left_child->end().node_id() < left_child->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(left_child);
    }
    
    if (right_child->end().node_id() < right_child->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(right_child);
    }
    
    if (left_child->start().node_id() > right_child->start().node_id()) {
        // Put them in a consistent orientation
        std::swap(left_child, right_child);
    }
    
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
    
    REQUIRE(left_child->start().node_id() == 2);
    REQUIRE(left_child->end().node_id() == 4);
    REQUIRE(right_child->start().node_id() == 4);
    REQUIRE(right_child->end().node_id() == 7);
    
    // Define the chromosome by telomere snarls (first and last)
    auto chromosome = make_pair(top_snarl, top_snarl);
    
    // Make a genome state for this genome, with the only top snarl being the
    // telomere snarls for the main chromosome
    GenomeState state(snarl_manager, &graph, {chromosome});
    
    // Make sure to get all the net graphs
    const NetGraph* top_graph = state.get_net_graph(top_snarl);
    const NetGraph* left_graph = state.get_net_graph(left_child);
    const NetGraph* right_graph = state.get_net_graph(right_child);
    
    SECTION("NetGraphs for snarls can be obtained") {
        REQUIRE(top_graph != nullptr);
        REQUIRE(left_graph != nullptr);
        REQUIRE(right_graph != nullptr);
    }
    
    SECTION("GenomeState starts empty") {
        REQUIRE(state.count_haplotypes(chromosome) == 0);
    }
    
    SECTION("A haplotype can be added") {
        // Define a haplotype across the entire graph, in three levels
        InsertHaplotypeCommand insert;
        
        // We fill these with handles from the appropriate net graphs.
        // TODO: really we could just use the backing graph. Would that be better???
        
        // For the top snarl we go 1, 2 (child chain), and 8
        insert.insertions.emplace(top_snarl, vector<vector<pair<handle_t, size_t>>>{{
            {top_graph->get_handle(1, false), 0},
            {top_graph->get_handle_from_inward_backing_handle(graph.get_handle(2, false)), 0},
            {top_graph->get_handle(8, false), 0}
        }});
        
        // For the left child snarl we go 2, 3, 4
        insert.insertions.emplace(left_child, vector<vector<pair<handle_t, size_t>>>{{
            {left_graph->get_handle(2, false), 0},
            {left_graph->get_handle(3, false), 0},
            {left_graph->get_handle(4, false), 0}
        }});
        
        // For the right child snarl we go 4, 5, 7
        insert.insertions.emplace(right_child, vector<vector<pair<handle_t, size_t>>>{{
            {right_graph->get_handle(4, false), 0},
            {right_graph->get_handle(5, false), 0},
            {right_graph->get_handle(7, false), 0}
        }});
        
        // Execute the command and get the undo command
        GenomeStateCommand* undo = state.execute(&insert);
        
        SECTION("The added haplotype is counted") {
            REQUIRE(state.count_haplotypes(chromosome) == 1);
        }
        
        SECTION("The haplotype can be traced") {
            // We trace out all the handles in the backing graph
            vector<handle_t> traced;
            
            state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                // Put every handle in the vector
                traced.push_back(visit);
            });
         
            // We should visit nodes 1, 2, 3, 4, 5, 7, 8 in that order
            REQUIRE(traced.size() == 7);
            REQUIRE(traced[0] == graph.get_handle(1, false));
            REQUIRE(traced[1] == graph.get_handle(2, false));
            REQUIRE(traced[2] == graph.get_handle(3, false));
            REQUIRE(traced[3] == graph.get_handle(4, false));
            REQUIRE(traced[4] == graph.get_handle(5, false));
            REQUIRE(traced[5] == graph.get_handle(7, false));
            REQUIRE(traced[6] == graph.get_handle(8, false));
        }
        
        delete undo;
    }
    
}

TEST_CASE("GenomeSate works on snarls with nontrivial child chains with backward snarls", "[genomestate]") {
    

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
    
    if (left_child->end().node_id() < left_child->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(left_child);
    }
    
    if (right_child->end().node_id() < right_child->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(right_child);
    }
    
    if (left_child->start().node_id() > right_child->start().node_id()) {
        // Put them in a consistent orientation
        std::swap(left_child, right_child);
    }
    
    REQUIRE(left_child->start().node_id() == 2);
    REQUIRE(left_child->end().node_id() == 4);
    
    // Make sure the right child is BACKWARD in the chain
    snarl_manager.flip(right_child);
    
    REQUIRE(right_child->start().node_id() == 7);
    REQUIRE(right_child->end().node_id() == 4);
    
    const Snarl* right_child_child = snarl_manager.children_of(right_child).at(0);
    
    if (right_child_child->end().node_id() < right_child_child->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(right_child_child);
    }
    
    REQUIRE(right_child_child->start().node_id() == 5);
    REQUIRE(right_child_child->end().node_id() == 6);
    
    // Define the chromosome by telomere snarls (first and last)
    auto chromosome = make_pair(top_snarl, top_snarl);
    
    // Make a genome state for this genome, with the only top snarl being the
    // telomere snarls for the main chromosome
    GenomeState state(snarl_manager, &graph, {chromosome});
    
    // Make sure to get all the net graphs
    const NetGraph* top_graph = state.get_net_graph(top_snarl);
    const NetGraph* left_graph = state.get_net_graph(left_child);
    const NetGraph* right_graph = state.get_net_graph(right_child);
    const NetGraph* right_child_graph = state.get_net_graph(right_child_child);
    
    SECTION("A haplotype can be added") {
        // Define a haplotype across the entire graph, in three levels
        InsertHaplotypeCommand insert;
        
        // We fill these with handles from the appropriate net graphs.
        // TODO: really we could just use the backing graph. Would that be better???
        
        // For the top snarl we go 1, 2 (child chain), and 8
        insert.insertions.emplace(top_snarl, vector<vector<pair<handle_t, size_t>>>{{
            {top_graph->get_handle(1, false), 0},
            {top_graph->get_handle_from_inward_backing_handle(graph.get_handle(2, false)), 0},
            {top_graph->get_handle(8, false), 0}
        }});
        
        // For the left child snarl we go 2, 3, 4
        insert.insertions.emplace(left_child, vector<vector<pair<handle_t, size_t>>>{{
            {left_graph->get_handle(2, false), 0},
            {left_graph->get_handle(3, false), 0},
            {left_graph->get_handle(4, false), 0}
        }});
        
        // For the right child snarl we go backward: 7 rev, 6 rev (child chain), 4 rev.
        insert.insertions.emplace(right_child, vector<vector<pair<handle_t, size_t>>>{{
            {right_graph->get_handle(7, true), 0},
            {right_graph->get_handle_from_inward_backing_handle(graph.get_handle(6, true)), 0},
            {right_graph->get_handle(4, true), 0}
        }});
        
        // For the right child's child, we go forward and take the insertion of 9
        insert.insertions.emplace(right_child_child, vector<vector<pair<handle_t, size_t>>>{{
            {right_child_graph->get_handle(5, false), 0},
            {right_child_graph->get_handle(9, false), 0},
            {right_child_graph->get_handle(6, false), 0}
        }});
        
        // Execute the command and get the undo command
        GenomeStateCommand* undo = state.execute(&insert);
        
        SECTION("The haplotype can be traced") {
            // We trace out all the handles in the backing graph
            vector<handle_t> traced;
            
            state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                // Put every handle in the vector
                traced.push_back(visit);
            });
         
            // We should visit nodes 1, 2, 3, 4, 5, 9, 6, 7, 8 in that order
            REQUIRE(traced.size() == 9);
            REQUIRE(traced[0] == graph.get_handle(1, false));
            REQUIRE(traced[1] == graph.get_handle(2, false));
            REQUIRE(traced[2] == graph.get_handle(3, false));
            REQUIRE(traced[3] == graph.get_handle(4, false));
            REQUIRE(traced[4] == graph.get_handle(5, false));
            REQUIRE(traced[5] == graph.get_handle(9, false));
            REQUIRE(traced[6] == graph.get_handle(6, false));
            REQUIRE(traced[7] == graph.get_handle(7, false));
            REQUIRE(traced[8] == graph.get_handle(8, false));
        }
        
        SECTION("Part of the haplotype can be swapped out") {
            ReplaceSnarlHaplotypeCommand replace;
            replace.snarl = right_child;
            replace.lane = state.count_haplotypes(right_child) - 1;
            replace.haplotype.push_back(graph.get_handle(7, true));
            replace.haplotype.push_back(graph.get_handle(4, true));
            
            GenomeStateCommand* unreplace = state.execute(&replace);
            
            SECTION("The haplotype can be traced") {
                // We trace out all the handles in the backing graph
                vector<handle_t> traced;
                
                state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                    // Put every handle in the vector
                    traced.push_back(visit);
                });
             
                // We should visit nodes 1, 2, 3, 4, 7, 8 in that order
                REQUIRE(traced.size() == 6);
                REQUIRE(traced[0] == graph.get_handle(1, false));
                REQUIRE(traced[1] == graph.get_handle(2, false));
                REQUIRE(traced[2] == graph.get_handle(3, false));
                REQUIRE(traced[3] == graph.get_handle(4, false));
                REQUIRE(traced[4] == graph.get_handle(7, false));
                REQUIRE(traced[5] == graph.get_handle(8, false));
            }
            
            SECTION("The swap can be reversed") {
                GenomeStateCommand* rereplace = state.execute(unreplace);
                
                SECTION("The haplotype can be traced") {
                    // We trace out all the handles in the backing graph
                    vector<handle_t> traced;
                    
                    state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                        // Put every handle in the vector
                        traced.push_back(visit);
                    });
                 
                    // We should visit nodes 1, 2, 3, 4, 5, 9, 6, 7, 8 in that order
                    REQUIRE(traced.size() == 9);
                    REQUIRE(traced[0] == graph.get_handle(1, false));
                    REQUIRE(traced[1] == graph.get_handle(2, false));
                    REQUIRE(traced[2] == graph.get_handle(3, false));
                    REQUIRE(traced[3] == graph.get_handle(4, false));
                    REQUIRE(traced[4] == graph.get_handle(5, false));
                    REQUIRE(traced[5] == graph.get_handle(9, false));
                    REQUIRE(traced[6] == graph.get_handle(6, false));
                    REQUIRE(traced[7] == graph.get_handle(7, false));
                    REQUIRE(traced[8] == graph.get_handle(8, false));
                }
                
                SECTION("The swap can be redone") {
                    delete state.execute(rereplace);
                    
                    SECTION("The haplotype can be traced") {
                        // We trace out all the handles in the backing graph
                        vector<handle_t> traced;
                        
                        state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                            // Put every handle in the vector
                            traced.push_back(visit);
                        });
                     
                        // We should visit nodes 1, 2, 3, 4, 7, 8 in that order
                        REQUIRE(traced.size() == 6);
                        REQUIRE(traced[0] == graph.get_handle(1, false));
                        REQUIRE(traced[1] == graph.get_handle(2, false));
                        REQUIRE(traced[2] == graph.get_handle(3, false));
                        REQUIRE(traced[3] == graph.get_handle(4, false));
                        REQUIRE(traced[4] == graph.get_handle(7, false));
                        REQUIRE(traced[5] == graph.get_handle(8, false));
                    }
                    
                }
                
                delete rereplace;
                
            }
            
            delete unreplace;
        }
        
        delete undo;
        
    }
    
    SECTION("A haplotype can be appended") {
        // Define a haplotype across the entire graph, all at once
        AppendHaplotypeCommand append;
        
        // We should visit nodes 1, 2, 3, 4, 5, 9, 6, 7, 8 in that order
        append.haplotype.push_back(graph.get_handle(1, false));
        append.haplotype.push_back(graph.get_handle(2, false));
        append.haplotype.push_back(graph.get_handle(3, false));
        append.haplotype.push_back(graph.get_handle(4, false));
        append.haplotype.push_back(graph.get_handle(5, false));
        append.haplotype.push_back(graph.get_handle(9, false));
        append.haplotype.push_back(graph.get_handle(6, false));
        append.haplotype.push_back(graph.get_handle(7, false));
        append.haplotype.push_back(graph.get_handle(8, false));
        
        // Execute the command and get the undo command
        GenomeStateCommand* undo = state.execute(&append);
        
        SECTION("The added haplotype is counted") {
            REQUIRE(state.count_haplotypes(chromosome) == 1);
        }
        
        SECTION("The haplotype can be traced") {
            
            // We trace out all the handles in the backing graph
            vector<handle_t> traced;
            
            state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                // Put every handle in the vector
                traced.push_back(visit);
            });
         
            // We should visit nodes 1, 2, 3, 4, 5, 9, 6, 7, 8 in that order
            REQUIRE(traced.size() == 9);
            REQUIRE(traced[0] == graph.get_handle(1, false));
            REQUIRE(traced[1] == graph.get_handle(2, false));
            REQUIRE(traced[2] == graph.get_handle(3, false));
            REQUIRE(traced[3] == graph.get_handle(4, false));
            REQUIRE(traced[4] == graph.get_handle(5, false));
            REQUIRE(traced[5] == graph.get_handle(9, false));
            REQUIRE(traced[6] == graph.get_handle(6, false));
            REQUIRE(traced[7] == graph.get_handle(7, false));
            REQUIRE(traced[8] == graph.get_handle(8, false));
        }
        
        delete undo;
        
    }
    
}

TEST_CASE("GenomeState can have haplotypes appended", "[genomestate]") {
    

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
    
    // Work out its snarls
    CactusSnarlFinder bubble_finder(graph);
    SnarlManager snarl_manager = bubble_finder.find_snarls();

    // Get the top snarl
    const Snarl* top_snarl = snarl_manager.top_level_snarls().at(0);
    
    if (top_snarl->end().node_id() < top_snarl->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(top_snarl);
    }
    
    // Make sure it's what we expect.
    REQUIRE(top_snarl->start().node_id() == 1);
    REQUIRE(top_snarl->end().node_id() == 8);
    
    // Get the middle snarl
    const Snarl* middle_snarl = snarl_manager.children_of(top_snarl).at(0);
    
    if (middle_snarl->end().node_id() < middle_snarl->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(middle_snarl);
    }
    
    // And the bottom snarl
    const Snarl* bottom_snarl = snarl_manager.children_of(middle_snarl).at(0);
    
    if (bottom_snarl->end().node_id() < bottom_snarl->start().node_id()) {
        // Put it a consistent way around
        snarl_manager.flip(bottom_snarl);
    }
    
    // Define the chromosome by telomere snarls (first and last)
    auto chromosome = make_pair(top_snarl, top_snarl);
    
    // Make a genome state for this genome, with the only top snarl being the
    // telomere snarls for the main chromosome
    GenomeState state(snarl_manager, &graph, {chromosome});
    
    SECTION("GenomeState starts empty") {
        REQUIRE(state.count_haplotypes(chromosome) == 0);
    }
    
    SECTION("A haplotype can be appended") {
        // Define a haplotype across the entire graph, all at once
        AppendHaplotypeCommand append;
        
        append.haplotype.push_back(graph.get_handle(1, false));
        append.haplotype.push_back(graph.get_handle(2, false));
        append.haplotype.push_back(graph.get_handle(3, false));
        append.haplotype.push_back(graph.get_handle(5, false));
        append.haplotype.push_back(graph.get_handle(7, false));
        append.haplotype.push_back(graph.get_handle(8, false));
        
        // Execute the command and get the undo command
        GenomeStateCommand* undo = state.execute(&append);
        
        SECTION("The added haplotype is counted") {
            REQUIRE(state.count_haplotypes(chromosome) == 1);
        }
        
        SECTION("The haplotype can be traced") {
            // We trace out all the handles in the backing graph
            vector<handle_t> traced;
            
            state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                // Put every handle in the vector
                traced.push_back(visit);
            });
         
            // We should visit nodes 1, 2, 3, 5, 7, 8 in that order
            REQUIRE(traced.size() == 6);
            REQUIRE(traced[0] == graph.get_handle(1, false));
            REQUIRE(traced[1] == graph.get_handle(2, false));
            REQUIRE(traced[2] == graph.get_handle(3, false));
            REQUIRE(traced[3] == graph.get_handle(5, false));
            REQUIRE(traced[4] == graph.get_handle(7, false));
            REQUIRE(traced[5] == graph.get_handle(8, false));
        }
        
        SECTION("Part of the haplotype can be replaced") {
            ReplaceSnarlHaplotypeCommand replace;
            replace.snarl = bottom_snarl;
            replace.lane = state.count_haplotypes(bottom_snarl) - 1;
            replace.haplotype.push_back(graph.get_handle(3, false));
            replace.haplotype.push_back(graph.get_handle(4, false));
            replace.haplotype.push_back(graph.get_handle(5, false));
            
            GenomeStateCommand* unreplace = state.execute(&replace);
            
            SECTION("The haplotype can be traced") {
                // We trace out all the handles in the backing graph
                vector<handle_t> traced;
                
                state.trace_haplotype(chromosome, 0, [&](const handle_t& visit) {
                    // Put every handle in the vector
                    traced.push_back(visit);
                });
             
                // We should visit nodes 1, 2, 3, 4, 5, 7, 8 in that order
                REQUIRE(traced.size() == 7);
                REQUIRE(traced[0] == graph.get_handle(1, false));
                REQUIRE(traced[1] == graph.get_handle(2, false));
                REQUIRE(traced[2] == graph.get_handle(3, false));
                REQUIRE(traced[3] == graph.get_handle(4, false));
                REQUIRE(traced[4] == graph.get_handle(5, false));
                REQUIRE(traced[5] == graph.get_handle(7, false));
                REQUIRE(traced[6] == graph.get_handle(8, false));
            }            
            
            delete unreplace;
        }
        
        SECTION("The haplotype can be removed again") {
            delete state.execute(undo);
            
            REQUIRE(state.count_haplotypes(chromosome) == 0);
        }
        
        delete undo;
    }
}
        
}
}

















