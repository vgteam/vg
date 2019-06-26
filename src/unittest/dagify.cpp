/// \file dagify.cpp
///  
/// unit tests for the handle based dagify algorithm
///

#include <iostream>
#include "json2pb.h"
#include "../vg.hpp"
#include "../split_strand_graph.hpp"
#include "../algorithms/dagify.hpp"
#include "../algorithms/is_acyclic.hpp"
#include "../algorithms/split_strands.hpp"
#include "random_graph.hpp"
#include "catch.hpp"

#include "sglib/hash_graph.hpp"

namespace vg {
namespace unittest {

    TEST_CASE("Dagify algorithm produces expected results for small test cases", "[algorithms][dagify]") {
        
        SECTION("Dagify can unroll a small loop between two nodes") {
            
            VG graph;
            
            handle_t n1 = graph.create_handle("AA");
            handle_t n2 = graph.create_handle("AC");
            handle_t n3 = graph.create_handle("AG");
            handle_t n4 = graph.create_handle("CC");
            
            graph.create_edge(n1, n2);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n2);
            
            VG dagified;
            
            size_t preserved_length = 1;
            
            unordered_map<id_t, id_t> trans = algorithms::dagify(&graph, &dagified, preserved_length);
            
            REQUIRE(algorithms::is_acyclic(&dagified));
            REQUIRE(dagified.get_node_count() == 6);
            
            handle_t d1, d2, d3, d4, d5, d6;
            bool f1 = false, f2 = false, f3 = false, f4 = false, f5 = false, f6 = false;
            
            dagified.for_each_handle([&](const handle_t& h) {
                if (dagified.get_sequence(h) == "AA") {
                    REQUIRE(!f1);
                    d1 = h;
                    f1 = true;
                }
                else if (dagified.get_sequence(h) == "AC") {
                    REQUIRE(!(f2 && f3));
                    if (f2) {
                        d3 = h;
                        f3 = true;
                    }
                    else {
                        d2 = h;
                        f2 = true;
                    }
                }
                else if (dagified.get_sequence(h) == "AG") {
                    REQUIRE(!(f4 && f5));
                    if (f4) {
                        d5 = h;
                        f5 = true;
                    }
                    else {
                        d4 = h;
                        f4 = true;
                    }
                }
                else if (dagified.get_sequence(h) == "CC") {
                    REQUIRE(!f6);
                    d6 = h;
                    f6 = true;
                }
                else {
                    REQUIRE(false);
                }
            });
            
            REQUIRE(f1);
            REQUIRE(f2);
            REQUIRE(f3);
            REQUIRE(f4);
            REQUIRE(f5);
            REQUIRE(f6);
            
            // all walks of length 2 should have been preserved in the new DAG, I enumerate them here
            vector<vector<handle_t>> walks;
            walks.emplace_back(); walks.back().push_back(n1); walks.back().push_back(n2);
            walks.emplace_back(); walks.back().push_back(n2); walks.back().push_back(n3);
            walks.emplace_back(); walks.back().push_back(n2); walks.back().push_back(n4);
            walks.emplace_back(); walks.back().push_back(n3); walks.back().push_back(n2);
            
            for (auto& walk : walks) {
                bool found = false;
                dagified.for_each_handle([&](const handle_t& h) {
                    // exhaustively enumerate all walks of the target length
                    vector<vector<handle_t>> local_walks;
                    local_walks.emplace_back(vector<handle_t>{h});
                    
                    while (local_walks.front().size() < walk.size() && !local_walks.empty()) {
                        vector<vector<handle_t>> next_local_walks;
                        
                        for (auto& prev_walk : local_walks) {
                            dagified.follow_edges(prev_walk.back(), false, [&](const handle_t& next) {
                                next_local_walks.push_back(prev_walk);
                                next_local_walks.back().push_back(next);
                            });
                        }
                        
                        local_walks = next_local_walks;
                    }
                    
                    for (auto& local_walk : local_walks) {
                        bool all_match = true;
                        for (size_t i = 0; i < walk.size(); i++) {
                            all_match = all_match && (graph.get_id(walk[i]) == trans[dagified.get_id(local_walk[i])] &&
                                                      graph.get_is_reverse(walk[i]) == dagified.get_is_reverse(local_walk[i]));
                        }
                        
                        if (all_match) {
                            found = true;
                            break;
                        }
                    }
                    
                    return !found;
                });
                
                REQUIRE(found);
            }
        }
        
        SECTION("Dagify can unroll a small loop between two nodes for a longer distance") {

            VG graph;

            handle_t n1 = graph.create_handle("AA");
            handle_t n2 = graph.create_handle("AC");
            handle_t n3 = graph.create_handle("AG");
            handle_t n4 = graph.create_handle("CC");

            graph.create_edge(n1, n2);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n2);

            VG dagified;

            size_t preserved_length = 5;

            unordered_map<id_t, id_t> trans = algorithms::dagify(&graph, &dagified, preserved_length);
            
            REQUIRE(algorithms::is_acyclic(&dagified));
            REQUIRE(dagified.get_node_count() == 8);

            handle_t d1, d2, d3, d4, d5, d6, d7, d8;
            bool f1 = false, f2 = false, f3 = false, f4 = false, f5 = false, f6 = false, f7 = false, f8 = false;

            dagified.for_each_handle([&](const handle_t& h) {
                if (dagified.get_sequence(h) == "AA") {
                    REQUIRE(!f1);
                    d1 = h;
                    f1 = true;
                }
                else if (dagified.get_sequence(h) == "AC") {
                    REQUIRE(!(f2 && f3 && f4));
                    if (f2 && f3) {
                        d4 = h;
                        f4 = true;
                    }
                    else if (f2) {
                        d3 = h;
                        f3 = true;
                    }
                    else {
                        d2 = h;
                        f2 = true;
                    }
                }
                else if (dagified.get_sequence(h) == "AG") {
                    REQUIRE(!(f5 && f6 && f7));
                    if (f5 && f6) {
                        d7 = h;
                        f7 = true;
                    }
                    else if (f5) {
                        d6 = h;
                        f6 = true;
                    }
                    else {
                        d5 = h;
                        f5 = true;
                    }
                }
                else if (dagified.get_sequence(h) == "CC") {
                    REQUIRE(!f8);
                    d8 = h;
                    f8 = true;
                }
                else {
                    REQUIRE(false);
                }
            });

            REQUIRE(f1);
            REQUIRE(f2);
            REQUIRE(f3);
            REQUIRE(f4);
            REQUIRE(f5);
            REQUIRE(f6);
            REQUIRE(f7);
            REQUIRE(f8);

            // all walks of length 3 should have been preserved in the new DAG, I enumerate them here
            vector<vector<handle_t>> walks;
            walks.emplace_back(); walks.back().push_back(n1); walks.back().push_back(n2); walks.back().push_back(n3);
            walks.emplace_back(); walks.back().push_back(n1); walks.back().push_back(n2); walks.back().push_back(n4);
            walks.emplace_back(); walks.back().push_back(n2); walks.back().push_back(n3); walks.back().push_back(n2);
            walks.emplace_back(); walks.back().push_back(n3); walks.back().push_back(n2); walks.back().push_back(n3);
            walks.emplace_back(); walks.back().push_back(n3); walks.back().push_back(n2); walks.back().push_back(n4);

            for (auto& walk : walks) {
                bool found = false;
                dagified.for_each_handle([&](const handle_t& h) {
                    // exhaustively enumerate all walks of the target length
                    vector<vector<handle_t>> local_walks;
                    local_walks.emplace_back(vector<handle_t>{h});
                    
                    while (local_walks.front().size() < walk.size() && !local_walks.empty()) {
                        vector<vector<handle_t>> next_local_walks;
                        
                        for (auto& prev_walk : local_walks) {
                            dagified.follow_edges(prev_walk.back(), false, [&](const handle_t& next) {
                                next_local_walks.push_back(prev_walk);
                                next_local_walks.back().push_back(next);
                            });
                        }
                        
                        local_walks = next_local_walks;
                    }
                    
                    for (auto& local_walk : local_walks) {
                        bool all_match = true;
                        for (size_t i = 0; i < walk.size(); i++) {
                            all_match = all_match && (graph.get_id(walk[i]) == trans[dagified.get_id(local_walk[i])] &&
                                                      graph.get_is_reverse(walk[i]) == dagified.get_is_reverse(local_walk[i]));
                        }
                        
                        if (all_match) {
                            found = true;
                            break;
                        }
                    }
                    
                    return !found;
                });
                
                REQUIRE(found);
            }
        }
        
        SECTION("Dagify can unroll a small loop between two nodes that has reversing edges") {
            
            VG graph;
            
            handle_t n1 = graph.create_handle("AA");
            handle_t n2 = graph.create_handle("AC");
            handle_t n3 = graph.create_handle("AG");
            handle_t n4 = graph.create_handle("CC");
            
            graph.create_edge(n1, n2);
            graph.create_edge(n2, graph.flip(n3));
            graph.create_edge(n2, graph.flip(n4));
            graph.create_edge(graph.flip(n3), n2);
            
            VG dagified;
            
            size_t preserved_length = 1;
            
            unordered_map<id_t, id_t> trans = algorithms::dagify(&graph, &dagified, preserved_length);
            
            REQUIRE(algorithms::is_acyclic(&dagified));
            REQUIRE(dagified.get_node_count() == 6);
            
            handle_t d1, d2, d3, d4, d5, d6;
            bool f1 = false, f2 = false, f3 = false, f4 = false, f5 = false, f6 = false;
            
            dagified.for_each_handle([&](const handle_t& h) {
                if (dagified.get_sequence(h) == "AA") {
                    REQUIRE(!f1);
                    d1 = h;
                    f1 = true;
                }
                else if (dagified.get_sequence(h) == "AC") {
                    REQUIRE(!(f2 && f3));
                    if (f2) {
                        d3 = h;
                        f3 = true;
                    }
                    else {
                        d2 = h;
                        f2 = true;
                    }
                }
                else if (dagified.get_sequence(h) == "AG") {
                    REQUIRE(!(f4 && f5));
                    if (f4) {
                        d5 = h;
                        f5 = true;
                    }
                    else {
                        d4 = h;
                        f4 = true;
                    }
                }
                else if (dagified.get_sequence(h) == "CC") {
                    REQUIRE(!f6);
                    d6 = h;
                    f6 = true;
                }
                else {
                    REQUIRE(false);
                }
            });
            
            REQUIRE(f1);
            REQUIRE(f2);
            REQUIRE(f3);
            REQUIRE(f4);
            REQUIRE(f5);
            REQUIRE(f6);
            
            // all walks of length 2 should have been preserved in the new DAG, I enumerate them here
            vector<vector<handle_t>> walks;
            walks.emplace_back(); walks.back().push_back(n1); walks.back().push_back(n2);
            walks.emplace_back(); walks.back().push_back(n2); walks.back().push_back(graph.flip(n3));
            walks.emplace_back(); walks.back().push_back(n2); walks.back().push_back(graph.flip(n4));
            walks.emplace_back(); walks.back().push_back(graph.flip(n3)); walks.back().push_back(n2);
            
            for (auto& walk : walks) {
                
                bool found = false;
                dagified.for_each_handle([&](const handle_t& h) {
                    
                    // go on the reverse for the nodes we reversed
                    handle_t start = h;
                    if (trans[dagified.get_id(start)] == graph.get_id(n3) ||
                        trans[dagified.get_id(start)] == graph.get_id(n4)) {
                        start = dagified.flip(start);
                    }
                    
                    // exhaustively enumerate all walks of the target length
                    vector<vector<handle_t>> local_walks;
                    local_walks.emplace_back(vector<handle_t>{start});
                    
                    while (local_walks.front().size() < walk.size() && !local_walks.empty()) {
                        vector<vector<handle_t>> next_local_walks;
                        
                        for (auto& prev_walk : local_walks) {
                            dagified.follow_edges(prev_walk.back(), false, [&](const handle_t& next) {
                                next_local_walks.push_back(prev_walk);
                                next_local_walks.back().push_back(next);
                            });
                        }
                        
                        local_walks = next_local_walks;
                    }
                    
                    for (auto& local_walk : local_walks) {
                        bool all_match = true;
                        for (size_t i = 0; i < walk.size(); i++) {
                            all_match = all_match && (graph.get_id(walk[i]) == trans[dagified.get_id(local_walk[i])] &&
                                                      graph.get_is_reverse(walk[i]) == dagified.get_is_reverse(local_walk[i]));
                        }
                        
                        if (all_match) {
                            found = true;
                            break;
                        }
                    }
                    
                    return !found;
                });
                
                REQUIRE(found);
            }
        }
    }
    
    TEST_CASE("Dagify algorithm produces expected results random test cases using both methods of making single stranded graphs", "[algorithms][dagify]") {
        
        size_t preserved_length = 15;
        size_t seq_size = 50;
        size_t variant_len = 5;
        size_t variant_count = 8;
        
        size_t num_trials = 1000;
        for (size_t i = 0; i < num_trials; ++i) {
            sglib::HashGraph graph;
            random_graph(seq_size, variant_len, variant_count, &graph);
            
            sglib::HashGraph direct_split;
            algorithms::split_strands(&graph, &direct_split);
            sglib::HashGraph direct_dagified;
            algorithms::dagify(&direct_split, &direct_dagified, preserved_length);

            StrandSplitGraph split(&graph);
            
            sglib::HashGraph dagified;
            algorithms::dagify(&split, &dagified, preserved_length);
            
            REQUIRE(algorithms::is_acyclic(&direct_dagified));
            REQUIRE(algorithms::is_acyclic(&dagified));
        }
    }
}
}
