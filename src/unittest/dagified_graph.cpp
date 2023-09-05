/// \file dagified_graph.cpp
///  
/// unit tests for the overlay based dagify algorithm
///

#include <iostream>
#include "vg/io/json2pb.h"
#include "../dagified_graph.hpp"
#include "../split_strand_graph.hpp"
#include "random_graph.hpp"
#include "catch.hpp"

#include "bdsg/hash_graph.hpp"

namespace vg {
namespace unittest {

    TEST_CASE("DagifiedGraph produces expected results for small test cases", "[dagify][overlay]") {
        
        SECTION("DagifiedGraph can unroll a small loop between two nodes") {
            
            bdsg::HashGraph graph;
            
            handle_t n1 = graph.create_handle("AA");
            handle_t n2 = graph.create_handle("AC");
            handle_t n3 = graph.create_handle("AG");
            handle_t n4 = graph.create_handle("CC");
            
            graph.create_edge(n1, n2);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n2);
            
            size_t preserved_length = 1;
            
            DagifiedGraph dagified(&graph, preserved_length);
            
            REQUIRE(handlealgs::is_acyclic(&dagified));
            REQUIRE(dagified.get_node_count() == 6);
            
            id_t min_id = numeric_limits<id_t>::max(), max_id = numeric_limits<id_t>::min();
            size_t count = 0;
            dagified.for_each_handle([&](const handle_t& h) {
                min_id = min(min_id, dagified.get_id(h));
                max_id = max(max_id, dagified.get_id(h));
                ++count;
                REQUIRE(graph.get_sequence(dagified.get_underlying_handle(h)) == dagified.get_sequence(h));
            });
            
            REQUIRE(min_id == dagified.min_node_id());
            REQUIRE(max_id == dagified.max_node_id());
            REQUIRE(count == dagified.get_node_count());
            
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
                            all_match = all_match && walk[i] == dagified.get_underlying_handle(local_walk[i]);
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
        
        SECTION("DagifiedGraph can unroll a small loop between two nodes for a longer distance") {

            bdsg::HashGraph graph;

            handle_t n1 = graph.create_handle("AA");
            handle_t n2 = graph.create_handle("AC");
            handle_t n3 = graph.create_handle("AG");
            handle_t n4 = graph.create_handle("CC");

            graph.create_edge(n1, n2);
            graph.create_edge(n2, n3);
            graph.create_edge(n2, n4);
            graph.create_edge(n3, n2);

            size_t preserved_length = 5;
            
            DagifiedGraph dagified(&graph, preserved_length);
            
            REQUIRE(handlealgs::is_acyclic(&dagified));
            REQUIRE(dagified.get_node_count() == 8);
            
            id_t min_id = numeric_limits<id_t>::max(), max_id = numeric_limits<id_t>::min();
            size_t count = 0;
            dagified.for_each_handle([&](const handle_t& h) {
                min_id = min(min_id, dagified.get_id(h));
                max_id = max(max_id, dagified.get_id(h));
                ++count;
                REQUIRE(graph.get_sequence(dagified.get_underlying_handle(h)) == dagified.get_sequence(h));
            });
            
            REQUIRE(min_id == dagified.min_node_id());
            REQUIRE(max_id == dagified.max_node_id());
            REQUIRE(count == dagified.get_node_count());

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
                            all_match = all_match && walk[i] == dagified.get_underlying_handle(local_walk[i]);
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
            
            bdsg::HashGraph graph;
            
            handle_t n1 = graph.create_handle("AA");
            handle_t n2 = graph.create_handle("AC");
            handle_t n3 = graph.create_handle("AG");
            handle_t n4 = graph.create_handle("CC");
            
            graph.create_edge(n1, n2);
            graph.create_edge(n2, graph.flip(n3));
            graph.create_edge(n2, graph.flip(n4));
            graph.create_edge(graph.flip(n3), n2);
                        
            size_t preserved_length = 1;
            
            DagifiedGraph dagified(&graph, preserved_length);
            
            REQUIRE(handlealgs::is_acyclic(&dagified));
            REQUIRE(dagified.get_node_count() == 6);
            
            id_t min_id = numeric_limits<id_t>::max(), max_id = numeric_limits<id_t>::min();
            size_t count = 0;
            dagified.for_each_handle([&](const handle_t& h) {
                min_id = min(min_id, dagified.get_id(h));
                max_id = max(max_id, dagified.get_id(h));
                ++count;
                REQUIRE(graph.get_sequence(dagified.get_underlying_handle(h)) == dagified.get_sequence(h));
            });
            
            REQUIRE(min_id == dagified.min_node_id());
            REQUIRE(max_id == dagified.max_node_id());
            REQUIRE(count == dagified.get_node_count());
            
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
                    if (graph.get_id(dagified.get_underlying_handle(start)) == graph.get_id(n3) ||
                        graph.get_id(dagified.get_underlying_handle(start)) == graph.get_id(n4)) {
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
                            all_match = all_match && walk[i] == dagified.get_underlying_handle(local_walk[i]);
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
    
    TEST_CASE("DagifiedGraph produces expected results random test cases", "[overlay][dagify]") {
        
        size_t preserved_length = 15;
        size_t seq_size = 50;
        size_t variant_len = 5;
        size_t variant_count = 8;
        
        size_t num_trials = 1000;
        for (size_t i = 0; i < num_trials; ++i) {
 
            bdsg::HashGraph graph;
            random_graph(seq_size, variant_len, variant_count, &graph);

            StrandSplitGraph split(&graph);
            
            DagifiedGraph dagified(&split, preserved_length);
            
            id_t min_id = numeric_limits<id_t>::max(), max_id = numeric_limits<id_t>::min();
            size_t count = 0;
            dagified.for_each_handle([&](const handle_t& h) {
                min_id = min(min_id, dagified.get_id(h));
                max_id = max(max_id, dagified.get_id(h));
                ++count;
                REQUIRE(split.get_sequence(dagified.get_underlying_handle(h)) == dagified.get_sequence(h));
            });
            
            REQUIRE(min_id == dagified.min_node_id());
            REQUIRE(max_id == dagified.max_node_id());
            REQUIRE(count == dagified.get_node_count());
            
            REQUIRE(handlealgs::is_acyclic(&dagified));
        }
    }
}
}
