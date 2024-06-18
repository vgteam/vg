/// \file packed_graph.cpp
///  
/// Unit tests for the PackedGraph class.
///

#include <iostream>
#include <sstream>

#include "vg/io/json2pb.h"
#include "../split_strand_graph.hpp"
#include "../subpath_overlay.hpp"
#include "../utility.hpp"
#include "../handle.hpp"
#include "random_graph.hpp"

#include "../vg.hpp"

#include "bdsg/hash_graph.hpp"

#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;
using namespace bdsg;

TEST_CASE("StrandSplitGraph overlay produces graphs that are single stranded", "[overlay][handle]") {
    
    int num_tests = 100;
    int seq_size = 500;
    int var_len = 10;
    int var_count = 30;
    
    for (int i = 0; i < num_tests; ++i) {
        
        bdsg::HashGraph graph;
        random_graph(seq_size, var_len, var_count, &graph);

        StrandSplitGraph split(&graph);
        REQUIRE(handlealgs::is_single_stranded(&split));
        
        REQUIRE(split.get_node_count() == 2 * graph.get_node_count());
        
        unordered_map<handle_t, int> node_count;
        
        int split_edge_trav_count = 0;
        int graph_edge_trav_count = 0;
        
        split.for_each_handle([&](const handle_t& h) {
            REQUIRE(split.get_sequence(h) == graph.get_sequence(split.get_underlying_handle(h)));
            
            split.follow_edges(h, false, [&](const handle_t& n) {
                REQUIRE(graph.has_edge(split.get_underlying_handle(h),
                                       split.get_underlying_handle(n)));
                split_edge_trav_count++;
            });
            split.follow_edges(h, true, [&](const handle_t& p) {
                REQUIRE(graph.has_edge(split.get_underlying_handle(p),
                                       split.get_underlying_handle(h)));
                split_edge_trav_count++;
            });
            
            node_count[split.get_underlying_handle(h)]++;
        });
        
        graph.for_each_handle([&](const handle_t& h) {
            REQUIRE(node_count[h] == 1);
            REQUIRE(node_count[graph.flip(h)] == 1);
            graph.follow_edges(h, false, [&](const handle_t& n) {
                graph_edge_trav_count++;
            });
            graph.follow_edges(h, true, [&](const handle_t& p) {
                graph_edge_trav_count++;
            });
        });
        
        REQUIRE(split_edge_trav_count == 2 * graph_edge_trav_count);
    }
}

TEST_CASE("SubpathOverlay works as expected", "[overlay][handle]") {
    
    
    HashGraph graph;
    auto h1 = graph.create_handle("A");
    auto h2 = graph.create_handle("C");
    auto h3 = graph.create_handle("G");
    auto h4 = graph.create_handle("T");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
    
    auto p = graph.create_path_handle("p");
    
    vector<step_handle_t> steps;
    steps.push_back(graph.append_step(p, h1));
    steps.push_back(graph.append_step(p, h2));
    steps.push_back(graph.append_step(p, h4));
    
    for (size_t i = 0; i < steps.size(); ++i) {
        for (size_t j = i; j < steps.size(); ++j) {
            
            SubpathOverlay subpath(&graph, steps[i], graph.get_next_step(steps[j]));
            
            REQUIRE(handlealgs::is_single_stranded(&subpath));
            
            REQUIRE(subpath.get_node_count() == j - i + 1);
            size_t manual_node_count = 0;
            subpath.for_each_handle([&](const handle_t& h) {
                ++manual_node_count;
            });
            REQUIRE(manual_node_count == subpath.get_node_count());
            
            auto nodes = handlealgs::lazier_topological_order(&subpath);
            
            for (size_t k = 0; k < nodes.size(); ++k) {
                auto s = graph.get_handle_of_step(steps[i + k]);
                REQUIRE(subpath.get_underlying_handle(nodes[k]) == s);
                REQUIRE(subpath.get_sequence(nodes[k]) == graph.get_sequence(s));
                
                REQUIRE(subpath.get_degree(nodes[k], true) == (k == 0 ? 0 : 1));
                REQUIRE(subpath.get_degree(nodes[k], false) == (k + 1 == nodes.size() ? 0 : 1));
            }
        }
    }
}


}
}
        
