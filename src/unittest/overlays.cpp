/// \file packed_graph.cpp
///  
/// Unit tests for the PackedGraph class.
///

#include <iostream>
#include <sstream>

#include "vg/io/json2pb.h"
#include "../split_strand_graph.hpp"
#include "../utility.hpp"
#include "random_graph.hpp"

#include "../vg.hpp"

#include "bdsg/hash_graph.hpp"

#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;

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
}
}
        
