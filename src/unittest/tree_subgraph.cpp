/// \file tree_subgraph.cpp
///  
/// Unit tests for the TreeSubgraph, which defines a subgraph with a tree topology.
///

#include <iostream>
#include <string>
#include "../json2pb.h"
#include "../tree_subgraph.hpp"
#include "../algorithms/dijkstra.hpp"
#include "catch.hpp"

#include <vg/vg.pb.h>
#include <sglib/hash_graph.hpp>

namespace vg {
namespace unittest {
using namespace std;



TEST_CASE("TreeSubgraph can be created and traversed", "[subgraph]") {
    
    // Make up a base graph
    sglib::HashGraph base;
    
    handle_t start = base.create_handle("GAT");
    handle_t middle = base.create_handle("TA");
    handle_t snp1 = base.create_handle("C");
    handle_t snp2 = base.create_handle("T");
    handle_t end = base.create_handle("A");
    
    base.create_edge(start, middle);
    base.create_edge(middle, snp1);
    base.create_edge(middle, snp2);
    base.create_edge(snp1, end);
    base.create_edge(snp2, end);
    
    // Make a tree of (prev tree item, base graph handle) pairs
    vector<pair<int64_t, handle_t>> tree;
    // Fill in the whole graph
    tree.emplace_back(-1, start);
    tree.emplace_back(0, middle);
    tree.emplace_back(1, snp1);
    tree.emplace_back(1, snp2);
    tree.emplace_back(2, end);
    tree.emplace_back(3, end);
    
    TreeSubgraph subgraph(&base, std::move(tree));
    tree.clear();
    
    // Make sure it works.
    
    // We need to get ahold of the root handle.
    handle_t root = subgraph.get_root();
    
    // Make sure it is correct
    REQUIRE(subgraph.get_sequence(root) == "GAT");
    REQUIRE(subgraph.get_degree(root, false) == 1);
    REQUIRE(subgraph.get_degree(root, true) == 0);
    
    // Make sure we have the right number of nodes
    REQUIRE(subgraph.get_node_count() == 6);
    
    // Traverse the graph
    unordered_map<handle_t, size_t> seen;
    algorithms::dijkstra(&subgraph, root, [&](const handle_t& reached, size_t distance) {
        seen[reached] = distance;
        return true;
    });
    
    // Make sure we see all of it
    REQUIRE(seen.size() == 6);
    
    unordered_map<handle_t, size_t> backing_counts;
    
    for (auto& kv : seen) {
        handle_t backing = subgraph.get_underlying_handle(kv.first);
        
        backing_counts[backing]++;
        
        if (backing == start) {
            // Distance 0 at start
            REQUIRE(kv.second == 0);
        } else if (backing == middle) {
            // Distance 0 at what we reach next also
            REQUIRE(kv.second == 0);
        } else if (backing == snp1 || backing == snp2) {
            // Distance 2 at the SNP
            REQUIRE(kv.second == 2);
        } else if (backing == end) {
            // Distance 3 after it
            REQUIRE(kv.second == 3);
        }
    }
    
    // Make sure we have the right stats on the backing handles
    REQUIRE(backing_counts.size() == 5);
    REQUIRE(backing_counts[start] == 1);
    REQUIRE(backing_counts[middle] == 1);
    REQUIRE(backing_counts[snp1] == 1);
    REQUIRE(backing_counts[snp2] == 1);
    REQUIRE(backing_counts[end] == 2);
    
}

}
}
        
