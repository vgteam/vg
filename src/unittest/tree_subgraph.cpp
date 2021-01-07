/// \file tree_subgraph.cpp
///  
/// Unit tests for the TreeSubgraph, which defines a subgraph with a tree topology.
///

#include <iostream>
#include <string>
#include "vg/io/json2pb.h"
#include "../tree_subgraph.hpp"
#include "catch.hpp"

#include <vg/vg.pb.h>
#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {
using namespace std;



TEST_CASE("TreeSubgraph can be created and traversed", "[subgraph]") {
    
    // Make up a base graph
    bdsg::HashGraph base;
    
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
    
    SECTION("works with no root trimming") {
    
        TreeSubgraph subgraph(&base, std::move(tree));
        tree.clear();
        
        // Make sure it works.
        
        // We need to get ahold of the root handle.
        vector<handle_t> order = subgraph.get_topological_order();
        const handle_t& root = order.at(0);
        
        REQUIRE(root == subgraph.get_root());
        
        SECTION("root node is correct") {
        
            // Make sure it is correct
            REQUIRE(subgraph.get_sequence(root) == "GAT");
            REQUIRE(subgraph.get_degree(root, false) == 1);
            REQUIRE(subgraph.get_degree(root, true) == 0);
            
            // Make sure it is numbered as we know the TreeSubgraph numbers things (i.e. by index in the tree)
            REQUIRE(subgraph.get_id(root) == 1);
            
        }
        
        SECTION("Size is correct") {
        
            // Make sure we have the right number of nodes
            REQUIRE(subgraph.get_node_count() == 6);
            
        }
        
        SECTION("Traversal from root is correct") {
        
            // Traverse the graph
            unordered_map<handle_t, size_t> seen;
            handlealgs::dijkstra(&subgraph, root, [&](const handle_t& reached, size_t distance) {
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
        
        SECTION("Root connectivity is correct") {
            vector<handle_t> neighbors;
            
            subgraph.follow_edges(root, true, [&](const handle_t& found) {
                neighbors.push_back(found);
            });
            
            // No left neighbors
            REQUIRE(neighbors.empty());
            
            neighbors.clear();
            
            subgraph.follow_edges(root, false, [&](const handle_t& found) {
                neighbors.push_back(found);
            });
            
            // One right neighbor
            REQUIRE(neighbors.size() == 1);
            // Should be at index 2 in the tree.
            REQUIRE(subgraph.get_id(neighbors[0]) == 2);
        }
        
        SECTION("Leaf connectivity is correct") {
            vector<handle_t> neighbors;
            
            handle_t leaf1 = subgraph.get_handle(5);
            
            subgraph.follow_edges(leaf1, true, [&](const handle_t& found) {
                neighbors.push_back(found);
            });
            
            // One left neighbor
            REQUIRE(neighbors.size() == 1);
            // Shouuld be at index 3 in the tree.
            REQUIRE(subgraph.get_id(neighbors[0]) == 3);
            
            neighbors.clear();
            
            subgraph.follow_edges(leaf1, false, [&](const handle_t& found) {
                neighbors.push_back(found);
            });
            
            // No right neighbors
            REQUIRE(neighbors.empty());
        }
    }
    
    SECTION("works with root trimming") {
        
        TreeSubgraph subgraph(&base, std::move(tree), 2);
        tree.clear();
        
        handle_t trimmed_root = subgraph.get_handle(1, false);
        
        REQUIRE(subgraph.get_sequence(trimmed_root) == "T");
        REQUIRE(subgraph.get_sequence(subgraph.flip(trimmed_root)) == "A");
        
        Path forward_path;
        json2pb(forward_path, R"(
            {"mapping": [
                {"position": {"node_id": 1}, "edit": [{"from_length": 1, "to_length": 1}]},
                {"position": {"node_id": 2}, "edit": [{"from_length": 2, "to_length": 2}]},
                {"position": {"node_id": 4}, "edit": [{"from_length": 1, "to_length": 1}]},
                {"position": {"node_id": 6}, "edit": [{"from_length": 1, "to_length": 1}]}
            ]}
        )");
        
        Path forward_translated = subgraph.translate_down(forward_path);
        
        // The last mapping should move in ID
        REQUIRE(forward_translated.mapping_size() == 4);
        REQUIRE(forward_translated.mapping(0).position().node_id() == 1);
        REQUIRE(forward_translated.mapping(1).position().node_id() == 2);
        REQUIRE(forward_translated.mapping(2).position().node_id() == 4);
        REQUIRE(forward_translated.mapping(3).position().node_id() == 5);
        // And the first should get an offset
        REQUIRE(forward_translated.mapping(0).position().offset() == 2);
        
        
        Path reverse_path;
        json2pb(reverse_path, R"(
            {"mapping": [
                {"position": {"node_id": 6, "is_reverse": true}, "edit": [{"from_length": 1, "to_length": 1}]},
                {"position": {"node_id": 4, "is_reverse": true}, "edit": [{"from_length": 1, "to_length": 1}]},
                {"position": {"node_id": 2, "is_reverse": true}, "edit": [{"from_length": 2, "to_length": 2}]},
                {"position": {"node_id": 1, "is_reverse": true}, "edit": [{"from_length": 1, "to_length": 1}]}
            ]}
        )");
        
        Path reverse_translated = subgraph.translate_down(reverse_path);
        
        // The first mapping should move in ID
        REQUIRE(reverse_translated.mapping_size() == 4);
        REQUIRE(reverse_translated.mapping(0).position().node_id() == 5);
        REQUIRE(reverse_translated.mapping(1).position().node_id() == 4);
        REQUIRE(reverse_translated.mapping(2).position().node_id() == 2);
        REQUIRE(reverse_translated.mapping(3).position().node_id() == 1);
        // And the last should not get an offset because it is in the same place as before on the reverse strand.
        REQUIRE(reverse_translated.mapping(3).position().offset() == 0);
        
    }
    
}

}
}
        
