/// \file three_edge_connected_components.cpp
///  
/// Unit tests for the three edge connected components algorithms
///

#include "catch.hpp"
#include "random_graph.hpp"

#include "../algorithms/three_edge_connected_components.hpp"

#include <structures/union_find.hpp>

#include <iostream>
#include <vector>
#include <map>



namespace vg {
namespace unittest {
using namespace std;

// We use this global adjacencies vector and functions to look at it as a "current" graph
static vector<vector<size_t>> adjacencies;
   
static void for_each_connected_node(size_t node, const function<void(size_t)>& iteratee) {
#ifdef debug
    cerr << "Asked for edges of node " << node << endl;
#endif
    for (auto& other : adjacencies.at(node)) {
        iteratee(other);
    }
}

// Represent the results as a union-find for checking, in this global union-find for the current graph.
static structures::UnionFind components(adjacencies.size(), true);

static void component_callback(const function<void(const function<void(size_t)>&)>& for_each_member) {
#ifdef debug
    cerr << "Got component" << endl;
#endif
    size_t first = 0;
    bool is_first = true;
    for_each_member([&](size_t member) {
#ifdef debug
        cerr << "Component contained " << member << endl;
#endif
        if (is_first) {
            // Find the first member of each component
            first = member;
            is_first = false;
        } else {
            // And union everything into it
            components.union_groups(first, member);
        }
    });
}

/// Compute 3-edge-connected-components in a slow but almost certainly right way
static structures::UnionFind brute_force_3ecc(size_t count, const function<void(size_t, const function<void(size_t)>&)>& get_connected) {
    // Make a big edge list
    // Only store edges one way
    vector<pair<size_t, size_t>> edges;
    for (size_t i = 0; i < count; i++) {
        get_connected(i, [&](size_t connected) {
            if (i < connected) {
                edges.emplace_back(i, connected);
            }
        });
    }
    
    // We need to count how many times each node is still connected to each other node.
    // Store by IDs in order.
    // We could use an unordered map but then we'd have to go find the pair hash...
    map<pair<size_t, size_t>, size_t> times_connected;
    
    // And we count how many subgraphs we look at
    size_t subgraph_count = 0;
    
    for (size_t r1 = 0; r1 < edges.size(); r1++) {
        // For each first edge we could remove
        for (size_t r2 = r1 + 1; r2 < edges.size(); r2++) {
            // For each second edge we could remove
            
            // Do connected components
            structures::UnionFind cc(count, true);
            for (size_t i = 0; i < edges.size(); i++) {
                if (i != r1 && i != r2) {
                    // Only count non-removed edges
                    cc.union_groups(edges[i].first, edges[i].second);
                }
            }
            
            for (auto& group : cc.all_groups()) {
                for (size_t i = 0; i < group.size(); i++) {
                    for (size_t j = i + 1; j < group.size(); j++) {
                        // For each pair of nodes in the same connected component this time, count that pair
                        times_connected[make_pair(min(group[i], group[j]), max(group[i], group[j]))]++;
                    }
                }
            }
            
            // Remember we did a subgraph
            subgraph_count++;
        }
    }
    
    // Now find pairs of nodes that were connected every time no matter what 2
    // edges we removed, and union them to get the 3eccs.

    structures::UnionFind to_return(count, true);
    for (auto& pair_and_count : times_connected) {
        if (pair_and_count.second == subgraph_count) {
            // This pair was connected every time
            to_return.union_groups(pair_and_count.first.first, pair_and_count.first.second);
        }
    }
    
    return to_return;
}

/// Determine if two union-finds are equal.
/// Can't be const because we need all_groups()
static bool uf_equal(structures::UnionFind& a, structures::UnionFind& b) {
    for (auto& group : a.all_groups()) {
        for (size_t i = 0; i < group.size(); i++) {
            if (b.group_size(group[i]) != group.size()) {
                return false;
            }
            if (b.find_group(group[i]) != b.find_group(group[0])) {
                return false;
            }
        }
    }
    return true;
}

TEST_CASE("3 edge connected components algorithms handle basic cases", "[3ecc][algorithms]") {
    
    SECTION("A cycle with an overlapping cycle has the overlapping cycle merged out") {
        adjacencies = {{1, 2}, {0, 2, 2}, {0, 1, 1}};
        components = structures::UnionFind(adjacencies.size(), true);
        
        SECTION("Works with Cactus") {
            algorithms::three_edge_connected_components_dense_cactus(adjacencies.size(), for_each_connected_node, component_callback);
            
            REQUIRE(components.all_groups().size() == 2);
            REQUIRE(components.group_size(0) == 1);
            REQUIRE(components.group_size(1) == 2);
            REQUIRE(components.group_size(2) == 2);
        }
        
        SECTION("Works with Tsin 2014") {
            algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
            REQUIRE(components.all_groups().size() == 2);
            REQUIRE(components.group_size(0) == 1);
            REQUIRE(components.group_size(1) == 2);
            REQUIRE(components.group_size(2) == 2);
        }
    }
    
    SECTION("An empty graph") {
        adjacencies = {};
        components = structures::UnionFind(adjacencies.size(), true);
    
        SECTION("Works with Cactus") {
            algorithms::three_edge_connected_components_dense_cactus(adjacencies.size(), for_each_connected_node, component_callback);
            REQUIRE(components.size() == 0);
        }
        SECTION("Works with Tsin 2014") {
            algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            REQUIRE(components.size() == 0);
        }
    }
    
    SECTION("A 4-node connected graph is one component") {
        adjacencies = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
        components = structures::UnionFind(adjacencies.size(), true);
        
        SECTION("Works with Cactus") {
            algorithms::three_edge_connected_components_dense_cactus(adjacencies.size(), for_each_connected_node, component_callback);
            
            REQUIRE(components.all_groups().size() == 1);
            REQUIRE(components.group_size(0) == 4);
            REQUIRE(components.group_size(1) == 4);
            REQUIRE(components.group_size(2) == 4);
            REQUIRE(components.group_size(3) == 4);
        }
        
        SECTION("Works with Tsin 2014") {
            algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
            REQUIRE(components.all_groups().size() == 1);
            REQUIRE(components.group_size(0) == 4);
            REQUIRE(components.group_size(1) == 4);
            REQUIRE(components.group_size(2) == 4);
            REQUIRE(components.group_size(3) == 4);
        }
    }
    
    SECTION("A graph with a bridge edge") {
        // This is two k-4s connected by a single edge
        adjacencies = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2, 7},
                       {5, 6, 7}, {4, 6, 7}, {4, 5, 7}, {4, 5, 6, 3}};
        components = structures::UnionFind(adjacencies.size(), true);
        
        SECTION("Works with Cactus") {
            algorithms::three_edge_connected_components_dense_cactus(adjacencies.size(), for_each_connected_node, component_callback);
            
            REQUIRE(components.all_groups().size() == 2);
            REQUIRE(components.group_size(0) == 4);
            REQUIRE(components.group_size(1) == 4);
            REQUIRE(components.group_size(2) == 4);
            REQUIRE(components.group_size(3) == 4);
            
            REQUIRE(components.group_size(4) == 4);
            REQUIRE(components.group_size(5) == 4);
            REQUIRE(components.group_size(6) == 4);
            REQUIRE(components.group_size(7) == 4);
        }
        
        SECTION("Works with Tsin 2014") {
            algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
            REQUIRE(components.all_groups().size() == 2);
            REQUIRE(components.group_size(0) == 4);
            REQUIRE(components.group_size(1) == 4);
            REQUIRE(components.group_size(2) == 4);
            REQUIRE(components.group_size(3) == 4);
            
            REQUIRE(components.group_size(4) == 4);
            REQUIRE(components.group_size(5) == 4);
            REQUIRE(components.group_size(6) == 4);
            REQUIRE(components.group_size(7) == 4);
        }
    }
    
    SECTION("The Tsin 2007 paper graph works in Tsin 2014") {
        // We reverse the adjacency lists because we reaverse them back to
        // front in our DFS, and we want to match the paper's example.
        // We add in a separate node 0 component so we can use the 1-based IDs in the paper.
        adjacencies = {
            {}, //0
            {10, 10, 2}, //1
            {1, 3, 8}, //2
            {5, 4, 2}, //3
            {6, 6, 3}, //4
            {6, 7, 6, 3}, //5
            {4, 4, 5, 5}, //6
            {5, 17, 11, 17, 8, 12}, //7
            {3, 7, 9}, //8
            {8, 10}, //9
            {1, 9, 1}, //10
            {12, 7, 17}, //11
            {16, 13, 7, 11}, //12
            {14, 15, 12, 16}, //13
            {15, 13, 16}, //14
            {13, 16, 14}, //15
            {13, 14, 15, 12}, //16
            {7, 11, 7} //17
        };
        components = structures::UnionFind(adjacencies.size(), true);
        
        algorithms::three_edge_connected_components_dense(adjacencies.size(), 3, for_each_connected_node, component_callback);
            
        REQUIRE(components.group_size(1) == 2);
        REQUIRE(components.group_size(10) == 2);
        REQUIRE(components.find_group(1) == components.find_group(10));
        
        REQUIRE(components.group_size(2) == 2);
        REQUIRE(components.group_size(8) == 2);
        REQUIRE(components.find_group(2) == components.find_group(8));
        
        REQUIRE(components.group_size(3) == 4);
        REQUIRE(components.group_size(4) == 4);
        REQUIRE(components.group_size(5) == 4);
        REQUIRE(components.group_size(6) == 4);
        REQUIRE(components.find_group(3) == components.find_group(4));
        REQUIRE(components.find_group(3) == components.find_group(5));
        REQUIRE(components.find_group(3) == components.find_group(6));
        
        REQUIRE(components.group_size(7) == 3);
        REQUIRE(components.group_size(11) == 3);
        REQUIRE(components.group_size(17) == 3);
        REQUIRE(components.find_group(7) == components.find_group(11));
        REQUIRE(components.find_group(7) == components.find_group(17));
        
        REQUIRE(components.group_size(9) == 1);
        
        REQUIRE(components.group_size(12) == 1);
        
        REQUIRE(components.group_size(13) == 4);
        REQUIRE(components.group_size(14) == 4);
        REQUIRE(components.group_size(15) == 4);
        REQUIRE(components.group_size(16) == 4);
        REQUIRE(components.find_group(13) == components.find_group(14));
        REQUIRE(components.find_group(13) == components.find_group(15));
        REQUIRE(components.find_group(13) == components.find_group(16));
    }
}

TEST_CASE("3ECC algorithms do not over-collapse an extra-edge triangle", "[3ecc][algorithms]") {
    
    adjacencies = {{2, 2, 1}, {2, 0}, {1, 0, 0}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    SECTION("Works with Cactus") {
        algorithms::three_edge_connected_components_dense_cactus(adjacencies.size(), for_each_connected_node, component_callback);
    
        // Only two things should merge.
        REQUIRE(components.all_groups().size() == 2);
    }
    
    SECTION("Works with Tsin 2014") {
        algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
        for (auto& group : components.all_groups()) {
            cerr << "Group:";
            for (auto& member : group) {
                cerr << " " << member;
            }
            cerr << endl;
        }
#endif
        
        // Only two things should merge.
        REQUIRE(components.all_groups().size() == 2);
        
    }
}

TEST_CASE("Tsin 2014 does not over-collapse in the presence of bridge edges", "[3ecc][algorithms]") {
    adjacencies = {{2, 1}, {3, 2, 0}, {1, 0}, {1}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
#endif
    
    // Node 0 should not merge
    REQUIRE(components.group_size(0) == 1);
    // In fact nothing should merge
    REQUIRE(components.all_groups().size() == adjacencies.size());
}

TEST_CASE("Tsin 2014 does not over-collapse in the presence of bridge edges with self loops", "[3ecc][algorithms]") {
    adjacencies = {{2, 1}, {3, 2, 0}, {1, 0}, {1, 3}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
#endif
    
    // Node 0 should not merge
    REQUIRE(components.group_size(0) == 1);
    // In fact nothing should merge
    REQUIRE(components.all_groups().size() == adjacencies.size());
}

TEST_CASE("Tsin 2014 does not over-collapse in the presence of bridge edges with multiple self loops", "[3ecc][algorithms]") {
    adjacencies = {{2, 1}, {3, 2, 0}, {1, 0}, {1, 3, 3}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
#endif
    
    // Node 0 should not merge
    REQUIRE(components.group_size(0) == 1);
    // In fact nothing should merge
    REQUIRE(components.all_groups().size() == adjacencies.size());
}

TEST_CASE("Tsin 2014 does not over-collapse in the presence of bridge sticks with self loops", "[3ecc][algorithms]") {
    adjacencies = {{2, 1}, {3, 2, 0}, {1, 0}, {1, 4}, {3, 5}, {4, 5}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
#endif
    
    // Node 0 should not merge
    REQUIRE(components.group_size(0) == 1);
    // In fact nothing should merge
    REQUIRE(components.all_groups().size() == adjacencies.size());
}

TEST_CASE("Tsin 2014 does not over-collapse in the presence of bridge edges with single-node loops", "[3ecc][algorithms]") {
    adjacencies = {{2, 1}, {3, 2, 0}, {1, 0}, {1, 4, 4}, {3, 3}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
#endif
    
    // Node 0 should not merge
    REQUIRE(components.group_size(0) == 1);
    // In fact nothing should merge
    REQUIRE(components.all_groups().size() == adjacencies.size());
}

TEST_CASE("Tsin 2014 does not over-collapse in the presence of bridge sticks with single-node loops", "[3ecc][algorithms]") {
    adjacencies = {{2, 1}, {3, 2, 0}, {1, 0}, {1, 4}, {3, 5}, {4, 6, 6}, {5, 5}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
#endif
    
    // Node 0 should not merge
    REQUIRE(components.group_size(0) == 1);
    // In fact nothing should merge
    REQUIRE(components.all_groups().size() == adjacencies.size());
}

TEST_CASE("Tsin 2014 does not over-collapse a run of edge pairs", "[3ecc][algorithms]") {
    adjacencies = {{1}, {0, 2, 2}, {1, 1, 3, 3}, {2, 2, 3, 4, 4}, {3, 3, 5, 5}, {4, 4, 5}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
#endif
    
    // Node 5 should not merge
    REQUIRE(components.group_size(5) == 1);
    // In fact nothing should merge
    REQUIRE(components.all_groups().size() == adjacencies.size());
}

TEST_CASE("Tsin 2014 handles a graph with self loops and extra-edge triangles", "[3ecc][algorithms]") {
    adjacencies = {{4, 1, 2}, {0}, {5, 0, 3}, {2, 4, 4}, {3, 3, 0}, {2, 5, 6}, {5}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
#ifdef debug
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
#endif
    
    // Only two things should merge.
    REQUIRE(components.all_groups().size() == 6);
}

TEST_CASE("Tsin 2014 works correctly on random graphs", "[3ecc][algorithms]") {
    
    for (size_t node_count = 2; node_count <= 20; node_count++) {
        for (size_t edge_count = 0; edge_count <= node_count * 3; edge_count += node_count/2) {
        
            for (size_t repeat = 0; repeat < 10; repeat++) {
        
        
#ifdef debug
                cerr << "Making random graph with " << node_count << " nodes and " << edge_count << " edges" << endl;
#endif
            
                // Make a random graph
                adjacencies = random_adjacency_list(node_count, edge_count);
                components = structures::UnionFind(adjacencies.size(), true);
            
#ifdef debug
                cerr << "Finding 3 edge connected components..." << endl;
#endif
            
                // Find 3ecc with Tsin 2014
                algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
                
#ifdef debug
                cerr << "Finding true 3 edge connected components..." << endl;
#endif
                
                // Find 3ecc manually
                structures::UnionFind truth(brute_force_3ecc(adjacencies.size(), for_each_connected_node));
                
                // Check
                bool correct = uf_equal(components, truth);
                
                if (!correct) {
                    // Dump instance
                    
#ifdef debug
                }
                {
#endif
                    
                    cerr << "Graph:" << endl;
                    for (size_t i = 0; i < adjacencies.size(); i++) {
                        cerr << i << ":";
                        for (auto& j : adjacencies[i]) {
                            cerr << " " << j;
                        }
                        cerr << endl;
                    }

                    cerr << "Our algorithm:" << endl;
                    for (auto& group : components.all_groups()) {
                        cerr << "Group:";
                        for (auto& member : group) {
                            cerr << " " << member;
                        }
                        cerr << endl;
                    }
                    
                    cerr << "Truth:" << endl;
                    for (auto& group : truth.all_groups()) {
                        cerr << "Group:";
                        for (auto& member : group) {
                            cerr << " " << member;
                        }
                        cerr << endl;
                    }
                }
                
                REQUIRE(correct);
            }
        }
    }
}

}
}
        
