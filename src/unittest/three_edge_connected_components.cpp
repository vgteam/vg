/// \file three_edge_connected_components.cpp
///  
/// Unit tests for the three edge connected components algorithms
///

#include <iostream>
#include <vector>
#include "../algorithms/three_edge_connected_components.hpp"
#include "catch.hpp"

#include <structures/union_find.hpp>


namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("3 edge connected components algorithms handle basic cases", "[3ecc][algorithms]") {
    vector<vector<size_t>> adjacencies;
   
    auto for_each_connected_node = [&](size_t node, const function<void(size_t)>& iteratee) {
#ifdef debug
        cerr << "Asked for edges of node " << node << endl;
#endif
        for (auto& other : adjacencies.at(node)) {
            iteratee(other);
        }
    };
    
    // Represent the results as a union-find for checking
    structures::UnionFind components(adjacencies.size(), true);
    
    auto component_callback = [&](const function<void(const function<void(size_t)>&)>& for_each_member) {
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
    };
    
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

TEST_CASE("Tsin 2014 handles a graph with self loops", "[3ecc][algorithms]") {
    vector<vector<size_t>> adjacencies;
   
    auto for_each_connected_node = [&](size_t node, const function<void(size_t)>& iteratee) {
#ifdef debug
        cerr << "Asked for edges of node " << node << endl;
#endif
        for (auto& other : adjacencies.at(node)) {
            iteratee(other);
        }
    };
    
    // Represent the results as a union-find for checking
    structures::UnionFind components(adjacencies.size(), true);
    
    auto component_callback = [&](const function<void(const function<void(size_t)>&)>& for_each_member) {
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
    };
    
    adjacencies = {{4, 1, 2}, {0}, {5, 0, 3}, {2, 4, 4}, {3, 3, 0}, {2, 5, 6}, {5}};
    components = structures::UnionFind(adjacencies.size(), true);
    
    algorithms::three_edge_connected_components_dense(adjacencies.size(), 0, for_each_connected_node, component_callback);
            
    for (auto& group : components.all_groups()) {
        cerr << "Group:";
        for (auto& member : group) {
            cerr << " " << member;
        }
        cerr << endl;
    }
    
    // Only two things should merge.
    REQUIRE(components.all_groups().size() == 6);
        
    
}

}
}
        
