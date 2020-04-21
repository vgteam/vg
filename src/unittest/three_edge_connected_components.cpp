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
        cerr << "Asked for edges of node " << node << endl;
        for (auto& other : adjacencies.at(node)) {
            iteratee(other);
        }
    };
    
    // Represent the results as a union-find for checking
    structures::UnionFind components(adjacencies.size(), true);
    
    auto component_callback = [&](const function<void(const function<void(size_t)>&)>& for_each_member) {
        cerr << "Got component" << endl;
        size_t first = 0;
        bool is_first = true;
        for_each_member([&](size_t member) {
            cerr << "Component contained " << member << endl;
            if (is_first) {
                // Find the first member of each component
                first = member;
            } else {
                // And union everything into it
                components.union_groups(first, member);
            }
        });
    };
    
    SECTION("Cactus algorithm handles an empty graph") {
        algorithms::three_edge_connected_components_dense_cactus(adjacencies.size(), for_each_connected_node, component_callback);
        REQUIRE(components.size() == 0);
    }
    
    SECTION("Cactus algorithm handles a 4-node connected graph as one component") {
        adjacencies = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
        components = structures::UnionFind(adjacencies.size(), true);
    
        algorithms::three_edge_connected_components_dense_cactus(adjacencies.size(), for_each_connected_node, component_callback);
        
        REQUIRE(components.size() == 4);
        REQUIRE(components.group_size(0) == 4);
        REQUIRE(components.group_size(1) == 4);
        REQUIRE(components.group_size(2) == 4);
        REQUIRE(components.group_size(3) == 4);
    }
}

}
}
        
