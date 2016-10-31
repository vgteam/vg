/**
 * unittest/constructor.cpp: test cases for the vg graph constructor tool
 */

#include "catch.hpp"
#include "../constructor.hpp"

#include "../path.hpp"

#include <vector>

namespace vg {
namespace unittest {

TEST_CASE( "An empty chunk with no variants can be constructed", "[constructor]" ) {
    Constructor constructor;
    
    auto result = constructor.construct_chunk("", "empty", std::vector<vcflib::Variant>());
    
    SECTION("the graph should have no elements") {
        REQUIRE(result.graph.node_size() == 0);
        REQUIRE(result.graph.edge_size() == 0);
        REQUIRE(result.left_ends.empty());
        REQUIRE(result.right_ends.empty());
    }

}

TEST_CASE( "A small linear chunk with no variants can be constructed", "[constructor]" ) {
    Constructor constructor;
    
    auto result = constructor.construct_chunk("GATTACA", "movie", std::vector<vcflib::Variant>());
    
    SECTION("the graph should have one node") {
        REQUIRE(result.graph.node_size() == 1);
        auto& node = result.graph.node(0);
        
        SECTION("the node should have the full sequence") {
            REQUIRE(node.sequence() == "GATTACA");
        }
        
        SECTION("the node should have the correct ID") {
            REQUIRE(node.id() == 1);
        }
        
        SECTION("the node should be the only exposed node on the left") {
            REQUIRE(result.left_ends.count(node.id()));
            REQUIRE(result.left_ends.size() == 1);
        }
        
        SECTION("the node should be the only exposed node on the right") {
            REQUIRE(result.right_ends.count(node.id()));
            REQUIRE(result.right_ends.size() == 1);
        }
    }
    
    SECTION("the graph should have no edges") {
        REQUIRE(result.graph.edge_size() == 0);
    }
    
    SECTION("the graph should have one path") {
        REQUIRE(result.graph.path_size() == 1);
        
        auto& path = result.graph.path(0);
        
        SECTION("the path should have the name we passed in") {
            REQUIRE(path.name() == "movie");
        }
        
        SECTION("the path should have one mapping") {
            REQUIRE(path.mapping_size() == 1);
            
            auto& mapping = path.mapping(0);
            
            SECTION("the mapping should be a perfect match") {
                REQUIRE(mapping_is_match(mapping));
            }
        }
    }
}

// TODO: divide a linear graph into multiple nodes

// TODO: graphs with variants

}
}
