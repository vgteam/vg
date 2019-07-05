//
//  \file next_path_position.cpp
//  
// Tests for the next_path_position algorithm.
//

#include "catch.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../json2pb.h"
#include "../algorithms/next_path_position.hpp"
#include <vg/vg.pb.h>
#include <stdio.h>

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("next_path_position can search out for several nodes on an XG", "[xg][algorithms]") {

    string graph_json = R"(
    {"node":[{"id":1,"sequence":"G"},
    {"id":2,"sequence":"A"},
    {"id":3,"sequence":"T"},
    {"id":4,"sequence":"T"},
    {"id":5,"sequence":"A"},
    {"id":6,"sequence":"C"},
    {"id":7,"sequence":"A"}],
    "edge":[{"from":1,"to":2},
    {"from":2,"to":3},
    {"from":3,"to":4},
    {"from":4,"to":5},
    {"from":5,"to":6},
    {"from":6,"to":7}],
    "path":[{"name":"x","mapping":[{"position":{"node_id":7}, "edit":[{"from_length":1,"to_length":1}]}]}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Build the xg index
    XG xg_index(proto_graph);
    
    SECTION("nothing is found when searching a zero distance") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, false, 0), 0);
        
        REQUIRE(id(result.first) == 0);
    }
    
    SECTION("nothing is found when searching a too-short distance") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, false, 0), 2);
        
        REQUIRE(id(result.first) == 0);
    }
    
    SECTION("nothing is found when searching a slightly too-short distance") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, false, 0), 5);
        REQUIRE(id(result.first) == 0);
    }
    
    SECTION("something is found at the right distance when searching an exactly sufficiently long distance") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, false, 0), 6);
        
        REQUIRE(id(result.first) == 7);
        REQUIRE(result.second == 6);
    }
    
    SECTION("something is found at the rigth distance when searching an overly-long distance") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, false, 0), 100);
        
        REQUIRE(id(result.first) == 7);
        REQUIRE(result.second == 6);
    }
    
    SECTION("nothing is found when searching a zero distance in reverse") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, true, 1), 0);
        
        REQUIRE(id(result.first) == 0);
    }
    
    SECTION("nothing is found when searching a too-short distance in reverse") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, true, 1), 2);
        
        REQUIRE(id(result.first) == 0);
    }
    
    SECTION("nothing is found when searching a slightly too-short distance in reverse") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, true, 1), 5);
        REQUIRE(id(result.first) == 0);
    }
    
    SECTION("something is found at the right distance when searching an exactly sufficiently long distance in reverse") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, true, 1), 6);
        
        REQUIRE(id(result.first) == 7);
        REQUIRE(result.second == 6);
    }
    
    SECTION("something is found at the right distance when searching an overly-long distance in reverse") {
        pair<pos_t, int64_t> result = algorithms::next_path_position(xg_index, make_pos_t(1, true, 1), 100);
        
        REQUIRE(id(result.first) == 7);
        REQUIRE(result.second == 6);
    }
}

}
}
