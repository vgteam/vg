/// \file cactus.cpp
///  
/// Unit tests for Cactus graphs and conversions.
///

#include <iostream>
#include <string>
#include "../io/json2graph.hpp"
#include "../cactus.hpp"
#include <bdsg/hash_graph.hpp>
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("We can convert a two-tailed graph to Cactus", "[cactus]") {

    string graph_json = R"(
    {"node":[{"sequence":"GT","id":7575},
    {"sequence":"TGTTAACAGCACAACATTTA","id":7580},
    {"sequence":"GGAAGTGTTTGCTACCAA","id":7576}],
    "path":[],
    "edge":[{"from":7575,"to":7580,"from_start":true},
    {"from":7575,"to":7576}]}
    )";

    bdsg::HashGraph graph;
    vg::io::json2graph(graph_json, &graph);

    // Make sure we can make a Cactus graph and get something out.
    auto cactusified = cactusify(graph);
    REQUIRE(cactusified.is_valid());
    
}

TEST_CASE("We can convert a hairpin graph to Cactus", "[cactus]") {

    // Here's a graph where only the left side of node 2 is dangling, and the right side of node 1 has a self loop.
    string graph_json = R"(
    {"node": [{"sequence": "A", "id": 1},
    {"sequence": "C", "id": 2}],
    "edge": [{"from": 2, "to": 1},
    {"from": 1, "to": 1, "to_end": true}]}
    )";

    bdsg::HashGraph graph;
    vg::io::json2graph(graph_json, &graph);

    // Make sure we can make a Cactus graph and get something out.
    auto cactusified = cactusify(graph);
    REQUIRE(cactusified.is_valid());
}

}
}
        
