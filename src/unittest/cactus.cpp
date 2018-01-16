/// \file cactus.cpp
///  
/// Unit tests for Cactus graphs and conversions.
///

#include <iostream>
#include <string>
#include "../json2pb.h"
#include "../cactus.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("We can convert a two-tailed graph to Cactus and back", "[cactus]") {
    
    VG graph;
    
    string graph_json = R"(
    {"node":[{"sequence":"GT","id":7575},
    {"sequence":"TGTTAACAGCACAACATTTA","id":7580},
    {"sequence":"GGAAGTGTTTGCTACCAA","id":7576}],
    "path":[],
    "edge":[{"from":7575,"to":7580,"from_start":true},
    {"from":7575,"to":7576}]}
    )";
    
    Graph g;
    json2pb(g, graph_json.c_str(), graph_json.size());
    graph.extend(g);

    // Make sure we can make a Cactus graph and get something out.    
    auto cactusified = cactusify(graph);
    REQUIRE(cactusified.is_valid());
    
}

}
}
        
