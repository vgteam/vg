// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

/// \file dijkstra.cpp
///  
/// Unit tests for the Dijkstra algorithm
///

#include <iostream>
#include <string>
#include "../handle.hpp"
#include "vg/io/json2pb.h"
#include "../vg.hpp"
#include "catch.hpp"

#include <vg/vg.pb.h>

#include <bdsg/hash_graph.hpp>


namespace vg {
namespace unittest {
using namespace std;

using bdsg::HashGraph;

TEST_CASE("Dijkstra search handles early stopping correctly", "[dijkstra][algorithms]") {
    
    HashGraph graph;
    
    handle_t start = graph.create_handle("GAT");
    handle_t middle = graph.create_handle("TA");
    handle_t snp1 = graph.create_handle("C");
    handle_t snp2 = graph.create_handle("T");
    handle_t end = graph.create_handle("A");
    
    graph.create_edge(start, middle);
    graph.create_edge(middle, snp1);
    graph.create_edge(middle, snp2);
    graph.create_edge(snp1, end);
    graph.create_edge(snp2, end);
    
    // Track what we reach and at what distance
    unordered_map<handle_t, size_t> seen;
    
    SECTION("Does not hit all handles at a given distance when told to stop at one") {
    
        handlealgs::dijkstra(&graph, start, [&](const handle_t& reached, size_t distance) {
            seen[reached] = distance;
            if (reached == snp1 || reached == snp2) {
                // Stop at either branch of the SNP
                return false;
            }
            return true;
        });
        
        // The handles up to the SNP and one branch of the SNP should be captured.
        REQUIRE(seen.size() == 3);
        REQUIRE(seen.count(start));
        REQUIRE(seen.count(middle));
        REQUIRE(seen.count(snp1) != seen.count(snp2));
        REQUIRE(!seen.count(end));
        
        // The distance calculation should be from the end of the start handle
        REQUIRE(seen.at(middle) == 0);
        
    }
    
    SECTION("Does hit all handles at a given distance when told to stop after") {
    
        handlealgs::dijkstra(&graph, start, [&](const handle_t& reached, size_t distance) {
            if (distance > 2) {
                // Stop after the SNP
                return false;
            }
            seen[reached] = distance;
            return true;
        });
        
        // The handles up to the SNP and one branch of the SNP should be captured.
        REQUIRE(seen.size() == 4);
        REQUIRE(seen.count(start));
        REQUIRE(seen.count(middle));
        REQUIRE(seen.count(snp1));
        REQUIRE(seen.count(snp2));
        REQUIRE(!seen.count(end));
        
        // The distance calculation should be from the end of the start handle
        REQUIRE(seen.at(snp1) == 2);
        REQUIRE(seen.at(snp2) == 2);
        
    }
    
}

TEST_CASE("Dijkstra search works on a particular problem graph", "[dijkstra][algorithms]") {

    string graph_json = R"(
{"node":[{"sequence":"A","id":"2454530"},{"sequence":"AGTGCTGGAGAGGATGTGGAGAAATAGGAAC","id":"2454529"},{"sequence":"C","id":"2454532"},{"sequence":"TTTTACACTGTTGGTGGGACTGTAAA","id":"2454533"},{"sequence":"A","id":"2454527"},{"sequence":"C","id":"2454528"},{"sequence":"G","id":"2454531"},{"sequence":"C","id":"2454534"},{"sequence":"T","id":"2454535"},{"sequence":"GGGTAATAA","id":"2454526"},{"sequence":"TAGTTCAACCATTGTGGAAGACTGTGGCAATT","id":"2454536"}],"edge":[{"from":"2454530","to":"2454532"},{"from":"2454530","to":"2454533"},{"from":"2454529","to":"2454530"},{"from":"2454529","to":"2454531"},{"from":"2454532","to":"2454533"},{"from":"2454533","to":"2454534"},{"from":"2454533","to":"2454535"},{"from":"2454527","to":"2454529"},{"from":"2454528","to":"2454529"},{"from":"2454531","to":"2454532"},{"from":"2454531","to":"2454533"},{"from":"2454534","to":"2454536"},{"from":"2454535","to":"2454536"},{"from":"2454526","to":"2454527"},{"from":"2454526","to":"2454528"}],"path":[{"name":"21","mapping":[{"position":{"node_id":"2454526"},"edit":[{"from_length":9,"to_length":9}],"rank":"3049077"},{"position":{"node_id":"2454528"},"edit":[{"from_length":1,"to_length":1}],"rank":"3049078"},{"position":{"node_id":"2454529"},"edit":[{"from_length":31,"to_length":31}],"rank":"3049079"},{"position":{"node_id":"2454531"},"edit":[{"from_length":1,"to_length":1}],"rank":"3049080"},{"position":{"node_id":"2454532"},"edit":[{"from_length":1,"to_length":1}],"rank":"3049081"},{"position":{"node_id":"2454533"},"edit":[{"from_length":26,"to_length":26}],"rank":"3049082"},{"position":{"node_id":"2454535"},"edit":[{"from_length":1,"to_length":1}],"rank":"3049083"},{"position":{"node_id":"2454536"},"edit":[{"from_length":32,"to_length":32}],"rank":"3049084"}]}]}    
    )";
    
    Graph g;
    json2pb(g, graph_json);
    
    // Wrap the graph in a HandleGraph
    VG graph(g);
    
    // Decide where to start
    handle_t start = graph.get_handle(2454536, true);
    
    // Track what we reach and at what distance
    unordered_map<handle_t, size_t> seen;
    
    
    handlealgs::dijkstra(&graph, start, [&](const handle_t& reached, size_t distance) {
        seen[reached] = distance;
        return true;
    });
    
    REQUIRE(seen.size() == graph.get_node_count());
        
}
    

   
}
}
        
