/// \file dijkstra.cpp
///  
/// Unit tests for the Dijkstra algorithm
///

#include <iostream>
#include <string>
#include "../algorithms/dijkstra.hpp"
#include "../handle.hpp"
#include "catch.hpp"

#include <sglib/hash_graph.hpp>


namespace vg {
namespace unittest {
using namespace std;

using sglib::HashGraph;

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
    
        algorithms::dijkstra(&graph, start, [&](const handle_t& reached, size_t distance) {
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
    
        algorithms::dijkstra(&graph, start, [&](const handle_t& reached, size_t distance) {
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


   
}
}
        
