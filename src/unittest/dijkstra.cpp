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
    SECTION("Shortest path across snarl is correct") {
        size_t i = 0;
        handlealgs::for_each_handle_in_shortest_path(&graph, start, end, [&](handle_t next, size_t distance) {
            cerr << graph.get_id(next) << " " << graph.get_is_reverse(next) << endl;
            if (i == 0) {
                REQUIRE(next == middle);
                REQUIRE(distance == 0);
            } else if (i == 1) {
                REQUIRE((next == snp1 || next == snp2));
                REQUIRE(distance == 2);
            } else {
                REQUIRE(false);
            }
            i++;
            return true;
        });
        REQUIRE( i == 2);
    }
    SECTION("Shortest path from within snarl is correct") {
        size_t i = 0;
        handlealgs::for_each_handle_in_shortest_path(&graph, middle, end, [&](handle_t next, size_t distance) {
            if (i == 0) {
                REQUIRE((next == snp1 || next == snp2));
                REQUIRE(distance == 0);
            } else {
                REQUIRE(false);
            }
            i++;
            return true;
        });
        REQUIRE( i == 1);
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
TEST_CASE( "Shortest path through chain with loop", "[dijkstra][algorithms]" ) {
        
    VG graph;

    handle_t n1 = graph.create_handle("GCA");
    handle_t n2 = graph.create_handle("T");
    handle_t n3 = graph.create_handle("GGCTGACTGA");
    handle_t n4 = graph.create_handle("CTGA");
    handle_t n5 = graph.create_handle("GCA");
    handle_t n6 = graph.create_handle("T");
    handle_t n7 = graph.create_handle("G");
    handle_t n8 = graph.create_handle("CTGA");
    handle_t n9 = graph.create_handle("GCA");
    handle_t n10 = graph.create_handle("T");
    handle_t n11 = graph.create_handle("G");
    handle_t n12 = graph.create_handle("CTGA");
    handle_t n13 = graph.create_handle("GCA");
 
    
    graph.create_edge(n1, n2);
    graph.create_edge(graph.flip(n2), n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(graph.flip(n3), n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    graph.create_edge(n5, n7);
    graph.create_edge(n6, n12);
    graph.create_edge(n7, n10);
    graph.create_edge(n7, n8);
    graph.create_edge(n7, graph.flip(n8));
    graph.create_edge(n7, n9);
    graph.create_edge(n9, n10);
    graph.create_edge(n10, n11);
    graph.create_edge(n11, n12);
    graph.create_edge(n12, n13);
    
    SECTION("Shortest path across chain is correct") {
        vector<pair<handle_t, size_t>> actual_path;
        actual_path.emplace_back(n4, 0);
        actual_path.emplace_back(n5, 4);
        actual_path.emplace_back(n7, 7);
        actual_path.emplace_back(n10, 8);
        size_t i = 0;
        handlealgs::for_each_handle_in_shortest_path(&graph, n2, n11, [&](handle_t next, size_t distance) {
            REQUIRE( i < 5);
            REQUIRE(next == actual_path[i].first);
            REQUIRE(distance == actual_path[i].second);
            i++;
            return true;
        });
        REQUIRE( i == 4);
    }
    SECTION("Shortest path taking loop is correct") {
        vector<pair<handle_t, size_t>> actual_path;
        actual_path.emplace_back(n4, 0);
        actual_path.emplace_back(n5, 4);
        actual_path.emplace_back(n7, 7);
        actual_path.emplace_back(n8, 8);
        actual_path.emplace_back(graph.flip(n7), 12);
        size_t i = 0;
        handlealgs::for_each_handle_in_shortest_path(&graph, n2, graph.flip(n5), [&](handle_t next, size_t distance) {
            cerr << "iterated on " << graph.get_id(next) << " " << graph.get_is_reverse(next) << endl;
            REQUIRE( i < 5);
            if (i == 3) {
                REQUIRE((next == actual_path[i].first || graph.flip(next) == actual_path[i].first));
            } else {
                REQUIRE(next == actual_path[i].first);
            }
            REQUIRE(distance == actual_path[i].second);
            i++;
            return true;
        });
        REQUIRE( i == 5);
    }
 }

    

   
}
}
        
