/// \file dinucleotide_machine.cpp
///  
/// Unit tests for the DinucleotideMachine
///

#include <iostream>
#include <random>

#include "../splice_region.hpp"
#include "catch.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {
using namespace std;

using bdsg::HashGraph;

TEST_CASE("SpliceRegion can detect a splice site on a trivial example",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h = graph.create_handle("AGTTGCAT");
    
    DinucleotideMachine machine;
    
    pos_t pos(graph.get_id(h), false, 3);
    bool search_left = false;
    int64_t search_dist = 2;
    vector<string> motifs{"GC", "CA", "GT"};
    
    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, motifs);
    
    auto m0 = splice_region.candidate_splice_sites(motifs[0]);
    auto m1 = splice_region.candidate_splice_sites(motifs[1]);
    auto m2 = splice_region.candidate_splice_sites(motifs[2]);
    
    REQUIRE(m0.size() == 1);
    REQUIRE(m1.size() == 1);
    REQUIRE(m2.size() == 0);
    
    REQUIRE(m0.front().first == pos_t(graph.get_id(h), false, 6));
    REQUIRE(m0.front().second == 3);
    REQUIRE(m1.front().first == pos_t(graph.get_id(h), false, 7));
    REQUIRE(m1.front().second == 4);
    
}

TEST_CASE("SpliceRegion can detect a splice site leftward",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h = graph.create_handle("AGTTGCA");
    
    DinucleotideMachine machine;
    
    pos_t pos(graph.get_id(h), false, 6);
    bool search_left = true;
    int64_t search_dist = 3;
    vector<string> motifs{"GC", "CA", "GT"};
    
    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, motifs);
    
    auto m0 = splice_region.candidate_splice_sites(motifs[0]);
    auto m1 = splice_region.candidate_splice_sites(motifs[1]);
    auto m2 = splice_region.candidate_splice_sites(motifs[2]);
    
    REQUIRE(m0.size() == 1);
    REQUIRE(m1.size() == 0);
    REQUIRE(m2.size() == 1);
    
    REQUIRE(m0.front().first == pos_t(graph.get_id(h), false, 4));
    REQUIRE(m0.front().second == 2);
    REQUIRE(m2.front().first == pos_t(graph.get_id(h), false, 1));
    REQUIRE(m2.front().second == 5);
    
}

TEST_CASE("SpliceRegion can detect a splice sites across node boundaries",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h0 = graph.create_handle("AAG");
    handle_t h1 = graph.create_handle("G");
    handle_t h2 = graph.create_handle("C");
    handle_t h3 = graph.create_handle("GTG");
    handle_t h4 = graph.create_handle("TTT");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h0, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    
    DinucleotideMachine machine;
    
    pos_t pos(graph.get_id(h0), false, 0);
    bool search_left = false;
    int64_t search_dist = 4;
    vector<string> motifs{"GC", "GG", "GT", "CG"};
    
    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, motifs);
    
    auto m0 = splice_region.candidate_splice_sites(motifs[0]);
    auto m1 = splice_region.candidate_splice_sites(motifs[1]);
    auto m2 = splice_region.candidate_splice_sites(motifs[2]);
    auto m3 = splice_region.candidate_splice_sites(motifs[3]);
    
    REQUIRE(m0.size() == 1);
    REQUIRE(m1.size() == 2);
    REQUIRE(m2.size() == 1);
    REQUIRE(m3.size() == 1);
    
    REQUIRE(m0.front().first == pos_t(graph.get_id(h2), false, 1));
    REQUIRE(m0.front().second == 4);
    REQUIRE(m1.front().first == pos_t(graph.get_id(h1), false, 1));
    REQUIRE(m1.front().second == 4);
    REQUIRE(m1.back().first == pos_t(graph.get_id(h3), false, 1));
    REQUIRE(m1.back().second == 5);
    REQUIRE(m2.front().first == pos_t(graph.get_id(h3), false, 2));
    REQUIRE(m2.front().second == 6);
    REQUIRE(m3.front().first == pos_t(graph.get_id(h3), false, 1));
    REQUIRE(m3.front().second == 5);
    
}



TEST_CASE("SpliceRegion can detect a splice sites across node boundaries going left",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h0 = graph.create_handle("GCG");
    handle_t h1 = graph.create_handle("G");
    handle_t h2 = graph.create_handle("C");
    handle_t h3 = graph.create_handle("AGCA");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h0, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    
    DinucleotideMachine machine;
    
    pos_t pos(graph.get_id(h3), false, 2);
    bool search_left = true;
    int64_t search_dist = 2;
    vector<string> motifs{"GC", "GG"};
    
    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, motifs);
    
    auto m0 = splice_region.candidate_splice_sites(motifs[0]);
    auto m1 = splice_region.candidate_splice_sites(motifs[1]);
    
    REQUIRE(m0.size() == 1);
    REQUIRE(m1.size() == 1);
    
    REQUIRE(m0.front().first == pos_t(graph.get_id(h0), false, 2));
    REQUIRE(m0.front().second == 4);
    REQUIRE(m1.front().first == pos_t(graph.get_id(h0), false, 2));
    REQUIRE(m1.front().second == 4);
    
}

}
}
        
