/// \file unittest/xdrop_aligner.cpp
///  
/// Unit tests for the XdropAligner class.
///

#include <iostream>
#include <string>
#include "../json2pb.h"
#include <vg/vg.pb.h>
#include "../vg.hpp"
#include "../xdrop_aligner.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("XdropAligner can compute an alignment with no MEMs", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 0, 40);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCTGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    
    vector<MaximalExactMatch> no_mems;
    
    aligner.align(aln, graph.graph, no_mems, false);
    
    // Make sure we got the right score
    REQUIRE(aln.score() == read.size());
    
    // Make sure we take the right path
    REQUIRE(aln.path().mapping_size() == 3);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
    REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
}

TEST_CASE("XdropAligner can compute an alignment with a MEM in the middle", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 0, 40);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCTGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    
    vector<MaximalExactMatch> fake_mems;
    fake_mems.emplace_back();
    // Claim a match on the "GT"
    fake_mems.back().begin = aln.sequence().begin() + 1;
    fake_mems.back().end = aln.sequence().begin() + 3;
    fake_mems.back().nodes.push_back(gcsa::Node::encode(n0->id(), 1, false));
    
    aligner.align(aln, graph.graph, fake_mems, false);
    
    // Make sure we got the right score
    REQUIRE(aln.score() == read.size());
    
    // Make sure we take the right path
    REQUIRE(aln.path().mapping_size() == 3);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
    REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
}

TEST_CASE("XdropAligner still incorrectly applies the full length bonus at only one end with a MEM", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 10, 40);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCTGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    
    vector<MaximalExactMatch> fake_mems;
    fake_mems.emplace_back();
    // Claim a match on the "GT"
    fake_mems.back().begin = aln.sequence().begin() + 1;
    fake_mems.back().end = aln.sequence().begin() + 3;
    fake_mems.back().nodes.push_back(gcsa::Node::encode(n0->id(), 1, false));
    
    aligner.align(aln, graph.graph, fake_mems, false);
    
    // Make sure we got the right score
    size_t expected_score = read.size() + 10 * 1;
    REQUIRE(aln.score() == expected_score);
}

TEST_CASE("XdropAligner still incorrectly applies the full length bonus at only one end with no MEM", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 10, 40);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCTGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    
    vector<MaximalExactMatch> no_mems;
    
    aligner.align(aln, graph.graph, no_mems, false);
    
    size_t expected_score = read.size() + 10 * 1;
    REQUIRE(aln.score() == expected_score);
}


   
}
}
        
