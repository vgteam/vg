/// \file aligner.cpp
///  
/// Unit tests for the basic methods of the Aligner class. See also:
/// pinned_alignment.cpp.
///

#include <iostream>
#include <string>
#include "../json2pb.h"
#include "../vg.pb.h"
#include "../gssw_aligner.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Aligner respects the full length bonus at both ends", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCTGAAGT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph);
    aligner_2.align(aln2, graph.graph);
    
    SECTION("bonus is collected at both ends") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Aligner respects the full length bonus for a single base read", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("G");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph);
    aligner_2.align(aln2, graph.graph);
    
    SECTION("bonus is collected twice even though both ends are one match") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
}

TEST_CASE("Aligner works when end bonus is granted to a match at the start of a node", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    // Make sure aligner runs
    aligner_1.align(aln1, graph.graph);
    aligner_2.align(aln2, graph.graph);
    
    SECTION("bonus is collected twice") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Full-length bonus can hold down the left end", "[aligner][alignment][mapping]") {
    VG graph;
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTGCTGAAGT");
    
    string read = string("AATGCTGAAGT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph);
    aligner_2.align(aln2, graph.graph);
    
    SECTION("left end is detatched without bonus") {
        REQUIRE(aln1.path().mapping_size() == 1);
        REQUIRE(aln1.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln1.path().mapping(0).position().offset() == 2);
        REQUIRE(aln1.path().mapping(0).edit_size() == 2);
        REQUIRE(aln1.path().mapping(0).edit(0).from_length() == 0);
        REQUIRE(aln1.path().mapping(0).edit(0).sequence() == "AA");
    }
    
    SECTION("left end is attached with bonus") {
        REQUIRE(aln2.path().mapping_size() == 1);
        REQUIRE(aln2.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln2.path().mapping(0).position().offset() == 0);
        REQUIRE(aln2.path().mapping(0).edit_size() == 3);
        REQUIRE(aln2.path().mapping(0).edit(0).from_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(0).to_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(0).sequence() == "");
    }
}

TEST_CASE("Full-length bonus can hold down the right end", "[aligner][alignment][mapping]") {
    VG graph;
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTGCTGAAGT");
    
    string read = string("AGTGCTGAAAT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph);
    aligner_2.align(aln2, graph.graph);
    
    SECTION("right end is detatched without bonus") {
        REQUIRE(aln1.path().mapping_size() == 1);
        REQUIRE(aln1.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln1.path().mapping(0).position().offset() == 0);
        REQUIRE(aln1.path().mapping(0).edit_size() == 2);
        REQUIRE(aln1.path().mapping(0).edit(1).from_length() == 0);
        REQUIRE(aln1.path().mapping(0).edit(1).sequence() == "AT");
    }
    
    SECTION("right end is attached with bonus") {
        REQUIRE(aln2.path().mapping_size() == 1);
        REQUIRE(aln2.path().mapping(0).position().node_id() == n0->id());
        REQUIRE(aln2.path().mapping(0).position().offset() == 0);
        REQUIRE(aln2.path().mapping(0).edit_size() == 3);
        REQUIRE(aln2.path().mapping(0).edit(2).from_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(2).to_length() == 1);
        REQUIRE(aln2.path().mapping(0).edit(2).sequence() == "");
    }
}

TEST_CASE("Full-length bonus can attach Ns", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("NNNNCTGANNN");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph);
    aligner_2.align(aln2, graph.graph);
    
    SECTION("bonused alignment ends in full-length match/mismatches") {
        REQUIRE(aln2.path().mapping_size() == 3);
        REQUIRE(mapping_from_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_to_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_from_length(aln2.path().mapping(2)) == 6);
        REQUIRE(mapping_to_length(aln2.path().mapping(2)) == 6);
    }
    
    SECTION("bonus is collected at both ends") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Full-length bonus can attach to Ns", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("NNNG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGANNN");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGTGCTGAAGT");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph);
    aligner_2.align(aln2, graph.graph);
    
    SECTION("bonused alignment ends in full-length match/mismatches") {
        REQUIRE(aln2.path().mapping_size() == 3);
        REQUIRE(mapping_from_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_to_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_from_length(aln2.path().mapping(2)) == 6);
        REQUIRE(mapping_to_length(aln2.path().mapping(2)) == 6);
    }
    
    SECTION("bonus is collected at both ends") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}

TEST_CASE("Full-length bonus can attach Ns to Ns", "[aligner][alignment][mapping]") {
    
    VG graph;
    
    Aligner aligner_1(1, 4, 6, 1, 0);
    Aligner aligner_2(1, 4, 6, 1, 10);
    
    Node* n0 = graph.create_node("NNNG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGANNN");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("NNNGCTGANNN");
    Alignment aln1, aln2;
    aln1.set_sequence(read);
    aln2.set_sequence(read);
    
    aligner_1.align(aln1, graph.graph);
    aligner_2.align(aln2, graph.graph);
    
    SECTION("bonused alignment ends in full-length match/mismatches") {
        REQUIRE(aln2.path().mapping_size() == 3);
        REQUIRE(mapping_from_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_to_length(aln2.path().mapping(0)) == 4);
        REQUIRE(mapping_from_length(aln2.path().mapping(2)) == 6);
        REQUIRE(mapping_to_length(aln2.path().mapping(2)) == 6);
    }
    
    SECTION("bonus is collected at both ends") {
        REQUIRE(aln2.score() == aln1.score() + 20);
    }
    
}
   
}
}
        
