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
    
    aligner.align(aln, graph, no_mems, false);
    
    // Make sure we got the right score
    REQUIRE(aln.score() == read.size());
    
    // Make sure we take the right path
    REQUIRE(aln.path().mapping_size() == 3);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
    REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
}

TEST_CASE("XdropAligner can compute an alignment with no MEMs in reverse mode", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 0, 40);
   
    // Build the graph in normal topological order
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAAGT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    // This should align to the reverse strand of the graph
    string read = string("AGTGCTGAAGT");
    Alignment aln;
    aln.set_sequence(read);
    
    vector<MaximalExactMatch> no_mems;
    
    // All this really changes is that we anchor on the last node instead of the first.
    aligner.align(aln, graph, no_mems, true);
    
    // Make sure we got the right score
    REQUIRE(aln.score() == read.size());
    
    // Make sure we take the right path
    REQUIRE(aln.path().mapping_size() == 3);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(0).position().is_reverse() == false);
    REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
    REQUIRE(aln.path().mapping(1).position().is_reverse() == false);
    REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
    REQUIRE(aln.path().mapping(2).position().is_reverse() == false);
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
    
    aligner.align(aln, graph, fake_mems, false);
    
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
    
    aligner.align(aln, graph, fake_mems, false);
    
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
    
    aligner.align(aln, graph, no_mems, false);
    
    size_t expected_score = read.size() + 10 * 1;
    REQUIRE(aln.score() == expected_score);
}

TEST_CASE("XdropAligner can be induced to pin with MEMs", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 0, 40);
    
    Node* n0 = graph.create_node("GAAAAAAAAAAAAAAAAAAAAA");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGATTACAT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("GATTACA");
    Alignment aln;
    aln.set_sequence(read);
    
    SECTION("The optimal alignment is foud without a MEM") {
        vector<MaximalExactMatch> no_mems;
        
        aligner.align(aln, graph, no_mems, false);
    
        // Make sure we got the right score
        REQUIRE(aln.score() == read.size());
        
        // Make sure we take the right path
        REQUIRE(aln.path().mapping_size() == 1);
        REQUIRE(aln.path().mapping(0).position().node_id() == n3->id());
    }
    
    SECTION("The suboptimal MEM-consistent alignment is foud with a MEM") {
    
        vector<MaximalExactMatch> fake_mems;
        fake_mems.emplace_back();
        // Claim a match on the initial "G" only
        fake_mems.back().begin = aln.sequence().begin();
        fake_mems.back().end = aln.sequence().begin() + 1;
        fake_mems.back().nodes.push_back(gcsa::Node::encode(n0->id(), 0, false));
        
        aligner.align(aln, graph, fake_mems, false);
        
        // The score will probably be terrible.
    
        // Make sure we land on the node the MEM was on, even though it is a terrible alignment.
        REQUIRE(aln.path().mapping_size() == 1);
        REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    }
}

TEST_CASE("XdropAligner can align pinned left", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 10, 40);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAACT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("ACT");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    aligner.align_pinned(aln, graph, true);
    
    // Make sure we got the right score.
    // Account for full length bonus, loss of a match, and gain of a mismatch.
    REQUIRE(aln.score() == read.size() + 10 - 1 - 4);
    
    // Make sure we take the right path
    REQUIRE(aln.path().mapping_size() == 1);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(0).position().offset() == 0);
    REQUIRE(aln.path().mapping(0).edit_size() == 3);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
    REQUIRE(aln.path().mapping(0).edit(1).from_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(1).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(1).sequence() == "C");
    REQUIRE(aln.path().mapping(0).edit(2).from_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(2).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(2).sequence() == "");
}

TEST_CASE("XdropAligner can align pinned right", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 10, 40);
    
    Node* n0 = graph.create_node("AGTG");
    Node* n1 = graph.create_node("C");
    Node* n2 = graph.create_node("A");
    Node* n3 = graph.create_node("TGAACT");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n3);
    
    string read = string("AGT");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned right, letting the graph compute a topological order
    aligner.align_pinned(aln, graph, false);
    
    // Make sure we got the right score.
    // Account for full length bonus, loss of a match, and gain of a mismatch.
    REQUIRE(aln.score() == read.size() + 10 - 1 - 4);
    
    // Make sure we take the right path
    REQUIRE(aln.path().mapping_size() == 1);
    REQUIRE(aln.path().mapping(0).position().node_id() == n3->id());
    REQUIRE(aln.path().mapping(0).position().offset() == 3);
    REQUIRE(aln.path().mapping(0).edit_size() == 3);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
    REQUIRE(aln.path().mapping(0).edit(1).from_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(1).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(1).sequence() == "G");
    REQUIRE(aln.path().mapping(0).edit(2).from_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(2).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(2).sequence() == "");
}


TEST_CASE("XdropAligner can align pinned left when that is a bad alignment", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 10, 40);
    
    Node* n0 = graph.create_node("TTAAGCTGAGGGAATAGTGCCTGGCATCGAGGAAAGCCTCTGA");
    
    string read = string("AGCTGAGGGAATAGTGCCTGGCATCGAGGAAAGCCTCTGA");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    aligner.align_pinned(aln, graph, true);
    
    // Make sure we got the right score.
    // Account for full length bonus, two extends, and one open
    REQUIRE(aln.score() == read.size() + 10 - 2 * 1 - 6);
    
    // Make sure we take the right path (leading 3 bp deletion)
    REQUIRE(aln.path().mapping_size() == 1);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(0).position().offset() == 0);
    REQUIRE(aln.path().mapping(0).edit_size() == 2);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == 3);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
    REQUIRE(aln.path().mapping(0).edit(1).from_length() == read.size());
    REQUIRE(aln.path().mapping(0).edit(1).to_length() == read.size());
    REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
}

TEST_CASE("XdropAligner can align pinned left with a leading insertion", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 10, 40);
    
    Node* n0 = graph.create_node("AAAGAGGTCAATAGCCAAAT");
    
    string read = string("GAAAGAGGTCAATAGCCAAAT");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    aligner.align_pinned(aln, graph, true);
    
    // Make sure we got the right score.
    // Account for full length bonus and one open, and the lack of a match on
    // the extra query base
    REQUIRE(aln.score() == read.size() - 1 + 10 - 6);
    
    // Make sure we take the right path (leading 1 bp insertion)
    REQUIRE(aln.path().mapping_size() == 1);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(0).position().offset() == 0);
    REQUIRE(aln.path().mapping(0).edit_size() == 2);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(0).sequence() == "G");
    REQUIRE(aln.path().mapping(0).edit(1).from_length() == read.size() - 1);
    REQUIRE(aln.path().mapping(0).edit(1).to_length() == read.size() - 1);
    REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
}

TEST_CASE("XdropAligner can align pinned right with a trailing insertion", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 10, 40);
    
    Node* n0 = graph.create_node("ATTTGGCTATTGACCTCTTT");
    
    string read = string("ATTTGGCTATTGACCTCTTTC");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned right, letting the graph compute a topological order
    aligner.align_pinned(aln, graph, false);
    
    // Make sure we got the right score.
    // Account for full length bonus and one open, and the lack of a match on
    // the extra query base
    REQUIRE(aln.score() == read.size() - 1 + 10 - 6);
    
    // Make sure we take the right path (trailing 1 bp insertion)
    REQUIRE(aln.path().mapping_size() == 1);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(0).position().offset() == 0);
    REQUIRE(aln.path().mapping(0).edit_size() == 2);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == read.size() - 1);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == read.size() - 1);
    REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
    REQUIRE(aln.path().mapping(0).edit(1).from_length() == 0);
    REQUIRE(aln.path().mapping(0).edit(1).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(1).sequence() == "C");

    // Make sure we got a rank set.
    REQUIRE(aln.path().mapping(0).rank() == 1);
}

TEST_CASE("XdropAligner can align pinned left when the entire read is an insertion", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    // Last parameter here is max gap length.
    XdropAligner aligner(1, 4, 6, 1, 10, 40);
    
    Node* n0 = graph.create_node("A");
    
    // Not even the full length bonus can save us
    string read = string("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    aligner.align_pinned(aln, graph, true);
    
    // Make sure we got the right score.
    // The whole sequence should just softclip.
    REQUIRE(aln.score() == 0);
    
    // The sequence should stay set
    REQUIRE(aln.sequence() == read);
    
    // Make sure we take the right path (whole read inserted)
    REQUIRE(aln.path().mapping_size() == 1);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(0).position().offset() == 0);
    REQUIRE(aln.path().mapping(0).edit_size() == 1);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == read.size());
    REQUIRE(aln.path().mapping(0).edit(0).sequence() == read);

    // Make sure we got a rank set.
    REQUIRE(aln.path().mapping(0).rank() == 1);
}


   
}
}
        
