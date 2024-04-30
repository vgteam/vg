/// \file unittest/xdrop_aligner.cpp
///  
/// Unit tests for the XdropAligner class.
///

#include <iostream>
#include <string>
#include "vg/io/json2pb.h"
#include "../alignment.hpp"
#include "../vg.hpp"
#include <vg/vg.pb.h>
#include "test_aligner.hpp"
#include "catch.hpp"
#include "bdsg/hash_graph.hpp"

namespace vg {
namespace unittest {
using namespace std;
using namespace vg::io;

TEST_CASE("XdropAligner can compute an alignment with no MEMs", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 0);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
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
    
    uint16_t max_gap_length =  40;
    aligner.align_xdrop(aln, graph, no_mems, false, max_gap_length);
    
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
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 0);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
   
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
    uint16_t max_gap_length =  40;
    aligner.align_xdrop(aln, graph, no_mems, true);
    
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
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 0);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
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
    
    uint16_t max_gap_length =  40;
    aligner.align_xdrop(aln, graph, fake_mems, false, max_gap_length);
    
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
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
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
    
    uint16_t max_gap_length =  40;
    aligner.align_xdrop(aln, graph, fake_mems, false, max_gap_length);
    
    // Make sure we got the right score
    size_t expected_score = read.size() + 10 * 1;
    REQUIRE(aln.score() == expected_score);
}

TEST_CASE("XdropAligner still incorrectly applies the full length bonus at only one end with no MEM", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
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
    
    uint16_t max_gap_length =  40;
    aligner.align_xdrop(aln, graph, no_mems, false, max_gap_length);
    
    size_t expected_score = read.size() + 10 * 1;
    REQUIRE(aln.score() == expected_score);
}

TEST_CASE("XdropAligner can be induced to pin with MEMs", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 0);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
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
        
        uint16_t max_gap_length =  40;
        aligner.align_xdrop(aln, graph, no_mems, false, max_gap_length);
    
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
        
        uint16_t max_gap_length =  40;
        aligner.align_xdrop(aln, graph, fake_mems, false, max_gap_length);
        
        // The score will probably be terrible.
    
        // Make sure we land on the node the MEM was on, even though it is a terrible alignment.
        REQUIRE(aln.path().mapping_size() == 1);
        REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    }
}

TEST_CASE("XdropAligner can align pinned left", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
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
    aligner.align_pinned(aln, graph, true, true);
    
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

TEST_CASE("XdropAligner can align pinned left across a large gap that occurrs immediately at the pin point, when the gap limit is high enough", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    auto h1 = graph.create_handle("CATCA");
    auto h2 = graph.create_handle("T");
    graph.create_edge(h1, h2);
    auto h3 = graph.create_handle("ATA");
    graph.create_edge(h2, h3);
    auto h4 = graph.create_handle("G");
    graph.create_edge(h3, h4);
    auto h5 = graph.create_handle("TGTTTGGAA");
    graph.create_edge(h4, h5);
    auto h6 = graph.create_handle("AATTACACCTTAAACTCAAGAGAATGAGGGTA");
    graph.create_edge(h5, h6);
    auto h7 = graph.create_handle("AAACATCAAATAATGTCTT");
    graph.create_edge(h6, h7);
    auto h8 = graph.create_handle("A");
    graph.create_edge(h7, h8);
    auto h9 = graph.create_handle("GTAGTATTTGGTGACTTAAATAGTTT");
    graph.create_edge(h8, h9);
    auto h10 = graph.create_handle("TGACCTATTGAACCCTCAAGAGTCCCTGGACC");
    graph.create_edge(h9, h10);
    auto h11 = graph.create_handle("A");
    graph.create_edge(h10, h11);
    auto h12 = graph.create_handle("TACTTTGTGAATTACTGTTTTA");
    graph.create_edge(h11, h12);
    auto h13 = graph.create_handle("G");
    graph.create_edge(h12, h13);
    auto h14 = graph.create_handle("TTGATGGACTA");
    graph.create_edge(h13, h14);
    auto h15 = graph.create_handle("ATAAATGCTATCAGATTCCTATTTTGTCAAGA");
    graph.create_edge(h14, h15);
    // Another branch off 12
    auto h16 = graph.create_handle("A");
    graph.create_edge(h12, h16);
    auto h17 = graph.create_handle("TTGATGGACTA");
    graph.create_edge(h16, h17);
    auto h18 = graph.create_handle("ATAAATGCTATCAGATTCCTATTTTGTCAAGA");
    graph.create_edge(h17, h18);
        
    string read = string("AAACATCAAATAATGTCTTAGTAGTATTTGGTGACTTAAATAGTTTTGACCTATTGAACCCTCAAGAG");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    // Make sure to allow big gaps.
    aligner.align_pinned(aln, graph, true, true, 60);
    
    // Make sure we got the right score.
    // Account for full length bonus, one gap open, 50 gap extends, and a whole read of matches.
    REQUIRE(aln.score() == read.size() + 10 - 6 - 50);
}

TEST_CASE("XdropAligner can align pinned left across an insertion with extra graph at the end", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    auto h1 = graph.create_handle("AATGAGAAAGGAAAAAGCTTTGGGAAAGTAGCTAAGCAGAGGATGCTTCTTAAAAAATGTCAAAAAAA");
        
    string read = string("CACATATGTCACTAGGAATGAGAAAGGAAAAAGCTTTGGGAAAGTAGCTAAGCAGAGGATGCTTCTTAAAAAATGTC");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align without xdrop
    aligner.align_pinned(aln, graph, true, false);
    auto no_xdrop_score = aln.score();
    
    // Align with xdrop
    aligner.align_pinned(aln, graph, true, true);
        
    // Score should be the same with and without xdrop.
    REQUIRE(aln.score() == no_xdrop_score);
    
    // Make sure we got the right score.
    // Account for full length bonus, one gap open, 15 gap extends, and the lack of matches for those.
    REQUIRE(aln.score() == read.size() + 10 - 6 - 15 - 16);
}

TEST_CASE("XdropAligner can align pinned left to a forking graph", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    Node* n0 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    Node* n1 = graph.create_node("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    Node* n2 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    Node* n3 = graph.create_node("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    Node* n4 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    Node* n5 = graph.create_node("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    Node* n6 = graph.create_node("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    Node* n7 = graph.create_node("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    Node* n8 = graph.create_node("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    Node* n9 = graph.create_node("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    
    graph.create_edge(n0, n1);
    graph.create_edge(n0, n3);
    graph.create_edge(n1, n2);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    graph.create_edge(n6, n7);
    graph.create_edge(n7, n8);
    graph.create_edge(n8, n9);
    
    string read = string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    aligner.align_pinned(aln, graph, true, true);
    
    // Make sure we got the right score.
    // Account for full length bonus.
    REQUIRE(aln.score() == read.size() + 10);
    
    // Make sure we take the right path
    REQUIRE(aln.path().mapping_size() == 3);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(0).position().offset() == 0);
    REQUIRE(aln.path().mapping(0).edit_size() == 1);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == 32);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == 32);
    REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
    REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
    REQUIRE(aln.path().mapping(1).position().offset() == 0);
    REQUIRE(aln.path().mapping(1).edit_size() == 1);
    REQUIRE(aln.path().mapping(1).edit(0).from_length() == 32);
    REQUIRE(aln.path().mapping(1).edit(0).to_length() == 32);
    REQUIRE(aln.path().mapping(1).edit(0).sequence() == "");
    REQUIRE(aln.path().mapping(2).position().node_id() == n2->id());
    REQUIRE(aln.path().mapping(2).position().offset() == 0);
    REQUIRE(aln.path().mapping(2).edit_size() == 1);
    REQUIRE(aln.path().mapping(2).edit(0).from_length() == 32);
    REQUIRE(aln.path().mapping(2).edit(0).to_length() == 32);
    REQUIRE(aln.path().mapping(2).edit(0).sequence() == "");
}

TEST_CASE("XdropAligner can align pinned right", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
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
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln, graph, false, true, false, max_gap_length);
    
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
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    Node* n0 = graph.create_node("TTAAGCTGAGGGAATAGTGCCTGGCATCGAGGAAAGCCTCTGA");
    
    string read = string("AGCTGAGGGAATAGTGCCTGGCATCGAGGAAAGCCTCTGA");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln, graph, true, true, false, max_gap_length);
    
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
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    Node* n0 = graph.create_node("AAAGAGGTCAATAGCCAAAT");
    
    string read = string("GAAAGAGGTCAATAGCCAAAT");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln, graph, true, true, false, max_gap_length);
    
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

TEST_CASE("XdropAligner can align pinned left with a leading deletion", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    Node* n0 = graph.create_node("GAAAGAGGTCAATAGCCAAAT");
    
    string read = string("AAAGAGGTCAATAGCCAAAT");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln, graph, true, true, false, max_gap_length);
    
    // Make sure we got the right score.
    // Account for full length bonus and one open
    REQUIRE(aln.score() == read.size() + 10 - 6);
    
    // Make sure we take the right path (leading 1 bp insertion)
    REQUIRE(aln.path().mapping_size() == 1);
    REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
    REQUIRE(aln.path().mapping(0).position().offset() == 0);
    REQUIRE(aln.path().mapping(0).edit_size() == 2);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
    REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
    REQUIRE(aln.path().mapping(0).edit(1).from_length() == read.size());
    REQUIRE(aln.path().mapping(0).edit(1).to_length() == read.size());
    REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
}

TEST_CASE("XdropAligner can align pinned right with a trailing insertion", "[xdrop][alignment][mapping]") {
    
    VG graph;
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    Node* n0 = graph.create_node("ATTTGGCTATTGACCTCTTT");
    
    string read = string("ATTTGGCTATTGACCTCTTTC");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned right, letting the graph compute a topological order
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln, graph, false, true, false, max_gap_length);
    
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
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    Node* n0 = graph.create_node("A");
    
    // Not even the full length bonus can save us
    string read = string("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    Alignment aln;
    aln.set_sequence(read);
    
    // Align pinned left, letting the graph compute a topological order
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln, graph, true, true, false, max_gap_length);
    
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

TEST_CASE("XdropAligner can select the best head and tail nodes automatically in pinned alignment", "[xdrop][alignment][mapping]") {
    
    bdsg::HashGraph graph;
    
    handle_t h1 = graph.create_handle("ATA");
    handle_t h2 = graph.create_handle("CGC");
    handle_t h3 = graph.create_handle("A");
    handle_t h4 = graph.create_handle("AGA");
    handle_t h5 = graph.create_handle("CTC");
    
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h3, h5);
    
    Alignment aln1;
    aln1.set_sequence("ATA");
    
    Alignment aln2;
    aln2.set_sequence("CGC");
    
    Alignment aln3;
    aln3.set_sequence("AGA");
    
    Alignment aln4;
    aln4.set_sequence("CTC");
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 5);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln1, graph, true, true, false, max_gap_length);
    aligner.align_pinned(aln2, graph, true, true, false, max_gap_length);
    aligner.align_pinned(aln3, graph, false, true, false, max_gap_length);
    aligner.align_pinned(aln4, graph, false, true, false, max_gap_length);
    
    REQUIRE(aln1.score() == 8);
    REQUIRE(aln2.score() == 8);
    REQUIRE(aln3.score() == 8);
    REQUIRE(aln4.score() == 8);
    
    REQUIRE(aln1.path().mapping_size() == 1);
    REQUIRE(aln2.path().mapping_size() == 1);
    REQUIRE(aln3.path().mapping_size() == 1);
    REQUIRE(aln4.path().mapping_size() == 1);
    
    REQUIRE(aln1.path().mapping(0).edit_size() == 1);
    REQUIRE(aln2.path().mapping(0).edit_size() == 1);
    REQUIRE(aln3.path().mapping(0).edit_size() == 1);
    REQUIRE(aln4.path().mapping(0).edit_size() == 1);
    
    REQUIRE(aln1.path().mapping(0).position().node_id() == graph.get_id(h1));
    REQUIRE(aln2.path().mapping(0).position().node_id() == graph.get_id(h2));
    REQUIRE(aln3.path().mapping(0).position().node_id() == graph.get_id(h4));
    REQUIRE(aln4.path().mapping(0).position().node_id() == graph.get_id(h5));
}

TEST_CASE("QualAdjXdropAligner can perform a quality-adjusted alignment without crashing", "[xdrop][alignment][mapping]") {

    bdsg::HashGraph graph;
    
    handle_t h1 = graph.create_handle("ATTTGGTACC");
    
    Alignment aln1;
    aln1.set_sequence("ATTTG");
    aln1.set_quality("HHHHH");
    alignment_quality_char_to_short(aln1);
    
    Alignment aln2;
    aln2.set_sequence("GTACC");
    aln2.set_quality("HHHHH");
    alignment_quality_char_to_short(aln2);
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 5);
    const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
    
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln1, graph, true, true, false, max_gap_length);
    aligner.align_pinned(aln2, graph, false, true, false, max_gap_length);
        
    REQUIRE(aln1.score() == 5 * 1 + 5);
    REQUIRE(aln1.path().mapping_size() == 1);
    REQUIRE(aln1.path().mapping(0).position().node_id() == graph.get_id(h1));
    REQUIRE(aln1.path().mapping(0).position().offset() == 0);
    REQUIRE(aln1.path().mapping(0).edit_size() == 1);
    REQUIRE(aln1.path().mapping(0).edit(0).from_length() == 5);
    REQUIRE(aln1.path().mapping(0).edit(0).to_length() == 5);
    REQUIRE(aln1.path().mapping(0).edit(0).sequence() == "");
    
    REQUIRE(aln2.score() == 5 * 1 + 5);
    REQUIRE(aln2.path().mapping_size() == 1);
    REQUIRE(aln2.path().mapping(0).position().node_id() == graph.get_id(h1));
    REQUIRE(aln2.path().mapping(0).position().offset() == 5);
    REQUIRE(aln2.path().mapping(0).edit_size() == 1);
    REQUIRE(aln2.path().mapping(0).edit(0).from_length() == 5);
    REQUIRE(aln2.path().mapping(0).edit(0).to_length() == 5);
    REQUIRE(aln2.path().mapping(0).edit(0).sequence() == "");
}



TEST_CASE("QualAdjXdropAligner will not penalize a low quality mismatch", "[xdrop][alignment][mapping]") {
    
    bdsg::HashGraph graph;
    
    handle_t h1 = graph.create_handle("ATTTGGTACC");
    
    Alignment aln;
    aln.set_sequence("ATCTG");
    aln.set_quality("HH#HH");
    alignment_quality_char_to_short(aln);
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 5);
    const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
    
    uint16_t max_gap_length = 40;
    aligner.align_pinned(aln, graph, true, true, false, max_gap_length);
    
    REQUIRE(aln.score() == 4 * 1 + 5);
    REQUIRE(aln.path().mapping_size() == 1);
    REQUIRE(aln.path().mapping(0).edit_size() == 3);
    REQUIRE(aln.path().mapping(0).edit(0).from_length() == 2);
    REQUIRE(aln.path().mapping(0).edit(0).to_length() == 2);
    REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
    REQUIRE(aln.path().mapping(0).edit(1).from_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(1).to_length() == 1);
    REQUIRE(aln.path().mapping(0).edit(1).sequence() == "C");
    REQUIRE(aln.path().mapping(0).edit(2).from_length() == 2);
    REQUIRE(aln.path().mapping(0).edit(2).to_length() == 2);
    REQUIRE(aln.path().mapping(0).edit(2).sequence() == "");
}

TEST_CASE("XdropAligner doesn't crash on a case where it is hard to find a seed", "[xdrop][alignment][mapping]") {
    
    string graph_json = R"({"edge": [{"from": "92345167", "to": "92345168"}, {"from": "92345182", "to": "92345183"}, {"from": "92345165", "to": "92345166"}, {"from": "92345177", "to": "92345178"}, {"from": "92345171", "to": "92345172"}, {"from": "92345161", "to": "92345162"}, {"from": "92345183", "to": "92345184"}, {"from": "92345181", "to": "92345182"}, {"from": "92345178", "to": "92345179"}, {"from": "92345166", "to": "92345167"}, {"from": "92345179", "to": "92345180"}, {"from": "92345173", "to": "92345174"}, {"from": "92345184", "to": "92345185"}, {"from": "92345169", "to": "92345170"}, {"from": "92345185", "to": "92345186"}, {"from": "92345160", "to": "92345161"}, {"from": "92345174", "to": "92345175"}, {"from": "92345162", "to": "92345163"}, {"from": "92345175", "to": "92345176"}, {"from": "92345168", "to": "92345169"}, {"from": "92345163", "to": "92345164"}, {"from": "92345172", "to": "92345173"}, {"from": "92345180", "to": "92345181"}, {"from": "92345176", "to": "92345177"}, {"from": "92345170", "to": "92345171"}, {"from": "92345164", "to": "92345165"}], "node": [{"id": "92345167", "sequence": "TTTATATATATATATTTATATATATATATTTA"}, {"id": "92345182", "sequence": "TATATATATTTATATATATATTTATATATATA"}, {"id": "92345165", "sequence": "ATATATATATATTTATATATATTTATATATTA"}, {"id": "92345177", "sequence": "TTTATATATATATTTATATATATATATTATAT"}, {"id": "92345171", "sequence": "TTATATATATATTTATATATATATTTATATAT"}, {"id": "92345161", "sequence": "ATATATTTATATATTTTTATATATTATATATT"}, {"id": "92345183", "sequence": "TTTATATATATTTATATATATATTTATATATA"}, {"id": "92345181", "sequence": "ATATATTATATATATATTTATATATATATTTA"}, {"id": "92345178", "sequence": "ATATATTTATATATATATTTATATATATATTT"}, {"id": "92345166", "sequence": "TTTATATATATTTATATATATATTTATATATA"}, {"id": "92345179", "sequence": "ATATATATATTTATATATATATTTATATATAT"}, {"id": "92345173", "sequence": "ATATTTATATATATATATTTATATATATATTT"}, {"id": "92345184", "sequence": "TATTTATATATATATTTATATATATTTATATA"}, {"id": "92345169", "sequence": "TTTATATATATATTTATATATATATTTATATA"}, {"id": "92345185", "sequence": "TATATTTATATATATATATATATATTTATATA"}, {"id": "92345160", "sequence": "ATTTATATATATATTTATATATATATTTATAT"}, {"id": "92345174", "sequence": "ATATATATATTTATATATATATTATTTATATA"}, {"id": "92345162", "sequence": "TATATATATATTTATATATTATATATATATTT"}, {"id": "92345175", "sequence": "TATATTTATATATATATTATATATATATTTAT"}, {"id": "92345168", "sequence": "TATATATATTTATATATATATTTATATATATA"}, {"id": "92345163", "sequence": "ATATATTTATATATATATTTATATATATTTAT"}, {"id": "92345172", "sequence": "ATATATATATATTTATATATATATTTATATAT"}, {"id": "92345180", "sequence": "ATTTATATATATATTTATATATATATTTATAT"}, {"id": "92345176", "sequence": "ATATATATATTATATATATATTTATATATATA"}, {"id": "92345170", "sequence": "TATATTTATATATATATATTATATATATATAT"}, {"id": "92345164", "sequence": "ATATATATTTATATATATTTATATATATATTT"}, {"id": "92345186", "sequence": "TATATTTATATATATTTATATATATATTTATA"}]})";
    
    Graph source;
    json2pb(source, graph_json.c_str(), graph_json.size());
    
    VG graph;
    graph.extend(source);
    
    Alignment aln;
    aln.set_sequence("CAGCACTTTGGGAGGCCAAGGTGGGTGGATCATCTGAGGTCAGGAGTTTGAGACCAGCCTGACCAACATGGTGAAATCCTGTCTCTACTGAAAATACTAAAATTAGCCAGGCGTGGCGGCCAGTGCCTGTAATCCCGGCTACTGGGGAGG");
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 10);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    aligner.align_xdrop(aln, graph, vector<MaximalExactMatch>(), false);
}

TEST_CASE("XdropAligner pinned alignment doesn't crash when aligning to a long stretch of mismatches",
          "[xdrop][alignment][mapping][pinned]") {
    
    bdsg::HashGraph graph;
    
    handle_t h0 = graph.create_handle("CCACCATCTTGTTCACTCTGGGGCCACAGACT");
    handle_t h1 = graph.create_handle("GTCTTTTTCTGGTCTCTGCTTCCCTGCTTCAT");
    handle_t h2 = graph.create_handle("CCTCCTTCTACTCTCTGCTTCCCTAGCGTGTG");
    handle_t h3 = graph.create_handle("GCCCAGATGGTCAGTCACAATCCTGACTCCAC");
    handle_t h4 = graph.create_handle("AGCAGTTTTGGGGTCAAGCCTGTAGACAGGAG");
    handle_t h5 = graph.create_handle("TTACTTATCATCTTTGAGTTTATTTAATTTTT");
    handle_t h6 = graph.create_handle("CAATGGGAGAACTAGATTGTCCAGTCTTGGCC");
    handle_t h7 = graph.create_handle("AAAAAAATGGTCTAGCTTTGAGTCATACTGTA");
    handle_t h8 = graph.create_handle("ATCATCTGTGGCTCAAAGGCAAGATCCTGCCC");
    handle_t h9 = graph.create_handle("ACTGTCCACTCGGCAGGGCTGTGGTGGGCACC");
    handle_t h10 = graph.create_handle("ACAAGGAGGAGTATTTCTTCTTCA");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    graph.create_edge(h5, h6);
    graph.create_edge(h6, h7);
    graph.create_edge(h7, h8);
    graph.create_edge(h8, h9);
    graph.create_edge(h9, h10);
    
    Alignment aln;
    aln.set_sequence("CAGATCCCTCGACCATCCGGTCAGGATACACAAAAGGACAGCAAAGGGGTTGAGAAGGGCTGAGGGGAGAAAAGCCAGGAAGCTGAGATCAGCAGAGGCCAAGCATAAAAACTGGGAGGATGCTACGAAGCTGCAGATGACAGCATCATTTTCTTGAAGAACATTCAAGGATTTGTCATAGTGGCTGGGCTTTCACTGATTGATTGAAGTCTACAAACAGCACTTCAATTGGTATCGGTCAAGTTCTTTAAGATTTAGGAAATTGATTGGAGCGGAAAATTGTAAGTTACAAAATTCGCACTGAAGTCCCATTAAAACCAC");
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 1, 1, 1, 0);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    aligner.align_pinned(aln, graph, false, true, 10);
}

TEST_CASE("XdropAligner pinned alignment doesn't crash when the optimal alignment column would be x-dropped if not for the full length bonus",
          "[xdrop][alignment][mapping][pinned]") {
    
    bdsg::HashGraph graph;
    
    handle_t h0 = graph.create_handle("AAGGG");
    
    
    Alignment aln;
    aln.set_sequence("AACGT");
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 9);
    const Aligner& aligner = *aligner_source.get_regular_aligner();
    
    aligner.align_pinned(aln, graph, true, true, 0);
    
    REQUIRE(aln.score() == 1 + 1 + 1 - 4 - 4 + 9);
}



TEST_CASE("QualAdjXdropAligner uses quality adjusted full length bonuses",
          "[xdrop][alignment][mapping][pinned]") {

    bdsg::HashGraph graph;
    
    handle_t h0 = graph.create_handle("AAGGG");
    
    Alignment aln1;
    aln1.set_sequence("AA");
    aln1.set_quality("H#");
    alignment_quality_char_to_short(aln1);
    
    Alignment aln2;
    aln2.set_sequence("GG");
    aln2.set_quality("#H");
    alignment_quality_char_to_short(aln2);
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 5);
    const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
    
    aligner.align_pinned(aln1, graph, true, true, 0);
    aligner.align_pinned(aln2, graph, false, true, 0);
        
    REQUIRE(aln1.score() == 1);
    REQUIRE(aln2.score() == 1);
}

TEST_CASE("QualAdjXdropAligner doesn't crash when a traceback goes through an X-dropped vector",
          "[xdrop][alignment][mapping][pinned]") {
    
    bdsg::HashGraph graph;
    
    handle_t h0 = graph.create_handle("ACCTACCAAATCACT");
    
    Alignment aln;
    aln.set_sequence("TCCTACCTAATCA");
    aln.set_quality(":F:F,FFF:::F,");
    alignment_quality_char_to_short(aln);
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 20);
    const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
    
    aligner.align_pinned(aln, graph, true, true, 0);
    
    REQUIRE(path_to_length(aln.path()) == aln.sequence().size());
}

TEST_CASE("XdropAligner can align across multiple nodes",
          "[xdrop][alignment][mapping][pinned]") {
    
    bdsg::HashGraph graph;
    
    handle_t h0 = graph.create_handle("TC");
    handle_t h1 = graph.create_handle("ATCATTATGTTTTTGAGTATAGTATAGTTACT");
    handle_t h2 = graph.create_handle("TGTTTTCACTGCTTTATATACCATTGTATAAA");
    handle_t h3 = graph.create_handle("TATACCACAATTTATCCATTCTAGTATTAATG");
    handle_t h4 = graph.create_handle("GACATTAGGTTTGTTGTTATTACAGTGTTCCT");
    handle_t h5 = graph.create_handle("GCAAGTAACCTTGCACATTTCTCCTGATACAC");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    
    Alignment aln;
    aln.set_sequence("CTTTAGGTAATTTTAAACAAACAAAAGCTATCTCAAAATTTTTTCACACTTTTCAAACCCCTGATCCCTCATTTACTCTTTGGAGATGCCCCACCTTATCTCTTTCCAAAAACTAAGGACAACTANNNNNNNNNCNNNNTNG");
    aln.set_quality("<0FI000F<I0FII7I77B0B00I077007<077<<I07FI0IB0FFB00FF0F<I00FIF77F<0II0<000<000<F<0B<F000000FB0BF<B7F0BBBF077000I0F00000000F000'''''''''0''''0'<");
    alignment_quality_char_to_short(aln);
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 5);
    const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
    
    aligner.align_pinned(aln, graph, true, true, 9);
    
    REQUIRE(path_to_length(aln.path()) == aln.sequence().size());
}

TEST_CASE("XdropAligner can align to a graph where some nodes are entirely unaligned",
          "[xdrop][alignment][mapping][pinned]") {
    
    bdsg::HashGraph graph;
    
    handle_t h0 = graph.create_handle("CTTCAGGAGGCAACATCAAGGTAACCCGAGTT");
    handle_t h1 = graph.create_handle("TAGGGATATTCTGCAAAATAACTGGCCTGTAA");
    handle_t h2 = graph.create_handle("ACCTCAAAATTAAGTCAAGATCATGGAACCAA");
    handle_t h3 = graph.create_handle("GAAAACACTGAGGAATTGTTCCACACTAAAGA");
    handle_t h4 = graph.create_handle("AGTGTGCAGAAACATGGCAACTCCATGTGGTT");
    handle_t h5 = graph.create_handle("CATGATTCTGAATCGGCA");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    
    Alignment aln;
    aln.set_sequence("CGGGGCCTTCCTCCAGTGAATATCTCTTCTAACATTGTCAATCATATCCAACAGAAGAAGAGATCGCTCCATTGGTACAAAAAGCTATAGCGCTGAATAATAGTTGGACGCTGACGGGCTATCTGAGTATGACAGAGCGCN");
    aln.set_quality("B7FFFFFIFIBFFFIII0FIIIIIF<IIBFIIIIIFF7F0BI7IFFF07FIIF00<0B<FB7BF00B007FF<0BI0BFF000000<000<0B07F00F<0B0B0F000007000000000000000000<000000F00'");
    alignment_quality_char_to_short(aln);
    
    TestAligner aligner_source;
    aligner_source.set_alignment_scores(1, 4, 6, 1, 5);
    const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
    
    aligner.align_pinned(aln, graph, false, true, 10);
    
    REQUIRE(path_to_length(aln.path()) == aln.sequence().size());
}

TEST_CASE("X-drop alignment applies quality penalities to the read sequence and not the ref", "[aligner][xdrop]") {
    
    // this case turned up as an unexpected instance of GSSW and dozeu producing different alignments
    
    bdsg::HashGraph graph;
    auto h1 = graph.create_handle("TGGTGCAGGAAAAGAAGATTGGA");
    auto h2 = graph.create_handle("AAAAGACCAAAACACACTTCTCCAGCAACACT");
    auto h3 = graph.create_handle("CCATCAGCTTTTTAAAATTTAGAGCTATCTGC");
    auto h4 = graph.create_handle("TAATTTTTTCCCTCTTCCTTCTCAATAAATGA");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    
    TestAligner aligner_source(0.417098);
    aligner_source.set_alignment_scores(1, 4, 6, 1, 5);
    const auto& aligner = *(aligner_source.get_qual_adj_aligner());
    
    Alignment aln;
    aln.set_sequence("CGGCGAGTGAAGGTGATTACCGGTTAGCAAAAAAACTAGCTTCTCAAGGACACCGCGAACTGCTCTGATTTGGAT");
    aln.set_quality("+,C########################################################################");
    alignment_quality_char_to_short(aln);
    
    aligner.align_pinned(aln, graph, true, true, 26);
    
    // should find a short alignment of 3 bases
    REQUIRE(aln.score() == 1);
    REQUIRE(path_from_length(aln.path()) == 3);
}


}
}
        
