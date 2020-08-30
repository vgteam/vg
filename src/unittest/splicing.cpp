/// \file dinucleotide_machine.cpp
///  
/// Unit tests for the DinucleotideMachine
///

#include <iostream>
#include <random>

#include "../splicing.hpp"
#include "catch.hpp"
#include "test_aligner.hpp"

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
    
    REQUIRE(m0.front().first == pos_t(graph.get_id(h), false, 4));
    REQUIRE(m0.front().second == 1);
    REQUIRE(m1.front().first == pos_t(graph.get_id(h), false, 5));
    REQUIRE(m1.front().second == 2);
    
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
    
    REQUIRE(m0.front().first == pos_t(graph.get_id(h), false, 6));
    REQUIRE(m0.front().second == 0);
    REQUIRE(m2.front().first == pos_t(graph.get_id(h), false, 3));
    REQUIRE(m2.front().second == 3);
    
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
    
    REQUIRE(m0.front().first == pos_t(graph.get_id(h0), false, 2));
    REQUIRE(m0.front().second == 2);
    REQUIRE(m1.front().first == pos_t(graph.get_id(h0), false, 2));
    REQUIRE(m1.front().second == 2);
    REQUIRE(m1.back().first == pos_t(graph.get_id(h1), false, 0));
    REQUIRE(m1.back().second == 3);
    REQUIRE(m2.front().first == pos_t(graph.get_id(h3), false, 0));
    REQUIRE(m2.front().second == 4);
    REQUIRE(m3.front().first == pos_t(graph.get_id(h2), false, 0));
    REQUIRE(m3.front().second == 3);
    
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
    
    REQUIRE(m0.front().first == pos_t(graph.get_id(h2), false, 1));
    REQUIRE(m0.front().second == 2);
    REQUIRE(m1.front().first == pos_t(graph.get_id(h1), false, 1));
    REQUIRE(m1.front().second == 2);
    
}

TEST_CASE("Softclip trimming works on a simple example",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h0 = graph.create_handle("TACCGATAGAC");
    
    Alignment aln;
    aln.set_sequence("ACCGATAGA");
    auto path = aln.mutable_path();
    auto mapping = path->add_mapping();
    auto position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h0));
    position->set_is_reverse(false);
    position->set_offset(1);
    auto edit = mapping->add_edit();
    edit->set_from_length(9);
    edit->set_to_length(9);
    aln.set_score(19);
    
    TestAligner test_aligner;
    test_aligner.set_alignment_scores(1, 4, 6, 1, 5);
    auto aligner = test_aligner.get_regular_aligner();
    
    auto trimmed_left = trimmed_end(aln, 2, false, graph, *aligner);
    auto trimmed_right = trimmed_end(aln, 2, true, graph, *aligner);
    
    REQUIRE(get<0>(trimmed_left) == pos_t(graph.get_id(h0), false, 3));
    REQUIRE(get<1>(trimmed_left) == 2);
    REQUIRE(get<2>(trimmed_left) == 7);
    REQUIRE(get<0>(trimmed_right) == pos_t(graph.get_id(h0), false, 8));
    REQUIRE(get<1>(trimmed_right) == 2);
    REQUIRE(get<2>(trimmed_right) == 7);
}

TEST_CASE("Softclip trimming can search for a trim point across multiple nodes",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h0 = graph.create_handle("GAC");
    handle_t h1 = graph.create_handle("T");
    handle_t h2 = graph.create_handle("TGG");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);
    
    Alignment aln;
    aln.set_sequence("ACTTG");
    auto path = aln.mutable_path();
    
    auto mapping = path->add_mapping();
    auto position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h0));
    position->set_is_reverse(false);
    position->set_offset(1);
    auto edit = mapping->add_edit();
    edit->set_from_length(2);
    edit->set_to_length(2);
    
    mapping = path->add_mapping();
    position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h1));
    position->set_is_reverse(false);
    position->set_offset(0);
    edit = mapping->add_edit();
    edit->set_from_length(1);
    edit->set_to_length(1);
    
    mapping = path->add_mapping();
    position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h2));
    position->set_is_reverse(false);
    position->set_offset(0);
    edit = mapping->add_edit();
    edit->set_from_length(2);
    edit->set_to_length(2);
    
    
    aln.set_score(15);
    
    TestAligner test_aligner;
    test_aligner.set_alignment_scores(1, 4, 6, 1, 5);
    auto aligner = test_aligner.get_regular_aligner();
    
    auto trimmed_left1 = trimmed_end(aln, 3, false, graph, *aligner);
    auto trimmed_left2 = trimmed_end(aln, 4, false, graph, *aligner);
    auto trimmed_right1 = trimmed_end(aln, 3, true, graph, *aligner);
    auto trimmed_right2 = trimmed_end(aln, 4, true, graph, *aligner);
    
    REQUIRE(get<0>(trimmed_left1) == pos_t(graph.get_id(h1), false, 1));
    REQUIRE(get<1>(trimmed_left1) == 3);
    REQUIRE(get<2>(trimmed_left1) == 8);
    
    REQUIRE(get<0>(trimmed_left2) == pos_t(graph.get_id(h2), false, 1));
    REQUIRE(get<1>(trimmed_left2) == 4);
    REQUIRE(get<2>(trimmed_left2) == 9);

    REQUIRE(get<0>(trimmed_right1) == pos_t(graph.get_id(h1), false, 0));
    REQUIRE(get<1>(trimmed_right1) == 3);
    REQUIRE(get<2>(trimmed_right1) == 8);

    REQUIRE(get<0>(trimmed_right2) == pos_t(graph.get_id(h0), false, 2));
    REQUIRE(get<1>(trimmed_right2) == 4);
    REQUIRE(get<2>(trimmed_right2) == 9);
}

TEST_CASE("Softclip trimming skips over softclips",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h0 = graph.create_handle("TACCGATAGAC");
    
    Alignment aln;
    aln.set_sequence("GGGGGGGACCGATAGATTTTTTT");
    auto path = aln.mutable_path();
    auto mapping = path->add_mapping();
    auto position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h0));
    position->set_is_reverse(false);
    position->set_offset(1);
    
    auto edit = mapping->add_edit();
    edit->set_from_length(0);
    edit->set_to_length(7);
    edit->set_sequence("GGGGGGG");
    
    edit = mapping->add_edit();
    edit->set_from_length(9);
    edit->set_to_length(9);
    
    edit = mapping->add_edit();
    edit->set_from_length(0);
    edit->set_to_length(7);
    edit->set_sequence("TTTTTTT");
    
    aln.set_score(9);
    
    TestAligner test_aligner;
    test_aligner.set_alignment_scores(1, 4, 6, 1, 5);
    auto aligner = test_aligner.get_regular_aligner();
    
    auto trimmed_left = trimmed_end(aln, 2, false, graph, *aligner);
    auto trimmed_right = trimmed_end(aln, 2, true, graph, *aligner);
    
    REQUIRE(get<0>(trimmed_left) == pos_t(graph.get_id(h0), false, 3));
    REQUIRE(get<1>(trimmed_left) == 9);
    REQUIRE(get<2>(trimmed_left) == 2);
    REQUIRE(get<0>(trimmed_right) == pos_t(graph.get_id(h0), false, 8));
    REQUIRE(get<1>(trimmed_right) == 9);
    REQUIRE(get<2>(trimmed_right) == 2);
}


TEST_CASE("Softclip trimming will go extra distance to avoid splitting an indel",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h0 = graph.create_handle("TACCGATAGAC");
    
    Alignment aln;
    aln.set_sequence("CGAGGTAG");
    auto path = aln.mutable_path();
    auto mapping = path->add_mapping();
    auto position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h0));
    position->set_is_reverse(false);
    position->set_offset(3);
    
    auto edit = mapping->add_edit();
    edit->set_from_length(3);
    edit->set_to_length(3);
    
    edit = mapping->add_edit();
    edit->set_from_length(0);
    edit->set_to_length(2);
    edit->set_sequence("GG");
    
    edit = mapping->add_edit();
    edit->set_from_length(3);
    edit->set_to_length(3);
    
    aln.set_score(9);
    
    TestAligner test_aligner;
    test_aligner.set_alignment_scores(1, 4, 6, 1, 5);
    auto aligner = test_aligner.get_regular_aligner();
    
    auto trimmed_left = trimmed_end(aln, 4, false, graph, *aligner);
    auto trimmed_right = trimmed_end(aln, 4, true, graph, *aligner);
    
    REQUIRE(get<0>(trimmed_left) == pos_t(graph.get_id(h0), false, 6));
    REQUIRE(get<1>(trimmed_left) == 5);
    REQUIRE(get<2>(trimmed_left) == 1);
    REQUIRE(get<0>(trimmed_right) == pos_t(graph.get_id(h0), false, 6));
    REQUIRE(get<1>(trimmed_right) == 5);
    REQUIRE(get<2>(trimmed_right) == 1);
}

TEST_CASE("Softclip trimming will trim the entire path if necessary",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h0 = graph.create_handle("GAC");
    handle_t h1 = graph.create_handle("T");
    handle_t h2 = graph.create_handle("TGG");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);
    
    Alignment aln;
    aln.set_sequence("AG");
    auto path = aln.mutable_path();
    
    auto mapping = path->add_mapping();
    auto position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h0));
    position->set_is_reverse(false);
    position->set_offset(1);
    auto edit = mapping->add_edit();
    edit->set_from_length(1);
    edit->set_to_length(1);
    
    edit = mapping->add_edit();
    edit->set_from_length(1);
    edit->set_to_length(0);
    
    mapping = path->add_mapping();
    position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h1));
    position->set_is_reverse(false);
    position->set_offset(0);
    edit = mapping->add_edit();
    edit->set_from_length(1);
    edit->set_to_length(0);
    
    mapping = path->add_mapping();
    position = mapping->mutable_position();
    position->set_node_id(graph.get_id(h2));
    position->set_is_reverse(false);
    position->set_offset(0);
    edit = mapping->add_edit();
    edit->set_from_length(1);
    edit->set_to_length(0);
    
    edit = mapping->add_edit();
    edit->set_from_length(1);
    edit->set_to_length(1);
    
    aln.set_score(4);
    
    TestAligner test_aligner;
    test_aligner.set_alignment_scores(1, 4, 6, 1, 5);
    auto aligner = test_aligner.get_regular_aligner();
    
    auto trimmed_left = trimmed_end(aln, 3, false, graph, *aligner);
    auto trimmed_right = trimmed_end(aln, 3, true, graph, *aligner);
    
    REQUIRE(get<0>(trimmed_left) == pos_t(graph.get_id(h2), false, 2));
    REQUIRE(get<1>(trimmed_left) == 2);
    REQUIRE(get<2>(trimmed_left) == 4);
    
    REQUIRE(get<0>(trimmed_right) == pos_t(graph.get_id(h0), false, 1));
    REQUIRE(get<1>(trimmed_right) == 2);
    REQUIRE(get<2>(trimmed_right) == 4);
    
}

}
}
        
