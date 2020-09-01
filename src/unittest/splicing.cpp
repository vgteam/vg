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

    TestAligner test_aligner;
    vector<tuple<string, string, double>> table;
    for (auto motif : motifs) {
        table.emplace_back(motif, motif, 1.0 / motifs.size());
    }
    SpliceMotifs splice_motifs(table, *test_aligner.get_regular_aligner());

    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, splice_motifs);

    auto m0 = splice_region.candidate_splice_sites(0);
    auto m1 = splice_region.candidate_splice_sites(1);
    auto m2 = splice_region.candidate_splice_sites(2);

    REQUIRE(m0.size() == 1);
    REQUIRE(m1.size() == 1);
    REQUIRE(m2.size() == 0);

    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m0.front())) == h);
    REQUIRE(get<1>(m0.front()) == 4);
    REQUIRE(get<2>(m0.front()) == 1);
    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m1.front())) == h);
    REQUIRE(get<1>(m1.front()) == 5);
    REQUIRE(get<2>(m1.front()) == 2);

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

    TestAligner test_aligner;
    vector<tuple<string, string, double>> table;
    for (auto motif : motifs) {
        table.emplace_back(motif, motif, 1.0 / motifs.size());
    }
    SpliceMotifs splice_motifs(table, *test_aligner.get_regular_aligner());

    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, splice_motifs);

    auto m0 = splice_region.candidate_splice_sites(0);
    auto m1 = splice_region.candidate_splice_sites(1);
    auto m2 = splice_region.candidate_splice_sites(2);

    REQUIRE(m0.size() == 1);
    REQUIRE(m1.size() == 0);
    REQUIRE(m2.size() == 1);

    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m0.front())) == h);
    REQUIRE(get<1>(m0.front()) == 6);
    REQUIRE(get<2>(m0.front()) == 0);
    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m2.front())) == h);
    REQUIRE(get<1>(m2.front()) == 3);
    REQUIRE(get<2>(m2.front()) == 3);

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

    TestAligner test_aligner;
    vector<tuple<string, string, double>> table;
    for (auto motif : motifs) {
        table.emplace_back(motif, motif, 1.0 / motifs.size());
    }
    SpliceMotifs splice_motifs(table, *test_aligner.get_regular_aligner());

    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, splice_motifs);

    auto m0 = splice_region.candidate_splice_sites(0);
    auto m1 = splice_region.candidate_splice_sites(1);
    auto m2 = splice_region.candidate_splice_sites(2);
    auto m3 = splice_region.candidate_splice_sites(3);

    REQUIRE(m0.size() == 1);
    REQUIRE(m1.size() == 2);
    REQUIRE(m2.size() == 1);
    REQUIRE(m3.size() == 1);

    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m0.front())) == h0);
    REQUIRE(get<1>(m0.front()) == 2);
    REQUIRE(get<2>(m0.front()) == 2);
    
    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m1.front())) == h0);
    REQUIRE(get<1>(m1.front()) == 2);
    REQUIRE(get<2>(m1.front()) == 2);
    
    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m1.back())) == h1);
    REQUIRE(get<1>(m1.back()) == 0);
    REQUIRE(get<2>(m1.back()) == 3);
    
    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m2.front())) == h3);
    REQUIRE(get<1>(m2.front()) == 0);
    REQUIRE(get<2>(m2.front()) == 4);
    
    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m3.front())) == h2);
    REQUIRE(get<1>(m3.front()) == 0);
    REQUIRE(get<2>(m3.front()) == 3);

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

    TestAligner test_aligner;
    vector<tuple<string, string, double>> table;
    for (auto motif : motifs) {
        table.emplace_back(motif, motif, 1.0 / motifs.size());
    }
    SpliceMotifs splice_motifs(table, *test_aligner.get_regular_aligner());

    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, splice_motifs);

    auto m0 = splice_region.candidate_splice_sites(0);
    auto m1 = splice_region.candidate_splice_sites(1);

    REQUIRE(m0.size() == 1);
    REQUIRE(m1.size() == 1);

    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m0.front())) == h2);
    REQUIRE(get<1>(m0.front()) == 1);
    REQUIRE(get<2>(m0.front()) == 2);
    
    REQUIRE(splice_region.get_subgraph().get_underlying_handle(get<0>(m1.front())) == h1);
    REQUIRE(get<1>(m1.front()) == 1);
    REQUIRE(get<2>(m1.front()) == 2);

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

TEST_CASE("JoinedSpliceGraph can correctly perform queries on a simple graph",
          "[splice]") {

    HashGraph graph;

    handle_t h0 = graph.create_handle("GAGACCC");
    handle_t h1 = graph.create_handle("TATACAT");
    handle_t h2 = graph.create_handle("TGGCCGG");

    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);

    pos_t left_pos(graph.get_id(h0), false, 1);
    pos_t right_pos(graph.get_id(h2), false, 6);
    int64_t search_dist = 4;

    IncrementalSubgraph left_subgraph(graph, left_pos, false, search_dist);
    IncrementalSubgraph right_subgraph(graph, right_pos, true, search_dist);

    handle_t left_splice_node = left_subgraph.handle_at_order(0);
    handle_t right_splice_node = right_subgraph.handle_at_order(0);

    JoinedSpliceGraph joined_graph(graph, left_subgraph, left_splice_node, 4,
                                   right_subgraph, right_splice_node, 3);

    // check generic invariants
    REQUIRE(joined_graph.has_node(joined_graph.get_id(joined_graph.left_seed_node())));
    REQUIRE(joined_graph.has_node(joined_graph.get_id(joined_graph.right_seed_node())));
    REQUIRE(joined_graph.left_seed_node() != joined_graph.right_seed_node());
    REQUIRE(joined_graph.get_sequence(joined_graph.flip(joined_graph.left_seed_node()))
            == reverse_complement(joined_graph.get_sequence(joined_graph.left_seed_node())));
    REQUIRE(joined_graph.get_sequence(joined_graph.flip(joined_graph.right_seed_node()))
            == reverse_complement(joined_graph.get_sequence(joined_graph.right_seed_node())));
    REQUIRE(!joined_graph.get_is_reverse(joined_graph.left_seed_node()));
    REQUIRE(!joined_graph.get_is_reverse(joined_graph.right_seed_node()));
    REQUIRE(!joined_graph.get_is_reverse(joined_graph.left_splice_node()));
    REQUIRE(!joined_graph.get_is_reverse(joined_graph.right_seed_node()));
    REQUIRE(joined_graph.get_is_reverse(joined_graph.flip(joined_graph.left_seed_node())));
    REQUIRE(joined_graph.get_is_reverse(joined_graph.flip(joined_graph.right_seed_node())));
    REQUIRE(joined_graph.get_sequence(joined_graph.left_seed_node()).size()
            == joined_graph.get_length(joined_graph.left_seed_node()));
    REQUIRE(joined_graph.get_sequence(joined_graph.right_seed_node()).size()
            == joined_graph.get_length(joined_graph.right_seed_node()));
    joined_graph.for_each_handle([&](const handle_t& handle) {
        for (handle_t h : {handle, joined_graph.flip(handle)}) {
            auto seq = joined_graph.get_sequence(h);
            for (size_t i = 0; i < seq.size(); ++i) {
                REQUIRE(joined_graph.get_base(h, i) == seq[i]);
                for (size_t j = 0; j < seq.size() + 1; ++j) {
                    REQUIRE(joined_graph.get_subsequence(h, i, j) == seq.substr(i, j));
                }
            }
            for (bool go_left : {false, true}) {
                size_t deg = 0;
                joined_graph.follow_edges(h, go_left, [&](const handle_t& n) {
                    ++deg;
                });
                REQUIRE(deg == joined_graph.get_degree(h, go_left));
            }
        }
    });

    // check for the specific topology we should have
    REQUIRE(joined_graph.get_node_count() == 2);
    bool found1 = false, found2 = false;
    size_t count = 0;
    joined_graph.for_each_handle([&](const handle_t& handle) {
        if (handle == joined_graph.left_seed_node()) {
            found1 = true;
            REQUIRE(joined_graph.get_sequence(handle) == "AGA");
        }
        else if (handle == joined_graph.right_seed_node()) {
            found2 = true;
            REQUIRE(joined_graph.get_sequence(handle) == "CCG");
        }
        ++count;
    });
    REQUIRE(found1);
    REQUIRE(found2);
    REQUIRE(count == 2);
    count = 0;
    found1 = found2 = false;

    joined_graph.follow_edges(joined_graph.left_seed_node(), true,
                              [&](const handle_t& next) {
        REQUIRE(false);
    });
    joined_graph.follow_edges(joined_graph.left_seed_node(), false,
                              [&](const handle_t& next) {
        REQUIRE(next == joined_graph.right_seed_node());
        count++;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(joined_graph.right_seed_node(), true,
                              [&](const handle_t& next) {
        REQUIRE(next == joined_graph.left_seed_node());
        count++;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(joined_graph.right_seed_node(), false,
                              [&](const handle_t& next) {
        REQUIRE(false);
    });

    joined_graph.follow_edges(joined_graph.flip(joined_graph.left_seed_node()), true,
                              [&](const handle_t& next) {
        REQUIRE(next == joined_graph.flip(joined_graph.right_seed_node()));
        count++;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(joined_graph.flip(joined_graph.left_seed_node()), false,
                              [&](const handle_t& next) {
        REQUIRE(false);
    });
    joined_graph.follow_edges(joined_graph.flip(joined_graph.right_seed_node()), true,
                              [&](const handle_t& next) {
        REQUIRE(false);
    });
    joined_graph.follow_edges(joined_graph.flip(joined_graph.right_seed_node()), false,
                              [&](const handle_t& next) {
        REQUIRE(next == joined_graph.flip(joined_graph.left_seed_node()));
        count++;
    });
    REQUIRE(count == 1);
    count = 0;
    
    {
        
        Path path;
        Position* p0 = path.add_mapping()->mutable_position();
        p0->set_node_id(joined_graph.get_id(joined_graph.left_seed_node()));
        p0->set_is_reverse(false);
        p0->set_offset(1);
        
        Position* p1 = path.add_mapping()->mutable_position();
        p1->set_node_id(joined_graph.get_id(joined_graph.right_seed_node()));
        p1->set_is_reverse(false);
        p1->set_offset(0);
        
        auto splice_idxs = joined_graph.translate_node_ids(path);
        
        REQUIRE(splice_idxs.first == 0);
        REQUIRE(splice_idxs.second == 1);
        
        REQUIRE(p0->node_id() == graph.get_id(h0));
        REQUIRE(p0->is_reverse() == false);
        REQUIRE(p0->offset() == 2);
        
        REQUIRE(p1->node_id() == graph.get_id(h2));
        REQUIRE(p1->is_reverse() == false);
        REQUIRE(p1->offset() == 3);
    }
    
    {
        Path path;
        Position* p0 = path.add_mapping()->mutable_position();
        p0->set_node_id(joined_graph.get_id(joined_graph.right_seed_node()));
        p0->set_is_reverse(true);
        p0->set_offset(2);
        
        Position* p1 = path.add_mapping()->mutable_position();
        p1->set_node_id(joined_graph.get_id(joined_graph.left_seed_node()));
        p1->set_is_reverse(true);
        p1->set_offset(0);
        
        auto splice_idxs = joined_graph.translate_node_ids(path);
        
        REQUIRE(splice_idxs.first == 1);
        REQUIRE(splice_idxs.second == 0);
        
        REQUIRE(p0->node_id() == graph.get_id(h2));
        REQUIRE(p0->is_reverse() == true);
        REQUIRE(p0->offset() == 3);
        
        REQUIRE(p1->node_id() == graph.get_id(h0));
        REQUIRE(p1->is_reverse() == true);
        REQUIRE(p1->offset() == 3);
        
    }
}

TEST_CASE("JoinedSpliceGraph can correctly perform queries a more complicated graph",
          "[splice]") {
    
    HashGraph graph;
    
    handle_t h0 = graph.create_handle("GAGACCC");
    handle_t h1 = graph.create_handle("TAT");
    handle_t h2 = graph.create_handle("TGG");
    handle_t h3 = graph.create_handle("CGCGAAA");
    handle_t h4 = graph.create_handle("GAC");
    handle_t h5 = graph.create_handle("AGA");
    handle_t h6 = graph.create_handle("CGATGAC");
    
    graph.create_edge(h0, h1);
    graph.create_edge(h0, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h3, h5);
    graph.create_edge(h4, h6);
    graph.create_edge(h5, h6);
    
    pos_t left_pos(graph.get_id(h0), false, 1);
    pos_t right_pos(graph.get_id(h6), false, 6);
    int64_t search_dist = 50;
    
    IncrementalSubgraph left_subgraph(graph, left_pos, false, search_dist);
    IncrementalSubgraph right_subgraph(graph, right_pos, true, search_dist);
    
    while (left_subgraph.is_extendable()) {
        auto g = left_subgraph.extend();
    }
    while (right_subgraph.is_extendable()) {
        auto g = right_subgraph.extend();
    }
    
    handle_t left_splice_node;
    handle_t right_splice_node;
    left_subgraph.for_each_handle([&](const handle_t& handle) {
        if (left_subgraph.get_underlying_handle(handle) == h1) {
            left_splice_node = handle;
        }
    });
    right_subgraph.for_each_handle([&](const handle_t& handle) {
        if (right_subgraph.get_underlying_handle(handle) == h3) {
            right_splice_node = handle;
        }
    });
        
    JoinedSpliceGraph joined_graph(graph, left_subgraph, left_splice_node, 2,
                                   right_subgraph, right_splice_node, 3);
    
    
    REQUIRE(joined_graph.get_node_count() == 6);
    bool found1 = false, found2 = false, found3 = false, found4 = false, found5 = false, found6 = false;
    handle_t t0, t1, t2, t3, t4, t5;
    joined_graph.for_each_handle([&](const handle_t& h) {
        if (joined_graph.get_sequence(h) == "AGACCC") {
            found1 = true;
            t0 = h;
        }
        else if (joined_graph.get_sequence(h) == "TA") {
            found2 = true;
            t1 = h;
        }
        else if (joined_graph.get_sequence(h) == "GAAA") {
            found3 = true;
            t2 = h;
        }
        else if (joined_graph.get_sequence(h) == "GAC") {
            found4 = true;
            t3 = h;
        }
        else if (joined_graph.get_sequence(h) == "AGA") {
            found5 = true;
            t4 = h;
        }
        else if (joined_graph.get_sequence(h) == "CGATGA") {
            found6 = true;
            t5 = h;
        }
        else{
            REQUIRE(false);
        }
    });
    REQUIRE(found1);
    REQUIRE(found2);
    REQUIRE(found3);
    REQUIRE(found4);
    REQUIRE(found5);
    REQUIRE(found6);
    found1 = found2 = found3 = found4 = found5 = found6 = false;
    
    
    size_t count = 0;
    joined_graph.follow_edges(t0, true, [&](const handle_t& n) {
        ++count;
    });
    REQUIRE(count == 0);
    count = 0;
    joined_graph.follow_edges(t0, false, [&](const handle_t& n) {
        REQUIRE(n == t1);
        ++count;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(t1, true, [&](const handle_t& n) {
        REQUIRE(n == t0);
        ++count;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(t1, false, [&](const handle_t& n) {
        REQUIRE(n == t2);
        ++count;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(t2, true, [&](const handle_t& n) {
        REQUIRE(n == t1);
        ++count;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(t2, false, [&](const handle_t& n) {
        bool correct = n == t3 || n == t4;
        REQUIRE(correct);
        ++count;
    });
    REQUIRE(count == 2);
    count = 0;
    joined_graph.follow_edges(t3, true, [&](const handle_t& n) {
        REQUIRE(n == t2);
        ++count;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(t3, false, [&](const handle_t& n) {
        REQUIRE(n == t5);
        ++count;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(t4, true, [&](const handle_t& n) {
        REQUIRE(n == t2);
        ++count;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(t4, false, [&](const handle_t& n) {
        REQUIRE(n == t5);
        ++count;
    });
    REQUIRE(count == 1);
    count = 0;
    joined_graph.follow_edges(t5, true, [&](const handle_t& n) {
        bool correct = n == t3 || n == t4;
        REQUIRE(correct);
        ++count;
    });
    REQUIRE(count == 2);
    count = 0;
    joined_graph.follow_edges(t5, false, [&](const handle_t& n) {
        ++count;
    });
    REQUIRE(count == 0);
    count = 0;
}

}
}
        
