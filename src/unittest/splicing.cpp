/// \file splicing.cpp
///  
/// Unit tests for the splicing functions
///

#include <iostream>
#include <random>

#include "../splicing.hpp"
#include "../multipath_mapper.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../build_index.hpp"
#include "xg.hpp"
#include "catch.hpp"
#include "test_aligner.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {
using namespace std;

class MultipathMapperSpliceTest : public MultipathMapper {
public:
    MultipathMapperSpliceTest(PathPositionHandleGraph* graph, gcsa::GCSA* gcsa_index,
                              gcsa::LCPArray* lcp_array, MinimumDistanceIndex* distance_index)
        : MultipathMapper(graph, gcsa_index, lcp_array, nullptr, nullptr, distance_index)
    {
        // nothing to do
    }
    ~MultipathMapperSpliceTest() = default;
    
    using MultipathMapper::find_spliced_alignments;
    using MultipathMapper::align_to_splice_candidates;
    using MultipathMapper::test_splice_candidates;
    using MultipathMapper::extract_cluster_graph;
    using MultipathMapper::splice_stats;
    using MultipathMapper::no_splice_log_odds;
};

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
    vector<double> weights{1.0};
    vector<pair<double, double>> params{{5.0, 2.0}};

    TestAligner test_aligner;
    vector<tuple<string, string, double>> table;
    for (auto motif : motifs) {
        table.emplace_back(motif, motif, 1.0 / motifs.size());
    }
    SpliceStats splice_stats(table, weights, params, *test_aligner.get_regular_aligner());

    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, splice_stats);

    auto m0 = splice_region.candidate_splice_sites(0);
    auto m1 = splice_region.candidate_splice_sites(2);
    auto m2 = splice_region.candidate_splice_sites(4);

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
    vector<double> weights{1.0};
    vector<pair<double, double>> params{{5.0, 2.0}};

    TestAligner test_aligner;
    vector<tuple<string, string, double>> table;
    for (auto motif : motifs) {
        table.emplace_back(motif, motif, 1.0 / motifs.size());
    }
    SpliceStats splice_stats(table, weights, params, *test_aligner.get_regular_aligner());

    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, splice_stats);

    auto m0 = splice_region.candidate_splice_sites(0);
    auto m1 = splice_region.candidate_splice_sites(2);
    auto m2 = splice_region.candidate_splice_sites(4);

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
    vector<double> weights{1.0};
    vector<pair<double, double>> params{{5.0, 2.0}};

    TestAligner test_aligner;
    vector<tuple<string, string, double>> table;
    for (auto motif : motifs) {
        table.emplace_back(motif, motif, 1.0 / motifs.size());
    }
    SpliceStats splice_stats(table, weights, params, *test_aligner.get_regular_aligner());

    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, splice_stats);

    auto m0 = splice_region.candidate_splice_sites(0);
    auto m1 = splice_region.candidate_splice_sites(2);
    auto m2 = splice_region.candidate_splice_sites(4);
    auto m3 = splice_region.candidate_splice_sites(6);

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
    vector<double> weights{1.0};
    vector<pair<double, double>> params{{5.0, 2.0}};

    TestAligner test_aligner;
    vector<tuple<string, string, double>> table;
    for (auto motif : motifs) {
        table.emplace_back(motif, motif, 1.0 / motifs.size());
    }
    SpliceStats splice_stats(table, weights, params, *test_aligner.get_regular_aligner());

    SpliceRegion splice_region(pos, search_left, search_dist, graph, machine, splice_stats);

    auto m0 = splice_region.candidate_splice_sites(0);
    auto m1 = splice_region.candidate_splice_sites(2);

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
    REQUIRE(joined_graph.min_link_length() == 6);
    REQUIRE(joined_graph.max_link_length() == 6);
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
    handle_t h5 = graph.create_handle("AGAT");
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
    REQUIRE(joined_graph.min_link_length() == 21);
    REQUIRE(joined_graph.max_link_length() == 22);
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
        else if (joined_graph.get_sequence(h) == "AGAT") {
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

TEST_CASE("trim_path trims paths correctly", "[splice][multipath]") {
    
    path_t path;
    
    auto m0 = path.add_mapping();
    m0->mutable_position()->set_node_id(1);
    m0->mutable_position()->set_is_reverse(false);
    m0->mutable_position()->set_offset(2);
    auto e0 = m0->add_edit();
    e0->set_from_length(4);
    e0->set_to_length(4);
    e0->set_sequence("ACGT");
    auto e1 = m0->add_edit();
    e1->set_from_length(0);
    e1->set_to_length(2);
    e1->set_sequence("AT");

    auto m1 = path.add_mapping();
    m1->mutable_position()->set_node_id(2);
    m1->mutable_position()->set_is_reverse(false);
    m1->mutable_position()->set_offset(0);
    auto e2 = m1->add_edit();
    e2->set_from_length(4);
    e2->set_to_length(0);
    e2->set_sequence("");
    auto e3 = m1->add_edit();
    e3->set_from_length(4);
    e3->set_to_length(4);
    e3->set_sequence("");
    
    SECTION("Beginning on left") {
        bool trimmed = trim_path(&path, true, 0, 0, 0);
        REQUIRE(!trimmed);
        REQUIRE(path.mapping_size() == 2);
        REQUIRE(path.mapping(0).edit_size() == 2);
        REQUIRE(path.mapping(1).edit_size() == 2);
    }
    
    SECTION("Beginning on right") {
        bool trimmed = trim_path(&path, false, 0, 0, 0);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 0);
    }
    
    SECTION("Between mappings on left") {
        bool trimmed = trim_path(&path, true, 1, 0, 0);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 1);
        REQUIRE(path.mapping(0).position().node_id() == 2);
        REQUIRE(path.mapping(0).position().offset() == 0);
    }
    
    SECTION("Between mappings on right") {
        bool trimmed = trim_path(&path, false, 1, 0, 0);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 1);
        REQUIRE(path.mapping(0).position().node_id() == 1);
        REQUIRE(path.mapping(0).position().offset() == 2);
    }
    
    SECTION("End on left") {
        bool trimmed = trim_path(&path, true, 2, 0, 0);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 0);
    }
    
    SECTION("End on right") {
        bool trimmed = trim_path(&path, false, 2, 0, 0);
        REQUIRE(!trimmed);
        REQUIRE(path.mapping_size() == 2);
        REQUIRE(path.mapping(0).edit_size() == 2);
        REQUIRE(path.mapping(1).edit_size() == 2);
    }
    
    SECTION("Between edits of first mapping on left") {
        bool trimmed = trim_path(&path, true, 0, 1, 0);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 2);
        REQUIRE(path.mapping(0).position().node_id() == 1);
        REQUIRE(path.mapping(0).position().offset() == 6);
        REQUIRE(path.mapping(0).edit_size() == 1);
        REQUIRE(path.mapping(0).edit(0).from_length() == 0);
        REQUIRE(path.mapping(0).edit(0).to_length() == 2);
        REQUIRE(path.mapping(0).edit(0).sequence() == "AT");
        REQUIRE(path.mapping(1).position().node_id() == 2);
        REQUIRE(path.mapping(1).position().offset() == 0);
        REQUIRE(path.mapping(1).edit_size() == 2);
    }
    
    SECTION("Between edits of first mapping on right") {
        bool trimmed = trim_path(&path, false, 0, 1, 0);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 1);
        REQUIRE(path.mapping(0).position().node_id() == 1);
        REQUIRE(path.mapping(0).position().offset() == 2);
        REQUIRE(path.mapping(0).edit_size() == 1);
        REQUIRE(path.mapping(0).edit(0).from_length() == 4);
        REQUIRE(path.mapping(0).edit(0).to_length() == 4);
        REQUIRE(path.mapping(0).edit(0).sequence() == "ACGT");
    }
    
    SECTION("Between edits of second mapping on left") {
        bool trimmed = trim_path(&path, true, 1, 1, 0);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 1);
        REQUIRE(path.mapping(0).position().node_id() == 2);
        REQUIRE(path.mapping(0).position().offset() == 4);
        REQUIRE(path.mapping(0).edit_size() == 1);
        REQUIRE(path.mapping(0).edit(0).from_length() == 4);
        REQUIRE(path.mapping(0).edit(0).to_length() == 4);
        REQUIRE(path.mapping(0).edit(0).sequence() == "");
    }
    
    SECTION("Between edits of second mapping on right") {
        bool trimmed = trim_path(&path, false, 1, 1, 0);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 2);
        REQUIRE(path.mapping(0).position().node_id() == 1);
        REQUIRE(path.mapping(0).position().offset() == 2);
        REQUIRE(path.mapping(0).edit_size() == 2);
        REQUIRE(path.mapping(1).position().node_id() == 2);
        REQUIRE(path.mapping(1).position().offset() == 0);
        REQUIRE(path.mapping(1).edit_size() == 1);
        REQUIRE(path.mapping(1).edit(0).from_length() == 4);
        REQUIRE(path.mapping(1).edit(0).to_length() == 0);
        REQUIRE(path.mapping(1).edit(0).sequence() == "");
    }
    
    SECTION("Within first edit of a mapping on left") {
        bool trimmed = trim_path(&path, true, 0, 0, 1);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 2);
        REQUIRE(path.mapping(0).position().node_id() == 1);
        REQUIRE(path.mapping(0).position().offset() == 3);
        REQUIRE(path.mapping(0).edit_size() == 2);
        REQUIRE(path.mapping(0).edit(0).from_length() == 3);
        REQUIRE(path.mapping(0).edit(0).to_length() == 3);
        REQUIRE(path.mapping(0).edit(0).sequence() == "CGT");
    }
    
    SECTION("Within first edit of a mapping on right") {
        bool trimmed = trim_path(&path, false, 0, 0, 1);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 1);
        REQUIRE(path.mapping(0).position().node_id() == 1);
        REQUIRE(path.mapping(0).position().offset() == 2);
        REQUIRE(path.mapping(0).edit_size() == 1);
        REQUIRE(path.mapping(0).edit(0).from_length() == 1);
        REQUIRE(path.mapping(0).edit(0).to_length() == 1);
        REQUIRE(path.mapping(0).edit(0).sequence() == "A");
    }
    
    SECTION("Within second edit of a mapping on left") {
        bool trimmed = trim_path(&path, true, 0, 1, 1);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 2);
        REQUIRE(path.mapping(0).position().node_id() == 1);
        REQUIRE(path.mapping(0).position().offset() == 6);
        REQUIRE(path.mapping(0).edit_size() == 1);
        REQUIRE(path.mapping(0).edit(0).from_length() == 0);
        REQUIRE(path.mapping(0).edit(0).to_length() == 1);
        REQUIRE(path.mapping(0).edit(0).sequence() == "T");
    }
    
    SECTION("Within second edit of a mapping on right") {
        bool trimmed = trim_path(&path, false, 0, 1, 1);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 1);
        REQUIRE(path.mapping(0).position().node_id() == 1);
        REQUIRE(path.mapping(0).position().offset() == 2);
        REQUIRE(path.mapping(0).edit_size() == 2);
        REQUIRE(path.mapping(0).edit(0).from_length() == 4);
        REQUIRE(path.mapping(0).edit(0).to_length() == 4);
        REQUIRE(path.mapping(0).edit(0).sequence() == "ACGT");
        REQUIRE(path.mapping(0).edit(1).from_length() == 0);
        REQUIRE(path.mapping(0).edit(1).to_length() == 1);
        REQUIRE(path.mapping(0).edit(1).sequence() == "A");
    }
    
    SECTION("Within an edit on the second mapping on left") {
        bool trimmed = trim_path(&path, true, 1, 1, 1);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 1);
        REQUIRE(path.mapping(0).position().node_id() == 2);
        REQUIRE(path.mapping(0).position().offset() == 5);
        REQUIRE(path.mapping(0).edit_size() == 1);
        REQUIRE(path.mapping(0).edit(0).from_length() == 3);
        REQUIRE(path.mapping(0).edit(0).to_length() == 3);
        REQUIRE(path.mapping(0).edit(0).sequence() == "");
    }
    
    SECTION("Within an edit on the second mapping on right") {
        bool trimmed = trim_path(&path, false, 1, 1, 1);
        REQUIRE(trimmed);
        REQUIRE(path.mapping_size() == 2);
        REQUIRE(path.mapping(1).position().node_id() == 2);
        REQUIRE(path.mapping(1).position().offset() == 0);
        REQUIRE(path.mapping(1).edit_size() == 2);
        REQUIRE(path.mapping(1).edit(1).from_length() == 1);
        REQUIRE(path.mapping(1).edit(1).to_length() == 1);
        REQUIRE(path.mapping(1).edit(1).sequence() == "");
    }
}

TEST_CASE("fuse_spliced_alignments produces the correct results",
          "[splice][multipath]") {

    bdsg::HashGraph graph;

    handle_t h0 = graph.create_handle("GATAAAA");
    handle_t h1 = graph.create_handle("AAAAAAA");
    handle_t h2 = graph.create_handle("AAATACA");

    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);

    TestAligner test_aligner;


    Alignment aln;
    aln.set_sequence("GATTACA");

    multipath_alignment_t left_mp_aln;

    left_mp_aln.set_sequence("GATTACA");
    auto s0 = left_mp_aln.add_subpath();
    auto m0 = s0->mutable_path()->add_mapping();
    m0->mutable_position()->set_node_id(graph.get_id(h0));
    m0->mutable_position()->set_is_reverse(false);
    m0->mutable_position()->set_offset(0);
    auto e0 = m0->add_edit();
    e0->set_from_length(3);
    e0->set_to_length(3);
    auto e1 = m0->add_edit();
    e1->set_from_length(0);
    e1->set_to_length(4);
    e1->set_sequence("TACA");
    s0->set_score(8);

    multipath_alignment_t right_mp_aln;

    auto s1 = right_mp_aln.add_subpath();
    auto m1 = s1->mutable_path()->add_mapping();
    m1->mutable_position()->set_node_id(graph.get_id(h2));
    m1->mutable_position()->set_is_reverse(false);
    m1->mutable_position()->set_offset(3);
    auto e2 = m1->add_edit();
    e2->set_from_length(0);
    e2->set_to_length(3);
    e2->set_sequence("GAT");
    auto e3 = m1->add_edit();
    e3->set_from_length(4);
    e3->set_to_length(4);
    s1->set_score(9);

    identify_start_subpaths(left_mp_aln);
    identify_start_subpaths(right_mp_aln);


    Alignment linker;
    auto m2 = linker.mutable_path()->add_mapping();
    m2->mutable_position()->set_node_id(graph.get_id(h0));
    m2->mutable_position()->set_is_reverse(false);
    m2->mutable_position()->set_offset(3);
    auto m3 = linker.mutable_path()->add_mapping();
    m3->mutable_position()->set_node_id(graph.get_id(h2));
    m3->mutable_position()->set_is_reverse(false);
    m3->mutable_position()->set_offset(3);

    int64_t left_bridge_point = 3;
    int64_t splice_idx = 1;
    int32_t splice_score = -2;

    
    // test for the spliced alignment we expect from this
    auto test_spliced_aln = [&](multipath_alignment_t& fused) {
        REQUIRE(fused.subpath_size() == 2);
        REQUIRE(fused.subpath(0).score() == 8);
        REQUIRE(fused.subpath(0).path().mapping_size() == 1);
        REQUIRE(fused.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h0));
        REQUIRE(fused.subpath(0).path().mapping(0).position().is_reverse() == false);
        REQUIRE(fused.subpath(0).path().mapping(0).position().offset() == 0);
        REQUIRE(fused.subpath(0).path().mapping(0).edit_size() == 1);
        REQUIRE(fused.subpath(0).path().mapping(0).edit(0).from_length() == 3);
        REQUIRE(fused.subpath(0).path().mapping(0).edit(0).to_length() == 3);
        REQUIRE(fused.subpath(0).path().mapping(0).edit(0).sequence() == "");
        REQUIRE(fused.subpath(0).next_size() == 0);
        REQUIRE(fused.subpath(0).connection_size() == 1);
        REQUIRE(fused.subpath(0).connection(0).next() == 1);
        REQUIRE(fused.subpath(0).connection(0).score() == splice_score);
        
        REQUIRE(fused.subpath(1).score() == 9);
        REQUIRE(fused.subpath(1).path().mapping_size() == 1);
        REQUIRE(fused.subpath(1).path().mapping(0).position().node_id() == graph.get_id(h2));
        REQUIRE(fused.subpath(1).path().mapping(0).position().is_reverse() == false);
        REQUIRE(fused.subpath(1).path().mapping(0).position().offset() == 3);
        REQUIRE(fused.subpath(1).path().mapping(0).edit_size() == 1);
        REQUIRE(fused.subpath(1).path().mapping(0).edit(0).from_length() == 4);
        REQUIRE(fused.subpath(1).path().mapping(0).edit(0).to_length() == 4);
        REQUIRE(fused.subpath(1).path().mapping(0).edit(0).sequence() == "");
        REQUIRE(fused.subpath(1).next_size() == 0);
        REQUIRE(fused.subpath(1).connection_size() == 0);
    };

    SECTION("Fusing works when linker is empty") {

        multipath_alignment_t fused = fuse_spliced_alignments(aln, move(left_mp_aln), move(right_mp_aln), left_bridge_point,
                                                              linker, splice_idx, splice_score,
                                                              *test_aligner.get_regular_aligner(), graph);
        test_spliced_aln(fused);
    }

    m2->mutable_position()->set_offset(2);
    auto e4 = m2->add_edit();
    e4->set_from_length(1);
    e4->set_to_length(1);
    linker.set_score(1);

    left_bridge_point = 2;

    SECTION("Fusing works when linker is only empty on the right side") {
        
        multipath_alignment_t fused = fuse_spliced_alignments(aln, move(left_mp_aln), move(right_mp_aln), left_bridge_point,
                                                              linker, splice_idx, splice_score,
                                                              *test_aligner.get_regular_aligner(), graph);
        test_spliced_aln(fused);
    }

    auto e5 = m3->add_edit();
    e5->set_from_length(1);
    e5->set_to_length(1);
    linker.set_score(2);

    SECTION("Fusing works when linker is non-empty on both sides") {

        multipath_alignment_t fused = fuse_spliced_alignments(aln, move(left_mp_aln), move(right_mp_aln), left_bridge_point,
                                                              linker, splice_idx, splice_score,
                                                              *test_aligner.get_regular_aligner(), graph);
        test_spliced_aln(fused);
    }


    m2->mutable_position()->set_offset(3);
    m2->mutable_edit()->Clear();
    linker.set_score(1);

    left_bridge_point = 3;

    SECTION("Fusing works when linker is only empty on the left side") {

        multipath_alignment_t fused = fuse_spliced_alignments(aln, move(left_mp_aln), move(right_mp_aln), left_bridge_point,
                                                              linker, splice_idx, splice_score,
                                                              *test_aligner.get_regular_aligner(), graph);
        test_spliced_aln(fused);
    }
}

TEST_CASE("fuse_spliced_alignments can handle multiple splice points",
          "[splice][mulipath]") {

    bdsg::HashGraph graph;

    handle_t h0 = graph.create_handle("GATAAAA");
    handle_t h1 = graph.create_handle("AAAAAAA");
    handle_t h2 = graph.create_handle("AAATACA");

    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);

    TestAligner test_aligner;


    Alignment aln;
    aln.set_sequence("GATTACA");

    multipath_alignment_t left_mp_aln;
    left_mp_aln.set_sequence("GATTACA");

    auto s0 = left_mp_aln.add_subpath();
    auto m0 = s0->mutable_path()->add_mapping();
    m0->mutable_position()->set_node_id(graph.get_id(h0));
    m0->mutable_position()->set_is_reverse(false);
    m0->mutable_position()->set_offset(0);
    auto e0 = m0->add_edit();
    e0->set_from_length(1);
    e0->set_to_length(1);

    s0->set_score(6);
    s0->add_next(1);
    s0->add_next(2);

    auto s1 = left_mp_aln.add_subpath();
    auto m1 = s1->mutable_path()->add_mapping();
    m1->mutable_position()->set_node_id(graph.get_id(h0));
    m1->mutable_position()->set_is_reverse(false);
    m1->mutable_position()->set_offset(1);
    auto e1 = m1->add_edit();
    e1->set_from_length(2);
    e1->set_to_length(2);
    auto e2 = m1->add_edit();
    e2->set_from_length(1);
    e2->set_to_length(1);
    e2->set_sequence("T");
    auto e3 = m1->add_edit();
    e3->set_from_length(0);
    e3->set_to_length(3);
    e3->set_sequence("ACA");

    s1->set_score(-2);

    auto s2 = left_mp_aln.add_subpath();
    auto m2 = s2->mutable_path()->add_mapping();
    m2->mutable_position()->set_node_id(graph.get_id(h0));
    m2->mutable_position()->set_is_reverse(false);
    m2->mutable_position()->set_offset(1);
    auto e4 = m2->add_edit();
    e4->set_from_length(2);
    e4->set_to_length(2);
    auto e5 = m2->add_edit();
    e5->set_from_length(0);
    e5->set_to_length(3);
    e5->set_sequence("TACA");

    s2->set_score(2);

    multipath_alignment_t right_mp_aln;
    right_mp_aln.set_sequence("GATTACA");

    auto s3 = right_mp_aln.add_subpath();
    auto m3 = s3->mutable_path()->add_mapping();
    m3->mutable_position()->set_node_id(graph.get_id(h2));
    m3->mutable_position()->set_is_reverse(false);
    m3->mutable_position()->set_offset(2);
    auto e6 = m3->add_edit();
    e6->set_from_length(0);
    e6->set_to_length(2);
    e6->set_sequence("GA");
    auto e7 = m3->add_edit();
    e7->set_from_length(1);
    e7->set_to_length(1);
    e7->set_sequence("T");
    auto e8 = m3->add_edit();
    e8->set_from_length(2);
    e8->set_to_length(2);

    s3->set_score(-2);
    s3->add_next(2);

    auto s4 = right_mp_aln.add_subpath();
    auto m4 = s4->mutable_path()->add_mapping();
    m4->mutable_position()->set_node_id(graph.get_id(h2));
    m4->mutable_position()->set_is_reverse(false);
    m4->mutable_position()->set_offset(3);
    auto e9 = m4->add_edit();
    e9->set_from_length(0);
    e9->set_to_length(3);
    e9->set_sequence("GAT");
    auto e10 = m4->add_edit();
    e10->set_from_length(2);
    e10->set_to_length(2);

    s4->set_score(2);
    s4->add_next(2);

    auto s5 = right_mp_aln.add_subpath();
    auto m5 = s5->mutable_path()->add_mapping();
    m5->mutable_position()->set_node_id(graph.get_id(h2));
    m5->mutable_position()->set_is_reverse(false);
    m5->mutable_position()->set_offset(5);
    auto e11 = m5->add_edit();
    e11->set_from_length(2);
    e11->set_to_length(2);

    s5->set_score(7);

    Alignment linker;
    auto m6 = linker.mutable_path()->add_mapping();
    m6->mutable_position()->set_node_id(graph.get_id(h0));
    m6->mutable_position()->set_is_reverse(false);
    m6->mutable_position()->set_offset(3);
    auto m7 = linker.mutable_path()->add_mapping();
    m7->mutable_position()->set_node_id(graph.get_id(h2));
    m7->mutable_position()->set_is_reverse(false);
    m7->mutable_position()->set_offset(3);

    int64_t left_bridge_point = 3;
    int64_t splice_idx = 1;
    int32_t splice_score = -2;

    identify_start_subpaths(left_mp_aln);
    identify_start_subpaths(right_mp_aln);

    multipath_alignment_t fused = fuse_spliced_alignments(aln, move(left_mp_aln), move(right_mp_aln), left_bridge_point,
                                                          linker, splice_idx, splice_score,
                                                          *test_aligner.get_regular_aligner(), graph);

    REQUIRE(fused.subpath_size() == 6);

    REQUIRE(fused.subpath(0).score() == 6);
    REQUIRE(fused.subpath(0).path().mapping_size() == 1);
    REQUIRE(fused.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h0));
    REQUIRE(fused.subpath(0).path().mapping(0).position().is_reverse() == false);
    REQUIRE(fused.subpath(0).path().mapping(0).position().offset() == 0);
    REQUIRE(fused.subpath(0).path().mapping(0).edit_size() == 1);
    REQUIRE(fused.subpath(0).path().mapping(0).edit(0).from_length() == 1);
    REQUIRE(fused.subpath(0).path().mapping(0).edit(0).to_length() == 1);
    REQUIRE(fused.subpath(0).path().mapping(0).edit(0).sequence() == "");
    REQUIRE(fused.subpath(0).next_size() == 2);
    REQUIRE(fused.subpath(0).next(0) == 1);
    REQUIRE(fused.subpath(0).next(1) == 2);
    REQUIRE(fused.subpath(0).connection_size() == 0);

    for (int i : {1, 2}) {
        REQUIRE(fused.subpath(i).score() == 2);
        REQUIRE(fused.subpath(i).path().mapping_size() == 1);
        REQUIRE(fused.subpath(i).path().mapping(0).position().node_id() == graph.get_id(h0));
        REQUIRE(fused.subpath(i).path().mapping(0).position().is_reverse() == false);
        REQUIRE(fused.subpath(i).path().mapping(0).position().offset() == 1);
        REQUIRE(fused.subpath(i).path().mapping(0).edit_size() == 1);
        REQUIRE(fused.subpath(i).path().mapping(0).edit(0).from_length() == 2);
        REQUIRE(fused.subpath(i).path().mapping(0).edit(0).to_length() == 2);
        REQUIRE(fused.subpath(i).path().mapping(0).edit(0).sequence() == "");
        REQUIRE(fused.subpath(i).next_size() == 0);
        REQUIRE(fused.subpath(i).connection_size() == 2);
        REQUIRE(fused.subpath(i).connection(0).next() == 3);
        REQUIRE(fused.subpath(i).connection(1).next() == 4);
    }

    for (int i : {3, 4}) {
        REQUIRE(fused.subpath(i).score() == 2);
        REQUIRE(fused.subpath(i).path().mapping_size() == 1);
        REQUIRE(fused.subpath(i).path().mapping(0).position().node_id() == graph.get_id(h2));
        REQUIRE(fused.subpath(i).path().mapping(0).position().is_reverse() == false);
        REQUIRE(fused.subpath(i).path().mapping(0).position().offset() == 3);
        REQUIRE(fused.subpath(i).path().mapping(0).edit_size() == 1);
        REQUIRE(fused.subpath(i).path().mapping(0).edit(0).from_length() == 2);
        REQUIRE(fused.subpath(i).path().mapping(0).edit(0).to_length() == 2);
        REQUIRE(fused.subpath(i).path().mapping(0).edit(0).sequence() == "");
        REQUIRE(fused.subpath(i).next_size() == 1);
        REQUIRE(fused.subpath(i).next(0) == 5);
        REQUIRE(fused.subpath(i).connection_size() == 0);
    }

    REQUIRE(fused.subpath(5).score() == 7);
    REQUIRE(fused.subpath(5).path().mapping_size() == 1);
    REQUIRE(fused.subpath(5).path().mapping(0).position().node_id() == graph.get_id(h2));
    REQUIRE(fused.subpath(5).path().mapping(0).position().is_reverse() == false);
    REQUIRE(fused.subpath(5).path().mapping(0).position().offset() == 5);
    REQUIRE(fused.subpath(5).path().mapping(0).edit_size() == 1);
    REQUIRE(fused.subpath(5).path().mapping(0).edit(0).from_length() == 2);
    REQUIRE(fused.subpath(5).path().mapping(0).edit(0).to_length() == 2);
    REQUIRE(fused.subpath(5).path().mapping(0).edit(0).sequence() == "");
    REQUIRE(fused.subpath(5).next_size() == 0);
    REQUIRE(fused.subpath(5).connection_size() == 0);
}

TEST_CASE("fuse_spliced_alignments can handle removing some subpaths", "[splice][multipath]") {
    
    bdsg::HashGraph graph;

    handle_t h0 = graph.create_handle("GATAAAA");
    handle_t h1 = graph.create_handle("AAAAAAA");
    handle_t h2 = graph.create_handle("AAATACA");

    graph.create_edge(h0, h1);
    graph.create_edge(h1, h2);

    TestAligner test_aligner;


    Alignment aln;
    aln.set_sequence("GATTACA");

    multipath_alignment_t left_mp_aln;

    left_mp_aln.set_sequence("GATTACA");
    auto s0 = left_mp_aln.add_subpath();
    auto m0 = s0->mutable_path()->add_mapping();
    m0->mutable_position()->set_node_id(graph.get_id(h0));
    m0->mutable_position()->set_is_reverse(false);
    m0->mutable_position()->set_offset(0);
    auto e0 = m0->add_edit();
    e0->set_from_length(3);
    e0->set_to_length(3);
    
    s0->set_score(8);
    s0->add_next(1);
    
    auto s1 = left_mp_aln.add_subpath();
    auto m1 = s1->mutable_path()->add_mapping();
    m1->mutable_position()->set_node_id(graph.get_id(h0));
    m1->mutable_position()->set_is_reverse(false);
    m1->mutable_position()->set_offset(3);
    auto e1 = m1->add_edit();
    e1->set_from_length(0);
    e1->set_to_length(4);
    e1->set_sequence("TACA");
    
    s1->set_score(0);

    multipath_alignment_t right_mp_aln;

    auto s2 = right_mp_aln.add_subpath();
    auto m2 = s2->mutable_path()->add_mapping();
    m2->mutable_position()->set_node_id(graph.get_id(h2));
    m2->mutable_position()->set_is_reverse(false);
    m2->mutable_position()->set_offset(3);
    auto e2 = m2->add_edit();
    e2->set_from_length(0);
    e2->set_to_length(3);
    e2->set_sequence("GAT");
    
    s2->set_score(0);
    s2->add_next(1);
    
    auto s3 = right_mp_aln.add_subpath();
    auto m3 = s3->mutable_path()->add_mapping();
    m3->mutable_position()->set_node_id(graph.get_id(h2));
    m3->mutable_position()->set_is_reverse(false);
    m3->mutable_position()->set_offset(3);
    auto e3 = m3->add_edit();
    e3->set_from_length(4);
    e3->set_to_length(4);
    
    s3->set_score(9);
    

    identify_start_subpaths(left_mp_aln);
    identify_start_subpaths(right_mp_aln);


    Alignment linker;
    auto m4 = linker.mutable_path()->add_mapping();
    m4->mutable_position()->set_node_id(graph.get_id(h0));
    m4->mutable_position()->set_is_reverse(false);
    m4->mutable_position()->set_offset(3);
    auto m5 = linker.mutable_path()->add_mapping();
    m5->mutable_position()->set_node_id(graph.get_id(h2));
    m5->mutable_position()->set_is_reverse(false);
    m5->mutable_position()->set_offset(3);

    int64_t left_bridge_point = 3;
    int64_t splice_idx = 1;
    int32_t splice_score = -2;
    
    multipath_alignment_t fused = fuse_spliced_alignments(aln, move(left_mp_aln), move(right_mp_aln), left_bridge_point,
                                                          linker, splice_idx, splice_score,
                                                          *test_aligner.get_regular_aligner(), graph);
        
    REQUIRE(fused.subpath_size() == 2);
    REQUIRE(fused.subpath(0).score() == 8);
    REQUIRE(fused.subpath(0).path().mapping_size() == 1);
    REQUIRE(fused.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h0));
    REQUIRE(fused.subpath(0).path().mapping(0).position().is_reverse() == false);
    REQUIRE(fused.subpath(0).path().mapping(0).position().offset() == 0);
    REQUIRE(fused.subpath(0).path().mapping(0).edit_size() == 1);
    REQUIRE(fused.subpath(0).path().mapping(0).edit(0).from_length() == 3);
    REQUIRE(fused.subpath(0).path().mapping(0).edit(0).to_length() == 3);
    REQUIRE(fused.subpath(0).path().mapping(0).edit(0).sequence() == "");
    REQUIRE(fused.subpath(0).next_size() == 0);
    REQUIRE(fused.subpath(0).connection_size() == 1);
    REQUIRE(fused.subpath(0).connection(0).next() == 1);
    REQUIRE(fused.subpath(0).connection(0).score() == splice_score);

    REQUIRE(fused.subpath(1).score() == 9);
    REQUIRE(fused.subpath(1).path().mapping_size() == 1);
    REQUIRE(fused.subpath(1).path().mapping(0).position().node_id() == graph.get_id(h2));
    REQUIRE(fused.subpath(1).path().mapping(0).position().is_reverse() == false);
    REQUIRE(fused.subpath(1).path().mapping(0).position().offset() == 3);
    REQUIRE(fused.subpath(1).path().mapping(0).edit_size() == 1);
    REQUIRE(fused.subpath(1).path().mapping(0).edit(0).from_length() == 4);
    REQUIRE(fused.subpath(1).path().mapping(0).edit(0).to_length() == 4);
    REQUIRE(fused.subpath(1).path().mapping(0).edit(0).sequence() == "");
    REQUIRE(fused.subpath(1).next_size() == 0);
    REQUIRE(fused.subpath(1).connection_size() == 0);
}

TEST_CASE("Seeds can be converted into multipath alignments", "[multipath][splice]") {
    
    bdsg::HashGraph graph;
    
    handle_t h1 = graph.create_handle("GCA");
    handle_t h2 = graph.create_handle("T");
    handle_t h3 = graph.create_handle("G");
    handle_t h4 = graph.create_handle("CTGA");
    handle_t h5 = graph.create_handle("T");
    handle_t h6 = graph.create_handle("A");
    handle_t h7 = graph.create_handle("TGAC");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    graph.create_edge(h4, h6);
    graph.create_edge(h5, h7);
    graph.create_edge(h6, h7);
    
    string read = string("GATCTGAATGC");
    Alignment aln;
    aln.set_sequence(read);
    
    pos_t hit_pos(graph.get_id(h1), false, 2);
    MaximalExactMatch mem;
    mem.begin = aln.sequence().begin() + 1;
    mem.end = aln.sequence().end() - 1;
    
    TestAligner test_aligner;
    auto& alnr = *test_aligner.get_regular_aligner();
    
    auto mp_aln = from_hit(aln, graph, hit_pos, mem, alnr);
    
    REQUIRE(mp_aln.sequence() == aln.sequence());
    REQUIRE(mp_aln.subpath_size() == 1);
    REQUIRE(mp_aln.subpath(0).score() == 9);
    REQUIRE(mp_aln.subpath(0).path().mapping_size() == 5);
    REQUIRE(mp_aln.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h1));
    REQUIRE(mp_aln.subpath(0).path().mapping(0).position().is_reverse() == false);
    REQUIRE(mp_aln.subpath(0).path().mapping(0).position().offset() == 2);
    REQUIRE(mp_aln.subpath(0).path().mapping(0).edit_size() == 2);
    REQUIRE(mp_aln.subpath(0).path().mapping(0).edit(0).from_length() == 0);
    REQUIRE(mp_aln.subpath(0).path().mapping(0).edit(0).to_length() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(0).edit(1).from_length() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(0).edit(1).to_length() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(1).position().node_id() == graph.get_id(h2));
    REQUIRE(mp_aln.subpath(0).path().mapping(1).position().is_reverse() == false);
    REQUIRE(mp_aln.subpath(0).path().mapping(1).position().offset() == 0);
    REQUIRE(mp_aln.subpath(0).path().mapping(1).edit_size() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(1).edit(0).from_length() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(1).edit(0).to_length() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(2).position().node_id() == graph.get_id(h4));
    REQUIRE(mp_aln.subpath(0).path().mapping(2).position().is_reverse() == false);
    REQUIRE(mp_aln.subpath(0).path().mapping(2).position().offset() == 0);
    REQUIRE(mp_aln.subpath(0).path().mapping(2).edit_size() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(2).edit(0).from_length() == 4);
    REQUIRE(mp_aln.subpath(0).path().mapping(2).edit(0).to_length() == 4);
    REQUIRE(mp_aln.subpath(0).path().mapping(3).position().node_id() == graph.get_id(h6));
    REQUIRE(mp_aln.subpath(0).path().mapping(3).position().is_reverse() == false);
    REQUIRE(mp_aln.subpath(0).path().mapping(3).position().offset() == 0);
    REQUIRE(mp_aln.subpath(0).path().mapping(3).edit_size() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(3).edit(0).from_length() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(3).edit(0).to_length() == 1);
    REQUIRE(mp_aln.subpath(0).path().mapping(4).position().node_id() == graph.get_id(h7));
    REQUIRE(mp_aln.subpath(0).path().mapping(4).position().is_reverse() == false);
    REQUIRE(mp_aln.subpath(0).path().mapping(4).position().offset() == 0);
    REQUIRE(mp_aln.subpath(0).path().mapping(4).edit_size() == 2);
    REQUIRE(mp_aln.subpath(0).path().mapping(4).edit(0).from_length() == 2);
    REQUIRE(mp_aln.subpath(0).path().mapping(4).edit(0).to_length() == 2);
    REQUIRE(mp_aln.subpath(0).path().mapping(4).edit(1).from_length() == 0);
    REQUIRE(mp_aln.subpath(0).path().mapping(4).edit(1).to_length() == 1);
}

TEST_CASE("MultipathMapper can make splice alignments from different candidates", "[multipath][splice]") {
    
    bdsg::HashGraph graph;
    
    handle_t h1 = graph.create_handle("CAAATTAGGGATGTGTAGATGATGATGT");
    handle_t h2 = graph.create_handle("TGCAGCATGCACAAGATCACCGATACCA");
    handle_t h3 = graph.create_handle("GGTAGAGACCAGTAGAGGAGCCCCTTTG");
    handle_t h4 = graph.create_handle("AGGATAAAGAAATGAGTAACCACGTACC");
    handle_t h5 = graph.create_handle("GAGTATTTAATATATGGATGGATTTTCA");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    
    path_handle_t ph = graph.create_path_handle("p");
    graph.append_step(ph, h1);
    graph.append_step(ph, h2);
    graph.append_step(ph, h3);
    graph.append_step(ph, h4);
    graph.append_step(ph, h5);
    
    // Make GCSA quiet
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Make pointers to fill in
    gcsa::GCSA* gcsaidx = nullptr;
    gcsa::LCPArray* lcpidx = nullptr;
    
    // Build the GCSA index
    build_gcsa_lcp(graph, gcsaidx, lcpidx, 16, 3);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    // Get snarls
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlManager snarl_manager = snarl_finder.find_snarls_parallel();
    
    // Build the distance index
    MinimumDistanceIndex distance_index(&xg_index, &snarl_manager);
    
    MultipathMapperSpliceTest mapper(&xg_index, gcsaidx, lcpidx, &distance_index);
    // slack this up a bit so it doesn't filter out good candidates in the small graph
    // otherwise, we'd need to calibrate, which i don't want to do in a test
    mapper.set_log_odds_against_splice(10.0 * 1.4);
    
    Alignment aln;
    aln.set_sequence(string("CAAATTAGGGATGTGTAGATGATGAT") + string("TATTTAATATATGGATGGATTTTCA"));
    
    multipath_alignment_t mp_aln;
    mp_aln.set_sequence(aln.sequence());
    auto s = mp_aln.add_subpath();
    s->set_score(26 + 5);
    auto p = s->mutable_path();
    auto m = p->add_mapping();
    m->mutable_position()->set_node_id(graph.get_id(h1));
    m->mutable_position()->set_is_reverse(false);
    m->mutable_position()->set_offset(0);
    auto e0 = m->add_edit();
    e0->set_from_length(26);
    e0->set_to_length(26);
    auto e1 = m->add_edit();
    e1->set_from_length(0);
    e1->set_to_length(25);
    e1->set_sequence("TATTTAATATATGGATGGATTTTCA");
    
    identify_start_subpaths(mp_aln);
    
    vector<multipath_alignment_t> mp_alns(1, mp_aln);
    vector<double> mults(1, 1.0);
    
    SECTION("From a MEM candidate") {

        MaximalExactMatch mem;
        mem.begin = aln.sequence().begin() + 26;
        mem.end = aln.sequence().end();
        mem.queried_count = 1;
        mem.match_count = 1;

        pos_t pos(graph.get_id(h5), false, 3);

        mem.nodes.push_back(gcsa::Node::encode(id(pos), offset(pos), is_rev(pos)));

        vector<MaximalExactMatch> mems(1, mem);
        vector<size_t> c_idx(1, 0);

        vector<MultipathMapper::clustergraph_t> cluster_graphs;

        mapper.find_spliced_alignments(aln, mp_alns, mults, c_idx, mems, cluster_graphs);

        REQUIRE(mp_alns.size() == 1);
        auto& m = mp_alns.front();
        REQUIRE(optimal_alignment_score(m) == aln.sequence().size() + 10 - mapper.no_splice_log_odds + mapper.splice_stats.intron_length_score(89));

    }
    
    SECTION("From a cluster candidate") {
        
        MaximalExactMatch mem1;
        mem1.begin = aln.sequence().begin() + 26;
        mem1.end = aln.sequence().begin() + 38;
        mem1.queried_count = 1;
        mem1.match_count = 1;
        mem1.nodes.push_back(gcsa::Node::encode(graph.get_id(h5), 3, false));
        
        MaximalExactMatch mem2;
        mem2.begin = aln.sequence().begin() + 39;
        mem2.end = aln.sequence().end();
        mem2.queried_count = 1;
        mem2.match_count = 1;
        mem2.nodes.push_back(gcsa::Node::encode(graph.get_id(h5), 16, false));
        
        vector<MaximalExactMatch> mems{mem1, mem2};
        vector<size_t> c_idx(1, 1);
        
        vector<MultipathMapper::clustergraph_t> cluster_graphs;
        cluster_graphs.emplace_back();
        auto& cluster_graph = cluster_graphs.back();
        get<1>(cluster_graph).first.emplace_back(&mems.front(), make_pos_t(mems.front().nodes.front()));
        get<1>(cluster_graph).first.emplace_back(&mems.back(), make_pos_t(mems.back().nodes.front()));
        get<1>(cluster_graph).second = 1.0;
        get<0>(cluster_graph) = mapper.extract_cluster_graph(aln, get<1>(cluster_graph)).first;
        get<2>(cluster_graph) = mem1.length() + mem2.length();
        
        mapper.find_spliced_alignments(aln, mp_alns, mults, c_idx, mems, cluster_graphs);
        
        REQUIRE(mp_alns.size() == 1);
        auto& m = mp_alns.front();
        REQUIRE(optimal_alignment_score(m) == aln.sequence().size() + 10 - mapper.no_splice_log_odds + mapper.splice_stats.intron_length_score(89));
        
    }
}

}
}
        
