///
/// \file position.cpp
///  
/// Unit tests for Position and pos_t manipulation
///

#include "catch.hpp"
#include "surjector.hpp"
#include "aligner.hpp"
#include "annotation.hpp"

#include "bdsg/hash_graph.hpp"
#include "bdsg/overlays/path_position_overlays.hpp"
#include <vg/vg.pb.h>

namespace vg {
namespace unittest {

class TestSurjector : public Surjector {
public:
    TestSurjector(const PathPositionHandleGraph* graph) : Surjector(graph) {}
    ~TestSurjector() = default;
    
    using Surjector::extract_overlapping_paths;
    using Surjector::filter_redundant_path_chunks;
    
};


TEST_CASE( "Spliced surject algorithm preserves deletions against the path", "[surject]" ) {
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GTCGT");
    handle_t h2 = graph.create_handle("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    handle_t h3 = graph.create_handle("TCCTTGC");
    handle_t h4 = graph.create_handle("A");
    handle_t h5 = graph.create_handle("T");
    handle_t h6 = graph.create_handle("GCCGA");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h3, h5);
    graph.create_edge(h4, h6);
    graph.create_edge(h5, h6);
    
    path_handle_t p = graph.create_path_handle("p");
    graph.append_step(p, h1);
    graph.append_step(p, h2);
    graph.append_step(p, h3);
    graph.append_step(p, h4);
    graph.append_step(p, h6);
    
    bdsg::PositionOverlay pos_graph(&graph);
    Surjector surjector(&pos_graph);
    
    vector<handle_t> read_path{h1, h3, h5, h6};
    
    Alignment read;
    string seq;
    Path* rpath = read.mutable_path();
    for (handle_t h : read_path) {
        Mapping* m = rpath->add_mapping();
        m->set_rank(rpath->mapping_size());
        m->mutable_position()->set_node_id(pos_graph.get_id(h));
        Edit* e = m->add_edit();
        e->set_from_length(pos_graph.get_length(h));
        e->set_to_length(pos_graph.get_length(h));
        
        seq += pos_graph.get_sequence(h);
    }
    read.set_sequence(seq);
    
    read.set_score(Aligner().scorer->score_contiguous_alignment(read));
    
    unordered_set<path_handle_t> paths{p};
    vector<Alignment> surjected_alns = surjector.surject(read, paths, true, true);
    REQUIRE(surjected_alns.size() == 1);
    auto& surjected = surjected_alns.front();
    
    vector<handle_t> surjected_path{h1, h3, h4, h6};
    
    REQUIRE(surjected.path().mapping_size() == read_path.size());
    
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(surjected.path().mapping(i).position().node_id() == graph.get_id(surjected_path[i]));
        REQUIRE(!surjected.path().mapping(i).position().is_reverse());
    }
    
    // should not penalize long deletions (assumed to be splices)
    REQUIRE(surjected.score() == read.score() - Aligner().scorer->mismatch - Aligner().scorer->match);
    REQUIRE(surjected.refpos_size() == 1);
    REQUIRE(surjected.refpos(0).name() == graph.get_path_name(p));
    REQUIRE(surjected.refpos(0).offset() == 0);
    
    
    
    Alignment rev_read;
    rev_read.set_sequence(reverse_complement(seq));
    
    Path* rev_rpath = rev_read.mutable_path();
    for (size_t i = 0; i < read_path.size(); ++i) {
        handle_t h = read_path[read_path.size() - i - 1];
        Mapping* m = rev_rpath->add_mapping();
        m->set_rank(rev_rpath->mapping_size());
        m->mutable_position()->set_node_id(pos_graph.get_id(h));
        m->mutable_position()->set_is_reverse(true);
        Edit* e = m->add_edit();
        e->set_from_length(pos_graph.get_length(h));
        e->set_to_length(pos_graph.get_length(h));
    }
    
    rev_read.set_score(Aligner().scorer->score_contiguous_alignment(rev_read));
    
    vector<Alignment> rev_surjected_alns = surjector.surject(rev_read, paths, true, true);
    REQUIRE(rev_surjected_alns.size() == 1);
    auto& rev_surjected = rev_surjected_alns.front();
    
    REQUIRE(rev_surjected.path().mapping_size() == read_path.size());
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(rev_surjected.path().mapping(i).position().node_id()
                == graph.get_id(surjected_path[surjected_path.size() - i - 1]));
        REQUIRE(rev_surjected.path().mapping(i).position().is_reverse());
    }
    
    // should not penalize long deletions (assumed to be splices)
    REQUIRE(rev_surjected.score() == read.score() - Aligner().scorer->mismatch - Aligner().scorer->match);
    REQUIRE(rev_surjected.refpos_size() == 1);
    REQUIRE(rev_surjected.refpos(0).name() == graph.get_path_name(p));
    REQUIRE(rev_surjected.refpos(0).offset() == 0);
    
}

TEST_CASE( "Spliced surject algorithm works when a read touches the same path in both orientations", "[surject]" ) {
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GGGGGGGGGGGGGGG");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("C");
    handle_t h4 = graph.create_handle("TTTTTTTTT");
    handle_t h5 = graph.create_handle("AAA");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, graph.flip(h4));
    graph.create_edge(h5, h4);
    graph.create_edge(h4, h4);
    graph.create_edge(graph.flip(h5), h1);
    
    path_handle_t p = graph.create_path_handle("p");
    graph.append_step(p, h1);
    graph.append_step(p, h2);
    graph.append_step(p, h4);
    graph.append_step(p, h4);
    graph.append_step(p, graph.flip(h4));
    graph.append_step(p, graph.flip(h5));
    
    bdsg::PositionOverlay pos_graph(&graph);
    Surjector surjector(&pos_graph);
    
    vector<handle_t> read_path{h1, h3, h4, h4};
    
    Alignment read;
    string seq;
    Path* rpath = read.mutable_path();
    for (handle_t h : read_path) {
        Mapping* m = rpath->add_mapping();
        m->set_rank(rpath->mapping_size());
        m->mutable_position()->set_node_id(pos_graph.get_id(h));
        Edit* e = m->add_edit();
        e->set_from_length(pos_graph.get_length(h));
        e->set_to_length(pos_graph.get_length(h));
        
        seq += pos_graph.get_sequence(h);
    }
    read.set_sequence(seq);
    
    read.set_score(Aligner().scorer->score_contiguous_alignment(read));
    
    unordered_set<path_handle_t> paths{p};
    vector<Alignment> surjected_alns = surjector.surject(read, paths, true, true);
    REQUIRE(surjected_alns.size() == 1);
    auto& surjected = surjected_alns.front();
    
    vector<handle_t> surjected_path{h1, h2, h4, h4};
    
    REQUIRE(surjected.path().mapping_size() == read_path.size());
    
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(surjected.path().mapping(i).position().node_id() == graph.get_id(surjected_path[i]));
        REQUIRE(surjected.path().mapping(i).position().is_reverse() == graph.get_is_reverse(surjected_path[i]));
    }
    
    Alignment rev_read;
    rev_read.set_sequence(reverse_complement(seq));
    
    Path* rev_rpath = rev_read.mutable_path();
    for (size_t i = 0; i < read_path.size(); ++i) {
        handle_t h = read_path[read_path.size() - i - 1];
        Mapping* m = rev_rpath->add_mapping();
        m->set_rank(rev_rpath->mapping_size());
        m->mutable_position()->set_node_id(pos_graph.get_id(h));
        m->mutable_position()->set_is_reverse(true);
        Edit* e = m->add_edit();
        e->set_from_length(pos_graph.get_length(h));
        e->set_to_length(pos_graph.get_length(h));
    }
    
    rev_read.set_score(Aligner().scorer->score_contiguous_alignment(rev_read));
    
    vector<Alignment> rev_surjected_alns = surjector.surject(rev_read, paths, true, true);
    REQUIRE(rev_surjected_alns.size() == 1);
    auto& rev_surjected = rev_surjected_alns.front();
    
    REQUIRE(rev_surjected.path().mapping_size() == read_path.size());
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(rev_surjected.path().mapping(i).position().node_id()
                == graph.get_id(surjected_path[surjected_path.size() - i - 1]));
        REQUIRE(rev_surjected.path().mapping(i).position().is_reverse()
                == !graph.get_is_reverse(surjected_path[surjected_path.size() - i - 1]));
    }
}

TEST_CASE("Path overlapping segments can be identified from multipath alignment",
          "[surject][multipath]"){
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GTCGT");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("T");
    handle_t h4 = graph.create_handle("TTAGAC");
    handle_t h5 = graph.create_handle("GCA");
    handle_t h6 = graph.create_handle("ATTAGACGCA");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    graph.create_edge(h5, h6);
    
    path_handle_t p = graph.create_path_handle("p");
    step_handle_t st0 = graph.append_step(p, h1);
    step_handle_t st1 = graph.append_step(p, h2);
    step_handle_t st2 = graph.append_step(p, h4);
    step_handle_t st3 = graph.append_step(p, h5);
    step_handle_t st4 = graph.append_step(p, h6);
    
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    
    multipath_alignment_t mp_aln;
    mp_aln.set_sequence("CGTATTAGACGC");
    
    auto s0 = mp_aln.add_subpath();
    auto m0 = s0->mutable_path()->add_mapping();
    m0->mutable_position()->set_node_id(graph.get_id(h1));
    m0->mutable_position()->set_is_reverse(false);
    m0->mutable_position()->set_offset(2);
    auto e0 = m0->add_edit();
    e0->set_from_length(3);
    e0->set_to_length(3);
    e0->set_sequence("");
    
    s0->add_next(1);
    s0->add_next(2);
    
    auto s1 = mp_aln.add_subpath();
    auto m1 = s1->mutable_path()->add_mapping();
    m1->mutable_position()->set_node_id(graph.get_id(h2));
    m1->mutable_position()->set_is_reverse(false);
    m1->mutable_position()->set_offset(0);
    auto e1 = m1->add_edit();
    e1->set_from_length(1);
    e1->set_to_length(1);
    e1->set_sequence("");
    
    s1->add_next(3);
    
    auto s2 = mp_aln.add_subpath();
    auto m2 = s2->mutable_path()->add_mapping();
    m2->mutable_position()->set_node_id(graph.get_id(h3));
    m2->mutable_position()->set_is_reverse(false);
    m2->mutable_position()->set_offset(0);
    auto e2 = m2->add_edit();
    e2->set_from_length(1);
    e2->set_to_length(1);
    e2->set_sequence("A");
    
    s2->add_next(3);
    
    auto s3 = mp_aln.add_subpath();
    auto m3 = s3->mutable_path()->add_mapping();
    m3->mutable_position()->set_node_id(graph.get_id(h4));
    m3->mutable_position()->set_is_reverse(false);
    m3->mutable_position()->set_offset(0);
    auto e3 = m3->add_edit();
    e3->set_from_length(6);
    e3->set_to_length(6);
    e3->set_sequence("");
    auto m4 = s3->mutable_path()->add_mapping();
    m4->mutable_position()->set_node_id(graph.get_id(h5));
    m4->mutable_position()->set_is_reverse(false);
    m4->mutable_position()->set_offset(0);
    auto e4 = m4->add_edit();
    e4->set_from_length(2);
    e4->set_to_length(2);
    e4->set_sequence("");
    
    identify_start_subpaths(mp_aln);
    
    unordered_set<path_handle_t> surjection_paths{p};
    
    function<int64_t(int64_t)> node_length = [&](int64_t node_id) {
        return int64_t(graph.get_length(graph.get_handle(node_id)));
    };
    
    SECTION("Forward strand of path") {
        
        unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>> connections;
        
        auto overlaps = surjector.extract_overlapping_paths(&pos_graph, mp_aln,
                                                            surjection_paths,
                                                            connections);
        
        auto fp = make_pair(p, false);
        
        REQUIRE(overlaps.count(fp));
        REQUIRE(overlaps.size() == 1);
        REQUIRE(connections.empty());
        
        auto& p_overlaps = overlaps[fp];
        
        REQUIRE(p_overlaps.first.size() == 1);
        REQUIRE(p_overlaps.first.front().first.first == mp_aln.sequence().begin());
        REQUIRE(p_overlaps.first.front().first.second == mp_aln.sequence().end());
        REQUIRE(p_overlaps.first.front().second.mapping_size() == 4);
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().node_id() == graph.get_id(h1));
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().is_reverse() == false);
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().offset() == 2);
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().node_id() == graph.get_id(h2));
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().is_reverse() == false);
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().offset() == 0);
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().node_id() == graph.get_id(h4));
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().is_reverse() == false);
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().offset() == 0);
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().node_id() == graph.get_id(h5));
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().is_reverse() == false);
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().offset() == 0);
        
        REQUIRE(p_overlaps.second.size() == 1);
        REQUIRE(p_overlaps.second.front().first == st0);
        REQUIRE(p_overlaps.second.front().second == st3);
    }
    
    SECTION("Reverse strand of path"){
        
        unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>> connections;
        
        multipath_alignment_t rev_mp_aln;
        rev_comp_multipath_alignment(mp_aln, node_length, rev_mp_aln);
        
        auto overlaps = surjector.extract_overlapping_paths(&pos_graph, rev_mp_aln,
                                                            surjection_paths,
                                                            connections);
        
        auto rp = make_pair(p, true);
            
        REQUIRE(overlaps.count(rp));
        REQUIRE(overlaps.size() == 1);
        REQUIRE(connections.empty());
        
        auto& p_overlaps = overlaps[rp];
        
        REQUIRE(p_overlaps.first.size() == 1);
        REQUIRE(p_overlaps.first.front().first.first == rev_mp_aln.sequence().begin());
        REQUIRE(p_overlaps.first.front().first.second == rev_mp_aln.sequence().end());
        REQUIRE(p_overlaps.first.front().second.mapping_size() == 4);
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().node_id() == graph.get_id(h5));
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().is_reverse() == true);
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().offset() == 1);
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().node_id() == graph.get_id(h4));
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().is_reverse() == true);
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().offset() == 0);
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().node_id() == graph.get_id(h2));
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().is_reverse() == true);
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().offset() == 0);
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().node_id() == graph.get_id(h1));
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().is_reverse() == true);
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().offset() == 0);
        
        REQUIRE(p_overlaps.second.size() == 1);
        REQUIRE(p_overlaps.second.front().first == st3);
        REQUIRE(p_overlaps.second.front().second == st0);
    }
    
    //ATTAGACGCA
    auto s4 = mp_aln.add_subpath();
    auto m5 = s4->mutable_path()->add_mapping();
    m5->mutable_position()->set_node_id(graph.get_id(h6));
    m5->mutable_position()->set_is_reverse(false);
    m5->mutable_position()->set_offset(0);
    auto e5 = m5->add_edit();
    e5->set_from_length(10);
    e5->set_to_length(10);
    e5->set_sequence("");
    
    auto s5 = mp_aln.add_subpath();
    auto m6 = s5->mutable_path()->add_mapping();
    m6->mutable_position()->set_node_id(graph.get_id(h6));
    m6->mutable_position()->set_is_reverse(false);
    m6->mutable_position()->set_offset(1);
    auto e6 = m6->add_edit();
    e6->set_from_length(9);
    e6->set_to_length(9);
    e6->set_sequence("");
    
    auto c0 = mp_aln.mutable_subpath(0)->add_connection();
    c0->set_next(4);
    c0->set_score(-2);
    auto c1 = mp_aln.mutable_subpath(1)->add_connection();
    c1->set_next(5);
    c1->set_score(-1);
    
    SECTION("Connections break segments and are recorded correctly") {
        
        unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>> connections;
        
        auto overlaps = surjector.extract_overlapping_paths(&pos_graph, mp_aln,
                                                            surjection_paths,
                                                            connections);
        
        auto fp = make_pair(p, false);
        
        REQUIRE(overlaps.count(fp));
        REQUIRE(overlaps.size() == 1);
        
        auto& p_overlaps = overlaps[fp];
        
        REQUIRE(p_overlaps.first.size() == 5);
        REQUIRE(p_overlaps.second.size() == 5);
        
        REQUIRE(connections.size() == 1);
        REQUIRE(connections.count(fp));
        
        auto& p_connections = connections[fp];
        
        REQUIRE(p_connections.size() == 2);
        
        for (auto& connection : p_connections) {
            if (get<2>(connection) == -2) {
                REQUIRE(p_overlaps.first[get<0>(connection)].second.mapping(0).position().node_id() == graph.get_id(h1));
                REQUIRE(p_overlaps.first[get<1>(connection)].second.mapping(0).position().node_id() == graph.get_id(h6));
                REQUIRE(p_overlaps.first[get<1>(connection)].second.mapping(0).position().offset() == 0);

            }
            else if (get<2>(connection) == -1) {
                REQUIRE(p_overlaps.first[get<0>(connection)].second.mapping(0).position().node_id() == graph.get_id(h2));
                REQUIRE(p_overlaps.first[get<1>(connection)].second.mapping(0).position().node_id() == graph.get_id(h6));
                REQUIRE(p_overlaps.first[get<1>(connection)].second.mapping(0).position().offset() == 1);
            }
            else {
                REQUIRE(false);
            }
        }
    }
}

TEST_CASE("Multipath alignments can be surjected", "[surject][multipath]") {

    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("T");
    handle_t h3 = graph.create_handle("TTAGAC");
    handle_t h4 = graph.create_handle("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    handle_t h5 = graph.create_handle("GCA");
    handle_t h6 = graph.create_handle("G");
    handle_t h7 = graph.create_handle("T");
    handle_t h8 = graph.create_handle("TTAGAC");
    
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h3, h5);
    graph.create_edge(h4, h5);
    graph.create_edge(h5, h6);
    graph.create_edge(h5, h7);
    graph.create_edge(h6, h8);
    graph.create_edge(h7, h8);
    
    path_handle_t p = graph.create_path_handle("p");
    step_handle_t st0 = graph.append_step(p, h1);
    step_handle_t st1 = graph.append_step(p, h3);
    step_handle_t st2 = graph.append_step(p, h4);
    step_handle_t st3 = graph.append_step(p, h5);
    step_handle_t st4 = graph.append_step(p, h6);
    step_handle_t st5 = graph.append_step(p, h8);
    
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    
    multipath_alignment_t mp_aln;
    mp_aln.set_sequence("TTTAGACGCTAGA");
    
    auto s0 = mp_aln.add_subpath();
    auto m0 = s0->mutable_path()->add_mapping();
    m0->mutable_position()->set_node_id(graph.get_id(h1));
    m0->mutable_position()->set_is_reverse(false);
    m0->mutable_position()->set_offset(0);
    auto e0 = m0->add_edit();
    e0->set_from_length(1);
    e0->set_to_length(1);
    e0->set_sequence("T");
    
    s0->set_score(1);
    s0->add_next(2);
    
    auto s1 = mp_aln.add_subpath();
    auto m1 = s1->mutable_path()->add_mapping();
    m1->mutable_position()->set_node_id(graph.get_id(h2));
    m1->mutable_position()->set_is_reverse(false);
    m1->mutable_position()->set_offset(0);
    auto e1 = m1->add_edit();
    e1->set_from_length(1);
    e1->set_to_length(1);
    e1->set_sequence("");
    
    s1->set_score(6);
    s1->add_next(2);
    
    auto s2 = mp_aln.add_subpath();
    auto m2 = s2->mutable_path()->add_mapping();
    m2->mutable_position()->set_node_id(graph.get_id(h3));
    m2->mutable_position()->set_is_reverse(false);
    m2->mutable_position()->set_offset(0);
    auto e2 = m2->add_edit();
    e2->set_from_length(6);
    e2->set_to_length(6);
    e2->set_sequence("");
    
    s2->set_score(6);
    s2->add_next(3);
    
    auto s3 = mp_aln.add_subpath();
    auto m3 = s3->mutable_path()->add_mapping();
    m3->mutable_position()->set_node_id(graph.get_id(h5));
    m3->mutable_position()->set_is_reverse(false);
    m3->mutable_position()->set_offset(0);
    auto e3 = m3->add_edit();
    e3->set_from_length(2);
    e3->set_to_length(2);
    e3->set_sequence("");
    
    s3->set_score(2);
    auto c0 = s3->add_connection();
    c0->set_next(4);
    c0->set_score(-1);
    
    auto s4 = mp_aln.add_subpath();
    auto m4 = s4->mutable_path()->add_mapping();
    m4->mutable_position()->set_node_id(graph.get_id(h8));
    m4->mutable_position()->set_is_reverse(false);
    m4->mutable_position()->set_offset(1);
    auto e4 = m4->add_edit();
    e4->set_from_length(4);
    e4->set_to_length(4);
    e4->set_sequence("");
    
    s4->set_score(9);
    
    identify_start_subpaths(mp_aln);
    
    multipath_alignment_t rc_mp_aln;
    rev_comp_multipath_alignment(mp_aln, [&](int64_t node_id) {
        return (int64_t) graph.get_length(graph.get_handle(node_id));
    }, rc_mp_aln);
    
    SECTION("With deletion splices allowed") {
        
        string path_name;
        int64_t path_pos;
        bool path_rev;
        unordered_set<path_handle_t> paths{p};
        vector<tuple<string, int64_t, bool>> positions;
        auto surjected_alns = surjector.surject(mp_aln, paths, positions, true, true);
        REQUIRE(surjected_alns.size() == 1);
        auto& surjected = surjected_alns.front();
        tie(path_name, path_pos, path_rev) = positions.front();
        
        REQUIRE(path_name == graph.get_path_name(p));
        REQUIRE(path_rev == false);
        REQUIRE(path_pos == 0);
        
        // normalize the representation
        merge_non_branching_subpaths(surjected);
        
        REQUIRE(surjected.sequence() == mp_aln.sequence());
        REQUIRE(surjected.subpath_size() == 2);
        REQUIRE(surjected.subpath(0).score() == 1 + 6 + 2);
        REQUIRE(surjected.subpath(0).path().mapping_size() == 3);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h1));
        REQUIRE(surjected.subpath(0).path().mapping(0).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(1).position().node_id() == graph.get_id(h3));
        REQUIRE(surjected.subpath(0).path().mapping(1).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(1).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(2).position().node_id() == graph.get_id(h5));
        REQUIRE(surjected.subpath(0).path().mapping(2).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(2).position().offset() == 0);
        REQUIRE(surjected.subpath(0).next_size() == 0);
        REQUIRE(surjected.subpath(0).connection_size() == 1);
        REQUIRE(surjected.subpath(0).connection(0).next() == 1);
        REQUIRE(surjected.subpath(0).connection(0).score() == -1);
        REQUIRE(surjected.subpath(1).score() == 9);
        REQUIRE(surjected.subpath(1).path().mapping_size() == 1);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().node_id() == graph.get_id(h8));
        REQUIRE(surjected.subpath(1).path().mapping(0).position().is_reverse() == false);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().offset() == 1);
        REQUIRE(surjected.subpath(1).next_size() == 0);
        REQUIRE(surjected.subpath(1).connection_size() == 0);
    }
    
    SECTION("Reverse with deletion splices allowed") {
        
        string path_name;
        int64_t path_pos;
        bool path_rev;
        unordered_set<path_handle_t> paths{p};
        vector<tuple<string, int64_t, bool>> positions;
        auto surjected_alns = surjector.surject(rc_mp_aln, paths, positions, true, true);
        REQUIRE(surjected_alns.size() == 1);
        auto& surjected = surjected_alns.front();
        tie(path_name, path_pos, path_rev) = positions.front();
        
        REQUIRE(path_name == graph.get_path_name(p));
        REQUIRE(path_rev == true);
        REQUIRE(path_pos == 0);
        
        // normalize the representation
        merge_non_branching_subpaths(surjected);
        
        REQUIRE(surjected.sequence() == rc_mp_aln.sequence());
        REQUIRE(surjected.subpath_size() == 2);
        REQUIRE(surjected.subpath(0).score() == 9);
        REQUIRE(surjected.subpath(0).path().mapping_size() == 1);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h8));
        REQUIRE(surjected.subpath(0).path().mapping(0).position().is_reverse() == true);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().offset() == 1);
        REQUIRE(surjected.subpath(0).connection_size() == 1);
        REQUIRE(surjected.subpath(0).connection(0).next() == 1);
        REQUIRE(surjected.subpath(0).connection(0).score() == -1);
        REQUIRE(surjected.subpath(0).score() == 1 + 6 + 2);
        REQUIRE(surjected.subpath(1).path().mapping_size() == 3);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().node_id() == graph.get_id(h5));
        REQUIRE(surjected.subpath(1).path().mapping(0).position().is_reverse() == true);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().offset() == 1);
        REQUIRE(surjected.subpath(1).path().mapping(1).position().node_id() == graph.get_id(h3));
        REQUIRE(surjected.subpath(1).path().mapping(1).position().is_reverse() == true);
        REQUIRE(surjected.subpath(1).path().mapping(1).position().offset() == 0);
        REQUIRE(surjected.subpath(1).path().mapping(2).position().node_id() == graph.get_id(h1));
        REQUIRE(surjected.subpath(1).path().mapping(2).position().is_reverse() == true);
        REQUIRE(surjected.subpath(1).path().mapping(2).position().offset() == 0);
        REQUIRE(surjected.subpath(1).next_size() == 0);
        REQUIRE(surjected.subpath(1).next_size() == 0);
        REQUIRE(surjected.subpath(1).connection_size() == 0);
    }
    
    SECTION("Without deletion splices allowed") {
        
        string path_name;
        int64_t path_pos;
        bool path_rev;
        unordered_set<path_handle_t> paths{p};
        vector<tuple<string, int64_t, bool>> positions;
        auto surjected_alns = surjector.surject(mp_aln, paths, positions, true, false);
        REQUIRE(surjected_alns.size() == 1);
        auto& surjected = surjected_alns.front();
        tie(path_name, path_pos, path_rev) = positions.front();
        
        REQUIRE(path_name == graph.get_path_name(p));
        REQUIRE(path_rev == false);
        REQUIRE(path_pos == 0);
        
        // normalize the representation
        merge_non_branching_subpaths(surjected);
        
        REQUIRE(surjected.sequence() == mp_aln.sequence());
        REQUIRE(surjected.subpath_size() == 2);
        REQUIRE(surjected.subpath(0).score() == 1 + 6 + 2 - 6 - 32);
        REQUIRE(surjected.subpath(0).path().mapping_size() == 4);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h1));
        REQUIRE(surjected.subpath(0).path().mapping(0).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(1).position().node_id() == graph.get_id(h3));
        REQUIRE(surjected.subpath(0).path().mapping(1).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(1).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(2).position().node_id() == graph.get_id(h4));
        REQUIRE(surjected.subpath(0).path().mapping(2).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(2).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(3).position().node_id() == graph.get_id(h5));
        REQUIRE(surjected.subpath(0).path().mapping(3).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(3).position().offset() == 0);
        REQUIRE(surjected.subpath(0).next_size() == 0);
        REQUIRE(surjected.subpath(0).connection_size() == 1);
        REQUIRE(surjected.subpath(0).connection(0).next() == 1);
        REQUIRE(surjected.subpath(0).connection(0).score() == -1);
        REQUIRE(surjected.subpath(1).score() == 9);
        REQUIRE(surjected.subpath(1).path().mapping_size() == 1);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().node_id() == graph.get_id(h8));
        REQUIRE(surjected.subpath(1).path().mapping(0).position().is_reverse() == false);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().offset() == 1);
        REQUIRE(surjected.subpath(1).next_size() == 0);
        REQUIRE(surjected.subpath(1).connection_size() == 0);
    }
}

TEST_CASE("Duplicate path chunks can be detected", "[surject][multipath]") {
    
    Alignment aln;
    aln.set_sequence("ACGT");
    
    Surjector::path_chunk_t chunk1;
    Surjector::path_chunk_t chunk2;
    Surjector::path_chunk_t chunk3;
    Surjector::path_chunk_t chunk4;
    
    chunk1.first.first = aln.sequence().begin();
    chunk1.first.second = aln.sequence().begin() + 4;
    auto m00 = chunk1.second.add_mapping();
    m00->mutable_position()->set_node_id(1);
    m00->mutable_position()->set_is_reverse(false);
    m00->mutable_position()->set_offset(2);
    auto e00 = m00->add_edit();
    e00->set_from_length(2);
    e00->set_to_length(2);
    auto m01 = chunk1.second.add_mapping();
    m01->mutable_position()->set_node_id(2);
    m01->mutable_position()->set_is_reverse(false);
    m01->mutable_position()->set_offset(0);
    auto e01 = m01->add_edit();
    e01->set_from_length(2);
    e01->set_to_length(2);
    
    chunk2.first.first = aln.sequence().begin();
    chunk2.first.second = aln.sequence().begin() + 2;
    auto m10 = chunk2.second.add_mapping();
    m10->mutable_position()->set_node_id(1);
    m10->mutable_position()->set_is_reverse(false);
    m10->mutable_position()->set_offset(2);
    auto e10 = m10->add_edit();
    e10->set_from_length(2);
    e10->set_to_length(2);
    
    chunk3.first.first = aln.sequence().begin() + 2;
    chunk3.first.second = aln.sequence().begin() + 4;
    auto m20 = chunk3.second.add_mapping();
    m20->mutable_position()->set_node_id(2);
    m20->mutable_position()->set_is_reverse(false);
    m20->mutable_position()->set_offset(0);
    auto e20 = m20->add_edit();
    e20->set_from_length(2);
    e20->set_to_length(2);
    
    chunk4.first.first = aln.sequence().begin();
    chunk4.first.second = aln.sequence().begin() + 4;
    auto m30 = chunk4.second.add_mapping();
    m30->mutable_position()->set_node_id(1);
    m30->mutable_position()->set_is_reverse(false);
    m30->mutable_position()->set_offset(2);
    auto e30 = m30->add_edit();
    e30->set_from_length(2);
    e30->set_to_length(2);
    auto m31 = chunk4.second.add_mapping();
    m31->mutable_position()->set_node_id(3);
    m31->mutable_position()->set_is_reverse(false);
    m31->mutable_position()->set_offset(0);
    auto e31 = m31->add_edit();
    e31->set_from_length(2);
    e31->set_to_length(2);
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("AAAC");
    handle_t h2 = graph.create_handle("GTGT");
    handle_t h3 = graph.create_handle("GTAC");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h1);
    
    path_handle_t p = graph.create_path_handle("path");
    
    step_handle_t s1 = graph.append_step(p, h1);
    step_handle_t s2 = graph.append_step(p, h2);
    step_handle_t s3 = graph.append_step(p, h1);
    step_handle_t s4 = graph.append_step(p, h3);
    
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    
    vector<Surjector::path_chunk_t> path_chunks{chunk1, chunk2, chunk3, chunk4};
    vector<pair<step_handle_t, step_handle_t>> ref_chunks;
    ref_chunks.emplace_back(s1, s2);
    ref_chunks.emplace_back(s1, s1);
    ref_chunks.emplace_back(s2, s2);
    ref_chunks.emplace_back(s3, s4);
    
    vector<tuple<size_t, size_t, int32_t>> connections;
    
    surjector.filter_redundant_path_chunks(false, path_chunks, ref_chunks, connections);
    
    REQUIRE(ref_chunks.size() == path_chunks.size());
    REQUIRE(path_chunks.size() == 2);
    
}

TEST_CASE("Supplementary alignments can be generated", "[surject]") {
    
    bdsg::HashGraph graph;

    path_handle_t p = graph.create_path_handle("p");

    handle_t h1 = graph.create_handle("GTCGT");
    graph.append_step(p, h1);

    handle_t prev = h1;
    for (size_t i = 0; i < 20; ++i) {
        handle_t h = graph.create_handle(string(64, 'A'));
        graph.create_edge(prev, h);
        graph.append_step(p, h);
        prev = h;
    }
    handle_t h2 = graph.create_handle("TCCTTGC");
    graph.create_edge(prev, h2);
    graph.append_step(p, h2);
    
    handle_t h3 = graph.create_handle("TGTC");

    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h3, h2);

    bdsg::PositionOverlay pos_graph(&graph);
    Surjector surjector(&pos_graph);
    surjector.report_supplementary = true;

    unordered_set<path_handle_t> paths{p};
    
    SECTION("In a single-path alignment") {
        // with and without an unaligned middle portion
        vector<vector<handle_t>> read_paths{{h1, h2}, {h1, h3, h2}};
        for (const auto& read_path : read_paths) {
            Alignment read;
            string seq;
            Path* rpath = read.mutable_path();
            for (handle_t h : read_path) {
                Mapping* m = rpath->add_mapping();
                m->set_rank(rpath->mapping_size());
                m->mutable_position()->set_node_id(pos_graph.get_id(h));
                Edit* e = m->add_edit();
                e->set_from_length(pos_graph.get_length(h));
                e->set_to_length(pos_graph.get_length(h));
                
                seq += pos_graph.get_sequence(h);
            }
            read.set_sequence(seq);
            
            read.set_score(Aligner().scorer->score_contiguous_alignment(read));
            
            vector<Alignment> surjected_alns = surjector.surject(read, paths, true, false);
        
            REQUIRE(surjected_alns.size() == 2);
            bool found1 = false, found2 = false;
            for (auto& aln : surjected_alns) {
                const auto& path = aln.path();
                REQUIRE(path.mapping_size() == 1);
                const auto& mapping = path.mapping(0);
                REQUIRE(mapping.edit_size() == 2); // match and soft-clip
                if (mapping.position().node_id() == graph.get_id(h1)) {
                    found1 = true;
                    const auto& clip = mapping.edit(1);
                    REQUIRE(clip.from_length() == 0);
                    REQUIRE(clip.sequence() == aln.sequence().substr(mapping.edit(0).to_length(), string::npos));
                }
                else if (mapping.position().node_id() == graph.get_id(h2)) {
                    found2 = true;
                    const auto& clip = mapping.edit(0);
                    REQUIRE(clip.from_length() == 0);
                    REQUIRE(clip.sequence() == aln.sequence().substr(0, aln.sequence().size() - mapping.edit(1).to_length()));
                }
            }
        
            REQUIRE(found1);
            REQUIRE(found2);
        
            bool suppl1 = is_supplementary(surjected_alns.front());
            bool suppl2 = is_supplementary(surjected_alns.back());
            REQUIRE(suppl1 != suppl2);
        }
    }

    SECTION("In a multipath alignment") {

        vector<handle_t> read_path{h1, h2};
        multipath_alignment_t mp_aln;
        string seq;
        for (handle_t h : read_path) {
            if (mp_aln.subpath_size() != 0) {
                mp_aln.mutable_subpath(mp_aln.subpath_size() - 1)->add_next(mp_aln.subpath_size());
            }
            auto subpath = mp_aln.add_subpath();
            subpath->set_score(pos_graph.get_length(h));
            auto mapping = subpath->mutable_path()->add_mapping();
            auto pos = mapping->mutable_position();
            pos->set_node_id(pos_graph.get_id(h));
            pos->set_is_reverse(pos_graph.get_is_reverse(h));
            pos->set_offset(0);
            auto edit = mapping->add_edit();
            edit->set_from_length(pos_graph.get_length(h));
            edit->set_to_length(pos_graph.get_length(h));
            seq += pos_graph.get_sequence(h);
        }
        mp_aln.add_start(0);
        mp_aln.set_sequence(seq);
        
        {
            // unspliced
            vector<tuple<string, int64_t, bool>> positions;
            vector<multipath_alignment_t> surjected = surjector.surject(mp_aln, paths, positions, true, false);
    
            REQUIRE(surjected.size() == 2);
            bool found1 = false, found2 = false;
            for (auto& surj : surjected) {
                REQUIRE(surj.subpath_size() == 1);
                const auto& path = surj.subpath(0).path();
                REQUIRE(path.mapping_size() == 1);
                const auto& mapping = path.mapping(0);
                REQUIRE(mapping.edit_size() == 2);
                const auto& pos = mapping.position();
                if (pos.node_id() == pos_graph.get_id(h1)) {
                    found1 = true;
                }
                else if (pos.node_id() == pos_graph.get_id(h2)) {
                    found2 = true;
                }
            }

            REQUIRE(found1);
            REQUIRE(found2);
        
            bool suppl1 = is_supplementary(surjected.front());
            bool suppl2 = is_supplementary(surjected.back());
            REQUIRE(suppl1 != suppl2);
        }

        {
            // spliced are not returned as supplementaries
            vector<tuple<string, int64_t, bool>> positions;
            vector<multipath_alignment_t> surjected = surjector.surject(mp_aln, paths, positions, true, true);
    
            REQUIRE(surjected.size() == 1);
            for (auto& surj : surjected) {
                REQUIRE(surj.subpath_size() == 2);
                for (const auto& subpath : surj.subpath()) {
                    REQUIRE(subpath.path().mapping_size() == 1);
                }
                REQUIRE(surj.subpath(0).path().mapping(0).position().node_id() == pos_graph.get_id(h1));
                REQUIRE(surj.subpath(1).path().mapping(0).position().node_id() == pos_graph.get_id(h2));
            }
        }
    }
}

TEST_CASE("Surjector promotes a mapped secondary when the primary fails to surject", "[surject][promote]") {

    // A minimal graph is enough; promote_secondary_if_primary_unmapped only
    // inspects the already-surjected alignments, it does not touch the graph.
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("ACGT");
    path_handle_t p = graph.create_path_handle("ref");
    graph.append_step(p, h1);
    bdsg::PositionOverlay pos_graph(&graph);

    Surjector surjector(&pos_graph);

    // Helper to build a surjected alignment. If mapped, give it one mapping so
    // path().mapping_size() > 0; if unmapped, leave the path empty.
    auto make_aln = [&](const string& name, bool secondary, bool mapped, int32_t score,
                        int32_t mapq, bool supplementary) {
        Alignment aln;
        aln.set_name(name);
        aln.set_sequence("ACGT");
        aln.set_is_secondary(secondary);
        aln.set_score(score);
        aln.set_mapping_quality(mapq);
        if (mapped) {
            auto* m = aln.mutable_path()->add_mapping();
            m->mutable_position()->set_node_id(graph.get_id(h1));
        }
        if (supplementary) {
            set_annotation<bool>(aln, "supplementary", true);
        }
        return aln;
    };

    SECTION("Disabled by default: no promotion even if primary is unmapped") {
        vector<Alignment> group;
        group.push_back(make_aln("r", false, false, 0, 0, false));   // unmapped primary
        group.push_back(make_aln("r", true, true, 30, 40, false));   // mapped secondary
        REQUIRE(surjector.promote_secondary_on_failed_surjection == false);
        bool promoted = surjector.promote_secondary_if_primary_unmapped(group);
        REQUIRE(promoted == false);
        // Group is unchanged: primary still unmapped, secondary still present.
        REQUIRE(group.size() == 2);
        REQUIRE(group[0].path().mapping_size() == 0);
        REQUIRE(group[0].is_secondary() == false);
    }

    surjector.promote_secondary_on_failed_surjection = true;

    SECTION("Promotes the single mapped secondary into the primary slot") {
        vector<Alignment> group;
        group.push_back(make_aln("r", false, false, 0, 0, false));   // unmapped primary
        group.push_back(make_aln("r", true, true, 30, 40, false));   // mapped secondary
        bool promoted = surjector.promote_secondary_if_primary_unmapped(group);
        REQUIRE(promoted == true);
        // Exactly one alignment remains: the promoted secondary as new primary.
        REQUIRE(group.size() == 1);
        REQUIRE(group[0].is_secondary() == false);
        REQUIRE(group[0].path().mapping_size() > 0);
        REQUIRE(group[0].score() == 30);
        REQUIRE(get_annotation<bool>(group[0], "promoted_from_secondary") == true);
    }

    SECTION("Promotes the highest-scoring mapped secondary") {
        vector<Alignment> group;
        group.push_back(make_aln("r", false, false, 0, 0, false));   // unmapped primary
        group.push_back(make_aln("r", true, true, 20, 30, false));   // lower-scoring secondary
        group.push_back(make_aln("r", true, true, 45, 35, false));   // higher-scoring secondary
        bool promoted = surjector.promote_secondary_if_primary_unmapped(group);
        REQUIRE(promoted == true);
        REQUIRE(group.front().is_secondary() == false);
        REQUIRE(group.front().score() == 45);
        // The other mapped secondary is retained (still marked secondary).
        REQUIRE(group.size() == 2);
        REQUIRE(group[1].is_secondary() == true);
        REQUIRE(group[1].score() == 20);
    }

    SECTION("No-op when the primary surjected successfully") {
        vector<Alignment> group;
        group.push_back(make_aln("r", false, true, 50, 60, false));  // mapped primary
        group.push_back(make_aln("r", true, true, 30, 40, false));   // mapped secondary
        bool promoted = surjector.promote_secondary_if_primary_unmapped(group);
        REQUIRE(promoted == false);
        REQUIRE(group.size() == 2);
        REQUIRE(group[0].is_secondary() == false);
        REQUIRE(group[0].score() == 50);
    }

    SECTION("No promotion when no secondary could be surjected") {
        vector<Alignment> group;
        group.push_back(make_aln("r", false, false, 0, 0, false));   // unmapped primary
        group.push_back(make_aln("r", true, false, 0, 0, false));    // unmapped secondary
        bool promoted = surjector.promote_secondary_if_primary_unmapped(group);
        REQUIRE(promoted == false);
        REQUIRE(group.size() == 2);
        REQUIRE(group[0].is_secondary() == false);
        REQUIRE(group[0].path().mapping_size() == 0);
    }

    SECTION("Supplementary alignments are never promoted") {
        vector<Alignment> group;
        group.push_back(make_aln("r", false, false, 0, 0, false));   // unmapped primary
        group.push_back(make_aln("r", true, true, 90, 60, true));    // mapped but supplementary
        bool promoted = surjector.promote_secondary_if_primary_unmapped(group);
        REQUIRE(promoted == false);
        REQUIRE(group.size() == 2);
        REQUIRE(group[0].is_secondary() == false);
        REQUIRE(group[0].path().mapping_size() == 0);
    }
}

TEST_CASE("Surjector promotes a secondary pair when the primary pair fails to surject",
          "[surject][promote]") {

    // As with the single-end promotion test, the graph is only needed to
    // construct a Surjector; promote_secondary_pair_if_primary_unmapped inspects
    // the already-surjected alignments and never touches the graph.
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("ACGT");
    path_handle_t p = graph.create_path_handle("ref");
    graph.append_step(p, h1);
    bdsg::PositionOverlay pos_graph(&graph);

    Surjector surjector(&pos_graph);

    // Build one mate of a surjected pair. Mapped mates get a single mapping so
    // path().mapping_size() > 0; unmapped mates leave the path empty.
    auto make_aln = [&](const string& name, bool secondary, bool mapped, int32_t score,
                        int32_t mapq) {
        Alignment aln;
        aln.set_name(name);
        aln.set_sequence("ACGT");
        aln.set_is_secondary(secondary);
        aln.set_score(score);
        aln.set_mapping_quality(mapq);
        if (mapped) {
            auto* m = aln.mutable_path()->add_mapping();
            m->mutable_position()->set_node_id(graph.get_id(h1));
        }
        return aln;
    };

    SECTION("Disabled by default: no promotion even if the primary pair is unmapped") {
        vector<Alignment> mate1, mate2;
        mate1.push_back(make_aln("r", false, false, 0, 0));   // unmapped primary mate 1
        mate2.push_back(make_aln("r", false, false, 0, 0));   // unmapped primary mate 2
        mate1.push_back(make_aln("r", true, true, 30, 40));   // mapped secondary mate 1
        mate2.push_back(make_aln("r", true, true, 30, 40));   // mapped secondary mate 2
        REQUIRE(surjector.promote_secondary_on_failed_surjection == false);
        bool promoted = surjector.promote_secondary_pair_if_primary_unmapped(mate1, mate2);
        REQUIRE(promoted == false);
        REQUIRE(mate1.size() == 2);
        REQUIRE(mate2.size() == 2);
        REQUIRE(mate1[0].path().mapping_size() == 0);
        REQUIRE(mate2[0].path().mapping_size() == 0);
        REQUIRE(mate1[0].is_secondary() == false);
        REQUIRE(mate2[0].is_secondary() == false);
    }

    surjector.promote_secondary_on_failed_surjection = true;

    SECTION("Promotes a fully-mapped secondary pair into the primary slot") {
        vector<Alignment> mate1, mate2;
        mate1.push_back(make_aln("r", false, false, 0, 0));   // unmapped primary pair
        mate2.push_back(make_aln("r", false, false, 0, 0));
        mate1.push_back(make_aln("r", true, true, 30, 40));   // mapped secondary pair
        mate2.push_back(make_aln("r", true, true, 25, 35));
        bool promoted = surjector.promote_secondary_pair_if_primary_unmapped(mate1, mate2);
        REQUIRE(promoted == true);
        // Vectors stay the same length and index-aligned.
        REQUIRE(mate1.size() == 2);
        REQUIRE(mate2.size() == 2);
        // The promoted pair is now the primary (index 0) for both mates.
        REQUIRE(mate1[0].is_secondary() == false);
        REQUIRE(mate2[0].is_secondary() == false);
        REQUIRE(mate1[0].path().mapping_size() > 0);
        REQUIRE(mate2[0].path().mapping_size() > 0);
        REQUIRE(mate1[0].score() == 30);
        REQUIRE(mate2[0].score() == 25);
        REQUIRE(get_annotation<bool>(mate1[0], "promoted_from_secondary") == true);
        REQUIRE(get_annotation<bool>(mate2[0], "promoted_from_secondary") == true);
        // The demoted former-primary pair is now secondary.
        REQUIRE(mate1[1].is_secondary() == true);
        REQUIRE(mate2[1].is_secondary() == true);
    }

    SECTION("Promotes the secondary pair with the highest summed mapped-mate score") {
        vector<Alignment> mate1, mate2;
        mate1.push_back(make_aln("r", false, false, 0, 0));   // unmapped primary pair
        mate2.push_back(make_aln("r", false, false, 0, 0));
        mate1.push_back(make_aln("r", true, true, 20, 30));   // secondary pair, summed 40
        mate2.push_back(make_aln("r", true, true, 20, 30));
        mate1.push_back(make_aln("r", true, true, 45, 35));   // secondary pair, summed 80
        mate2.push_back(make_aln("r", true, true, 35, 25));
        bool promoted = surjector.promote_secondary_pair_if_primary_unmapped(mate1, mate2);
        REQUIRE(promoted == true);
        REQUIRE(mate1[0].is_secondary() == false);
        REQUIRE(mate2[0].is_secondary() == false);
        // The higher summed-score pair (45 + 35 = 80) wins.
        REQUIRE(mate1[0].score() == 45);
        REQUIRE(mate2[0].score() == 35);
    }

    SECTION("Prefers a fully-mapped pair over a half-mapped pair with equal summed score") {
        vector<Alignment> mate1, mate2;
        mate1.push_back(make_aln("r", false, false, 0, 0));   // unmapped primary pair
        mate2.push_back(make_aln("r", false, false, 0, 0));
        // Half-mapped secondary: only mate 1 mapped, summed score 50.
        mate1.push_back(make_aln("r", true, true, 50, 40));
        mate2.push_back(make_aln("r", true, false, 0, 0));
        // Fully-mapped secondary: both mapped, summed score 50 (30 + 20).
        mate1.push_back(make_aln("r", true, true, 30, 20));
        mate2.push_back(make_aln("r", true, true, 20, 20));
        bool promoted = surjector.promote_secondary_pair_if_primary_unmapped(mate1, mate2);
        REQUIRE(promoted == true);
        // The fully-mapped pair is chosen despite the tie in summed score.
        REQUIRE(mate1[0].path().mapping_size() > 0);
        REQUIRE(mate2[0].path().mapping_size() > 0);
        REQUIRE(mate1[0].score() == 30);
        REQUIRE(mate2[0].score() == 20);
    }

    SECTION("Promotes a half-mapped secondary pair when no fully-mapped pair exists") {
        vector<Alignment> mate1, mate2;
        mate1.push_back(make_aln("r", false, false, 0, 0));   // unmapped primary pair
        mate2.push_back(make_aln("r", false, false, 0, 0));
        // Only mate 2 of the secondary surjected.
        mate1.push_back(make_aln("r", true, false, 0, 0));
        mate2.push_back(make_aln("r", true, true, 40, 30));
        bool promoted = surjector.promote_secondary_pair_if_primary_unmapped(mate1, mate2);
        REQUIRE(promoted == true);
        REQUIRE(mate1[0].is_secondary() == false);
        REQUIRE(mate2[0].is_secondary() == false);
        // The mapped mate carries the alignment; its partner stays unmapped.
        REQUIRE(mate1[0].path().mapping_size() == 0);
        REQUIRE(mate2[0].path().mapping_size() > 0);
        REQUIRE(mate2[0].score() == 40);
    }

    SECTION("No-op when either mate of the primary pair surjected") {
        vector<Alignment> mate1, mate2;
        // Primary pair is half-mapped: a valid SAM state, left untouched.
        mate1.push_back(make_aln("r", false, true, 50, 60));
        mate2.push_back(make_aln("r", false, false, 0, 0));
        mate1.push_back(make_aln("r", true, true, 30, 40));
        mate2.push_back(make_aln("r", true, true, 30, 40));
        bool promoted = surjector.promote_secondary_pair_if_primary_unmapped(mate1, mate2);
        REQUIRE(promoted == false);
        REQUIRE(mate1[0].is_secondary() == false);
        REQUIRE(mate1[0].score() == 50);
    }

    SECTION("No promotion when no secondary pair has a mapped mate") {
        vector<Alignment> mate1, mate2;
        mate1.push_back(make_aln("r", false, false, 0, 0));   // unmapped primary pair
        mate2.push_back(make_aln("r", false, false, 0, 0));
        mate1.push_back(make_aln("r", true, false, 0, 0));    // unmapped secondary pair
        mate2.push_back(make_aln("r", true, false, 0, 0));
        bool promoted = surjector.promote_secondary_pair_if_primary_unmapped(mate1, mate2);
        REQUIRE(promoted == false);
        REQUIRE(mate1[0].is_secondary() == false);
        REQUIRE(mate2[0].is_secondary() == false);
        REQUIRE(mate1[0].path().mapping_size() == 0);
        REQUIRE(mate2[0].path().mapping_size() == 0);
    }

    SECTION("No-op when the mate vectors are not index-aligned") {
        vector<Alignment> mate1, mate2;
        mate1.push_back(make_aln("r", false, false, 0, 0));   // unmapped primary pair
        mate2.push_back(make_aln("r", false, false, 0, 0));
        mate1.push_back(make_aln("r", true, true, 30, 40));   // secondary only present for mate 1
        bool promoted = surjector.promote_secondary_pair_if_primary_unmapped(mate1, mate2);
        REQUIRE(promoted == false);
        REQUIRE(mate1.size() == 2);
        REQUIRE(mate2.size() == 1);
    }
}

}
}
