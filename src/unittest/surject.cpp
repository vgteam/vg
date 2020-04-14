///
/// \file position.cpp
///  
/// Unit tests for Position and pos_t manipulation
///

#include "catch.hpp"
#include "surjector.hpp"
#include "aligner.hpp"

#include "bdsg/hash_graph.hpp"
#include "bdsg/overlays/path_position_overlays.hpp"
#include <vg/vg.pb.h>

namespace vg {
namespace unittest {

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
    
    read.set_score(Aligner().score_ungapped_alignment(read));
    
    set<string> path_names{pos_graph.get_path_name(p)};
    Alignment surjected = surjector.surject(read, path_names, true, true);
    
    vector<handle_t> surjected_path{h1, h3, h4, h6};
    
    REQUIRE(surjected.path().mapping_size() == read_path.size());
    
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(surjected.path().mapping(i).position().node_id() == graph.get_id(surjected_path[i]));
        REQUIRE(!surjected.path().mapping(i).position().is_reverse());
    }
    
    // should not penalize long deletions (assumed to be splices)
    REQUIRE(surjected.score() == read.score() - Aligner().mismatch - Aligner().match);
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
    
    rev_read.set_score(Aligner().score_ungapped_alignment(rev_read));
    
    Alignment rev_surjected = surjector.surject(rev_read, path_names, true, true);
    
    REQUIRE(rev_surjected.path().mapping_size() == read_path.size());
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(rev_surjected.path().mapping(i).position().node_id()
                == graph.get_id(surjected_path[surjected_path.size() - i - 1]));
        REQUIRE(rev_surjected.path().mapping(i).position().is_reverse());
    }
    
    // should not penalize long deletions (assumed to be splices)
    REQUIRE(rev_surjected.score() == read.score() - Aligner().mismatch - Aligner().match);
    REQUIRE(rev_surjected.refpos_size() == 1);
    REQUIRE(rev_surjected.refpos(0).name() == graph.get_path_name(p));
    REQUIRE(rev_surjected.refpos(0).offset() == 0);
    
}

}
}
