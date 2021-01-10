/// \file deletion_aligner.cpp
///  
/// unit tests for the DeletionAligner
///

#include <iostream>
#include "path.hpp"
#include "deletion_aligner.hpp"
#include "random_graph.hpp"
#include "catch.hpp"

#include "bdsg/hash_graph.hpp"

namespace vg {
namespace unittest {

TEST_CASE("Deletion aligner finds optimal deletions", "[aligner][deletionaligner]") {

    bdsg::HashGraph graph;
    
    auto check_aln = [&](const Alignment& aln,
                         vector<handle_t> correct_path) {
        REQUIRE(aln.path().mapping_size() == correct_path.size());
        int total_length = 0;
        for (size_t i = 0; i < correct_path.size(); ++i) {
            auto h = correct_path[i];
            const auto& m = aln.path().mapping(i);
            REQUIRE(m.position().node_id() == graph.get_id(h));
            REQUIRE(m.position().is_reverse() == graph.get_is_reverse(h));
            REQUIRE(m.position().offset() == 0);
            REQUIRE(mapping_from_length(m) == graph.get_length(h));
            REQUIRE(mapping_to_length(m) == 0);
            total_length += graph.get_length(h);
        }
        REQUIRE(aln.score() == (total_length ? -total_length - 5 : 0));
    };
    
    // scoring bubbles with powers of two makes ordered
    // scores act like binary bits in counting order
    auto h1 = graph.create_handle("AA");// +/- 1
    auto h2 = graph.create_handle("A");
    auto h3 = graph.create_handle("AAA");
    auto h4 = graph.create_handle("A");
    auto h5 = graph.create_handle("AAA");// +/- 2
    auto h6 = graph.create_handle("A");
    auto h7 = graph.create_handle("AAAA"); // +/- 4
    auto h8 = graph.create_handle("AA");
    auto h9 = graph.create_handle("A");
    auto h10 = graph.create_handle("AAAAAAAAA"); // +/- 8
    
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h3, h5);
    graph.create_edge(h4, h6);
    graph.create_edge(h5, h6);
    graph.create_edge(h6, h7);
    graph.create_edge(h6, h8);
    graph.create_edge(h7, h8);
    graph.create_edge(h8, h9);
    graph.create_edge(h8, h10);
    
    DeletionAligner aligner(6, 1);
    
    Alignment aln;
    
    SECTION("Single traceback works") {
        aligner.align(aln, graph);
        check_aln(aln, {h2, h3, h4, h6, h8, h9});
    }
    
    SECTION("Multi traceback works") {
        
        int num_alts = 15;
        vector<Alignment> alts;
        aligner.align_multi(aln, alts, graph, num_alts);
        
        
        vector<vector<handle_t>> corrects{
            {h2, h3, h4, h6, h8, h9},
            {h1, h3, h4, h6, h8, h9},
            {h2, h3, h5, h6, h8, h9},
            {h1, h3, h5, h6, h8, h9},
            {h2, h3, h4, h6, h7, h8, h9},
            {h1, h3, h4, h6, h7, h8, h9},
            {h2, h3, h5, h6, h7, h8, h9},
            {h1, h3, h5, h6, h7, h8, h9},
            {h2, h3, h4, h6, h8, h10},
            {h1, h3, h4, h6, h8, h10},
            {h2, h3, h5, h6, h8, h10},
            {h1, h3, h5, h6, h8, h10},
            {h2, h3, h4, h6, h7, h8, h10},
            {h1, h3, h4, h6, h7, h8, h10},
            {h2, h3, h5, h6, h7, h8, h10}
        };
        
        REQUIRE(corrects.size() == num_alts);
        REQUIRE(alts.size() == num_alts);
        for (size_t i = 0; i < num_alts; ++i) {
            check_aln(alts[i], corrects[i]);
        }
        
    }
}

}
}
