///
/// \file masker.cpp
///
/// Unit tests for the Masker class
///

#include "catch.hpp"
#include "../masker.hpp"
#include "../handle.hpp"
#include "../gbwt_helper.hpp"
#include "../gbwtgraph_helper.hpp"

#include <bdsg/hash_graph.hpp>
#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/utils.h>
#include <gbwtgraph/path_cover.h>

namespace vg {
namespace unittest {

TEST_CASE( "Masker produces expected results on simple paths", "[masker]" ) {
    
    SECTION( "Simple masking task") {
        
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("ACGT");
        handle_t h2 = graph.create_handle("TTA");
        handle_t h3 = graph.create_handle("GTAC");
        handle_t h4 = graph.create_handle("AA");
        handle_t h5 = graph.create_handle("G");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        
        path_handle_t p = graph.create_path_handle("p");
        
        graph.append_step(p, h1);
        graph.append_step(p, h2);
        graph.append_step(p, h3);
        graph.append_step(p, h4);
        graph.append_step(p, h5);
        
        Masker masker(graph);
        
        std::vector<std::tuple<std::string, size_t, size_t>> intervals{
            {"p", 4, 13}
        };
        
        masker.mask_sequences(intervals);
        
        REQUIRE(graph.get_sequence(h1) == "ACGT");
        REQUIRE(graph.get_sequence(h2) == "NNN");
        REQUIRE(graph.get_sequence(h3) == "NNNN");
        REQUIRE(graph.get_sequence(h4) == "NN");
        REQUIRE(graph.get_sequence(h5) == "G");
    }
    
    SECTION( "Mask on reverse strand") {
        
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("ACGT");
        handle_t h2 = graph.create_handle("TTA");
        handle_t h3 = graph.create_handle("GTAC");
        handle_t h4 = graph.create_handle("AA");
        handle_t h5 = graph.create_handle("G");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        
        path_handle_t p = graph.create_path_handle("p");
        
        graph.append_step(p, graph.flip(h5));
        graph.append_step(p, graph.flip(h4));
        graph.append_step(p, graph.flip(h3));
        graph.append_step(p, graph.flip(h2));
        graph.append_step(p, graph.flip(h1));
        
        Masker masker(graph);
        
        std::vector<std::tuple<std::string, size_t, size_t>> intervals{
            {"p", 4, 11}
        };
        
        masker.mask_sequences(intervals);
        
        REQUIRE(graph.get_sequence(h1) == "ACGN");
        REQUIRE(graph.get_sequence(h2) == "NNN");
        REQUIRE(graph.get_sequence(h3) == "NNNC");
        REQUIRE(graph.get_sequence(h4) == "AA");
        REQUIRE(graph.get_sequence(h5) == "G");
    }

    SECTION( "Mask on a single node") {
        
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("ACGTACGT");
        
        path_handle_t p = graph.create_path_handle("p");
        
        graph.append_step(p, h1);
        
        Masker masker(graph);
        
        std::vector<std::tuple<std::string, size_t, size_t>> intervals{
            {"p", 1, 2},
            {"p", 5, 7}
        };
        
        masker.mask_sequences(intervals);
        
        REQUIRE(graph.get_sequence(h1) == "ANGTANNT");
    }
    
    SECTION( "Mask across a clipped path" ) {
        
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("ACGTACGT");
        handle_t h2 = graph.create_handle("ATCTCTAA");
        handle_t h3 = graph.create_handle("G");
        handle_t h4 = graph.create_handle("GGACTCA");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        
        std::string pname1 = graph.create_path_name(PathSense::GENERIC,
                                                    PathMetadata::NO_SAMPLE_NAME,
                                                    "p",
                                                    PathMetadata::NO_HAPLOTYPE,
                                                    PathMetadata::NO_PHASE_BLOCK,
                                                    handlegraph::subrange_t(0, PathMetadata::NO_END_POSITION));
        std::string pname2 = graph.create_path_name(PathSense::GENERIC,
                                                    PathMetadata::NO_SAMPLE_NAME,
                                                    "p",
                                                    PathMetadata::NO_HAPLOTYPE,
                                                    PathMetadata::NO_PHASE_BLOCK,
                                                    handlegraph::subrange_t(17, PathMetadata::NO_END_POSITION));
        
        path_handle_t p1 = graph.create_path_handle(pname1);
        path_handle_t p2 = graph.create_path_handle(pname2);
        
        graph.append_step(p1, h1);
        graph.append_step(p1, h2);
        
        graph.append_step(p2, h4);
        
        Masker masker(graph);
        
        std::vector<std::tuple<std::string, size_t, size_t>> intervals{
            {"p", 12, 20},
        };
        
        masker.mask_sequences(intervals);
        
        REQUIRE(graph.get_sequence(h1) == "ACGTACGT");
        REQUIRE(graph.get_sequence(h2) == "ATCTNNNN");
        REQUIRE(graph.get_sequence(h3) == "G");
        REQUIRE(graph.get_sequence(h4) == "NNNCTCA");
    }
    
    SECTION( "Masker works on GBWTGraph" ) {
        
        std::vector<gbwt::vector_type> paths;
        
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("ACGT");
        handle_t h2 = graph.create_handle("TTA");
        handle_t h3 = graph.create_handle("GTAC");
        handle_t h4 = graph.create_handle("AA");
        handle_t h5 = graph.create_handle("G");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        
        std::string ref = graph.create_path_name(PathSense::GENERIC,
                                                 PathMetadata::NO_SAMPLE_NAME,
                                                 "chr",
                                                 PathMetadata::NO_HAPLOTYPE,
                                                 PathMetadata::NO_PHASE_BLOCK,
                                                 PathMetadata::NO_SUBRANGE);
        
        auto p = graph.create_path_handle(ref);
        
        graph.append_step(p, h1);
        graph.append_step(p, h2);
        graph.append_step(p, h3);
        graph.append_step(p, h4);
        graph.append_step(p, h5);
        
        // parameters surmised from from get_gbwt() source code
        gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        gbwt::GBWTBuilder gbwt_builder(8, 20);
                
        gbwtgraph::MetadataBuilder metadata;
        // stolen from gbwtgraph::GFAParsingParameters::DEFAULT_XXXXX
        metadata.add_path_name_format(".*", "C", PathSense::GENERIC);
        
        auto jobs = gbwtgraph::gbwt_construction_jobs(graph, 10000);
        std::vector<std::vector<path_handle_t>> assignments = assign_paths(graph, jobs, &metadata, nullptr);
        assert(assignments.size() == 1);
        assert(assignments.front().size() == 1);
        gbwtgraph::insert_paths(graph, assignments.front(), gbwt_builder, 0, false);
        gbwt_builder.index.addMetadata();
        gbwt_builder.index.metadata = metadata.get_metadata();
                
        gbwt_builder.finish();
        gbwt::GBWT gbwt_index(gbwt_builder.index);
        
        gbwtgraph::GBWTGraph gbwt_graph(gbwt_index, graph);
        
        Masker masker(gbwt_graph);
        
        std::vector<std::tuple<std::string, size_t, size_t>> intervals{
            {ref, 4, 13}
        };
        
        masker.mask_sequences(intervals);
        
        REQUIRE(gbwt_graph.get_sequence(h1) == "ACGT");
        REQUIRE(gbwt_graph.get_sequence(h2) == "NNN");
        REQUIRE(gbwt_graph.get_sequence(h3) == "NNNN");
        REQUIRE(gbwt_graph.get_sequence(h4) == "NN");
        REQUIRE(gbwt_graph.get_sequence(h5) == "G");
        REQUIRE(gbwt_graph.get_sequence(graph.flip(h1)) == "ACGT");
        REQUIRE(gbwt_graph.get_sequence(graph.flip(h2)) == "NNN");
        REQUIRE(gbwt_graph.get_sequence(graph.flip(h3)) == "NNNN");
        REQUIRE(gbwt_graph.get_sequence(graph.flip(h4)) == "NN");
        REQUIRE(gbwt_graph.get_sequence(graph.flip(h5)) == "C");
    }
}

TEST_CASE( "Masker produces expected results on paths with variation", "[masker]" ) {
    
    SECTION( "Simple bubble" ) {
        
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("ACGT");
        handle_t h2 = graph.create_handle("TTA");
        handle_t h3 = graph.create_handle("GTAC");
        handle_t h4 = graph.create_handle("AA");
        handle_t h5 = graph.create_handle("G");
        handle_t h6 = graph.create_handle("CG");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h2, h6);
        graph.create_edge(h3, h4);
        graph.create_edge(h6, h4);
        graph.create_edge(h4, h5);
        
        path_handle_t p = graph.create_path_handle("p");
        
        graph.append_step(p, h1);
        graph.append_step(p, h2);
        graph.append_step(p, h3);
        graph.append_step(p, h4);
        graph.append_step(p, h5);
        
        Masker masker(graph);
        
        std::vector<std::tuple<std::string, size_t, size_t>> intervals{
            {"p", 4, 13}
        };
        
        masker.mask_sequences(intervals);
        
        REQUIRE(graph.get_sequence(h1) == "ACGT");
        REQUIRE(graph.get_sequence(h2) == "NNN");
        REQUIRE(graph.get_sequence(h3) == "NNNN");
        REQUIRE(graph.get_sequence(h4) == "NN");
        REQUIRE(graph.get_sequence(h5) == "G");
        REQUIRE(graph.get_sequence(h6) == "NN");
    }

    SECTION( "Nested bubbles" ) {
        
        bdsg::HashGraph graph;
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");
        handle_t h6 = graph.create_handle("A");
        handle_t h7 = graph.create_handle("A");
        handle_t h8 = graph.create_handle("A");
        handle_t h9 = graph.create_handle("A");
        handle_t h10 = graph.create_handle("A");
        handle_t h11 = graph.create_handle("A");
        handle_t h12 = graph.create_handle("A");
        handle_t h13 = graph.create_handle("A");
        handle_t h14 = graph.create_handle("A");
        handle_t h15 = graph.create_handle("A");
        handle_t h16 = graph.create_handle("A");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h1, h3);
        graph.create_edge(h2, h7);
        graph.create_edge(h3, h4);
        graph.create_edge(h3, h5);
        graph.create_edge(h4, h6);
        graph.create_edge(h5, h6);
        graph.create_edge(h6, h7);
        graph.create_edge(h7, h8);
        graph.create_edge(h7, h9);
        graph.create_edge(h8, h16);
        graph.create_edge(h9, h10);
        graph.create_edge(h9, h11);
        graph.create_edge(h10, h12);
        graph.create_edge(h11, h12);
        graph.create_edge(h12, h13);
        graph.create_edge(h12, h14);
        graph.create_edge(h13, h15);
        graph.create_edge(h14, h15);
        graph.create_edge(h15, h16);
        
        // some extra bubbles to make sure we get the right snarl decomposition
        handle_t h17 = graph.create_handle("A");
        handle_t h18 = graph.create_handle("A");
        handle_t h19 = graph.create_handle("A");
        handle_t h20 = graph.create_handle("A");
        handle_t h21 = graph.create_handle("A");
        handle_t h22 = graph.create_handle("A");
        handle_t h23 = graph.create_handle("A");
        handle_t h24 = graph.create_handle("A");
        handle_t h25 = graph.create_handle("A");
        
        graph.create_edge(h16, h17);
        graph.create_edge(h16, h18);
        graph.create_edge(h17, h19);
        graph.create_edge(h18, h19);
        graph.create_edge(h19, h20);
        graph.create_edge(h19, h21);
        graph.create_edge(h20, h22);
        graph.create_edge(h21, h22);
        graph.create_edge(h22, h23);
        graph.create_edge(h22, h24);
        graph.create_edge(h23, h25);
        graph.create_edge(h24, h25);
        
        path_handle_t p = graph.create_path_handle("p");
        
        graph.append_step(p, h1);
        graph.append_step(p, h3);
        graph.append_step(p, h4);
        graph.append_step(p, h6);
        graph.append_step(p, h7);
        graph.append_step(p, h9);
        graph.append_step(p, h10);
        graph.append_step(p, h12);
        graph.append_step(p, h13);
        graph.append_step(p, h15);
        graph.append_step(p, h16);
        
        Masker masker(graph);
        
        std::vector<std::tuple<std::string, size_t, size_t>> intervals{
            {"p", 1, 11}
        };
        
        masker.mask_sequences(intervals);
        
        REQUIRE(graph.get_sequence(h1) == "A");
        REQUIRE(graph.get_sequence(h2) == "A");
        REQUIRE(graph.get_sequence(h3) == "N");
        REQUIRE(graph.get_sequence(h4) == "N");
        REQUIRE(graph.get_sequence(h5) == "N");
        REQUIRE(graph.get_sequence(h6) == "N");
        REQUIRE(graph.get_sequence(h7) == "N");
        REQUIRE(graph.get_sequence(h8) == "N");
        REQUIRE(graph.get_sequence(h9) == "N");
        REQUIRE(graph.get_sequence(h10) == "N");
        REQUIRE(graph.get_sequence(h11) == "N");
        REQUIRE(graph.get_sequence(h12) == "N");
        REQUIRE(graph.get_sequence(h13) == "N");
        REQUIRE(graph.get_sequence(h14) == "N");
        REQUIRE(graph.get_sequence(h15) == "N");
        REQUIRE(graph.get_sequence(h16) == "N");
    }
}

}
}
