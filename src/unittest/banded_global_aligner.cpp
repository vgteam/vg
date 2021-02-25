//
//  banded_global_aligner.cpp
//  
//  Unit tests for banded global aligner module.
//

#include <stdio.h>
#include "catch.hpp"
#include "test_aligner.hpp"
#include "vg.hpp"
#include "path.hpp"
#include "banded_global_aligner.hpp"
#include "vg/io/json2pb.h"
#include "bdsg/hash_graph.hpp"

using namespace google::protobuf;
using namespace vg::io;
namespace vg {
    namespace unittest {
        
        TEST_CASE( "Banded global aligner produces correct alignments with all types of edits",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner produces correct alignment when read matches exactly") {
                
                VG graph;
                
                TestAligner aligner_source;
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
                
                int band_width = 1;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when read matches across a doubly-reversing edge") {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n1, n0, true, true);
                graph.create_edge(n2, n0, true, true);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AGTGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when there is a mismatch" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("CCCAGTT");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CCCAGTGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 6);
                REQUIRE(path.mapping(0).edit(0).to_length() == 6);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 1);
                REQUIRE(path.mapping(0).edit(1).to_length() == 1);
                REQUIRE(path.mapping(0).edit(1).sequence() == string("G"));
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when a there is a single base deletion" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("CCCAGATG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CCCAGTGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 5);
                REQUIRE(path.mapping(0).edit(0).to_length() == 5);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 1);
                REQUIRE(path.mapping(0).edit(1).to_length() == 0);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(2).from_length() == 2);
                REQUIRE(path.mapping(0).edit(2).to_length() == 2);
                REQUIRE(path.mapping(0).edit(2).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when a there is a multi-base deletion"  ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AACCCAGATG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AACCCAGGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);;
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 7);
                REQUIRE(path.mapping(0).edit(0).to_length() == 7);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 2);
                REQUIRE(path.mapping(0).edit(1).to_length() == 0);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(2).from_length() == 1);
                REQUIRE(path.mapping(0).edit(2).to_length() == 1);
                REQUIRE(path.mapping(0).edit(2).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when a there is a single base insertion")  {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AACCCAGAGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 7);
                REQUIRE(path.mapping(0).edit(0).to_length() == 7);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 0);
                REQUIRE(path.mapping(0).edit(1).to_length() == 1);
                REQUIRE(path.mapping(0).edit(1).sequence() == string("A"));
                
                REQUIRE(path.mapping(0).edit(2).from_length() == 1);
                REQUIRE(path.mapping(0).edit(2).to_length() == 1);
                REQUIRE(path.mapping(0).edit(2).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when a there is a multi-base insertion" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AACCCAGATGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 7);
                REQUIRE(path.mapping(0).edit(0).to_length() == 7);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 0);
                REQUIRE(path.mapping(0).edit(1).to_length() == 2);
                REQUIRE(path.mapping(0).edit(1).sequence() == string("AT"));
                
                REQUIRE(path.mapping(0).edit(2).from_length() == 1);
                REQUIRE(path.mapping(0).edit(2).to_length() == 1);
                REQUIRE(path.mapping(0).edit(2).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when a there is a deletion across a node boundary" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AACCCAGG");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("AT");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AACCCAGATGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 7);
                REQUIRE(path.mapping(0).edit(0).to_length() == 7);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 1);
                REQUIRE(path.mapping(0).edit(1).to_length() == 0);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 0);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(1).from_length() == 1);
                REQUIRE(path.mapping(1).edit(1).to_length() == 1);
                REQUIRE(path.mapping(1).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }

            
            SECTION( "Banded global aligner produces correct alignment when it begins with a single base insertion") {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AACCCAGG");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("AT");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("TAACCCAGGCATGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 0);
                REQUIRE(path.mapping(0).edit(0).to_length() == 1);
                REQUIRE(path.mapping(0).edit(0).sequence() == string("T"));
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 8);
                REQUIRE(path.mapping(0).edit(1).to_length() == 8);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 2);
                REQUIRE(path.mapping(1).edit(0).to_length() == 2);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
            }
            
            SECTION( "Banded global aligner produces correct alignment when it begins with a multi-base insertion" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AACCCAGG");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("AT");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("TGAACCCAGGCATGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 0);
                REQUIRE(path.mapping(0).edit(0).to_length() == 2);
                REQUIRE(path.mapping(0).edit(0).sequence() == string("TG"));
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 8);
                REQUIRE(path.mapping(0).edit(1).to_length() == 8);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 2);
                REQUIRE(path.mapping(1).edit(0).to_length() == 2);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when it begins with a single base deletion" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AACCCAGG");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("AT");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("ACCCAGGCATGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);;
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 1);
                REQUIRE(path.mapping(0).edit(0).to_length() == 0);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 7);
                REQUIRE(path.mapping(0).edit(1).to_length() == 7);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 2);
                REQUIRE(path.mapping(1).edit(0).to_length() == 2);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when it begins with a multi-base deletion" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AACCCAGG");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("AT");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CCCAGGCATGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 2;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);;
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 2);
                REQUIRE(path.mapping(0).edit(0).to_length() == 0);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 6);
                REQUIRE(path.mapping(0).edit(1).to_length() == 6);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 2);
                REQUIRE(path.mapping(1).edit(0).to_length() == 2);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when it begins with a deletion that crosses a node boundary" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AA");
                Node* n1 = graph.create_node("CCCAGGCA");
                Node* n2 = graph.create_node("GTGCTATA");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CCAGGCATGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 3;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 2);
                REQUIRE(path.mapping(0).edit(0).to_length() == 0);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 0);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(1).from_length() == 7);
                REQUIRE(path.mapping(1).edit(1).to_length() == 7);
                REQUIRE(path.mapping(1).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when it begins with a deletion that ends at a node boundary",
                    "[alignment][banded][mapping]" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("A");
                Node* n1 = graph.create_node("C");
                
                graph.create_edge(n0, n1);
                
                string read = string("C");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_padding = 1;
                
                aligner.align_global_banded(aln, graph, band_padding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 1);
                REQUIRE(path.mapping(0).edit(0).to_length() == 0);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
            }
            
            SECTION( "Banded global aligner produces correct alignment when it begins with a deletion that ends at a node boundary",
                    "[alignment][banded][mapping]" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AATG");
                Node* n1 = graph.create_node("C");
                
                graph.create_edge(n0, n1);
                
                string read = string("C");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_padding = 1;
                
                aligner.align_global_banded(aln, graph, band_padding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 0);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
            }
        }
        
        TEST_CASE( "Banded global aligner produces correct alignments with different graph structures",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner produces correct alignment when there are paths of different lengths in the graph" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("TG");
                Node* n1 = graph.create_node("TGGC");
                Node* n2 = graph.create_node("AAA");
                Node* n3 = graph.create_node("AGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("TGTGGCAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 2);
                REQUIRE(path.mapping(0).edit(0).to_length() == 2);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 4);
                REQUIRE(path.mapping(1).edit(0).to_length() == 4);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 3);
                REQUIRE(path.mapping(2).edit(0).to_length() == 3);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment against a singleton graph" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("C");
                
                string read = string("A");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_padding = 1;
                
                aligner.align_global_banded(aln, graph, band_padding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 1);
                REQUIRE(path.mapping(0).edit(0).to_length() == 1);
                REQUIRE(path.mapping(0).edit(0).sequence() == "A");
                
            }
            
            SECTION( "Banded global aligner can align to a graph that is not topologically sorted" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("CTAG");
                Node* n1 = graph.create_node("T");
                Node* n2 = graph.create_node("CC");
                Node* n3 = graph.create_node("GTA");
                
                graph.create_edge(n2, n0);
                graph.create_edge(n2, n1);
                graph.create_edge(n0, n3);
                graph.create_edge(n1, n3);

                
                string read = string("CCCTAGGTA");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_padding = 1;
                
                aligner.align_global_banded(aln, graph, band_padding);
                
                const Path& path = aln.path();
                                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(3).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 3);
                REQUIRE(path.mapping(1).position().node_id() == 1);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 2);
                REQUIRE(path.mapping(0).edit(0).to_length() == 2);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 4);
                REQUIRE(path.mapping(1).edit(0).to_length() == 4);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 3);
                REQUIRE(path.mapping(2).edit(0).to_length() == 3);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when a node is masked" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("TG");
                Node* n1 = graph.create_node("TGGC");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("AGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("TGTGGCAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 2);
                REQUIRE(path.mapping(0).edit(0).to_length() == 2);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 4);
                REQUIRE(path.mapping(1).edit(0).to_length() == 4);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 3);
                REQUIRE(path.mapping(2).edit(0).to_length() == 3);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when there are multiple source nodes" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("CGCC");
                Node* n2 = graph.create_node("ATTA");
                
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n2);
                
                string read = string("AGTGATTA");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 3);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 4);
                REQUIRE(path.mapping(1).edit(0).to_length() == 4);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when there are multiple sink nodes" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                
                string read = string("AGTGC");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when there are multiple source and multiple sink nodes" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("CGCC");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("TGGG");
                Node* n4 = graph.create_node("GTTA");
                
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                
                string read = string("AGTGTTGGG");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 3);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 4);
                REQUIRE(path.mapping(2).edit(0).to_length() == 4);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment on a separated graph" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("CGCC");
                Node* n2 = graph.create_node("GTTA");
                Node* n3 = graph.create_node("TGGG");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n2, n3);
                
                string read = string("AGTGCGCC");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 4);
                REQUIRE(path.mapping(1).edit(0).to_length() == 4);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
            }
            
        }
        
        TEST_CASE( "Banded global aligner produces correct alignments in edge cases of banding",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner produces correct alignment when the band is wider than the read" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("CAGGA");
                Node* n1 = graph.create_node("AA");
                Node* n2 = graph.create_node("GTGTATA");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CAGGAAATGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 20;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 5);
                REQUIRE(path.mapping(0).edit(0).to_length() == 5);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 2);
                REQUIRE(path.mapping(1).edit(0).to_length() == 2);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct alignment when the band is one base wide" ) {
                
                VG graph;
                
                
                
                TestAligner aligner_source;
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
                
                int band_width = 0;
                bool permissive_banding = false;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner does not have overflow errors when scores are in the 100s" ) {
                
                VG graph;
                
                // make the alignment take a long gap that would overflow an 8 bit integer
                Node* n0 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
                Node* n1 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
                Node* n2 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
                Node* n3 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
                Node* n4 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
                Node* n5 = graph.create_node("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n1, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                
                
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                string read = "C";
                Alignment aln;
                aln.set_sequence(read);
                
                int band_padding = 1;
                aligner.align_global_banded(aln, graph, band_padding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // score would have triggered overflow
                REQUIRE(aln.score() < numeric_limits<int8_t>::min());
                
                // mostly looking for it to not explode on this example, don't bother checking path
            }
        }
        
        TEST_CASE( "Banded global aligner produces correct alignments with permissive banding option",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner runs with permissive banding" ) {
                
                VG graph;
                
                TestAligner aligner_source;
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
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }

            SECTION( "Banded global aligner with permissive banding can complete alignment outside of the padded band" ) {
                
                VG graph;
                
                
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("CTGGTGTAGTA");
                Node* n2 = graph.create_node("AATTCCCCCGG");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AGTGAATTGGTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 0;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 3);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 4);
                REQUIRE(path.mapping(1).edit(0).to_length() == 4);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(1).from_length() == 5);
                REQUIRE(path.mapping(1).edit(1).to_length() == 0);
                REQUIRE(path.mapping(1).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(2).from_length() == 2);
                REQUIRE(path.mapping(1).edit(2).to_length() == 2);
                REQUIRE(path.mapping(1).edit(2).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
        }
        
        TEST_CASE( "Banded global aligner produces correct alignments with base quality adjustments",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner runs with base quality adjusted alignments" ) {
                
                VG graph;
                
                QualAdjAligner aligner = QualAdjAligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AGTGCTGAAGT");
                string qual = string("HHHHHHHHHHH");
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner produces correct base quality adjusted alignments with variable base qualities" ) {
                
                VG graph;
                
                QualAdjAligner aligner = QualAdjAligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AGTGCTGAAGT");
                string qual = string("HHDD><<9861");
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 4);
                REQUIRE(path.mapping(0).edit(0).to_length() == 4);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner increases score when adjusting for low quality mismatches" ) {
                
                VG graph;
                
                QualAdjAligner aligner = QualAdjAligner();
                
                Node* n0 = graph.create_node("AGTGCTGAAGT");
                
                string read = string("AGTGCTCAAGT");
                string qual_full = string("HHHHHHHHHHH");
                string qual_reduced = string("HHHHHH$HHHH");
                
                Alignment aln_full;
                aln_full.set_sequence(read);
                aln_full.set_quality(qual_full);
                alignment_quality_char_to_short(aln_full);
                
                Alignment aln_reduced;
                aln_reduced.set_sequence(read);
                aln_reduced.set_quality(qual_reduced);
                alignment_quality_char_to_short(aln_reduced);
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln_full, graph, band_width, permissive_banding);
                aligner.align_global_banded(aln_reduced, graph, band_width, permissive_banding);
                
                const Path& path_full = aln_full.path();
                const Path& path_reduced = aln_reduced.path();
                
                // is a global alignment
                REQUIRE(path_full.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path_full.mapping(path_full.mapping_size() - 1)) == graph.graph.node(path_full.mapping(path_full.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                REQUIRE(path_reduced.mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(path_reduced.mapping(path_reduced.mapping_size() - 1)) == graph.graph.node(path_reduced.mapping(path_reduced.mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(path_full.mapping(0).position().node_id() == 1);
                
                REQUIRE(path_reduced.mapping(0).position().node_id() == 1);
                
                // has corrects edits
                REQUIRE(path_full.mapping(0).edit(0).from_length() == 6);
                REQUIRE(path_full.mapping(0).edit(0).to_length() == 6);
                REQUIRE(path_full.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path_full.mapping(0).edit(1).from_length() == 1);
                REQUIRE(path_full.mapping(0).edit(1).to_length() == 1);
                REQUIRE(path_full.mapping(0).edit(1).sequence() == string("C"));
                
                REQUIRE(path_full.mapping(0).edit(2).from_length() == 4);
                REQUIRE(path_full.mapping(0).edit(2).to_length() == 4);
                REQUIRE(path_full.mapping(0).edit(2).sequence().empty());
                
                REQUIRE(path_reduced.mapping(0).edit(0).from_length() == 6);
                REQUIRE(path_reduced.mapping(0).edit(0).to_length() == 6);
                REQUIRE(path_reduced.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path_reduced.mapping(0).edit(1).from_length() == 1);
                REQUIRE(path_reduced.mapping(0).edit(1).to_length() == 1);
                REQUIRE(path_reduced.mapping(0).edit(1).sequence() == string("C"));
                
                REQUIRE(path_reduced.mapping(0).edit(2).from_length() == 4);
                REQUIRE(path_reduced.mapping(0).edit(2).to_length() == 4);
                REQUIRE(path_reduced.mapping(0).edit(2).sequence().empty());
                
                // scores behave as expected
                REQUIRE(aln_reduced.score() > aln_full.score());
            }
        }
        
        TEST_CASE( "Banded global aligner produces correct alternate alignments in different traceback scenarios",
                  "[alignment][multialignment][banded][mapping]" ) {
            
            
            SECTION( "Banded global aligner returns alternate alignments in descending score order" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("ACGTAGTCTGAA");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("TT");
                Node* n3 = graph.create_node("TGACGTACGTTA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("ACGTAGTCTGACATGACGTACGTTA");
                Alignment aln;
                aln.set_sequence(read);
                
                
                int max_multi_alns = 20;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                const Path& path = aln.path();
                
                
                int64_t prev_score = numeric_limits<int64_t>::max();
                for (Alignment& alt_aln : multi_alns) {
                    const Path& path = alt_aln.path();
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // score is in descending order
                    REQUIRE(aln.score() <= prev_score);
                    prev_score = aln.score();
                }
            }
            
            SECTION( "Banded global aligner stores the optimal alignment in both the main Alignment object and the first position in the return vector in multiple alignment" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("ACGTAGTCTGAA");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("TT");
                Node* n3 = graph.create_node("TGACGTACGTTA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("ACGTAGTCTGACATGACGTACGTTA");
                Alignment aln;
                aln.set_sequence(read);
                
                int max_multi_alns = 20;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                REQUIRE(aln.sequence() == multi_alns[0].sequence());
                REQUIRE(aln.score() == multi_alns[0].score());
                REQUIRE(aln.identity() == multi_alns[0].identity());
                REQUIRE(aln.path().mapping_size() == multi_alns[0].path().mapping_size());
                for (int i = 0; i < aln.path().mapping_size(); i++) {
                    REQUIRE(aln.path().mapping(i).position().node_id() == multi_alns[0].path().mapping(i).position().node_id());
                    REQUIRE(aln.path().mapping(i).position().offset() == multi_alns[0].path().mapping(i).position().offset());
                    REQUIRE(aln.path().mapping(i).edit_size() == multi_alns[0].path().mapping(i).edit_size());
                    for (int j = 0; j < aln.path().mapping(i).edit_size(); j++) {
                        REQUIRE(aln.path().mapping(i).edit(j).from_length() == multi_alns[0].path().mapping(i).edit(j).from_length());
                        REQUIRE(aln.path().mapping(i).edit(j).to_length() == multi_alns[0].path().mapping(i).edit(j).to_length());
                        REQUIRE(aln.path().mapping(i).edit(j).sequence() == multi_alns[0].path().mapping(i).edit(j).sequence());
                    }
                }
            }
            
            SECTION( "Banded global aligner can identify both optimal alignments when there is a deletion in a homodimer" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("ACGTAGTCTGAA");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("TT");
                Node* n3 = graph.create_node("TGACGTACGTTA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("ACGTAGTCTGACATGACGTACGTTA");
                Alignment aln;
                aln.set_sequence(read);
                
                int max_multi_alns = 2;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    const Path& path = alt_aln.path();
                    
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // follows correct path
                    is_first_opt = is_first_opt && (path.mapping(0).position().node_id() == 1);
                    is_first_opt = is_first_opt && (path.mapping(1).position().node_id() == 2);
                    is_first_opt = is_first_opt && (path.mapping(2).position().node_id() == 4);
                    
                    // has corrects edits
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).from_length() == 11);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).to_length() == 11);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(0).edit(1).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(1).to_length() == 0);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(1).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).from_length() == 2);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).to_length() == 2);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).from_length() == 12);
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).to_length() == 12);
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).sequence().empty());
                    
                    // follows correct path
                    is_second_opt = is_second_opt && (path.mapping(0).position().node_id() == 1);
                    is_second_opt = is_second_opt && (path.mapping(1).position().node_id() == 2);
                    is_second_opt = is_second_opt && (path.mapping(2).position().node_id() == 4);
                    
                    // has corrects edits
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).from_length() == 10);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).to_length() == 10);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(0).edit(1).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(1).to_length() == 0);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(1).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(0).edit(2).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(2).to_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(2).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).from_length() == 2);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).to_length() == 2);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).from_length() == 12);
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).to_length() == 12);
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).sequence().empty());
                    
                    found_first_opt = found_first_opt || is_first_opt;
                    found_second_opt = found_second_opt || is_second_opt;
                }
                
                REQUIRE(found_first_opt);
                REQUIRE(found_second_opt);
            }
            
            SECTION( "Banded global aligner can identify alternate alignments that follow a different node sequence than the optimal alignment" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("ACGTAGTCTGAA");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("TGACGTACGTTA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("ACGTAGTCTGAACTGACGTACGTTA");
                Alignment aln;
                aln.set_sequence(read);
                
                int max_multi_alns = 20;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                bool took_alternate_path = false;
                for (Alignment& alt_aln : multi_alns) {
                    const Path& path = alt_aln.path();
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    if (path.mapping(1).position().node_id() != aln.path().mapping(1).position().node_id()) {
                        took_alternate_path = true;
                        break;
                    }
                }
                
                REQUIRE(took_alternate_path);
            }
            
            SECTION( "Banded global aligner can identify an alternate alignment that branches from another alternate alignment at a node boundary" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AAAAAAAA");
                Node* n1 = graph.create_node("GGG");
                Node* n2 = graph.create_node("G");
                Node* n3 = graph.create_node("C");
                Node* n4 = graph.create_node("T");
                Node* n5 = graph.create_node("G");
                Node* n6 = graph.create_node("AAAAAA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n5);
                graph.create_edge(n4, n5);
                graph.create_edge(n1, n6);
                graph.create_edge(n5, n6);
                
                string read = string("AAAAAAAAGGGAAAAAA");
                Alignment aln;
                aln.set_sequence(read);
                
                int max_multi_alns = 3;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    
                    const Path& path = alt_aln.path();
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // follows correct path
                    is_first_opt = is_first_opt && (path.mapping(0).position().node_id() == 1);
                    is_first_opt = is_first_opt && (path.mapping(1).position().node_id() == 3);
                    is_first_opt = is_first_opt && (path.mapping(2).position().node_id() == 4);
                    is_first_opt = is_first_opt && (path.mapping(3).position().node_id() == 6);
                    is_first_opt = is_first_opt && (path.mapping(4).position().node_id() == 7);
                    
                    // has corrects edits
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).from_length() == 8);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).to_length() == 8);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).to_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).to_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).sequence() == "G");
                    
                    is_first_opt = is_first_opt && (path.mapping(3).edit(0).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(3).edit(0).to_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(3).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(4).edit(0).from_length() == 6);
                    is_first_opt = is_first_opt && (path.mapping(4).edit(0).to_length() == 6);
                    is_first_opt = is_first_opt && (path.mapping(4).edit(0).sequence().empty());
                    
                    // follows correct path
                    is_second_opt = is_second_opt && (path.mapping(0).position().node_id() == 1);
                    is_second_opt = is_second_opt && (path.mapping(1).position().node_id() == 3);
                    is_second_opt = is_second_opt && (path.mapping(2).position().node_id() == 5);
                    is_second_opt = is_second_opt && (path.mapping(3).position().node_id() == 6);
                    is_second_opt = is_second_opt && (path.mapping(4).position().node_id() == 7);
                    
                    // has corrects edits
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).from_length() == 8);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).to_length() == 8);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).to_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).to_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).sequence() == "G");
                    
                    is_second_opt = is_second_opt && (path.mapping(3).edit(0).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(3).edit(0).to_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(3).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(4).edit(0).from_length() == 6);
                    is_second_opt = is_second_opt && (path.mapping(4).edit(0).to_length() == 6);
                    is_second_opt = is_second_opt && (path.mapping(4).edit(0).sequence().empty());
                    
                    found_first_opt = found_first_opt || is_first_opt;
                    found_second_opt = found_second_opt || is_second_opt;
                    
                }
                
                REQUIRE(found_first_opt);
                REQUIRE(found_second_opt);
            }
            
            
            
            SECTION( "Banded global aligner can identify four alignments with equal scores" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("A");
                Node* n1 = graph.create_node("T");
                Node* n2 = graph.create_node("G");
                Node* n3 = graph.create_node("A");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("C");
                Alignment aln;
                aln.set_sequence(read);
                
                int max_multi_alns = 4;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                bool found_third_opt = false;
                bool found_fourth_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    bool is_third_opt = true;
                    bool is_fourth_opt = true;
                    
                    const Path& path = alt_aln.path();
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // follows correct path
                    REQUIRE(path.mapping(0).position().node_id() == n0->id());
                    is_first_opt = is_first_opt && (path.mapping(1).position().node_id() == n1->id());
                    is_second_opt = is_second_opt && (path.mapping(1).position().node_id() == n1->id());
                    is_third_opt = is_third_opt && (path.mapping(1).position().node_id() == n2->id());
                    is_fourth_opt = is_fourth_opt && (path.mapping(1).position().node_id() == n2->id());
                    REQUIRE(path.mapping(2).position().node_id() == n3->id());
                    
                    // has corrects edit
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).to_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).sequence() == "C");
                    
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).to_length() == 0);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_third_opt = is_third_opt && (path.mapping(0).edit(0).from_length() == 1);
                    is_third_opt = is_third_opt && (path.mapping(0).edit(0).to_length() == 1);
                    is_third_opt = is_third_opt && (path.mapping(0).edit(0).sequence() == "C");
                    
                    is_fourth_opt = is_fourth_opt && (path.mapping(0).edit(0).from_length() == 1);
                    is_fourth_opt = is_fourth_opt && (path.mapping(0).edit(0).to_length() == 0);
                    is_fourth_opt = is_fourth_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                    REQUIRE(path.mapping(1).edit(0).to_length() == 0);
                    REQUIRE(path.mapping(1).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).to_length() == 0);
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).to_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).sequence() == "C");
                    
                    is_third_opt = is_third_opt && (path.mapping(2).edit(0).from_length() == 1);
                    is_third_opt = is_third_opt && (path.mapping(2).edit(0).to_length() == 0);
                    is_third_opt = is_third_opt && (path.mapping(2).edit(0).sequence().empty());
                    
                    is_fourth_opt = is_fourth_opt && (path.mapping(2).edit(0).from_length() == 1);
                    is_fourth_opt = is_fourth_opt && (path.mapping(2).edit(0).to_length() == 1);
                    is_fourth_opt = is_fourth_opt && (path.mapping(2).edit(0).sequence() == "C");
                    
                    REQUIRE(alt_aln.score() == -11);
                    
                    found_first_opt = found_first_opt || is_first_opt;
                    found_second_opt = found_second_opt || is_second_opt;
                    found_third_opt = found_third_opt || is_third_opt;
                    found_fourth_opt = found_fourth_opt || is_fourth_opt;
                }
                
                REQUIRE(found_first_opt);
                REQUIRE(found_second_opt);
                REQUIRE(found_third_opt);
                REQUIRE(found_fourth_opt);
            }
            
            SECTION( "Banded global aligner can identify an alternate alignment that branches from another alternate alignment inside a node sequence" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AAAAAAAAAA");
                Node* n1 = graph.create_node("CGGC");
                Node* n2 = graph.create_node("CGGT");
                Node* n3 = graph.create_node("AAAAAAAAAA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAAAAAAAAACGGGCAAAAAAAAAA");
                Alignment aln;
                aln.set_sequence(read);
                
                int max_multi_alns = 10;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    
                    const Path& path = alt_aln.path();
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // follows correct path
                    is_first_opt = is_first_opt && (path.mapping(0).position().node_id() == 1);
                    is_first_opt = is_first_opt && (path.mapping(1).position().node_id() == 3);
                    is_first_opt = is_first_opt && (path.mapping(2).position().node_id() == 4);
                    
                    // has corrects edits
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).from_length() == 10);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).to_length() == 10);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    if (path.mapping(1).edit_size() >= 4) {
                        is_first_opt = is_first_opt && (path.mapping(1).edit(0).from_length() == 1);
                        is_first_opt = is_first_opt && (path.mapping(1).edit(0).to_length() == 1);
                        is_first_opt = is_first_opt && (path.mapping(1).edit(0).sequence().empty());
                        
                        is_first_opt = is_first_opt && (path.mapping(1).edit(1).from_length() == 0);
                        is_first_opt = is_first_opt && (path.mapping(1).edit(1).to_length() == 1);
                        is_first_opt = is_first_opt && (path.mapping(1).edit(1).sequence() == "G");
                        
                        is_first_opt = is_first_opt && (path.mapping(1).edit(2).from_length() == 2);
                        is_first_opt = is_first_opt && (path.mapping(1).edit(2).to_length() == 2);
                        is_first_opt = is_first_opt && (path.mapping(1).edit(2).sequence().empty());
                        
                        is_first_opt = is_first_opt && (path.mapping(1).edit(3).from_length() == 1);
                        is_first_opt = is_first_opt && (path.mapping(1).edit(3).to_length() == 1);
                        is_first_opt = is_first_opt && (path.mapping(1).edit(3).sequence() == "C");
                    }
                    else {
                        is_first_opt = false;
                    }
                    
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).from_length() == 10);
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).to_length() == 10);
                    is_first_opt = is_first_opt && (path.mapping(2).edit(0).sequence().empty());
                    
                    // follows correct path
                    is_second_opt = is_second_opt && (path.mapping(0).position().node_id() == 1);
                    is_second_opt = is_second_opt && (path.mapping(1).position().node_id() == 3);
                    is_second_opt = is_second_opt && (path.mapping(2).position().node_id() == 4);
                    
                    // has corrects edits
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).from_length() == 10);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).to_length() == 10);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    if (path.mapping(1).edit_size() >= 4) {
                        is_second_opt = is_second_opt && (path.mapping(1).edit(0).from_length() == 2);
                        is_second_opt = is_second_opt && (path.mapping(1).edit(0).to_length() == 2);
                        is_second_opt = is_second_opt && (path.mapping(1).edit(0).sequence().empty());
                        
                        is_second_opt = is_second_opt && (path.mapping(1).edit(1).from_length() == 0);
                        is_second_opt = is_second_opt && (path.mapping(1).edit(1).to_length() == 1);
                        is_second_opt = is_second_opt && (path.mapping(1).edit(1).sequence() == "G");
                        
                        is_second_opt = is_second_opt && (path.mapping(1).edit(2).from_length() == 1);
                        is_second_opt = is_second_opt && (path.mapping(1).edit(2).to_length() == 1);
                        is_second_opt = is_second_opt && (path.mapping(1).edit(2).sequence().empty());
                        
                        is_second_opt = is_second_opt && (path.mapping(1).edit(3).from_length() == 1);
                        is_second_opt = is_second_opt && (path.mapping(1).edit(3).to_length() == 1);
                        is_second_opt = is_second_opt && (path.mapping(1).edit(3).sequence() == "C");
                    }
                    else {
                        is_second_opt = false;
                    }
                    
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).from_length() == 10);
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).to_length() == 10);
                    is_second_opt = is_second_opt && (path.mapping(2).edit(0).sequence().empty());
                    
                    found_first_opt = found_first_opt || is_first_opt;
                    found_second_opt = found_second_opt || is_second_opt;
                    
                }
                
                REQUIRE(found_first_opt);
                REQUIRE(found_second_opt);
            }
            
            
            
            SECTION( "Banded global aligner can identify both of two lead deletions of equal length" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("C");
                Node* n1 = graph.create_node("T");
                Node* n2 = graph.create_node("GA");
                
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n2);
                
                string read = string("A");
                Alignment aln;
                aln.set_sequence(read);
                
                int max_multi_alns = 2;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    const Path& path = alt_aln.path();
                    
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // follows correct path
                    is_first_opt = is_first_opt && (path.mapping(0).position().node_id() == n0->id());
                    is_first_opt = is_first_opt && (path.mapping(1).position().node_id() == n2->id());
                    
                    // has corrects edits
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).to_length() == 0);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).to_length() == 0);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(1).edit(1).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(1).to_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(1).sequence().empty());
                    
                    // follows correct path
                    is_second_opt = is_second_opt && (path.mapping(0).position().node_id() == n1->id());
                    is_second_opt = is_second_opt && (path.mapping(1).position().node_id() == n2->id());
                    
                    // has corrects edits
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).to_length() == 0);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).to_length() == 0);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(1).edit(1).from_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(1).to_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(1).sequence().empty());
                    
                    found_first_opt = found_first_opt || is_first_opt;
                    found_second_opt = found_second_opt || is_second_opt;
                }
                
                REQUIRE(found_first_opt);
                REQUIRE(found_second_opt);
            }
            
            SECTION( "Global banded aligner does not produce duplicate alternate alignments" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                // low complexity sequences to ensure many alternate alignments
                Node* n0 = graph.create_node("CCCCCCCCCTCCCCCCCCCCTCCCCCCCCCCGACCCCCCCCCCC");
                Node* n1 = graph.create_node("CCCCCCCCCCACCCCCCCCCCACCCCCCCCCCTCCCA");
                Node* n2 = graph.create_node("CCCCCACCCCCCCCGTCCCCCCCCCCCA");
                Node* n3 = graph.create_node("CCCCCCCCCCCCGCCCCCCCCCCGCCCCCCCCC");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = "CCCCCCCCCTCCCCCCCCCCTCCCCCCCCCCGACCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCACCCCCCCCCCTCCCACCCCCCCCCCCCGCCCCCCCCCCGCCCCCCCCC";
                Alignment aln;
                aln.set_sequence(read);
                
                int max_multi_alns = 2000;
                int band_padding = 1;
                bool permissive_banding = true;
                vector<Alignment> multi_alns;
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                unordered_set<string> alns_seen;
                for (Alignment& alt_aln : multi_alns) {
                    string aln_string = hash_alignment(alt_aln);
                    
                    REQUIRE(alns_seen.count(aln_string) == 0);
                    alns_seen.insert(aln_string);
                }
            }
        }
        
        TEST_CASE( "Banded global aligner works with Ns",
                  "[alignment][banded][mapping]"  ) {
            
            SECTION( "Banded global aligner can align Ns to letters" ) {
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("NNNGCTGANNN");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                SECTION("alignment ends in full-length matches/mismatches") {
                    REQUIRE(aln.path().mapping_size() == 3);
                    REQUIRE(mapping_from_length(aln.path().mapping(0)) == 4);
                    REQUIRE(mapping_to_length(aln.path().mapping(0)) == 4);
                    REQUIRE(mapping_from_length(aln.path().mapping(2)) == 6);
                    REQUIRE(mapping_to_length(aln.path().mapping(2)) == 6);
                }
                
                
            }
            
            SECTION( "Banded global aligner can align letters to Ns" ) {
                VG graph;
                
                TestAligner aligner_source;
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
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                SECTION("alignment ends in full-length matches/mismatches") {
                    REQUIRE(aln.path().mapping_size() == 3);
                    REQUIRE(mapping_from_length(aln.path().mapping(0)) == 4);
                    REQUIRE(mapping_to_length(aln.path().mapping(0)) == 4);
                    REQUIRE(mapping_from_length(aln.path().mapping(2)) == 6);
                    REQUIRE(mapping_to_length(aln.path().mapping(2)) == 6);
                }
                
                
            }
            
            SECTION( "Banded global aligner can align Ns to Ns" ) {
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("NNNG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGANNN");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("NNNGCTGANNN");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, graph, band_width, permissive_banding);
                
                SECTION("alignment ends in full-length matches/mismatches") {
                    REQUIRE(aln.path().mapping_size() == 3);
                    REQUIRE(mapping_from_length(aln.path().mapping(0)) == 4);
                    REQUIRE(mapping_to_length(aln.path().mapping(0)) == 4);
                    REQUIRE(mapping_from_length(aln.path().mapping(2)) == 6);
                    REQUIRE(mapping_to_length(aln.path().mapping(2)) == 6);
                }
            }
        }

        
        TEST_CASE( "Banded global aligner can align to graphs with empty nodes",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner can align to a graph with an empty source node") {
                
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("");
                Node* n1 = graph.create_node("CT");
                Node* n2 = graph.create_node("TA");
                Node* n3 = graph.create_node("GA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "CTGA";
                string qual = "HHHH";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
                
                aln.Clear();
                read = "AGA";
                qual = "HHH";
                
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n2->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(1).from_length() == 1);
                REQUIRE(aln.path().mapping(1).edit(1).to_length() == 1);
                REQUIRE(aln.path().mapping(1).edit(1).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
                
            }
            
            SECTION( "Banded global aligner can align to a graph with an empty sink node") {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("GA");
                Node* n1 = graph.create_node("CT");
                Node* n2 = graph.create_node("TA");
                Node* n3 = graph.create_node("");

                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "GACT";
                string qual = "HHHH";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner can align to a graph with an empty source and sink node") {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("");
                Node* n1 = graph.create_node("CT");
                Node* n2 = graph.create_node("TA");
                Node* n3 = graph.create_node("");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "CT";
                string qual = "HH";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner can align to a graph with a chained empty source and sink nodes") {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("");
                Node* n1 = graph.create_node("");
                Node* n2 = graph.create_node("CT");
                Node* n3 = graph.create_node("TA");
                Node* n4 = graph.create_node("");
                Node* n5 = graph.create_node("");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "CT";
                string qual = "HH";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n2->id());
                REQUIRE(aln.path().mapping(3).position().node_id() == n4->id());
                REQUIRE(aln.path().mapping(4).position().node_id() == n5->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(3).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(3).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(3).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(4).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(4).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(4).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner can align to a graph with an empty nodes that is both a source and sink") {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("");
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "CT";
                string qual = "HH";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(0).edit(0).sequence() == read);
            }
            
            SECTION( "Banded global aligner can align to a graph with both empty and non-empty sources and sinks") {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("");
                Node* n1 = graph.create_node("GA");
                Node* n2 = graph.create_node("CT");
                Node* n3 = graph.create_node("GA");
                Node* n4 = graph.create_node("");
                
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n2);
                graph.create_edge(n2, n3);
                graph.create_edge(n2, n4);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "CTGA";
                string qual = "HHHH";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n2->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
                
                
                aln.Clear();
                read = "GACT";
                qual = "HHHH";
                
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n2->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n4->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner can align to a graph with empty interior nodes") {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("GA");
                Node* n1 = graph.create_node("");
                Node* n2 = graph.create_node("");
                Node* n3 = graph.create_node("CT");
                Node* n4 = graph.create_node("TA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n3);
                graph.create_edge(n1, n2);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "GATA";
                string qual = "HHHH";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n2->id());
                REQUIRE(aln.path().mapping(3).position().node_id() == n4->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(3).edit(0).from_length() == 2);
                REQUIRE(aln.path().mapping(3).edit(0).to_length() == 2);
                REQUIRE(aln.path().mapping(3).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner can align to an empty graph" ) {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("");
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "G";
                string qual = "H";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 1);
                REQUIRE(aln.path().mapping(0).edit(0).sequence() == read);
            }
            
            
            SECTION( "Banded global aligner can align to an empty graph of more than one node" ) {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("");
                Node* n1 = graph.create_node("");
                
                graph.create_edge(n0, n1);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "G";
                string qual = "H";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 0);
                
                REQUIRE((aln.path().mapping(0).edit(0).to_length() == 1) != (aln.path().mapping(1).edit(0).to_length() == 1));
                REQUIRE((aln.path().mapping(0).edit(0).to_length() == 0) != (aln.path().mapping(1).edit(0).to_length() == 0));
                
                REQUIRE((aln.path().mapping(0).edit(0).sequence() == read) != (aln.path().mapping(1).edit(0).sequence() == read));
                REQUIRE((aln.path().mapping(0).edit(0).sequence().empty()) != (aln.path().mapping(1).edit(0).sequence().empty()));
            }
            
            SECTION( "Banded global aligner can align to a graph with empty and non-empty paths" ) {
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("");
                Node* n1 = graph.create_node("");
                Node* n2 = graph.create_node("TCA");
                
                graph.create_edge(n0, n1);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "G";
                string qual = "H";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 0);
                
                REQUIRE((aln.path().mapping(0).edit(0).to_length() == 1) != (aln.path().mapping(1).edit(0).to_length() == 1));
                REQUIRE((aln.path().mapping(0).edit(0).to_length() == 0) != (aln.path().mapping(1).edit(0).to_length() == 0));
                
                REQUIRE((aln.path().mapping(0).edit(0).sequence() == read) != (aln.path().mapping(1).edit(0).sequence() == read));
                REQUIRE((aln.path().mapping(0).edit(0).sequence().empty()) != (aln.path().mapping(1).edit(0).sequence().empty()));
                
                
                aln.Clear();
                read = "TCA";
                qual = "HHH";
                
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n2->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 3);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 3);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
            }
            
            SECTION( "Banded global aligner can find multi-alignments over empty sink nodes" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("GA");
                Node* n1 = graph.create_node("");
                Node* n2 = graph.create_node("");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                
                bool permissive_banding = true;
                int band_padding = 1;
                int max_multi_alns = 2;
                vector<Alignment> multi_alns;
                
                string read = "GA";
                string qual = "HH";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    const Path& path = alt_aln.path();
                    
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // follows correct path
                    is_first_opt = is_first_opt && (path.mapping(0).position().node_id() == n0->id());
                    is_first_opt = is_first_opt && (path.mapping(1).position().node_id() == n1->id());
                    
                    // has corrects edits
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).from_length() == 2);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).to_length() == 2);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).from_length() == 0);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).to_length() == 0);
                    is_first_opt = is_first_opt && (path.mapping(1).edit(0).sequence().empty());
                    
                    // follows correct path
                    is_second_opt = is_second_opt && (path.mapping(0).position().node_id() == n0->id());
                    is_second_opt = is_second_opt && (path.mapping(1).position().node_id() == n2->id());
                    
                    // has corrects edits
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).from_length() == 2);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).to_length() == 2);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).from_length() == 0);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).to_length() == 0);
                    is_second_opt = is_second_opt && (path.mapping(1).edit(0).sequence().empty());
                    
                    found_first_opt = found_first_opt || is_first_opt;
                    found_second_opt = found_second_opt || is_second_opt;
                }
                
                REQUIRE(found_first_opt);
                REQUIRE(found_second_opt);
                
            }
            
            // note: alternate alignments over different paths of empty nodes from the same non-empty node
            // to either the same non-empty node or a source are not supported
            
            SECTION( "Banded global aligner can find multi-alignments over empty and non-empty separated components" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("C");
                Node* n1 = graph.create_node("TT");
                Node* n2 = graph.create_node("");
                
                bool permissive_banding = true;
                int band_padding = 1;
                int max_multi_alns = 3;
                vector<Alignment> multi_alns;
                
                string read = "C";
                string qual = "H";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                REQUIRE(multi_alns.size() <= 3);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                bool found_third_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    
                    const Path& path = alt_aln.path();
                    
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    bool is_third_opt = true;
                    
                    // is a global alignment
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == graph.graph.node(path.mapping(path.mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // follows correct path
                    is_first_opt = is_first_opt && (path.mapping(0).position().node_id() == n0->id());
                    
                    // has corrects edits
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).from_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).to_length() == 1);
                    is_first_opt = is_first_opt && (path.mapping(0).edit(0).sequence().empty());
                    
                    // follows correct path
                    is_second_opt = is_second_opt && (path.mapping(0).position().node_id() == n2->id());
                    
                    // has corrects edits
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).from_length() == 0);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).to_length() == 1);
                    is_second_opt = is_second_opt && (path.mapping(0).edit(0).sequence() == read);
                    
                    // follows correct path
                    is_third_opt = is_third_opt && (path.mapping(0).position().node_id() == n1->id());
                    
                    // has corrects edits
                    is_third_opt = is_third_opt && (path.mapping(0).edit(0).from_length() == 1);
                    is_third_opt = is_third_opt && (path.mapping(0).edit(1).from_length() == 1);
                    is_third_opt = is_third_opt && (path.mapping(0).edit(0).to_length() == 0 || path.mapping(0).edit(0).to_length() == 1);
                    is_third_opt = is_third_opt && (path.mapping(0).edit(1).to_length() == 0 || path.mapping(0).edit(1).to_length() == 1);
                    is_third_opt = is_third_opt && (path.mapping(0).edit(0).to_length() == 0 != path.mapping(0).edit(1).to_length() == 0);
                    
                    found_first_opt = found_first_opt || is_first_opt;
                    found_second_opt = found_second_opt || is_second_opt;
                    found_third_opt = found_third_opt || is_third_opt;
                }
                
                REQUIRE(found_first_opt);
                REQUIRE(found_second_opt);
                REQUIRE(found_third_opt);
            }
        }
        
        TEST_CASE( "Banded global aligner can align empty reads",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner can align an empty read to a graph with only one path") {
                
                VG graph;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                Node* n0 = graph.create_node("A");
                Node* n1 = graph.create_node("G");
                Node* n2 = graph.create_node("T");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n1, n2);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "";
                string qual = "";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n2->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
                
            }
            
            SECTION( "Banded global aligner can align an empty read to a graph with multiple paths") {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("A");
                Node* n1 = graph.create_node("G");
                Node* n2 = graph.create_node("TC");
                Node* n3 = graph.create_node("C");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "";
                string qual = "";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
                
                REQUIRE(aln.score() == -8);
            }
            
            SECTION( "Banded global aligner can align an empty read to a graph with empty nodes" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("");
                Node* n1 = graph.create_node("G");
                Node* n2 = graph.create_node("TC");
                Node* n3 = graph.create_node("");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "";
                string qual = "";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                REQUIRE(aln.path().mapping(1).position().node_id() == n1->id());
                REQUIRE(aln.path().mapping(2).position().node_id() == n3->id());
                REQUIRE(aln.path().mapping(3).position().node_id() == n4->id());
                REQUIRE(aln.path().mapping(4).position().node_id() == n5->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(2).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(2).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(3).edit(0).from_length() == 1);
                REQUIRE(aln.path().mapping(3).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(3).edit(0).sequence().empty());
                
                REQUIRE(aln.path().mapping(4).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(4).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(4).edit(0).sequence().empty());
                
                REQUIRE(aln.score() == -7);
            }
            
            SECTION( "Banded global aligner can align an empty read to an empty graph" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("");
                
                bool permissive_banding = true;
                int band_padding = 1;
                
                string read = "";
                string qual = "";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                aligner.align_global_banded(aln, graph, band_padding, permissive_banding);
                
                // is a global alignment
                REQUIRE(aln.path().mapping(0).position().offset() == 0);
                REQUIRE(mapping_from_length(aln.path().mapping(aln.path().mapping_size() - 1)) == graph.graph.node(aln.path().mapping(aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                
                // follows correct path
                REQUIRE(aln.path().mapping(0).position().node_id() == n0->id());
                
                // has corrects edits
                REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
                REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
                
                REQUIRE(aln.score() == 0);
            }
            
            SECTION( "Banded global aligner obtain multi-alignments from an empty read" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                Node* n0 = graph.create_node("A");
                Node* n1 = graph.create_node("GG");
                Node* n2 = graph.create_node("TC");
                Node* n3 = graph.create_node("CA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                bool permissive_banding = true;
                int band_padding = 1;
                int max_multi_alns = 2;
                
                string read = "";
                string qual = "";
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                vector<Alignment> multi_alns;
                aligner.align_global_banded_multi(aln, multi_alns, graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                bool found_first_opt = false, found_second_opt = false;
                
                for (auto& alt_aln : multi_alns) {
                    
                    // is a global alignment
                    REQUIRE(alt_aln.path().mapping(0).position().offset() == 0);
                    REQUIRE(mapping_from_length(alt_aln.path().mapping(alt_aln.path().mapping_size() - 1)) == graph.graph.node(alt_aln.path().mapping(alt_aln.path().mapping_size() - 1).position().node_id() - 1).sequence().length());
                    
                    // follows correct path
                    REQUIRE(alt_aln.path().mapping(0).position().node_id() == n0->id());
                    bool is_first_opt = alt_aln.path().mapping(1).position().node_id() == n1->id();
                    bool is_second_opt = alt_aln.path().mapping(1).position().node_id() == n2->id();
                    REQUIRE(alt_aln.path().mapping(2).position().node_id() == n3->id());
                    
                    // has corrects edits
                    REQUIRE(alt_aln.path().mapping(0).edit(0).from_length() == 1);
                    REQUIRE(alt_aln.path().mapping(0).edit(0).to_length() == 0);
                    REQUIRE(alt_aln.path().mapping(0).edit(0).sequence().empty());
                    
                    REQUIRE(alt_aln.path().mapping(1).edit(0).from_length() == 2);
                    REQUIRE(alt_aln.path().mapping(1).edit(0).to_length() == 0);
                    REQUIRE(alt_aln.path().mapping(1).edit(0).sequence().empty());
                    
                    REQUIRE(alt_aln.path().mapping(2).edit(0).from_length() == 2);
                    REQUIRE(alt_aln.path().mapping(2).edit(0).to_length() == 0);
                    REQUIRE(alt_aln.path().mapping(2).edit(0).sequence().empty());
                    
                    REQUIRE(alt_aln.score() == -10);
                    
                    found_first_opt = found_first_opt || is_first_opt;
                    found_second_opt = found_second_opt || is_second_opt;
                }
                
                REQUIRE(found_first_opt);
                REQUIRE(found_second_opt);
            }
            
            
            
            SECTION( "Banded global aligner does not produce empty edits when there is an insertion an empty node") {
                string graph_json = R"({"edge": [{"to_end": true, "from_start": true, "to": 22, "from": 20}, {"to": 26, "from": 20}, {"to": 24, "from": 20}, {"to_end": true, "from_start": true, "to": 26, "from": 4}, {"to_end": true, "from_start": true, "to": 24, "from": 4}], "node": [{"sequence": "C", "id": 24}, {"sequence": "GAGA", "id": 20}, {"sequence": "T", "id": 26}, {"sequence": "GGAGTCT", "id": 4}, {"id": 22}]})";
                
                Graph graph;
                json2pb(graph, graph_json.c_str(), graph_json.size());
                VG vg_graph(graph);
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                string read = "TTTTGATGGAGGCC";
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                bool permissive_banding = true;
                aligner.align_global_banded(aln, vg_graph, band_width, permissive_banding);
                
                const Path& p = aln.path();
                for (int i = 0; i < p.mapping_size();i++) {
                    const auto& m = p.mapping(i);
                    for (int j = 0; j < m.edit_size(); j++) {
                        const auto& e = m.edit(j);
                        bool empty = e.to_length() == 0 && e.from_length() == 0;
                        REQUIRE(!empty);
                    }
                }
            }
        }
    
        TEST_CASE( "Banded global aligner doesn't crash when the worst possible score is just on the edge of needing a larger int size",
                  "[alignment][banded][mapping]" ) {
            
            bdsg::HashGraph graph;
            
            handle_t h0 = graph.create_handle("", 68181350);
            handle_t h1 = graph.create_handle("G", 68181343);
            handle_t h2 = graph.create_handle("TGAGTGG", 68181344);
            handle_t h3 = graph.create_handle("CTTTGGTTCCCGGCTGAGGTGGAGTGGGCTGA", 68181345);
            handle_t h4 = graph.create_handle("GGACTAGACTGAGCCCTCGGACATGGAGGTGG", 68181346);
            handle_t h5 = graph.create_handle("GGATGGGGCAGACTCATCCCATTCTTGACCAA", 68181347);
            handle_t h6 = graph.create_handle("GCCCTTGTTCTGCTCCCTTCCCAG", 68181348);
            handle_t h7 = graph.create_handle("", 68181349);
            
            graph.create_edge(h0, h1);
            graph.create_edge(h0, h7);
            graph.create_edge(h1, h2);
            graph.create_edge(h2, h3);
            graph.create_edge(h3, h4);
            graph.create_edge(h4, h5);
            graph.create_edge(h5, h6);
            graph.create_edge(h6, h7);
            
            string sequence = "AA";
            Alignment aln;
            aln.set_sequence(sequence);
            
            TestAligner aligner_source;
            aligner_source.set_alignment_scores(1, 1, 1, 1, 0);
            const Aligner& aligner = *aligner_source.get_regular_aligner();
            
            aligner.align_global_banded(aln, graph, 2, true);
            
            REQUIRE(aln.path().mapping_size() == 2);
            REQUIRE(aln.path().mapping(0).position().node_id() == graph.get_id(h0));
            REQUIRE(aln.path().mapping(1).position().node_id() == graph.get_id(h7));
        }
    
        TEST_CASE("Try to recreate a memory access bug", "[alignment][banded][mapping][memory]") {
            
            // note: this never crashed, but the bug shows up on valgrind
            
            bdsg::HashGraph graph;
            
            handle_t h0 = graph.create_handle("T");
            
            string sequence = "CTCATTCCCGGAACCTTGAAATGGAGCT";
            string qual = "DCDD=2DECBEC=F@E?BEFEEFECED<";
            Alignment aln;
            aln.set_sequence(sequence);
            aln.set_quality(string_quality_short_to_char(qual));
            
            TestAligner aligner_source;
            const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
            
            aligner.align_global_banded(aln, graph, 1, true);
        }
    }
}






