//
//  banded_global_aligner.cpp
//  
//  Unit tests for banded global aligner module.
//

#include <stdio.h>
#include "catch.hpp"
#include "gssw_aligner.hpp"
#include "vg.hpp"
#include "path.hpp"
#include "banded_global_aligner.hpp"
#include "json2pb.h"

using namespace google::protobuf;

namespace vg {
    namespace unittest {
        
        
        TEST_CASE( "Banded global aligner produces correct alignments with all types of edits",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner produces correct alignment when read matches exactly") {
                
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
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                Aligner aligner;
                
                Node* n0 = graph.create_node("A");
                Node* n1 = graph.create_node("C");
                
                graph.create_edge(n0, n1);
                
                string read = string("C");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_padding = 1;
                
                aligner.align_global_banded(aln, graph.graph, band_padding);
                
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

            SECTION( "Global banded aligner doesn't explode" ) {

                string json_graph = "{\"node\": [{\"sequence\": \"GT\", \"id\": 86051}, {\"sequence\": \"G\", \"id\": 86052}, {\"sequence\": \"A\", \"id\": 86053}, {\"sequence\": \"TGCCCAGTTCTAGACAATT\", \"id\": 86054}, {\"sequence\": \"T\", \"id\": 86055}, {\"sequence\": \"A\", \"id\": 86056}, {\"sequence\": \"CAA\", \"id\": 86057}, {\"sequence\": \"C\", \"id\": 86058}, {\"sequence\": \"A\", \"id\": 86059}, {\"sequence\": \"C\", \"id\": 86060}, {\"sequence\": \"T\", \"id\": 86061}, {\"sequence\": \"C\", \"id\": 86062}, {\"sequence\": \"CAACAAGCTCAGTGCATTTGCCACTT\", \"id\": 86063}, {\"sequence\": \"GATAAGTGCTGAGCTGAATGGTTTGC\", \"id\": 86064}, {\"sequence\": \"AGGTGTTACCTCAAATTCATGCCCAT\", \"id\": 86065}, {\"sequence\": \"CATGCTTAATAAGGTCAAAACC\", \"id\": 86066}], \"edge\": [{\"from\": 86051, \"to\": 86052}, {\"from\": 86051, \"to\": 86053}, {\"from\": 86052, \"to\": 86054}, {\"from\": 86053, \"to\": 86054}, {\"from\": 86054, \"to\": 86055}, {\"from\": 86054, \"to\": 86056}, {\"from\": 86055, \"to\": 86057}, {\"from\": 86056, \"to\": 86057}, {\"from\": 86057, \"to\": 86058}, {\"from\": 86057, \"to\": 86059}, {\"from\": 86058, \"to\": 86060}, {\"from\": 86059, \"to\": 86060}, {\"from\": 86060, \"to\": 86061}, {\"from\": 86060, \"to\": 86062}, {\"from\": 86061, \"to\": 86063}, {\"from\": 86062, \"to\": 86063}, {\"from\": 86063, \"to\": 86064}, {\"from\": 86064, \"to\": 86065}, {\"from\": 86065, \"to\": 86066}]}";

                Graph fail_graph;
                json2pb(fail_graph, json_graph.c_str(), json_graph.size());

                VG graph;
                graph.extend(fail_graph);

                Aligner aligner;
                
                string read = "CTATTGTTTCTGATGGGGTCTCAGTGGTGCCCAAAGCAACAGATCCCCATGGCAGAGGGA";
                Alignment aln;
                aln.set_sequence(read);
                
                int band_padding = 1;
                aligner.align_global_banded(aln, graph.graph, band_padding);
                
                const Path& path = aln.path();
                
                // is a global alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                
            }
        }
        
        TEST_CASE( "Banded global aligner produces correct alignments with different graph structures",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner produces correct alignment when there are paths of different lengths in the graph" ) {
                
                VG graph;
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
            
            SECTION( "Banded global aligner produces correct alignment against a singleton graph",
                    "[alignment][banded][mapping]" ) {
                
                VG graph;
                
                Aligner aligner;
                
                Node* n0 = graph.create_node("C");
                
                string read = string("A");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_padding = 1;
                
                aligner.align_global_banded(aln, graph.graph, band_padding);
                
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
            
            SECTION( "Banded global aligner produces correct alignment when a node is masked" ) {
                
                VG graph;
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("CGCC");
                Node* n2 = graph.create_node("ATTA");
                
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n2);
                
                string read = string("AGTGATTA");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
                Node* n0 = graph.create_node("AGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                
                string read = string("AGTGC");
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 0;
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
        }
        
        TEST_CASE( "Banded global aligner produces correct alignments with permissive banding option",
                  "[alignment][banded][mapping]" ) {
            
            SECTION( "Banded global aligner runs with permissive banding" ) {
                
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
                Alignment aln;
                aln.set_sequence(read);
                
                int band_width = 1;
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width,
                                                                         true);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                QualAdjAligner aligner = QualAdjAligner();
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width,
                                                                         true);
                
                
                banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width,
                                                                         true,
                                                                         true);
                
                
                banded_aligner.align(aligner.adjusted_score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                BandedGlobalAligner banded_aligner = BandedGlobalAligner(aln,
                                                                         graph.graph,
                                                                         band_width,
                                                                         true,
                                                                         true);
                
                
                banded_aligner.align(aligner.adjusted_score_matrix, aligner.nt_table, aligner.gap_open,
                                     aligner.gap_extension);
                
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
                
                BandedGlobalAligner banded_aligner_full = BandedGlobalAligner(aln_full,
                                                                              graph.graph,
                                                                              band_width,
                                                                              true,
                                                                              true);
                
                BandedGlobalAligner banded_aligner_reduced = BandedGlobalAligner(aln_reduced,
                                                                                 graph.graph,
                                                                                 band_width,
                                                                                 true,
                                                                                 true);
                
                
                banded_aligner_full.align(aligner.adjusted_score_matrix, aligner.nt_table, aligner.scaled_gap_open,
                                          aligner.scaled_gap_extension);
                
                banded_aligner_reduced.align(aligner.adjusted_score_matrix, aligner.nt_table, aligner.scaled_gap_open,
                                             aligner.scaled_gap_extension);
                
                
                
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
                  "[alignment][banded][mapping]" ) {
            
            
            SECTION( "Banded global aligner returns alternate alignments in descending score order" ) {
                
                VG graph;
                
                Aligner aligner = Aligner();
                
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
                
                aligner.align_global_banded_multi(aln, multi_alns, graph.graph, max_multi_alns,
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
                
                Aligner aligner;
                
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
                
                aligner.align_global_banded_multi(aln, multi_alns, graph.graph, max_multi_alns,
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
                
                Aligner aligner;
                
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
                
                aligner.align_global_banded_multi(aln, multi_alns, graph.graph, max_multi_alns,
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
                
                Aligner aligner;
                
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
                
                aligner.align_global_banded_multi(aln, multi_alns, graph.graph, max_multi_alns,
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
                
                Aligner aligner;
                
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
                
                aligner.align_global_banded_multi(aln, multi_alns, graph.graph, max_multi_alns,
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
            
            
            
            SECTION( "Banded global aligner can identify an alternate alignment that branches from another alternate alignment inside a node sequence" ) {
                
                VG graph;
                
                Aligner aligner;
                
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
                
                aligner.align_global_banded_multi(aln, multi_alns, graph.graph, max_multi_alns,
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
            
            SECTION( "Global banded aligner does not produce duplicate alternate alignments" ) {
                
                VG graph;
                
                Aligner aligner;
                
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
                
                aligner.align_global_banded_multi(aln, multi_alns, graph.graph, max_multi_alns,
                                                  band_padding, permissive_banding);
                
                unordered_set<string> alns_seen;
                for (Alignment& alt_aln : multi_alns) {
                    string aln_string = hash_alignment(alt_aln);
                    
                    REQUIRE(alns_seen.count(aln_string) == 0);
                    alns_seen.insert(aln_string);
                }
            }
        }
    }
}






