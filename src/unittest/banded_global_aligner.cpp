//
//  banded_global_aligner.cpp
//  
//
//  Created by Jordan Eizenga on 7/26/16.
//
//

#include <stdio.h>
#include "catch.hpp"
#include "gssw_aligner.hpp"
#include "vg.hpp"
#include "banded_global_aligner.hpp"

using namespace google::protobuf;

namespace vg {
    namespace unittest {
        
        TEST_CASE( "Banded global alignment works when read matches exactly",
                   "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when there is a mismatch",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when there are paths of different lenths in the graph",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a node is masked",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a single base deletion",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a multi-base deletion",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a single base insertion",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a multi-base insertion",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a deletion across a node boundary",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a single base lead insertion",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a multi-base lead insertion",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a single base lead deletion",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there is a multi-base lead deletion",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when a there lead deletion across a node boundary",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when band is wider than read",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works when band one base wide",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works with permissive banding",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width,
                                                                                   true);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment with permissive banding can complete alignment outside of the padded band",
                   "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width,
                                                                                   true);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works with multiple source nodes",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works with multiple sink nodes",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works with multiple source and multiple sink nodes",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment works with a separated graph",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width);
            
            
            banded_aligner.align(aligner.score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
            
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment runs with quality adjusted alignments",
                  "[alignment][banded][mapping]" ) {
            
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
            BandedGlobalAlignmentGraph banded_aligner = BandedGlobalAlignmentGraph(aln,
                                                                                   graph.graph,
                                                                                   band_width,
                                                                                   true,
                                                                                   true);
            
            
            banded_aligner.align(aligner.adjusted_score_matrix, aligner.nt_table, aligner.gap_open,
                                 aligner.gap_extension);
            
            const Path& path = aln.path();
                        
            // is a global alignment
            REQUIRE(path.mapping(0).position().offset() == 0);
            
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
        
        TEST_CASE( "Banded global alignment increases score when adjusting for low quality mismatches",
                  "[alignment][banded][mapping]" ) {
            
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
            
            BandedGlobalAlignmentGraph banded_aligner_full = BandedGlobalAlignmentGraph(aln_full,
                                                                                        graph.graph,
                                                                                        band_width,
                                                                                        true,
                                                                                        true);
            
            BandedGlobalAlignmentGraph banded_aligner_reduced = BandedGlobalAlignmentGraph(aln_reduced,
                                                                                           graph.graph,
                                                                                           band_width,
                                                                                           true,
                                                                                           true);
            
            
            banded_aligner_full.align(aligner.adjusted_score_matrix, aligner.nt_table, aligner.gap_open,
                                      aligner.gap_extension);
            
            banded_aligner_reduced.align(aligner.adjusted_score_matrix, aligner.nt_table, aligner.gap_open,
                                         aligner.gap_extension);
            
            const Path& path_full = aln_full.path();
            const Path& path_reduced = aln_reduced.path();
            
            // is a global alignment
            REQUIRE(path_full.mapping(0).position().offset() == 0);
            REQUIRE(path_reduced.mapping(0).position().offset() == 0);
            
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
}
