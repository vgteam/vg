//
// pinned_alignment.cpp
//  
// Unit tests for pinned alignment and pinned multi-alignment functions in Aligner
//

#include <stdio.h>
#include <unordered_set>
#include "alignment.hpp"
#include "test_aligner.hpp"
#include "vg.hpp"
#include "path.hpp"
#include "catch.hpp"
#include "vg/io/json2pb.h"
#include "bdsg/hash_graph.hpp"

using namespace google::protobuf;
using namespace vg::io;

namespace vg {
    namespace unittest {
        
        TEST_CASE( "Pinned alignment produces correct alignments with different types of edits",
                  "[alignment][pinned][mapping]" ) {
            
            SECTION( "Pinned alignment produces correct alignment when read matches exactly") {
                
                VG graph;
                
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
                
                Node* pinned_node = n3;
                bool pin_left = false;
                                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
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
            
            SECTION( "Pinned alignment produces same alignment for an exact match regardless of left or right pinning") {
                
                VG graph;
                
                Node* n0 = graph.create_node("GGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("GGTGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
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
            
            SECTION( "Right-pinned alignment soft clips off the left end when there is a mismatch at final base") {
                
                VG graph;
                
                Node* n0 = graph.create_node("CGTG");
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
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, 0);
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).position().offset() == 1);
                
                REQUIRE(path.mapping(0).edit(0).to_length() == 1);
                REQUIRE(path.mapping(0).edit(0).sequence() == "A");
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 3);
                REQUIRE(path.mapping(0).edit(1).to_length() == 3);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Left-pinned alignment soft clips off the right end when there is a mismatch at final base") {
                
                VG graph;
                
                Node* n0 = graph.create_node("TGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGC");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("TGTGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
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
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 5);
                REQUIRE(path.mapping(2).edit(0).to_length() == 5);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "T");
            }
            
            SECTION( "Right-pinned alignment accepts mismatch on the right end when there is a mismatch at final base") {
                
                VG graph;
                
                Node* n0 = graph.create_node("CGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGC");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CGTGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
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
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 5);
                REQUIRE(path.mapping(2).edit(0).to_length() == 5);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 1);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "T");
            }
            
            SECTION( "Left-pinned alignment accepts mismatch on the right end when there is a mismatch at final base") {
                
                VG graph;
                
                Node* n0 = graph.create_node("CGTG");
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
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 1);
                REQUIRE(path.mapping(0).edit(0).to_length() == 1);
                REQUIRE(path.mapping(0).edit(0).sequence() == "A");
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 3);
                REQUIRE(path.mapping(0).edit(1).to_length() == 3);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Right-pinned alignment accepts a match of N to N on the right end") {
                
                VG graph;
                
                Node* n0 = graph.create_node("CGTG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGANNN");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CGTGCTGANNN");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has correct mapping characteristics
                // TODO: the N-N matches are produced as separate edits from the normal matches.
                REQUIRE(mapping_from_length(path.mapping(0)) == 4);
                REQUIRE(mapping_from_length(path.mapping(0)) == 4);
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(mapping_from_length(path.mapping(2)) == 6);
                REQUIRE(mapping_from_length(path.mapping(2)) == 6);
            }
            
            SECTION( "Left-pinned alignment accepts a match of N to N on the left end") {
                
                VG graph;
                
                Node* n0 = graph.create_node("NNNG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("NNNGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has correct mapping characteristics
                REQUIRE(mapping_from_length(path.mapping(0)) == 4);
                REQUIRE(mapping_from_length(path.mapping(0)) == 4);
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(mapping_from_length(path.mapping(2)) == 6);
                REQUIRE(mapping_from_length(path.mapping(2)) == 6);
            }
            
            SECTION( "Pinned alignment produces correct alignment when there is a mismatch" ) {
                
                VG graph;
                
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
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
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
            
            SECTION( "Pinned alignment produces correct alignment when a there is a deletion"  ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGATG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 8);
                REQUIRE(path.mapping(0).edit(0).to_length() == 8);
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
            
            SECTION( "Pinned alignment produces correct alignment when a there is an insertion" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGATGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 8);
                REQUIRE(path.mapping(0).edit(0).to_length() == 8);
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
            
            SECTION( "Pinned alignment produces correct alignment when a there is a deletion across a node boundary" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAAACCCAGC");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("AT");
                Node* n3 = graph.create_node("TGAAGTAA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAAACCCAGTTGAAGTAA");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 3);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
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
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 8);
                REQUIRE(path.mapping(2).edit(0).to_length() == 8);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Pinned alignment produces correct alignment when a there is an N match" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAAACCCAGC");
                Node* n1 = graph.create_node("CA");
                Node* n2 = graph.create_node("TT");
                Node* n3 = graph.create_node("TGAAGTAA");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAAACCCAGCNATGAAGTAA");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 10);
                REQUIRE(path.mapping(0).edit(0).to_length() == 10);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence() == "N");
                
                REQUIRE(path.mapping(1).edit(1).from_length() == 1);
                REQUIRE(path.mapping(1).edit(1).to_length() == 1);
                REQUIRE(path.mapping(1).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 8);
                REQUIRE(path.mapping(2).edit(0).to_length() == 8);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Pinned alignment produces a correct right-pinned null alignment when there is no match" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAAAA");
                
                string read = string("CCC");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned null alignment
                REQUIRE(path.mapping(0).position().offset() == n0->sequence().length());
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 0);
                REQUIRE(path.mapping(0).edit(0).to_length() == 3);
                REQUIRE(path.mapping(0).edit(0).sequence() == aln.sequence());
            }
            
            SECTION( "Pinned alignment produces a correct left-pinned null alignment when there is no match" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAAAA");
                
                string read = string("CCC");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned null alignment
                REQUIRE(path.mapping(0).position().offset() == 0);
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 0);
                REQUIRE(path.mapping(0).edit(0).to_length() == 3);
                REQUIRE(path.mapping(0).edit(0).sequence() == aln.sequence());
            }
        }
        
        TEST_CASE( "Pinned alignment produces correct alignments in edge cases involving the pinned end",
                  "[alignment][pinned][mapping]" ) {
        
            SECTION( "Right-pinned alignment produces correct alignment when a there is an insertion on the pinned end" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGGCTGAAGTA");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 0);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "A");
            }
            
            SECTION( "Left-pinned alignment produces correct alignment when a there is an insertion on the pinned end" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CAAACCCAGGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 0);
                REQUIRE(path.mapping(0).edit(0).to_length() == 1);
                REQUIRE(path.mapping(0).edit(0).sequence() == "C");
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 9);
                REQUIRE(path.mapping(0).edit(1).to_length() == 9);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
            }
            
            SECTION( "Right-pinned alignment produces correct alignment when there is a deletion on the pinned end" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGACTGGATAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("TGACTGGATAA");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 11);
                REQUIRE(path.mapping(0).edit(0).to_length() == 11);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 2);
                REQUIRE(path.mapping(0).edit(1).to_length() == 0);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
            }
            
            SECTION( "Right-pinned alignment produces correct alignment when there is a deletion of multiple nodes on the pinned end" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCTGAACGTAGAGGCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("C");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n1, n2);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCTGAACGTAGAGGCAGG");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 3);
                REQUIRE(path.mapping(3).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 21);
                REQUIRE(path.mapping(0).edit(0).to_length() == 21);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 0);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 1);
                REQUIRE(path.mapping(2).edit(0).to_length() == 0);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(3).edit(0).from_length() == 1);
                REQUIRE(path.mapping(3).edit(0).to_length() == 0);
                REQUIRE(path.mapping(3).edit(0).sequence().empty());
            }
            
            SECTION( "Left-pinned alignment produces correct alignment when there is a deletion on the pinned end" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("TGACTGGATAAGT");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("AAACCCAGG");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("ACTGGATAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 2);
                REQUIRE(path.mapping(0).edit(0).to_length() == 0);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 11);
                REQUIRE(path.mapping(0).edit(1).to_length() == 11);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
            }
            
            SECTION( "Left-pinned alignment produces correct alignment when there is a deletion of multiple nodes on the pinned end" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("A");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("CAACCTGAACGTAGAGGCAGG");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n1, n2);
                graph.create_edge(n2, n3);
                
                string read = string("CAACCTGAACGTAGAGGCAGG");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 3);
                REQUIRE(path.mapping(3).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 1);
                REQUIRE(path.mapping(0).edit(0).to_length() == 0);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 0);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 1);
                REQUIRE(path.mapping(2).edit(0).to_length() == 0);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(3).edit(0).from_length() == 21);
                REQUIRE(path.mapping(3).edit(0).to_length() == 21);
                REQUIRE(path.mapping(3).edit(0).sequence().empty());
            }
            
            SECTION( "Right-pinned alignment produces correct alignment when there is a deletion on the pinned end followed by an N match" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGACTGGATAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("TGACTGGATAN");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 10);
                REQUIRE(path.mapping(0).edit(0).to_length() == 10);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                bool option_1 = true;
                bool option_2 = true;
                
                option_1 = option_1 && path.mapping(0).edit(1).from_length() == 1;
                option_1 = option_1 && path.mapping(0).edit(1).to_length() == 1;
                option_1 = option_1 && path.mapping(0).edit(1).sequence() == "N";
                
                option_1 = option_1 && path.mapping(0).edit(2).from_length() == 2;
                option_1 = option_1 && path.mapping(0).edit(2).to_length() == 0;
                option_1 = option_1 && path.mapping(0).edit(2).sequence().empty();
                
                option_2 = option_2 && path.mapping(0).edit(1).from_length() == 2;
                option_2 = option_2 && path.mapping(0).edit(1).to_length() == 0;
                option_2 = option_2 && path.mapping(0).edit(1).sequence().empty();
                
                option_2 = option_2 && path.mapping(0).edit(2).from_length() == 1;
                option_2 = option_2 && path.mapping(0).edit(2).to_length() == 1;
                option_2 = option_2 && path.mapping(0).edit(2).sequence() == "N";
                
                bool found_correct_option = option_1 || option_2;
                REQUIRE(found_correct_option);
                
            }
            
            SECTION( "Left-pinned alignment produces correct alignment when there is a deletion on the pinned end followed by an N match" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("TGACTGGATAAGT");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("AAACCCAGG");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("NCTGGATAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                
                // has corrects edits
                
                bool option_1 = true;
                bool option_2 = true;
                
                option_1 = option_1 && path.mapping(0).edit(0).from_length() == 1;
                option_1 = option_1 && path.mapping(0).edit(0).to_length() == 1;
                option_1 = option_1 && path.mapping(0).edit(0).sequence() == "N";
                
                option_1 = option_1 && path.mapping(0).edit(1).from_length() == 2;
                option_1 = option_1 && path.mapping(0).edit(1).to_length() == 0;
                option_1 = option_1 && path.mapping(0).edit(1).sequence().empty();
                
                option_2 = option_2 && path.mapping(0).edit(0).from_length() == 2;
                option_2 = option_2 && path.mapping(0).edit(0).to_length() == 0;
                option_2 = option_2 && path.mapping(0).edit(0).sequence().empty();
                
                option_2 = option_2 && path.mapping(0).edit(1).from_length() == 1;
                option_2 = option_2 && path.mapping(0).edit(1).to_length() == 1;
                option_2 = option_2 && path.mapping(0).edit(1).sequence() == "N";
                
                bool found_correct_option = option_1 || option_2;
                REQUIRE(found_correct_option);
                
            }
            
            SECTION( "Pinned alignment can correctly choose between multiple sink nodes" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("T");
                Node* n4 = graph.create_node("G");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n0, n3);
                graph.create_edge(n0, n4);
                
                string read = string("AAACCCAGGA");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n2;
                bool pin_left = false;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 3);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
            }
        }
        
        TEST_CASE ( "Pinned alignment utilizes the pinned end full-length alignment bonus correctly",
                   "[alignment][pinned][mapping]" ) {
        
            SECTION( "Right-pinned alignment will take a mismatch near the left tail to gain a full-length alignment bonus" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AATCCCAGGCTGAAGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                int8_t full_length_bonus = 3;
                                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, full_length_bonus);
                const Aligner& aligner = *aligner_source.get_regular_aligner();

                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 2);
                REQUIRE(path.mapping(0).edit(0).to_length() == 2);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 1);
                REQUIRE(path.mapping(0).edit(1).to_length() == 1);
                REQUIRE(path.mapping(0).edit(1).sequence() == "T");
                
                REQUIRE(path.mapping(0).edit(2).from_length() == 6);
                REQUIRE(path.mapping(0).edit(2).to_length() == 6);
                REQUIRE(path.mapping(0).edit(2).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 6);
                REQUIRE(path.mapping(2).edit(0).to_length() == 6);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                // score is correct
                REQUIRE(aln.score() == 2 - 4 + 6 + 1 + 6 + full_length_bonus);
            }
            
            SECTION( "Left-pinned alignment will take a mismatch near the right tail to gain a full-length alignment bonus" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGGCTGACGT");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                int8_t full_length_bonus = 3;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, full_length_bonus);
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 3);
                REQUIRE(path.mapping(2).edit(0).to_length() == 3);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 1);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "C");
                
                REQUIRE(path.mapping(2).edit(2).from_length() == 2);
                REQUIRE(path.mapping(2).edit(2).to_length() == 2);
                REQUIRE(path.mapping(2).edit(2).sequence().empty());
                
                // score is correct
                REQUIRE(aln.score() == 9 + 1 + 3 - 4 + 2 + full_length_bonus);
            }
            
            SECTION( "Pinned alignment will take a mismatch on the tail to gain a full-length alignment bonus if it is greater than the mismatch penalty" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGGCTGAAGG");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                int8_t full_length_bonus = 5;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, full_length_bonus);
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                aligner.align_pinned(aln, graph, pin_left);
                                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 5);
                REQUIRE(path.mapping(2).edit(0).to_length() == 5);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 1);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "G");
                
                // score is correct
                REQUIRE(aln.score() == 9 + 1 + 5 - 4 + full_length_bonus);
            }
            
            SECTION( "Pinned alignment will not take a mismatch on the tail to gain a full-length alignment bonus if it is less than the mismatch penalty" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGGCTGAAGG");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                int8_t full_length_bonus = 2;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, full_length_bonus);
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 5);
                REQUIRE(path.mapping(2).edit(0).to_length() == 5);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 0);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "G");
                
                // score is correct
                REQUIRE(aln.score() == 9 + 1 + 5);
            }
        }
        
        TEST_CASE( "Pinned alignment behaves correctly with quality adjusted alignments",
                  "[alignment][pinned][mapping]" ) {
        
            SECTION( "Quality adjusted left-pinned alignment runs" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGC");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGGCTGAAGG");
                string qual = string("HHHHHHHHHHHHHHHH");
                
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                int8_t full_length_bonus = 5;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, full_length_bonus);
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 5);
                REQUIRE(path.mapping(2).edit(0).to_length() == 5);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 1);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "G");
            }
            
            SECTION( "Quality adjusted right-pinned alignment runs" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGC");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGGCTGAAGG");
                string qual = string("HHHHHHHHHHHHHHHH");

                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                int8_t full_length_bonus = 5;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, full_length_bonus);
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 5);
                REQUIRE(path.mapping(2).edit(0).to_length() == 5);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 1);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "G");
            }
            
            SECTION( "Quality adjusted pinned alignment runs with non-trivial quality scores" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGC");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("AAACCCAGGCTGAAGG");
                string qual = string("HHHDDD<<<88800((");
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                int8_t full_length_bonus = 5;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, full_length_bonus);
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 9);
                REQUIRE(path.mapping(0).edit(0).to_length() == 9);
                REQUIRE(path.mapping(0).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 5);
                REQUIRE(path.mapping(2).edit(0).to_length() == 5);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 1);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "G");
            }
            
            SECTION( "Quality adjusted pinned alignment runs with all features activated" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("AAACCCAGG");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("TGAAGT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = string("CAACCCAGGCTGAAGG");
                string qual = string("HHHDDD<<<88800HH");
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                int8_t full_length_bonus = 3;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, full_length_bonus);
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                aligner.align_pinned(aln, graph, pin_left);
                
                const Path& path = aln.path();
                
                // is a pinned alignment
                if (pin_left) {
                    REQUIRE(path.mapping(0).position().offset() == 0);
                    REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                }
                else {
                    REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                    REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                }
                
                // follows correct path
                REQUIRE(path.mapping(0).position().node_id() == 1);
                REQUIRE(path.mapping(1).position().node_id() == 2);
                REQUIRE(path.mapping(2).position().node_id() == 4);
                
                // has corrects edits
                REQUIRE(path.mapping(0).edit(0).from_length() == 1);
                REQUIRE(path.mapping(0).edit(0).to_length() == 1);
                REQUIRE(path.mapping(0).edit(0).sequence() == "C");
                
                REQUIRE(path.mapping(0).edit(1).from_length() == 8);
                REQUIRE(path.mapping(0).edit(1).to_length() == 8);
                REQUIRE(path.mapping(0).edit(1).sequence().empty());
                
                REQUIRE(path.mapping(1).edit(0).from_length() == 1);
                REQUIRE(path.mapping(1).edit(0).to_length() == 1);
                REQUIRE(path.mapping(1).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(0).from_length() == 5);
                REQUIRE(path.mapping(2).edit(0).to_length() == 5);
                REQUIRE(path.mapping(2).edit(0).sequence().empty());
                
                REQUIRE(path.mapping(2).edit(1).from_length() == 0);
                REQUIRE(path.mapping(2).edit(1).to_length() == 1);
                REQUIRE(path.mapping(2).edit(1).sequence() == "G");
            }
        }
        
        TEST_CASE( "Pinned multi-alignment produces correct alternate alignments in different traceback scenarios",
                  "[alignment][multialignment][pinned][mapping]" ) {
            
            SECTION( "Pinned multi-alignment returns alignments in descending score order" ) {
                
                VG graph;
                
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
                
                Node* pinned_node = n3;
                bool pin_left = false;
                int max_multi_alns = 20;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
                int64_t prev_score = numeric_limits<int64_t>::max();
                for (Alignment& alt_aln : multi_alns) {
                    const Path& path = alt_aln.path();
                    
                    // is a pinned alignment
                    if (pin_left) {
                        REQUIRE(path.mapping(0).position().offset() == 0);
                        REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                    }
                    else {
                        REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                        REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                    }
                    
                    // score is in descending order
                    REQUIRE(aln.score() <= prev_score);
                    prev_score = aln.score();
                }
            }
            
            SECTION( "Pinned multi-alignment stores the optimal alignment in both the main Alignment object and the first position in the return vector" ) {
                
                VG graph;
                
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
                
                Node* pinned_node = n3;
                bool pin_left = false;
                int max_multi_alns = 2;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
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
            
            SECTION( "Pinned multi-alignment identifies both optimal alignments when there is a deletion in a homodimer" ) {
                
                VG graph;
                
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
                
                Node* pinned_node = n0;
                bool pin_left = true;
                int max_multi_alns = 2;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    const Path& path = alt_aln.path();
                    
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    
                    // is a pinned alignment
                    if (pin_left) {
                        REQUIRE(path.mapping(0).position().offset() == 0);
                        REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                    }
                    else {
                        REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                        REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                    }
                    
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
            
            SECTION( "Pinned multi-alignment can identify alternate alignments that follow a different node sequence than the optimal alignment" ) {
                
                VG graph;
                
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
                
                Node* pinned_node = n0;
                bool pin_left = true;
                int max_multi_alns = 20;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
                bool took_alternate_path = false;
                for (Alignment& alt_aln : multi_alns) {
                    const Path& path = alt_aln.path();
                    
                    // is a pinned alignment
                    if (pin_left) {
                        REQUIRE(path.mapping(0).position().offset() == 0);
                        REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                    }
                    else {
                        REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                        REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                    }
                    
                    if (path.mapping(1).position().node_id() != aln.path().mapping(1).position().node_id()) {
                        took_alternate_path = true;
                        break;
                    }
                }
                
                REQUIRE(took_alternate_path);
            }
            
            SECTION( "Pinned multi-alignment will not return alternates when none score positively" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("CA");
                
                string read = string("A");
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n0;
                bool pin_left = false;
                int max_multi_alns = 100;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
                REQUIRE(multi_alns.size() == 1);
            }
            
            SECTION( "Pinned multi-alignment can identify an alternate alignment that branches from another alternate alignment at a node boundary" ) {
                
                VG graph;
                
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
                
                Node* pinned_node = n6;
                bool pin_left = false;
                int max_multi_alns = 3;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    
                    const Path& path = alt_aln.path();
                    
                    // is a pinned alignment
                    if (pin_left) {
                        REQUIRE(path.mapping(0).position().offset() == 0);
                        REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                    }
                    else {
                        REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                        REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                    }
                    
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
            
            SECTION( "Pinned multi-alignment can identify an alternate alignment that branches from another alternate alignment inside a node sequence" ) {
                
                VG graph;
                
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
                
                Node* pinned_node = n3;
                bool pin_left = false;
                int max_multi_alns = 10;
                
                TestAligner aligner_source;
                const Aligner& aligner = *aligner_source.get_regular_aligner();
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
                bool found_first_opt = false;
                bool found_second_opt = false;
                for (Alignment& alt_aln : multi_alns) {
                    bool is_first_opt = true;
                    bool is_second_opt = true;
                    
                    const Path& path = alt_aln.path();
                    
                    // is a pinned alignment
                    if (pin_left) {
                        REQUIRE(path.mapping(0).position().offset() == 0);
                        REQUIRE(path.mapping(0).position().node_id() == pinned_node->id());
                    }
                    else {
                        REQUIRE(mapping_from_length(path.mapping(path.mapping_size() - 1)) == pinned_node->sequence().length());
                        REQUIRE(path.mapping(path.mapping_size() - 1).position().node_id() == pinned_node->id());
                    }
                    
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
            
            SECTION( "Pinned multi-alignment can take alternate matches over node boundaries without incurring a known bug" ) {
                
                VG graph;
                
                Node* n0 = graph.create_node("T");
                Node* n1 = graph.create_node("C");
                Node* n2 = graph.create_node("A");
                Node* n3 = graph.create_node("CCCTGCTAGTCTGGAGTTGATCAAGGAACCTGTCT");
                
                graph.create_edge(n0, n1);
                graph.create_edge(n0, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                
                string read = "CCCGG";
                string qual = "<<<''";
                Alignment aln;
                aln.set_sequence(read);
                aln.set_quality(qual);
                alignment_quality_char_to_short(aln);
                
                Node* pinned_node = n0;
                bool pin_left = true;
                int max_multi_alns = 100;
                
                TestAligner aligner_source;
                const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
                
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
                for (Alignment& alt_aln : multi_alns) {
                    for (size_t i = 0; i < alt_aln.path().mapping_size(); i++) {
                        REQUIRE(alt_aln.path().mapping(i).position().offset() >= 0);
                        REQUIRE(path_to_length(alt_aln.path()) == read.size());
                    }
                }
            }
            
            
            SECTION( "Pinned multi-alignment does not produce duplicate alignments" ) {
                
                VG graph;
                
                TestAligner aligner_source;
                aligner_source.set_alignment_scores(1, 4, 6, 1, 0);
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
                
                string read = "CCCCCCCTCCCCCCCCCCTCCCCCCCCCCGACCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCACCCCCCCCCCTCCCACCCCCCCCCCCCGCCCCCCCCCCGCCCCCCCCC";
                Alignment aln;
                aln.set_sequence(read);
                
                Node* pinned_node = n3;
                bool pin_left = false;
                int max_multi_alns = 5000;
                
                vector<Alignment> multi_alns;
                aligner.align_pinned_multi(aln, multi_alns, graph, pin_left, max_multi_alns);
                
                unordered_set<string> alns_seen;
                for (Alignment& alt_aln : multi_alns) {
                    string aln_string = hash_alignment(alt_aln);
                    
                    REQUIRE(alns_seen.count(aln_string) == 0);
                    alns_seen.insert(aln_string);
                }
            }
        }
        
        TEST_CASE("Quality adjusted alignment scores full length bonuses correctly",
                  "[alignment][pinned][mapping]" ) {
            
            bdsg::HashGraph graph;
            
            handle_t h = graph.create_handle("ACGTAGTCTGAA");
            
            Alignment aln1;
            aln1.set_sequence("ACGT");
            aln1.set_quality("HHH#");
            alignment_quality_char_to_short(aln1);
            
            Alignment aln2;
            aln2.set_sequence("TGAA");
            aln2.set_quality("#HHH");
            alignment_quality_char_to_short(aln2);
            
            TestAligner aligner_source;
            const QualAdjAligner& aligner = *aligner_source.get_qual_adj_aligner();
            aligner.align_pinned(aln1, graph, true, false);
            aligner.align_pinned(aln2, graph, false, false);
            
            REQUIRE(aln1.score() == 3);
            REQUIRE(aln2.score() == 3);
        }

        TEST_CASE("Pinned alignment doesn't produce invalid alignments",
                  "[alignment][pinned][mapping]" ) {
            
            std::string read_string = "AAGTGGACTGCATTAGAGGGAGCAGATTTGAAACACTGTTTTTGTGGAATTTGCAAGTGGAGATTTCAAGCGCTTTGGGGCCAAAGGCAGAGAAAGGAAATATCTTCGTATAAAAACTAGACAGAATCATTCTCAGAAACTGCTGCGTGATGTGCGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCGGTTTGGAAACACTCTGTTTGTAAAGTCTGCACGTGGATATTTTGACCACTTAGAGGCCTTCGTTGGAAACGGGTTTTTTGCATGTAAGGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTTGTGTGCATTCAACTCACAGAGTTGAACGTTCCCTTAGACAGAGCAGATTTGAAACACTCTATTTGTGCAATTTGAAAGTGTAGATTTCAAGCGCTTTAAGGTCAGTGGCAGAAAAGGAAATATCATCGTTTCAAAACTAGACAGAATCATTCCCACAAACTGCGTTGTGATGTGTTCGTTCAACTCACAGAGTTTAACCTTTCTGTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGTAAGTGGATATTCTGACATCTAGTGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTATTTCCTTGTTTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAACGGTAGAAAAGGAAATATCTTCGTAATAAAGACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACTTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTGAGTCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATAGTAGAAAAGGAAATATCTTCGTAGAAAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGTGTTCAACTCACAGAGTTTAACCTTACTTTTCATAGAGCAGTTAGTAAACACTCTGTTTATAAAGTCTGCAAGTGGATATTCAGACCCCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTATGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACTTTCATTTAGAGGGAGCAGATTTGAAACACTGTTTTTGTGGAATTTGCAAGTGGAGATTTCAAGCGCTTTGGGGCCAAAGGCAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAAGATTCTCAGAAACTGCTGCGTGATGTGCGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCTAAGGTTTGGAAACACTCTGTTTGTAAAGTCTGCACGTGGATATTTTGACCACTTAGAGGCCTTCGTTGGAAACGGGTTTTTTGCATGTAAGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTTGTGTGCATTCAACTCACAGAGTTGAACGTTCCCTTAGACAGAGCAGATTTGAAACACTCTATTTGTGCAATTTGCAAGTGTAGATTTCAAGCGCTTTAAGGTCAGTGGCAGAAAAGGATATATCTTCGTTTCAAAACTAGACAGAATCATTCCCACAAACTGCGTTGTGATGTGTTCGTTCAACTCACAGAGTTTAACCTTTCTGTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGTAAGTGGATATTCTGACATCTTGTGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACAATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAAAAGGAAATATCTTCGTATAAAGACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACTTTTCTATTCATAGACCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAGGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACACACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATAGTAGAAAAGGGAATATCTTCGTAGAAAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGTGTTCAACTTACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGTAAACACTCTGTTTATAAAGTCTGCAAGTGGATATTCAGACCCCTTTGAGCCTTCGTTGGAAACGGGATTTCTTCATATTATGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACTTTCATTTAGAGGGAGCAGATTTGAAACACTGTTTTTGTGGAATTTGCAAGTGGAGATTTCAAGCGCTTTGGGGCCAAAGGCAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAATCATTCTCAGAAACTGCTGCGTGATGTGCGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCGGTTTGGAAACACTCTGTTTGTAAAGTCTGCACGTGGATATTTTGACCACTTAGAGGCCTTCGTTGGAAACGGGTTTTTTGCATGTAAGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTTGTGTGCATTCAACTCACAGAGTTGAACGTTCCCTTAGACAGAGCAGATTTGAAACACTCTATTTGTGCAATTTGCAAGTGTAGATTTCAAGCGCTTTAAGGTCAGTGGCAGAAAAGGATATATCTTCGTTTCAAAACTAGACAGAATCATTCCCACAAACTGCGTTGTGATGTGTTCGTTCAACTCACAGAGTTTAACCTTTCTGTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGTAAGTGGATATTCCGACATCTTGTGGCCTTCGTTGGAAACGGGATTTCTTCATATCTGGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACAATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTGAATGGTAGAAAAGGAAATATCTTCGTATAAAGATAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAGCTTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAGGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGACACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACACACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGAGGAGATTTAGCCGCTTTGAGGTCAATAGTAGAAAAGGGAATATCTTCGTAGAAAAACTAGACAGAATGATTCTCAGAAACTCTTTGTGATGTGTGTGTTCAACTTACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGTAAACACTCTGTTTATAAAGTCTGCAAGTGGATATTCAGACCCCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTATGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACTTTCATTTAGAGGGAGCAGATTTGAAACACTGTTTTTGTGGAATTTGCAAGTGGAGATTTCAAGCGCTTTGGGGCCAAAGGCAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAATCATTCTCAGAAACTGCTGCGTGATGTGCGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCGGTTTGGAAACACTCTGTTTGTAAAGTCTGCACGTGGATATTTTGACCACTTAGAGGCCTTCGTTGGAAACGGGTTTTTTGCATGTAAGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTTGTGTGCATTCAACTCACAGAGTTGAACGTTCCCTTAGACAGAGCAGATTTGAAACACTCTATTTTGTGCAATTTGCAAGTGTAGATTTCAAGCGCTTTAAGGTCAGTGGCAGAAAAGGAAATATCTTCGTTTCAAAACTAGACAGAATCATTCCCACAAACTGCGTTGTGATGTGTTCGTTCAACTCACAGAGTTTAACCTTTCTGTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGTAAGTGATATTCTGACATCTTGTGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTGTGTGTATTCAACTCACAGAGTTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAAAAGGAAATATCTTCGTATAAAGACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACTTTTCTTTTCATAAAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACAGAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAATAGGAAATATCTTCCTATAGTAACTGACAGAATGATCTCAGAAGACTCCTTTGGTGATGGTGGCGTTCAACTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAAACACCCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATAGTAGAAAAGGAAATATCTTCGTAGAAAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGTGTTCAACTCACAGAGTTTAACCTGCTTTTCATAGAGCAGTTAGTAAACACTCTGTTTATAAAGTCTGCAAGTGGATATTTAGACCCCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTATGCTAGACAGAAGACGTTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACTTTCATTTAGAGAGAGCAGATTTGAAACACTGTTTTTGTGGAATTTGCAAGTGGAGATTTCAAGCGCATTGGGGCCAAAGGCAGAAAAGGAAATATCTTCGTATAAAAACTACACAGAATCATTCTCAGAAACTGCTGCGTGATGTGTGCGTTCAACTCTCTGAGTTTAACTTTTCTTTTCATTCAGCGGTTTGGAAACACTCTGTTTGTAAAGTCTGCACGTAGATATTTTGACCACTTAGAGGCCTCGTTGGAAACGGGTTTTTTTCATGTAAGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTTGTGTGCATTCAACTCACAAGAGAACGTTCCCTTAGACAGAGCAGATTTGAAACACTCTATTTGTGCAATTGGCAAGTGTAGATTTCAAACGCTTTAAGGTCAAAGGCAGAAAAGGAAATATCTTCGTTTCAAAACTAGACAGAATCATTCCCACAAACTGCGTTGTGATGTGTTCGTTCAACTCACAGAGTTTAACCTTTCTGTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGTAAGTGGATATTCTGACATCTTGTGGCCTTCGTTGGAAACGGAATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACAGAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTCGAGGTCAATGGTAGAATAGGT";
            std::string graph_json = R"(

{
  "node": [
    {
      "id": "56",
      "sequence": "C"
    },
    {
      "id": "35",
      "sequence": "AAAAACTAGACAGAATCATTCTCAGAAACTGCTGCGTGATGTGTGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCGGTTTGGAAACACTCTGTTTGTAAAGTCTGCACGTGGATATTTTGACCACTTAGAGGCCTTCGTTGGAAACGGGTTTTTTTCATGTAAGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTTGTGTGCATTCAACTCACAGAGTTGAACGTTCCCTTAGACAGAGCAGATTTGAAACACTCTATTTGTGCAATTGGCAAGTGTAGATTTCAAGCGCTTTAA"
    },
    {
      "id": "60",
      "sequence": "AGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAAAAGGAAATATCTTCGTATAAAGACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACTTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGTATATCCAGATCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTATGC"
    },
    {
      "id": "67",
      "sequence": "TCT"
    },
    {
      "id": "73",
      "sequence": "AAACACTCTGTTT"
    },
    {
      "id": "115",
      "sequence": "C"
    },
    {
      "id": "112",
      "sequence": "GA"
    },
    {
      "id": "86",
      "sequence": "CCTTCGTTGGAAAC"
    },
    {
      "id": "168",
      "sequence": "A"
    },
    {
      "id": "12",
      "sequence": "G"
    },
    {
      "id": "75",
      "sequence": "TAAAGTCTGCA"
    },
    {
      "id": "23",
      "sequence": "CTTTCTTTTCATAGAGCAGTTAGGAGACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTATGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACTTTCATTTAGAGAGAGCAGATTGCAAACACTGTTTTTGTGGAATTTGCAAGTGGAGATTTCAAGCGCTTTGGGGCCAAAGGCAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAATCATTCTCAGAAACTGCAGCGTGATGTGTGCATTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCTGTTTGGAAACACTCTGTTTGTAAAGTCTGCACGTGGATATTTTGACCACTTAGAG"
    },
    {
      "id": "111",
      "sequence": "A"
    },
    {
      "id": "41",
      "sequence": "CTTCGT"
    },
    {
      "id": "68",
      "sequence": "T"
    },
    {
      "id": "82",
      "sequence": "CTT"
    },
    {
      "id": "130",
      "sequence": "GGCAGAAAAGGAAATATCTTCGT"
    },
    {
      "id": "125",
      "sequence": "AA"
    },
    {
      "id": "77",
      "sequence": "GTGGATATT"
    },
    {
      "id": "172",
      "sequence": "T"
    },
    {
      "id": "71",
      "sequence": "G"
    },
    {
      "id": "66",
      "sequence": "GCTAGACAGAAGAATTC"
    },
    {
      "id": "103",
      "sequence": "C"
    },
    {
      "id": "59",
      "sequence": "G"
    },
    {
      "id": "26",
      "sequence": "A"
    },
    {
      "id": "127",
      "sequence": "T"
    },
    {
      "id": "116",
      "sequence": "T"
    },
    {
      "id": "100",
      "sequence": "GCTAGACAGAAGAATTC"
    },
    {
      "id": "79",
      "sequence": "GAC"
    },
    {
      "id": "141",
      "sequence": "T"
    },
    {
      "id": "135",
      "sequence": "C"
    },
    {
      "id": "138",
      "sequence": "AAACTGC"
    },
    {
      "id": "107",
      "sequence": "G"
    },
    {
      "id": "46",
      "sequence": "C"
    },
    {
      "id": "57",
      "sequence": "GTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACAGAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTCGAGGTCAATGGTAGAATAGGTAATATCTTCCTATAGAAACTAGACAGAATAATTCTCAGAAACTC"
    },
    {
      "id": "152",
      "sequence": "AGC"
    },
    {
      "id": "170",
      "sequence": "C"
    },
    {
      "id": "129",
      "sequence": "T"
    },
    {
      "id": "78",
      "sequence": "TT"
    },
    {
      "id": "133",
      "sequence": "C"
    },
    {
      "id": "72",
      "sequence": "G"
    },
    {
      "id": "1",
      "sequence": ""
    },
    {
      "id": "137",
      "sequence": "C"
    },
    {
      "id": "22",
      "sequence": "G"
    },
    {
      "id": "154",
      "sequence": "GTT"
    },
    {
      "id": "33",
      "sequence": "GTCAATGGCAGAAAAGGAAATATCTTCGTTTCAAAACTAGACAGAATCATTCCCACAAACTGCGTTGTGATGTGTTCGTTCAACTCACAGAGTTGAACCTTTCTGTTCATAGAGCAGTTAGGAAACATTCTGTTTGTAAAGTCTGTAAGTGGATATTCTCACATCTTGTGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTCACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACA"
    },
    {
      "id": "40",
      "sequence": "G"
    },
    {
      "id": "113",
      "sequence": "C"
    },
    {
      "id": "165",
      "sequence": "G"
    },
    {
      "id": "142",
      "sequence": "C"
    },
    {
      "id": "5",
      "sequence": "TTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTT"
    },
    {
      "id": "55",
      "sequence": "TTTGTGATGTGTGCGTTCAACTCACAAAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGG"
    },
    {
      "id": "114",
      "sequence": "AGAGCAGATTTGAAACACT"
    },
    {
      "id": "136",
      "sequence": "CA"
    },
    {
      "id": "117",
      "sequence": "A"
    },
    {
      "id": "45",
      "sequence": "AAATATCTTCCTATAGAAACTAGACAGAAAGATTCTCATAAACTCCTTTGTGATGTGTGTGTTCAACTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAGACACTCTGTTTGTAAAGTCTGCAAGT"
    },
    {
      "id": "145",
      "sequence": "A"
    },
    {
      "id": "158",
      "sequence": "A"
    },
    {
      "id": "28",
      "sequence": "C"
    },
    {
      "id": "148",
      "sequence": "TTTCT"
    },
    {
      "id": "92",
      "sequence": "T"
    },
    {
      "id": "36",
      "sequence": "A"
    },
    {
      "id": "118",
      "sequence": "TTTGTG"
    },
    {
      "id": "162",
      "sequence": "TGAC"
    },
    {
      "id": "84",
      "sequence": "GA"
    },
    {
      "id": "7",
      "sequence": "AGAGCAG"
    },
    {
      "id": "25",
      "sequence": "TTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAATAGGAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCATAAACTCCTTTGTGATGAGTGCGTTCAACTCACAGAGTTTAA"
    },
    {
      "id": "95",
      "sequence": "G"
    },
    {
      "id": "93",
      "sequence": "G"
    },
    {
      "id": "18",
      "sequence": "G"
    },
    {
      "id": "147",
      "sequence": "C"
    },
    {
      "id": "157",
      "sequence": "T"
    },
    {
      "id": "16",
      "sequence": "A"
    },
    {
      "id": "19",
      "sequence": "TGTGTATTCAACTCACAGAGTTGAACGATCCTTTACTGAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCATGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAATAGGAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACCTTTCTTTTCA"
    },
    {
      "id": "44",
      "sequence": "A"
    },
    {
      "id": "31",
      "sequence": "AGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAAAAGGAAATATCTTCGTATAAAGACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACTTTTCTTTTCATAGAGCAGTTAGGAAATGCTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACAGAGAGCAGACTTGAAACACTGTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAATAGGAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCATAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAAAGACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCCTTG"
    },
    {
      "id": "146",
      "sequence": "CAGAGTTTAAC"
    },
    {
      "id": "74",
      "sequence": "G"
    },
    {
      "id": "61",
      "sequence": "G"
    },
    {
      "id": "29",
      "sequence": "CGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAAAAGGAAATATCTTCGTATAAAGACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACTTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTATGCTAGACAGAATAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACAGAGAGCAGACTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAATAGGAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCACCTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTTGAGGCCTTCTTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTC"
    },
    {
      "id": "159",
      "sequence": "A"
    },
    {
      "id": "101",
      "sequence": "C"
    },
    {
      "id": "105",
      "sequence": "C"
    },
    {
      "id": "17",
      "sequence": "AGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATACTATGATAGACAGAATAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAGCAGACTTGATATACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAATAGGAAATATCTTCCTATAGAAACTAGACAGAATGATTCTCATAAACTCCATTGTGATGTGTGCGTTCAACTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAGACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCCCTTTGAGGCCTTCGATGGAAACGGGATTTCTTCATATTATGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACTTTCATTTAGAGAGAGCAGATTTG"
    },
    {
      "id": "166",
      "sequence": "T"
    },
    {
      "id": "89",
      "sequence": "T"
    },
    {
      "id": "80",
      "sequence": "G"
    },
    {
      "id": "51",
      "sequence": "TTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATC"
    },
    {
      "id": "143",
      "sequence": "C"
    },
    {
      "id": "48",
      "sequence": "C"
    },
    {
      "id": "15",
      "sequence": "AACACTGTTTTTGTGGAATTTGCAAGTGGAGATTTCAAGCGCTTTGGGGCCAAAGGCAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAATCATTGTCAGAAACTGCTGCGTGATGTGTGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCGGTTTGGAAACAGTCTGTTTGTAAATTCTGCACGTGGATATTTTGACCACTTAGAGGCCTTCGTTGGAAACGGGATTTTTTCATGTAAGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTTGTGTGCATTCAACTCACAGAGTTGAACGTTCCCTTAGACAGAGCAGATTTGAAACAGTCTATTTGTGCAATTTGCAAGTGTAGATTTCAAGCGCTTTAAGGTCAATGGCAGAAAAGGAAATATCTTCGTTTCAAAACTAGACAGAATCATTCCCACAAACTGCGTTGTGATGTGTTCGTTCAACTCACAGAGTTTAACCTTTCT"
    },
    {
      "id": "97",
      "sequence": "A"
    },
    {
      "id": "134",
      "sequence": "AAAACTAGACAGAATCATTC"
    },
    {
      "id": "110",
      "sequence": "TT"
    },
    {
      "id": "30",
      "sequence": "AGGCCTTCGTTGGAAACGGGCTTTCTTCATATTCTCCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACGATCCTTTACAGAGAGCAGACTTGAAACACTCTTTCTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATAGTAGAAAAGGAAATATCTTCGTAGAAAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGTGTTCAACTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGTAAACACTCTGTTTATAAAGTCTGCAAGTGGATATTCAGACCCCTTTGAGGCCTTCGTTGGAAACGGGATTTCTTCATATTATGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACTTTCATTTAGAGAGAGCAGATTTGAAACACTGTTTTTGTGGAATTTCCAAGGGAGATTTCAAGCGCTTTGTGGCCAAAGGCAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAATCATTCTCAGAAACTGCTGCGTGATGTGTGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCGGTTTGGAAACACTCTGTTTGTAAGTCTGCACGTGGATATTTTGACCACTTAGAGGCCTTCGTTGGAAACGGGTTTTTTTCATGTAAGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTTGTGTGCATTCAACTCACAGAGTTGAACGTTCCCTTAGACAGAGCAGATTTGAAACACTCTATTTGTGCAATTGGCAAGTGTAGATTTCAAGCGCTTTAAGGTCAATGGCAGAAAAGGAAATATCTTCGTTTCAAAACTTGACAGAATCATTCCCACAAACTGCGTTGTGATGTGTTCGTTCAACTCACAGAGTTTAACCTTTCTGTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGTAAGTGGATATTCTGACATCTTGTGGCCTT"
    },
    {
      "id": "6",
      "sequence": "ACTTGAAACACTCTTT"
    },
    {
      "id": "164",
      "sequence": "CTT"
    },
    {
      "id": "153",
      "sequence": "A"
    },
    {
      "id": "64",
      "sequence": "CA"
    },
    {
      "id": "90",
      "sequence": "TTT"
    },
    {
      "id": "139",
      "sequence": "GTT"
    },
    {
      "id": "4",
      "sequence": "C"
    },
    {
      "id": "13",
      "sequence": "TTCATAGA"
    },
    {
      "id": "104",
      "sequence": "ATTCAACT"
    },
    {
      "id": "52",
      "sequence": "G"
    },
    {
      "id": "43",
      "sequence": "GATATTCAGACCTCTTTGAGG"
    },
    {
      "id": "11",
      "sequence": "CAGTTAGGAAACACTCTGTTTGTAAAGTCTGTAAGTGGATATTCTGACATCTTGTGGCCTTCGTTGGAAACGGGATTTCTTCATATTCTGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTG"
    },
    {
      "id": "69",
      "sequence": "GTT"
    },
    {
      "id": "171",
      "sequence": "T"
    },
    {
      "id": "85",
      "sequence": "GG"
    },
    {
      "id": "119",
      "sequence": "C"
    },
    {
      "id": "39",
      "sequence": "GGAAACGGGATTTCTTCATATTATGCTAGACAGAAGAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTGACAGAGTTGAACTTTCATTTAGAGAGAGCAGATTTGAAACACTGTTTTT"
    },
    {
      "id": "126",
      "sequence": "GG"
    },
    {
      "id": "108",
      "sequence": "TTC"
    },
    {
      "id": "156",
      "sequence": "GGAAACACTCTGTTTGTAAAGTCTG"
    },
    {
      "id": "2",
      "sequence": "T"
    },
    {
      "id": "10",
      "sequence": "T"
    },
    {
      "id": "27",
      "sequence": "GTAA"
    },
    {
      "id": "124",
      "sequence": "AGATTTCAAGCGCTTT"
    },
    {
      "id": "144",
      "sequence": "TTCAACTC"
    },
    {
      "id": "20",
      "sequence": "G"
    },
    {
      "id": "81",
      "sequence": "A"
    },
    {
      "id": "9",
      "sequence": "GTGTATTCAACTCACAGAGTTGAACGATCCTTTACA"
    },
    {
      "id": "109",
      "sequence": "CC"
    },
    {
      "id": "161",
      "sequence": "C"
    },
    {
      "id": "88",
      "sequence": "GG"
    },
    {
      "id": "120",
      "sequence": "AATT"
    },
    {
      "id": "24",
      "sequence": "G"
    },
    {
      "id": "8",
      "sequence": "G"
    },
    {
      "id": "37",
      "sequence": "TGGAATTTGCAAGTGGAGATTTCAAGCGCTTTGGGGCCAAAGGCAGAAAAGGAAATATCTTCGTA"
    },
    {
      "id": "83",
      "sequence": "A"
    },
    {
      "id": "99",
      "sequence": "G"
    },
    {
      "id": "121",
      "sequence": "T"
    },
    {
      "id": "14",
      "sequence": "G"
    },
    {
      "id": "174",
      "sequence": "AT"
    },
    {
      "id": "123",
      "sequence": "T"
    },
    {
      "id": "32",
      "sequence": "C"
    },
    {
      "id": "151",
      "sequence": "AG"
    },
    {
      "id": "54",
      "sequence": "A"
    },
    {
      "id": "63",
      "sequence": "GTAACTTCCTTGTGTT"
    },
    {
      "id": "91",
      "sequence": "T"
    },
    {
      "id": "62",
      "sequence": "GTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAG"
    },
    {
      "id": "150",
      "sequence": "TTCAT"
    },
    {
      "id": "122",
      "sequence": "GCAAGTG"
    },
    {
      "id": "58",
      "sequence": "AGACAGAATAATTCTCA"
    },
    {
      "id": "173",
      "sequence": "C"
    },
    {
      "id": "98",
      "sequence": "A"
    },
    {
      "id": "76",
      "sequence": "C"
    },
    {
      "id": "34",
      "sequence": "G"
    },
    {
      "id": "50",
      "sequence": "C"
    },
    {
      "id": "167",
      "sequence": "GGCCTTCGTTGGAAACGGG"
    },
    {
      "id": "42",
      "sequence": "T"
    },
    {
      "id": "87",
      "sequence": "T"
    },
    {
      "id": "132",
      "sequence": "T"
    },
    {
      "id": "140",
      "sequence": "GTGATGTGT"
    },
    {
      "id": "169",
      "sequence": "TTT"
    },
    {
      "id": "160",
      "sequence": "GTGGATATT"
    },
    {
      "id": "49",
      "sequence": "TTTACACAGAGCAGACTT"
    },
    {
      "id": "106",
      "sequence": "ACAGAGTTGAAC"
    },
    {
      "id": "94",
      "sequence": "CAT"
    },
    {
      "id": "102",
      "sequence": "CAGTAACTTCCTTGTGTTGTGTG"
    },
    {
      "id": "128",
      "sequence": "CAA"
    },
    {
      "id": "70",
      "sequence": "T"
    },
    {
      "id": "21",
      "sequence": "CCTTCTTTGGAAACGGGTTTTTTTCATGTAAGGCTAGACAGAAGAATTCCCAGTAACTTCCTTGTGTT"
    },
    {
      "id": "38",
      "sequence": "G"
    },
    {
      "id": "163",
      "sequence": "AT"
    },
    {
      "id": "131",
      "sequence": "T"
    },
    {
      "id": "53",
      "sequence": "TATTCAGACCTCTTTGAGGCCTTC"
    },
    {
      "id": "47",
      "sequence": "AAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAATAG"
    },
    {
      "id": "175",
      "sequence": "A"
    },
    {
      "id": "3",
      "sequence": "GAGGTCAATGGTAGAATAGG"
    },
    {
      "id": "96",
      "sequence": "T"
    },
    {
      "id": "149",
      "sequence": "G"
    },
    {
      "id": "155",
      "sequence": "A"
    },
    {
      "id": "65",
      "sequence": "T"
    }
  ],
  "edge": [
    {
      "from": "56",
      "from_start": true,
      "to": "57",
      "to_end": true
    },
    {
      "from": "35",
      "from_start": true,
      "to": "36",
      "to_end": true
    },
    {
      "from": "60",
      "from_start": true,
      "to": "61",
      "to_end": true
    },
    {
      "from": "67",
      "from_start": true,
      "to": "68",
      "to_end": true
    },
    {
      "from": "73",
      "to": "74"
    },
    {
      "from": "115",
      "to": "116"
    },
    {
      "from": "112",
      "to": "113"
    },
    {
      "from": "86",
      "to": "87"
    },
    {
      "from": "168",
      "to": "169"
    },
    {
      "from": "12",
      "from_start": true,
      "to": "13",
      "to_end": true
    },
    {
      "from": "75",
      "to": "76"
    },
    {
      "from": "23",
      "from_start": true,
      "to": "24",
      "to_end": true
    },
    {
      "from": "111",
      "to": "112"
    },
    {
      "from": "41",
      "from_start": true,
      "to": "42",
      "to_end": true
    },
    {
      "from": "68",
      "from_start": true,
      "to": "175",
      "to_end": true
    },
    {
      "from": "82",
      "to": "83"
    },
    {
      "from": "130",
      "to": "131"
    },
    {
      "from": "125",
      "to": "126"
    },
    {
      "from": "77",
      "to": "78"
    },
    {
      "from": "172",
      "to": "173"
    },
    {
      "from": "71",
      "to": "72"
    },
    {
      "from": "66",
      "from_start": true,
      "to": "67",
      "to_end": true
    },
    {
      "from": "103",
      "to": "104"
    },
    {
      "from": "59",
      "from_start": true,
      "to": "60",
      "to_end": true
    },
    {
      "from": "26",
      "from_start": true,
      "to": "27",
      "to_end": true
    },
    {
      "from": "127",
      "to": "128"
    },
    {
      "from": "116",
      "to": "117"
    },
    {
      "from": "100",
      "to": "101"
    },
    {
      "from": "79",
      "to": "80"
    },
    {
      "from": "141",
      "to": "142"
    },
    {
      "from": "135",
      "to": "136"
    },
    {
      "from": "138",
      "to": "139"
    },
    {
      "from": "107",
      "to": "108"
    },
    {
      "from": "46",
      "from_start": true,
      "to": "47",
      "to_end": true
    },
    {
      "from": "57",
      "from_start": true,
      "to": "58",
      "to_end": true
    },
    {
      "from": "152",
      "to": "153"
    },
    {
      "from": "170",
      "to": "171"
    },
    {
      "from": "129",
      "to": "130"
    },
    {
      "from": "78",
      "to": "79"
    },
    {
      "from": "133",
      "to": "134"
    },
    {
      "from": "72",
      "to": "73"
    },
    {
      "from": "1",
      "from_start": true,
      "to": "2",
      "to_end": true
    },
    {
      "from": "137",
      "to": "138"
    },
    {
      "from": "22",
      "from_start": true,
      "to": "23",
      "to_end": true
    },
    {
      "from": "154",
      "to": "155"
    },
    {
      "from": "33",
      "from_start": true,
      "to": "34",
      "to_end": true
    },
    {
      "from": "40",
      "from_start": true,
      "to": "41",
      "to_end": true
    },
    {
      "from": "113",
      "to": "114"
    },
    {
      "from": "165",
      "to": "166"
    },
    {
      "from": "142",
      "to": "143"
    },
    {
      "from": "5",
      "from_start": true,
      "to": "6",
      "to_end": true
    },
    {
      "from": "55",
      "from_start": true,
      "to": "56",
      "to_end": true
    },
    {
      "from": "114",
      "to": "115"
    },
    {
      "from": "136",
      "to": "137"
    },
    {
      "from": "117",
      "to": "118"
    },
    {
      "from": "45",
      "from_start": true,
      "to": "46",
      "to_end": true
    },
    {
      "from": "145",
      "to": "146"
    },
    {
      "from": "158",
      "to": "159"
    },
    {
      "from": "28",
      "from_start": true,
      "to": "29",
      "to_end": true
    },
    {
      "from": "148",
      "to": "149"
    },
    {
      "from": "92",
      "to": "93"
    },
    {
      "from": "36",
      "from_start": true,
      "to": "37",
      "to_end": true
    },
    {
      "from": "118",
      "to": "119"
    },
    {
      "from": "162",
      "to": "163"
    },
    {
      "from": "84",
      "to": "85"
    },
    {
      "from": "7",
      "from_start": true,
      "to": "8",
      "to_end": true
    },
    {
      "from": "25",
      "from_start": true,
      "to": "26",
      "to_end": true
    },
    {
      "from": "95",
      "to": "96"
    },
    {
      "from": "93",
      "to": "94"
    },
    {
      "from": "18",
      "from_start": true,
      "to": "19",
      "to_end": true
    },
    {
      "from": "147",
      "to": "148"
    },
    {
      "from": "157",
      "to": "158"
    },
    {
      "from": "16",
      "from_start": true,
      "to": "17",
      "to_end": true
    },
    {
      "from": "19",
      "from_start": true,
      "to": "20",
      "to_end": true
    },
    {
      "from": "44",
      "from_start": true,
      "to": "45",
      "to_end": true
    },
    {
      "from": "31",
      "from_start": true,
      "to": "32",
      "to_end": true
    },
    {
      "from": "146",
      "to": "147"
    },
    {
      "from": "74",
      "to": "75"
    },
    {
      "from": "61",
      "from_start": true,
      "to": "62",
      "to_end": true
    },
    {
      "from": "29",
      "from_start": true,
      "to": "30",
      "to_end": true
    },
    {
      "from": "159",
      "to": "160"
    },
    {
      "from": "101",
      "to": "102"
    },
    {
      "from": "105",
      "to": "106"
    },
    {
      "from": "17",
      "from_start": true,
      "to": "18",
      "to_end": true
    },
    {
      "from": "166",
      "to": "167"
    },
    {
      "from": "89",
      "to": "90"
    },
    {
      "from": "80",
      "to": "81"
    },
    {
      "from": "51",
      "from_start": true,
      "to": "52",
      "to_end": true
    },
    {
      "from": "143",
      "to": "144"
    },
    {
      "from": "48",
      "from_start": true,
      "to": "49",
      "to_end": true
    },
    {
      "from": "15",
      "from_start": true,
      "to": "16",
      "to_end": true
    },
    {
      "from": "97",
      "to": "98"
    },
    {
      "from": "134",
      "to": "135"
    },
    {
      "from": "110",
      "to": "111"
    },
    {
      "from": "30",
      "from_start": true,
      "to": "31",
      "to_end": true
    },
    {
      "from": "6",
      "from_start": true,
      "to": "7",
      "to_end": true
    },
    {
      "from": "164",
      "to": "165"
    },
    {
      "from": "153",
      "to": "154"
    },
    {
      "from": "64",
      "from_start": true,
      "to": "65",
      "to_end": true
    },
    {
      "from": "90",
      "to": "91"
    },
    {
      "from": "139",
      "to": "140"
    },
    {
      "from": "4",
      "from_start": true,
      "to": "5",
      "to_end": true
    },
    {
      "from": "13",
      "from_start": true,
      "to": "14",
      "to_end": true
    },
    {
      "from": "104",
      "to": "105"
    },
    {
      "from": "52",
      "from_start": true,
      "to": "53",
      "to_end": true
    },
    {
      "from": "43",
      "from_start": true,
      "to": "44",
      "to_end": true
    },
    {
      "from": "11",
      "from_start": true,
      "to": "12",
      "to_end": true
    },
    {
      "from": "69",
      "to": "70"
    },
    {
      "from": "171",
      "to": "172"
    },
    {
      "from": "85",
      "to": "86"
    },
    {
      "from": "119",
      "to": "120"
    },
    {
      "from": "39",
      "from_start": true,
      "to": "40",
      "to_end": true
    },
    {
      "from": "126",
      "to": "127"
    },
    {
      "from": "108",
      "to": "109"
    },
    {
      "from": "156",
      "to": "157"
    },
    {
      "from": "2",
      "from_start": true,
      "to": "3",
      "to_end": true
    },
    {
      "from": "10",
      "from_start": true,
      "to": "11",
      "to_end": true
    },
    {
      "from": "27",
      "from_start": true,
      "to": "28",
      "to_end": true
    },
    {
      "from": "124",
      "to": "125"
    },
    {
      "from": "144",
      "to": "145"
    },
    {
      "from": "20",
      "from_start": true,
      "to": "21",
      "to_end": true
    },
    {
      "from": "81",
      "to": "82"
    },
    {
      "from": "9",
      "from_start": true,
      "to": "10",
      "to_end": true
    },
    {
      "from": "109",
      "to": "110"
    },
    {
      "from": "161",
      "to": "162"
    },
    {
      "from": "88",
      "to": "89"
    },
    {
      "from": "120",
      "to": "121"
    },
    {
      "from": "24",
      "from_start": true,
      "to": "25",
      "to_end": true
    },
    {
      "from": "8",
      "from_start": true,
      "to": "9",
      "to_end": true
    },
    {
      "from": "37",
      "from_start": true,
      "to": "38",
      "to_end": true
    },
    {
      "from": "83",
      "to": "84"
    },
    {
      "from": "99",
      "to": "100"
    },
    {
      "from": "121",
      "to": "122"
    },
    {
      "from": "14",
      "from_start": true,
      "to": "15",
      "to_end": true
    },
    {
      "from": "174",
      "to": "175"
    },
    {
      "from": "123",
      "to": "124"
    },
    {
      "from": "32",
      "from_start": true,
      "to": "33",
      "to_end": true
    },
    {
      "from": "151",
      "to": "152"
    },
    {
      "from": "54",
      "from_start": true,
      "to": "55",
      "to_end": true
    },
    {
      "from": "63",
      "from_start": true,
      "to": "64",
      "to_end": true
    },
    {
      "from": "91",
      "to": "92"
    },
    {
      "from": "62",
      "from_start": true,
      "to": "63",
      "to_end": true
    },
    {
      "from": "150",
      "to": "151"
    },
    {
      "from": "122",
      "to": "123"
    },
    {
      "from": "58",
      "from_start": true,
      "to": "59",
      "to_end": true
    },
    {
      "from": "173",
      "to": "174"
    },
    {
      "from": "98",
      "to": "99"
    },
    {
      "from": "76",
      "to": "77"
    },
    {
      "from": "34",
      "from_start": true,
      "to": "35",
      "to_end": true
    },
    {
      "from": "50",
      "from_start": true,
      "to": "51",
      "to_end": true
    },
    {
      "from": "167",
      "to": "168"
    },
    {
      "from": "42",
      "from_start": true,
      "to": "43",
      "to_end": true
    },
    {
      "from": "87",
      "to": "88"
    },
    {
      "from": "132",
      "to": "133"
    },
    {
      "from": "140",
      "to": "141"
    },
    {
      "from": "169",
      "to": "170"
    },
    {
      "from": "160",
      "to": "161"
    },
    {
      "from": "49",
      "from_start": true,
      "to": "50",
      "to_end": true
    },
    {
      "from": "106",
      "to": "107"
    },
    {
      "from": "94",
      "to": "95"
    },
    {
      "from": "102",
      "to": "103"
    },
    {
      "from": "128",
      "to": "129"
    },
    {
      "from": "70",
      "to": "71"
    },
    {
      "from": "21",
      "from_start": true,
      "to": "22",
      "to_end": true
    },
    {
      "from": "38",
      "from_start": true,
      "to": "39",
      "to_end": true
    },
    {
      "from": "163",
      "to": "164"
    },
    {
      "from": "131",
      "to": "132"
    },
    {
      "from": "53",
      "from_start": true,
      "to": "54",
      "to_end": true
    },
    {
      "from": "47",
      "from_start": true,
      "to": "48",
      "to_end": true
    },
    {
      "from": "3",
      "from_start": true,
      "to": "4",
      "to_end": true
    },
    {
      "from": "96",
      "to": "97"
    },
    {
      "from": "149",
      "to": "150"
    },
    {
      "from": "155",
      "to": "156"
    },
    {
      "from": "65",
      "from_start": true,
      "to": "66",
      "to_end": true
    }
  ]
}

            )";

            Graph chunk;
            json2pb(chunk, graph_json.c_str(), graph_json.size());
            vg::VG graph(chunk);

            Alignment aln;
            aln.set_sequence(read_string);
            
            TestAligner aligner_source;
            const Aligner& aligner = *aligner_source.get_regular_aligner();
            aligner.align_pinned(aln, graph, false, true, 150);

            // Check before simplification
            REQUIRE(alignment_is_valid(aln, &graph));

            *(aln.mutable_path()) = simplify(aln.path());
            
            // Check after simplification
            REQUIRE(alignment_is_valid(aln, &graph));
        }
    }
}

