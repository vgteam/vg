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
    }
}

