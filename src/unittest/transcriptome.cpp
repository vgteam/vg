/// \file transcriptome.cpp
///  
/// unit tests for transcriptome class
///

#include <stdio.h>
#include <iostream>

#include "gbwt/dynamic_gbwt.h"
#include "bdsg/packed_graph.hpp"

#include "../transcriptome.hpp"

#include "catch.hpp"

namespace vg {
    namespace unittest {

        vector<vector<uint64_t> > transcript_paths_to_int_vectors(const vector<CompletedTranscriptPath> & transcript_paths) {

            vector<vector<uint64_t> > int_vectors;
            int_vectors.reserve(transcript_paths.size());

            for (auto & transcript_path: transcript_paths) {

                int_vectors.emplace_back(vector<uint64_t>());
                int_vectors.back().reserve(transcript_path.path.size());

                for (auto & handle: transcript_path.path) {

                    int_vectors.back().emplace_back(bdsg::as_integer(handle));
                }
            }

            sort(int_vectors.begin(), int_vectors.end());
            return int_vectors;
        }

        vector<string> transcript_paths_to_sequences(const vector<CompletedTranscriptPath> & transcript_paths, const MutablePathDeletableHandleGraph & graph) {

            vector<string> sequences;
            sequences.reserve(transcript_paths.size());

            for (auto & transcript_path: transcript_paths) {

                sequences.emplace_back("");

                for (auto & handle: transcript_path.path) {

                    sequences.back() += graph.get_sequence(handle);
                }
            }

            sort(sequences.begin(), sequences.end());
            return sequences;
        }

        TEST_CASE("Transcriptome can add splice-junctions and project transcripts", "[transcriptome]") {
         
            unique_ptr<MutablePathDeletableHandleGraph> graph(new bdsg::PackedGraph);
               
            handle_t node1 = graph->create_handle("AAAA");
            handle_t node2 = graph->create_handle("CC");
            handle_t node3 = graph->create_handle("G");
            handle_t node4 = graph->create_handle("TTTTTTTT");
            handle_t node5 = graph->create_handle("CCC");
            handle_t node6 = graph->create_handle("AAAA");
            
            graph->create_edge(node1, node2);
            graph->create_edge(node1, node3);
            graph->create_edge(node2, node4);
            graph->create_edge(node3, node4);
            graph->create_edge(node4, node5);
            graph->create_edge(node4, node6);
            graph->create_edge(node5, node6);

            path_handle_t path1 = graph->create_path_handle("path1");
            graph->append_step(path1, node1);
            graph->append_step(path1, node2);
            graph->append_step(path1, node4);
            graph->append_step(path1, node5);
            graph->append_step(path1, node6);

            path_handle_t path2 = graph->create_path_handle("path2");
            graph->append_step(path2, node1);
            graph->append_step(path2, node2);
            graph->append_step(path2, node4);
            graph->append_step(path2, node6);

            path_handle_t path3 = graph->create_path_handle("path3");
            graph->append_step(path3, node1);
            graph->append_step(path3, node3);
            graph->append_step(path3, node4);
            graph->append_step(path3, node5);
            graph->append_step(path3, node6);

            path_handle_t path4 = graph->create_path_handle("path4");
            graph->append_step(path4, node1);
            graph->append_step(path4, node3);
            graph->append_step(path4, node4);
            graph->append_step(path4, node6);

            Transcriptome transcriptome(move(graph));
            REQUIRE(graph == nullptr);

            transcriptome.num_threads = 2;

            REQUIRE(transcriptome.reference_transcript_paths().size() == 0);
            REQUIRE(transcriptome.haplotype_transcript_paths().size() == 0);

            REQUIRE(transcriptome.graph().get_node_count() == 6);
            REQUIRE(transcriptome.graph().get_edge_count() == 7);
            REQUIRE(transcriptome.graph().get_path_count() == 4);

            REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path1")) == 5);
            REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path2")) == 4);
            REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path3")) == 5);
            REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path4")) == 4);

            unique_ptr<gbwt::GBWT> empty_haplotype_index(new gbwt::GBWT());

            stringstream transcript_stream;
            transcript_stream << "# header 1" << endl;
            transcript_stream << "# header\t2\t" << endl;
            transcript_stream << "path1\t.\texon\t2\t7\t.\t+\t.\ttranscript_id \"transcript1\"; exon_number 1;" << endl;
            transcript_stream << "path1\t.\texon\t9\t10\t.\t+\t.\ttranscript_id \"transcript1\"; exon_number 2;" << endl;
            transcript_stream << "path1\t.\texon\t19\t21\t.\t+\t.\ttranscript_id \"transcript1\"; exon_number 3;" << endl;
            transcript_stream << "path1\t.\texon\t2\t7\t.\t-\t.\texon_number 1; transcript_id \"transcript2\";" << endl;
            transcript_stream << "path1\t.\texon\t16\t21\t.\t-\t.\texon_number 2; transcript_id \"transcript2\";" << endl;
            transcript_stream << "path1\t.\texon\t9\t11\t.\t+\t.\ttranscript_id \"transcript3\";" << endl;
            transcript_stream << "path1\t.\texon\t18\t21\t.\t+\t.\ttranscript_id \"transcript3\";" << endl;

            SECTION("Transcriptome can add splice-junctions and reference transcript paths") {

                transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream}), empty_haplotype_index, false, false);

                REQUIRE(transcriptome.transcript_paths().size() == 3);
                REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                REQUIRE(transcriptome.haplotype_transcript_paths().size() == 0);

                REQUIRE(transcriptome.sort_compact_nodes());

                REQUIRE(transcriptome.graph().get_node_count() == 13);
                REQUIRE(transcriptome.graph().get_edge_count() == 18);
                REQUIRE(transcriptome.graph().get_path_count() == 4);

                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path1")) == 12);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path2")) == 10);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path3")) == 12);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path4")) == 10);

                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(1)) == "A");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(2)) == "AAA");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(3)) == "CC");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(4)) == "G");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(5)) == "T");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(6)) == "T");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(7)) == "TT");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(8)) == "T");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(9)) == "TTT");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(10)) == "C");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(11)) == "CC");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(12)) == "A");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(13)) == "AAA");

                auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                REQUIRE(int_ref_transcript_paths.at(1) == vector<uint64_t>({14, 16, 24, 26}));
                REQUIRE(int_ref_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 7, 5}));

                auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                REQUIRE(seq_ref_transcript_paths.front() == "AAACCTTTAAA");
                REQUIRE(seq_ref_transcript_paths.at(1) == "TTTAAAA");
                REQUIRE(seq_ref_transcript_paths.back() == "TTTTGGAGGTTT");

                transcript_stream.clear();
                transcript_stream.seekg(0,ios::beg);

                SECTION("Transcriptome can project transcripts onto embedded haplotype paths") {

                    transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *empty_haplotype_index, true);
                    
                    REQUIRE(transcriptome.transcript_paths().size() == 5);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 4);

                    auto int_ref_transcript_paths2 = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());
                    REQUIRE(int_ref_transcript_paths == int_ref_transcript_paths2);
                    
                    auto int_hap_transcript_paths = transcript_paths_to_int_vectors(transcriptome.haplotype_transcript_paths());

                    REQUIRE(int_hap_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(1) == vector<uint64_t>({4, 8, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(2) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_hap_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 9, 5}));

                    auto seq_hap_transcript_paths = transcript_paths_to_sequences(transcriptome.haplotype_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_hap_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(1) == "AAAGTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(2) == "TTTAAAA");
                    REQUIRE(seq_hap_transcript_paths.back() == "TTTTGGACTTT");
                }

                SECTION("Transcriptome can project transcripts onto embedded haplotype paths and not collapse redundant paths") {

                    transcriptome.path_collapse_type = "no";
                    
                    transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *empty_haplotype_index, true);

                    REQUIRE(transcriptome.transcript_paths().size() == 8);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 5);

                    auto int_ref_transcript_paths2 = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());
                    REQUIRE(int_ref_transcript_paths == int_ref_transcript_paths2);

                    auto int_hap_transcript_paths = transcript_paths_to_int_vectors(transcriptome.haplotype_transcript_paths());

                    REQUIRE(int_hap_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(1) == vector<uint64_t>({4, 8, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(2) == vector<uint64_t>({4, 8, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(3) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_hap_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 9, 5}));

                    auto seq_hap_transcript_paths = transcript_paths_to_sequences(transcriptome.haplotype_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_hap_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(1) == "AAAGTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(2) == "AAAGTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(3) == "TTTAAAA");
                    REQUIRE(seq_hap_transcript_paths.back() == "TTTTGGACTTT");
                }

                SECTION("Transcriptome can remove non-transcribed nodes") {

                    transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *empty_haplotype_index, true);

                    REQUIRE(transcriptome.transcript_paths().size() == 5);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 4);

                    transcriptome.remove_non_transcribed_nodes();

                    REQUIRE(transcriptome.graph().get_node_count() == 9);
                    REQUIRE(transcriptome.graph().get_edge_count() == 11);
                    REQUIRE(transcriptome.graph().get_path_count() == 0);   

                    auto int_ref_transcript_paths2 = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());
                    REQUIRE(int_ref_transcript_paths == int_ref_transcript_paths2);

                    auto int_hap_transcript_paths = transcript_paths_to_int_vectors(transcriptome.haplotype_transcript_paths());

                    REQUIRE(int_hap_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(1) == vector<uint64_t>({4, 8, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(2) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_hap_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 9, 5}));             

                    auto seq_hap_transcript_paths = transcript_paths_to_sequences(transcriptome.haplotype_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_hap_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(1) == "AAAGTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(2) == "TTTAAAA");
                    REQUIRE(seq_hap_transcript_paths.back() == "TTTTGGACTTT");
                }

                SECTION("Transcriptome can chop long nodes") {

                    transcriptome.chop_nodes(2);

                    REQUIRE(transcriptome.sort_compact_nodes());

                    REQUIRE(transcriptome.graph().get_node_count() == 16);
                    REQUIRE(transcriptome.graph().get_edge_count() == 21);
                    REQUIRE(transcriptome.graph().get_path_count() == 4);

                    REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path1")) == 15);
                    REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path2")) == 13);
                    REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path3")) == 15);
                    REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path4")) == 13);

                    auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                    REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({4, 6, 8, 12, 16, 30, 32}));
                    REQUIRE(int_ref_transcript_paths.at(1) == vector<uint64_t>({16, 18, 28, 30, 32}));
                    REQUIRE(int_ref_transcript_paths.back() == vector<uint64_t>({33, 31, 29, 27, 13, 9, 7, 5}));

                    auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_ref_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_ref_transcript_paths.at(1) == "TTTAAAA");
                    REQUIRE(seq_ref_transcript_paths.back() == "TTTTGGAGGTTT");
                }

                SECTION("Transcriptome can add reference transcript paths to graph") {

                    transcriptome.embed_transcript_paths(true, false);
                    REQUIRE(transcriptome.graph().get_path_count() == 7);            
      
                    transcriptome.embed_transcript_paths(false, true);
                    REQUIRE(transcriptome.graph().get_path_count() == 7);
                }

                SECTION("Transcriptome can add haplotype transcript paths to graph") {

                    transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *empty_haplotype_index, true);

                    transcriptome.embed_transcript_paths(false, true);
                    REQUIRE(transcriptome.graph().get_path_count() == 8);  
                }

                gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
                gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(11, true)));

                gbwt::vector_type gbwt_thread_1(10);
                gbwt::vector_type gbwt_thread_2(12);
   
                gbwt_thread_1[0] = gbwt::Node::encode(1, false);
                gbwt_thread_1[1] = gbwt::Node::encode(2, false);
                gbwt_thread_1[2] = gbwt::Node::encode(4, false);
                gbwt_thread_1[3] = gbwt::Node::encode(5, false);
                gbwt_thread_1[4] = gbwt::Node::encode(6, false);
                gbwt_thread_1[5] = gbwt::Node::encode(7, false);
                gbwt_thread_1[6] = gbwt::Node::encode(8, false);
                gbwt_thread_1[7] = gbwt::Node::encode(9, false);
                gbwt_thread_1[8] = gbwt::Node::encode(12, false);
                gbwt_thread_1[9] = gbwt::Node::encode(13, false);

                gbwt_thread_2[0] = gbwt::Node::encode(1, false);
                gbwt_thread_2[1] = gbwt::Node::encode(2, false);
                gbwt_thread_2[2] = gbwt::Node::encode(4, false);
                gbwt_thread_2[3] = gbwt::Node::encode(5, false);
                gbwt_thread_2[4] = gbwt::Node::encode(6, false);
                gbwt_thread_2[5] = gbwt::Node::encode(7, false);
                gbwt_thread_2[6] = gbwt::Node::encode(8, false);
                gbwt_thread_2[7] = gbwt::Node::encode(9, false);
                gbwt_thread_2[8] = gbwt::Node::encode(10, false);
                gbwt_thread_2[9] = gbwt::Node::encode(11, false);
                gbwt_thread_2[10] = gbwt::Node::encode(12, false);
                gbwt_thread_2[11] = gbwt::Node::encode(13, false);

                gbwt_builder.insert(gbwt_thread_1, true);
                gbwt_builder.insert(gbwt_thread_2, true);

                gbwt_builder.finish();

                std::stringstream gbwt_stream;
                gbwt_builder.index.serialize(gbwt_stream);

                unique_ptr<gbwt::GBWT> haplotype_index(new gbwt::GBWT());
                
                haplotype_index->load(gbwt_stream);
                REQUIRE(haplotype_index->bidirectional());   

                SECTION("Transcriptome can project transcripts onto GBWT haplotypes") {

                    transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *haplotype_index, false);

                    REQUIRE(transcriptome.transcript_paths().size() == 5);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 3);

                    auto int_hap_transcript_paths = transcript_paths_to_int_vectors(transcriptome.haplotype_transcript_paths());

                    REQUIRE(int_hap_transcript_paths.front() == vector<uint64_t>({4, 8, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(1) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_hap_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 9, 5}));

                    auto seq_hap_transcript_paths = transcript_paths_to_sequences(transcriptome.haplotype_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_hap_transcript_paths.front() == "AAAGTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(1) == "TTTAAAA");
                    REQUIRE(seq_hap_transcript_paths.back() == "TTTTGGACTTT");
                }

                SECTION("Transcriptome can project transcripts onto embedded paths and GBWT haplotypes") {

                    transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *haplotype_index, true);

                    REQUIRE(transcriptome.transcript_paths().size() == 5);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 4);

                    auto int_hap_transcript_paths = transcript_paths_to_int_vectors(transcriptome.haplotype_transcript_paths());

                    REQUIRE(int_hap_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(1) == vector<uint64_t>({4, 8, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(2) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_hap_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 9, 5}));

                    auto seq_hap_transcript_paths = transcript_paths_to_sequences(transcriptome.haplotype_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_hap_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(1) == "AAAGTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(2) == "TTTAAAA");
                    REQUIRE(seq_hap_transcript_paths.back() == "TTTTGGACTTT");
                }

                SECTION("Transcriptome can project transcripts onto embedded paths and GBWT haplotypes and collapse redundant paths across haplotypes") {

                    transcript_stream << "path1\t.\texon\t2\t7\t.\t-\t.\texon_number 1; transcript_id \"transcript4\";" << endl;
                    transcript_stream << "path1\t.\texon\t16\t21\t.\t-\t.\texon_number 2; transcript_id \"transcript4\";" << endl;

                    transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *haplotype_index, true);

                    REQUIRE(transcriptome.transcript_paths().size() == 6);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 5);

                    auto int_hap_transcript_paths = transcript_paths_to_int_vectors(transcriptome.haplotype_transcript_paths());

                    REQUIRE(int_hap_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(1) == vector<uint64_t>({4, 8, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(2) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_hap_transcript_paths.at(3) == vector<uint64_t>({27, 25, 23, 11, 9, 5}));
                    REQUIRE(int_hap_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 9, 5}));

                    auto seq_hap_transcript_paths = transcript_paths_to_sequences(transcriptome.haplotype_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_hap_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(1) == "AAAGTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(2) == "TTTAAAA");
                    REQUIRE(seq_hap_transcript_paths.at(3) == "TTTTGGACTTT");
                    REQUIRE(seq_hap_transcript_paths.back() == "TTTTGGACTTT");
                }

                SECTION("Transcriptome can project transcripts onto embedded paths and GBWT haplotypes and collapse redundant paths across haplotypes and transcripts") {

                    transcript_stream << "path1\t.\texon\t2\t7\t.\t-\t.\texon_number 1; transcript_id \"transcript4\";" << endl;
                    transcript_stream << "path1\t.\texon\t16\t21\t.\t-\t.\texon_number 2; transcript_id \"transcript4\";" << endl;

                    transcriptome.path_collapse_type = "all";

                    transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *haplotype_index, true);

                    REQUIRE(transcriptome.transcript_paths().size() == 5);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 4);

                    auto int_hap_transcript_paths = transcript_paths_to_int_vectors(transcriptome.haplotype_transcript_paths());

                    REQUIRE(int_hap_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(1) == vector<uint64_t>({4, 8, 10, 14, 26}));
                    REQUIRE(int_hap_transcript_paths.at(2) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_hap_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 9, 5}));

                    auto seq_hap_transcript_paths = transcript_paths_to_sequences(transcriptome.haplotype_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_hap_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(1) == "AAAGTTTAAA");
                    REQUIRE(seq_hap_transcript_paths.at(2) == "TTTAAAA");
                    REQUIRE(seq_hap_transcript_paths.back() == "TTTTGGACTTT");
                }
            }

            SECTION("Transcriptome can parse gff3 file and reorder reverse exons") {

                stringstream transcript_stream2;
                transcript_stream2 << "path1\t.\texon\t2\t7\t.\t+\t.\ttranscript_id=transcript1;exon_number=0" << endl;
                transcript_stream2 << "path1\t.\texon\t9\t10\t.\t+\t.\ttranscript_id=transcript1;exon_number=1" << endl;
                transcript_stream2 << "path1\t.\texon\t19\t21\t.\t+\t.\ttranscript_id=transcript1;exon_number=2" << endl;
                transcript_stream2 << "path1\t.\texon\t16\t21\t.\t-\t.\tID=exon:transcript2:0;transcript_id=transcript2;" << endl;
                transcript_stream2 << "path1\t.\texon\t2\t7\t.\t-\t.\tID=exon:transcript2:1;transcript_id=transcript2" << endl;
                transcript_stream2 << "path1\t.\texon\t9\t11\t.\t+\t.\ttranscript_id=transcript3" << endl;
                transcript_stream2 << "path1\t.\texon\t18\t21\t.\t+\t.\ttranscript_id=transcript3" << endl;

                transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream2}), empty_haplotype_index, false, false);

                REQUIRE(transcriptome.transcript_paths().size() == 3);
                REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                REQUIRE(transcriptome.haplotype_transcript_paths().size() == 0);

                REQUIRE(transcriptome.sort_compact_nodes());

                REQUIRE(transcriptome.graph().get_node_count() == 13);
                REQUIRE(transcriptome.graph().get_edge_count() == 18);
                REQUIRE(transcriptome.graph().get_path_count() == 4);

                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path1")) == 12);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path2")) == 10);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path3")) == 12);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path4")) == 10);

                auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                REQUIRE(int_ref_transcript_paths.at(1) == vector<uint64_t>({14, 16, 24, 26}));
                REQUIRE(int_ref_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 7, 5}));

                auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                REQUIRE(seq_ref_transcript_paths.front() == "AAACCTTTAAA");
                REQUIRE(seq_ref_transcript_paths.at(1) == "TTTAAAA");
                REQUIRE(seq_ref_transcript_paths.back() == "TTTTGGAGGTTT");
            }
            
            SECTION("Transcriptome can parse gff3 file and exclude transcripts with incorrect exon order") {

                stringstream transcript_stream2;
                transcript_stream2 << "path1\t.\texon\t9\t10\t.\t+\t.\ttranscript_id=transcript1;exon_number=2" << endl;
                transcript_stream2 << "path1\t.\texon\t2\t7\t.\t+\t.\ttranscript_id=transcript1;exon_number=1" << endl;
                transcript_stream2 << "path1\t.\texon\t19\t21\t.\t+\t.\ttranscript_id=transcript1;exon_number=3" << endl;
                transcript_stream2 << "path1\t.\texon\t16\t21\t.\t-\t.\tID=exon:transcript2:0;transcript_id=transcript2;" << endl;
                transcript_stream2 << "path1\t.\texon\t2\t7\t.\t-\t.\tID=exon:transcript2:1;transcript_id=transcript2" << endl;
                transcript_stream2 << "path1\t.\texon\t8\t9\t.\t-\t.\tID=exon:transcript2:2;transcript_id=transcript2" << endl;
                transcript_stream2 << "path1\t.\texon\t9\t11\t.\t+\t.\texon_number=1;transcript_id=transcript3" << endl;
                transcript_stream2 << "path1\t.\texon\t18\t21\t.\t+\t.\texon_number=2;transcript_id=transcript3" << endl;
                transcript_stream2 << "path1\t.\texon\t8\t10\t.\t+\t.\ttranscript_id=transcript4" << endl;
                transcript_stream2 << "path1\t.\texon\t3\t6\t.\t+\t.\ttranscript_id=transcript4" << endl;

                transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream2}), empty_haplotype_index, false, false);

                REQUIRE(transcriptome.transcript_paths().size() == 1);
                REQUIRE(transcriptome.reference_transcript_paths().size() == 1);

                REQUIRE(transcriptome.sort_compact_nodes());

                REQUIRE(transcriptome.graph().get_node_count() == 8);
                REQUIRE(transcriptome.graph().get_edge_count() == 10);
                REQUIRE(transcriptome.graph().get_path_count() == 4);

                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path1")) == 7);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path2")) == 6);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path3")) == 7);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path4")) == 6);

                auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({10, 16}));

                auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                REQUIRE(seq_ref_transcript_paths.front() == "TTTAAAA");
            }

            SECTION("Transcriptome can parse gff3 file and exclude transcripts with overlapping exons") {

                stringstream transcript_stream2;
                transcript_stream2 << "path1\t.\texon\t2\t7\t.\t+\t.\ttranscript_id=transcript1;exon_number=0" << endl;
                transcript_stream2 << "path1\t.\texon\t7\t10\t.\t+\t.\ttranscript_id=transcript1;exon_number=1" << endl;
                transcript_stream2 << "path1\t.\texon\t19\t21\t.\t+\t.\ttranscript_id=transcript1;exon_number=2" << endl;
                transcript_stream2 << "path1\t.\texon\t16\t21\t.\t-\t.\tID=exon:transcript2:0;transcript_id=transcript2;" << endl;
                transcript_stream2 << "path1\t.\texon\t2\t17\t.\t-\t.\tID=exon:transcript2:1;transcript_id=transcript2" << endl;
                transcript_stream2 << "path1\t.\texon\t9\t11\t.\t+\t.\ttranscript_id=transcript3" << endl;
                transcript_stream2 << "path1\t.\texon\t18\t21\t.\t+\t.\ttranscript_id=transcript3" << endl;

                transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream2}), empty_haplotype_index, false, false);

                REQUIRE(transcriptome.transcript_paths().size() == 1);
                REQUIRE(transcriptome.reference_transcript_paths().size() == 1);

                REQUIRE(transcriptome.sort_compact_nodes());

                REQUIRE(transcriptome.graph().get_node_count() == 8);
                REQUIRE(transcriptome.graph().get_edge_count() == 10);
                REQUIRE(transcriptome.graph().get_path_count() == 4);

                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path1")) == 7);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path2")) == 6);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path3")) == 7);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path4")) == 6);

                auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({10, 16}));

                auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                REQUIRE(seq_ref_transcript_paths.front() == "TTTAAAA");
            }

            SECTION("Transcriptome can add splice-junctions and update GBWT threads") {

                gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
                gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(6, true)));

                gbwt::vector_type gbwt_thread_1(4);
                gbwt::vector_type gbwt_thread_2(5);
   
                gbwt_thread_1[0] = gbwt::Node::encode(1, false);
                gbwt_thread_1[1] = gbwt::Node::encode(3, false);
                gbwt_thread_1[2] = gbwt::Node::encode(4, false);
                gbwt_thread_1[3] = gbwt::Node::encode(6, false);

                gbwt_thread_2[0] = gbwt::Node::encode(1, false);
                gbwt_thread_2[1] = gbwt::Node::encode(3, false);
                gbwt_thread_2[2] = gbwt::Node::encode(4, false);
                gbwt_thread_2[3] = gbwt::Node::encode(5, false);
                gbwt_thread_2[4] = gbwt::Node::encode(6, false);

                gbwt_builder.insert(gbwt_thread_1, true);
                gbwt_builder.insert(gbwt_thread_2, true);

                gbwt_builder.finish();

                std::stringstream gbwt_stream;
                gbwt_builder.index.serialize(gbwt_stream);

                unique_ptr<gbwt::GBWT> haplotype_index(new gbwt::GBWT());
                
                haplotype_index->load(gbwt_stream);
                REQUIRE(haplotype_index->bidirectional());

                transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream}), haplotype_index, false, true);

                transcript_stream.clear();
                transcript_stream.seekg(0,ios::beg);

                transcriptome.add_haplotype_transcripts(vector<istream *>({&transcript_stream}), *haplotype_index, false);

                REQUIRE(transcriptome.transcript_paths().size() == 5);
                REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                REQUIRE(transcriptome.haplotype_transcript_paths().size() == 3);

                REQUIRE(transcriptome.sort_compact_nodes());

                auto int_hap_transcript_paths = transcript_paths_to_int_vectors(transcriptome.haplotype_transcript_paths());

                REQUIRE(int_hap_transcript_paths.front() == vector<uint64_t>({4, 8, 10, 14, 26}));
                REQUIRE(int_hap_transcript_paths.at(1) == vector<uint64_t>({14, 16, 24, 26}));
                REQUIRE(int_hap_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 9, 5}));

                auto seq_hap_transcript_paths = transcript_paths_to_sequences(transcriptome.haplotype_transcript_paths(), transcriptome.graph());

                REQUIRE(seq_hap_transcript_paths.front() == "AAAGTTTAAA");
                REQUIRE(seq_hap_transcript_paths.at(1) == "TTTAAAA");
                REQUIRE(seq_hap_transcript_paths.back() == "TTTTGGACTTT");
            }

            SECTION("Transcriptome can use GBWT haplotype transcript annotations") {

                gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
                gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(6, true)));

                gbwt_builder.index.addMetadata();

                gbwt::vector_type gbwt_thread_1(5);
                gbwt::vector_type gbwt_thread_2(4);
                gbwt::vector_type gbwt_thread_3(3);
                gbwt::vector_type gbwt_thread_4(1);

                gbwt_thread_1[0] = gbwt::Node::encode(1, false);
                gbwt_thread_1[1] = gbwt::Node::encode(2, false);
                gbwt_thread_1[2] = gbwt::Node::encode(4, false);
                gbwt_thread_1[3] = gbwt::Node::encode(5, false);
                gbwt_thread_1[4] = gbwt::Node::encode(6, false);

                gbwt_thread_2[0] = gbwt::Node::encode(1, false);
                gbwt_thread_2[1] = gbwt::Node::encode(2, false);
                gbwt_thread_2[2] = gbwt::Node::encode(4, false);
                gbwt_thread_2[3] = gbwt::Node::encode(6, false);

                gbwt_thread_3[0] = gbwt::Node::encode(1, false);
                gbwt_thread_3[1] = gbwt::Node::encode(2, false);
                gbwt_thread_3[2] = gbwt::Node::encode(4, false);

                gbwt_thread_4[0] = gbwt::Node::encode(6, false);

                gbwt_builder.insert(gbwt_thread_1, true);
                gbwt_builder.index.metadata.addPath(0, 0, 0, 0);

                gbwt_builder.insert(gbwt_thread_2, true);
                gbwt_builder.index.metadata.addPath(0, 1, 0, 0);

                gbwt_builder.insert(gbwt_thread_3, true);
                gbwt_builder.index.metadata.addPath(1, 2, 0, 0);

                gbwt_builder.insert(gbwt_thread_4, true);
                gbwt_builder.index.metadata.addPath(1, 2, 0, 17);

                gbwt_builder.index.metadata.addSamples(vector<string>({"sample1", "sample2"}));
                gbwt_builder.index.metadata.addContigs(vector<string>({"path1", "path2", "path3"}));

                gbwt_builder.finish();

                std::stringstream gbwt_stream;
                gbwt_builder.index.serialize(gbwt_stream);

                unique_ptr<gbwt::GBWT> haplotype_index(new gbwt::GBWT());
                
                haplotype_index->load(gbwt_stream);
                REQUIRE(haplotype_index->bidirectional());

                stringstream transcript_stream2;
                transcript_stream2 << "sample1#0#path1\t.\texon\t2\t7\t.\t+\t.\ttranscript_id \"transcript1\"; exon_number 1;" << endl;
                transcript_stream2 << "sample1#0#path1\t.\texon\t9\t10\t.\t+\t.\ttranscript_id \"transcript1\"; exon_number 2;" << endl;
                transcript_stream2 << "sample1#0#path1\t.\texon\t19\t21\t.\t+\t.\ttranscript_id \"transcript1\"; exon_number 3;" << endl;
                transcript_stream2 << "sample1#0#path1\t.\texon\t2\t7\t.\t-\t.\texon_number 1; transcript_id \"transcript2\";" << endl;
                transcript_stream2 << "sample1#0#path1\t.\texon\t16\t21\t.\t-\t.\texon_number 2; transcript_id \"transcript2\";" << endl;
                transcript_stream2 << "sample1#0#path1\t.\texon\t9\t11\t.\t+\t.\ttranscript_id \"transcript3\";" << endl;
                transcript_stream2 << "sample1#0#path1\t.\texon\t18\t21\t.\t+\t.\ttranscript_id \"transcript3\";" << endl;

                SECTION("Transcriptome can add splice-junctions and reference transcript paths using GBWT haplotypes") {

                    transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream2}), haplotype_index, true, false);

                    REQUIRE(transcriptome.transcript_paths().size() == 3);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 0);
                    
                    REQUIRE(transcriptome.sort_compact_nodes());

                    REQUIRE(transcriptome.graph().get_node_count() == 13);
                    REQUIRE(transcriptome.graph().get_edge_count() == 18);
                    REQUIRE(transcriptome.graph().get_path_count() == 4);

                    REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path1")) == 12);
                    REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path2")) == 10);
                    REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path3")) == 12);
                    REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path4")) == 10);

                    auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                    REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_ref_transcript_paths.at(1) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_ref_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 7, 5}));

                    auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_ref_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_ref_transcript_paths.at(1) == "TTTAAAA");
                    REQUIRE(seq_ref_transcript_paths.back() == "TTTTGGAGGTTT");
                }
         
                SECTION("Transcriptome can collapse redundant GBWT annotated reference transcript paths across haplotypes") {

                    stringstream transcript_stream3;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t2\t7\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t9\t10\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t16\t18\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t9\t11\t.\t+\t.\ttranscript_id \"transcript4\";" << endl;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t15\t18\t.\t+\t.\ttranscript_id \"transcript4\";" << endl;

                    transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream2, &transcript_stream3}), haplotype_index, true, false);

                    REQUIRE(transcriptome.transcript_paths().size() == 4);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 4);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 0);

                    REQUIRE(transcriptome.sort_compact_nodes());

                    auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                    REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_ref_transcript_paths.at(1) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_ref_transcript_paths.at(2) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_ref_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 7, 5}));

                    auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_ref_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_ref_transcript_paths.at(1) == "TTTAAAA");
                    REQUIRE(seq_ref_transcript_paths.at(2) == "TTTAAAA");
                    REQUIRE(seq_ref_transcript_paths.back() == "TTTTGGAGGTTT");
                }

                SECTION("Transcriptome can collapse redundant GBWT annotated reference transcript paths across haplotypes and transcripts") {

                    stringstream transcript_stream3;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t2\t7\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t9\t10\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t16\t18\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t9\t11\t.\t+\t.\ttranscript_id \"transcript4\";" << endl;
                    transcript_stream3 << "sample1#0#path2\t.\texon\t15\t18\t.\t+\t.\ttranscript_id \"transcript4\";" << endl;

                    transcriptome.path_collapse_type = "all";

                    transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream2, &transcript_stream3}), haplotype_index, true, false);

                    REQUIRE(transcriptome.transcript_paths().size() == 3);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 3);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 0);

                    REQUIRE(transcriptome.sort_compact_nodes());

                    auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                    REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_ref_transcript_paths.at(1) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_ref_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 7, 5}));

                    auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_ref_transcript_paths.front() == "AAACCTTTAAA");
                    REQUIRE(seq_ref_transcript_paths.at(1) == "TTTAAAA");
                    REQUIRE(seq_ref_transcript_paths.back() == "TTTTGGAGGTTT");
                }

                SECTION("Transcriptome can construct GBWT annotated reference transcript paths spanning multiple broken contigs") {

                    transcript_stream2 << "sample2#0#path3\t.\texon\t9\t10\t.\t-\t.\ttranscript_id \"transcript4\";" << endl;
                    transcript_stream2 << "sample2#0#path3\t.\texon\t19\t21\t.\t-\t.\ttranscript_id \"transcript4\";" << endl;
                    transcript_stream2 << "sample2#0#path3\t.\texon\t9\t11\t.\t+\t.\ttranscript_id \"transcript5\";" << endl;
                    transcript_stream2 << "sample2#0#path3\t.\texon\t16\t21\t.\t+\t.\ttranscript_id \"transcript5\";" << endl;
                    transcript_stream2 << "sample2#0#path3\t.\texon\t19\t21\t.\t+\t.\ttranscript_id \"transcript6\";" << endl;

                    transcriptome.add_reference_transcripts(vector<istream *>({&transcript_stream2}), haplotype_index, true, false);

                    REQUIRE(transcriptome.transcript_paths().size() == 5);
                    REQUIRE(transcriptome.reference_transcript_paths().size() == 5);
                    REQUIRE(transcriptome.haplotype_transcript_paths().size() == 0);

                    REQUIRE(transcriptome.sort_compact_nodes());

                    auto int_ref_transcript_paths = transcript_paths_to_int_vectors(transcriptome.reference_transcript_paths());

                    REQUIRE(int_ref_transcript_paths.front() == vector<uint64_t>({4, 6, 10, 14, 26}));
                    REQUIRE(int_ref_transcript_paths.at(1) == vector<uint64_t>({14, 16, 24, 26}));
                    REQUIRE(int_ref_transcript_paths.at(2) == vector<uint64_t>({26}));
                    REQUIRE(int_ref_transcript_paths.at(3) == vector<uint64_t>({27, 15}));
                    REQUIRE(int_ref_transcript_paths.back() == vector<uint64_t>({27, 25, 23, 11, 7, 5}));

                    auto seq_ref_transcript_paths = transcript_paths_to_sequences(transcriptome.reference_transcript_paths(), transcriptome.graph());

                    REQUIRE(seq_ref_transcript_paths.front() == "AAA");
                    REQUIRE(seq_ref_transcript_paths.at(1) == "AAACCTTTAAA");
                    REQUIRE(seq_ref_transcript_paths.at(2) == "TTTAA");
                    REQUIRE(seq_ref_transcript_paths.at(3) == "TTTAAAA");
                    REQUIRE(seq_ref_transcript_paths.back() == "TTTTGGAGGTTT");
                }
            }

            SECTION("Transcriptome can parse intron BED file and add splice-junctions") {

                stringstream intron_stream;
                intron_stream << "# header 1" << endl;
                intron_stream << "# header\t2\t" << endl;
                intron_stream << "path1\t7\t8" << endl;
                intron_stream << "path1\t10\t18\texon\t0" << endl;
                intron_stream << "path1\t7\t15\texon\t0\t-" << endl;
                intron_stream << "path1\t11\t17\texon\t0\t+\t.\t." << endl;

                transcriptome.add_intron_splice_junctions(vector<istream *>({&intron_stream}), empty_haplotype_index, false);

                REQUIRE(transcriptome.transcript_paths().empty());

                REQUIRE(transcriptome.sort_compact_nodes());

                REQUIRE(transcriptome.graph().get_node_count() == 12);
                REQUIRE(transcriptome.graph().get_edge_count() == 17);
                REQUIRE(transcriptome.graph().get_path_count() == 4);

                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path1")) == 11);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path2")) == 9);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path3")) == 11);
                REQUIRE(transcriptome.graph().get_step_count(transcriptome.graph().get_path_handle("path4")) == 9);

                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(1)) == "AAAA");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(2)) == "CC");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(3)) == "G");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(4)) == "T");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(5)) == "T");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(6)) == "TT");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(7)) == "T");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(8)) == "TTT");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(9)) == "C");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(10)) == "CC");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(11)) == "A");
                REQUIRE(transcriptome.graph().get_sequence(transcriptome.graph().get_handle(12)) == "AAA");             
            }
        }
    }
}
