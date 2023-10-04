/// \file minimizer_mapper.cpp
///  
/// unit tests for the minimizer mapper

#include <iostream>
#include "vg/io/json2pb.h"
#include "../io/json2graph.hpp"
#include <vg/vg.pb.h>
#include "../minimizer_mapper.hpp"
#include "../build_index.hpp"
#include "../integrated_snarl_finder.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

// We define a child class to expose protected stuff for testing
class TestMinimizerMapper : public MinimizerMapper {
public:
    TestMinimizerMapper(
        gbwtgraph::GBWTGraph gbwt_graph,
        gbwtgraph::DefaultMinimizerIndex minimizer_index,
        SnarlDistanceIndex* distance_index,
        PathPositionHandleGraph* handle_graph) 
        : MinimizerMapper(gbwt_graph, minimizer_index, distance_index, nullptr, handle_graph){};
    using MinimizerMapper::MinimizerMapper;
    using MinimizerMapper::Minimizer;
    using MinimizerMapper::fragment_length_distr;
    using MinimizerMapper::faster_cap;
    using MinimizerMapper::with_dagified_local_graph;
    using MinimizerMapper::align_sequence_between;
};

TEST_CASE("Fragment length distribution gets reasonable value", "[giraffe][mapping]") {

        vector<int64_t> distances { 27, 69, 76, 88, 107, 114, 119, 121, 124, 124, 125, 125, 125, 125, 126, 126, 126, 127, 127, 127, 127, 127, 127, 
            134, 134, 134, 134, 135, 135, 136, 136, 136, 136, 136, 136, 136, 136, 137, 137, 137, 138, 139, 139, 139, 140, 140, 140, 140, 140, 140, 
            144, 144, 144, 145, 145, 145, 145, 145, 145, 145, 145, 145, 146, 146, 146, 146, 146, 147, 147, 147, 147, 147, 148, 148, 148, 148, 148, 
            154, 154, 154, 154, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 156, 156, 156, 156, 156, 156, 156, 156, 
            160, 160, 160, 160, 160, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 162, 162, 162, 162, 163, 163, 163, 
            167, 167, 167, 167, 167, 167, 167, 168, 168, 168, 168, 169, 169, 169, 169, 169, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 171, 
            174, 174, 174, 174, 174, 175, 175, 175, 175, 175, 175, 176, 176, 176, 176, 176, 176, 176, 177, 178, 178, 178, 178, 178, 178, 178, 178,
            182, 182, 182, 182, 182, 182, 182, 182, 183, 183, 183, 183, 183, 183, 183, 183, 183, 184, 184, 184, 184, 184, 185, 185, 185, 185, 185,
            189, 189, 189, 189, 189, 189, 189, 189, 189, 189, 189, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 191, 191, 191, 
            195, 195, 195, 195, 196, 196, 196, 196, 196, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 198, 198, 198, 198, 198, 198, 198, 
            202, 202, 202, 202, 202, 202, 202, 202, 202, 203, 203, 203, 203, 203, 203, 203, 203, 204, 204, 204, 204, 204, 204, 204, 204, 204, 
            207, 207, 208, 208, 208, 208, 208, 208, 208, 208, 208, 209, 209, 209, 209, 209, 209, 210, 210, 210, 210, 210, 210, 210, 211, 211, 
            214, 214, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 216, 216, 216, 216, 216, 216, 216, 216, 216, 217, 217, 217, 217, 
            220, 220, 220, 220, 220, 220, 220, 220, 220, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 222, 222, 222, 222, 222, 
            226, 226, 226, 226, 226, 226, 226, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 228, 228, 228, 229, 229, 229, 229, 
            233, 233, 233, 233, 233, 233, 233, 233, 234, 234, 234, 234, 234, 234, 234, 234, 234, 234, 234, 235, 235, 235, 235, 235, 236, 236, 
            238, 238, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 240, 240, 240, 240, 240, 240, 240, 241, 241, 
            244, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 246, 246, 246, 246, 246, 246, 246, 246, 246, 246, 246, 247, 247, 247, 
            249, 250, 250, 250, 250, 250, 250, 250, 250, 251, 251, 251, 251, 251, 251, 251, 252, 252, 252, 252, 252, 252, 252, 252, 252, 252, 
            255, 255, 255, 255, 255, 255, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 257, 257, 257, 257, 257, 257, 257, 257, 
            261, 261, 261, 261, 261, 262, 262, 262, 262, 262, 262, 262, 262, 262, 262, 262, 263, 263, 263, 263, 264, 264, 264, 264, 264, 264, 
            268, 268, 268, 268, 268, 269, 269, 269, 269, 269, 269, 269, 269, 269, 269, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 
            274, 275, 275, 275, 276, 276, 276, 276, 276, 276, 276, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 278, 278, 278, 278, 
            281, 281, 281, 281, 281, 281, 281, 282, 282, 282, 282, 282, 282, 282, 282, 283, 283, 283, 283, 283, 283, 283, 283, 284, 284, 284, 
            288, 288, 288, 288, 288, 289, 289, 289, 289, 289, 289, 289, 290, 290, 290, 291, 291, 291, 291, 291, 291, 291, 291, 291, 292, 292, 
            295, 296, 296, 296, 296, 296, 296, 296, 296, 297, 297, 297, 297, 297, 297, 297, 298, 298, 298, 298, 298, 298, 298, 298, 298, 299,
            302, 302, 302, 302, 302, 303, 303, 303, 303, 303, 303, 303, 303, 304, 304, 304, 304, 304, 304, 304, 305, 305, 305, 305, 306, 306,
            312, 312, 312, 312, 312, 313, 313, 313, 313, 313, 314, 314, 314, 314, 314, 315, 315, 315, 315, 315, 315, 316, 316, 316, 316, 316,
            320, 320, 321, 321, 321, 321, 321, 321, 321, 322, 322, 322, 322, 322, 322, 323, 323, 323, 323, 323, 324, 324, 324, 324, 325, 325,
            332, 333, 333, 333, 333, 333, 333, 333, 333, 333, 334, 334, 334, 334, 335, 335, 335, 335, 335, 336, 336, 336, 336, 336, 337, 337,
            342, 343, 344, 344, 344, 344, 344, 344, 345, 345, 346, 346, 346, 346, 346, 347, 347, 347, 347, 348, 348, 348, 348, 348, 348, 348,
            356, 356, 357, 357, 358, 358, 359, 359, 359, 359, 361, 361, 361, 361, 362, 362, 362, 363, 363, 363, 363, 364, 364, 364, 364, 365,
            374, 374, 374, 375, 376, 376, 376, 377, 377, 377, 378, 378, 379, 379, 379, 379, 380, 380, 380, 380, 381, 381, 382, 382, 383, 384,
            397, 398, 398, 398, 399, 399, 399, 399, 399, 399, 400, 400, 402, 402, 402, 403, 403, 403, 403, 404, 404, 404, 404, 406, 407, 407,
            463, 465, 466, 466, 467, 470, 471, 473, 474, 479, 479, 480, 481, 482, 483, 485, 491, 493, 493, 495, 496, 497, 498, 501, 507, 511,
            512, 512, 515, 519, 521, 521, 521, 523, 524, 527, 531, 535, 537, 541, 542, 550, 556, 557, 557, 559, 561, 569, 572, 573, 575, 579,
            580, 585, 589, 622, 627, 647, 667, 715, 757, 780, 1302, 4224, 10520, 16912, 17605, 17773, 18141, 18234, 19908, 22806, 26071, 
            33167, 39460, 39642, 55666, 59773, 63297, 74729, 82293, 84261, 103051, 125968, 126638, 133620, 134000, 156120, 156834, 158566, 159945,
            163617, 168131, 170576, 186151, 187063, 196981, 199264, 205006, 211618, 214498, 232698, 241394, 242280, 253799, 254397, 257196, 261598,
            316979, 329457, 654834,
            332, 333, 333, 333, 333, 333, 333, 333, 333, 333, 334, 334, 334, 334, 335, 335, 335, 335, 335, 336, 336, 336, 336, 336, 337, 337,
            342, 343, 344, 344, 344, 344, 344, 344, 345, 345, 346, 346, 346, 346, 346, 347, 347, 347, 347, 348, 348, 348, 348, 348, 348, 348,
            356, 356, 357, 357, 358, 358, 359, 359, 359, 359, 361, 361, 361, 361, 362, 362, 362, 363, 363, 363, 363, 364, 364, 364, 364, 365,
            274, 275, 275, 275, 276, 276, 276, 276, 276, 276, 276, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 278, 278, 278, 278 };


        //Make an empty minimizer mapper just to get the fragment length distr
        gbwtgraph::GBWTGraph gbwt_graph;
        gbwt::GBWT gbwt;
        gbwt_graph.set_gbwt(gbwt);
        gbwtgraph::DefaultMinimizerIndex minimizer_index;
        SnarlDistanceIndex distance_index;
        PathPositionHandleGraph* handle_graph;
        TestMinimizerMapper test_mapper (gbwt_graph, minimizer_index, &distance_index, handle_graph);
        for (int64_t dist : distances) {
            if (dist <= test_mapper.max_fragment_length) {
                test_mapper.fragment_length_distr.register_fragment_length(dist);
            }
        }

        SECTION("Make sure we have a reasonable distribution") {
            //Distribution should ignore outliers
            REQUIRE(test_mapper.fragment_length_distr.std_dev() <= 400);
        }
}

/// Cover a sequence of all Gs in minimizers.
static void cover_in_minimizers(const std::string sequence, int core_width, int flank_width, int stride, std::vector<TestMinimizerMapper::Minimizer>& minimizers, std::vector<size_t>& minimizers_explored) {

    string min_seq;
    for (int i = 0; i < core_width; i++) {
        min_seq.push_back('G');
    }
    auto encoded = gbwtgraph::DefaultMinimizerIndex::key_type::encode(min_seq);
    
    for (int core_start = 0; core_start + core_width < sequence.size(); core_start += stride) {
        minimizers_explored.push_back(minimizers.size());
        minimizers.emplace_back();
        TestMinimizerMapper::Minimizer& m = minimizers.back();
        
        if (core_start <= flank_width) {
            // Partial left flank
            m.agglomeration_start = 0;
            m.agglomeration_length = core_width + flank_width + core_start;
        } else if (sequence.size() - core_start - core_width <= flank_width) {
            // Partial right flank
            m.agglomeration_start = core_start - flank_width;
            m.agglomeration_length = sequence.size() - m.agglomeration_start - 1;
        } else {
            // Full flanks
            m.agglomeration_start = core_start - flank_width;
            m.agglomeration_length = core_width + flank_width * 2;
        }
        
        // We need to set the key and its hash
        m.value.key = encoded;
        m.value.hash = m.value.key.hash();
        m.value.offset = core_start;
        m.value.is_reverse = false;
        m.length = core_width;
        
        // We knowe the occurrences won't be used.
        m.occs = nullptr;
        m.hits = 1;
        m.score = 1;
    }

}

TEST_CASE("Mapping quality cap cannot be confused by excessive Gs", "[giraffe][mapping]") {
    string sequence;
    string quality;
    for (size_t i = 0; i < 150; i++) {
        sequence.push_back('G');
        quality.push_back((char)0x1E);
    }
    
    // Cover the read in 25bp cores with 10bp flanks on each side
    int core_width = 25;
    int flank_width = 10;
    vector<TestMinimizerMapper::Minimizer> minimizers;
    // They are all going to be explored
    vector<size_t> minimizers_explored;
    
    cover_in_minimizers(sequence, core_width, flank_width, 1, minimizers, minimizers_explored);
    
    // Compute the MAPQ cap
    double cap = TestMinimizerMapper::faster_cap(minimizers, minimizers_explored, sequence, quality);
    
    // The MAPQ cap should not be infinite.
    REQUIRE(!isinf(cap));
}

TEST_CASE("Mapping quality cap cannot be confused by fuzzing with high base qualities", "[giraffe][mapping]") {
    string sequence;
    string quality;
    for (size_t i = 0; i < 100; i++) {
        sequence.push_back('G');
        quality.push_back((char)(unsigned char)60);
    }
    
    TestMinimizerMapper::Minimizer minimizer_template;
    minimizer_template.value.is_reverse = false;
    
    minimizer_template.hits = 229;
    // We knowe the occurrences won't be used.
    minimizer_template.occs = nullptr;
    minimizer_template.score = 1;
    
    for (size_t try_number = 0; try_number < 100000; try_number++) {
    
        vector<TestMinimizerMapper::Minimizer> minimizers;
        // They are all going to be explored
        vector<size_t> minimizers_explored;
        
        size_t minimizer_count = rand() % 100 + 5;
        
        for (size_t i = 0; i < minimizer_count; i++) {
            // Generate a random and not very realistic agglomeration
            size_t core_width = rand() % std::min(sequence.size()/2 - 1, (size_t)32 - 1) + 1;
            size_t run_length = rand() % std::min(sequence.size() - core_width, (size_t)32 - core_width);
            size_t flank_width = rand() % 10;
            size_t core_start = rand() % (sequence.size() - core_width - run_length);
                
            string min_seq;
            for (int i = 0; i < core_width; i++) {
                min_seq.push_back('G');
            }
            auto encoded = gbwtgraph::DefaultMinimizerIndex::key_type::encode(min_seq);
        
            minimizers_explored.push_back(minimizers.size());
            minimizers.emplace_back();
            TestMinimizerMapper::Minimizer& m = minimizers.back();
            m = minimizer_template;
            
            // Now clip the agglomeration to the read.
            m.agglomeration_start = core_start;
            m.agglomeration_length = core_width + run_length + flank_width * 2;
            if (flank_width > m.agglomeration_start) {
                m.agglomeration_length -= (flank_width - m.agglomeration_start);
                m.agglomeration_start = 0;
            } else {
                m.agglomeration_start -= flank_width;
            }
            if (m.agglomeration_start + m.agglomeration_length > sequence.size()) {
                m.agglomeration_length = sequence.size() - m.agglomeration_start;
            }
           
            // We need to set the key and its hash
            m.value.key = encoded;
            m.value.hash = m.value.key.hash();
            m.value.offset = core_start;
            m.value.is_reverse = false;
            m.length = core_width;
            
            m.hits = 229;
            // We knowe the occurrences won't be used.
            m.occs = nullptr;
            m.score = 1;
        }
        
        // Compute the MAPQ cap
        double cap = TestMinimizerMapper::faster_cap(minimizers, minimizers_explored, sequence, quality);
        
        // The MAPQ cap should not be infinite.
        REQUIRE(!isinf(cap));
    }
}

TEST_CASE("MinimizerMapper can map against subgraphs between points", "[giraffe][mapping]") {

        Aligner aligner;
        HashGraph graph;
        
        // We have a real path with a mismatch
        auto h1 = graph.create_handle("AAAAGAT");
        auto h2 = graph.create_handle("TG");
        graph.create_edge(h1, h2);
        // This node is backward
        auto h3 = graph.create_handle("AAAAAAAAATG");
        graph.create_edge(h2, graph.flip(h3));
        // And we have a dangling tip that is a better matchn
        auto h4 = graph.create_handle("TA");
        graph.create_edge(h2, h4);
        auto h5 = graph.create_handle("CA");
        graph.create_edge(h4, h5);
        
        
        Alignment aln;
        aln.set_sequence("GATTACA");
        
        // Left anchor should be on start
        pos_t left_anchor {graph.get_id(h1), false, 4};
        // Right anchor should be past end
        pos_t right_anchor {graph.get_id(h3), true, 2};
        
        TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, 20, &graph, &aligner, aln);
        
        // Make sure we get the right alignment
        REQUIRE(aln.path().mapping_size() == 3);
        REQUIRE(aln.path().mapping(0).position().node_id() == graph.get_id(h1));
        REQUIRE(aln.path().mapping(0).position().is_reverse() == graph.get_is_reverse(h1));
        REQUIRE(aln.path().mapping(0).position().offset() == offset(left_anchor));
        REQUIRE(aln.path().mapping(1).position().node_id() == graph.get_id(h2));
        REQUIRE(aln.path().mapping(1).position().is_reverse() == graph.get_is_reverse(h2));
        REQUIRE(aln.path().mapping(1).position().offset() == 0);
        REQUIRE(aln.path().mapping(2).position().node_id() == graph.get_id(h3));
        REQUIRE(aln.path().mapping(2).position().is_reverse() == !graph.get_is_reverse(h3));
        REQUIRE(aln.path().mapping(2).position().offset() == 0);
}

TEST_CASE("MinimizerMapper can map against subgraphs between abutting points", "[giraffe][mapping]") {

        Aligner aligner;
        HashGraph graph;
        
        // We have a big node
        auto h1 = graph.create_handle("AAAAGAT");
        auto h2 = graph.create_handle("TG");
        graph.create_edge(h1, h2);
        
        Alignment aln;
        aln.set_sequence("A");
        
        SECTION("Abutting points on same node") {
            // Left anchor should be on start
            pos_t left_anchor {graph.get_id(h1), false, 3};
            // Right anchor should be past end
            pos_t right_anchor {graph.get_id(h1), false, 3};
            
            TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, 20, &graph, &aligner, aln);
            
            // Make sure we get the right alignment
            REQUIRE(aln.path().mapping_size() == 1);
            REQUIRE(aln.path().mapping(0).position().node_id() == graph.get_id(h1));
            REQUIRE(aln.path().mapping(0).position().is_reverse() == graph.get_is_reverse(h1));
            REQUIRE(aln.path().mapping(0).position().offset() == offset(left_anchor));
            REQUIRE(aln.path().mapping(0).edit_size() == 1);
            REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
            REQUIRE(aln.path().mapping(0).edit(0).to_length() == 1);
            REQUIRE(aln.path().mapping(0).edit(0).sequence() == "A");
        }
        
        SECTION("Abutting points on different nodes") {
            // Left anchor should be on start
            pos_t left_anchor {graph.get_id(h1), false, 7};
            // Right anchor should be past end
            pos_t right_anchor {graph.get_id(h2), false, 0};
            
            TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, 20, &graph, &aligner, aln);
            
            // Make sure we get the right alignment
            REQUIRE(aln.path().mapping_size() == 1);
            REQUIRE(aln.path().mapping(0).position().node_id() == graph.get_id(h1));
            REQUIRE(aln.path().mapping(0).position().is_reverse() == graph.get_is_reverse(h1));
            REQUIRE(aln.path().mapping(0).position().offset() == offset(left_anchor));
            REQUIRE(aln.path().mapping(0).edit_size() == 1);
            REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
            REQUIRE(aln.path().mapping(0).edit(0).to_length() == 1);
            REQUIRE(aln.path().mapping(0).edit(0).sequence() == "A");
        }
}

TEST_CASE("MinimizerMapper can map an empty string between odd points", "[giraffe][mapping]") {

        Aligner aligner;
        
        string graph_json = R"({
            "edge": [
                {"from": "55511923", "to": "55511925"},
                {"from": "55511923", "to": "55511924"},
                {"from": "55511921", "to": "55511924"},
                {"from": "55511921", "to": "55511922"},
                {"from": "55511922", "to": "55511923"},
                {"from": "55511922", "to": "55511924"},
                {"from": "55511924", "to": "55511925"}
            ],
            "node": [
                {"id": "55511923", "sequence": "T"},
                {"id": "55511921", "sequence": "TTCCTT"},
                {"id": "55511922", "sequence": "CC"},
                {"id": "55511924", "sequence": "TC"},
                {"id": "55511925", "sequence": "CTTCCTTCC"}
            ]
        })";
        
        // TODO: Write a json_to_handle_graph
        vg::Graph proto_graph;
        json2pb(proto_graph, graph_json.c_str(), graph_json.size());
        auto graph = vg::VG(proto_graph);
        
        Alignment aln;
        aln.set_sequence("");
        
        pos_t left_anchor {55511921, false, 5}; // This is on the final base of the node
        pos_t right_anchor {55511925, false, 6};
        
        TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, 20, &graph, &aligner, aln);
        
        // Make sure we get the right alignment. We should see the last base of '21 and go '21 to '24 to '25 and delete everything
        REQUIRE(aln.path().mapping_size() == 3);
        REQUIRE(aln.path().mapping(0).position().node_id() == 55511921);
        REQUIRE(aln.path().mapping(0).position().is_reverse() == false);
        REQUIRE(aln.path().mapping(0).position().offset() == 5);
        REQUIRE(aln.path().mapping(1).position().node_id() == 55511924);
        REQUIRE(aln.path().mapping(1).position().is_reverse() == false);
        REQUIRE(aln.path().mapping(1).position().offset() == 0);
        REQUIRE(aln.path().mapping(2).position().node_id() == 55511925);
        REQUIRE(aln.path().mapping(2).position().is_reverse() == false);
        REQUIRE(aln.path().mapping(2).position().offset() == 0);
}

TEST_CASE("MinimizerMapper can align a reverse strand string to the middle of a node", "[giraffe][mapping]") {

        Aligner aligner;
        
        string graph_json = R"({
            "node": [
                {"id": "48732576", "sequence": "GACAGCTTGTGGCCTTCGTTGGAAACGGGATTTCTTCATACTATGCTAGACAGAAGAATACTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATGGTTTACACAGAGCAGATTTGAAACACTCTTTTTGTGGAATTAGCAAGTGGAGATTTCAGCCGCTTTGAGGTCAATGGTAGAAAAGGAAATATCTTCGTATAAAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTTTGAGGCCTTCGTTGGAAACGGGTTTTTTTCATATAAGGCTAGACAGAAGAATTCCTAGTAATTTCCTTGTGTTGTGTGTGTTCAACTCACAGAGTTGAACTTTCATTTACACAGAGCAGATTTGAAACACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAAGCGCTTTGAGACCAAAGGCAGAAAAGGATATATCTTCGTATAAAAACTAGACAGAATCATTCTCAGAAAATGCTCTGCGATGTGTGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCAGTTTGGAAACAATCTGTTTGTAAAGTCTGCACGTGGATAATTTGACCACTTAGAGGCCTTCGTTGGAAACGGGTTTTTTTCATGTAAGGCTAGACACAATTCTCAGTAACTTCCTTGTGTTGTGTGTATTCAACTCACAGAGTTGAACGATCCTTTACACAGAGCATACTTGGAACACTCTTTTTGTGGAAGTTGCAAGTGGAGATTTCAGCCGCTTTGAAGTCAAAGGTAGAAAAGGAAATATCTTCCTATAAAAACTAGACAGAATGATTCTCAGAAACTCCTTTGTGATGTGTGCATTCAACTCACAGAGTTTAACCTTTCTTTTCATAGAGCAGTTAGGAAACACTCTGTTTGTAAAGTCTGCAAGTGGATATTCAGACCTCTT"}
            ]
        })";
        
        vg::VG graph;
        vg::io::json2graph(graph_json, &graph);
        
        Alignment aln;
        aln.set_sequence("CAAATTCCACAAAAAGAGTGTTACAAGTCTGCTCTGTGTAAAGGATCGTTCAACTCTGGGAGTTGAATACACACAACACGCGGAAGTTACTGAGAATTCTTCTGTCTAGCCTTACATGAAAAAAACCCGTTTCCAACGAAGGCCTCAAAGAGGTCAAAATATCCACTTGCAGACTTTACAAACAGAGTGTTTCCTAACTACTCTATGAATAGAAAGGTTAAACTCTGTGAGATGAACACACACATCACAAAGGAGTTTCTGAGAATCATTCTGTCTAGTTTTTATAGGAAGATATTTCCTTTTCTACCATTGACCTCAAAGCGGCTGAAATCTCCACTTGCAAATTCCTCAAAAAGAGTGTTTCAAGTCTGCTCTGTGTAAAGGATCGTCAACTCTGTGAGTTGAATACACACAACACGCGGAAGTTACTGAGAATTCTTCTGTCTAGCATAGTATGAAGAAATCCCGTTTCCAACGAAGGCCTCAAAGAGGTCTGAATATCCACTTGCAGAGTTTACAAACAGAGTGTTTCCTAACTGCTCTATGAAAAGAAAGGTTAAACTCTGTGAGTTGAACGCACACATCACAAAGAAGTTTCTGAGAATCATCTGTCTAGTTTTTATACGAAGATATTTCCTTTTCTACCATTGACCTCAAAGCGGCTGAAATCTCCACTTGCAAATTCCACAAAAAGAGTGTTT");


        pos_t left_anchor {48732576, true, 193};
        pos_t right_anchor {48732576, true, 893};
        
        TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 800, 50, &graph, &aligner, aln);

        // We demand a positive-score alignment
        REQUIRE(aln.score() > 0);
}

TEST_CASE("MinimizerMapper can align a long tail", "[giraffe][mapping]") {

        Aligner aligner;
        
        string graph_json = R"(
            {"edge": [{"from": "28131", "to": "28132"}, {"from": "28132", "to": "28133"}, {"from": "28130", "to": "28131"}, {"from": "28129", "to": "28130"}, {"from": "28128", "to": "28129"}], "node": [{"id": "28131", "sequence": "GAATTATGATCAAATGGAATCGAATGTAATCATCATCAAATGGAATCAAAAATAACCATCATCAATTGGTATTGAATGGAATTGTCATCAAATGGAATTCAAAGGAATCATCATCAAATGGAACCGAATGGAATCCTCATTGAATGGAAATGAAAGGGGTCATCATCTAATGGAATCGCATGGAATCATCATCAAATGGAATCGAATGGAATCATCATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATGGAATCATCTTCTAATGGAAAGGAATGGAATCATCGCATAGAATCGAATGGAATTATCATCGAATGGAATCGAATGGTATCAACACCAAACGGAAAAAAACGGAATTATCGAATGGAATCGAAGAGAATCTTCGAACGGACCCGAATGGAATCATCTAATGGAATGGAATGGAATAATCCATGGACTCGAATGCAATCATCATCGAATGGAATCGAATGGAATCATCGAATGGACTCGAATGGAATAATCATTGAACGGAATCGAATGGAATCATCATCGGATGGAAACGAATGGAATCATCATCGAATGGAAATGAAAGGAGTCATCATCTAATGGAATTGCATGGAATCATCATAAAATGGAATCGAATGGAATCAACATCAAATGGAATCAAATGGAATCATTGAACGGAATTGAATGGAATCGTCATCGAATGAATTGACTGCAATCATCGAATGGTCTCGAATGGAATCATCTTCAAATGGAATGGAATGGAATCATCGCATAGAATCGAATGGAATTATCATCGAATGGAATCGAATGGAATCAACATCAAACGGAAAAAAACAGAATTATCGTATGGAATCGAAGAGAATCATCGAGTGGACCCGAATGGAATCATCTAATGGAATGGAATGGAATAATCCATGGACTCGAATGCAATCATCATCGAATGGAATCGAATGGAATCATCGAATGGACTCGAATGGAATAATCAT"}, {"id": "28132", "sequence": "TGAACGGAATCGAATGGAATCATCATCGGATGGAAATGAATGGAATCATCATCGAATGGAATCGAATAGAATTATGGAATGAAATCCAGTGTGATCATCATCGAATGGACCCGAATGGAATCATCATCCAACGGAAGCTAATGGAATCAACATCGAATGAATCGAATGGAAACACCATCGAATTGAAACGAATGGAATTATCATGAAATTGAAATGGATGGACTCATCATCGAATGGATTCGAATGGAATCATCGAATAAAATTGATTGAAATCATCATCCAATGGAATCGAATGGTATCATTGAATGGAATCGAATGGAATCATCATCAGATGGAAATGAATGGAATCGTCATAGAATGGAATCGAATGGATTCATTGAATGGAATCAGATGGAATCATCGAATGGACTGGAATGGAATCATTGAATGGACTCGAAAGGGATCATGATTGAATGGAATTGAATGGAATCATCGAATGGTCTCGATTGGAATCATTATCAAATGGAATCGAATGGAATCATCGAATAGAATCGAATGGAACAATCATCGAATGTACTCAAATGGAATTATCCTCAAATGGAATCGAATGGAATTATCGAATGCAATCGAATGGAATTATCGAATGCAATCGAATAGAATCATCGAATGGACTCGAATGGAATCATCGAATGGAATGGAATGGAACAGTCAATGAACACGAATGGAATCATCATTGAATGGAATCTAATGGAATCATCGAGTGGAATCGAATGGAATTATGATCAAATGGAATCGAATGTAATCATCATCAAATGGAATCAAAAATAACCATCATCAATTGCTATTGAATGGAATTGTCATCAAATGGAATTCAAAGGAATCATCATCAAATGGAACCGAATGGAATCCTCATTGAATGGAAATGAAAGGGGTCATCATCTAATGGAATCGCATGGAATCATCATCAAATGGAATCGAATGGAATCATCATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAAT"}, {"id": "28133", "sequence": "GAATTGAATGCAATCATCGAATGGTCTCGAATGGAATCATCTTCTAATGGAAAGGAATGGAATCATCGCATAGAATCGAATGGAATTATCATCGAATGGAATCGAATGGTATCAACACCAAACGGAAAAAAACGGAATTATCGAATGGAATCGAAGAGAATCTTCGAACGGACCCGAATGGAATCATCTAATGGAATGGAATGGAATAATCCATGGACTCGAATGCAATCATCATCGAATGGAATCGAATGGAATCATCGAATGGACTCGAATGGAATAATCATTGAACGGAATCGAATGGAATCATCATCGGATGGAAACGAATGGAATCATCATCGAATGGAAATGAAAGGAGTCATCATCTAATGCAATTGCATGGAATCATCATCAAATAGAATCGAATGGAATCAACATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAATGAATTGACTGCAATCATCGAATGGTCTCGAATGGAATCATCTTCAAATGGAATGGAATGGAATCATCGCATAGAATCGAATGGAATTATCATCGAATGGAATCGAATGGAATCAACAACAAACGGAAAAAAACGGAATTATCGAATGGAATCGAAGAGAATCATCGAATGGACCCGAATGGAATCATCTAATGGAATGGAATGGAATAATCCATGGACTCGAATGCAATCATCATCGAATGGAATCGAATGGAATCATCGAATGGACTCGAATGGAATAATCATTGAACGGAATCGAATGGAATCATCATCGGATGGAAATGAATGGAATCATCATCGAATGGAATCGAATAGAATTATGGAATGAAATCCAGTGTGGTCATCATCGAATGGACCCGAATGGAATCATCATCCAACGGAAGCTAATGGAATCAACATCGAATGAATCAAATGGAAACACCATCGAATTGAAACGAATGGAATTATCATGAAATTGAAACGGATGGACTCATCATCGAATGGATTCGAATGGAATCATCGAATAAAATTGATTGAAA"}, {"id": "28130", "sequence": "ATCATCGAATGGTCTCGAATGGAATCATCTTCTAATGGAAAGGAATGGAATCATCGCATAGAATCGAATGGAATTATCATCGAATGGAATCGAATGGTATCAACACCAAACGGAAAAAAACGGAATTATCGAAAGGAATCGAAGAGAATCTTCGAACGGACCCGAATGGAATCATCTAATGGAATGGAATGGAATAATCCATGGACTCGAATGCAATCATCATCGAATGGAATCGAATGGAATCATCGAATGGACTCGAATGGAATAATCATTGAACGGAATCGAATGGAATCATCATCGGATGGAAATGAATGGAATCATCATCGAATGGAATCGAATAGAATTATGGAATGAAATCCAGTGTGATCATCATCGAATGGACCCGAATGGAATCATCATCCAACAGAAGCTAATGGAATCAACATCGAATGAATCGAATGGAAACACCATCGAATTCAAACGAATGGAATTACCATGAAATTGAAATGGATGGACTCATCATCGAATGGATTCGGATGGAATCATCGAATAAAATTGATTGAAATCATCATCGAATGGAATCGAATGGTATCATTGAATGGAATCGAATGGAATCATCATCAGATGGAAATGAATGGAATCGTCATAGAATGGAATCGAATGGATTCATTGAATGGAATCAGATGGAATCATCGAATGGACTGGAATGGAATCATTGAATGGACTCGAAAGGGATCATGATTGAATGGAATTGAATGGAATCATCGAATGGTCTCGATTGGAATCATTATGAAATGGAATCGAATGGAATCACCGAATAGAATCGAATGGAACAATCATCGAATGGACTCAAATGGAATTATCCTCAAATGGAATCGAATGGAATTATCAAATGCAATCGAATGGAATTATCGAATGCAATCGAATAGAATCATCGAATGGACTCGAATGGAATCATCGAATGGAATGGAATGGAACAGTCAATGAACTCGAATGGAATCATCATTGAATGGAATCGAATGTAATCATCCAGTGGAATCGAATG"}, {"id": "28129", "sequence": "CTCGATTGGAATCATTATCAAATGGAATCGAATGGAATCACCGAATAGAATCGAATGGAACAATCATCGAATGGACTCAAATGGAATTATCCTCAAATGGAATCGAATGGAATTATCGAATGCAATCGAATGGAATTATCGAATGCAATCGAATAGAATCATCGAATGGACTCGAATGGAATCATCGAATGGAATGGAATGGAACAGTCAATGAACACGAATGGAATCATCATTGAATGGAATCGAATGGAATCATCGAGTGGAATCGAATGGAATTATGATCAAATGGAATCGAATGTAATCATCATCAAATGGAATCAAAAATAACCATCATCAATTGGTATTGAATGGAATTGTCATCAAATGGAATTCAAAGGAATCATCATCAAATGGAACCGAATGGAATCCTCATTGAATGGAAATGAAAGGGGTCATCATCTAATGGAATCGCATGGAATCATCACCAAATGGAATCGAATGGAATCATCATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATGGAATCATCTTCTAATGGAAAGGAATGGAATCATCGCATAGAATCGAATGGAATTATCATCGAATGGAATCGAATGGTATCAACACCAAACGGAAAAAAACGGAATTATCGAATGGAATCGAAGAGAATCTTCGAACGGACCCGAATGGAATCATCTAATGGAATGGAATGGAATAATCCATGGACTCGAATGCAATCATCATCGAATGGAATCGAATGGAATCATCGAATGGACTCGAATGGAATAATCATTGAACGGAATCGAATGGAATCATCATCGGATGGAAACGAATGGAATCATCATCGAATGGAAATGAAAGGAGTCATCATCTAATGCAATTGCATGGAATCATCATCAAATGGAATCGAATGGAATCAACATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATCGTCATCGAATGAATTGACTGCA"}, {"id": "28128", "sequence": "ATCCATGGACTCGAATGCAATCATCATCGAATGGAATCGAATGGAATCATCAAATGGACTCGAATGGAATAATCATTGAACGGAATCGAATGGAATCATCATCGGATGGAAACGAATGGAATCATCATCGAATGGAAATGAAAGGAGTCATCATCTAATGGAATTGCATGGAATCATCATAAAATGGAATCGAATGGAATCAATATCAAATGGAATCAAATGGAATCATTGAACGGAATTGAATGGAATCGTCATCGAATGAATTGACTGCAATCATCGAATGGTCTCGAATGGAATCATCTTCAAATGGAATGGAATGGAATCATCGCATAGAATCGAATGGAATTATCATCGAATGGAATCGAATGGAATCAACATCAAACGGAAAAAAACGGAATTATCGAATGGAATCGAAGAGAATCATCGAATGGACCCGAATGGAATCATCTAATGGAATGGAATGGAATAATCCATGGACTCGAATGCAATCATCATCGAATGGAATCGAATGGAATCATCGAATGGACTCGAATGGAATAATCATTGAACGGAATCGAATGGAATCATCATCGGATGGAAATGAATGGAATCATCATCGAATGGAATCGAATAGAATTATGGAATGAAATCCAGTGTGATCATCATCGAATGGACCCGAATGGAATCATCATCCAACGGAAGCTAATGGAATCAACATCGAATGAATCGAATGGAAACACCATCGAATTGAAACGAATGGAATTATCATGAAATTGAAATGGATGGACTCATCATCGAATGGATTCGAATGGAATCATCGAATAAAATTGATTGAAATCATCATCGAATGGAATCGAATGGTATCATTGAATGGAATCGAATGGAATCATCATCAGATGGAAATGAATGGAATCGTCATAGAATGGAATCGAATGGATTCATTGAATGGAATCAGATGGAATCATCGAATGGACTGGAATGGAATCATTGAATGGACTCGAAAGGGATCATGATTGAATGGAATTGAATGGAATCATCGAATGGT"}]}
        )";
        
        vg::VG graph;
        vg::io::json2graph(graph_json, &graph);
        
        Alignment aln;
        aln.set_sequence("TGGATGATGATTCCATTTGGGTCCATTCGATGATGATCACACTGGATTTCATTCCATAATTCTATTCGATTCCATTCGATGATGATTCCATACATTTCCATCCGATGATGATTCCATTCGATTCCGTTCAATGATTATTCCATTCGAGTCCATTCGATGATTCCATTCGATTCCATTCGATGATGATTGCATTCGAGTCCATGGATTATTCCATTCCATTCCATTAGGTGATTCCATTCGGGTCCGTTCGAAGATTCTCTTCGATTCCATTCGATAATTCCGTTTTTTTCCGTTTGATGTTGATTCCATTCGACTCCATTCGATGATAATTCCACTCGATTCTATGCGATGATTCCATTCCATTCCATTTGAAGATGATTCCATTCGAGACCATTCGATGATTGCATTCAATTCATTCGATGACGATTCCATTCAATTCCGTTCAATGATTCCATTTGATTCCATTTGATGTTGATTCCATTCGATTCCATTTTATGATGATTCCATGCAATTCCATTAGATGATGACTCCTTTCATTTCCATTCGATGATGATTCCATTCGGTTCCATTTGATGATGATTCCTTTGAATTCCGTTTGATGACAATTCCATTCAATACCAATTGATGATGGTTATTTTTGATTCCATTTGATGAGGATTACATTCGATTCCATTGGATCATAATTCCATTCGATTCCACTCGATGATTCCATTCGATTCCATTCAATGATGATTCCATTCGAGTTCATTGACTGTTCCATTCCATTCCATTCGATGATTCCATTCGAGTCCATTCGATGATTCTATTCGATTGCATTCGATAATTCCATTCGATTGCATTCGATAATTCCCTTCGATTCCATTTGAGGATAATTCCATTTGAGTCCATTCGATGATTGTTCCATTCGATTCTATTCGGTGATTCCATTCGATTCCATTTGATAATGATTCCAATCGAGACCATTCGATGATTCCATTCAATTCCATTCAACAATGATTCCATTCGAGTCCATTCAATGATTCCATTCCAGTCCATTCGATGATTCCATCTGACTCCATTCAATGAATCCATTCGATTCCATTCTATGACGATTCCATTCATTTCCATCTGATGATGATTCCATTCGATCCCATCCAATGACACCATTCGATTCCATTCGATGATGATTTCAATCAATTTTATTCGATGATTCCATTCGAATCCATTCGATGATGGGTCCATCCATTTCAATTTCATGATAATTCCATTCGTTTCAATTCGATGGTTTTTCCATTCGATTCATTCGATGTTGATTCCATTAGCTTCCGTTGGATGATGATTCCATTCGGGTCCATTCGATGATGATCACACTGGATTTCATTCCATAATTCTATTCGATTCCATTCGATGATGATTCCATTCATTTCCATCCGATGATGATTCCATTCGATTCCGTTCAATGATTATTCCATTCGAGTCCATTCGATGATTCCATTCGATTCCATTCGATGATGATTGCATTCGAGTCCATGGATTATTCCATTCCATTCCATTAGATGATTCCATTCGGGTCCGTTCGAAGATTCTCTTCGATTCCATTCGATAATCCCGTTTTTTTCCGTTTGATATTGATACCATTCGATTCCATTCAATGATAATTCCATTCGATTCTATGCGATGATTCCATTCCATTCCATTGGAAGATGATTCCATTCGAGACCATTCGATGATTGCATTCAATTCATTCGATGACGATTCCATTCAATTCCGTTCAATGATTCCATTTGATTCCATTTGATGTTGATTCCATTCGATTCCATTTGATGATGATTCCATGCAATTCCATTAGATGATGACTCCTTTCATTTCCATTCGATGATGATTCCATTCGTTTCCATCCGAAGATGATTCCATTCGATTCCGTTCAATGATTATTCCATTCGAGTCCATTCGATGATTCCATTCGATTCTATACGATGATGATTGCATTCGAGTCCGTGGATCATTCCATTCAATTCCATTAGATTATTCCATTCGAGTCCATTCGATGATTCTCTTCGATTACATTCGACGATGATTGCATTCGAGTCCATGGATTATTCCATTCCATTCCATTAGATGATTCCATTCGGGTCCATTCGATGATTCTCTTCGATTCCATTCGATAATTCCGTTTTTTTCCGTTTGATGTTGATTCCATTCGATTCCATTCGATGATAATTCCATTCGATTCTATGCGATGATTCCATTCCATTCCATTTGAAGATGATTCCATTCGAGACCATTCGATGATTGCATTCAATTCATTCGATGACGATTCCATTCAATTCCGTTCAATGATTCCATTTGATTCCATTTGATGTTGATTCCATTCGATTCCATTTTATGATGATTCAATGCAATTCCATTAGATGATGACTCCTTTCATTTACATTCGATGATGATTCCATTCGTTTCCATCCGATGATGATTCCATTCGATTCTCTTCAATGCTTATTCCATTCGAGTCCATTCGATGATTCCATTCGATTCCATTCGATGATGATTGCATTCGAGTCCATGGATTATTCCATTCAATTCCATTAGATGATTCCATTCGGGTCCGTTCGAAGATTCTCTTCGATTCCATTCGATAATTCCGTTTTTCTCCGTTTGGTGTTGATACCATTCGATTCCATTCGATGATAATTCCTTTCGATTCTATGCGATGATTCCATTCCTTTCCATTAGAAGACGATTCCATTCGAGACCATTCGATGATTGCATTCAATTCATTCGATGACGATTCCATTCAATTCTGTTCAATGATTCCATCAGATTCCATTTGATGATGATTCCATTCGATTCCATTTGATGATGATTCCATGCGATTCCATTAGATGATGACCCCTTTCATTTCCATTCAATGAGGATTCCATTCGGTTCCATTTCATGATGTTTCCTTTGAATTCCATTTGATGACAATTCCATTCAATACCAATTGATGATGGTTATTTTTGATTCCATTTGATGATGATTACATTCGATTCCATTTGATCATAATTCCATTCGATTCCACTCGATGATTCCATTCGATTCCATTCAATGATGATTCCATTCGAGTTCATTGACTGTTCCATTCCATTCCATTCGATGATTCCATTCGAGTCCATTCGATGATTCTATTCGATTGCATTCGATAATTCCATTCGATTGCATTCGATAATTCCATTCGATTCCATTGGAGGATAATTCCATTTGAGTCCATTCGATGATTGTTCCATTCGATTCTATTCGGTGATTCCATTCGATTCCATTTGATAATGATTCCAATCGAGACCATTCGATGATTCCATTCAATTCCATTCAATAATGATCCCTTTCGAGTCCATTCAATGATTCCATTCCAGTCCATTCGATGATTCCATCTGATTCCATTCAATGAATCCATTCGATTCCATTCTATGACGATTCCATTCATTTCCATCTGATGATGATTACATTCGATCCCATTCAATGACACCATTAGATTCCATTCGATGATGATTTCAATCAATTTTATTCGATGATTCCATTCGAATCCATTCGATGATGGGTCCATCCATTTCAATTTCATGATAATTCCATTCGTTTCAATTCGATGGTGTTTCCATTCGATTCATTCGATGTTGATTCCATTAGCTTCCGTTGGATGATGATTCCATTCGGGTACATTCGATGATGATCACACTGGATTTCATTCCATAATTCTATTCGATTCCATTCGATGATGATTCCATTCATTTCCATCCGATGATGATTCCATTCGATTCCGTTCAATGATTATTCCATTCGAGTCCATTCGATGATTCCATTCGATTCCATTCGATGATGATTGCATTCGAGTCCATGGATTATTCCATTCCATTCCATTAGATGATTCCATTCGGGTCCGTTCGAAGATTCTCTTCGATTCCATTCGATAATTCCGTTTTTTTCCGTTTGATGTTGATACCATTCGATTCCATTCGATGATAATTC");


        pos_t left_anchor {28132, true, 892};
        
        TestMinimizerMapper::align_sequence_between(left_anchor, empty_pos_t(), 5000, 500, &graph, &aligner, aln);

        // We demand a positive-score alignment
        REQUIRE(aln.score() > 0);
        // We demand not having a very long softclip at the end
        REQUIRE(aln.path().mapping_size() > 0);
        auto& last_mapping = aln.path().mapping(aln.path().mapping_size() - 1);
        REQUIRE(last_mapping.edit_size() > 0);
        auto& last_edit = last_mapping.edit(last_mapping.edit_size() - 1);
        REQUIRE(last_edit.to_length() <= std::max(10, last_edit.from_length()));
}

TEST_CASE("MinimizerMapper can align a long left tail", "[giraffe][mapping]") {

        Aligner aligner;
        
        string graph_json = R"(
{
  "edge": [
    {
      "from": "56",
      "to": "58"
    },
    {
      "from": "56",
      "to": "57"
    },
    {
      "from": "35",
      "to": "36"
    },
    {
      "from": "60",
      "to": "61"
    },
    {
      "from": "220",
      "to": "221"
    },
    {
      "from": "220",
      "to": "222"
    },
    {
      "from": "308",
      "to": "309"
    },
    {
      "from": "308",
      "to": "310"
    },
    {
      "from": "67",
      "from_start": true,
      "to": "69"
    },
    {
      "from": "215",
      "to": "217"
    },
    {
      "from": "73",
      "to": "75"
    },
    {
      "from": "319",
      "to": "320"
    },
    {
      "from": "251",
      "to": "252"
    },
    {
      "from": "251",
      "to": "253"
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
      "from": "348",
      "to": "349"
    },
    {
      "from": "185",
      "to": "186"
    },
    {
      "from": "185",
      "to": "187"
    },
    {
      "from": "365",
      "to": "367"
    },
    {
      "from": "333",
      "to": "334"
    },
    {
      "from": "86",
      "to": "87"
    },
    {
      "from": "168",
      "to": "170"
    },
    {
      "from": "364",
      "to": "365"
    },
    {
      "from": "364",
      "to": "366"
    },
    {
      "from": "207",
      "to": "209"
    },
    {
      "from": "263",
      "to": "264"
    },
    {
      "from": "263",
      "to": "265"
    },
    {
      "from": "242",
      "to": "244"
    },
    {
      "from": "183",
      "to": "185"
    },
    {
      "from": "376",
      "to": "377"
    },
    {
      "from": "224",
      "to": "226"
    },
    {
      "from": "177",
      "to": "178"
    },
    {
      "from": "12",
      "to": "13"
    },
    {
      "from": "75",
      "to": "76",
      "to_end": true
    },
    {
      "from": "75",
      "to": "77"
    },
    {
      "from": "111",
      "to": "113"
    },
    {
      "from": "23",
      "to": "24"
    },
    {
      "from": "264",
      "to": "266"
    },
    {
      "from": "41",
      "to": "42"
    },
    {
      "from": "41",
      "to": "43"
    },
    {
      "from": "68",
      "to": "69"
    },
    {
      "from": "82",
      "from_start": true,
      "to": "84"
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
      "from": "125",
      "to": "127"
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
      "from": "339",
      "to": "340"
    },
    {
      "from": "66",
      "to": "68"
    },
    {
      "from": "66",
      "to": "67",
      "to_end": true
    },
    {
      "from": "103",
      "to": "104"
    },
    {
      "from": "280",
      "to": "281"
    },
    {
      "from": "59",
      "to": "61"
    },
    {
      "from": "208",
      "to": "209"
    },
    {
      "from": "336",
      "to": "337"
    },
    {
      "from": "26",
      "to": "27"
    },
    {
      "from": "358",
      "to": "359"
    },
    {
      "from": "358",
      "to": "360"
    },
    {
      "from": "366",
      "to": "367"
    },
    {
      "from": "211",
      "to": "212"
    },
    {
      "from": "343",
      "to": "344"
    },
    {
      "from": "343",
      "to": "345"
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
      "from": "116",
      "to": "118"
    },
    {
      "from": "100",
      "to": "101"
    },
    {
      "from": "230",
      "to": "232"
    },
    {
      "from": "279",
      "to": "281"
    },
    {
      "from": "79",
      "from_start": true,
      "to": "81"
    },
    {
      "from": "195",
      "to": "197"
    },
    {
      "from": "374",
      "to": "375"
    },
    {
      "from": "141",
      "to": "143"
    },
    {
      "from": "278",
      "to": "280"
    },
    {
      "from": "278",
      "to": "279"
    },
    {
      "from": "135",
      "to": "137"
    },
    {
      "from": "138",
      "to": "140"
    },
    {
      "from": "222",
      "to": "223"
    },
    {
      "from": "107",
      "to": "109"
    },
    {
      "from": "107",
      "to": "108"
    },
    {
      "from": "46",
      "to": "47"
    },
    {
      "from": "276",
      "to": "278"
    },
    {
      "from": "295",
      "to": "296"
    },
    {
      "from": "57",
      "to": "58"
    },
    {
      "from": "381",
      "to": "383"
    },
    {
      "from": "247",
      "to": "248"
    },
    {
      "from": "152",
      "to": "153"
    },
    {
      "from": "152",
      "to": "154"
    },
    {
      "from": "170",
      "to": "171"
    },
    {
      "from": "170",
      "to": "172"
    },
    {
      "from": "129",
      "to": "131"
    },
    {
      "from": "250",
      "to": "251"
    },
    {
      "from": "238",
      "to": "239"
    },
    {
      "from": "238",
      "to": "240"
    },
    {
      "from": "78",
      "to": "80"
    },
    {
      "from": "78",
      "to": "79",
      "to_end": true
    },
    {
      "from": "133",
      "to": "134"
    },
    {
      "from": "258",
      "to": "260"
    },
    {
      "from": "72",
      "to": "73"
    },
    {
      "from": "72",
      "to": "74"
    },
    {
      "from": "184",
      "to": "185"
    },
    {
      "from": "252",
      "to": "254"
    },
    {
      "from": "1",
      "to": "2"
    },
    {
      "from": "1",
      "to": "6"
    },
    {
      "from": "137",
      "to": "138"
    },
    {
      "from": "137",
      "to": "139"
    },
    {
      "from": "154",
      "to": "155"
    },
    {
      "from": "154",
      "to": "156"
    },
    {
      "from": "22",
      "from_start": true,
      "to": "24"
    },
    {
      "from": "313",
      "to": "314"
    },
    {
      "from": "237",
      "to": "238"
    },
    {
      "from": "206",
      "to": "207"
    },
    {
      "from": "206",
      "to": "208"
    },
    {
      "from": "288",
      "to": "290"
    },
    {
      "from": "270",
      "to": "272"
    },
    {
      "from": "354",
      "to": "355"
    },
    {
      "from": "299",
      "to": "301"
    },
    {
      "from": "299",
      "to": "300"
    },
    {
      "from": "33",
      "to": "35"
    },
    {
      "from": "33",
      "to": "34",
      "to_end": true
    },
    {
      "from": "345",
      "to": "346"
    },
    {
      "from": "40",
      "to": "41"
    },
    {
      "from": "231",
      "to": "232"
    },
    {
      "from": "113",
      "to": "114"
    },
    {
      "from": "113",
      "to": "115"
    },
    {
      "from": "245",
      "to": "246"
    },
    {
      "from": "254",
      "to": "256"
    },
    {
      "from": "254",
      "to": "255"
    },
    {
      "from": "283",
      "to": "284"
    },
    {
      "from": "165",
      "to": "167"
    },
    {
      "from": "309",
      "to": "313"
    },
    {
      "from": "142",
      "to": "143"
    },
    {
      "from": "5",
      "to": "6"
    },
    {
      "from": "114",
      "to": "116"
    },
    {
      "from": "55",
      "to": "56"
    },
    {
      "from": "265",
      "to": "266"
    },
    {
      "from": "325",
      "to": "326"
    },
    {
      "from": "136",
      "to": "137"
    },
    {
      "from": "117",
      "to": "119"
    },
    {
      "from": "45",
      "from_start": true,
      "to": "47"
    },
    {
      "from": "145",
      "to": "146"
    },
    {
      "from": "282",
      "to": "284"
    },
    {
      "from": "337",
      "to": "339"
    },
    {
      "from": "337",
      "to": "338"
    },
    {
      "from": "342",
      "to": "343"
    },
    {
      "from": "275",
      "to": "277"
    },
    {
      "from": "275",
      "to": "276"
    },
    {
      "from": "363",
      "to": "364"
    },
    {
      "from": "378",
      "to": "380"
    },
    {
      "from": "351",
      "to": "352"
    },
    {
      "from": "158",
      "to": "160"
    },
    {
      "from": "218",
      "to": "220"
    },
    {
      "from": "176",
      "to": "178"
    },
    {
      "from": "176",
      "to": "177"
    },
    {
      "from": "28",
      "from_start": true,
      "to": "30"
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
      "from": "92",
      "to": "94"
    },
    {
      "from": "36",
      "to": "37"
    },
    {
      "from": "36",
      "to": "38"
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
      "from": "84",
      "to": "86"
    },
    {
      "from": "7",
      "from_start": true,
      "to": "8"
    },
    {
      "from": "25",
      "from_start": true,
      "to": "27"
    },
    {
      "from": "203",
      "to": "204"
    },
    {
      "from": "203",
      "to": "205"
    },
    {
      "from": "95",
      "to": "96"
    },
    {
      "from": "95",
      "to": "97"
    },
    {
      "from": "292",
      "to": "293"
    },
    {
      "from": "353",
      "to": "355"
    },
    {
      "from": "232",
      "to": "233"
    },
    {
      "from": "232",
      "to": "234"
    },
    {
      "from": "93",
      "to": "95"
    },
    {
      "from": "296",
      "to": "298"
    },
    {
      "from": "296",
      "to": "297"
    },
    {
      "from": "304",
      "to": "305"
    },
    {
      "from": "18",
      "to": "19"
    },
    {
      "from": "240",
      "to": "241"
    },
    {
      "from": "147",
      "to": "149"
    },
    {
      "from": "157",
      "to": "158"
    },
    {
      "from": "157",
      "to": "159"
    },
    {
      "from": "16",
      "to": "18"
    },
    {
      "from": "16",
      "to": "17",
      "to_end": true
    },
    {
      "from": "370",
      "to": "372"
    },
    {
      "from": "341",
      "to": "343"
    },
    {
      "from": "287",
      "to": "289"
    },
    {
      "from": "287",
      "to": "288"
    },
    {
      "from": "349",
      "to": "350"
    },
    {
      "from": "349",
      "to": "351"
    },
    {
      "from": "19",
      "to": "20"
    },
    {
      "from": "19",
      "to": "21"
    },
    {
      "from": "44",
      "to": "46"
    },
    {
      "from": "44",
      "to": "45",
      "to_end": true
    },
    {
      "from": "368",
      "to": "369"
    },
    {
      "from": "217",
      "to": "219"
    },
    {
      "from": "217",
      "to": "218"
    },
    {
      "from": "31",
      "from_start": true,
      "to": "33"
    },
    {
      "from": "266",
      "to": "267"
    },
    {
      "from": "266",
      "to": "268"
    },
    {
      "from": "146",
      "to": "147"
    },
    {
      "from": "146",
      "to": "148"
    },
    {
      "from": "74",
      "to": "75"
    },
    {
      "from": "61",
      "to": "63"
    },
    {
      "from": "61",
      "to": "62"
    },
    {
      "from": "29",
      "to": "30"
    },
    {
      "from": "380",
      "to": "382"
    },
    {
      "from": "380",
      "to": "381"
    },
    {
      "from": "212",
      "to": "213"
    },
    {
      "from": "212",
      "to": "214"
    },
    {
      "from": "303",
      "to": "305"
    },
    {
      "from": "228",
      "from_start": true,
      "to": "229"
    },
    {
      "from": "159",
      "to": "160"
    },
    {
      "from": "193",
      "to": "194"
    },
    {
      "from": "226",
      "to": "227"
    },
    {
      "from": "226",
      "to": "228",
      "to_end": true
    },
    {
      "from": "101",
      "to": "102"
    },
    {
      "from": "101",
      "to": "103"
    },
    {
      "from": "360",
      "to": "361"
    },
    {
      "from": "223",
      "to": "225"
    },
    {
      "from": "223",
      "to": "224"
    },
    {
      "from": "105",
      "to": "107"
    },
    {
      "from": "285",
      "to": "287"
    },
    {
      "from": "17",
      "from_start": true,
      "to": "19"
    },
    {
      "from": "271",
      "to": "272"
    },
    {
      "from": "335",
      "to": "337"
    },
    {
      "from": "198",
      "to": "200"
    },
    {
      "from": "166",
      "to": "167"
    },
    {
      "from": "214",
      "to": "215"
    },
    {
      "from": "214",
      "to": "216"
    },
    {
      "from": "331",
      "to": "332"
    },
    {
      "from": "331",
      "to": "333"
    },
    {
      "from": "80",
      "to": "81"
    },
    {
      "from": "51",
      "from_start": true,
      "to": "53"
    },
    {
      "from": "89",
      "to": "90"
    },
    {
      "from": "274",
      "to": "275"
    },
    {
      "from": "246",
      "to": "247"
    },
    {
      "from": "246",
      "to": "248"
    },
    {
      "from": "143",
      "to": "144"
    },
    {
      "from": "143",
      "to": "145"
    },
    {
      "from": "48",
      "from_start": true,
      "to": "50"
    },
    {
      "from": "15",
      "to": "16"
    },
    {
      "from": "97",
      "to": "98"
    },
    {
      "from": "330",
      "to": "331"
    },
    {
      "from": "284",
      "to": "286"
    },
    {
      "from": "284",
      "to": "285"
    },
    {
      "from": "134",
      "to": "136"
    },
    {
      "from": "134",
      "to": "135"
    },
    {
      "from": "110",
      "to": "112"
    },
    {
      "from": "110",
      "to": "111"
    },
    {
      "from": "30",
      "to": "31",
      "to_end": true
    },
    {
      "from": "30",
      "to": "32"
    },
    {
      "from": "6",
      "to": "8"
    },
    {
      "from": "6",
      "to": "7",
      "to_end": true
    },
    {
      "from": "234",
      "to": "235"
    },
    {
      "from": "219",
      "to": "220"
    },
    {
      "from": "367",
      "to": "368"
    },
    {
      "from": "367",
      "to": "369"
    },
    {
      "from": "272",
      "to": "273"
    },
    {
      "from": "272",
      "to": "274"
    },
    {
      "from": "182",
      "to": "183"
    },
    {
      "from": "182",
      "to": "184"
    },
    {
      "from": "253",
      "to": "254"
    },
    {
      "from": "153",
      "to": "156"
    },
    {
      "from": "186",
      "to": "188"
    },
    {
      "from": "164",
      "to": "165"
    },
    {
      "from": "164",
      "to": "166"
    },
    {
      "from": "64",
      "to": "66"
    },
    {
      "from": "64",
      "to": "65"
    },
    {
      "from": "267",
      "to": "269"
    },
    {
      "from": "90",
      "to": "91"
    },
    {
      "from": "90",
      "to": "92"
    },
    {
      "from": "139",
      "to": "140"
    },
    {
      "from": "4",
      "from_start": true,
      "to": "5"
    },
    {
      "from": "359",
      "to": "361"
    },
    {
      "from": "13",
      "to": "14",
      "to_end": true
    },
    {
      "from": "13",
      "to": "15"
    },
    {
      "from": "104",
      "to": "106"
    },
    {
      "from": "104",
      "to": "105"
    },
    {
      "from": "316",
      "to": "317"
    },
    {
      "from": "328",
      "to": "329"
    },
    {
      "from": "328",
      "to": "330"
    },
    {
      "from": "52",
      "to": "53"
    },
    {
      "from": "179",
      "to": "180"
    },
    {
      "from": "179",
      "to": "181"
    },
    {
      "from": "369",
      "to": "370"
    },
    {
      "from": "369",
      "to": "371"
    },
    {
      "from": "356",
      "to": "358"
    },
    {
      "from": "300",
      "to": "302"
    },
    {
      "from": "43",
      "to": "44"
    },
    {
      "from": "11",
      "from_start": true,
      "to": "13"
    },
    {
      "from": "69",
      "to": "70"
    },
    {
      "from": "69",
      "to": "71"
    },
    {
      "from": "171",
      "to": "173"
    },
    {
      "from": "302",
      "to": "303"
    },
    {
      "from": "302",
      "to": "304"
    },
    {
      "from": "85",
      "to": "87"
    },
    {
      "from": "119",
      "to": "120"
    },
    {
      "from": "119",
      "to": "121"
    },
    {
      "from": "39",
      "to": "41"
    },
    {
      "from": "39",
      "to": "40"
    },
    {
      "from": "216",
      "to": "217"
    },
    {
      "from": "126",
      "to": "128"
    },
    {
      "from": "108",
      "to": "110"
    },
    {
      "from": "382",
      "to": "383"
    },
    {
      "from": "156",
      "to": "157"
    },
    {
      "from": "124",
      "to": "125"
    },
    {
      "from": "27",
      "to": "29"
    },
    {
      "from": "27",
      "to": "28",
      "to_end": true
    },
    {
      "from": "10",
      "to": "12"
    },
    {
      "from": "10",
      "to": "11",
      "to_end": true
    },
    {
      "from": "261",
      "to": "263"
    },
    {
      "from": "307",
      "to": "308"
    },
    {
      "from": "2",
      "to": "4",
      "to_end": true
    },
    {
      "from": "2",
      "to": "3"
    },
    {
      "from": "144",
      "to": "146"
    },
    {
      "from": "273",
      "to": "275"
    },
    {
      "from": "257",
      "to": "259"
    },
    {
      "from": "257",
      "to": "258"
    },
    {
      "from": "352",
      "to": "353"
    },
    {
      "from": "352",
      "to": "354"
    },
    {
      "from": "312",
      "to": "314"
    },
    {
      "from": "200",
      "to": "201"
    },
    {
      "from": "200",
      "to": "202"
    },
    {
      "from": "81",
      "to": "82",
      "to_end": true
    },
    {
      "from": "81",
      "to": "83"
    },
    {
      "from": "20",
      "to": "21"
    },
    {
      "from": "290",
      "to": "292"
    },
    {
      "from": "290",
      "to": "291"
    },
    {
      "from": "340",
      "to": "341"
    },
    {
      "from": "340",
      "to": "342"
    },
    {
      "from": "187",
      "to": "188"
    },
    {
      "from": "213",
      "to": "214"
    },
    {
      "from": "329",
      "to": "331"
    },
    {
      "from": "9",
      "to": "10"
    },
    {
      "from": "346",
      "to": "348"
    },
    {
      "from": "346",
      "to": "347"
    },
    {
      "from": "189",
      "to": "191"
    },
    {
      "from": "344",
      "to": "346"
    },
    {
      "from": "227",
      "to": "229"
    },
    {
      "from": "294",
      "to": "296"
    },
    {
      "from": "109",
      "to": "110"
    },
    {
      "from": "161",
      "to": "163"
    },
    {
      "from": "249",
      "to": "251"
    },
    {
      "from": "372",
      "to": "374"
    },
    {
      "from": "372",
      "to": "373"
    },
    {
      "from": "241",
      "to": "242"
    },
    {
      "from": "241",
      "to": "243"
    },
    {
      "from": "88",
      "from_start": true,
      "to": "90"
    },
    {
      "from": "209",
      "to": "211"
    },
    {
      "from": "209",
      "to": "210"
    },
    {
      "from": "236",
      "to": "238"
    },
    {
      "from": "120",
      "to": "122"
    },
    {
      "from": "323",
      "to": "324"
    },
    {
      "from": "323",
      "to": "325"
    },
    {
      "from": "260",
      "to": "261"
    },
    {
      "from": "260",
      "to": "262"
    },
    {
      "from": "297",
      "to": "299"
    },
    {
      "from": "24",
      "to": "26"
    },
    {
      "from": "24",
      "to": "25",
      "to_end": true
    },
    {
      "from": "8",
      "to": "9"
    },
    {
      "from": "8",
      "to": "10"
    },
    {
      "from": "37",
      "to": "39"
    },
    {
      "from": "83",
      "to": "84"
    },
    {
      "from": "190",
      "to": "191"
    },
    {
      "from": "201",
      "to": "203"
    },
    {
      "from": "99",
      "to": "101"
    },
    {
      "from": "121",
      "to": "122"
    },
    {
      "from": "311",
      "to": "313"
    },
    {
      "from": "281",
      "to": "282"
    },
    {
      "from": "281",
      "to": "283"
    },
    {
      "from": "14",
      "from_start": true,
      "to": "16"
    },
    {
      "from": "314",
      "to": "315"
    },
    {
      "from": "314",
      "to": "316"
    },
    {
      "from": "357",
      "to": "358"
    },
    {
      "from": "334",
      "to": "335"
    },
    {
      "from": "334",
      "to": "336"
    },
    {
      "from": "174",
      "to": "176"
    },
    {
      "from": "322",
      "to": "323"
    },
    {
      "from": "269",
      "to": "270"
    },
    {
      "from": "269",
      "to": "271"
    },
    {
      "from": "315",
      "to": "317"
    },
    {
      "from": "123",
      "to": "125"
    },
    {
      "from": "305",
      "to": "306"
    },
    {
      "from": "305",
      "to": "307"
    },
    {
      "from": "268",
      "to": "269"
    },
    {
      "from": "32",
      "to": "33"
    },
    {
      "from": "197",
      "to": "199"
    },
    {
      "from": "197",
      "to": "198"
    },
    {
      "from": "233",
      "to": "235"
    },
    {
      "from": "196",
      "to": "197"
    },
    {
      "from": "262",
      "to": "263"
    },
    {
      "from": "320",
      "to": "322"
    },
    {
      "from": "320",
      "to": "321"
    },
    {
      "from": "324",
      "to": "326"
    },
    {
      "from": "210",
      "to": "212"
    },
    {
      "from": "151",
      "to": "152"
    },
    {
      "from": "239",
      "to": "241"
    },
    {
      "from": "63",
      "to": "64"
    },
    {
      "from": "54",
      "to": "55"
    },
    {
      "from": "54",
      "to": "56"
    },
    {
      "from": "191",
      "to": "193"
    },
    {
      "from": "191",
      "to": "192"
    },
    {
      "from": "91",
      "to": "92"
    },
    {
      "from": "244",
      "to": "246"
    },
    {
      "from": "244",
      "to": "245"
    },
    {
      "from": "205",
      "to": "206"
    },
    {
      "from": "62",
      "to": "64"
    },
    {
      "from": "150",
      "to": "152"
    },
    {
      "from": "327",
      "to": "328"
    },
    {
      "from": "122",
      "to": "124"
    },
    {
      "from": "122",
      "to": "123"
    },
    {
      "from": "58",
      "to": "59"
    },
    {
      "from": "58",
      "to": "60"
    },
    {
      "from": "199",
      "to": "200"
    },
    {
      "from": "173",
      "to": "174"
    },
    {
      "from": "173",
      "to": "175"
    },
    {
      "from": "256",
      "to": "257"
    },
    {
      "from": "188",
      "to": "189"
    },
    {
      "from": "188",
      "to": "190"
    },
    {
      "from": "277",
      "to": "278"
    },
    {
      "from": "361",
      "to": "362"
    },
    {
      "from": "361",
      "to": "363"
    },
    {
      "from": "98",
      "to": "100"
    },
    {
      "from": "98",
      "to": "99"
    },
    {
      "from": "355",
      "to": "357"
    },
    {
      "from": "355",
      "to": "356"
    },
    {
      "from": "235",
      "to": "237"
    },
    {
      "from": "235",
      "to": "236"
    },
    {
      "from": "204",
      "to": "206"
    },
    {
      "from": "377",
      "to": "379"
    },
    {
      "from": "377",
      "to": "378"
    },
    {
      "from": "310",
      "to": "311"
    },
    {
      "from": "310",
      "to": "312"
    },
    {
      "from": "321",
      "to": "323"
    },
    {
      "from": "371",
      "to": "372"
    },
    {
      "from": "76",
      "from_start": true,
      "to": "78"
    },
    {
      "from": "34",
      "from_start": true,
      "to": "36"
    },
    {
      "from": "318",
      "to": "320"
    },
    {
      "from": "243",
      "to": "244"
    },
    {
      "from": "50",
      "to": "52"
    },
    {
      "from": "50",
      "to": "51",
      "to_end": true
    },
    {
      "from": "194",
      "to": "196"
    },
    {
      "from": "194",
      "to": "195"
    },
    {
      "from": "167",
      "to": "169",
      "to_end": true
    },
    {
      "from": "167",
      "to": "168"
    },
    {
      "from": "301",
      "to": "302"
    },
    {
      "from": "317",
      "to": "319"
    },
    {
      "from": "317",
      "to": "318"
    },
    {
      "from": "132",
      "to": "134"
    },
    {
      "from": "140",
      "to": "142"
    },
    {
      "from": "140",
      "to": "141"
    },
    {
      "from": "202",
      "to": "203"
    },
    {
      "from": "248",
      "to": "250"
    },
    {
      "from": "248",
      "to": "249"
    },
    {
      "from": "169",
      "from_start": true,
      "to": "170"
    },
    {
      "from": "42",
      "to": "44"
    },
    {
      "from": "180",
      "to": "182"
    },
    {
      "from": "255",
      "to": "257"
    },
    {
      "from": "160",
      "to": "161"
    },
    {
      "from": "160",
      "to": "162"
    },
    {
      "from": "87",
      "to": "88",
      "to_end": true
    },
    {
      "from": "87",
      "to": "89"
    },
    {
      "from": "289",
      "to": "290"
    },
    {
      "from": "49",
      "to": "50"
    },
    {
      "from": "291",
      "to": "293"
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
      "from": "225",
      "to": "226"
    },
    {
      "from": "128",
      "to": "130"
    },
    {
      "from": "128",
      "to": "129"
    },
    {
      "from": "347",
      "to": "349"
    },
    {
      "from": "259",
      "to": "260"
    },
    {
      "from": "350",
      "to": "352"
    },
    {
      "from": "379",
      "to": "380"
    },
    {
      "from": "375",
      "to": "377"
    },
    {
      "from": "375",
      "to": "376"
    },
    {
      "from": "21",
      "to": "22",
      "to_end": true
    },
    {
      "from": "21",
      "to": "23"
    },
    {
      "from": "229",
      "to": "231"
    },
    {
      "from": "229",
      "to": "230"
    },
    {
      "from": "38",
      "to": "39"
    },
    {
      "from": "163",
      "to": "179"
    },
    {
      "from": "163",
      "to": "164"
    },
    {
      "from": "332",
      "to": "334"
    },
    {
      "from": "131",
      "to": "132"
    },
    {
      "from": "131",
      "to": "133"
    },
    {
      "from": "102",
      "to": "104"
    },
    {
      "from": "192",
      "to": "194"
    },
    {
      "from": "70",
      "to": "72"
    },
    {
      "from": "326",
      "to": "328"
    },
    {
      "from": "326",
      "to": "327"
    },
    {
      "from": "221",
      "to": "223"
    },
    {
      "from": "373",
      "to": "375"
    },
    {
      "from": "53",
      "to": "56"
    },
    {
      "from": "53",
      "to": "54"
    },
    {
      "from": "362",
      "to": "364"
    },
    {
      "from": "47",
      "to": "48",
      "to_end": true
    },
    {
      "from": "47",
      "to": "49"
    },
    {
      "from": "175",
      "to": "176"
    },
    {
      "from": "286",
      "to": "287"
    },
    {
      "from": "338",
      "to": "340"
    },
    {
      "from": "178",
      "to": "179"
    },
    {
      "from": "3",
      "to": "5"
    },
    {
      "from": "96",
      "to": "98"
    },
    {
      "from": "306",
      "to": "308"
    },
    {
      "from": "149",
      "to": "151"
    },
    {
      "from": "149",
      "to": "150"
    },
    {
      "from": "155",
      "to": "157"
    },
    {
      "from": "181",
      "to": "182"
    },
    {
      "from": "65",
      "to": "66"
    },
    {
      "from": "293",
      "to": "294"
    },
    {
      "from": "293",
      "to": "295"
    },
    {
      "from": "298",
      "to": "299"
    }
  ],
  "node": [
    {
      "id": "56",
      "sequence": "GTGTAGTGGAGTGAAGTGGGTTCGACTGGAATGGAATTGAACGGAATGGAATGGAATTTAATGGAATGGAATGGAATGGAATGGAA"
    },
    {
      "id": "35",
      "sequence": "A"
    },
    {
      "id": "60",
      "sequence": "C"
    },
    {
      "id": "220",
      "sequence": "TGGAGTGAAGTTGAATGAAAGAATGGAATGGAATGGAGTGGA"
    },
    {
      "id": "308",
      "sequence": "TGGAATGGAATGGAAT"
    },
    {
      "id": "67",
      "sequence": "G"
    },
    {
      "id": "215",
      "sequence": "G"
    },
    {
      "id": "73",
      "sequence": "G"
    },
    {
      "id": "319",
      "sequence": "A"
    },
    {
      "id": "251",
      "sequence": "A"
    },
    {
      "id": "115",
      "sequence": "TG"
    },
    {
      "id": "112",
      "sequence": "A"
    },
    {
      "id": "348",
      "sequence": "GA"
    },
    {
      "id": "185",
      "sequence": "TGGAATT"
    },
    {
      "id": "365",
      "sequence": "A"
    },
    {
      "id": "333",
      "sequence": "AC"
    },
    {
      "id": "86",
      "sequence": "G"
    },
    {
      "id": "168",
      "sequence": "T"
    },
    {
      "id": "364",
      "sequence": "TGGAATGGA"
    },
    {
      "id": "207",
      "sequence": "A"
    },
    {
      "id": "263",
      "sequence": "TGGA"
    },
    {
      "id": "242",
      "sequence": "G"
    },
    {
      "id": "183",
      "sequence": "T"
    },
    {
      "id": "376",
      "sequence": "GGAAT"
    },
    {
      "id": "224",
      "sequence": "C"
    },
    {
      "id": "177",
      "sequence": "TCCAT"
    },
    {
      "id": "12",
      "sequence": "A"
    },
    {
      "id": "75",
      "sequence": "AATGTAATGGCATGAAATAGAATGGAATGGAATGGAGTGGAATGGAGTGGAGTAGAATGGAATGGAGCGGAATGGATTGAAGTGGAGTGGAATGCAATGGAGTGGAATGGAGTGGAGAGAAACGGAACGGAATGGATTCCTGTGGAAAGAATGAATTGGAATGCATTGGAGTGGATTGGAGAGGAATGGAGTGGAGGGCAATGGAAA"
    },
    {
      "id": "111",
      "sequence": "T"
    },
    {
      "id": "23",
      "sequence": "A"
    },
    {
      "id": "264",
      "sequence": "G"
    },
    {
      "id": "41",
      "sequence": "AGTGGAGTGGAATGGAATGGAGTGATATGGAATGGAGTGGAATGGAATGGCATCGAATGGAATGAAATAGAAGGGAATGGAATGGAATGGAA"
    },
    {
      "id": "68",
      "sequence": "A"
    },
    {
      "id": "82",
      "sequence": "C"
    },
    {
      "id": "130",
      "sequence": "AT"
    },
    {
      "id": "125",
      "sequence": "TGGAATGGA"
    },
    {
      "id": "77",
      "sequence": "C"
    },
    {
      "id": "172",
      "sequence": "A"
    },
    {
      "id": "71",
      "sequence": "G"
    },
    {
      "id": "339",
      "sequence": "A"
    },
    {
      "id": "66",
      "sequence": "AGTGGAATGGAATGGAA"
    },
    {
      "id": "103",
      "sequence": "A"
    },
    {
      "id": "280",
      "sequence": "A"
    },
    {
      "id": "59",
      "sequence": "T"
    },
    {
      "id": "208",
      "sequence": "G"
    },
    {
      "id": "336",
      "sequence": "C"
    },
    {
      "id": "26",
      "sequence": "G"
    },
    {
      "id": "358",
      "sequence": "TGGAAT"
    },
    {
      "id": "366",
      "sequence": "G"
    },
    {
      "id": "211",
      "sequence": "A"
    },
    {
      "id": "343",
      "sequence": "TGGA"
    },
    {
      "id": "127",
      "sequence": "AT"
    },
    {
      "id": "116",
      "sequence": "G"
    },
    {
      "id": "100",
      "sequence": "ATA"
    },
    {
      "id": "230",
      "sequence": "G"
    },
    {
      "id": "279",
      "sequence": "G"
    },
    {
      "id": "79",
      "sequence": "T"
    },
    {
      "id": "195",
      "sequence": "G"
    },
    {
      "id": "374",
      "sequence": "T"
    },
    {
      "id": "141",
      "sequence": "T"
    },
    {
      "id": "278",
      "sequence": "A"
    },
    {
      "id": "135",
      "sequence": "A"
    },
    {
      "id": "138",
      "sequence": "A"
    },
    {
      "id": "222",
      "sequence": "A"
    },
    {
      "id": "107",
      "sequence": "AATGGAG"
    },
    {
      "id": "46",
      "sequence": "G"
    },
    {
      "id": "276",
      "sequence": "G"
    },
    {
      "id": "295",
      "sequence": "T"
    },
    {
      "id": "57",
      "sequence": "TGGAA"
    },
    {
      "id": "381",
      "sequence": "GGA"
    },
    {
      "id": "247",
      "sequence": "GGAAT"
    },
    {
      "id": "152",
      "sequence": "GG"
    },
    {
      "id": "170",
      "sequence": "GA"
    },
    {
      "id": "129",
      "sequence": "CA"
    },
    {
      "id": "250",
      "sequence": "G"
    },
    {
      "id": "238",
      "sequence": "GAATGGAATG"
    },
    {
      "id": "78",
      "sequence": "GAGAGGAATGGAACAGAGTGGAATGGAGTTGAGTGGAGTGGGATAGATTGGAGTGTAATGGAGTTTAGTGGAGAGGAATGGAATAGAGTGGAATGGAGTTG"
    },
    {
      "id": "133",
      "sequence": "G"
    },
    {
      "id": "258",
      "sequence": "G"
    },
    {
      "id": "72",
      "sequence": "GATTGGAATGGAATGAAGTG"
    },
    {
      "id": "184",
      "sequence": "A"
    },
    {
      "id": "252",
      "sequence": "G"
    },
    {
      "id": "1",
      "sequence": "ATGGAGTGGTGTGAAATGAAAAGGAATGGAATGGAATGGAATGGATTGGAAAAGAATGGAATGGAGGGGAATGGAATGGAATGGAAGGGACTGGAATGGCTTCGAGTGGAGTGTAGTGGAATGGAGTGGAATAGAATGGAAAGGAGTGGAATGGAATCGAATGAGTGGAACGGAATGGAATGCAATGGAATGGAATGGAATGGAATGTAGTGGAGCAGAGTGGAATGGAATGGAATGGAATATAGAGTAGTGGAATGGAATGGAATGGAATGCAATGGAATGGA"
    },
    {
      "id": "137",
      "sequence": "TG"
    },
    {
      "id": "154",
      "sequence": "CTGGGA"
    },
    {
      "id": "22",
      "sequence": "C"
    },
    {
      "id": "313",
      "sequence": "A"
    },
    {
      "id": "237",
      "sequence": "G"
    },
    {
      "id": "206",
      "sequence": "TG"
    },
    {
      "id": "288",
      "sequence": "GG"
    },
    {
      "id": "270",
      "sequence": "A"
    },
    {
      "id": "354",
      "sequence": "C"
    },
    {
      "id": "299",
      "sequence": "AATGGAATGGAATGGAATGGAATGGAATGGAA"
    },
    {
      "id": "33",
      "sequence": "GGAGTGGAATGGATTGGAGAGGAGTGGAGTACATTGGAATGGAGTGGAATGGAGTGAAGTGCAATGGAATGGAATGGAATGAGTGGAGTGGAATGGAATGGAGTGGAACGGAGTGGAGGGGAATGGAATGGAGTGGAAAGGAATGGAGTGGAATGGATTGGAGTGGAGTGGAGTCGAATGGAATGGAGTGAAATGGAGTGGAGCGTAATTGAATGGAAAGGTGTGGAGTTGAGTGGAATGGAA"
    },
    {
      "id": "345",
      "sequence": "A"
    },
    {
      "id": "40",
      "sequence": "AGTGGAGTGGAATGG"
    },
    {
      "id": "231",
      "sequence": "A"
    },
    {
      "id": "113",
      "sequence": "G"
    },
    {
      "id": "245",
      "sequence": "TGAAA"
    },
    {
      "id": "254",
      "sequence": "TG"
    },
    {
      "id": "283",
      "sequence": "G"
    },
    {
      "id": "165",
      "sequence": "G"
    },
    {
      "id": "309",
      "sequence": "TCC"
    },
    {
      "id": "142",
      "sequence": "A"
    },
    {
      "id": "5",
      "sequence": "GAATGGAATGGAATGCAATGGAATGGA"
    },
    {
      "id": "114",
      "sequence": "CA"
    },
    {
      "id": "55",
      "sequence": "ATGGAATGGAATGGA"
    },
    {
      "id": "265",
      "sequence": "A"
    },
    {
      "id": "325",
      "sequence": "GGAATG"
    },
    {
      "id": "136",
      "sequence": "G"
    },
    {
      "id": "117",
      "sequence": "T"
    },
    {
      "id": "45",
      "sequence": "T"
    },
    {
      "id": "145",
      "sequence": "G"
    },
    {
      "id": "282",
      "sequence": "A"
    },
    {
      "id": "337",
      "sequence": "ATGGAATGGA"
    },
    {
      "id": "342",
      "sequence": "A"
    },
    {
      "id": "275",
      "sequence": "G"
    },
    {
      "id": "363",
      "sequence": "A"
    },
    {
      "id": "378",
      "sequence": "C"
    },
    {
      "id": "351",
      "sequence": "A"
    },
    {
      "id": "158",
      "sequence": "TT"
    },
    {
      "id": "218",
      "sequence": "G"
    },
    {
      "id": "176",
      "sequence": "T"
    },
    {
      "id": "28",
      "sequence": "C"
    },
    {
      "id": "148",
      "sequence": "C"
    },
    {
      "id": "92",
      "sequence": "TGGAATT"
    },
    {
      "id": "36",
      "sequence": "GGACTG"
    },
    {
      "id": "118",
      "sequence": "AGACTG"
    },
    {
      "id": "162",
      "sequence": "A"
    },
    {
      "id": "84",
      "sequence": "AGTGGAATAGAGTGGAATGTAATATAACGGTGTGTAGTGGAATGGAATGCAATGGAATGAAATGGAATGAAATAAAAAGGAATGGAACTAAGTGTAGTGGAGTGGAATGTAATTGAGTGGAGTGGAATGGAATAAATTGGAATGGAATGCATTGGAGTGGAGTGGAGGTGAGTGGAAGGGAATGGATCGGAATGGAACGGACGGGAATGGATTGGAATGGAATGGAGGGGAATGGAATGGCATGGAATGGATTTGAATGTAAT"
    },
    {
      "id": "7",
      "sequence": "TCCATTCCATTTCATTCCATTCCAT"
    },
    {
      "id": "25",
      "sequence": "T"
    },
    {
      "id": "203",
      "sequence": "TGGA"
    },
    {
      "id": "95",
      "sequence": "AGTGGA"
    },
    {
      "id": "292",
      "sequence": "A"
    },
    {
      "id": "353",
      "sequence": "T"
    },
    {
      "id": "232",
      "sequence": "TG"
    },
    {
      "id": "93",
      "sequence": "G"
    },
    {
      "id": "296",
      "sequence": "G"
    },
    {
      "id": "304",
      "sequence": "G"
    },
    {
      "id": "18",
      "sequence": "A"
    },
    {
      "id": "240",
      "sequence": "G"
    },
    {
      "id": "147",
      "sequence": "A"
    },
    {
      "id": "157",
      "sequence": "T"
    },
    {
      "id": "16",
      "sequence": "GATTGGAGAGGAATGGATTGGAGTGGAATCGACTGGAGTGGAATGGAAAGGATTGGAGTGGACAGGAATGGAATGAAGTGGATTGGAGTGGAGTGGAACAGAATGGAACGGAGTGCAGTGGAGTAGAATGGAATGGAGTGGAACGGAATGGAGTGGAAGAGAATGGAGTGGGGCAGAGTGGAGTGGACTCGAATGGAATGGAATGGAGTGGAATGGATTGGAACGAAATGGGAAGGAATGGATTGGAGTGGAATAGAATGGAGTGGGATGGAATGAAGTGGAATGGAATGGAGAGGAGTGGAG"
    },
    {
      "id": "370",
      "sequence": "A"
    },
    {
      "id": "341",
      "sequence": "G"
    },
    {
      "id": "287",
      "sequence": "TGGAAT"
    },
    {
      "id": "349",
      "sequence": "GGA"
    },
    {
      "id": "19",
      "sequence": "GGAATAGAATGGAGTGAAATACAGTAGAGTGGAATGGAATGGAATGTAGTGGAGAGGAATGGAATTGAATGGAATGGAATTCAGAGGAATGAAGTGGAGTGGAGTGGAATGGAATGGA"
    },
    {
      "id": "44",
      "sequence": "GGAATGGAGTGGAGCGGAATGGAATGGAATGGAATGCAATGGAATGGAGTGGAGTGGAATGGAATGGAATGCAAAGGAATGGACTGGAACGGAGTGGAGTGGAGCGGAATGTAATGGAGACGATTGGGGTAGAAAGGAACGGAATGGAATGGAGTGGAGTGGAATGGAGTTGAGTGGATTGCAATGGAAAGGAATGGAATGGAGTGATATGGAATGGTGAGGAAGGGAGTGGATTGGAAAGGAATGGAGAGCAACGAATTGGAGTGGAGTGGATTGGAATGGAATGTAGAGGAACTGAACGGAAAGGAGTGGATTGAAATGGAATGGAATGGAACAGAATGGAAAGGAACATAAAGAAATGGAATGGAATGCAATGGAGTGGGGTGGAGGTTAATGGAATAGAGTGGAGAGGAATAGAATGGAATGGAAAAGAAT"
    },
    {
      "id": "368",
      "sequence": "G"
    },
    {
      "id": "217",
      "sequence": "TGGAATGGA"
    },
    {
      "id": "31",
      "sequence": "C"
    },
    {
      "id": "266",
      "sequence": "TGGA"
    },
    {
      "id": "146",
      "sequence": "GAATTC"
    },
    {
      "id": "74",
      "sequence": "C"
    },
    {
      "id": "61",
      "sequence": "AAT"
    },
    {
      "id": "29",
      "sequence": "A"
    },
    {
      "id": "380",
      "sequence": "AATGGAAT"
    },
    {
      "id": "212",
      "sequence": "T"
    },
    {
      "id": "303",
      "sequence": "A"
    },
    {
      "id": "228",
      "sequence": "TGGA"
    },
    {
      "id": "159",
      "sequence": "GG"
    },
    {
      "id": "193",
      "sequence": "AT"
    },
    {
      "id": "226",
      "sequence": "GGAAT"
    },
    {
      "id": "101",
      "sequence": "ATG"
    },
    {
      "id": "360",
      "sequence": "G"
    },
    {
      "id": "223",
      "sequence": "TGGAATGGAATGGAA"
    },
    {
      "id": "105",
      "sequence": "C"
    },
    {
      "id": "285",
      "sequence": "G"
    },
    {
      "id": "17",
      "sequence": "A"
    },
    {
      "id": "271",
      "sequence": "G"
    },
    {
      "id": "335",
      "sequence": "A"
    },
    {
      "id": "198",
      "sequence": "T"
    },
    {
      "id": "166",
      "sequence": "A"
    },
    {
      "id": "214",
      "sequence": "GGA"
    },
    {
      "id": "331",
      "sequence": "TGGGAAAGAATGGAATGGAGTGC"
    },
    {
      "id": "80",
      "sequence": "G"
    },
    {
      "id": "51",
      "sequence": "T"
    },
    {
      "id": "89",
      "sequence": "T"
    },
    {
      "id": "274",
      "sequence": "GA"
    },
    {
      "id": "246",
      "sequence": "GGAATGGAATGGAATGGAATGGAAT"
    },
    {
      "id": "143",
      "sequence": "GGAAT"
    },
    {
      "id": "48",
      "sequence": "C"
    },
    {
      "id": "15",
      "sequence": "G"
    },
    {
      "id": "97",
      "sequence": "CT"
    },
    {
      "id": "330",
      "sequence": "A"
    },
    {
      "id": "284",
      "sequence": "A"
    },
    {
      "id": "134",
      "sequence": "TGGA"
    },
    {
      "id": "110",
      "sequence": "GGAGTGG"
    },
    {
      "id": "30",
      "sequence": "AGTGGAATAGAATGGAATGGAGACGAATTGAATGGATTGACTTGAATGGAGTGGAATAAAGTCCAGTGGAATGGAAAGGAGAGGAATGGGA"
    },
    {
      "id": "6",
      "sequence": "ATGGAGTGGA"
    },
    {
      "id": "234",
      "sequence": "G"
    },
    {
      "id": "219",
      "sequence": "A"
    },
    {
      "id": "367",
      "sequence": "TGGAATGGAATGGAATG"
    },
    {
      "id": "272",
      "sequence": "TGGAATGGAATGGA"
    },
    {
      "id": "182",
      "sequence": "GA"
    },
    {
      "id": "253",
      "sequence": "A"
    },
    {
      "id": "153",
      "sequence": "AATTCC"
    },
    {
      "id": "186",
      "sequence": "TA"
    },
    {
      "id": "164",
      "sequence": "GGAATGGA"
    },
    {
      "id": "64",
      "sequence": "CGATGGGGGG"
    },
    {
      "id": "267",
      "sequence": "G"
    },
    {
      "id": "90",
      "sequence": "GTGGAGTGAAGTGGAGTGTAGAGGAGTCGAGTGGATGGGACTGGAATGGAATGGAGTGGAAAGGTGTGGAGTGGAAAGGAATGGA"
    },
    {
      "id": "139",
      "sequence": "T"
    },
    {
      "id": "4",
      "sequence": "C"
    },
    {
      "id": "359",
      "sequence": "A"
    },
    {
      "id": "13",
      "sequence": "AGTAGAGTGGAGTGAAATGTTGTGGAGTGGAGTGGAATGGAGTAAAATGGAATGGAATGAAGTGGAGTGGAATGGAATGGAGTGGAATGTAACGGAGT"
    },
    {
      "id": "104",
      "sequence": "AATG"
    },
    {
      "id": "316",
      "sequence": "A"
    },
    {
      "id": "328",
      "sequence": "GGA"
    },
    {
      "id": "52",
      "sequence": "G"
    },
    {
      "id": "179",
      "sequence": "GGAAT"
    },
    {
      "id": "369",
      "sequence": "A"
    },
    {
      "id": "356",
      "sequence": "G"
    },
    {
      "id": "300",
      "sequence": "T"
    },
    {
      "id": "43",
      "sequence": "C"
    },
    {
      "id": "11",
      "sequence": "A"
    },
    {
      "id": "69",
      "sequence": "GGAAA"
    },
    {
      "id": "171",
      "sequence": "G"
    },
    {
      "id": "302",
      "sequence": "GGA"
    },
    {
      "id": "85",
      "sequence": "T"
    },
    {
      "id": "119",
      "sequence": "GT"
    },
    {
      "id": "39",
      "sequence": "AATGCAATGGAGTGGAATGGATTGAAGTGGAATGGAATGGAGTGGAGTGGAGAGGAATGGAATGGAGTGGAATGCAGTGG"
    },
    {
      "id": "216",
      "sequence": "A"
    },
    {
      "id": "126",
      "sequence": "GG"
    },
    {
      "id": "108",
      "sequence": "A"
    },
    {
      "id": "382",
      "sequence": "TCC"
    },
    {
      "id": "156",
      "sequence": "A"
    },
    {
      "id": "124",
      "sequence": "G"
    },
    {
      "id": "27",
      "sequence": "AATGGAATGGAGTAGCATAGAATGAAATGGAATGGAGTGGGGTGGAGTGGAGTGGAATTGACTGGAGTGGTATAGAATGCAATGGAATGGAGAGGAGGGCAGTGGAGTGGAGTGGGGTC"
    },
    {
      "id": "10",
      "sequence": "AGGTATGGAGTGGAGGGGAGTGGATTGGAGTGGAGAGGAATGGAGTGGAATCTTGTTCAATGGAGTGGAATATAATGGAATCAAGTGGAGTGGAATGGATTGGAGTGGAGTGGAATGGAGTGGAGTGGAGAGGAATGGAATGGAGTGGAATGCAGTGGAGTGGAGTGGAATGGAGGGCAGTGGAATGGAATGGATAGGAGTGGAGTGGAGAGGACTGGACTTGTGTGGAATGGAATGGAATGGAATGGAGTGGGATTGAGAGGAGTGGAGTGGAGTAGAATGGATTGCACTGGAATGGAATGGAATGGAATTCAGTTGAATGGAATAGATTGGAATGGAACGGAGTTCAATGGAATGGAGAGTAATGAAGTGGAGTGGAGAGGAGTGGAATGGAATGGAGTGGAATGGAGTGGAGTGGAATGGAATAAAGTGGAATGGAGTGGATTGGAACGGAATGGAATGGAATGGATTCAAGTGGTGTGGGTGGAATGGAATGAAATGGAATGGAGTGGACAGAAGTGGAGTGGAATGCATTGGAATGGAGTGGCTTCGAATGGTGTCGGTGGAATGGAAGGAAATGAAATGGAGTGAAGTGGAATGGAGTGGAATGCAATTGTTTGGAGTGGTGTGGAGAT"
    },
    {
      "id": "261",
      "sequence": "A"
    },
    {
      "id": "307",
      "sequence": "A"
    },
    {
      "id": "2",
      "sequence": "ATGGAGTGGAAT"
    },
    {
      "id": "144",
      "sequence": "T"
    },
    {
      "id": "273",
      "sequence": "AT"
    },
    {
      "id": "257",
      "sequence": "AATG"
    },
    {
      "id": "352",
      "sequence": "TGGAA"
    },
    {
      "id": "312",
      "sequence": "TG"
    },
    {
      "id": "200",
      "sequence": "TT"
    },
    {
      "id": "81",
      "sequence": "ATAGATTGGAATGGAATGGAATGCAATCGAATGGATTGGAATGGAATGGAATGGAATGGAAATGAGTGGAGTGGAGTGAAATGGAATGCAGTTCAATGGAGGGGAGAGAAATGGAAAGGAATGGAATGGAATGAGGCGGTGTGAAATGAAATGCAGTGGAATTGAATAGAGTGGAATGGAATGGATTGGAGGGGATTGGAATGGAATGGAGTTGAATGGAATATAGTGTAATGGAATG"
    },
    {
      "id": "20",
      "sequence": "ATGGA"
    },
    {
      "id": "290",
      "sequence": "AAT"
    },
    {
      "id": "340",
      "sequence": "TGGAATGGA"
    },
    {
      "id": "187",
      "sequence": "CC"
    },
    {
      "id": "213",
      "sequence": "GGAATTGACTGGAATGGAATGGAGCGGAAAGCAGTGGAGT"
    },
    {
      "id": "329",
      "sequence": "T"
    },
    {
      "id": "9",
      "sequence": "TGGAG"
    },
    {
      "id": "346",
      "sequence": "TGGA"
    },
    {
      "id": "189",
      "sequence": "T"
    },
    {
      "id": "344",
      "sequence": "G"
    },
    {
      "id": "227",
      "sequence": "GGTG"
    },
    {
      "id": "294",
      "sequence": "A"
    },
    {
      "id": "109",
      "sequence": "T"
    },
    {
      "id": "161",
      "sequence": "G"
    },
    {
      "id": "249",
      "sequence": "C"
    },
    {
      "id": "383",
      "sequence": "ATGGAA"
    },
    {
      "id": "372",
      "sequence": "TGGAA"
    },
    {
      "id": "241",
      "sequence": "A"
    },
    {
      "id": "88",
      "sequence": "T"
    },
    {
      "id": "209",
      "sequence": "A"
    },
    {
      "id": "236",
      "sequence": "A"
    },
    {
      "id": "120",
      "sequence": "T"
    },
    {
      "id": "323",
      "sequence": "GGAATGGAATGGAAT"
    },
    {
      "id": "260",
      "sequence": "A"
    },
    {
      "id": "297",
      "sequence": "G"
    },
    {
      "id": "24",
      "sequence": "TGGAATGGAATGGAATCTAATGGAAAGGAATGGAATGGAAAGGACTGGAGTTGAAAGGAATTGAGAGGAATGAAATGGACTAGAATGTCATGGAATGGAATGGAATGTAGTGGATTTCAATGGAATGTAATAGAATAGAGTGGAATGTAGTTGTGTGGAGTGCAGTGGAATGGAAAGTTGTGGATTGGGGTGGAGGGGAATGGTGTGGAAAGAATGGAGTGCAGTGGAGTGGAATGGAGGGTAGTGGAGTGGAATGGAAAGGAATAGAATCGAAACGAATTGAATGGAATGGAATGCAGAAGACAGGAGTGGAGTGGAATTGATTGGAGTGGAATGTAGCGGAGTGGAGTGGATTGGAATGGAATGCAAAGGAATGGAATGGAAACGAGTACAATGGAATGGAAAGGAACGGAATGAAGTGGGGTGGAGTGGAATGGAATGGAGTGGAATGCAGTTGAGTAAAGTGGATTGGAATGGAATGTAGTGGAATG"
    },
    {
      "id": "8",
      "sequence": "G"
    },
    {
      "id": "37",
      "sequence": "C"
    },
    {
      "id": "83",
      "sequence": "C"
    },
    {
      "id": "190",
      "sequence": "A"
    },
    {
      "id": "201",
      "sequence": "GAG"
    },
    {
      "id": "99",
      "sequence": "GGC"
    },
    {
      "id": "121",
      "sequence": "C"
    },
    {
      "id": "311",
      "sequence": "A"
    },
    {
      "id": "281",
      "sequence": "TG"
    },
    {
      "id": "14",
      "sequence": "T"
    },
    {
      "id": "314",
      "sequence": "TG"
    },
    {
      "id": "357",
      "sequence": "A"
    },
    {
      "id": "334",
      "sequence": "TGG"
    },
    {
      "id": "174",
      "sequence": "G"
    },
    {
      "id": "322",
      "sequence": "A"
    },
    {
      "id": "269",
      "sequence": "TGGA"
    },
    {
      "id": "315",
      "sequence": "G"
    },
    {
      "id": "123",
      "sequence": "A"
    },
    {
      "id": "305",
      "sequence": "TGGAATGGA"
    },
    {
      "id": "268",
      "sequence": "A"
    },
    {
      "id": "32",
      "sequence": "A"
    },
    {
      "id": "197",
      "sequence": "TGGAATGGA"
    },
    {
      "id": "233",
      "sequence": "T"
    },
    {
      "id": "196",
      "sequence": "A"
    },
    {
      "id": "262",
      "sequence": "G"
    },
    {
      "id": "320",
      "sequence": "AA"
    },
    {
      "id": "324",
      "sequence": "A"
    },
    {
      "id": "210",
      "sequence": "G"
    },
    {
      "id": "151",
      "sequence": "AT"
    },
    {
      "id": "239",
      "sequence": "C"
    },
    {
      "id": "63",
      "sequence": "G"
    },
    {
      "id": "54",
      "sequence": "ATGGA"
    },
    {
      "id": "191",
      "sequence": "TGGA"
    },
    {
      "id": "91",
      "sequence": "ATGGAATGGAGTCGTG"
    },
    {
      "id": "244",
      "sequence": "TGGAAT"
    },
    {
      "id": "205",
      "sequence": "A"
    },
    {
      "id": "62",
      "sequence": "T"
    },
    {
      "id": "150",
      "sequence": "GA"
    },
    {
      "id": "327",
      "sequence": "TCCAT"
    },
    {
      "id": "122",
      "sequence": "GA"
    },
    {
      "id": "58",
      "sequence": "ATGAATA"
    },
    {
      "id": "199",
      "sequence": "A"
    },
    {
      "id": "173",
      "sequence": "TGGA"
    },
    {
      "id": "256",
      "sequence": "A"
    },
    {
      "id": "188",
      "sequence": "ATGGA"
    },
    {
      "id": "277",
      "sequence": "C"
    },
    {
      "id": "361",
      "sequence": "GA"
    },
    {
      "id": "98",
      "sequence": "GGAAT"
    },
    {
      "id": "355",
      "sequence": "GGA"
    },
    {
      "id": "235",
      "sequence": "AATGGAAT"
    },
    {
      "id": "204",
      "sequence": "G"
    },
    {
      "id": "377",
      "sequence": "G"
    },
    {
      "id": "310",
      "sequence": "GG"
    },
    {
      "id": "321",
      "sequence": "T"
    },
    {
      "id": "371",
      "sequence": "G"
    },
    {
      "id": "76",
      "sequence": "C"
    },
    {
      "id": "34",
      "sequence": "A"
    },
    {
      "id": "318",
      "sequence": "G"
    },
    {
      "id": "243",
      "sequence": "A"
    },
    {
      "id": "50",
      "sequence": "CAGAGTAGAGTGGAGTGAGGACGACTGGATGGTAATTGAAAGGAATGGAATGGAACGGAGTTGAATGGAATGGAGAGGAATGCAATGGAATGGAGTGGAATGGAATGGAGTGGAGTGGAGTGGAGTTGAATAGAATGTACTGGAATGGCATGGAATGGAATGGAATGGAATGGAGTGGAGTGGAATGGAGTGGAGGGGAGACAAACGGAATGGAATGGAATGGAGGGGAGGGGAGTGAAGTGGAATGTAAACCAGTGG"
    },
    {
      "id": "194",
      "sequence": "GGA"
    },
    {
      "id": "167",
      "sequence": "TTGAATGGAATGGAATGGAAT"
    },
    {
      "id": "301",
      "sequence": "G"
    },
    {
      "id": "317",
      "sequence": "AATG"
    },
    {
      "id": "132",
      "sequence": "A"
    },
    {
      "id": "140",
      "sequence": "AA"
    },
    {
      "id": "202",
      "sequence": "CCA"
    },
    {
      "id": "248",
      "sequence": "GGAATGGAATGGAATG"
    },
    {
      "id": "169",
      "sequence": "C"
    },
    {
      "id": "42",
      "sequence": "T"
    },
    {
      "id": "180",
      "sequence": "A"
    },
    {
      "id": "255",
      "sequence": "G"
    },
    {
      "id": "160",
      "sequence": "A"
    },
    {
      "id": "87",
      "sequence": "GAGGGGAAAGAAATTGAGTGGAATTGAGTGG"
    },
    {
      "id": "289",
      "sequence": "TT"
    },
    {
      "id": "49",
      "sequence": "A"
    },
    {
      "id": "291",
      "sequence": "G"
    },
    {
      "id": "106",
      "sequence": "G"
    },
    {
      "id": "94",
      "sequence": "A"
    },
    {
      "id": "225",
      "sequence": "T"
    },
    {
      "id": "128",
      "sequence": "GGA"
    },
    {
      "id": "347",
      "sequence": "AT"
    },
    {
      "id": "259",
      "sequence": "A"
    },
    {
      "id": "350",
      "sequence": "G"
    },
    {
      "id": "379",
      "sequence": "G"
    },
    {
      "id": "375",
      "sequence": "GGAAT"
    },
    {
      "id": "21",
      "sequence": "GTAGAATGGAATGGAATGAAATGGAATGGATTGGAGTGCAGGGGAGCAGAATGCAATGGAAAGGAGTGAA"
    },
    {
      "id": "229",
      "sequence": "TGGA"
    },
    {
      "id": "38",
      "sequence": "G"
    },
    {
      "id": "163",
      "sequence": "TGGAAT"
    },
    {
      "id": "332",
      "sequence": "CA"
    },
    {
      "id": "131",
      "sequence": "GGA"
    },
    {
      "id": "102",
      "sequence": "G"
    },
    {
      "id": "192",
      "sequence": "CG"
    },
    {
      "id": "70",
      "sequence": "T"
    },
    {
      "id": "326",
      "sequence": "GAAT"
    },
    {
      "id": "221",
      "sequence": "G"
    },
    {
      "id": "373",
      "sequence": "C"
    },
    {
      "id": "53",
      "sequence": "ATGTAGTGGAGTGAAGTGGATTGGAATGGAATATAGTGGAATTGAATGGAATGGAGTGGAATGCAATTTACCGAAATGGAAAGGAACGGAATGGAGTAAAGTTGAGTGGAATGGAATTGAGTGGAGTGGTATGGAATGGAATGGAATGGAATGGA"
    },
    {
      "id": "362",
      "sequence": "C"
    },
    {
      "id": "47",
      "sequence": "TC"
    },
    {
      "id": "175",
      "sequence": "A"
    },
    {
      "id": "286",
      "sequence": "A"
    },
    {
      "id": "338",
      "sequence": "G"
    },
    {
      "id": "178",
      "sequence": "GGAATGGAAT"
    },
    {
      "id": "3",
      "sequence": "A"
    },
    {
      "id": "96",
      "sequence": "AA"
    },
    {
      "id": "306",
      "sequence": "G"
    },
    {
      "id": "149",
      "sequence": "ATGGAATGGAATGGAATGGA"
    },
    {
      "id": "155",
      "sequence": "G"
    },
    {
      "id": "181",
      "sequence": "G"
    },
    {
      "id": "65",
      "sequence": "G"
    },
    {
      "id": "293",
      "sequence": "GAA"
    },
    {
      "id": "298",
      "sequence": "C"
    }
  ]
}
        
        )";
        
        vg::VG graph;
        vg::io::json2graph(graph_json, &graph);
        
        Alignment aln;
        aln.set_sequence("GGAATGCAATGGAAAGAAATGGAATGGAATGGAATGAAAAGGAATGGAATGGAAAGAAGTGCAGTGGAGTGGAATGGAATTGAGTGAAATGGAATGGAAAGGAAATGGAATGGAGTGCAGTGGAGTGGAGTGGGGTCGAGTGGAATGGAATTGAACGGAATGGAATGGAATTTAATGGAATGGAATGGAAAGGATTGGAATGGAATGGAACAGAATTCTATGGAGTGGAATCGAATGGAATGGAAACGAAAGGATTGGAATGGAAAGGAAAGGAACGGATTTGCCTGGAATGGTTTGGAATGGAATGCAGTGGAACGCATTGGAGTGGAATGGAATGGAGTGGAATGGATTGGAGTGGAGTCTAATTCAATGGAGTGGAATGGAGTGGAATGGAATGGAATGGAATGGATTCCTGTGGAAAGAATATTAATGGAATGGATTGGAGTGGAATGGAGAAGAATGGCGTGGAGTGAAATGGAATGGAGAGCAATGGAATTGAGTGGAATGGAGTTAAGTGCTGTGGAATAGATTAGAGTGCAATGGAGCTTAGGGGAGTGCAGTGGAATGGAGTGGAATAGATTTGAATGTATTGGAATGAAATGGAATAGAAAGAAATGGAATGGAATGGAAAGAAATGGAATGGAATAGAATGGAATGCTATTGAGTGGAGTGGAGTTGGTTCGAGTGGATGGGGATGAAATGGAATGAAATGGATAGTAATAGAATAGAATAAAATGGAAATGAGGGGAGTGGAGTGAAATGGAAGGCAGTCGATTGGAGTGCAGTAGAATGGAATGGAATGGAATGACTTGGTGTGGAATGAAATGGAGTGGAATTGATGGAGTGGAATGGAATGGATTGGAATGGACTGGAATCGATTGAAGTGGAATGGAATAGAGTGGAATGTATTGGAACGGAGTGTATTGGAATGGAACGCAATGGAAAATGATGGAATGAAATAGAAAGGAATGGAACTAAGTTTAGTGCATTGGAATGGAGTTGAGTGGATTGGAAAGGAATAAAAGGGAATGGAATGCAATGGAGTGGAGTGGAGTGGAGCGGAAGGGAATGGAACGGAATGGAATGGAGTGGAATGGAATGGAGTGGAATGGAATGGCATGGAATGGATTGGAATTGAATGGAGTGGAATGGAATTGACTGGAATGGAATGGAGTGGAAAGCAATGGAGTGGAGTGGAACGGAGGAGGGGTCGAGTAGATGGGAATGGAATGGAATGGAGTGGAGTGGAATAGAGTGGAATGGAGAGGAGTGGTGTGGAGTGTAATGGATTGAGTAGAGAGAAATAGAATGGAATGGAATGGAATGGAATGCAATGGAATTCAATTGAATTCAATATAATGAAATAGAATGGAGAGGATGGGAATTAACTAGAGTGGAATGGAGTGGAATGAGTGGAGTGGAATGGAATGGAATCGAATTAAGCGGGATGTAATGGAATAGAATGCATTGAAATGGAATGGATTGGACGGGACTGGAATGGAATTGAGAGGAGAAAAGCAGAATTGAATGGCATTGAATAGAGTGGAATGCAGTGCATTGGGGTGGAGTGGAATGGAACGGAATGGAGTGAAGTTGAAGGGAACGGAATGCAATGGAATGCAATGGAATGGAATGGAATGGAATGGAATGGAATCCAGTGGAGTGGAATGGAATGGAATGTAAAGGAATGGAATGGAATGGTGTGGAGTGGAATAGAATGGAAGGGAATGCAGTGGAACGGAATGGAATGCAATGGAATGGAATGGAGTGGGGTGGAGTGGAATGGAATTAAGTGGACTGGAATATAATGAAATGGAATGGAGTGGAGTCGAGTGGAGACTGGTCGAGTGGAATGGAATGGAATGGAGTGGAGTGTAAAGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATTCCATGGAATGGAATGGAATTCCATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATTCCATGGAATTCCATGGAATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATTGCATGGAATGGAATTCCATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCCATGGAATGGAATGGAATGGAATGGAATTCC");


        pos_t right_anchor {383, false, 0};
        
        TestMinimizerMapper::align_sequence_between(empty_pos_t(), right_anchor, 5000, 500, &graph, &aligner, aln);

        // We demand a positive-score alignment
        REQUIRE(aln.score() > 0);
        // We demand not having a very long softclip at the start
        REQUIRE(aln.path().mapping_size() > 0);
        auto& first_mapping = aln.path().mapping(0);
        REQUIRE(first_mapping.edit_size() > 0);
        auto& first_edit = first_mapping.edit(0);
        REQUIRE(first_edit.to_length() <= std::max(10, first_edit.from_length()));
}

TEST_CASE("MinimizerMapper can extract a strand-split dagified local graph without extraneous tips", "[giraffe][mapping]") {
    // Make the graph that was causing trouble (it's just a stick)
    std::string graph_json = R"(
        {
            "edge": [{"from": "60245280", "to": "60245281"},
                     {"from": "60245283", "to": "60245284"},
                     {"from": "60245282", "to": "60245283"},
                     {"from": "60245277", "to": "60245278"},
                     {"from": "60245279", "to": "60245280"},
                     {"from": "60245284", "to": "60245285"},
                     {"from": "60245281", "to": "60245282"},
                     {"from": "60245278", "to": "60245279"}],
            "node": [{"id": "60245280", "sequence": "GATTACAGATTACA"},
                     {"id": "60245283", "sequence": "GATTACAGATTACA"},
                     {"id": "60245282", "sequence": "GATTACAGATTACA"},
                     {"id": "60245277", "sequence": "GATTACAGATTACA"},
                     {"id": "60245285", "sequence": "GATTACAGATTACA"},
                     {"id": "60245279", "sequence": "GATTACAGATTACA"},
                     {"id": "60245284", "sequence": "GATTACAGATTACA"},
                     {"id": "60245281", "sequence": "GATTACAGATTACA"},
                     {"id": "60245278", "sequence": "GATTACAGATTACA"}]
        }
    )";
    vg::Graph graph_chunk;
    json2pb(graph_chunk, graph_json.c_str(), graph_json.size());
    vg::VG graph(graph_chunk);
    
    TestMinimizerMapper::with_dagified_local_graph(make_pos_t(60245283, false, 10), empty_pos_t(), 50, graph, [&](DeletableHandleGraph& dagified_graph, const std::function<std::pair<nid_t, bool>(const handle_t&)>& dagified_handle_to_base) {
        // The graph started as a stick
        // We strand-split it to two disconnected sticks, and then dagify from the one start node in the one orientation, so it should go back to being one stick, with 2 tips.
        auto tip_handles = handlegraph::algorithms::find_tips(&dagified_graph);
#ifdef debug
        for (auto& h : tip_handles) {
            // Dump all the tips for debugging
            auto original = dagified_handle_to_base(h);
            std::cerr << "Found tip handle " << dagified_graph.get_id(h) << (dagified_graph.get_is_reverse(h) ? "-" : "+") << " representing " << original.first << (original.second ? "-" : "+") << std::endl;
        }
#endif
        for (auto& h : tip_handles) {
            auto original = dagified_handle_to_base(h);
            if (!dagified_graph.get_is_reverse(h)) {
                // Any head must correspond to the anchoring node
                REQUIRE(original.first == 60245283);
                REQUIRE(original.second == false);
            }
        }
        // There should be that head and also some tail where we ran out of search bases.
        REQUIRE(tip_handles.size() == 2);
    });
}



}

}

