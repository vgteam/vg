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
    using MinimizerMapper::to_anchor;
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

TEST_CASE("MinimizerMapper can map with an initial deletion", "[giraffe][mapping][right_tail]") {

        Aligner aligner;
        
        string graph_json = R"({
            "edge": [
                {"from": "1", "to": "2"},
                {"from": "1", "to": "3"}
            ],
            "node": [
                {"id": "1", "sequence": "T"},
                {"id": "2", "sequence": "GATTACA"},
                {"id": "3", "sequence": "CATTAG"}
            ]
        })";
        
        // TODO: Write a json_to_handle_graph
        vg::Graph proto_graph;
        json2pb(proto_graph, graph_json.c_str(), graph_json.size());
        auto graph = vg::VG(proto_graph);
        
        Alignment aln;
        aln.set_sequence("CATTAG");
        
        pos_t left_anchor {1, false, 0}; // This includes the base on node 1
        pos_t right_anchor = empty_pos_t();
        
        TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, 20, &graph, &aligner, aln);

        // Make sure we get the right alignment. We should have a 1bp deletion and then the matching node.
        REQUIRE(aln.path().mapping_size() == 2);
        REQUIRE(aln.path().mapping(0).position().node_id() == 1);
        REQUIRE(aln.path().mapping(0).position().is_reverse() == false);
        REQUIRE(aln.path().mapping(0).position().offset() == 0);
        REQUIRE(aln.path().mapping(0).edit_size() == 1);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 1);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
        REQUIRE(aln.path().mapping(1).position().node_id() == 3);
        REQUIRE(aln.path().mapping(1).position().is_reverse() == false);
        REQUIRE(aln.path().mapping(1).position().offset() == 0);
        REQUIRE(aln.path().mapping(1).edit_size() == 1);
        REQUIRE(aln.path().mapping(1).edit(0).from_length() == 6);
        REQUIRE(aln.path().mapping(1).edit(0).to_length() == 6);
        REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
}

TEST_CASE("MinimizerMapper can map with an initial deletion on a multi-base node", "[giraffe][mapping][right_tail]") {

        Aligner aligner;
        
        string graph_json = R"({
            "edge": [
                {"from": "1", "to": "2"},
                {"from": "1", "to": "3"}
            ],
            "node": [
                {"id": "1", "sequence": "TATA"},
                {"id": "2", "sequence": "GATTACA"},
                {"id": "3", "sequence": "CATTAG"}
            ]
        })";
        
        // TODO: Write a json_to_handle_graph
        vg::Graph proto_graph;
        json2pb(proto_graph, graph_json.c_str(), graph_json.size());
        auto graph = vg::VG(proto_graph);
        
        Alignment aln;
        aln.set_sequence("CATTAG");
        
        pos_t left_anchor {1, false, 3}; // This includes the last base on node 1
        pos_t right_anchor = empty_pos_t();
        
        TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, 20, &graph, &aligner, aln);

        // Make sure we get the right alignment. We should have a 1bp deletion and then the matching node.
        REQUIRE(aln.path().mapping_size() == 2);
        REQUIRE(aln.path().mapping(0).position().node_id() == 1);
        REQUIRE(aln.path().mapping(0).position().is_reverse() == false);
        REQUIRE(aln.path().mapping(0).position().offset() == 3);
        REQUIRE(aln.path().mapping(0).edit_size() == 1);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 1);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
        REQUIRE(aln.path().mapping(1).position().node_id() == 3);
        REQUIRE(aln.path().mapping(1).position().is_reverse() == false);
        REQUIRE(aln.path().mapping(1).position().offset() == 0);
        REQUIRE(aln.path().mapping(1).edit_size() == 1);
        REQUIRE(aln.path().mapping(1).edit(0).from_length() == 6);
        REQUIRE(aln.path().mapping(1).edit(0).to_length() == 6);
        REQUIRE(aln.path().mapping(1).edit(0).sequence().empty());
}

TEST_CASE("MinimizerMapper can map right off the past-the-end base", "[giraffe][mapping][right_tail]") {

        Aligner aligner;
        
        string graph_json = R"({
            "edge": [
                {"from": "1", "to": "2"},
                {"from": "1", "to": "3"}
            ],
            "node": [
                {"id": "1", "sequence": "T"},
                {"id": "2", "sequence": "GATTACA"},
                {"id": "3", "sequence": "CATTAG"}
            ]
        })";
        
        // TODO: Write a json_to_handle_graph
        vg::Graph proto_graph;
        json2pb(proto_graph, graph_json.c_str(), graph_json.size());
        auto graph = vg::VG(proto_graph);
        
        Alignment aln;
        aln.set_sequence("CATTAG");
        
        pos_t left_anchor {1, false, 1}; // This is the past-end position
        pos_t right_anchor = empty_pos_t();
        
        TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, 20, &graph, &aligner, aln);

        // Make sure we get the right alignment. We should pick the matching node and use it. 
        REQUIRE(aln.path().mapping_size() == 1);
        REQUIRE(aln.path().mapping(0).position().node_id() == 3);
        REQUIRE(aln.path().mapping(0).position().is_reverse() == false);
        REQUIRE(aln.path().mapping(0).position().offset() == 0);
        REQUIRE(aln.path().mapping(0).edit_size() == 1);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 6);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 6);
        REQUIRE(aln.path().mapping(0).edit(0).sequence().empty());
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

TEST_CASE("MinimizerMapper can make correct anchors from minimizers and their zip codes", "[giraffe][mapping]") {
    Alignment aln;
    aln.set_sequence("AAAAAAAAAA"); // 10 bp

    // I only need a linear graph to test translation (ignoring running off the ends).
    // TODO: Test trimmign back from node ends.
    VG graph;

    Node* n1 = graph.create_node("AAAAAAAAAA");

    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);

    for (bool graph_reverse_strand : {false, true}) {
        // Try the read running both forward and backward along the graph.
        
        for (bool anchor_a_reverse : {false, true}) {
            for (bool anchor_b_reverse : {false, true}) {
                // Try all combinations of first and second hit minimizer
                // orientations relative to the read.
        
                // These are graph positions for each minimizer hit. They are first read
                // bases for forward-read-strand minimizers, and last read bases for
                // reverse-read-strand minimizers, and they always point in the read's
                // forward direction.
                std::vector<pos_t> graph_positions;

                // These are read positions for each minimizer hit, in the form of an
                // anchoring base on the read's forward strand, and an orientation from
                // that anchoring base for the minimizer sequence's orientation/where the
                // rest of the minimizer sequence falls in the read. 
                //
                // False is read forward (minimizer occurrence is here and to the right),
                // true is read reverse (minimizer occurrence is here and to the left,
                // minimal sequence is from the read's reverse strand).
                std::vector<std::pair<size_t, bool>> read_positions;

                // These are the minimizer lengths
                std::vector<size_t> lengths;

                if (anchor_a_reverse) {
                    // Have a 3bp hit at the start of the read and graph. It is anchored at its
                    // final location in the read.
                    graph_positions.emplace_back(1, graph_reverse_strand, 2);
                    read_positions.emplace_back(2, true);
                    lengths.emplace_back(3);
                } else {
                    // Have a 3bp hit at the start of the read and graph. It is anchored at its
                    // start location in the read.
                    graph_positions.emplace_back(1, graph_reverse_strand, 0);
                    read_positions.emplace_back(0, false);
                    lengths.emplace_back(3);
                }
                
                if (anchor_b_reverse) {
                    // Have another 3bp hit at the end, with the graph and read still going in
                    // the same direction, but with the minimizer on the other strand of the
                    // read.
                    //
                    // It is anchored at its final location in the read, but the position is
                    // still on the forward strand of the graph, since the read is still going
                    // forward along the graph node. 
                    graph_positions.emplace_back(1, graph_reverse_strand, 9);
                    read_positions.emplace_back(9, true);
                    lengths.emplace_back(3);
                } else {
                    // Have another 3bp hit at the end, anchored at its start location in the read.
                    graph_positions.emplace_back(1, graph_reverse_strand, 7);
                    read_positions.emplace_back(7, false);
                    lengths.emplace_back(3);
                }

                // Add a middle anchor overlapping the left one
                graph_positions.emplace_back(1, graph_reverse_strand, 1);
                read_positions.emplace_back(1, false);
                lengths.emplace_back(3);

                // Add a middle anchor actually in the middle, abutting the left one, and shorter
                graph_positions.emplace_back(1, graph_reverse_strand, 3);
                read_positions.emplace_back(3, false);
                lengths.emplace_back(2);


                vector<MinimizerMapper::Minimizer> minimizers;
                vector<SnarlDistanceIndexClusterer::Seed> seeds;
                for (size_t i = 0; i < read_positions.size(); i++) {
                    // Make a minimizer
                    minimizers.emplace_back();
                    minimizers.back().length = lengths.at(i);
                    minimizers.back().value.offset = read_positions.at(i).first;
                    minimizers.back().value.is_reverse = read_positions.at(i).second;

                    // Make a zipcode for its graph position
                    ZipCode zipcode;
                    zipcode.fill_in_zipcode(distance_index, graph_positions.at(i));

                    // Make a seed attaching that graph position to its minimizer.
                    seeds.push_back({ graph_positions.at(i), i, zipcode});
                }
                VectorView<MinimizerMapper::Minimizer> minimizer_vector (minimizers);

                // Make and check the zip code tree
                ZipCodeForest zip_forest;
                zip_forest.fill_in_forest(seeds, minimizer_vector, distance_index, 10);
                REQUIRE(zip_forest.trees.size() == 1);

                // Make an aligner for scoring
                Aligner aligner;

                // Make the anchors
                std::vector<algorithms::Anchor> anchors;
                for (size_t i = 0; i < seeds.size(); i++) {
#ifdef debug
                    std::cerr << "Anchor " << i << ":" << std::endl;
#endif
                    anchors.push_back(TestMinimizerMapper::to_anchor(aln, minimizers, seeds, i, graph, &aligner));

                    // Make sure the anchor is right.
                    // It needs to start at the right place in the read.
                    REQUIRE(anchors.back().read_start() == minimizers.at(seeds.at(i).source).forward_offset());
                    // Sinve the minimizers are all within single nodes here, the anchor should be as long as the minimizer.
                    REQUIRE(anchors.back().length() == minimizers.at(seeds.at(i).source).length);
                }

                // For each form anchor and to anchor, remember the read and graph distances.
                std::unordered_map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> all_transitions;

                // Set up to get all the transitions between anchors in the zip code tree
                auto transition_iterator = algorithms::zip_tree_transition_iterator(seeds, zip_forest.trees.at(0), std::numeric_limits<size_t>::max());
                // And get them
                transition_iterator(anchors, distance_index, graph, std::numeric_limits<size_t>::max(), [&](size_t from_anchor, size_t to_anchor, size_t read_distance, size_t graph_distance) {
                    // And for each of them, remember them
#ifdef debug
                    std::cerr << "From anchor " << from_anchor << " to anchor " << to_anchor << " we cross " << read_distance << " bp of read and " << graph_distance << " bp of graph" << std::endl;
#endif
                    all_transitions.emplace(std::make_pair(from_anchor, to_anchor), std::make_pair(read_distance, graph_distance));
                });

                // Make sure we got the right transitions for these anchors
                // AAAAAAAAAA
                // XXX----YYY
                //   01234
                REQUIRE(all_transitions.at(std::make_pair(0, 1)).first == 4);
                REQUIRE(all_transitions.at(std::make_pair(0, 1)).second == 4);

                // AAAAAAAAAA
                // -XXX---YYY
                //    0123
                REQUIRE(all_transitions.at(std::make_pair(2, 1)).first == 3);
                REQUIRE(all_transitions.at(std::make_pair(2, 1)).second == 3);

                // AAAAAAAAAA
                // ---XX--YYY
                //     012
                REQUIRE(all_transitions.at(std::make_pair(3, 1)).first == 2);
                REQUIRE(all_transitions.at(std::make_pair(3, 1)).second == 2);

                // AAAAAAAAAA
                // XXXYY-----
                //   0
                REQUIRE(all_transitions.at(std::make_pair(0, 3)).first == 0);
                REQUIRE(all_transitions.at(std::make_pair(0, 3)).second == 0);

                // We shouldn't see any extra transitions, like between overlapping anchors.
                REQUIRE(all_transitions.size() == 4);
            }
        }
    }
}



}

}

