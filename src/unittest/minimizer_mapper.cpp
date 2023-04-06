/// \file minimizer_mapper.cpp
///  
/// unit tests for the minimizer mapper

#include <iostream>
#include "vg/io/json2pb.h"
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
        : MinimizerMapper(gbwt_graph, minimizer_index, distance_index, handle_graph){};
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

TEST_CASE("Mapping quality cap cannot be confused by excessively high base qualities", "[giraffe][mapping]") {
    string sequence;
    string quality;
    for (size_t i = 0; i < 250; i++) {
        sequence.push_back('G');
        quality.push_back((char)(unsigned char)0xFF);
    }
    
    vector<TestMinimizerMapper::Minimizer> minimizers;
    // They are all going to be explored
    vector<size_t> minimizers_explored;
    
    TestMinimizerMapper::Minimizer minimizer_template;
    minimizer_template.value.is_reverse = false;
    
    minimizer_template.hits = 229;
    // We knowe the occurrences won't be used.
    minimizer_template.occs = nullptr;
    minimizer_template.score = 1;
    
    for (size_t try_number = 0; try_number < 1000; try_number++) {
    
        for (size_t i = 0; i < 10; i++) {
            size_t core_width = rand() % std::min(sequence.size()/2, (size_t)32);
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
        
        TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, &graph, &aligner, aln);
        
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
        
        TestMinimizerMapper::align_sequence_between(left_anchor, right_anchor, 100, &graph, &aligner, aln);
        
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

