/// \file minimizer_mapper.cpp
///  
/// unit tests for the minimizer mapper

#include <iostream>
#include "json2pb.h"
#include <vg/vg.pb.h>
#include "../minimizer_mapper.hpp"
#include "../build_index.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

// We define a child class to expose protected stuff for testing
class TestMinimizerMapper : public MinimizerMapper {
public:
    using MinimizerMapper::MinimizerMapper;
    using MinimizerMapper::score_extension_group;
    using MinimizerMapper::window_breaking_quality;
    using MinimizerMapper::Minimizer;
};

TEST_CASE("MinimizerMapper::score_extension_group works", "[giraffe][mapping]") {

    // Define an Alignment to score extensions against.
    Alignment aln;
    aln.set_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    
    // Define a bunch of fake GaplessExtensions to chain up and score
    vector<GaplessExtension> to_score;
    
    SECTION("Score of no gapless extensions is 0") {
        REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 0);
    }
    
    SECTION("One-base extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 1;
        to_score.back().read_interval.second = 2;
        to_score.back().score = 1;
        
        SECTION("Score of a 1-base gapless extension is 1") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 1);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 2;
        to_score.back().read_interval.second = 3;
        to_score.back().score = 1;
        
        SECTION("Score of two adjacent 1-base gapless extensions is 2") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 2);
        }
        
    }
    
    SECTION("Longer extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 1;
        to_score.back().read_interval.second = 10;
        to_score.back().score = 9;
        
        SECTION("Score of one 9-base extension is 9") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 9);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 10;
        to_score.back().read_interval.second = 19;
        to_score.back().score = 9;
        
        SECTION("Score of two 9-base extensions abutting is 18") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 18);
        }
        
        to_score.back().read_interval.first++;
        to_score.back().read_interval.second++;
        
        SECTION("Score of two 9-base extensions separated by a 1-base gap is 12") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 12);
        }
        
        to_score.back().read_interval.first++;
        to_score.back().read_interval.second++;
        
        SECTION("Score of two 9-base extensions separated by a 2-base gap is 11") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 11);
        }
        
        to_score.back().read_interval.first -= 3;
        to_score.back().read_interval.second -= 3;
        
        SECTION("Score of two 9-base extensions with a 1-base ovealap is 12") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 11);
        }
    }
    
    SECTION("Many possibly overlapping extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 0;
        to_score.back().read_interval.second = 1;
        to_score.back().score = 1;
    
        for (size_t i = 0; i < 35; i++) {
        
            to_score.emplace_back();
            to_score.back().read_interval.first = i + 1;
            to_score.back().read_interval.second = i + 1 + 30;
            to_score.back().score = 30;
        
        }
    
        
        
        SECTION("Score of one 1-base extensions and 2 30-base extensions is 61") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 61);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 28;
        to_score.back().read_interval.second = 28 + 45;
        to_score.back().score = 45;
        
        // Sort by read interval as is required
        std::sort(to_score.begin(), to_score.end(), [](const GaplessExtension& a, const GaplessExtension& b) {
            return (a.read_interval.first < b.read_interval.first) ||
                (a.read_interval.first == b.read_interval.first && a.read_interval.second < b.read_interval.second);
        });
        
        SECTION("Score of one 1-base extension, a 30-base extension, a backtrack of 4, and a 45-base extension is correct") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == (31 + 45 - 6 - 3));
        }
        
        
        for (size_t i = 3; i < 29; i++) {
             to_score.emplace_back();
            to_score.back().read_interval.first = i;
            to_score.back().read_interval.second = i + 15;
            to_score.back().score = 15;
        }
        
        // Sort by read interval as is required
        std::sort(to_score.begin(), to_score.end(), [](const GaplessExtension& a, const GaplessExtension& b) {
            return (a.read_interval.first < b.read_interval.first) ||
                (a.read_interval.first == b.read_interval.first && a.read_interval.second < b.read_interval.second);
        });
        
        SECTION("Score of one 1-base extension, a 30-base extension, a backtrack of 4, and a 45-base extension is not distracted") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == (31 + 45 - 6 - 3));
        }
    
    }
}


}

}
