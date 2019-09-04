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

// We define a child class to expose all the protected stuff for testing
class TestMinimizerMapper : public MinimizerMapper {
public:
    using MinimizerMapper::MinimizerMapper;
    using MinimizerMapper::score_extension_group;
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
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 11;
        to_score.back().read_interval.second = 20;
        to_score.back().score = 9;
        
        SECTION("Score of two 9-base extensions separated by a 1-base gap is 12") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 12);
        }
        
    }
}


}

}
