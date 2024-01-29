/// \file funnel.cpp
///  
/// Unit tests for the Funnel class.
///

#include <iostream>
#include <vector>

#include "../funnel.hpp"

#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;
    
TEST_CASE("Funnel tracks tags correctly through merge_group", "[funnel]") {
    
    Funnel funnel;

    funnel.stage("seed");
    funnel.introduce(3);

    funnel.tag(1, Funnel::State::CORRECT, 0, 10);
    funnel.tag(2, Funnel::State::PLACED, 100, 110);

    std::vector<size_t> seeds_to_merge {0, 1, 2};

    funnel.stage("tree");
    funnel.merge_group(seeds_to_merge.begin(), seeds_to_merge.end());

    funnel.stage("fragment");
    funnel.introduce();
    funnel.also_merge_group(2, seeds_to_merge.begin(), seeds_to_merge.end());
    funnel.also_relevant(1, 0);

    std::vector<size_t> fragments_to_merge {0};

    funnel.stage("chain");
    funnel.merge_group(fragments_to_merge.begin(), fragments_to_merge.end());

    REQUIRE(funnel.last_correct_stage() == "chain"); 

}
}
}
        
