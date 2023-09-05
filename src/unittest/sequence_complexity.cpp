///
/// \file sequence_complexity.cpp
///  
/// Unit tests for SeqComplexity
///

#include "catch.hpp"
#include "../sequence_complexity.hpp"
#include <iostream>

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE( "Sequence complexity can be identified", "[complexity]" ) {
    
    string seq = "ACGTACGTACGTACGT";
    SeqComplexity<16> complexity(seq);
    for (int order = 1; order <= 16; ++order) {
        if (order % 4 == 0 && order < seq.size()) {
            REQUIRE(complexity.p_value(order) < .01);
        }
        else {
            REQUIRE(complexity.p_value(order) > .95);
        }
    }
}

}
}
