/**
 * unittest/distributions.cpp: test cases for distributions.hpp
 */

#include "catch.hpp"
#include "distributions.hpp"

namespace vg {
namespace unittest {

TEST_CASE( "Factorials are computed", "[distributions][factorial]" ) {
    REQUIRE(factorial_ln(0) == Approx(log(1)).epsilon(1E-10));
    REQUIRE(factorial_ln(1) == Approx(log(1)).epsilon(1E-10));
    REQUIRE(factorial_ln(2) == Approx(log(2)).epsilon(1E-10));
    REQUIRE(factorial_ln(3) == Approx(log(6)).epsilon(1E-10));
    REQUIRE(factorial_ln(10) == Approx(log(3628800)).epsilon(1E-10));
}

}
}
