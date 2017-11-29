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

TEST_CASE( "Ambiguous multinomial works", "[distributions][multinomial]" ) {

    vector<double> probs{0.5, 0.25, 0.25};
    unordered_map<vector<bool>, int> obs{{{true, false, false}, 1}, {{false, true, true}, 2}};

    REQUIRE(multinomial_censored_sampling_prob_ln(probs, obs) == Approx(log(0.125)).epsilon(1E-10));

}

}
}
