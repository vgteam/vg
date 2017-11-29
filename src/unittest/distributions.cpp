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

    // Define some category probabilities
    vector<double> probs{0.5, 0.25, 0.25};
    
    // And an ambiguous observation: 1 in category 1, 2 split between 2 and 3.
    unordered_map<vector<bool>, int> obs{{{true, false, false}, 1}, {{false, true, true}, 2}};
    
    // Define the 3 real possible cases matching the observations
    vector<int> case1{1, 2, 0};
    vector<int> case2{1, 1, 1};
    vector<int> case3{1, 0, 2};
    
    // Sum up their probabilities
    auto case_total = logprob_add(logprob_add(multinomial_sampling_prob_ln(probs, case1),
        multinomial_sampling_prob_ln(probs, case2)), multinomial_sampling_prob_ln(probs, case3));

    // Make sure we get the same answer for the probability of matching those constraints.
    REQUIRE(multinomial_censored_sampling_prob_ln(probs, obs) == Approx(case_total).epsilon(1E-10));

}

}
}
