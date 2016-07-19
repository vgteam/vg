#ifndef VG_DISTRIBUTIONS_H
#define VG_DISTRIBUTIONS_H

// distributions.hpp: contains functions for probability distributions used in genotyping.
// Also includes a fancy cacheing factorial system from Freebayes


#include <map>
#include <cmath>

#include "utility.hpp"

namespace vg {

using namespace std;

// We use this slightly nonstandard type for our math. We wrap it so it's easy
// to change later.
using real_t = long double;

/**
 * Calculate the natural log of the gamma function of the given argument.
 */
inline real_t gamma_ln(real_t x) {

    real_t cofactors[] = {76.18009173, 
        -86.50532033,
        24.01409822,
        -1.231739516,
        0.120858003E-2,
        -0.536382E-5};    

    real_t x1 = x - 1.0;
    real_t tmp = x1 + 5.5;
    tmp -= (x1 + 0.5) * log(tmp);
    real_t ser = 1.0;
    for (int j=0; j<=5; j++) {
        x1 += 1.0;
        ser += cofactors[j]/x1;
    }
    real_t y =  (-1.0 * tmp + log(2.50662827465 * ser));

    return y;
}

/**
 * Calculate the natural log of the factorial of the given integer. TODO:
 * replace with a cache or giant lookup table from Freebayes.
 */
inline real_t factorial_ln(int n) {
    if (n < 0) {
        return (long double)-1.0;
    }
    else if (n == 0) {
        return (long double)0.0;
    }
    else {
        return gamma_ln(n + 1.0);
    }
}

/**
 * Raise a log probability to a power
 */
inline real_t pow_ln(real_t m, int n) {
    return m * n;
}

/**
 * Compute the number of ways to select k items from a collection of n
 * distinguishable items, ignoring order. Returns the natural log of the
 * (integer) result.
 */
inline real_t choose_ln(int n, int k) {
    return factorial_ln(n) - (factorial_ln(k) + factorial_ln(n - k));
}

/**
 * Get the probability for sampling the counts in obs from a set of categories
 * weighted by the probabilities in probs. Works for both double and real_t
 * probabilities. Also works for binomials.
 */
template <typename ProbIn>
real_t multinomial_sampling_prob_ln(const vector<ProbIn>& probs, const vector<int>& obs) {
    vector<real_t> factorials;
    vector<real_t> probsPowObs;
    factorials.resize(obs.size());
    transform(obs.begin(), obs.end(), factorials.begin(), factorial_ln);
    typename vector<ProbIn>::const_iterator p = probs.begin();
    vector<int>::const_iterator o = obs.begin();
    for (; p != probs.end() && o != obs.end(); ++p, ++o) {
        probsPowObs.push_back(pow_ln(log(*p), *o));
    }
    // Use the collection sum defined in utility.hpp
    return factorial_ln(sum(obs)) - sum(factorials) + sum(probsPowObs);
}

/**
 * Compute the probability of having the given number of successes or fewer in
 * the given number of trials, with the given success probability. Returns the
 * resulting log probability.
 */
template <typename ProbIn>
real_t binomial_cmf_ln(ProbIn success_logprob, size_t trials, size_t successes) {
    // Compute log probabilities for all cases
    vector<real_t> case_logprobs;
    
    if(successes > trials) {
        return prob_to_logprob(0);
    }
    
    for(size_t considered_successes = 0; considered_successes <= successes; considered_successes++) {
        // For every number of successes up to this one, add in the probability.
        case_logprobs.push_back(choose_ln(trials, considered_successes) +
            success_logprob * considered_successes +
            logprob_invert(success_logprob) * (trials - considered_successes));
    }
    
    // Sum up all those per-case probabilities
    return logprob_sum(case_logprobs);
}



}

#endif
