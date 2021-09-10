#ifndef VG_STATISTICS_HPP_INCLUDED
#define VG_STATISTICS_HPP_INCLUDED

/**
 * \file statistics.hpp
 *
 * Defines a range of statistical functions
 *
 */

#include <cmath>
#include <algorithm>

#include "utility.hpp"

namespace vg {

using namespace std;
 
double median(std::vector<int> &v);
// Online mean-variance computation with wellfords algorithm (pass 0's to 1st 3 params to start)
void wellford_update(size_t& count, double& mean, double& M2, double new_val);
pair<double, double> wellford_mean_var(size_t count, double mean, double M2, bool sample_variance = false);


template<typename T>
double stdev(const T& v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    return std::sqrt(sq_sum / v.size());
}


/// The standard normal cumulative distribution function
double Phi(double x);

/// Inverse CDF of a standard normal distribution. Must have 0 < quantile < 1.
double Phi_inv(double quantile);

/*
 * Return the log of the sum of two log-transformed values without taking them
 * out of log space.
 */
inline double add_log(double log_x, double log_y) {
    return log_x > log_y ? log_x + log(1.0 + exp(log_y - log_x)) : log_y + log(1.0 + exp(log_x - log_y));
}

/*
 * Return the log of the difference of two log-transformed values without taking
 * them out of log space.
 */
inline double subtract_log(double log_x, double log_y) {
    return log_x + log(1.0 - exp(log_y - log_x));
}

/**
 * Convert a number ln to the same number log 10.
 */
inline double ln_to_log10(double ln) {
    return ln / log(10);
}

/**
 * Convert a number log 10 to the same number ln.
 */
inline double log10_to_ln(double l10) {
    return l10 * log(10);
}

/**
 * Given the log10 of a value, retunr the log10 of (that value plus one).
 */
inline double log10_add_one(double x) {
    return log10(pow(10, x) + 1);
}

/**
 * Return the log of the sum of two log10-transformed values without taking them
 * out of log space.
 */
inline double add_log10(double i, double j) {
    if (i < j) {
        return log10_add_one(j - i) + ((j - i < 10) ? i : j);
    }
    return log10_add_one(i - j) + ((i - j <= 10) ? j : i);
}

/**
 * Assume that we have n independent random events that occur with probability
 * p each (p is interpreted as a real number between 0 at 0 and 1 at its
 * maximum value). Return an approximate probability for at least one event
 * occurring as a phred score.
 *
 * n must be <= MAX_AT_LEAST_ONE_EVENTS.
 */
double phred_for_at_least_one(size_t p, size_t n);

/**
 * Assume that we have n independent random events that occur with probability
 * p each (p is interpreted as a real number between 0 at 0 and 1 at its
 * maximum value). Return an approximate probability for at least one event
 * occurring as a raw probability.
 *
 * n must be <= MAX_AT_LEAST_ONE_EVENTS.
 */
double prob_for_at_least_one(size_t p, size_t n);

/// How many events should we allow in phred_for_at_least_one and
/// prob_for_at_least_one? Should be >= our longest supported minimizer index.
constexpr static size_t MAX_AT_LEAST_ONE_EVENTS = 32;
/// Use this many bits for approximate probabilities.
constexpr static size_t AT_LEAST_ONE_PRECISION = 8; 

// normal pdf, from http://stackoverflow.com/a/10848293/238609
template <typename T>
T normal_pdf(T x, T m = 0.0, T s = 1.0)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T z = (x - m) / s;
    
    return inv_sqrt_2pi / s * std::exp(T(-0.5) * z * z);
}


/// Convert a probability to a natural log probability.
inline double prob_to_logprob(double prob) {
    return log(prob);
}
/// Convert natural log probability to a probability
inline double logprob_to_prob(double logprob) {
    return exp(logprob);
}
/// Add two probabilities (expressed as logprobs) together and return the result
/// as a logprob.
inline double logprob_add(double logprob1, double logprob2) {
    // Pull out the larger one to avoid underflows
    double pulled_out = max(logprob1, logprob2);
    return pulled_out + prob_to_logprob(logprob_to_prob(logprob1 - pulled_out) + logprob_to_prob(logprob2 - pulled_out));
}
/// Invert a logprob, and get the probability of its opposite.
inline double logprob_invert(double logprob) {
    return prob_to_logprob(1.0 - logprob_to_prob(logprob));
}

/// Convert 8-bit Phred quality score to probability of wrongness, using a lookup table.
double phred_to_prob(uint8_t phred);

/// Convert floating point Phred quality score to probability of wrongness.
inline double phred_to_prob(double phred) {
    return pow(10, -phred / 10);
}

/// Convert probability of wrongness to integer Phred quality score.
inline double prob_to_phred(double prob) {
    return -10.0 * log10(prob);
}

/// Convert a Phred quality score directly to a natural log probability of wrongness.
inline double phred_to_logprob(int phred) {
    return (-((double)phred) / 10) / log10(exp(1.0));
}

/// Convert a natural log probability of wrongness directly to a Phred quality score.
inline double logprob_to_phred(double logprob ) {
    return -10.0 * logprob * log10(exp(1.0));
}

/// Take the geometric mean of two logprobs
inline double logprob_geometric_mean(double lnprob1, double lnprob2) {
    return log(sqrt(exp(lnprob1 + lnprob2)));
}

/// Take the geometric mean of two phred-encoded probabilities
inline double phred_geometric_mean(double phred1, double phred2) {
    return prob_to_phred(sqrt(phred_to_prob(phred1 + phred2)));
}

/// Add two probabilities (expressed as phred scores) together and return the result
/// as a phred score.
inline double phred_add(double phred1, double phred2) {
    return logprob_to_phred(logprob_add(phred_to_logprob(phred1), phred_to_logprob(phred2)));
}



/**
 * Compute the sum of the values in a collection, where the values are log
 * probabilities and the result is the log of the total probability. Items must
 * be convertible to/from doubles for math.
 */
template<typename Collection>
typename Collection::value_type logprob_sum(const Collection& collection) {
    
    // Set up an alias
    using Item = typename Collection::value_type;
    
    // Pull out the minimum value
    auto min_iterator = min_element(begin(collection), end(collection));
    
    if(min_iterator == end(collection)) {
        // Nothing there, p = 0
        return Item(prob_to_logprob(0));
    }
    
    auto check_iterator = begin(collection);
    ++check_iterator;
    if(check_iterator == end(collection)) {
        // We only have a single element anyway. We don't want to subtract it
        // out because we'll get 0s.
        return *min_iterator;
    }
    
    // Pull this much out of every logprob.
    Item pulled_out = *min_iterator;
    
    if(logprob_to_prob(pulled_out) == 0) {
        // Can't divide by 0!
        // TODO: fix this in selection
        pulled_out = prob_to_logprob(1);
    }
    
    Item total(0);
    for(auto& to_add : collection) {
        // Sum up all the scaled probabilities.
        total += logprob_to_prob(to_add - pulled_out);
    }
    
    // Re-log and re-scale
    return pulled_out + prob_to_logprob(total);
}

/**
 * Compute the sum of the values in a collection, represented by an iterator
 * range, where the values are Phred scores and the result is the Phred score
 * of the total probability. Items must be convertible to/from doubles for
 * math.
 */
template<typename Iterator>
typename std::iterator_traits<Iterator>::value_type phred_sum(const Iterator& begin_it, const Iterator& end_it) {
    
    // Set up an alias for the type we're operating on
    using Item = typename std::iterator_traits<Iterator>::value_type;
    
    // Pull out the minimum probability
    auto min_iterator = max_element(begin_it, end_it);
    
    if (min_iterator == end_it) {
        // Nothing there, p = 0
        return Item(logprob_to_phred(prob_to_logprob(0)));
    }
    
    auto check_iterator = begin_it;
    ++check_iterator;
    if (check_iterator == end_it) {
        // We only have a single element anyway. We don't want to subtract it
        // out because we'll get 0s.
        return *min_iterator;
    }
    
    // Pull this much out of every logprob.
    double pulled_out = phred_to_logprob(*min_iterator);
    
    if (logprob_to_prob(pulled_out) == 0) {
        // Can't divide by 0!
        // TODO: fix this in selection
        pulled_out = prob_to_logprob(1);
    }
    
    double total(0);
    for(auto to_add_it = begin_it; to_add_it != end_it; ++to_add_it) {
        // Sum up all the scaled probabilities.
        total += logprob_to_prob(phred_to_logprob(*to_add_it) - pulled_out);
    }
    
    // Re-log and re-scale
    return Item(logprob_to_phred(pulled_out + prob_to_logprob(total)));
}

/**
 * Compute the sum of the values in a collection, where the values are Phred
 * scores and the result is the Phred score of the total probability. Items
 * must be convertible to/from doubles for math.
 */
template<typename Collection>
typename Collection::value_type phred_sum(const Collection& collection) {
    return phred_sum(begin(collection), end(collection));
}


double slope(const std::vector<double>& x, const std::vector<double>& y);
double fit_zipf(const vector<double>& y);

/// Returns the MLE rate parameter for the distribution of (shape) iid exponential RVs
double fit_fixed_shape_max_exponential(const vector<double>& x, double shape, double tolerance = 1e-8);

/// Returns the MLE estimate for the number of iid exponential RVs the data are maxima of
double fit_fixed_rate_max_exponential(const vector<double>& x, double rate, double tolerance = 1e-8);

/// Returns the MLE rate and shape parameters of a max exponential
pair<double, double> fit_max_exponential(const vector<double>& x, double tolerance = 1e-8);

// TODO: I'm eliminating this algorithm because it is approx non-identifiable for large values of shape
///// Returns the MLE rate, shape, and location parameters  of an offset max exponential
//tuple<double, double, double> fit_offset_max_exponential(const vector<double>& x, double tolerance = 1e-8);

/// Return the CDF of a max exponential with the given parameters
inline double max_exponential_cdf(double x, double rate, double shape, double location = 0.0) {
    return x > location ? pow(1.0 - exp(-(x - location) * rate), shape) : 0.0;
}

/// The log likelihood of a max exponential with the given parameters on the given data
double max_exponential_log_likelihood(const vector<double>& x, double rate, double shape,
                                      double location = 0.0);

/// Returns an estimate of the rate and shape parameters of a Weibull distribution
pair<double, double> fit_weibull(const vector<double>& x);

/// Returns an estimate of the rate, shape, and location (minimum value) of a 3-parameter Weibull distribution
tuple<double, double, double> fit_offset_weibull(const vector<double>& x,
                                                 double tolerance = 1e-8);

/// Return the CDF of a max exponential with the given parameters
inline double weibull_cdf(double x, double scale, double shape, double location = 0.0) {
    return x > location ? 1.0 - exp(-pow((x - location) / scale, shape)) : 0.0;
}

/// Returns the log likelihood of some data generated by a Weibull distribution
double weibull_log_likelihood(const vector<double>& x, double scale, double shape,
                              double location = 0.0);

/// Returns a local maximum of a function within an interval
double golden_section_search(const function<double(double)>& f, double x_min, double x_max,
                             double tolerance = 1e-8);

/// A shitty set of linear algebra functions

vector<vector<double>> transpose(const vector<vector<double>>& A);

vector<vector<double>> matrix_multiply(const vector<vector<double>>& A,
                                       const vector<vector<double>>& B);

vector<double> matrix_multiply(const vector<vector<double>>& A,
                               const vector<double>& b);

vector<vector<double>> matrix_invert(const vector<vector<double>>& A);


/// Returns the coefficients of a regression (does not automatically compute constant)
vector<double> regress(const vector<vector<double>>& X, vector<double>& y);

/*
 *********************************
 * Code ported over from FreeBayes
 *********************************
 */


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
 * Compute the number of ways to select k_1, k_2, ... k_i items into i buckets
 * from a collection of n distinguishable items, ignoring order. All of the
 * items have to go into the buckets, so all k_i must sum to n. To compute
 * choose you have to call this function with a 2-element vector, to represent
 * the chosen and not-chosen buckets. Returns the natural log of the (integer)
 * result.
 *
 * TODO: Turns out we don't actually need this for the ambiguous multinomial
 * after all.
 */
inline real_t multinomial_choose_ln(int n, vector<int> k) {
    // We use the product-of-binomial-coefficients approach from
    // <https://en.wikipedia.org/wiki/Multinomial_theorem#Multinomial_coefficients>
    real_t product_of_binomials_ln = 0;
    
    // We sum up the bucket sizes as we go
    int bucket_sum = 0;
    
    for (auto& bucket_size : k) {
        // Increment the size of what we choose from
        bucket_sum += bucket_size;
        // Choose this many
        product_of_binomials_ln += choose_ln(bucket_sum, bucket_size);
    }
    
    // Make sure they actually gave us a proper decomposition of the items into
    // buckets.
    assert(bucket_sum == n);
    
    return product_of_binomials_ln;
}

/**
 * Compute the log probability of a Poisson-distributed process: observed events
 * in an interval where expected events happen on average.
 */
inline real_t poisson_prob_ln(int observed, real_t expected) {
    return log(expected) * (real_t) observed - expected - factorial_ln(observed);
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

/**
 * Get the log probability for sampling the given value from a geometric
 * distribution with the given success log probability. The geometric
 * distribution is the distribution of the number of trials, with a given
 * success probability, required to observe a single success.
 */
template <typename ProbIn>
real_t geometric_sampling_prob_ln(ProbIn success_logprob, size_t trials) {
    return logprob_invert(success_logprob) * (trials - 1) + success_logprob;
}

/**
 * Given a split of items across a certain number of categories, as ints between
 * the two given bidirectional iterators, advance to the next split and return
 * true. If there is no next split, leave the collection unchanged and return
 * false.
 */
template<typename Iter>
bool advance_split(Iter start, Iter end) {
    if (start == end) {
        // Base case: we hit the end. No more possible splits.
#ifdef debug
        cerr << "Not advancing empty split" << endl;
#endif
        return false;
    } else {
        
#ifdef debug
        cerr << "Trying to advance split with " << *start << " items in first category" << endl;
#endif
        
        // Try advancing what comes after us.
        auto next = start;
        ++next;
        if (advance_split(next, end)) {
            // It worked.
#ifdef debug
            cerr << "Advanced child split" << endl;
#endif
            return true;
        }
        
#ifdef debug
        cerr << "Could not advance child split" << endl;
#endif
        
        // If that didn't work, try moving an item from here to what comes after us.
        // We also need to reset what comes after us to its initial state of everything in the first category.
        // This is easy because we know everything in what comes after us has made its way to the end.
        if (*start != 0 && next != end) {
            // We have something to move
            
            // Do the reset so everything after us is in the first category
            // after us.
            auto next_to_last = end;
            --next_to_last;
            
#ifdef debug
            cerr << "Want to move " << *next_to_last << " items to next which has " << *next << " and also move one from start which has " << *start << endl;
#endif
            
            if (next_to_last != next) {
                (*next) += *next_to_last;
                *next_to_last = 0;
            }
            
            (*start)--;
            (*next)++;
            
#ifdef debug
            cerr << "Reset child split to have " << *next << " items vs. our " << *start << endl;
#endif
            
            return true;
        }
        
        // If that didn't work, we're out of stuff to do.
#ifdef debug
        cerr << "Could not advance or reset child split" << endl;
#endif
        
        
        return false;
    }
}


/**
 * Get the log probability for sampling any actual set of category counts that
 * is consistent with the constraints specified by obs, using the per-category
 * probabilities defined in probs.
 *
 * Obs maps from a vector of per-category flags (called a "class") to a number
 * of items that might be in any of the flagged categories.
 *
 * For example, if there are two equally likely categories, and one item flagged
 * as potentially from either category, the probability of sampling a set of
 * category counts consistent with that constraint is 1. If instead there are
 * three equally likely categories, and one item flagged as potentially from two
 * of the three but not the third, the probability of sampling a set of category
 * counts consistent with that constraint is 2/3.
 */
template<typename ProbIn>
real_t multinomial_censored_sampling_prob_ln(const vector<ProbIn>& probs, const unordered_map<vector<bool>, int>& obs) {
    // We fill this with logprobs for all the different cases and then sum them
    // up.
    vector<real_t> case_logprobs;
    
    // We have a state. We advance this state until we can't anymore.
    //
    // The state is, for each ambiguity class, a vector of length equal to
    // the number of set bits in the valence, and sum equal to the number of
    // reads int he category.
    //
    // We start with all the reads in the first spot in each class, and
    // advance/reset until we have iterated over all combinations of category
    // assignments for all classes.
    unordered_map<vector<bool>, vector<int>> splits_by_class;
    
    // Prepare the state
    for (auto& kv : obs) {
        // For each input class
        
        if (kv.second == 0) {
            // No reads are in this class, so we can skip it.
            continue;
        }
        
        // Work out if it actually matches any categories
        bool has_any_categories = false;
        for (const auto& bit : kv.first) {
            if (bit) {
                has_any_categories = true;
                break;
            }
        }
        if (!has_any_categories) {
            // There are reads and they match nothing.
            // So this case is impossible.
            return prob_to_logprob(0);
        }
        
        // Otherwise there are reads and they match something.
        
        // For each class, find the vector we will use to describe its read-to-
        // category assignments.
        auto& class_state = splits_by_class[kv.first];
        
        for (const auto& bit : kv.first) {
            // Allocate a spot for each set bit
            if (bit) {
                class_state.push_back(0);
            }
        }
        
        // Drop all the reads in the first category
        class_state.front() = kv.second;
    }
    
    if (splits_by_class.empty()) {
        // There are no classes with any reads.
        // P(nothing happened) = 1.
        return prob_to_logprob(1);
    }
    
    // Now we loop over all the combinations of class states using a stack thing.
    list<decltype(splits_by_class)::iterator> stack;
    
    // And maintain this vector of category counts for the state we are in. We
    // incrementally update it so we aren't always looping over all the classes
    // to rebuild it.
    vector<int> category_counts(probs.size());
    
    // We have a function to add in the contribution of a class's state
    auto add_class_state = [&](const pair<vector<bool>, vector<int>>& class_state) {
        auto count_it = class_state.second.begin();
        for (size_t i = 0; i < category_counts.size(); i++) {
            // For each category
            if (class_state.first.at(i)) {
                // If this ambiguity class touches it
                
                assert(count_it != class_state.second.end());
                
                // Add in the state's count
                category_counts.at(i) += *count_it;
                
                // And consume that state entry
                ++count_it;
            }
        }
        assert(count_it == class_state.second.end());
    };
    
    // And a function to back it out again
    auto remove_class_state = [&](const pair<vector<bool>, vector<int>>& class_state) {
        auto count_it = class_state.second.begin();
        for (size_t i = 0; i < category_counts.size(); i++) {
            // For each category
            if (class_state.first.at(i)) {
                // If this ambiguity class touches it
                
                assert(count_it != class_state.second.end());
                
                // Back out the state's count
                category_counts.at(i) -= *count_it;
                
                // And consume that state entry
                ++count_it;
            }
        }
        assert(count_it == class_state.second.end());
    };
    
    for (auto it = splits_by_class.begin(); it != splits_by_class.end(); ++it) {
        // Populate the stack with everything
        stack.push_back(it);
        
        // And make sure the category counts are up to date.
        add_class_state(*it);
    }
    
    while (!stack.empty()) {
        // Emit the current state
        
#ifdef debug
        cerr << "Category counts:" << endl;
        for (auto& count : category_counts) {
            cerr << count << endl;
        }
#endif
        
        auto case_logprob = multinomial_sampling_prob_ln(probs, category_counts);
        
#ifdef debug
        cerr << "Case probability: " << logprob_to_prob(case_logprob) << endl;
#endif
        
        // Put in the logprob for this case.
        case_logprobs.push_back(case_logprob);
        
        while (!stack.empty()) {
            // See if we can advance what's at the bottom of the stack
            // First clear it out of the category counts.
            remove_class_state(*stack.back());
            if (advance_split(stack.back()->second.begin(), stack.back()->second.end())) {
                // We advanced it successfully.
                
#ifdef debug
                cerr << "Advanced class at stack depth " << stack.size() - 1 << endl;
#endif
                
                // Put it back in the category counts
                add_class_state(*stack.back());
                
                // We finally found something to advance, so stop ascending the stack.
                break;
            } else {
                
#ifdef debug
                cerr << "Could not advanced class at stack depth " << stack.size() - 1 << endl;
#endif
                
                // Pop off the back of the stack.
                stack.pop_back();
                
                // Keep looking up
            }
        }
        
        if (!stack.empty()) {
            // We found *something* to advance and haven't finished.
            
            // Now fill in the whole stack again with the first split for every category.
            auto it = stack.back();
            ++it;
            
            while (it != splits_by_class.end()) {
                // Reset the split to all 0s except for the first entry.
                for (auto& entry : it->second) {
                    entry = 0;
                }
                it->second.front() = obs.at(it->first);
                
                // Populate the stack with the next class
                stack.push_back(it);
                
                // And make sure the category counts are up to date.
                add_class_state(*it);
                
#ifdef debug
                cerr << "Reset class at stack depth " << stack.size() - 1 << endl;
#endif
                
                // Look for the next class
                ++it;
            }
        }
        
        // Otherwise we have finished looping over everything and so we should leave the stack empty.
        
    }
    
    // Sum up all those per-case probabilities
    return logprob_sum(case_logprobs);
}


// These handy sampling distribution implementations (uniform_real_distribution
// and normal_distribution) matching the C++ <random> API are copied and
// adapted from code at https://stackoverflow.com/a/34962942
template<typename T = double>
class uniform_real_distribution {
public:
    typedef T result_type;
    
    uniform_real_distribution(T _a = 0.0, T _b = 1.0) : m_a(_a), m_b(_b) {
        // Nothing to do!
    }
    
    void reset() {
        // Also nothing to do!
    }
    
    template<class Generator>
    T operator()(Generator &_g) {
        double dScale = (m_b - m_a) / ((T)(_g.max() - _g.min()) + (T)1);
        return (_g() - _g.min()) * dScale  + m_a;
    }
    
    T a() const {
        return m_a;
    }
    
    T b() const {
        return m_b;
    }
    
protected:
    T m_a;
    T m_b;
};

template<typename T = double>
class normal_distribution {
public:
    typedef T result_type;
    
    normal_distribution(T _mean = 0.0, T _stddev = 1.0) : m_mean(_mean), m_stddev(_stddev) {
        // Nothing to do!
    }
    
    void reset() {
        m_distU1.reset();
    }
    
    template<class Generator>
    T operator()(Generator &_g) {
        // Use Box-Muller algorithm
        const double pi = 3.14159265358979323846264338327950288419716939937511;
        double u1 = m_distU1(_g);
        double u2 = m_distU1(_g);
        double r = sqrt(-2.0 * log(u1));
        return m_mean + m_stddev * r * sin(2.0 * pi * u2);
    }
    
    T mean() const {
        return m_mean;
    }
    T stddev() const {
        return m_stddev;
    }
    
protected:
    T m_mean;
    T m_stddev;
    vg::uniform_real_distribution<T> m_distU1;
};

template<typename T = double>
class truncated_normal_distribution {
public:
    typedef T result_type;
    
    truncated_normal_distribution(T _mu = 0.0,
                                  T _sigma = 1.0,
                                  T _a = -numeric_limits<double>::max() / 2.0,
                                  T _b = numeric_limits<double>::max() / 2.0)
        : m_mu(_mu), m_sigma(_sigma), m_alpha((_a - _mu) / _sigma), m_beta((_b - _mu) / _sigma)
    {
        assert(m_sigma > 0.0);
        assert(m_alpha < m_beta);
    }
    
    void reset() {
        m_distU1.reset();
    }
    
    template<class Generator>
    T operator()(Generator &_g) {
        T u = m_distU1(_g);
        return m_mu + m_sigma * Phi_inv(u * Phi(m_beta) + (1.0 - u) * Phi(m_alpha));
    }
    
    T mean() const {
        return m_mu + m_sigma * A();
    }
    
    T stddev() const {
        T b = (m_alpha * normal_pdf(m_alpha) - m_beta * normal_pdf(m_beta)) / Z();
        T a = A();
        return m_sigma * std::sqrt(1.0 + b - a * a);
    }
    
    T density(T x) const {
        T z = (x - m_mu) / m_sigma;
        if (z >= m_alpha && z <= m_beta) {
            return normal_pdf(z) / (m_sigma * Z());
        }
        else {
            return 0.0;
        }
    }
    
    T cumul(T x) const {
        T z = (x - m_mu) / m_sigma;
        if (z >= m_alpha && z <= m_beta) {
            return (Phi(z) - Phi(m_alpha)) / Z();
        }
        else if (z > m_beta) {
            return 1.0;
        }
        else {
            return 0.0;
        }
    }
    
protected:
    T Z() const {
        return Phi(m_beta) - Phi(m_alpha);
    }
    T A() const {
        return (normal_pdf(m_alpha) - normal_pdf(m_beta)) / Z();
    }
    
    T m_mu;
    T m_sigma;
    T m_alpha;
    T m_beta;
    vg::uniform_real_distribution<T> m_distU1;
};

/// We use this widerer to widen the output of a PRNG that generates only
/// numbers in a smaller range so they cover a wider int type.
template<typename PRNG, typename OutType>
class WideningPRNG {
public:
    
    using result_type = OutType;
    
    WideningPRNG(PRNG& to_widen) : base(to_widen) {
        // Nothing to do
    }
    
    OutType min() const {
        return numeric_limits<OutType>::min();
    }
    
    OutType max() const {
        return numeric_limits<OutType>::max();
    }
    
    /// Generate a random number filling all OutType's bits with random bits.
    OutType operator()() {
        static_assert(is_unsigned<OutType>::value, "OutType must be an unsigned int type");
        static_assert(numeric_limits<long long unsigned int>::digits >= numeric_limits<OutType>::digits, "OutType is too wide to bit count in");
        
        // Work out how wide the base PRNG range is, in total (counting the inclusive bounds)
        // We assume that this will fit in the wider type we are trying to fill.
        // Otherwise you don't need this class.
        OutType base_range = (OutType) base.max() - (OutType) base.min() + 1;
        
        // The base range must be 1 bit at least. Otherwise we will make no progress.
        assert(base_range >= 2);
        
        // Count unset leading bits with this useful compiler builtin.
        // Hope your compiler has it.
        auto unused_bits = __builtin_clzll(base_range);
        auto used_bits = numeric_limits<long long unsigned int>::digits - unused_bits;
        
        // Get just that max used bit
        OutType used_bit = 1 << (used_bits - 1);
        
        // If a generated number from the RNG has this bit flag in it, it has
        // passed the largest complete power of 2 it can make and needs to be
        // rerolled.
        OutType reroll_flag = 0;
        
        if (base_range > used_bit) {
            // We don't cover a power of 2 exactly.
            // If the high bit is set we're going to need to reroll.
            reroll_flag = used_bit;
            // Lop off a bit for computing how many bits we actually got.
            used_bit = used_bit >> 1;
            used_bits--;
        }
        
        assert(used_bits > 0);
        
        OutType generated = 0;
        int generated_bits = 0;
        
        while (generated_bits < numeric_limits<OutType>::digits) {
            // Shift what's there up to make room for new bits.
            generated = generated << used_bits;
            
            OutType new_bits;
            do {
                // Generate bits until we're below the largest power of 2 we can generate above, if any.
                new_bits = (OutType) base() - (OutType) base.min();
            } while (!(new_bits & reroll_flag) && reroll_flag);
            
            // Add in the new bits and record that they are there.
            generated |= new_bits;
            generated_bits += used_bits;
        }
        
        // Now we have a full type full of bits.
        return generated;
    }
    
protected:
    PRNG& base;
};


/// This uniform_int_distribution implementation is based on the
/// UniformRealDistribution from https://stackoverflow.com/a/34962942
template<typename T = int>
class uniform_int_distribution {
public:
    typedef T result_type;
    
    uniform_int_distribution(T _a = 0, T _b = numeric_limits<T>::max()) : m_a(_a), m_b(_b) {
        // Make sure inclusive bounds are valid
        assert(_b >= _a);
    }
    
    void reset() {
        // Also nothing to do!
    }
    
    
    template<class Generator>
    T operator()(Generator &_g) {
        
#ifdef debug
        cerr << "Source range " << _g.min() << " to " << _g.max() << endl;
        cerr << "Dest range " << m_a << " to " << m_b << endl;
#endif
        
        // Define an unsigned widest type to work in
        using WorkType = typename make_unsigned<typename common_type<typename Generator::result_type, T>::type>::type;
        
        // How big are the source and destination ranges?
        // Since they are so big and inclusive we can't always hold their real sizes, so hold size-1
        WorkType source_range_size_minus_1 = (WorkType) _g.max() - (WorkType) _g.min();
        WorkType dest_range_size_minus_1 = (WorkType) m_b - (WorkType) m_a;
        
        if (source_range_size_minus_1 >= dest_range_size_minus_1) {
            // The generator's result is going to be wide enough
            return generate_from_wide_generator(_g);
        } else {
            // The hard way is generating a bigger range from a smaller range.
            // Wrap the generator in something to widen it to our work type
            // and recurse.
            WideningPRNG<Generator, WorkType> widened(_g);
            
            // Generate with that, which had better be wide enough
            return generate_from_wide_generator(widened);
        }
    }
            
    T a() const {
        return m_a;
    }
            
    T b() const {
        return m_b;
    }
        
protected:
        
    /// Generate a result when we know the generator will produce a result on a
    /// range as big as or bigger than ours.
    template<class Generator>
    T generate_from_wide_generator(Generator &_g) {
        // Jordan's strategy: discard anything above the highest multiple of your range, then mod down to your range.

#ifdef debug
        cerr << "Source range " << _g.min() << " to " << _g.max() << endl;
        cerr << "Dest range " << m_a << " to " << m_b << endl;
#endif

        // Define an unsigned widest type to work in
        using WorkType = typename make_unsigned<typename common_type<typename Generator::result_type, T>::type>::type;

        // How big are the source and destination ranges?
        // Since they are so big and inclusive we can't always hold their real sizes, so hold size-1
        WorkType source_range_size_minus_1 = (WorkType) _g.max() - (WorkType) _g.min();
        WorkType dest_range_size_minus_1 = (WorkType) m_b - (WorkType) m_a;

        // We must be generating a smaller range from a bigger rnage here.
        assert(source_range_size_minus_1 >= dest_range_size_minus_1);

        if (dest_range_size_minus_1 == source_range_size_minus_1) {
            // Ranges are the same size. No real work to do.
            return (WorkType) _g() - (WorkType) _g.min() + (WorkType) m_a;
        }

        // Otherwise the ranges differ in size. Which means the dest range must
        // be smaller. Which means the dest range's real size is representable.
        WorkType dest_range_size = dest_range_size_minus_1 + 1;

        // Find how many numbers we have to clip off of the top of the source
        // range so the rest can be covered by tiled destination ranges.
        WorkType remainder = source_range_size_minus_1 % dest_range_size;
        // Change the remainder from source_range_size_minus_1 to the remainder for the actual source range size
        remainder = (remainder + 1) % dest_range_size;

        if (remainder == 0) {
            // We perfectly tiled the source range
            return ((WorkType) _g() - (WorkType) _g.min()) % dest_range_size + (WorkType) m_a;
        }

        // Otherwise there are some values we need to reject

        // Sample a value until we get one that isn't too close to the top of the range.
        WorkType sampled;
        do {
            sampled = (WorkType) _g();
        } while (_g.max() - sampled < remainder);

        // Convert to destination range.
        return (sampled - (WorkType) _g.min()) % dest_range_size + m_a;
    }


    T m_a;
    T m_b;
};
        
/// We provide a partial discrete_distribution implementation that is just the parts we need
template<typename T = int>
class discrete_distribution {
public:
    typedef T result_type;
    typedef double param_type;
    
    template<class InputIt>
    discrete_distribution(InputIt first, InputIt last) : m_weights{first, last} {
        // We can't use an empty weights vector
        assert(!m_weights.empty());
        // Compute partial sums
        std::partial_sum(m_weights.begin(), m_weights.end(), std::back_inserter(m_sums));
    }
    
    discrete_distribution(initializer_list<double> weights = {1}) : discrete_distribution(weights.begin(), weights.end()) {
        // Nothing to do
    }
    
    void reset() {
        // Also nothing to do!
    }
    
    template<class Generator>
    T operator()(Generator &_g) {
    
        // Set up to generate a double from 0 to max weight
        vg::uniform_real_distribution<double> backing_dist(0, m_sums.back());
        // Do it and find which cumumative sum is greater than it
        auto winning_iterator = std::lower_bound(m_sums.begin(), m_sums.end(), backing_dist(_g));
        
        // Find its category number and return that.
        return winning_iterator - m_sums.begin();
    
    }
        
    protected:
    // If we ever want to implement the params stuff we need the weights stored.
    vector<double> m_weights;
    vector<double> m_sums;
    
};
    
// ewen's allele sampling distribution.  for use in genotype prior (as in freebayes)
// gives Pr(a1, ...,an;theta) where ai is the number of sampled haplotypes (out of n) that
// have i different alleles at a given locus. theta is the population mutation rate.
// ex: for a single diploid genotype, a={2,0} = heterozygous: 2 alleles occur once.
//                                    a={0,1} = homozygous: 1 allele occurs twice.
//
// https://en.wikipedia.org/wiki/Ewens%27s_sampling_formula
// https://github.com/ekg/freebayes/blob/master/src/Ewens.cpp#L17
inline real_t ewens_af_prob_ln(const vector<int>& a, real_t theta) {

    // first term (wrt formula as stated on wikipedia)
    // n! / (theta * (theta + 1) * ... (theta + n - 1))
    real_t term1_num_ln = factorial_ln(a.size());
    real_t term1_denom_ln = 0.;
    for (int i = 0; i < a.size(); ++i) {
        term1_denom_ln += log(theta + i);
    }
    real_t term1_ln = term1_num_ln - term1_denom_ln;

    // second term
    // prod [ (theta^aj) / (j^aj * aj!)
    real_t term2_ln = 0.;
    for (int j = 0; j < a.size(); ++j) {
        real_t num = log(pow(theta, a[j]));
        real_t denom = log(pow(1. + j, a[j]) + factorial_ln(a[j]));
        term2_ln += num - denom;
    }

    return term1_ln + term2_ln;
}


}

#endif
