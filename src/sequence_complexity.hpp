/**
 * \file sequence_complexity.hpp
 *
 * Defines and implements an algorithm to identify low-complexity sequences
 *
 */
#ifndef VG_SEQUENCE_COMPLEXITY_HPP_INCLUDED
#define VG_SEQUENCE_COMPLEXITY_HPP_INCLUDED

#include <string>
#include <cmath>
#include <iostream>

namespace vg {

using namespace std;

//#define debug_seq_complexity

/*
 * Struct to compute the complexity of sequence at different orders
 */
template<int MaxOrder = 1>
struct SeqComplexity {
    
    SeqComplexity(string::const_iterator begin, string::const_iterator end);
    SeqComplexity(const string& seq);
    
    // The approximate p-value of the n-th order correlation between
    // nucleotides. Only valid for 1 <= order <= MaxOrder template param.
    double p_value(int order) const;
    
    // The fraction of pairs that are repeats at this order.
    double repetitiveness(int order) const;
    
private:
    
    int len;
    int matches[MaxOrder];
};

/*
 * Template implementations
 */

template<int MaxOrder>
SeqComplexity<MaxOrder>::SeqComplexity(const string& seq) : SeqComplexity(seq.begin(), seq.end()) {
    
}

template<int MaxOrder>
SeqComplexity<MaxOrder>::SeqComplexity(string::const_iterator begin, string::const_iterator end) {

    len = end - begin;
    
    for (int i = 0; i < MaxOrder; ++i) {
        matches[i] = 0;
    }
    
    for (int i = 1; i < len; ++i) {
        for (int j = max<int>(0, i - MaxOrder); j < i; ++j) {
            matches[i - j - 1] += (*(begin + j) == *(begin + i));
        }
    }
#ifdef debug_seq_complexity
    cerr << "match table for seq of length " << len << ":" << endl;
    for (int i = 1; i < len; ++i) {
        cerr << i << ": " << matches[i - 1] << endl;
    }
#endif
}

// TODO: have GC bias instead of uniform random? maybe not since repetition
// hurts alignment uncertainty even in biased sequencesc

template<int MaxOrder>
double SeqComplexity<MaxOrder>::p_value(int order) const {
    if (order < len && order + 8 > len) {
        // exact binomial CDF
        // TODO: flip sum if is smaller?
        double x = 1.0;
        double y = len - order;
        double term = pow(0.75, y);
        double accum = 0.0;
        for (int i = 0, k = matches[order - 1]; i < k; ++i) {
            accum += term;
            // the ratio of successive terms in the binomial distr pmf
            term *= 0.333333333333333333 * y / x;
            x += 1.0;
            y -= 1.0;
        }
        
        return 1.0 - accum;
    }
    else if (order < len) {
        // normal approximation to binomial
        static const double root_pq = sqrt(0.25 * 0.75);
        static const double root_1_2 = sqrt(0.5);
        double z = (double(matches[order - 1]) - double(len - order) * 0.25) / (sqrt(len - order) * root_pq);
        // the normal CDF
        return 1.0 - 0.5 * erfc(-root_1_2 * z);
    }
    else {
        return 1.0;
    }
}

template<int MaxOrder>
double SeqComplexity<MaxOrder>::repetitiveness(int order) const {
    return len > order ? double(matches[order - 1]) / double(len - order) : 0.0;
}


}

#endif
