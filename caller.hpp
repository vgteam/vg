#ifndef CALLER_H
#define CALLER_H

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include "vg.pb.h"
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "pileup.hpp"

namespace vg {

using namespace std;

// Super simple variant caller, for now written to get bakeoff evaluation bootstrapped.
// Idea: Idependently process Pileup records, using simple thresholding to make calls.
// Output stored as GAM. 
class Caller {
public:

    static const int Default_buffer_size;
    static const int Default_min_depth;
    static const double Default_min_frac;
    static const double Default_min_likelihood;
    static const int Default_context;
    static const char Default_default_quality;
    
    Caller(int buffer_size = Default_buffer_size,
           int min_depth = Default_min_depth,
           double min_frac = Default_min_frac,
           double _min_likelihood = Default_min_likelihood, 
           int context = Default_context,
           int default_quality = Default_default_quality) :
        _buffer_size(buffer_size),
        _min_depth(min_depth),
        _min_frac(min_frac),
        _min_log_likelihood(safe_log(_min_likelihood)),
        _context(context),
        _default_quality(default_quality) {
    }
        
    // copy constructor
    Caller(const Caller& other) : _alignments(other._alignments),
                                  _buffer_size(other._buffer_size),
                                  _min_depth(other._min_depth),
                                  _min_frac(other._min_frac),
                                  _min_log_likelihood(other._min_log_likelihood),
                                  _context(other._context),
                                  _default_quality(other._default_quality) {
    }

    // move constructor
    Caller(Caller&& other) : _alignments(other._alignments),
                             _buffer_size(other._buffer_size),
                             _min_depth(other._min_depth),
                             _min_frac(other._min_frac),
                             _min_log_likelihood(other._min_log_likelihood),
                             _context(other._context),
                             _default_quality(other._default_quality) {
        other._alignments.clear();
    }

    // copy assignment operator
    Caller& operator=(const Caller& other) {
        Caller tmp(other);
        *this = move(tmp);
        return *this;
    }

    // move assignment operator
    Caller& operator=(Caller&& other) {
        _buffer_size = other._buffer_size;
        _min_depth = other._min_depth;
        _min_frac = other._min_frac;
        _min_log_likelihood = other._min_log_likelihood;
        _context = other._context;
        _default_quality = other._default_quality;
        swap(_alignments, other._alignments);
        other._alignments.clear();
        return *this;
    }

    // delete contents of table
    ~Caller() {
        clear();
    }
    void clear();

    // buffer of output alignments (how we're storing variation for now)
    vector<Alignment> _alignments;
    // maximum size _alignments can get to before writing out output stream
    int _buffer_size;
    // minimum depth of pileup to call variants on
    int _min_depth;
    // minimum fraction of bases in pileup that nucleotide must have to be snp
    double _min_frac;
    // minimum log likelihood for a snp to be called
    double _min_log_likelihood;
    // vg surject doesn't work when alignments conly contain snps.  try adding
    // some matching context on each side to see if it helps.  this specifies
    // context length (on each side)
    int _context;
    // if we don't have a mapping quality for a read position, use this
    char _default_quality;

    // write GAM to JSON
    void to_json(ostream& out);
    // write GAM to protobuf
    void write(ostream& out);
    // write out the _alignments buffer
    void flush_buffer(ostream& out, bool json);

    // call every position in the node pileup
    void call_node_pileup(const NodePileup& pileup, ostream& out, bool json);
    // call position at given base
    void call_base_pileup(const NodePileup& np, int64_t offset);
    
    // Find the top-two bases in a pileup, along with their counts
    void compute_top_frequencies(const BasePileup& bp,
                                 const vector<pair<int, int> >& base_offsets,
                                 char& top_base, int& top_count,
                                 char& second_base, int& second_count);
    
    // compute genotype for a position with maximum likelihood
    pair<char, char> ml_snp_genotype(const BasePileup& bp,
                                     const vector<pair<int, int> >& base_offsets,
                                     char top_base, char second_base);
    // compute likelihood of a genotype samtools-style
    double genotype_log_likelihood(const BasePileup& pb,
                                   const vector<pair<int, int> >& base_offsets,
                                   double g, char first, char second);

    // create a snp from a pileup. The snp is stored as an Alignment
    // in _alignments.  When possible, the snp is padded by _context
    // exact matches on each side, to help vg mod -i merge it. 
    void create_snp(const NodePileup& np, int64_t offset, char base);

    // convert nucleotide base into integer
    static int nidx(char c) {
        c = ::toupper(c);
        switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        }
        return 4;
    }

    // and back
    static char idxn(int i) {
        return "ACGTN"[i];
    }

    // log function that tries to avoid 0s
    static double safe_log(double v) {
        return ::log(std::max(v, numeric_limits<double>::epsilon()));
    }

    // tranform ascii phred into probability of error
    static double phred2prob(char q) {
        //q -= 33; dont think we need this after all
        assert(q >= 0);
        // error prob = 10^(-PHRED/10)
        return pow(10., -(double)q / (10.));
    }
};



}

#endif
