#ifndef CALLER_H
#define CALLER_H

#include <iostream>
#include <algorithm>
#include <functional>
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

    Caller(int buffer_size = 1000, int min_depth = 5, double min_frac = 0.3) :
        _buffer_size(buffer_size),
        _min_depth(min_depth),
        _min_frac(min_frac){
    }
        
    // copy constructor
    Caller(const Caller& other) : _alignments(other._alignments),
                                  _buffer_size(other._buffer_size),
                                  _min_depth(other._min_depth),
                                  _min_frac(other._min_frac){
    }

    // move constructor
    Caller(Caller&& other) : _alignments(other._alignments),
                             _buffer_size(other._buffer_size),
                             _min_depth(other._min_depth),
                             _min_frac(other._min_frac){
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
    // create a snp from a pileup.
    void create_snp(int64_t node_id, int64_t offset, char base);

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
};



}

#endif
