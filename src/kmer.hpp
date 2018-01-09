#ifndef VG_KMER_HPP_INCLUDED
#define VG_KMER_HPP_INCLUDED

#include "vg.pb.h"
#include <iostream>
#include "json2pb.h"
#include "handle.hpp"
#include "position.hpp"

/** \file 
 * Functions for working with `kmers_t`'s in HandleGraphs.
 */

namespace vg {

using namespace std;

/// Stores a kmer in the context of a graph.
/// The context may be used to double the kmer length
/// of a deBruijn graph encoded in these,
/// which is done in GCSA2's index construction.
struct kmer_t {
    kmer_t(const string& s,
           const pos_t& b,
           const pos_t& e,
           const handle_t& c)
        : seq(s), begin(b), end(e), curr(c) { };
    /// the kmer
    string seq;
    /// our start position
    pos_t begin;
    /// Used in construction
    pos_t end; /// one past the (current) end of the kmer
    handle_t curr; /// the next handle we extend into
    /// These are filled when construction is complete
    vector<pos_t> prev_pos; // previous positions and their chars
    vector<pos_t> next_pos;
    vector<char> prev_char; // next positions and their chars
    vector<char> next_char;
};

/// Iterate over all the kmers in the graph, running lambda on each
void for_each_kmer(const HandleGraph& graph, size_t k,
                   const function<void(const kmer_t&)>& lambda,
                   id_t head_id = 0, id_t tail_id = 0);

/// Print a kmer_t to a stream.
ostream& operator<<(ostream& out, const kmer_t& kmer);

}

#endif
