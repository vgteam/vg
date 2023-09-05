#ifndef VG_KMER_HPP_INCLUDED
#define VG_KMER_HPP_INCLUDED

#include <vg/vg.pb.h>
#include <iostream>
#include <atomic>
#include "vg/io/json2pb.h"
#include "handle.hpp"
#include "position.hpp"
#include "gcsa/gcsa.h"

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

/**
 * Exception that indicates that a limit on disk size has been exceeded
 */
class SizeLimitExceededException : public std::exception {
public:
    
    SizeLimitExceededException() noexcept = default;
    ~SizeLimitExceededException() noexcept = default;
    
    const char* what() const noexcept;
private:
    
    static const string msg;
};

/// Iterate over all the kmers in the graph, running lambda on each
/// If the stop flag is included, stop execution if it ever evaluates to true
void for_each_kmer(const HandleGraph& graph, size_t k,
                   const function<void(const kmer_t&)>& lambda,
                   id_t head_id = 0, id_t tail_id = 0,
                   atomic<int>* stop_flag = nullptr);

/// Print a kmer_t to a stream.
ostream& operator<<(ostream& out, const kmer_t& kmer);

/// Convert the kmer_t to a set of gcsa2 binary kmers which are exposed via a callback.
void kmer_to_gcsa_kmers(const kmer_t& kmer, const gcsa::Alphabet& alpha, const function<void(const gcsa::KMer&)>& lambda);

/// Encode the chars into the gcsa2 byte
gcsa::byte_type encode_chars(const vector<char>& chars, const gcsa::Alphabet& alpha);

/**
 * Write GCSA2 formatted binary KMers to the given ostream.
 * size_limit is the maximum size of the kmer file in bytes. When the function
 * returns, size_limit is the size of the kmer file in bytes.
 */
void write_gcsa_kmers(const HandleGraph& graph, int kmer_size, ostream& out, size_t& size_limit, id_t head_id, id_t tail_id);

/// Open a tempfile and write the kmers to it. The calling context should remove it
/// with temp_file::remove(). In the case that the size limit is exceeded, throws a
/// SizeLimitExceededException and deletes the temp file.
string write_gcsa_kmers_to_tmpfile(const HandleGraph& graph, int kmer_size, size_t& size_limit, id_t head_id, id_t tail_id,
                                   const string& base_file_name = "vg-kmers-tmp-");

}

#endif
