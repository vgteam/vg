#ifndef VG_KFF_HPP_INCLUDED
#define VG_KFF_HPP_INCLUDED

/** \file 
 * Tools for working with the Kmer File Format (KFF).
 */

#include <deque>
#include <mutex>

#include <kff_io.hpp>

#include <gbwtgraph/minimizer.h>

namespace vg {

//------------------------------------------------------------------------------

/// A mapping of character values from KFF encoding to minimizer index encoding.
struct kff_recoding_t {
    uint8_t data[4];
};

/// Returns the number of bytes required for a kmer in KFF format.
inline size_t kff_bytes(size_t k) {
    return (k + 3) / 4;
}

/// Returns `true` if the encoding is trivial (0, 1, 2, 3).
bool kff_is_trivial(const uint8_t* encoding);

/// Inverts the KFF encoding into a packed -> char table.
std::string kff_invert(const uint8_t* encoding);

/// Returns a recoding for the given encoding.
kff_recoding_t kff_recoding(const uint8_t* encoding);

/// Parses a big-endian integer from KFF data.
uint64_t kff_parse(const uint8_t* data, size_t bytes);

//------------------------------------------------------------------------------

/// Encodes a kmer in KFF format according to the given encoding.
/// Non-ACGT characters are encoded as 0s.
std::vector<uint8_t> kff_encode(const std::string& kmer, const uint8_t* encoding);

/// Decodes a kmer in KFF format according to the given encoding.
std::string kff_decode(const uint8_t* kmer, size_t k, const std::string& decoding);

//------------------------------------------------------------------------------

/// Recodes a kmer from a minimizer index in KFF format according to the given encoding.
std::vector<uint8_t> kff_recode(gbwtgraph::Key64::value_type kmer, size_t k, const uint8_t* encoding);

/// Recodes a KFF kmer in the minimizer index format according to the given encoding.
/// Will fail silently if `k` is too large or `recoding` is not from `kff_recoding()`.
gbwtgraph::Key64::value_type kff_recode(const uint8_t* kmer, size_t k, kff_recoding_t recoding);

/// Recodes a KFF kmer in the minimizer index format, assuming that the encoding is
/// the same. Will fail silently if `k` or `bytes` is too large.
gbwtgraph::Key64::value_type kff_recode_trivial(const uint8_t* kmer, size_t k, size_t bytes);

/// Recodes `n` KFF kmers in the minimizer index format according to the given encoding.
/// Will fail silently if `k` is too large or `recoding` is not from `kff_recoding()`.
std::vector<gbwtgraph::Key64::value_type> kff_recode(const uint8_t* kmers, size_t n, size_t k, kff_recoding_t recoding);

//------------------------------------------------------------------------------

/// Returns the reverse complement of a KFF kmer.
std::vector<uint8_t> kff_reverse_complement(const uint8_t* kmer, size_t k, const uint8_t* encoding);

/// Returns the reverse complement of a minimizer index kmer.
inline gbwtgraph::Key64::value_type minimizer_reverse_complement(gbwtgraph::Key64::value_type kmer, size_t k) {
    return gbwtgraph::Key64(kmer).reverse_complement(k).get_key();
}

//------------------------------------------------------------------------------

/**
 * A wrapper over `Kff_reader` that allows reading kmers safely from multiple threads.
 */
class ParallelKFFReader {
public:
    typedef gbwtgraph::Key64::value_type kmer_type;

    /// Creates a new reader for the given file. Throws `std::runtime_error` if the
    /// file cannot be opened or the data is unacceptable.
    ParallelKFFReader(const std::string& filename);

    /// Reads the next `n` kmers and counts from the file. This can be called safely
    /// from multiple threads. If the returned vector contains fewer than `n` kmers
    /// this indicates that the reader has reached the end of the file.
    std::vector<std::pair<kmer_type, size_t>> read(size_t n);

    /// KFF reader.
    Kff_reader reader;

    /// Buffer from unused kmers from the latest block.
    std::deque<std::pair<kmer_type, size_t>> buffer;

    /// Mutex for accessing `reader` and `buffer`.
    std::mutex mtx;

    /// Length of the kmers.
    size_t k;

    /// Maximum number of kmers per block.
    size_t max_kmers_per_block;

    /// Number of bytes reserved for each kmer count.
    size_t data_bytes;

    /// Encoding used for the kmers.
    std::uint8_t encoding[4];

    /// Recoding from KFF kmers to minimizer index kmers.
    kff_recoding_t recoding;

private:
    ParallelKFFReader(const ParallelKFFReader&) = delete;
    ParallelKFFReader& operator= (const ParallelKFFReader&) = delete;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_KFF_HPP_INCLUDED
