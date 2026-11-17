#ifndef VG_DEDUP_INTERNAL_HPP_INCLUDED
#define VG_DEDUP_INTERNAL_HPP_INCLUDED

/**
 * \file dedup_internal.hpp
 *
 * Internal tools for content-aware node chopping (see dedup.hpp).
 */

#include "dedup.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace vg {

//------------------------------------------------------------------------------

/// A k-mer with a locally minimal hash. A candidate position for chopping.
struct KmerAnchor {
    /// Starting offset of the k-mer in the sequence.
    size_t offset;
    /// The k-mer encoded using 2 bits/base, with the last base in the low bits.
    std::uint64_t kmer;
};

/**
 * Returns the k-mers of the sequence whose hash is strictly smaller than the
 * hash of every other k-mer starting at most `radius` bases away, in sequence
 * order. Hashing is case insensitive and only considers the bases
 * A, C, G and T.
 *
 * Similar to finding forward strand minimizers with window length
 * `2 * radius + 1`.
 *
 * @pre `kmer_length` is between 1 and 32.
 */
std::vector<KmerAnchor> find_anchors(std::string_view sequence, size_t kmer_length, size_t radius);

//------------------------------------------------------------------------------

/// Returns the fingerprint of the entire sequence, as computed by
/// `SubstringFingerprinter`.
std::uint64_t sequence_fingerprint(std::string_view sequence);

/**
 * A precomputed table for quickly computing the fingerprint of any substring
 * `s` of length `n`, up to a maximum length.
 *
 * The fingerprint is the polynomial `s[0] * B^(n-1) + ... + s[n-1]` modulo
 * the prime `2^61 - 1`, for a fixed base `B`. Two different strings of length
 * at most `n` have the same fingerprint with probability at most `n / 2^61`.
 */
class SubstringFingerprinter {
public:
    /// Builds the table for substrings of length up to `max_length`.
    SubstringFingerprinter(std::string_view sequence, size_t max_length);

    /// Returns the fingerprint of the substring [begin, end).
    std::uint64_t fingerprint(size_t begin, size_t end) const;

private:
    /// The fingerprint of each prefix of the sequence, by its length.
    std::vector<std::uint64_t> prefix_fingerprints;
    /// `B^i` modulo the prime, for each `i` up to the maximum length.
    std::vector<std::uint64_t> base_powers;
};

//------------------------------------------------------------------------------

/// A set of candidate positions for chopping. The set consists of candidates
/// returned by `find_anchors()`, with filler positions every `max_node_length`
/// bases if the gap would be longer than that.
class CutCandidateSet {
public:
    /// Finds the candidate positions in the sequence.
    CutCandidateSet(std::string_view sequence, const ContentChopParameters& parameters);

    /// Returns the candidate positions, including sentinels at the start and
    /// end of the sequence.
    const std::vector<size_t>& offsets() const { return this->cut_offsets; }

    /// Return the maximum length of a segment.
    size_t max_segment_length() const { return this->max_length; }

    /// Returns the fingerprint of the segment between candidates `from` and
    /// `to`, given as indexes into `offsets()`.
    ///
    /// @pre `from < to`, and the candidates are at most the maximum segment
    /// length apart.
    std::uint64_t segment_fingerprint(size_t from, size_t to) const {
        return this->fingerprinter.fingerprint(this->cut_offsets[from], this->cut_offsets[to]);
    }

    /// Calls `visit(from, to)` for each segment, where `from < to`
    /// are indexes into `offsets()` at most `max_segment_length()` apart.
    template<typename Visitor>
    void for_each_segment(const Visitor& visit) const {
        for (size_t from = 0; from < this->cut_offsets.size(); from++) {
            for (size_t to = from + 1;
                 to < this->cut_offsets.size() && this->cut_offsets[to] - this->cut_offsets[from] <= this->max_length;
                 to++) {
                visit(from, to);
            }
        }
    }

private:
    std::vector<size_t> cut_offsets;
    size_t max_length;
    SubstringFingerprinter fingerprinter;
};

//------------------------------------------------------------------------------

/// A map from segment fingerprints to their occurrence counts.
class SegmentCountTable {
public:
    /// Builds the mapping from an arbitrary list of fingerprints.
    /// The list will be sorted using multiple OpenMP threads.
    explicit SegmentCountTable(std::vector<std::uint64_t> fingerprints);

    /// Returns the number of occurrences of the segment with this fingerprint.
    size_t count(std::uint64_t fingerprint) const;

    /// Return the number of distinct fingerprints.
    size_t size() const { return this->entries.size(); }

private:
    /// A fingerprint and its number of occurrences.
    struct Entry {
        std::uint64_t fingerprint;
        size_t count;
    };

    /// One entry for each distinct fingerprint, in sorted order.
    std::vector<Entry> entries;
};

/**
 * Selects chopping positions from `candidates`, in sorted order. Returned
 * offsets do not include the sentinels at the start and end of the sequence.
 *
 * Assumes that `counts` has a non-zero count for every segment. The cost of
 * a segment is `piece_weight + length / occurrences`. Since the number of
 * occurrences is based on fingerprints, rare collisions may lead to
 * suboptimal choices.
 *
 * NOTE: Dynamic programming prioritizes the leftmost minimum-cost chopping.
 * In particular, if we have a long node with unique sequence, the first piece
 * may be below max length, while subsequent pieces are max length. Chopping
 * with vg mod --chop does the opposite.
 */
std::vector<size_t> choose_cuts(
    const CutCandidateSet& candidates,
    const SegmentCountTable& counts,
    double piece_weight
);

//------------------------------------------------------------------------------

}

#endif
