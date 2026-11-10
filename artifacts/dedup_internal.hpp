#ifndef VG_DEDUP_INTERNAL_HPP_INCLUDED
#define VG_DEDUP_INTERNAL_HPP_INCLUDED

/**
 * \file dedup_internal.hpp
 * The parts of the content-aware chopping planner, for its implementation and its tests.
 *
 * Callers of plan_content_chopping() need only dedup.hpp. This header declares the tools that
 * decide where to cut each long node: find_anchors() and CutCandidateSet pick the candidate cuts,
 * SegmentCountTable counts the candidate segments, and choose_cuts() chooses the cuts in a node.
 */

#include "dedup.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace vg {

//------------------------------------------------------------------------------

/// A k-mer whose hash is the strict minimum of the hashes of the k-mers near it.
struct KmerAnchor {
    /// The offset of the first base of the k-mer in the sequence.
    size_t offset;
    /// The k-mer with 2 bits per base, the last base in the lowest bits.
    uint64_t kmer;
};

/**
 * Return the k-mers of the sequence whose hash is smaller than the hash of every other k-mer
 * starting at most `radius` bases away, in sequence order.
 *
 * The k-mers are case-insensitive and contain only the bases A, C, G and T. A k-mer near an end
 * of the sequence is compared only with the k-mers inside the sequence. A k-mer that occurs
 * twice within `radius` bases is never an anchor, so a homopolymer run, or a tandem repeat with
 * a period of at most `radius`, has no anchors inside it.
 *
 * Whether a k-mer is an anchor depends only on the k-mers starting at most `radius` bases from
 * it. If two sequences share a stretch that holds all of them, the k-mer is an anchor in both
 * sequences or in neither.
 *
 * @pre `kmer_length` is between 1 and 32.
 */
std::vector<KmerAnchor> find_anchors(std::string_view sequence, size_t kmer_length, size_t radius);

//------------------------------------------------------------------------------

/**
 * Return the fingerprint of the sequence. It equals the fingerprint that a
 * SubstringFingerprinter gives for any substring with the same bases.
 */
uint64_t sequence_fingerprint(std::string_view sequence);

/**
 * A SubstringFingerprinter is a table of the fingerprints of the prefixes of one sequence, which
 * gives the fingerprint of any substring up to a maximum length in constant time.
 *
 * The fingerprint of a string `s` of length `n` is the polynomial `s[0] * B^(n-1) + ... + s[n-1]`
 * modulo the prime `2^61 - 1`, for a fixed base `B`. The characters are their byte values. Equal
 * strings have equal fingerprints, and two different strings of length at most `n` have the same
 * fingerprint with probability at most `n / 2^61`. Case matters.
 */
class SubstringFingerprinter {
public:
    /// Build the table for the sequence. Substrings up to `max_length` bases long can be
    /// fingerprinted.
    SubstringFingerprinter(std::string_view sequence, size_t max_length);

    /// Return the fingerprint of the substring [begin, end).
    /// @pre `begin <= end`, `end` is at most the sequence length, and `end - begin` is at most
    /// the maximum length.
    uint64_t fingerprint(size_t begin, size_t end) const;

private:
    /// The fingerprint of each prefix of the sequence, by its length.
    std::vector<uint64_t> prefix_fingerprints;
    /// `B^i` modulo the prime, for each `i` up to the maximum length.
    std::vector<uint64_t> base_powers;
};

//------------------------------------------------------------------------------

/**
 * A CutCandidateSet is the set of offsets where a node's sequence may be cut, with the
 * fingerprints of the segments between them. A segment is the sequence between two candidate
 * cuts at most the maximum node length apart.
 *
 * The candidates are the start and the end of the sequence, the offsets of its k-mer anchors,
 * and filler offsets. Wherever two consecutive candidates would be more than the maximum piece
 * length apart, a filler goes at that length after the first of them. So the gap between
 * consecutive candidates is at most the maximum node length, and a chopping into pieces of
 * that length always exists.
 */
class CutCandidateSet {
public:
    /// Find the candidate cuts in the sequence.
    CutCandidateSet(std::string_view sequence, const ContentChopParameters& parameters);

    /// Return the candidate cuts, in increasing order. The first is 0 and the last is the length
    /// of the sequence. An empty sequence has just one candidate.
    const std::vector<size_t>& offsets() const { return cut_offsets; }

    /// Return the maximum length of a segment.
    size_t max_segment_length() const { return max_length; }

    /// Return the fingerprint of the segment between candidates `from` and `to`, given as
    /// indexes into offsets().
    /// @pre `from < to`, and the candidates are at most the maximum segment length apart.
    uint64_t segment_fingerprint(size_t from, size_t to) const {
        return fingerprinter.fingerprint(cut_offsets[from], cut_offsets[to]);
    }

    /// Call `visit(from, to)` for each segment, where `from < to` are indexes into offsets().
    template<typename Visitor>
    void for_each_segment(const Visitor& visit) const {
        for (size_t from = 0; from < cut_offsets.size(); from++) {
            for (size_t to = from + 1;
                 to < cut_offsets.size() && cut_offsets[to] - cut_offsets[from] <= max_length;
                 to++) {
                visit(from, to);
            }
        }
    }

private:
    std::vector<size_t> cut_offsets;
    size_t max_length;
    /// Gives the fingerprints of the segments.
    SubstringFingerprinter fingerprinter;
};

//------------------------------------------------------------------------------

/**
 * A SegmentCountTable is a table of how many times each segment occurs, keyed by the
 * fingerprint of its sequence.
 */
class SegmentCountTable {
public:
    /// Count the fingerprints, one for each occurrence of a segment.
    explicit SegmentCountTable(std::vector<uint64_t> fingerprints);

    /// Return the number of occurrences of the segment with this fingerprint, or 0 if there are
    /// none.
    size_t count(uint64_t fingerprint) const;

    /// Return the number of distinct fingerprints.
    size_t size() const { return entries.size(); }

private:
    /// A fingerprint and its number of occurrences.
    struct Entry {
        uint64_t fingerprint;
        size_t count;
    };

    /// One entry for each distinct fingerprint, in increasing order of fingerprint.
    std::vector<Entry> entries;
};

/**
 * Return the offsets at which to cut a node, in increasing order, without 0 and the length of
 * the node.
 *
 * Choose the chopping into segments of the candidate set that has the least total cost. A
 * segment of length `L` that occurs `c` times costs `piece_weight + L / c`, where `c` is at
 * least 1. Among choppings of equal cost, prefer longer pieces at the end of the node.
 */
std::vector<size_t> choose_cuts(const CutCandidateSet& candidates, const SegmentCountTable& counts,
                                double piece_weight);

//------------------------------------------------------------------------------

}

#endif
