#ifndef VG_DIPLOID_PAIR_SCORER_HPP_INCLUDED
#define VG_DIPLOID_PAIR_SCORER_HPP_INCLUDED

/** \file diploid_pair_scorer.hpp
 *
 * Scores of pairs of haplotype sequences in a subchain as the two haplotypes of
 * a diploid sample, which diploid sampling uses to choose the pair.
 */

#include "../src/recombinator_haplotypes.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace vg {

//------------------------------------------------------------------------------

/**
 * A judge that scores pairs of haplotype sequences of one subchain as the two
 * haplotypes of a diploid sample.
 *
 * A pair holds 0, 1 or 2 copies of each kmer. An absent kmer scores 1, 0 or -1
 * for those copies, a present kmer -1, 0 or 1, and a heterozygous kmer 1 for
 * one copy and 0 otherwise. Frequent and ignored kmers score 0. The score of the
 * pair is the sum over the kmers, which is what
 * `KmerPresenceMatrix::diploid_score()` returns.
 *
 * The judge decodes the presence row of each sequence once and keeps a mask of
 * the kmers in each scored class. It then scores a pair with five popcounts per
 * 64 kmers.
 */
class DiploidPairScorer {
public:
    /// Decodes the presence rows of the given sequences of the subchain, and the
    /// classes of its kmers. `kmer_types` has an entry for each kmer of the
    /// subchain.
    DiploidPairScorer(
        const Haplotypes::Subchain& subchain, const std::vector<size_t>& sequences,
        const std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>& kmer_types
    );

    /// Returns the score of the pair of sequences at the given positions in the
    /// list given to the constructor.
    std::int64_t score(size_t first, size_t second) const;

private:
    /// The kmers in each scored class among 64 consecutive kmers. Kmer `i` of
    /// the subchain is bit `i % 64` of word `i / 64`.
    struct ClassMasks {
        std::uint64_t absent = 0;
        std::uint64_t heterozygous = 0;
        std::uint64_t present = 0;
    };

    /// The number of 64-bit words in a presence row.
    size_t words;

    /// The class masks of each word.
    std::vector<ClassMasks> classes;

    /// The presence rows of the sequences, one after another, in the bit order
    /// of `ClassMasks`.
    std::vector<std::uint64_t> rows;
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_DIPLOID_PAIR_SCORER_HPP_INCLUDED
