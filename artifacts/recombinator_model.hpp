#ifndef VG_RECOMBINATOR_MODEL_HPP_INCLUDED
#define VG_RECOMBINATOR_MODEL_HPP_INCLUDED

/** \file recombinator_model.hpp
 *
 * A model of the kmer counts in a sample and a hidden Markov model over pairs of
 * haplotypes. Haplotype sampling uses them to choose the first two haplotypes in
 * each subchain (`Recombinator::Parameters::model_pair`).
 */

#include "../src/recombinator_haplotypes.hpp"

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace vg::pair_model {

//------------------------------------------------------------------------------

/**
 * The count of a kmer in the sample's reads, given the number of copies n of the
 * kmer in the two haplotypes of the sample. The count is negative binomial with
 * mean `mu * n` and size `r1 * n` when n > 0, and with mean `eps` and size `r0`
 * when n = 0; the latter covers sequencing errors and kmers from the rest of the
 * genome. Counts below `MIN_COUNT` are reported as 0, because kmer counters drop
 * them, and counts of at least `MAX_COUNT` are reported as `MAX_COUNT`.
 */
struct CountModel {
    /// Counts below this are reported as 0.
    constexpr static size_t MIN_COUNT = 2;

    /// Counts of at least this are reported as this.
    constexpr static size_t MAX_COUNT = 255;

    /// The fit is a mixture over copy numbers 0 to this.
    constexpr static size_t FIT_COPIES = 8;

    /// Mean count per copy.
    double mu = 0.0;

    /// Mean count with no copies.
    double eps = 0.0;

    /// Negative binomial size per copy.
    double r1 = 0.0;

    /// Negative binomial size with no copies.
    double r0 = 0.0;

    /// Mixture weights of copy numbers 0 to `FIT_COPIES` in the fit.
    std::vector<double> pi;

    /// Return the size and the mean of the negative binomial distribution of the
    /// count given the number of copies.
    std::pair<double, double> nb_parameters(size_t copies) const;

    /// Return the log-probability of the reported count given the number of copies.
    double log_probability(size_t count, size_t copies) const;

    /**
     * Fit the model by maximum likelihood to a histogram of reported counts, where
     * `histogram[y]` is the number of kmers with count `y` for `y <= MAX_COUNT`. The
     * mixture weights are fitted with the other parameters. The coverage is the
     * expected count of a kmer with two copies. Counts above `FIT_COPIES / 2` times
     * the coverage (at most `MAX_COUNT - 1`) are only known to be above it, so that
     * repeats are explained by copy numbers and not by the tail of the distribution
     * for no copies. The mean per copy is fitted between 1/8 and 1 times the coverage.
     */
    static CountModel fit(const std::vector<double>& histogram, double coverage);
};

//------------------------------------------------------------------------------

/// The largest number of copies of a kmer in a pair of haplotypes.
constexpr size_t PAIR_COPIES = 2 * Haplotypes::Subchain::MAX_COPIES;

/// With this probability, a count is a random draw from the histogram instead of a
/// draw from the count model. This limits the effect of any single kmer on a score.
constexpr double CONTAMINATION = 0.001;

/// Log-probabilities of reported counts given the copies in a pair of haplotypes:
/// `table[y][m]` for count `y` (0 to `CountModel::MAX_COUNT`) and `m` copies (0 to
/// `PAIR_COPIES`).
typedef std::vector<std::array<double, PAIR_COPIES + 1>> CountTable;

/// Return the table of log-probabilities of the count model, with `CONTAMINATION`
/// from the histogram.
CountTable count_table(const CountModel& model, const std::vector<double>& histogram);

/// A kmer count that marks a kmer the pair scores ignore.
constexpr size_t IGNORED_COUNT = static_cast<size_t>(-1);

/// The number of best pairs in a subchain that the candidates are taken from.
constexpr size_t CANDIDATE_PAIRS = 400;

/// The largest number of candidate groups in a subchain.
constexpr size_t CANDIDATES = 32;

/**
 * The candidate groups of haplotypes in one subchain and the scores of their pairs.
 * The sequences in a group have identical kmer presence and copies. For each group,
 * `members` and `haplotypes` have one entry per sequence of the group.
 */
struct Candidates {
    /// Sequence offsets in the subchain for each group.
    std::vector<std::vector<size_t>> members;

    /// For each group, an identifier of the haplotype each member belongs to.
    /// The Viterbi algorithm uses them to follow haplotypes between subchains.
    std::vector<std::vector<std::uint64_t>> haplotypes;

    /// The score of the pair of groups `(i, j)` at `i * size() + j`: the
    /// log-likelihood of the kmer counts given the copies in the pair.
    std::vector<double> scores;

    /// The length of the subchain in bp, which the transition probabilities use.
    double length = 0.0;

    /// Return the number of groups.
    size_t size() const { return this->members.size(); }
};

/// The kmers with more than one copy in a sequence.
typedef Haplotypes::Subchain::kmer_copies_type KmerCopies;

/**
 * Score the pairs of haplotypes in a subchain of the given length. Each sequence is
 * given as a row of kmer presence bits (as in `Haplotypes::Subchain::presence_row()`),
 * its `KmerCopies`, and the identifier of the haplotype it belongs to. The counts are
 * the sample's counts of the kmers, at most `CountModel::MAX_COUNT`, or
 * `IGNORED_COUNT`. Only the usable sequences become candidates.
 *
 * Sequences with identical kmer presence and copies form a group. The candidates are
 * the first `CANDIDATES` groups to appear in the best `CANDIDATE_PAIRS` pairs, ranked
 * with the copies clipped to one, and they are returned in that order. Pairs with
 * equal scores are ranked by the first sequences of their groups. The scores of the
 * candidate pairs use the copies.
 */
Candidates score_pairs(
    const std::vector<std::vector<std::uint64_t>>& rows,
    const std::vector<KmerCopies>& copies,
    const std::vector<std::uint64_t>& haplotype_ids,
    const std::vector<bool>& usable,
    const std::vector<size_t>& counts,
    const CountTable& table,
    double length
);

/**
 * Return the most likely path of ordered pairs of candidate groups along a chain of
 * subchains, as (first group, second group) in each subchain. In each pair, each
 * member keeps its haplotype into the next subchain with probability 1 - r and
 * otherwise moves to a candidate group chosen uniformly at random, where
 * r = 1 - exp(-rho * d) for the mean length d of the two subchains. A member whose
 * haplotype is not among the candidates in the next subchain also moves uniformly.
 * Pair scores are divided by `tau`. Throws `std::invalid_argument` if a subchain has
 * no candidates.
 */
std::vector<std::pair<size_t, size_t>> viterbi(const std::vector<Candidates>& subchains, double tau, double rho);

/**
 * Choose a sequence offset for each member of the pair in each subchain, given the
 * path from `viterbi()`. Each member keeps the haplotype it had in the previous
 * subchain when its group has a sequence of it. Otherwise it takes the first
 * sequence in the group, except that when both members are in the same group, the
 * second one then avoids the sequence of the first if the group has another
 * sequence. The second member avoids the sequence of the first only when it does not
 * keep its own haplotype, so both members can take the same sequence.
 */
std::vector<std::pair<size_t, size_t>> choose_sequences(
    const std::vector<Candidates>& subchains, const std::vector<std::pair<size_t, size_t>>& path
);

//------------------------------------------------------------------------------

} // namespace vg::pair_model

#endif // VG_RECOMBINATOR_MODEL_HPP_INCLUDED
