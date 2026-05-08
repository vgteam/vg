#ifndef VG_LOGGED_GAP_ALIGNMENT_SCORER_HPP_INCLUDED
#define VG_LOGGED_GAP_ALIGNMENT_SCORER_HPP_INCLUDED

#include <cstdint>
#include <vector>

#include <vg/vg.pb.h>

#include "alignment_scorer.hpp"

namespace vg {

/**
 * Scorer implementing the minimap2 long-indel rescoring formula
 * (Heng Li, Bioinformatics 2021, doi:10.1093/bioinformatics/btab705).
 *
 * Constructed from a "standard" alignment that defines the scoring scheme.
 * minimap2 only ever applies the scoring scheme back to that same alignment,
 * but we allow applying it to others.
 *
 * The case of applying the scoring scheme back to the standard alignment is
 * handled specially for speed; you MUST NOT apply the scorer to a different
 * alignmnet at the same address as the standard alignment. (An easy way to
 * ensure this is to ensure that the standard alignment outlives this object
 * and is not modified.)
 *
 * Hides everything about the empirical minimap2 formula inside the class.
 */
class LoggedGapAlignmentScorer : public AlignmentScorer {
public:
    /// Construct a LoggedGapAlignmentScorer using the scoring scheme defined
    /// by `standard`.
    ///
    /// `standard` must not be deallocated or modified between construction and
    /// any call to `score_alignment()` on any alignment at that address,
    /// because pointer identity is used to recognize that operations are on
    /// the standard alignment.
    explicit LoggedGapAlignmentScorer(const Alignment& standard);

    /// Score an alignment under this scheme. Guaranteed to be O(1) when aln is
    /// at the address of standard.
    int32_t score_alignment(const Alignment& aln) const override;

    /// Compute the "log base" that can be used to interpret scores
    /// probabilistically.
    double recover_log_base(double gc_content, double tol = 1e-12) const override;

    // Because all these fields are const, we can't just run a static member
    // that writes into each, so we need to do some backflips to fille
    // everything in with just initializer list expressions.

    // We expose the precomputed statistics about the standard alignment for
    // the user of the class to read off. These are used to initialize the
    // other members and MUST appear first here.

    /// Number of matches in the standard alignment
    const size_t matches;
    /// Number of mismatches in the standard alignment
    const size_t mismatches;
    /// Lengths of all gaps in the standard alignment
    const std::vector<size_t> gap_lengths;

private:

    /// Estimate of the divergence (what minimap2 calls "d").
    /// This is used to initialize mismatch, and depends on the
    /// statistics, and so MUST appear between them.
    const double divergence;

public:

    /// Marginal per-base match score under this scheme. Always 1.0.
    const double match = 1.0;
    /// Marginal per-base mismatch (and per-gap-open) score under this scheme.
    const double mismatch;

private:

    /// Build a LoggedGapAlignmentScorer with precomputed operation counts
    LoggedGapAlignmentScorer(const Alignment& standard, std::tuple<size_t, size_t, std::vector<size_t>>&& operation_counts);

    /// Score a precomputed set of edit counts under this scheme.
    int32_t score_from_counts(size_t matches, size_t mismatches,
                              const std::vector<size_t>& gap_lengths) const;

    /// Address where the standard alignment is, so we can recognize it later.
    const Alignment* standard_address;

    /// Walk an alignment and count matches, mismatches, and per-gap lengths.
    static void count_alignment_operations(const Alignment& aln,
                                           size_t& matches,
                                           size_t& mismatches,
                                           std::vector<size_t>& gap_lengths);
    
    /// Walk an alignment and count matches, mismatches, and per-gap lengths.
    /// Returns the result in a format usable to construct a LoggedGapAlignmentScorer
    static std::tuple<size_t, size_t, std::vector<size_t>> count_alignment_operations(const Alignment& aln);

    /// Helper function to compute the divergence estimate from edit information.
    static double compute_divergence(size_t matches, size_t mismatches,
                                     const std::vector<size_t>& gap_lengths);
};

} // namespace vg

#endif
