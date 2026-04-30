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
 * Constructed from a reference alignment that defines the scoring scheme;
 * derives d at construction and exposes the marginal per-base match/mismatch
 * scores (in fractional precision) through which a MappingQualityCalculator
 * can be built. Can score the reference alignment (for free, via a cached
 * count) or any other alignment (which costs one path walk and applies the
 * formula with the reference's d).
 *
 * All knowledge of the formula's coefficients lives in this class.
 */
class LoggedGapAlignmentScorer : public AlignmentScorer {
public:
    /// Walks `reference` once, counts edits, derives d.
    /// `reference` must outlive any call to `score_alignment(reference)`
    /// since pointer identity is used to skip the second walk.
    explicit LoggedGapAlignmentScorer(const Alignment& reference);

    /// Score an alignment under this scheme. If `&aln == &reference_aln`,
    /// returns the cached value without re-walking the path.
    int32_t score_alignment(const Alignment& aln) const override;

    void fill_substitution_matrix(double out[16]) const override;

    /// Score a precomputed set of edit counts under this scheme.
    int32_t score_from_counts(size_t matches, size_t mismatches,
                              const std::vector<size_t>& gap_lengths) const;

    /// Walk an alignment and tally matches, mismatches, and per-gap lengths.
    /// Pure path traversal; no scheme parameter required.
    static void count_alignment_operations(const Alignment& aln,
                                           size_t& matches,
                                           size_t& mismatches,
                                           std::vector<size_t>& gap_lengths);

    /// Marginal per-base match score under this scheme. Always 1.0.
    double match = 1.0;
    /// Marginal per-base mismatch (and per-gap-open) score under this scheme.
    double mismatch = 0.0;

    /// Get the d parameter of the scheme (read-only). Public for diagnostics.
    double get_d() const { return d; }

private:
    /// max(0.02, (mm + go) / (m + mm + go)). Recovered at construction.
    double d;

    /// Identity and cached counts of the reference alignment, so
    /// `score_alignment(reference)` doesn't walk the path twice.
    const Alignment* reference_aln;
    size_t reference_matches;
    size_t reference_mismatches;
    std::vector<size_t> reference_gap_lengths;

    static double recover_d(size_t matches, size_t mismatches,
                            const std::vector<size_t>& gap_lengths);
};

} // namespace vg

#endif
