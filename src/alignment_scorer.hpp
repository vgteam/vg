#ifndef VG_ALIGNMENT_SCORER_HPP_INCLUDED
#define VG_ALIGNMENT_SCORER_HPP_INCLUDED

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include <vg/vg.pb.h>

#include "handle.hpp"
#include "path.hpp"
#include "types.hpp"

namespace vg {

/// Defaults previously declared in aligner.hpp; duplicated here so callers can
/// build a scorer without pulling in the GSSW machinery.
static constexpr int8_t default_scorer_match = 1;
static constexpr int8_t default_scorer_mismatch = 4;
static constexpr int8_t default_scorer_score_matrix[16] = {
     default_scorer_match,    -default_scorer_mismatch, -default_scorer_mismatch, -default_scorer_mismatch,
    -default_scorer_mismatch,  default_scorer_match,    -default_scorer_mismatch, -default_scorer_mismatch,
    -default_scorer_mismatch, -default_scorer_mismatch,  default_scorer_match,    -default_scorer_mismatch,
    -default_scorer_mismatch, -default_scorer_mismatch, -default_scorer_mismatch,  default_scorer_match
};
static constexpr int8_t default_scorer_gap_open = 6;
static constexpr int8_t default_scorer_gap_extension = 1;
static constexpr int8_t default_scorer_full_length_bonus = 5;

/// Score a gap with the given open and extension scores. Lives here (rather
/// than aligner.hpp) so AlignmentScorer can use it without a circular include.
int32_t score_gap(size_t gap_length, int32_t gap_open, int32_t gap_extension);

/**
 * Root contract every scorer satisfies. Just enough virtual surface for
 * MappingQualityCalculator to compute log_base from any scorer.
 */
class AlignmentScorer {
public:
    virtual ~AlignmentScorer() = default;

    /// Compute a single integer score for the alignment under this scorer's scheme.
    virtual int32_t score_alignment(const Alignment& aln) const = 0;

    /// Fill a 4x4 substitution matrix in double precision. Used at MQ calc
    /// construction time to recover log_base.
    virtual void fill_substitution_matrix(double out[16]) const = 0;

    /// Bisects to find the log base under which the partition function over
    /// the given 4x4 substitution matrix and inferred (gc-content-driven)
    /// nucleotide frequencies equals 1. Static so that scorers needing
    /// log_base internally (QualAdjAlignmentScorer, while building its
    /// quality-adjusted tables) can reuse the same routine the
    /// MappingQualityCalculator does.
    static double recover_log_base(const double matrix[16], double gc_content, double tol = 1e-12);

private:
    static bool verify_valid_log_odds_score_matrix(const double matrix[16], const double nt_freqs[4]);
    static double alignment_score_partition_function(double lambda, const double matrix[16], const double nt_freqs[4]);
};

/**
 * Scorer surface that GSSW DP and surjection rescoring need: per-edit and
 * per-partial-path scoring functions. Matrix-style scorers implement this;
 * the logged-gap scorer does not.
 */
class EditAlignmentScorer : public AlignmentScorer {
public:
    /// Compute the score of an exact match in the given alignment, from the
    /// given offset, of the given length.
    virtual int32_t score_exact_match(const Alignment& aln, size_t read_offset, size_t length) const = 0;
    /// Compute the score of an exact match of the given sequence with the given qualities.
    /// Qualities may be ignored by some implementations.
    virtual int32_t score_exact_match(const std::string& sequence, const std::string& base_quality) const = 0;
    /// Compute the score of an exact match of the given range of sequence with the given qualities.
    /// Qualities may be ignored by some implementations.
    virtual int32_t score_exact_match(std::string::const_iterator seq_begin, std::string::const_iterator seq_end,
                                      std::string::const_iterator base_qual_begin) const = 0;
    /// Compute the score of a mismatch of the given range of sequence with the given qualities.
    /// Return is signed and almost certainly negative.
    virtual int32_t score_mismatch(std::string::const_iterator seq_begin, std::string::const_iterator seq_end,
                                   std::string::const_iterator base_qual_begin) const = 0;

    virtual int32_t score_full_length_bonus(bool left_side, std::string::const_iterator seq_begin,
                                            std::string::const_iterator seq_end,
                                            std::string::const_iterator base_qual_begin) const = 0;
    virtual int32_t score_full_length_bonus(bool left_side, const Alignment& alignment) const = 0;

    /// Compute the score of a path against the given range of subsequence with the given qualities.
    virtual int32_t score_partial_alignment(const Alignment& alignment, const HandleGraph& graph,
                                            const path_t& path,
                                            std::string::const_iterator seq_begin,
                                            bool no_read_end_scoring = false) const = 0;

    /// Score a gap with this scorer's gap_open / gap_extension parameters.
    virtual int32_t score_gap(size_t gap_length) const;

    /// Use the score values in the scorer to score the given alignment,
    /// scoring gaps caused by jumping between nodes using a custom gap length
    /// estimation function.
    virtual int32_t score_discontiguous_alignment(const Alignment& aln,
        const std::function<size_t(pos_t, pos_t, size_t)>& estimate_distance,
        bool allow_left_bonus = true,
        bool allow_right_bonus = true) const;

    /// Score the given alignment assuming there are no gaps between Mappings.
    virtual int32_t score_contiguous_alignment(const Alignment& aln,
                                               bool allow_left_bonus = true,
                                               bool allow_right_bonus = true) const;

    /// Return the score with full-length bonuses removed.
    virtual int32_t remove_bonuses(const Alignment& aln, bool pinned = false, bool pin_left = false) const;

    /// Default scorer-level score_alignment is a contiguous-alignment score
    /// with both bonuses included.
    int32_t score_alignment(const Alignment& aln) const override;

    /// The longest gap detectable from a read position without soft-clipping.
    size_t longest_detectable_gap(const Alignment& alignment, const std::string::const_iterator& read_pos) const;
    size_t longest_detectable_gap(size_t read_length, size_t read_pos) const;
    size_t longest_detectable_gap(const Alignment& alignment) const;
    size_t longest_detectable_gap(size_t read_length) const;

    // Direct-access score parameters. Concrete subclasses fill these in.
    int8_t  match = 0;
    int8_t  mismatch = 0;
    int8_t  gap_open = 0;
    int8_t  gap_extension = 0;
    int8_t  full_length_bonus = 0;
};

/**
 * Concrete scorer that holds an int8_t score matrix and ignores read base
 * qualities. Owns its score matrix and nt_table; their layout matches what
 * GSSW expects (5x5, N-padded).
 */
class MatrixAlignmentScorer : public EditAlignmentScorer {
public:
    /// Build from a 4x4 substitution matrix and gap parameters. Allocates
    /// owned 5x5 N-padded `score_matrix` and `nt_table`.
    MatrixAlignmentScorer(const int8_t* score_matrix_4x4 = default_scorer_score_matrix,
                          int8_t gap_open = default_scorer_gap_open,
                          int8_t gap_extension = default_scorer_gap_extension,
                          int8_t full_length_bonus = default_scorer_full_length_bonus);

    ~MatrixAlignmentScorer() override;

    // Non-copyable, non-movable to keep ownership semantics simple.
    MatrixAlignmentScorer(const MatrixAlignmentScorer&) = delete;
    MatrixAlignmentScorer& operator=(const MatrixAlignmentScorer&) = delete;

    int32_t score_exact_match(const Alignment& aln, size_t read_offset, size_t length) const override;
    int32_t score_exact_match(const std::string& sequence, const std::string& base_quality) const override;
    int32_t score_exact_match(std::string::const_iterator seq_begin, std::string::const_iterator seq_end,
                              std::string::const_iterator base_qual_begin) const override;
    /// Length-only convenience for callers that lack a sequence.
    int32_t score_exact_match(const std::string& sequence) const;
    int32_t score_exact_match(std::string::const_iterator seq_begin, std::string::const_iterator seq_end) const;

    int32_t score_mismatch(std::string::const_iterator seq_begin, std::string::const_iterator seq_end,
                           std::string::const_iterator base_qual_begin) const override;
    /// Length-only mismatch score (returns a signed, usually negative value).
    int32_t score_mismatch(size_t length) const;

    int32_t score_full_length_bonus(bool left_side, std::string::const_iterator seq_begin,
                                    std::string::const_iterator seq_end,
                                    std::string::const_iterator base_qual_begin) const override;
    int32_t score_full_length_bonus(bool left_side, const Alignment& alignment) const override;

    int32_t score_partial_alignment(const Alignment& alignment, const HandleGraph& graph,
                                    const path_t& path,
                                    std::string::const_iterator seq_begin,
                                    bool no_read_end_scoring = false) const override;

    void fill_substitution_matrix(double out[16]) const override;

    /// 5x5 N-padded matrix in row-major order. Owned. Consumed directly by
    /// GSSW and by score_partial_alignment in the qual-adjusted subclass.
    int8_t* score_matrix = nullptr;
    /// Char->int nt translation. Owned.
    int8_t* nt_table = nullptr;
    /// Original 4x4 substitution matrix in row-major order. Owned. Used to
    /// compute log_base so the qual-adjusted subclass can keep a separate
    /// (overwritten) score_matrix while still exposing the un-quality-adjusted
    /// substitution scheme for MAPQ.
    int8_t* substitution_matrix = nullptr;
};

/**
 * Quality-adjusted matrix scorer. Owns a quality-indexed score matrix and a
 * per-quality full-length bonus table, and overrides the quality-aware
 * scoring functions accordingly. The int8_t scalar fields inherited from
 * MatrixAlignmentScorer (match, mismatch, ...) keep their original
 * un-quality-adjusted values for log_base recovery and identity correction;
 * they are not consulted by this class's own scoring routines.
 */
class QualAdjAlignmentScorer : public MatrixAlignmentScorer {
public:
    QualAdjAlignmentScorer(const int8_t* score_matrix_4x4 = default_scorer_score_matrix,
                           int8_t gap_open = default_scorer_gap_open,
                           int8_t gap_extension = default_scorer_gap_extension,
                           int8_t full_length_bonus = default_scorer_full_length_bonus,
                           double gc_content_for_qual_adj = 0.5);

    ~QualAdjAlignmentScorer() override;

    int32_t score_exact_match(const Alignment& aln, size_t read_offset, size_t length) const override;
    int32_t score_exact_match(const std::string& sequence, const std::string& base_quality) const override;
    int32_t score_exact_match(std::string::const_iterator seq_begin, std::string::const_iterator seq_end,
                              std::string::const_iterator base_qual_begin) const override;
    int32_t score_mismatch(std::string::const_iterator seq_begin, std::string::const_iterator seq_end,
                           std::string::const_iterator base_qual_begin) const override;
    int32_t score_full_length_bonus(bool left_side, std::string::const_iterator seq_begin,
                                    std::string::const_iterator seq_end,
                                    std::string::const_iterator base_qual_begin) const override;
    int32_t score_full_length_bonus(bool left_side, const Alignment& alignment) const override;
    int32_t score_partial_alignment(const Alignment& alignment, const HandleGraph& graph,
                                    const path_t& path,
                                    std::string::const_iterator seq_begin,
                                    bool no_read_end_scoring = false) const override;

    /// Per-quality full-length bonus table. Owned.
    int8_t* qual_adj_full_length_bonuses = nullptr;

protected:
    int8_t* qual_adjusted_matrix(const int8_t* score_matrix_4x4, double gc_content,
                                 double log_base, uint32_t max_qual) const;
    int8_t* qual_adjusted_bonuses(int8_t base_full_length_bonus, double log_base,
                                  uint32_t max_qual) const;
};

} // namespace vg

#endif
