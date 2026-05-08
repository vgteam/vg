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

/// Default scoring parameters.
static constexpr int8_t default_match = 1;
static constexpr int8_t default_mismatch = 4;
static constexpr int8_t default_score_matrix[16] = {
     default_match,    -default_mismatch, -default_mismatch, -default_mismatch,
    -default_mismatch,  default_match,    -default_mismatch, -default_mismatch,
    -default_mismatch, -default_mismatch,  default_match,    -default_mismatch,
    -default_mismatch, -default_mismatch, -default_mismatch,  default_match
};
static constexpr int8_t default_gap_open = 6;
static constexpr int8_t default_gap_extension = 1;
static constexpr int8_t default_full_length_bonus = 5;
static constexpr double default_gc_content = 0.5;

/// Score a gap with the given open and extension scores.
int32_t score_gap(size_t gap_length, int32_t gap_open, int32_t gap_extension);

/**
 * Base interface for computing scores for alignments under some (not
 * necessarily affine-gap) scoring scheme.
 */
class AlignmentScorer {
public:
    virtual ~AlignmentScorer() = default;

    /// Compute a single integer score for the alignment under this scorer's scheme.
    virtual int32_t score_alignment(const Alignment& aln) const = 0;

    /// Get the log base value used to probabilistically interpret scores.
    /// Implementations probably need to have a stores GC content to implement this.
    virtual double get_log_base() const = 0;

protected:

    /// Bisects to find the log base under which the partition function over
    /// the given 4x4 substitution matrix and nucleotide frequencies (from GC
    /// content) equals 1.
    ///
    /// Useful for implementing the recover_log_base operation for a derived
    /// class, if that class can produce a scoring matrix.
    static double recover_log_base(const double matrix[16], double gc_content, double tol = 1e-12);

private:
    /// Check a score matrix to make sure it has negative scores for random sequence.
    static bool verify_valid_log_odds_score_matrix(const double matrix[16], const double nt_freqs[4]);
    /// Partition function used for recovering the log base for a scoring scheme.
    /// TODO: What is a partition function exactly?
    static double alignment_score_partition_function(double lambda, const double matrix[16], const double nt_freqs[4]);
};

/**
 * An abstract AlignmentScorer that can score individual edits independently.
 */
class EditAlignmentScorer : public AlignmentScorer {
public:

    EditAlignmentScorer(int8_t match = default_match,
                        int8_t mismatch = default_mismatch,
                        int8_t gap_open = default_gap_open,
                        int8_t gap_extension = default_gap_extension,
                        int8_t full_length_bonus = default_full_length_bonus);
   
    ////
    // Implement AlignmentScorer
    ////

    /// Score an alignment as a contiguous alignment, including full-length bonuses.
    int32_t score_alignment(const Alignment& aln) const override;
    
    ////
    // Provide extra tools for working with pieces of alignments
    ////

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

    /// The longest gap detectable from a read position without soft-clipping.
    size_t longest_detectable_gap(const Alignment& alignment, const std::string::const_iterator& read_pos) const;
    size_t longest_detectable_gap(size_t read_length, size_t read_pos) const;
    size_t longest_detectable_gap(const Alignment& alignment) const;
    size_t longest_detectable_gap(size_t read_length) const;

    ////
    // Expose individual generic operation scores.
    ////

    int8_t match;
    int8_t mismatch;
    int8_t gap_open;
    int8_t gap_extension;
    int8_t full_length_bonus;
};

/**
 * An alignment scorer that scores bassed on a score matrix.
 *
 * Takes 4x4 matrices for construction, but internally uses 5x5 N-padded,
 * GSSW-compatible matrices.
 */
class MatrixAlignmentScorer : public EditAlignmentScorer {
public:
    /// Build from a 4x4 substitution matrix and gap parameters.
    /// 
    /// Fills in the generic base match and mismatch values from the matrix.
    MatrixAlignmentScorer(const int8_t* score_matrix_4x4 = default_score_matrix,
                          int8_t gap_open = default_gap_open,
                          int8_t gap_extension = default_gap_extension,
                          int8_t full_length_bonus = default_full_length_bonus,
                          double gc_content = default_gc_content);

    ~MatrixAlignmentScorer() override;

    // Non-copyable, non-movable to keep ownership semantics simple.
    MatrixAlignmentScorer(const MatrixAlignmentScorer&) = delete;
    MatrixAlignmentScorer& operator=(const MatrixAlignmentScorer&) = delete;

    /// Retriece the log base value that was computed at construction.
    double get_log_base() const override;

    int32_t score_exact_match(const Alignment& aln, size_t read_offset, size_t length) const override;
    int32_t score_exact_match(const std::string& sequence, const std::string& base_quality) const override;
    int32_t score_exact_match(std::string::const_iterator seq_begin, std::string::const_iterator seq_end,
                              std::string::const_iterator base_qual_begin) const override;
    /// Score an exact match of unspecified sequence and the given length.
    int32_t score_exact_match(const std::string& sequence) const;
    int32_t score_exact_match(std::string::const_iterator seq_begin, std::string::const_iterator seq_end) const;

    int32_t score_mismatch(std::string::const_iterator seq_begin, std::string::const_iterator seq_end,
                           std::string::const_iterator base_qual_begin) const override;
    /// Score a mismatch of unspecified sequence of the given length.
    /// This usually produces a negative score.
    int32_t score_mismatch(size_t length) const;

    int32_t score_full_length_bonus(bool left_side, std::string::const_iterator seq_begin,
                                    std::string::const_iterator seq_end,
                                    std::string::const_iterator base_qual_begin) const override;
    int32_t score_full_length_bonus(bool left_side, const Alignment& alignment) const override;

    int32_t score_partial_alignment(const Alignment& alignment, const HandleGraph& graph,
                                    const path_t& path,
                                    std::string::const_iterator seq_begin,
                                    bool no_read_end_scoring = false) const override;

    /// 5x5 GSSW-style N-padded matrix in row-major order, possibly with multiple levels.
    /// Owned by this class.
    ///
    /// This class makes this a 5x5x1 matrix, but subclasses can make this have
    /// more dimensions.
    int8_t* score_matrix;
    /// Char->int nt translation. Owned by this class.
    int8_t* nt_table;

protected:
    /// We remember our computed log base.
    double log_base;

};

/**
 * AlignmentScorer that scores alignments using quality-adjusted scoring.
 *
 * Handles adjusting a 4x4 input scoring matrix for base quality.
 */
class QualAdjAlignmentScorer : public MatrixAlignmentScorer {
public:
    QualAdjAlignmentScorer(const int8_t* score_matrix_4x4 = default_score_matrix,
                           int8_t gap_open = default_gap_open,
                           int8_t gap_extension = default_gap_extension,
                           int8_t full_length_bonus = default_full_length_bonus,
                           double gc_content = default_gc_content);

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

    /// Per-quality full-length bonus table. Owned by this class.
    int8_t* qual_adj_full_length_bonuses = nullptr;

protected:
    int8_t* qual_adjusted_matrix(const int8_t* score_matrix_4x4, double gc_content,
                                 double log_base, uint32_t max_qual) const;
    int8_t* qual_adjusted_bonuses(int8_t base_full_length_bonus, double log_base,
                                  uint32_t max_qual) const;
};

} // namespace vg

#endif
