#ifndef VG_STRING_ALIGNER_HPP_INCLUDED
#define VG_STRING_ALIGNER_HPP_INCLUDED

/** \file
 * A global aligner of two strings that returns the operations of the alignment, and a function
 * that moves the gaps of an alignment as far left as they go.
 */

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace vg {

using namespace std;

/// The costs of a global alignment of two strings, and the limit on the work to find it.
///
/// A match costs nothing. A gap of `length` bases costs `gap_open + gap_extend * length`, or,
/// with the second piece, `long_gap_open + long_gap_extend * length` if that is less, so that a
/// long gap costs less than the short gaps and mismatches that could stand in for it.
struct StringAlignerParams {
    /// The cost of a mismatch.
    int32_t mismatch = 4;
    /// The cost of opening a gap.
    int32_t gap_open = 6;
    /// The cost of each base of a gap.
    int32_t gap_extend = 2;
    /// Whether a gap can take the second piece of the gap cost.
    bool long_gaps = true;
    /// The cost of opening a gap, in the second piece.
    int32_t long_gap_open = 24;
    /// The cost of each base of a gap, in the second piece.
    int32_t long_gap_extend = 1;
    /// The number of diagonals on each side of those between the two corners that the first band
    /// holds. Each wider band doubles it.
    size_t first_band = 16;
    /// The most cells a band that is traced back may hold. The aligner keeps one byte for each.
    size_t max_cells = (size_t)1 << 27;
    /// The most cells a band that is filled without tracing may hold. The aligner keeps only two
    /// rows of it, so this limits time, not memory.
    size_t max_work = (size_t)1 << 30;
    /// Whether the alignment is found with the wavefront algorithm, which needs no band and has
    /// none of the limits above: its time grows with the square of the alignment's cost and
    /// its memory with the cost, whatever the lengths of the strings.
    bool wavefront = false;
    /// With the wavefront algorithm and not negative: the alignment is given up, and nothing
    /// returned, once every alignment is known to cost more than this. An alignment that
    /// costs a little more can still be returned.
    int64_t max_cost = -1;

    /// Return the cost of a gap of `length` bases.
    int64_t gap_cost(size_t length) const;
};

/// What an alignment does with a run of bases.
enum class AlignmentOp : uint8_t {
    /// The bases of the two strings are the same.
    Match,
    /// Each base of the first string is replaced by a different base of the second.
    Mismatch,
    /// The bases are in the second string alone.
    Insertion,
    /// The bases are in the first string alone.
    Deletion
};

/// A run of one operation of an alignment.
struct AlignmentRun {
    AlignmentOp op;
    /// The number of bases, which is not 0.
    size_t length;

    bool operator==(const AlignmentRun& other) const = default;
};

/// A global alignment of two strings.
struct StringAlignment {
    /// The cost of the alignment.
    int64_t cost = 0;
    /// The runs of the alignment, from the start of both strings to their ends. No two neighbours
    /// have the same operation.
    vector<AlignmentRun> runs;
    /// The number of diagonals on each side of those between the two corners that the last band
    /// held, or 0 if no band was needed.
    size_t band = 0;
};

/**
 * Align all of `first` to all of `second` at the least cost.
 *
 * The aligner fills a band of diagonals around those between the two corners of the alignment
 * matrix, keeping one byte for each cell, so its memory is the length of `first` times the number
 * of diagonals in the band. An alignment that leaves the band has a gap that reaches past the
 * band's edge and a gap that comes back, and the two cost at least a known amount. If the best
 * alignment in the band costs no more than that, no alignment costs less, and it is returned.
 * Otherwise the aligner doubles the width of the band and fills it again.
 *
 * Where doubling the band finds nothing cheaper, and the cost is still more than leaving the
 * band could cost, the aligner finds the least cost another way. It fills, keeping two rows and
 * no trace, the narrowest band that an alignment of that cost cannot leave. It then returns the
 * first alignment of the least cost that a traced band holds.
 *
 * Bases are compared as characters.
 *
 * @return Nothing if the alignment of least cost is not found within `params.max_cells` and
 * `params.max_work`.
 */
optional<StringAlignment> align_global(const string& first, const string& second,
                                       const StringAlignerParams& params = StringAlignerParams());

/**
 * Align all of `first` to all of `second` at the least cost among the alignments that stay
 * inside a region guided by `corners`.
 *
 * Each corner is a number of bases of `first` and of `second` that the alignment is expected to
 * pass near, in order. The first cell and the last are corners too, and a corner is moved forward
 * where it would go back from the one before. The region is, for each two corners that follow
 * each other, the band of the rows from the one to the other and of the diagonals from the one's
 * to the other's, each widened by `margin`. So a stretch whose two strings have about the same
 * length is a band about one diagonal, and only a stretch with a long gap is wide. As in
 * `align_global`, the bases the two share at their starts and at their ends are matched first,
 * and the corners are taken in what lies between. The aligner keeps a byte for each cell of the
 * region.
 *
 * @return Nothing if the region has more cells than `params.max_cells`.
 */
optional<StringAlignment> align_guided(const string& first, const string& second,
                                       const vector<pair<size_t, size_t>>& corners,
                                       size_t margin, const StringAlignerParams& params);

/**
 * Move each gap of `runs`, an alignment of `first` to `second`, as far left as it goes.
 *
 * A gap moves one base left where the base before it is matched and is the same as the last base
 * of the gap. The alignment spells the same two strings and has the same numbers of each
 * operation, so two alignments that differ only in where a gap sits in a repeat come out the
 * same. Gaps of one kind that come to touch are joined.
 *
 * @param keep_first Leave the first column of the alignment in place, so that a gap that could
 * move to the start of the strings keeps one base before it.
 */
void left_align_gaps(const string& first, const string& second, vector<AlignmentRun>& runs,
                     bool keep_first = false);

}

#endif
