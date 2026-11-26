#ifndef VG_READ_CONSENSUS_HPP_INCLUDED
#define VG_READ_CONSENSUS_HPP_INCLUDED

/** \file
 * The consensus of reads that span one window, built by aligning each read to a backbone.
 */

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "string_aligner.hpp"

namespace vg {

using namespace std;

/// How a consensus of reads is built.
struct ReadConsensusParams {
    /// The most reads the consensus is built from. With more, those whose lengths are nearest
    /// the median length are used.
    size_t max_reads = 12;
    /// The fewest reads a consensus is built from.
    size_t min_reads = 3;
    /// The most times the reads are aligned to a backbone.
    static constexpr size_t max_rounds = 4;
    /// The most inserted strings that are compared with each other to pick the one an insertion
    /// takes. With more, those whose lengths are nearest the median length are compared.
    static constexpr size_t max_medoid = 8;
    /// Whether the consensus starts from the read nearest the others, rather than from the
    /// backbone, where the backbone is further from the reads than that read is.
    bool seed_from_reads = true;
    /// Whether the consensus keeps the backbone's length of a run of one base wherever it
    /// differs from the backbone only in the length of that run.
    bool keep_homopolymer_lengths = true;
    /// The fewest bases a run of one base has in the backbone for the consensus to keep its
    /// length. A shorter run takes the length the reads give it.
    size_t homopolymer_min_run = 9;
    /// The costs and limits of each alignment of a read to a backbone.
    StringAlignerParams aligner{.wavefront = true};
};

/**
 * Return the consensus of `reads`, which all span one window from end to end in one orientation,
 * given `backbone`, the sequence now called across the window.
 *
 * Each round aligns every read to the current backbone and moves the gaps left. Each base of
 * the backbone becomes whichever of a base or a deletion most reads have there, and stays as
 * it is on a tie. Between two bases of the backbone an insertion is made where more than half
 * the reads insert bases, and its bases are those of the medoid of the strings the reads
 * insert: the one with the least summed edit distance to the others. The result is the next
 * round's backbone, until a round changes nothing or `params.max_rounds` rounds are done.
 *
 * The first round is from `backbone`. `backbone` is returned where that round changes nothing,
 * or only what `params.keep_homopolymer_lengths` undoes, and `backbone` is no shorter than the
 * shortest of the reads used and no longer than the longest. Otherwise, with
 * `params.seed_from_reads`, the rounds start again from the medoid read where it is nearer
 * the reads than `backbone` is.
 *
 * The first round stops aligning reads once the rest cannot change the backbone, where the
 * backbone lies within the lengths of the reads (see `consensus_round`), so the backbone can be
 * returned where a read it did not align has no alignment within the limits of
 * `params.aligner`.
 *
 * @return Nothing if there are fewer than `params.min_reads` reads, or an alignment is not found
 * within the limits of `params.aligner`.
 */
optional<string> read_consensus(const vector<string>& reads, const string& backbone,
                                const ReadConsensusParams& params = ReadConsensusParams());

/**
 * Return one round of `read_consensus`: `backbone` with each base and each gap between bases set
 * to what most of `reads` have there.
 *
 * With `settle`, return `backbone` as soon as the reads aligned so far show that the rest cannot
 * change it: each base of the backbone leads every other vote by at least the number of reads
 * left, and no gap has inserted bases in so many reads that the rest could make them more than
 * half. That is the round's result where every read aligns. Where one of the rest has no
 * alignment within the limits of `params.aligner`, the full round gives nothing, and the round
 * with `settle` can still give `backbone`.
 *
 * @return Nothing if an alignment it makes is not found within the limits of `params.aligner`.
 */
optional<string> consensus_round(const vector<string>& reads, const string& backbone,
                                 const ReadConsensusParams& params, bool settle);

/**
 * Return `sequence` with each difference from `backbone` that only changes the length of a run of
 * one base undone.
 *
 * The two are aligned and the gaps moved left. A gap whose bases are all one base, beside a
 * matched base of the backbone that is the same, changes only the length of a run, and the result
 * has the backbone's bases there if the backbone's run has at least `min_run` bases.
 *
 * @return Nothing if the alignment is not found within the limits of `aligner`.
 */
optional<string> keep_homopolymer_lengths(const string& sequence, const string& backbone,
                                          const StringAlignerParams& aligner =
                                              StringAlignerParams(),
                                          size_t min_run = 1);

/// Return the least number of single-base substitutions, insertions and deletions that turn
/// `first` into `second`, or nothing if it is not found.
optional<int64_t> edit_distance(const string& first, const string& second);

}

#endif
