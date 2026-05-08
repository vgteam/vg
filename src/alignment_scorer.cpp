#include "alignment_scorer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "alignment.hpp"
#include "crash.hpp"
#include "gssw.h"

namespace vg {

using std::function;
using std::max;
using std::min;
using std::string;
using std::vector;

int32_t score_gap(size_t gap_length, int32_t gap_open, int32_t gap_extension) {
    return gap_length ? -(gap_open + (gap_length - 1) * gap_extension) : 0;
}

// ----- AlignmentScorer -----

double AlignmentScorer::recover_log_base(const double matrix[16], double gc_content, double tol) {
    // convert gc content into base-wise frequencies
    double nt_freqs[4];
    nt_freqs[0] = 0.5 * (1 - gc_content);
    nt_freqs[1] = 0.5 * gc_content;
    nt_freqs[2] = 0.5 * gc_content;
    nt_freqs[3] = 0.5 * (1 - gc_content);

    if (!verify_valid_log_odds_score_matrix(matrix, nt_freqs)) {
        // Dump the offending matrix
        for (size_t x = 0; x < 4; x++) {
            std::cerr << "error:[AlignmentScorer]";
            for (size_t y = 0; y < 4; y++) {
                std::cerr << " " << matrix[y * 4 + x];
            }
            std::cerr << std::endl;
        }
        std::cerr << "error:[AlignmentScorer] Score matrix is invalid. Must have a negative expected score against random sequence." << std::endl;
        // TODO: Use new logging stuff.
        std::exit(1);
    }

    // searching for a positive value (because it's a base of a logarithm)
    double lower_bound;
    double upper_bound;

    // arbitrary starting point greater than zero
    double lambda = 1.0;
    // exponential search for a window containing lambda where total probability is 1
    double partition = alignment_score_partition_function(lambda, matrix, nt_freqs);
    if (partition < 1.0) {
        lower_bound = lambda;
        while (partition <= 1.0) {
            lower_bound = lambda;
            lambda *= 2.0;
            partition = alignment_score_partition_function(lambda, matrix, nt_freqs);
        }
        upper_bound = lambda;
    } else {
        upper_bound = lambda;
        while (partition >= 1.0) {
            upper_bound = lambda;
            lambda /= 2.0;
            partition = alignment_score_partition_function(lambda, matrix, nt_freqs);
        }
        lower_bound = lambda;
    }

    // bisect to find a log base where total probability is 1
    while (upper_bound / lower_bound - 1.0 > tol) {
        lambda = 0.5 * (lower_bound + upper_bound);
        if (alignment_score_partition_function(lambda, matrix, nt_freqs) < 1.0) {
            lower_bound = lambda;
        } else {
            upper_bound = lambda;
        }
    }
    return 0.5 * (lower_bound + upper_bound);
}

bool AlignmentScorer::verify_valid_log_odds_score_matrix(const double matrix[16], const double nt_freqs[4]) {
    bool contains_positive_score = false;
    for (int i = 0; i < 16; ++i) {
        if (matrix[i] > 0) {
            contains_positive_score = true;
            break;
        }
    }
    if (!contains_positive_score) return false;

    double expected_score = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            expected_score += nt_freqs[i] * nt_freqs[j] * matrix[i * 4 + j];
        }
    }
    return expected_score < 0.0;
}

double AlignmentScorer::alignment_score_partition_function(double lambda, const double matrix[16], const double nt_freqs[4]) {
    double partition = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            partition += nt_freqs[i] * nt_freqs[j] * std::exp(lambda * matrix[i * 4 + j]);
        }
    }
    if (std::isnan(partition)) {
        std::cerr << "error:[AlignmentScorer] overflow in log-odds base recovery." << std::endl;
        std::exit(1);
    }
    return partition;
}

// ----- EditAlignmentScorer -----

EditAlignmentScorer::EditAlignmentScorer(
    int8_t match,
    int8_t mismatch,
    int8_t gap_open,
    int8_t gap_extension,
    int8_t full_length_bonus
) : AlignmentScorer(),
    match(match), mismatch(mismatch),
    gap_open(gap_open), gap_extension(gap_extension),
    full_length_bonus(full_length_bonus)
{
    // Nothing to do!
}

int32_t EditAlignmentScorer::score_gap(size_t gap_length) const {
    return vg::score_gap(gap_length, gap_open, gap_extension);
}

int32_t EditAlignmentScorer::score_alignment(const Alignment& aln) const {
    return score_contiguous_alignment(aln, true, true);
}

int32_t EditAlignmentScorer::score_discontiguous_alignment(const Alignment& aln,
    const function<size_t(pos_t, pos_t, size_t)>& estimate_distance,
    bool allow_left_bonus, bool allow_right_bonus) const {

    int score = 0;
    int read_offset = 0;
    auto& path = aln.path();

    // We keep track of whether the last edit was a deletion for coalescing
    // adjacent deletions across node boundaries
    bool last_was_deletion = false;

    for (int i = 0; i < path.mapping_size(); ++i) {
        // For each mapping
        auto& mapping = path.mapping(i);
        for (int j = 0; j < mapping.edit_size(); ++j) {
            // For each edit in the mapping
            auto& edit = mapping.edit(j);

            // Score the edit according to its type
            if (edit_is_match(edit)) {
                score += score_exact_match(aln, read_offset, edit.to_length());
                last_was_deletion = false;
            } else if (edit_is_sub(edit)) {
                score += score_mismatch(aln.sequence().begin() + read_offset,
                                        aln.sequence().begin() + read_offset + edit.to_length(),
                                        aln.quality().begin() + read_offset);
                last_was_deletion = false;
            } else if (edit_is_deletion(edit)) {
                if (last_was_deletion) {
                    // No need to charge a gap open
                    score -= edit.from_length() * gap_extension;
                } else {
                    // We need a gap open
                    score -= edit.from_length() ? gap_open + (edit.from_length() - 1) * gap_extension : 0;
                }

                if (edit.from_length()) {
                    // We already charged a gap open
                    last_was_deletion = true;
                }
                // If there's a 0-length deletion, leave the last_was_deletion flag unchanged.
            } else if (edit_is_insertion(edit) && !((i == 0 && j == 0) ||
                                                    (i == path.mapping_size() - 1 && j == mapping.edit_size() - 1))) {
                // todo how do we score this qual adjusted?
                score -= edit.to_length() ? gap_open + (edit.to_length() - 1) * gap_extension : 0;
                last_was_deletion = false;
                // No need to track if the last edit was an insertion because
                // insertions will be all together in a single edit at a point.
            } else {
                // Edit has no score effect. Probably a softclip.
                last_was_deletion = false;
            }
            read_offset += edit.to_length();
        }
        // score any intervening gaps in mappings using approximate distances
        if (i + 1 < path.mapping_size()) {
            // what is the distance between the last position of this mapping
            // and the first of the next
            Position last_pos = mapping.position();
            last_pos.set_offset(last_pos.offset() + mapping_from_length(mapping));
            Position next_pos = path.mapping(i + 1).position();
            // Estimate the distance
            int dist = estimate_distance(make_pos_t(last_pos), make_pos_t(next_pos), aln.sequence().size());
            if (dist > 0) {
                // If it's nonzero, score it as a deletion gap
                score -= gap_open + (dist - 1) * gap_extension;
            }
        }
    }

    // We should report any bonuses used in the DP in the final score
    if (allow_left_bonus && !softclip_start(aln)) {
        score += score_full_length_bonus(true, aln);
    }
    if (allow_right_bonus && !softclip_end(aln)) {
        score += score_full_length_bonus(false, aln);
    }

    return score;
}

int32_t EditAlignmentScorer::score_contiguous_alignment(const Alignment& aln,
                                                        bool allow_left_bonus,
                                                        bool allow_right_bonus) const {
    return score_discontiguous_alignment(aln, [](pos_t, pos_t, size_t) { return (size_t) 0; },
                                         allow_left_bonus, allow_right_bonus);
}

int32_t EditAlignmentScorer::remove_bonuses(const Alignment& aln, bool pinned, bool pin_left) const {
    int32_t score = aln.score();
    if (softclip_start(aln) == 0 && !(pinned && pin_left)) {
        // No softclip at the start, and a left end bonus was applied.
        score -= score_full_length_bonus(true, aln);
    }
    if (softclip_end(aln) == 0 && !(pinned && !pin_left)) {
        // No softclip at the end, and a right end bonus was applied.
        score -= score_full_length_bonus(false, aln);
    }
    return score;
}

size_t EditAlignmentScorer::longest_detectable_gap(const Alignment& alignment, const string::const_iterator& read_pos) const {
    return longest_detectable_gap(alignment.sequence().size(), read_pos - alignment.sequence().begin());
}

size_t EditAlignmentScorer::longest_detectable_gap(size_t read_length, size_t read_pos) const {
    // algebraic solution for when score is > 0 assuming perfect match other than gap
    assert(read_length >= read_pos);
    int64_t overhang_length = min(read_pos, read_length - read_pos);
    int64_t numer = match * overhang_length + full_length_bonus;
    int64_t gap_length = (numer - gap_open) / gap_extension + 1;
    return gap_length >= 0 && overhang_length > 0 ? gap_length : 0;
}

size_t EditAlignmentScorer::longest_detectable_gap(const Alignment& alignment) const {
    // longest detectable gap across entire read is in the middle
    return longest_detectable_gap(alignment.sequence().size(), alignment.sequence().size() / 2);
}

size_t EditAlignmentScorer::longest_detectable_gap(size_t read_length) const {
    return longest_detectable_gap(read_length, read_length / 2);
}

// ----- MatrixAlignmentScorer -----

MatrixAlignmentScorer::MatrixAlignmentScorer(
    const int8_t* score_matrix_4x4,
    int8_t gap_open,
    int8_t gap_extension,
    int8_t full_length_bonus,
    double gc_content
) : EditAlignmentScorer(score_matrix_4x4[0], -score_matrix_4x4[1], gap_open, gap_extension, full_length_bonus),
    score_matrix((int8_t*) std::malloc(sizeof(int8_t) * 25)),
    nt_table(gssw_create_nt_table())
{
    // Fill in the 5x5 score matrix, which adds in the 5th row and column of 0s
    // for N matches like GSSW wants
    crash_unless(score_matrix != nullptr);
    for (size_t i = 0, j = 0; i < 25; ++i) {
        if (i % 5 == 4 || i / 5 == 4) {
            score_matrix[i] = 0;
        } else {
            score_matrix[i] = score_matrix_4x4[j];
            ++j;
        }
    }

    // Also make a double-type score matrix and get the log base. We can't
    // easily do this in an initializer.
    double double_matrix[16];
    for (int i = 0; i < 16; ++i) {
        double_matrix[i] = static_cast<double>(score_matrix_4x4[i]);
    }
    // Recover a log base from that
    log_base = AlignmentScorer::recover_log_base(double_matrix, gc_content);
}

MatrixAlignmentScorer::~MatrixAlignmentScorer() {
    if (nt_table) std::free(nt_table);
    if (score_matrix) std::free(score_matrix);
}

int32_t MatrixAlignmentScorer::score_exact_match(const Alignment&, size_t, size_t length) const {
    return match * length;
}

int32_t MatrixAlignmentScorer::score_exact_match(const string& sequence) const {
    return match * sequence.length();
}

int32_t MatrixAlignmentScorer::score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end) const {
    return match * (seq_end - seq_begin);
}

int32_t MatrixAlignmentScorer::score_exact_match(const string& sequence, const string&) const {
    return score_exact_match(sequence);
}

int32_t MatrixAlignmentScorer::score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end,
                                                 string::const_iterator) const {
    return score_exact_match(seq_begin, seq_end);
}

int32_t MatrixAlignmentScorer::score_mismatch(string::const_iterator seq_begin, string::const_iterator seq_end,
                                              string::const_iterator) const {
    return -mismatch * (seq_end - seq_begin);
}

int32_t MatrixAlignmentScorer::score_mismatch(size_t length) const {
    return -match * length;
}

int32_t MatrixAlignmentScorer::score_full_length_bonus(bool, string::const_iterator,
                                                       string::const_iterator,
                                                       string::const_iterator) const {
    return full_length_bonus;
}

int32_t MatrixAlignmentScorer::score_full_length_bonus(bool, const Alignment&) const {
    return full_length_bonus;
}

int32_t MatrixAlignmentScorer::score_partial_alignment(const Alignment& alignment, const HandleGraph& graph,
                                                       const path_t& path,
                                                       string::const_iterator seq_begin,
                                                       bool no_read_end_scoring) const {
    int32_t score = 0;
    string::const_iterator read_pos = seq_begin;
    bool in_deletion = false;
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        const auto& mapping = path.mapping(i);
        for (size_t j = 0; j < mapping.edit_size(); ++j) {
            const auto& edit = mapping.edit(j);
            if (edit.from_length() > 0) {
                if (edit.to_length() > 0) {
                    if (edit.sequence().empty()) {
                        // match
                        score += match * edit.from_length();
                    } else {
                        // mismatch
                        score -= mismatch * edit.from_length();
                    }
                    // apply full length bonus
                    if (read_pos == alignment.sequence().begin() && !no_read_end_scoring) {
                        score += score_full_length_bonus(true, alignment);
                    }
                    if (read_pos + edit.to_length() == alignment.sequence().end() && !no_read_end_scoring) {
                        score += score_full_length_bonus(false, alignment);
                    }
                    in_deletion = false;
                } else if (in_deletion) {
                    score -= edit.from_length() * gap_extension;
                } else {
                    // deletion
                    score -= gap_open + (edit.from_length() - 1) * gap_extension;
                    in_deletion = true;
                }
            } else if (edit.to_length() > 0) {
                // don't score soft clips if scoring read ends
                if (no_read_end_scoring ||
                    (read_pos != alignment.sequence().begin() &&
                     read_pos + edit.to_length() != alignment.sequence().end())) {
                    // insert
                    score -= gap_open + (edit.to_length() - 1) * gap_extension;
                }
                in_deletion = false;
            }
            read_pos += edit.to_length();
        }
    }
    return score;
}

double MatrixAlignmentScorer::get_log_base() const {
    // We have our log_base precomputed.
    return log_base;
}

// ----- QualAdjAlignmentScorer -----

QualAdjAlignmentScorer::QualAdjAlignmentScorer(const int8_t* score_matrix_4x4,
                                               int8_t gap_open,
                                               int8_t gap_extension,
                                               int8_t full_length_bonus,
                                               double gc_content)
    : MatrixAlignmentScorer(score_matrix_4x4, gap_open, gap_extension, full_length_bonus, gc_content) {

    constexpr uint32_t max_base_qual = 255;

    // Replace the 5x5 N-padded matrix with a quality-indexed 5x5xQ table.
    std::free(score_matrix);
    score_matrix = qual_adjusted_matrix(score_matrix_4x4, gc_content, log_base, max_base_qual);
    qual_adj_full_length_bonuses = qual_adjusted_bonuses(full_length_bonus, log_base, max_base_qual);
}

QualAdjAlignmentScorer::~QualAdjAlignmentScorer() {
    if (qual_adj_full_length_bonuses) std::free(qual_adj_full_length_bonuses);
}

int8_t* QualAdjAlignmentScorer::qual_adjusted_matrix(const int8_t* score_matrix_4x4,
                                                     double gc_content,
                                                     double log_base,
                                                     uint32_t max_qual) const {
    // TODO: duplicative with recover_log_base() in building nt_freqs
    double nt_freqs[4];
    nt_freqs[0] = 0.5 * (1 - gc_content);
    nt_freqs[1] = 0.5 * gc_content;
    nt_freqs[2] = 0.5 * gc_content;
    nt_freqs[3] = 0.5 * (1 - gc_content);

    // recover the emission probabilities of the align state of the HMM
    double align_prob[16];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            align_prob[i * 4 + j] = std::exp(log_base * score_matrix_4x4[i * 4 + j]) * nt_freqs[i] * nt_freqs[j];
        }
    }

    // compute the sum of the emission probabilities under a base error
    double align_complement_prob[16];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            align_complement_prob[i * 4 + j] = 0.0;
            for (int k = 0; k < 4; ++k) {
                if (k != j) {
                    align_complement_prob[i * 4 + j] += align_prob[i * 4 + k];
                }
            }
        }
    }

    // quality score of random guessing
    int lowest_meaningful_qual = std::ceil(-10.0 * std::log10(0.75));

    int8_t* qual_adj_mat = (int8_t*) std::malloc(25 * (max_qual + 1) * sizeof(int8_t));
    crash_unless(qual_adj_mat != nullptr);
    // compute the adjusted alignment scores for each quality level
    for (uint32_t q = 0; q <= max_qual; ++q) {
        double err = std::pow(10.0, -((double) q) / 10.0);
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                int8_t score;
                if (i == 4 || j == 4 || (int) q < lowest_meaningful_qual) {
                    score = 0;
                } else {
                    score = std::round(std::log(((1.0 - err) * align_prob[i * 4 + j] + (err / 3.0) * align_complement_prob[i * 4 + j])
                                                / (nt_freqs[i] * ((1.0 - err) * nt_freqs[j] + (err / 3.0) * (1.0 - nt_freqs[j])))) / log_base);
                }
                qual_adj_mat[q * 25 + i * 5 + j] = std::round(score);
            }
        }
    }
    return qual_adj_mat;
}

int8_t* QualAdjAlignmentScorer::qual_adjusted_bonuses(int8_t base_full_length_bonus,
                                                      double log_base,
                                                      uint32_t max_qual) const {
    double p_full_len = std::exp(log_base * base_full_length_bonus) / (1.0 + std::exp(log_base * base_full_length_bonus));

    int8_t* qual_adj_bonuses = (int8_t*) std::calloc(max_qual + 1, sizeof(int8_t));
    crash_unless(qual_adj_bonuses != nullptr);

    int lowest_meaningful_qual = std::ceil(-10.0 * std::log10(0.75));
    // hack because i want the minimum qual value from illumina (2) to have zero score, but phred
    // values are spaced out in a way to approximate this singularity well
    ++lowest_meaningful_qual;

    for (uint32_t q = lowest_meaningful_qual; q <= max_qual; ++q) {
        double err = std::pow(10.0, -((double) q) / 10.0);
        double score = std::log(((1.0 - err * 4.0 / 3.0) * p_full_len + (err * 4.0 / 3.0) * (1.0 - p_full_len)) / (1.0 - p_full_len)) / log_base;
        qual_adj_bonuses[q] = std::round(score);
    }
    return qual_adj_bonuses;
}

int32_t QualAdjAlignmentScorer::score_exact_match(const Alignment& aln, size_t read_offset, size_t length) const {
    auto& sequence = aln.sequence();
    auto& base_quality = aln.quality();
    int32_t score = 0;
    for (size_t i = 0; i < length; ++i) {
        score += score_matrix[25 * (uint8_t) base_quality[read_offset + i] + 6 * nt_table[(uint8_t) sequence[read_offset + i]]];
    }
    return score;
}

int32_t QualAdjAlignmentScorer::score_exact_match(const string& sequence, const string& base_quality) const {
    int32_t score = 0;
    for (size_t i = 0; i < sequence.length(); ++i) {
        score += score_matrix[25 * (uint8_t) base_quality[i] + 6 * nt_table[(uint8_t) sequence[i]]];
    }
    return score;
}

int32_t QualAdjAlignmentScorer::score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end,
                                                  string::const_iterator base_qual_begin) const {
    int32_t score = 0;
    auto seq_iter = seq_begin;
    auto qual_iter = base_qual_begin;
    for (; seq_iter != seq_end; ++seq_iter, ++qual_iter) {
        score += score_matrix[25 * (uint8_t) (*qual_iter) + 6 * nt_table[(uint8_t) (*seq_iter)]];
    }
    return score;
}

int32_t QualAdjAlignmentScorer::score_mismatch(string::const_iterator seq_begin, string::const_iterator seq_end,
                                               string::const_iterator base_qual_begin) const {
    int32_t score = 0;
    auto seq_iter = seq_begin;
    auto qual_iter = base_qual_begin;
    for (; seq_iter != seq_end; ++seq_iter, ++qual_iter) {
        score += score_matrix[25 * (uint8_t) (*qual_iter) + 1];
    }
    return score;
}

int32_t QualAdjAlignmentScorer::score_full_length_bonus(bool left_side, string::const_iterator seq_begin,
                                                        string::const_iterator seq_end,
                                                        string::const_iterator base_qual_begin) const {
    if (seq_begin != seq_end) {
        return qual_adj_full_length_bonuses[(uint8_t) (left_side ? *base_qual_begin
                                                                 : *(base_qual_begin + (seq_end - seq_begin) - 1))];
    }
    return 0;
}

int32_t QualAdjAlignmentScorer::score_full_length_bonus(bool left_side, const Alignment& alignment) const {
    return score_full_length_bonus(left_side, alignment.sequence().begin(), alignment.sequence().end(),
                                   alignment.quality().begin());
}

int32_t QualAdjAlignmentScorer::score_partial_alignment(const Alignment& alignment, const HandleGraph& graph,
                                                        const path_t& path,
                                                        string::const_iterator seq_begin,
                                                        bool no_read_end_scoring) const {
    int32_t score = 0;
    string::const_iterator read_pos = seq_begin;
    string::const_iterator qual_pos = alignment.quality().begin() + (seq_begin - alignment.sequence().begin());

    bool in_deletion = false;
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        const auto& mapping = path.mapping(i);

        // get the sequence of this node on the proper strand
        string node_seq = graph.get_sequence(graph.get_handle(mapping.position().node_id(),
                                                              mapping.position().is_reverse()));
        string::const_iterator ref_pos = node_seq.begin() + mapping.position().offset();

        for (size_t j = 0; j < mapping.edit_size(); ++j) {
            const auto& edit = mapping.edit(j);
            if (edit.from_length() > 0) {
                if (edit.to_length() > 0) {
                    auto siter = read_pos;
                    auto riter = ref_pos;
                    auto qiter = qual_pos;
                    for (; siter != read_pos + edit.to_length(); ++siter, ++qiter, ++riter) {
                        score += score_matrix[25 * (uint8_t) (*qiter) + 5 * nt_table[(uint8_t) (*riter)] + nt_table[(uint8_t) (*siter)]];
                    }

                    // apply full length bonus
                    if (read_pos == alignment.sequence().begin() && !no_read_end_scoring) {
                        score += score_full_length_bonus(true, alignment);
                    }
                    if (read_pos + edit.to_length() == alignment.sequence().end() && !no_read_end_scoring) {
                        score += score_full_length_bonus(false, alignment);
                    }
                    in_deletion = false;
                } else if (in_deletion) {
                    score -= edit.from_length() * gap_extension;
                } else {
                    // deletion
                    score -= gap_open + (edit.from_length() - 1) * gap_extension;
                    in_deletion = true;
                }
            } else if (edit.to_length() > 0) {
                // don't score soft clips if read end scoring
                if (no_read_end_scoring ||
                    (read_pos != alignment.sequence().begin() &&
                     read_pos + edit.to_length() != alignment.sequence().end())) {
                    score -= gap_open + (edit.to_length() - 1) * gap_extension;
                }
                in_deletion = false;
            }
            read_pos += edit.to_length();
            qual_pos += edit.to_length();
            ref_pos += edit.from_length();
        }
    }
    return score;
}

} // namespace vg
