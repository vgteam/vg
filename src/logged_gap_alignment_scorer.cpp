#include "logged_gap_alignment_scorer.hpp"

#include <algorithm>
#include <cmath>

namespace vg {

LoggedGapAlignmentScorer::LoggedGapAlignmentScorer(const Alignment& standard, double gc_content)
    : LoggedGapAlignmentScorer(standard, count_alignment_operations(standard), gc_content) {
    
    // Nothing to do!
}

int32_t LoggedGapAlignmentScorer::score_alignment(const Alignment& aln) const {
    if (&aln == standard_address) {
        // Just score from the stored statistics
        return score_from_counts(matches, mismatches, gap_lengths);
    }
    // Otherwise, recompute and score from new statistics
    size_t m, mm;
    std::vector<size_t> gaps;
    count_alignment_operations(aln, m, mm, gaps);
    return score_from_counts(m, mm, gaps);
}

double LoggedGapAlignmentScorer::get_log_base() const {
    // TODO: This is the same approach as for the matrix-based scorer, but the
    // members have different heritage so we can't just re-use the code...
    return log_base;
}

LoggedGapAlignmentScorer::LoggedGapAlignmentScorer(
    const Alignment& standard,
    std::tuple<size_t, size_t, std::vector<size_t>>&& operation_counts,
    double gc_content
) : 
    // Very carefully initialize everything from only what preceeds it in the class definition
    matches(std::get<0>(operation_counts)),
    mismatches(std::get<1>(operation_counts)),
    gap_lengths(std::move(std::get<2>(operation_counts))),
    divergence(compute_divergence(matches, mismatches, gap_lengths)),
    mismatch(-1.0 / (2.0 * divergence)), 
    standard_address(&standard)
{
    // We have to go back and fill in log_base which is too hard to do in an initializer
    
    // Synthesize a uniform substitution matrix from the marginal scores.
    double matrix[16];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            matrix[i * 4 + j] = (i == j) ? match : mismatch;
        }
    }
    // Recover a log base from that
    log_base = AlignmentScorer::recover_log_base(matrix, gc_content);
}

int32_t LoggedGapAlignmentScorer::score_from_counts(size_t matches, size_t mismatches,
                                                    const std::vector<size_t>& gap_lengths) const {
    double non_match_penalty = static_cast<double>(mismatches + gap_lengths.size()) / (2.0 * divergence);
    double indel_penalty = 0.0;
    for (size_t gap_length : gap_lengths) {
        indel_penalty += std::log2(1.0 + gap_length);
    }
    return std::round(static_cast<double>(matches) - non_match_penalty - indel_penalty);
}

void LoggedGapAlignmentScorer::count_alignment_operations(const Alignment& aln,
                                                          size_t& matches,
                                                          size_t& mismatches,
                                                          std::vector<size_t>& gap_lengths) {
    matches = 0;
    mismatches = 0;
    gap_lengths.clear();

    enum class EditType { MATCH, MISMATCH, INS, DEL, COMPLEX, NONE };
    EditType prev_type = EditType::NONE;
    size_t current_gap_length = 0;

    auto finish_gap = [&]() {
        if (current_gap_length > 0) {
            gap_lengths.push_back(current_gap_length);
            current_gap_length = 0;
        }
    };

    for (size_t i = 0; i < (size_t) aln.path().mapping_size(); ++i) {
        auto& mapping = aln.path().mapping(i);
        for (size_t j = 0; j < (size_t) mapping.edit_size(); ++j) {
            auto& edit = mapping.edit(j);
            if (edit.from_length() == edit.to_length() && edit.from_length() > 0) {
                finish_gap();
                if (edit.sequence().empty()) {
                    matches += edit.from_length();
                    prev_type = EditType::MATCH;
                } else {
                    mismatches += edit.from_length();
                    prev_type = EditType::MISMATCH;
                }
            } else if (edit.from_length() == 0 && edit.to_length() > 0) {
                if (prev_type != EditType::INS) finish_gap();
                current_gap_length += edit.to_length();
                prev_type = EditType::INS;
            } else if (edit.from_length() > 0 && edit.to_length() == 0) {
                if (prev_type != EditType::DEL) finish_gap();
                current_gap_length += edit.from_length();
                prev_type = EditType::DEL;
            } else {
                finish_gap();
                mismatches += std::max(edit.from_length(), edit.to_length());
                prev_type = EditType::COMPLEX;
            }
        }
    }
    finish_gap();
}

std::tuple<size_t, size_t, std::vector<size_t>> LoggedGapAlignmentScorer::count_alignment_operations(const Alignment& aln) {
    // Make the right shape of tuple
    std::tuple<size_t, size_t, std::vector<size_t>> to_return;
    // Fill it in
    count_alignment_operations(aln, std::get<0>(to_return), std::get<1>(to_return), std::get<2>(to_return));
    // Return it
    return to_return;
}

double LoggedGapAlignmentScorer::compute_divergence(size_t matches, size_t mismatches,
                                                    const std::vector<size_t>& gap_lengths) {
    double total = static_cast<double>(matches + mismatches + gap_lengths.size());
    if (total == 0.0) return 0.02;
    return std::max(0.02, static_cast<double>(mismatches + gap_lengths.size()) / total);
}

} // namespace vg
