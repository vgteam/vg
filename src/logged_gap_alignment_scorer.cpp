#include "logged_gap_alignment_scorer.hpp"

#include <algorithm>
#include <cmath>

namespace vg {

LoggedGapAlignmentScorer::LoggedGapAlignmentScorer(const Alignment& reference)
    : reference_aln(&reference) {
    count_alignment_operations(reference, reference_matches, reference_mismatches, reference_gap_lengths);
    d = recover_d(reference_matches, reference_mismatches, reference_gap_lengths);
    match = 1.0;
    mismatch = -1.0 / (2.0 * d);
}

double LoggedGapAlignmentScorer::recover_d(size_t matches, size_t mismatches,
                                           const std::vector<size_t>& gap_lengths) {
    double total = static_cast<double>(matches + mismatches + gap_lengths.size());
    if (total == 0.0) return 0.02;
    return std::max(0.02, static_cast<double>(mismatches + gap_lengths.size()) / total);
}

int32_t LoggedGapAlignmentScorer::score_from_counts(size_t matches, size_t mismatches,
                                                    const std::vector<size_t>& gap_lengths) const {
    double non_match_penalty = static_cast<double>(mismatches + gap_lengths.size()) / (2.0 * d);
    double indel_penalty = 0.0;
    for (size_t gap_length : gap_lengths) {
        indel_penalty += std::log2(1.0 + gap_length);
    }
    return std::round(static_cast<double>(matches) - non_match_penalty - indel_penalty);
}

int32_t LoggedGapAlignmentScorer::score_alignment(const Alignment& aln) const {
    if (&aln == reference_aln) {
        return score_from_counts(reference_matches, reference_mismatches, reference_gap_lengths);
    }
    size_t m, mm;
    std::vector<size_t> gaps;
    count_alignment_operations(aln, m, mm, gaps);
    return score_from_counts(m, mm, gaps);
}

void LoggedGapAlignmentScorer::fill_substitution_matrix(double out[16]) const {
    // Synthesize a uniform substitution matrix from the marginal scores.
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            out[i * 4 + j] = (i == j) ? match : mismatch;
        }
    }
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

} // namespace vg
