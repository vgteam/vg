#include "mapping_quality_calculator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

#include "alignment.hpp"
#include "statistics.hpp"

namespace vg {

using std::log;
using std::log10;
using std::max;
using std::min;
using std::pair;
using std::pow;
using std::round;
using std::sqrt;
using std::vector;

/// Phred scale conversion factor.
static const double quality_scale_factor = 10.0 / std::log(10.0);

double MappingQualityCalculator::maximum_mapping_quality_exact(const vector<double>& scaled_scores,
                                                               size_t* max_idx_out,
                                                               const vector<double>* multiplicities) {
    // work in log transformed values to avoid risk of overflow
    double log_sum_exp = std::numeric_limits<double>::lowest();
    double to_score = std::numeric_limits<double>::lowest();

    // go in reverse order because this has fewer numerical problems when the scores are sorted (as usual)
    for (int64_t i = scaled_scores.size() - 1; i >= 0; --i) {
        // get the value of one copy of the score and check if it's the max
        double score = scaled_scores.at(i);
        if (max_idx_out && score >= to_score) {
            // Since we are going in reverse order, make sure to break ties in favor of the earlier item.
            *max_idx_out = i;
            to_score = score;
        }

        // add all copies of the score
        if (multiplicities && multiplicities->at(i) > 1.0) {
            score += log(multiplicities->at(i));
        }

        // accumulate the sum of all score
        log_sum_exp = add_log(log_sum_exp, score);
    }

    // if necessary, assume a null alignment of 0.0 for comparison since this is local
    if (scaled_scores.size() == 1) {
        if (multiplicities && multiplicities->at(0) <= 1.0) {
            log_sum_exp = add_log(log_sum_exp, 0.0);
        } else if (!multiplicities) {
            log_sum_exp = add_log(log_sum_exp, 0.0);
        }
    }

    if (!max_idx_out) {
        to_score = scaled_scores.empty() ? 0.0 : scaled_scores.front();
    }

    double direct_mapq = -quality_scale_factor * subtract_log(0.0, to_score - log_sum_exp);
    return std::isinf(direct_mapq) ? (double) std::numeric_limits<int32_t>::max() : direct_mapq;
}

double MappingQualityCalculator::maximum_mapping_quality_approx(const vector<double>& scaled_scores,
                                                                size_t* max_idx_out,
                                                                const vector<double>* multiplicities) {
    assert(!scaled_scores.empty());

    // TODO: this isn't very well-named now that it also supports computing non-maximum
    // mapping qualities

    // determine the maximum score and the count of the next highest score
    double max_score = scaled_scores.at(0);
    size_t max_idx = 0;

    // we start with the possibility of a null score of 0.0
    double next_score = 0.0;
    double next_count = 1.0;

    if (multiplicities) {
        if (multiplicities->at(0) > 1.0) {
            // there are extra copies of this one, so we'll init with those
            next_score = max_score;
            next_count = multiplicities->at(0) - 1.0;
        }
    }

    for (size_t i = 1; i < scaled_scores.size(); ++i) {
        double score = scaled_scores.at(i);
        if (score > max_score) {
            if (multiplicities && multiplicities->at(i) > 1.0) {
                // there are extra counts of the new highest score due to multiplicity
                next_score = score;
                next_count = multiplicities->at(i) - 1.0;
            } else if (next_score == max_score) {
                // the next highest was the same score as the old max, so we can
                // add its count back in
                next_count += 1.0;
            } else {
                // the old max score is now the second highest
                next_score = max_score;
                next_count = multiplicities ? multiplicities->at(max_idx) : 1.0;
            }
            max_score = score;
            max_idx = i;
        } else if (score > next_score) {
            // the new score is the second highest
            next_score = score;
            next_count = multiplicities ? multiplicities->at(i) : 1.0;
        } else if (score == next_score) {
            // the new score ties the second highest, so we combine their counts
            next_count += multiplicities ? multiplicities->at(i) : 1.0;
        }
    }

    // record the index of the highest score
    if (max_idx_out) {
        // we're either returning the mapping quality of whichever was the best, or we're
        // returning the mapping quality of the first, which also is the best
        *max_idx_out = max_idx;
    }
    if (max_idx_out || max_idx == 0) {
        return max(0.0, quality_scale_factor * (max_score - next_score - (next_count > 1.0 ? log(next_count) : 0.0)));
    } else {
        // we're returning the mapping quality of the first, which is not the best. the approximation
        // gets complicated here, so lets just fall back on the exact computation
        return maximum_mapping_quality_exact(scaled_scores, nullptr, multiplicities);
    }
}

double MappingQualityCalculator::first_mapping_quality_exact(const vector<double>& scaled_scores,
                                                             const vector<double>* multiplicities) {
    return maximum_mapping_quality_exact(scaled_scores, nullptr, multiplicities);
}

double MappingQualityCalculator::first_mapping_quality_approx(const vector<double>& scaled_scores,
                                                              const vector<double>* multiplicities) {
    return maximum_mapping_quality_approx(scaled_scores, nullptr, multiplicities);
}

double MappingQualityCalculator::group_mapping_quality_exact(const vector<double>& scaled_scores,
                                                             const vector<size_t>& group,
                                                             const vector<double>* multiplicities) const {
    // work in log transformed values to avoid risk of overflow
    double total_log_sum_exp = std::numeric_limits<double>::lowest();
    double non_group_log_sum_exp = std::numeric_limits<double>::lowest();

    // go in reverse order because this has fewer numerical problems when the scores are sorted (as usual)
    int64_t group_idx = group.size() - 1;
    for (int64_t i = scaled_scores.size() - 1; i >= 0; --i) {

        // the score of one alignment
        double score = scaled_scores.at(i);

        // the score all the multiples of this score combined
        double multiple_score = score;
        if (multiplicities && multiplicities->at(i) > 1.0) {
            multiple_score += log(multiplicities->at(i));
        }

        total_log_sum_exp = add_log(total_log_sum_exp, multiple_score);

        if (group_idx >= 0 && i == (int64_t) group[group_idx]) {
            // this is the next index in the group
            group_idx--;
            if (multiplicities && multiplicities->at(i) > 1.0) {
                // there's some remaining multiples of this score that don't get added into the group
                non_group_log_sum_exp = add_log(non_group_log_sum_exp,
                                                score + log(multiplicities->at(i) - 1.0));
            }
        } else {
            // this index is not part of the group
            non_group_log_sum_exp = add_log(non_group_log_sum_exp, multiple_score);
        }
    }

    if (scaled_scores.size() == 1) {
        if ((multiplicities && multiplicities->at(0) <= 1.0) || !multiplicities) {
            // assume a null alignment of 0.0 for comparison since this is local
            // TODO: repetitive, do I need to be this careful to not deref a null?
            non_group_log_sum_exp = add_log(non_group_log_sum_exp, 0.0);
            total_log_sum_exp = add_log(total_log_sum_exp, 0.0);
        }
    }

    double direct_mapq = quality_scale_factor * (total_log_sum_exp - non_group_log_sum_exp);
    return (std::isinf(direct_mapq) || direct_mapq > std::numeric_limits<int32_t>::max()) ?
           (double) std::numeric_limits<int32_t>::max() : direct_mapq;
}

vector<double> MappingQualityCalculator::all_mapping_qualities_exact(const vector<double>& scaled_scores,
                                                                     const vector<double>* multiplicities) const {
    vector<double> mapping_qualities(scaled_scores.size());

    // iterate backwards for improved numerical performance in sorted scores
    double log_denom = 0.0;
    for (int64_t i = scaled_scores.size() - 1; i >= 0; --i) {
        double score = scaled_scores[i];
        if (multiplicities && (*multiplicities)[i] != 1.0) {
            score += log((*multiplicities)[i]);
        }
        log_denom = add_log(log_denom, score);
    }
    // compute the mapping qualities
    for (size_t i = 0; i < scaled_scores.size(); ++i) {
        double log_prob_error = log10(1.0 - std::exp(scaled_scores[i] - log_denom));
        if (std::isnormal(log_prob_error) || log_prob_error == 0.0) {
            mapping_qualities[i] = -10.0 * log_prob_error;
        } else {
            mapping_qualities[i] = (double) std::numeric_limits<int32_t>::max();
        }
    }
    return mapping_qualities;
}

void MappingQualityCalculator::compute_mapping_quality(vector<Alignment>& alignments,
                                                       int max_mapping_quality,
                                                       bool fast_approximation,
                                                       double cluster_mq,
                                                       bool use_cluster_mq,
                                                       int overlap_count,
                                                       double mq_estimate,
                                                       double maybe_mq_threshold,
                                                       double identity_weight) const {
    assert(log_base > 0.0);
    if (alignments.empty()) return;

    vector<double> scaled_scores(alignments.size());
    for (size_t i = 0; i < alignments.size(); ++i) {
        scaled_scores[i] = log_base * alignments[i].score();
    }

    double mapping_quality;
    size_t max_idx;
    if (!fast_approximation) {
        mapping_quality = maximum_mapping_quality_exact(scaled_scores, &max_idx);
    } else {
        mapping_quality = maximum_mapping_quality_approx(scaled_scores, &max_idx);
    }

    if (use_cluster_mq) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(cluster_mq + mapping_quality)));
    }
    if (overlap_count) {
        mapping_quality -= quality_scale_factor * log(overlap_count);
    }

    auto& max_aln = alignments.at(max_idx);
    int l = max(alignment_to_length(max_aln), alignment_from_length(max_aln));
    double identity = 1. - (double)(l * rep_match - max_aln.score()) / (rep_match + rep_mismatch) / l;

    mapping_quality /= 2;
    mapping_quality *= pow(identity, identity_weight);

    if (mq_estimate < maybe_mq_threshold && mq_estimate < mapping_quality) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(mq_estimate + mapping_quality)));
    }
    if (mapping_quality > max_mapping_quality) {
        mapping_quality = max_mapping_quality;
    }
    if (alignments[max_idx].score() == 0) {
        mapping_quality = 0;
    }

    alignments[max_idx].set_mapping_quality(max(0, (int32_t) round(mapping_quality)));
    for (size_t i = 1; i < alignments.size(); ++i) {
        alignments[0].add_secondary_score(alignments[i].score());
    }
}

void MappingQualityCalculator::compute_paired_mapping_quality(pair<vector<Alignment>, vector<Alignment>>& alignment_pairs,
                                                              const vector<double>& frag_weights,
                                                              int max_mapping_quality1,
                                                              int max_mapping_quality2,
                                                              bool fast_approximation,
                                                              double cluster_mq,
                                                              bool use_cluster_mq,
                                                              int overlap_count1,
                                                              int overlap_count2,
                                                              double mq_estimate1,
                                                              double mq_estimate2,
                                                              double maybe_mq_threshold,
                                                              double identity_weight) const {
    assert(log_base > 0.0);

    size_t size = min(alignment_pairs.first.size(), alignment_pairs.second.size());
    if (size == 0) return;

    vector<double> scaled_scores(size);
    for (size_t i = 0; i < size; ++i) {
        scaled_scores[i] = log_base * (alignment_pairs.first[i].score() + alignment_pairs.second[i].score());
        // + frag_weights[i]);
        // ^^^ we could also incorporate the fragment weights, but this does not seem to help performance in the current form
    }

    size_t max_idx;
    double mapping_quality;
    if (!fast_approximation) {
        mapping_quality = maximum_mapping_quality_exact(scaled_scores, &max_idx);
    } else {
        mapping_quality = maximum_mapping_quality_approx(scaled_scores, &max_idx);
    }
    if (use_cluster_mq) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(cluster_mq + mapping_quality)));
    }

    double mapping_quality1 = mapping_quality;
    double mapping_quality2 = mapping_quality;
    if (overlap_count1) mapping_quality1 -= quality_scale_factor * log(overlap_count1);
    if (overlap_count2) mapping_quality2 -= quality_scale_factor * log(overlap_count2);

    auto& max_aln1 = alignment_pairs.first.at(max_idx);
    int len1 = max(alignment_to_length(max_aln1), alignment_from_length(max_aln1));
    double identity1 = 1. - (double)(len1 * rep_match - max_aln1.score()) / (rep_match + rep_mismatch) / len1;
    auto& max_aln2 = alignment_pairs.second.at(max_idx);
    int len2 = max(alignment_to_length(max_aln2), alignment_from_length(max_aln2));
    double identity2 = 1. - (double)(len2 * rep_match - max_aln2.score()) / (rep_match + rep_mismatch) / len2;

    mapping_quality1 /= 2;
    mapping_quality2 /= 2;
    mapping_quality1 *= pow(identity1, identity_weight);
    mapping_quality2 *= pow(identity2, identity_weight);

    double mq_estimate = min(mq_estimate1, mq_estimate2);
    if (mq_estimate < maybe_mq_threshold && mq_estimate < mapping_quality1) {
        mapping_quality1 = prob_to_phred(sqrt(phred_to_prob(mq_estimate + mapping_quality1)));
    }
    if (mq_estimate < maybe_mq_threshold && mq_estimate < mapping_quality2) {
        mapping_quality2 = prob_to_phred(sqrt(phred_to_prob(mq_estimate + mapping_quality2)));
    }

    if (mapping_quality1 > max_mapping_quality1) mapping_quality1 = max_mapping_quality1;
    if (mapping_quality2 > max_mapping_quality2) mapping_quality2 = max_mapping_quality2;

    if (alignment_pairs.first[max_idx].score() == 0) mapping_quality1 = 0;
    if (alignment_pairs.second[max_idx].score() == 0) mapping_quality2 = 0;

    int32_t final_mq = max(0, (int32_t) round(min(mapping_quality1, mapping_quality2)));
    alignment_pairs.first[max_idx].set_mapping_quality(final_mq);
    alignment_pairs.second[max_idx].set_mapping_quality(final_mq);

    for (size_t i = 1; i < alignment_pairs.first.size(); ++i) {
        alignment_pairs.first[0].add_secondary_score(alignment_pairs.first[i].score());
    }
    for (size_t i = 1; i < alignment_pairs.second.size(); ++i) {
        alignment_pairs.second[0].add_secondary_score(alignment_pairs.second[i].score());
    }
}

int32_t MappingQualityCalculator::compute_max_mapping_quality(const vector<double>& scores, bool fast_approximation,
                                                              const vector<double>* multiplicities) const {
    vector<double> scaled_scores(scores.size());
    for (size_t i = 0; i < scores.size(); ++i) {
        scaled_scores[i] = log_base * scores[i];
    }
    size_t idx;
    return (int32_t) (fast_approximation ? maximum_mapping_quality_approx(scaled_scores, &idx, multiplicities)
                                         : maximum_mapping_quality_exact(scaled_scores, &idx, multiplicities));
}

int32_t MappingQualityCalculator::compute_first_mapping_quality(const vector<double>& scores, bool fast_approximation,
                                                                const vector<double>* multiplicities) const {
    vector<double> scaled_scores(scores.size());
    for (size_t i = 0; i < scores.size(); ++i) {
        scaled_scores[i] = log_base * scores[i];
    }
    return (int32_t) (fast_approximation ? first_mapping_quality_approx(scaled_scores, multiplicities)
                                         : first_mapping_quality_exact(scaled_scores, multiplicities));
}

int32_t MappingQualityCalculator::compute_group_mapping_quality(const vector<double>& scores,
                                                                const vector<size_t>& group,
                                                                const vector<double>* multiplicities) const {
    // make a non-const local version in case we need to sort it
    vector<size_t> non_const_group;
    const vector<size_t>* grp_ptr = &group;

    // ensure that group is in sorted order as following function expects
    if (!std::is_sorted(group.begin(), group.end())) {
        non_const_group = group;
        std::sort(non_const_group.begin(), non_const_group.end());
        grp_ptr = &non_const_group;
    }

    vector<double> scaled_scores(scores.size(), 0.0);
    for (size_t i = 0; i < scores.size(); ++i) {
        scaled_scores[i] = log_base * scores[i];
    }
    return group_mapping_quality_exact(scaled_scores, *grp_ptr, multiplicities);
}

vector<int32_t> MappingQualityCalculator::compute_all_mapping_qualities(const vector<double>& scores,
                                                                        const vector<double>* multiplicities) const {
    vector<double> scaled_scores(scores.size(), 0.0);
    for (size_t i = 0; i < scores.size(); ++i) {
        scaled_scores[i] = log_base * scores[i];
    }
    vector<double> double_mapqs = all_mapping_qualities_exact(scaled_scores, multiplicities);
    vector<int32_t> to_return(double_mapqs.size(), 0);
    for (size_t i = 0; i < to_return.size(); ++i) {
        to_return[i] = double_mapqs[i];
    }
    return to_return;
}

double MappingQualityCalculator::mapping_quality_score_diff(double mapping_quality) const {
    return mapping_quality / (quality_scale_factor * log_base);
}

double MappingQualityCalculator::score_to_unnormalized_likelihood_ln(double score) const {
    // Log base needs to be set, or this can't work.
    assert(log_base != 0);
    // Likelihood is proportional to e^(lambda * score), so ln is just the exponent.
    return log_base * score;
}

double MappingQualityCalculator::estimate_next_best_score(int length, double min_diffs) const {
    return ((length - min_diffs) * rep_match - min_diffs * rep_mismatch);
}

double MappingQualityCalculator::max_possible_mapping_quality(int length) const {
    double max_score = log_base * length * rep_match;
    vector<double> v = { max_score };
    size_t max_idx;
    return maximum_mapping_quality_approx(v, &max_idx);
}

double MappingQualityCalculator::estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs) const {
    double max_score = log_base * ((length - min_diffs) * rep_match - min_diffs * rep_mismatch);
    double next_max_score = log_base * ((length - next_min_diffs) * rep_match - next_min_diffs * rep_mismatch);
    vector<double> v = { max_score, next_max_score };
    size_t max_idx;
    return maximum_mapping_quality_approx(v, &max_idx);
}

} // namespace vg
