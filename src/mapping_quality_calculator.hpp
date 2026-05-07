#ifndef VG_MAPPING_QUALITY_CALCULATOR_HPP_INCLUDED
#define VG_MAPPING_QUALITY_CALCULATOR_HPP_INCLUDED

#include <cstdint>
#include <utility>
#include <vector>

#include <vg/vg.pb.h>

#include "alignment_scorer.hpp"

namespace vg {

/**
 * Widget for computing mapping qualities from collections of alignment scores.
 *
 * Constructable from any AlignmentScorer that also exposes arithmetic `match`
 * and `mismatch` members. No reference to the scorer is retained after
 * construction.
 *
 * Responsible for caching the computed log base value.
 */
class MappingQualityCalculator {
public:
    template<typename ScorerT>
    MappingQualityCalculator(const ScorerT& scorer, double gc_content)
        : gc_content(gc_content),
          rep_match(static_cast<double>(scorer.match)),
          rep_mismatch(static_cast<double>(scorer.mismatch)),
          log_base(scorer.recover_log_base(gc_content)) {}

    /// Stores -10 * log_10(P_err) in alignment mapping_quality field. P_err
    /// is the probability that the alignment is not the correct one.
    void compute_mapping_quality(std::vector<Alignment>& alignments,
                                 int max_mapping_quality,
                                 bool fast_approximation,
                                 double cluster_mq,
                                 bool use_cluster_mq,
                                 int overlap_count,
                                 double mq_estimate,
                                 double maybe_mq_threshold,
                                 double identity_weight) const;

    /// Same for paired reads. Mapping qualities are stored in both alignments.
    void compute_paired_mapping_quality(std::pair<std::vector<Alignment>, std::vector<Alignment>>& alignment_pairs,
                                        const std::vector<double>& frag_weights,
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
                                        double identity_weight) const;

    /// Computes mapping quality for the optimal score in a vector of scores.
    int32_t compute_max_mapping_quality(const std::vector<double>& scores, bool fast_approximation,
                                        const std::vector<double>* multiplicities = nullptr) const;

    /// Computes mapping quality for the first score in a vector of scores.
    int32_t compute_first_mapping_quality(const std::vector<double>& scores, bool fast_approximation,
                                          const std::vector<double>* multiplicities = nullptr) const;

    /// Computes mapping quality for a group of scores.
    int32_t compute_group_mapping_quality(const std::vector<double>& scores, const std::vector<size_t>& group,
                                          const std::vector<double>* multiplicities = nullptr) const;

    /// Computes mapping quality for all of a vector of scores.
    std::vector<int32_t> compute_all_mapping_qualities(const std::vector<double>& scores,
                                                       const std::vector<double>* multiplicities = nullptr) const;

    /// Difference between optimal and second-best alignment scores that would
    /// produce the given mapping quality under the fast approximation.
    double mapping_quality_score_diff(double mapping_quality) const;

    /// Convert a score to an unnormalized log likelihood for the sequence.
    double score_to_unnormalized_likelihood_ln(double score) const;

    double max_possible_mapping_quality(int length) const;
    double estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs) const;
    double estimate_next_best_score(int length, double min_diffs) const;

    double get_log_base() const { return log_base; }
    double get_gc_content() const { return gc_content; }

    // Static helpers (no internal MQ state needed).

    /// Given a nonempty vector of nonnegative scaled alignment scores,
    /// compute the mapping quality of the maximal score. Sets *max_idx_out
    /// (if non-null) to the index of that score.
    static double maximum_mapping_quality_exact(const std::vector<double>& scaled_scores,
                                                size_t* max_idx_out,
                                                const std::vector<double>* multiplicities = nullptr);
    static double maximum_mapping_quality_approx(const std::vector<double>& scaled_scores,
                                                 size_t* max_idx_out,
                                                 const std::vector<double>* multiplicities = nullptr);
    static double first_mapping_quality_exact(const std::vector<double>& scaled_scores,
                                              const std::vector<double>* multiplicities = nullptr);
    static double first_mapping_quality_approx(const std::vector<double>& scaled_scores,
                                               const std::vector<double>* multiplicities = nullptr);

private:
    double group_mapping_quality_exact(const std::vector<double>& scaled_scores, const std::vector<size_t>& group,
                                       const std::vector<double>* multiplicities = nullptr) const;
    std::vector<double> all_mapping_qualities_exact(const std::vector<double>& scaled_scores,
                                                    const std::vector<double>* multiplicities = nullptr) const;

    double gc_content;
    double rep_match;
    double rep_mismatch;
    double log_base;
};

} // namespace vg

#endif
