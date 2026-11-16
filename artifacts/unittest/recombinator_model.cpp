/** \file
 *
 * Unit tests for the count model and the pair hidden Markov model in
 * recombinator_model.cpp.
 */

#include "../recombinator_model.hpp"

#include "../../src/unittest/catch.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

// A count model with the given parameters and no mixture weights.
pair_model::CountModel make_model(double mu, double eps, double r1, double r0) {
    pair_model::CountModel model;
    model.mu = mu; model.eps = eps; model.r1 = r1; model.r0 = r0;
    return model;
}

// The expected histogram of reported counts for `total` kmers with copy number
// `n` in proportion `weights[n]`.
std::vector<double> expected_histogram(const pair_model::CountModel& model, const std::vector<double>& weights, double total) {
    std::vector<double> histogram(pair_model::CountModel::MAX_COUNT + 1, 0.0);
    for (size_t y = 0; y <= pair_model::CountModel::MAX_COUNT; y++) {
        if (y > 0 && y < pair_model::CountModel::MIN_COUNT) {
            continue; // Reported as 0.
        }
        for (size_t n = 0; n < weights.size(); n++) {
            histogram[y] += total * weights[n] * std::exp(model.log_probability(y, n));
        }
    }
    return histogram;
}

// Kmer presence rows for sequences given as strings of `0` and `1`.
std::vector<std::vector<std::uint64_t>> presence_rows(const std::vector<std::string>& rows) {
    std::vector<std::vector<std::uint64_t>> result;
    for (const std::string& row : rows) {
        result.emplace_back((row.size() + 63) / 64, 0);
        for (size_t j = 0; j < row.size(); j++) {
            if (row[j] == '1') {
                result.back()[j / 64] |= std::uint64_t(1) << (j % 64);
            }
        }
    }
    return result;
}

// Candidates where group `g` has a single member: sequence `g` of haplotype `haplotypes[g]`.
pair_model::Candidates single_members(const std::vector<std::uint64_t>& haplotypes, const std::vector<double>& scores, double length) {
    pair_model::Candidates result;
    for (size_t g = 0; g < haplotypes.size(); g++) {
        result.members.push_back({ g });
        result.haplotypes.push_back({ haplotypes[g] });
    }
    result.scores = scores;
    result.length = length;
    return result;
}

// Returns the group that contains the sequence.
size_t group_of(const pair_model::Candidates& candidates, size_t sequence) {
    for (size_t g = 0; g < candidates.size(); g++) {
        for (size_t member : candidates.members[g]) {
            if (member == sequence) {
                return g;
            }
        }
    }
    return candidates.size();
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Kmer count model", "[haplotypes][recombinator]") {
    pair_model::CountModel truth = make_model(15.0, 0.5, 20.0, 0.5);

    SECTION("counts below the minimum are reported as zero") {
        for (size_t copies : { size_t(0), size_t(1), size_t(2) }) {
            REQUIRE(truth.log_probability(1, copies) == truth.log_probability(0, copies));
        }
    }

    SECTION("the probabilities of reported counts sum to one") {
        // With 20 copies, the mean is above the largest count.
        for (size_t copies : { size_t(0), size_t(1), size_t(2), size_t(5), size_t(20) }) {
            double total = std::exp(truth.log_probability(0, copies));
            for (size_t y = pair_model::CountModel::MIN_COUNT; y <= pair_model::CountModel::MAX_COUNT; y++) {
                total += std::exp(truth.log_probability(y, copies));
            }
            REQUIRE(total == Approx(1.0).epsilon(1e-6));
        }
    }

    SECTION("the fit recovers the parameters") {
        // With weight on copy numbers above 4, twice the copies at half the mean
        // would need more copies than the fit has.
        std::vector<double> weights { 0.45, 0.15, 0.25, 0.05, 0.04, 0.03, 0.03 };
        std::vector<double> histogram = expected_histogram(truth, weights, 1e6);
        pair_model::CountModel fitted = pair_model::CountModel::fit(histogram, 30.0);
        REQUIRE(fitted.mu == Approx(truth.mu).epsilon(0.01));
        REQUIRE(fitted.r1 == Approx(truth.r1).epsilon(0.05));
        REQUIRE(fitted.pi.size() == pair_model::CountModel::FIT_COPIES + 1);
        for (size_t n = 0; n < weights.size(); n++) {
            REQUIRE(fitted.pi[n] == Approx(weights[n]).margin(0.01));
        }
    }
}

TEST_CASE("Pair scores from the count model", "[haplotypes][recombinator]") {
    pair_model::CountModel model = make_model(15.0, 0.5, 20.0, 0.5);
    std::vector<double> flat(pair_model::CountModel::MAX_COUNT + 1, 1.0);
    std::vector<double> table = pair_model::count_table(model, flat);

    // The sample's counts for the pair of sequences with the given rows and copies.
    auto counts_for = [&](const std::vector<std::string>& rows, const std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>& copies,
                          size_t first, size_t second) -> std::vector<size_t> {
        std::vector<size_t> counts(rows.front().size(), 0);
        for (size_t sequence : { first, second }) {
            for (size_t j = 0; j < counts.size(); j++) {
                size_t c = (rows[sequence][j] == '1' ? 1 : 0);
                for (auto entry : copies[sequence]) {
                    if (entry.first == j) {
                        c = entry.second;
                    }
                }
                counts[j] += 15 * c;
            }
        }
        return counts;
    };

    SECTION("the best pair is the sample's pair") {
        std::vector<std::string> rows {
            "11000011",
            "10100101",
            "01010110",
            "10100101", // Same as sequence 1.
            "00001111"
        };
        std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> copies(rows.size());
        std::vector<bool> usable(rows.size(), true);
        std::vector<size_t> counts = counts_for(rows, copies, 1, 2);
        pair_model::Candidates candidates = pair_model::score_pairs(presence_rows(rows), copies, usable, counts, table);
        REQUIRE(candidates.size() == 4);
        REQUIRE(candidates.members[group_of(candidates, 1)] == std::vector<size_t>({ 1, 3 }));
        size_t best = std::max_element(candidates.scores.begin(), candidates.scores.end()) - candidates.scores.begin();
        size_t a = best / candidates.size(), b = best % candidates.size();
        std::vector<size_t> pair { std::min(a, b), std::max(a, b) };
        std::vector<size_t> expected { group_of(candidates, 1), group_of(candidates, 2) };
        std::sort(expected.begin(), expected.end());
        REQUIRE(pair == expected);
        REQUIRE(candidates.scores[a * candidates.size() + b] == Approx(candidates.scores[b * candidates.size() + a]));
    }

    SECTION("copies separate sequences with the same kmers") {
        // Sequences 0 and 1 have the same kmers, but sequence 1 has two copies of kmer 0.
        std::vector<std::string> rows {
            "1100",
            "1100",
            "1011"
        };
        std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> copies { {}, { { 0, 2 } }, {} };
        std::vector<bool> usable(rows.size(), true);
        std::vector<size_t> counts = counts_for(rows, copies, 1, 2);
        pair_model::Candidates candidates = pair_model::score_pairs(presence_rows(rows), copies, usable, counts, table);
        REQUIRE(candidates.size() == 3);
        size_t g0 = group_of(candidates, 0), g1 = group_of(candidates, 1), g2 = group_of(candidates, 2);
        REQUIRE(g0 != g1);
        size_t K = candidates.size();
        REQUIRE(candidates.scores[g1 * K + g2] > candidates.scores[g0 * K + g2]);
    }

    SECTION("each pair score is the sum of the table over the kmers") {
        // Kmers 0, 2, and 3 are carried by most sequences, and sequences 1 and 3
        // have multiple copies of a kmer. Kmers 2 and 4 are ignored.
        std::vector<std::string> rows {
            "1100111",
            "1010110",
            "0110101",
            "1111011",
            "0000111"
        };
        std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> copies { {}, { { 0, 3 } }, {}, { { 4, 2 } }, {} };
        std::vector<bool> usable(rows.size(), true);
        std::vector<size_t> counts = counts_for(rows, copies, 1, 3);
        counts[2] = counts[4] = pair_model::IGNORED_COUNT;
        pair_model::Candidates candidates = pair_model::score_pairs(presence_rows(rows), copies, usable, counts, table);
        REQUIRE(candidates.size() == rows.size());
        auto copies_of = [&](size_t sequence, size_t kmer) -> size_t {
            size_t result = (rows[sequence][kmer] == '1' ? 1 : 0);
            for (auto entry : copies[sequence]) {
                if (entry.first == kmer) {
                    result = entry.second;
                }
            }
            return result;
        };
        size_t K = candidates.size(), W = pair_model::PAIR_COPIES + 1;
        for (size_t i = 0; i < K; i++) {
            for (size_t k = 0; k < K; k++) {
                size_t a = candidates.members[i].front(), b = candidates.members[k].front();
                double expected = 0.0;
                for (size_t j = 0; j < counts.size(); j++) {
                    if (counts[j] != pair_model::IGNORED_COUNT) {
                        expected += table[counts[j] * W + std::min(copies_of(a, j) + copies_of(b, j), pair_model::PAIR_COPIES)];
                    }
                }
                REQUIRE(candidates.scores[i * K + k] == Approx(expected));
            }
        }
    }

    SECTION("candidates are the groups in the best pairs") {
        // Each sequence has its own kmer, and the sample has sequences 37 and 38.
        size_t n = 40;
        std::vector<std::string> rows(n, std::string(n, '0'));
        for (size_t i = 0; i < n; i++) {
            rows[i][i] = '1';
        }
        std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> copies(n);
        std::vector<bool> usable(n, true);
        pair_model::Candidates candidates = pair_model::score_pairs(presence_rows(rows), copies, usable, counts_for(rows, copies, 37, 38), table);
        REQUIRE(candidates.size() == pair_model::CANDIDATES);
        REQUIRE(group_of(candidates, 37) < candidates.size());
        REQUIRE(group_of(candidates, 38) < candidates.size());
        size_t best = std::max_element(candidates.scores.begin(), candidates.scores.end()) - candidates.scores.begin();
        std::set<size_t> found { best / candidates.size(), best % candidates.size() };
        std::set<size_t> expected { group_of(candidates, 37), group_of(candidates, 38) };
        REQUIRE(found == expected);
    }

    SECTION("ignored kmers and unusable sequences") {
        std::vector<std::string> rows {
            "1100",
            "0011",
            "1111"
        };
        std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> copies(rows.size());
        std::vector<size_t> counts = counts_for(rows, copies, 0, 0);
        std::vector<bool> usable { true, false, true };
        pair_model::Candidates candidates = pair_model::score_pairs(presence_rows(rows), copies, usable, counts, table);
        REQUIRE(group_of(candidates, 1) == candidates.size());

        // With every kmer ignored, all pairs have the same score.
        std::vector<size_t> ignored(counts.size(), pair_model::IGNORED_COUNT);
        usable.assign(rows.size(), true);
        candidates = pair_model::score_pairs(presence_rows(rows), copies, usable, ignored, table);
        for (double score : candidates.scores) {
            REQUIRE(score == 0.0);
        }
    }
}

TEST_CASE("Viterbi path over pairs of haplotypes", "[haplotypes][recombinator]") {
    // Haplotypes 10, 11, 12 in both subchains. The first subchain prefers (10, 11), and
    // the second prefers (12, 11) slightly over (10, 11).
    std::vector<double> first(9, -100.0), second(9, -100.0);
    first[0 * 3 + 1] = 0.0;
    second[0 * 3 + 1] = 0.0;
    second[2 * 3 + 1] = 0.1;
    std::vector<pair_model::Candidates> subchains {
        single_members({ 10, 11, 12 }, first, 1000.0),
        single_members({ 10, 11, 12 }, second, 1000.0)
    };

    SECTION("a member keeps its haplotype when switching is rare") {
        auto path = pair_model::viterbi(subchains, 1.0, 3e-6);
        REQUIRE(path.size() == 2);
        REQUIRE(path[0] == std::make_pair(size_t(0), size_t(1)));
        REQUIRE(path[1] == std::make_pair(size_t(0), size_t(1)));
    }

    SECTION("a member switches when switching is free") {
        auto path = pair_model::viterbi(subchains, 1.0, 1.0);
        REQUIRE(path[0] == std::make_pair(size_t(0), size_t(1)));
        REQUIRE(path[1] == std::make_pair(size_t(2), size_t(1)));
    }

    SECTION("a member follows its haplotype to another group") {
        // In the second subchain, haplotype 10 is in group 2, and the pair (12, 11)
        // scores slightly better than (10, 11).
        std::vector<double> scores(9, -100.0);
        scores[0 * 3 + 1] = 0.1;
        scores[2 * 3 + 1] = 0.0;
        subchains[1] = single_members({ 12, 11, 10 }, scores, 1000.0);
        auto path = pair_model::viterbi(subchains, 1.0, 3e-6);
        REQUIRE(path[1] == std::make_pair(size_t(2), size_t(1)));
    }

    SECTION("subchains with different numbers of candidates") {
        // The middle subchain has four groups in another order, and it prefers (12, 13)
        // slightly over (10, 11). With rare switching, the members keep 10 and 11.
        std::vector<double> middle(16, -100.0);
        middle[3 * 4 + 0] = 0.0;
        middle[1 * 4 + 2] = 0.1;
        subchains[1] = single_members({ 11, 12, 13, 10 }, middle, 1000.0);
        subchains.push_back(single_members({ 10, 11, 12 }, first, 1000.0));
        auto path = pair_model::viterbi(subchains, 1.0, 3e-6);
        REQUIRE(path.size() == 3);
        REQUIRE(path[0] == std::make_pair(size_t(0), size_t(1)));
        REQUIRE(path[1] == std::make_pair(size_t(3), size_t(0)));
        REQUIRE(path[2] == std::make_pair(size_t(0), size_t(1)));
        REQUIRE(pair_model::viterbi(subchains, 1.0, 1.0)[1] == std::make_pair(size_t(1), size_t(2)));
    }

    SECTION("tau divides the pair scores") {
        // Switching the first member to haplotype 12 costs about 6.9.
        std::vector<double> scores(9, -100.0);
        scores[0 * 3 + 1] = 0.0;
        scores[2 * 3 + 1] = 10.0;
        subchains[1] = single_members({ 10, 11, 12 }, scores, 1000.0);
        REQUIRE(pair_model::viterbi(subchains, 1.0, 3e-6)[1] == std::make_pair(size_t(2), size_t(1)));
        REQUIRE(pair_model::viterbi(subchains, 3.0, 3e-6)[1] == std::make_pair(size_t(0), size_t(1)));
    }
}

TEST_CASE("Sequences for the model pair", "[haplotypes][recombinator]") {
    pair_model::Candidates a, b;
    a.members = { { 4, 5 } };
    a.haplotypes = { { 100, 101 } };
    b.members = { { 7, 8, 9 } };
    b.haplotypes = { { 102, 101, 100 } };
    std::vector<pair_model::Candidates> subchains { a, b };

    SECTION("both members in one group take different sequences, then keep their haplotypes") {
        auto chosen = pair_model::choose_sequences(subchains, { { 0, 0 }, { 0, 0 } });
        REQUIRE(chosen[0] == std::make_pair(size_t(4), size_t(5)));
        REQUIRE(chosen[1] == std::make_pair(size_t(9), size_t(8)));
    }

    SECTION("members that keep their haplotypes may take the same sequence") {
        // Both members take haplotype 100 in the first subchain and keep it in the
        // second, although its group also has haplotype 101.
        subchains[0].members = { { 4 } };
        subchains[0].haplotypes = { { 100 } };
        auto chosen = pair_model::choose_sequences(subchains, { { 0, 0 }, { 0, 0 } });
        REQUIRE(chosen[0] == std::make_pair(size_t(4), size_t(4)));
        REQUIRE(chosen[1] == std::make_pair(size_t(9), size_t(9)));
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg
