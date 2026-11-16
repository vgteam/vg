#include "recombinator_model.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <unordered_map>

namespace vg {

namespace pair_model {

//------------------------------------------------------------------------------

// Numerical constants.

constexpr size_t CountModel::MIN_COUNT;
constexpr size_t CountModel::MAX_COUNT;
constexpr size_t CountModel::FIT_COPIES;

//------------------------------------------------------------------------------

namespace {

constexpr size_t TABLE_WIDTH = PAIR_COPIES + 1;

// Return log(exp(a) + exp(b)). Both may be -inf, which `add_log()` does not allow.
double log_add(double a, double b) {
    if (a < b) {
        std::swap(a, b);
    }
    if (b == -std::numeric_limits<double>::infinity()) {
        return a;
    }
    return a + std::log1p(std::exp(b - a));
}

// log P(Y = y) for Y ~ NB(size r, mean m).
double nb_log_pmf(double y, double r, double m) {
    double log_rm = std::log(r + m);
    return std::lgamma(y + r) - std::lgamma(r) - std::lgamma(y + 1.0) + r * (std::log(r) - log_rm) + y * (std::log(m) - log_rm);
}

// log P(Y >= from) for Y ~ NB(size r, mean m).
double nb_log_tail(size_t from, double r, double m) {
    double mode = (r > 1.0 ? (r - 1.0) * m / r : 0.0);
    if (mode < double(from)) {
        // The terms decrease from `from` on, so we sum them relative to the first.
        double q = m / (r + m), relative = 1.0, sum = 1.0;
        for (size_t k = from; k < from + 1000000; k++) {
            relative *= (double(k) + r) / (double(k) + 1.0) * q;
            sum += relative;
            if (relative < 1e-17 * sum) {
                break;
            }
        }
        return nb_log_pmf(double(from), r, m) + std::log(sum);
    }
    double head = 0.0;
    for (size_t k = 0; k < from; k++) {
        head += std::exp(nb_log_pmf(double(k), r, m));
    }
    return std::log1p(-std::min(head, 1.0));
}

// Minimize `f` within the box [lo, hi] by BFGS over a logistic transformation of
// the box, with central-difference gradients, starting from `x0` in the interior of
// the box. Return the minimizer.
std::vector<double> minimize_in_box(
    const std::function<double(const std::vector<double>&)>& f,
    const std::vector<double>& x0, const std::vector<double>& lo, const std::vector<double>& hi
) {
    size_t d = x0.size();
    auto to_x = [&](const std::vector<double>& z) -> std::vector<double> {
        std::vector<double> x(d);
        for (size_t i = 0; i < d; i++) {
            x[i] = lo[i] + (hi[i] - lo[i]) / (1.0 + std::exp(-z[i]));
        }
        return x;
    };
    auto value = [&](const std::vector<double>& z) -> double {
        return f(to_x(z));
    };
    auto gradient = [&](const std::vector<double>& z) -> std::vector<double> {
        std::vector<double> g(d), t = z;
        for (size_t i = 0; i < d; i++) {
            constexpr double h = 1e-5;
            t[i] = z[i] + h;
            double above = value(t);
            t[i] = z[i] - h;
            double below = value(t);
            t[i] = z[i];
            g[i] = (above - below) / (2.0 * h);
        }
        return g;
    };
    auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) -> double {
        double result = 0.0;
        for (size_t i = 0; i < d; i++) {
            result += a[i] * b[i];
        }
        return result;
    };

    std::vector<double> z(d);
    for (size_t i = 0; i < d; i++) {
        double u = (x0[i] - lo[i]) / (hi[i] - lo[i]);
        z[i] = std::log(u / (1.0 - u));
    }
    double fz = value(z);
    std::vector<double> g = gradient(z);
    std::vector<double> H(d * d, 0.0); // Approximate inverse Hessian.
    for (size_t i = 0; i < d; i++) {
        H[i * d + i] = 1.0;
    }
    bool scaled = false;
    size_t flat = 0;
    for (size_t iteration = 0; iteration < 5000; iteration++) {
        std::vector<double> step(d, 0.0);
        for (size_t i = 0; i < d; i++) {
            for (size_t j = 0; j < d; j++) {
                step[i] -= H[i * d + j] * g[j];
            }
        }
        double slope = dot(g, step);
        if (!(slope < 0.0)) {
            // Not a descent direction, so we restart from the gradient.
            std::fill(H.begin(), H.end(), 0.0);
            for (size_t i = 0; i < d; i++) {
                H[i * d + i] = 1.0;
                step[i] = -g[i];
            }
            scaled = false;
            slope = -dot(g, g);
            if (!(slope < 0.0)) {
                break;
            }
        }

        // Backtracking line search.
        std::vector<double> next(d);
        double f_next = fz, t = 1.0;
        bool found = false;
        for (size_t k = 0; k < 80; k++, t *= 0.5) {
            for (size_t i = 0; i < d; i++) {
                next[i] = z[i] + t * step[i];
            }
            f_next = value(next);
            if (f_next <= fz + 1e-4 * t * slope) {
                found = true;
                break;
            }
        }
        if (!found) {
            break;
        }

        // BFGS update of the inverse Hessian.
        std::vector<double> g_next = gradient(next);
        std::vector<double> s(d), y(d);
        for (size_t i = 0; i < d; i++) {
            s[i] = next[i] - z[i];
            y[i] = g_next[i] - g[i];
        }
        double sy = dot(s, y);
        if (sy > 1e-12 * std::sqrt(dot(s, s) * dot(y, y))) {
            if (!scaled) {
                double c = sy / dot(y, y);
                for (size_t i = 0; i < d; i++) {
                    H[i * d + i] = c;
                }
                scaled = true;
            }
            double rho = 1.0 / sy;
            std::vector<double> Hy(d, 0.0);
            for (size_t i = 0; i < d; i++) {
                for (size_t j = 0; j < d; j++) {
                    Hy[i] += H[i * d + j] * y[j];
                }
            }
            double yHy = dot(y, Hy);
            for (size_t i = 0; i < d; i++) {
                for (size_t j = 0; j < d; j++) {
                    H[i * d + j] += -rho * (s[i] * Hy[j] + Hy[i] * s[j]) + (rho * rho * yHy + rho) * s[i] * s[j];
                }
            }
        }

        // Stop after five steps that do not decrease the value.
        double decrease = fz - f_next;
        z = next; g = g_next; fz = f_next;
        flat = (decrease <= 1e-13 * std::abs(fz) ? flat + 1 : 0);
        if (flat >= 5) {
            break;
        }
    }
    return to_x(z);
}

std::uint64_t mix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Return the transition probabilities (Ka x Kb) of a member of the pair from the
// groups of `a` to the groups of `b` with move probability r.
std::vector<double> transition(const Candidates& a, const Candidates& b, double r) {
    size_t Ka = a.size(), Kb = b.size();
    std::unordered_map<std::uint64_t, size_t> where;
    for (size_t gb = 0; gb < Kb; gb++) {
        for (std::uint64_t x : b.haplotypes[gb]) {
            where[x] = gb;
        }
    }
    std::vector<double> A(Ka * Kb, 0.0);
    for (size_t ga = 0; ga < Ka; ga++) {
        double share = 1.0 / double(a.haplotypes[ga].size());
        for (std::uint64_t x : a.haplotypes[ga]) {
            auto iter = where.find(x);
            if (iter != where.end()) {
                A[ga * Kb + iter->second] += share;
            }
        }
        double kept = 0.0;
        for (size_t gb = 0; gb < Kb; gb++) {
            kept += A[ga * Kb + gb];
        }
        double leak = 1.0 - kept;
        for (size_t gb = 0; gb < Kb; gb++) {
            A[ga * Kb + gb] = (1.0 - r) * (A[ga * Kb + gb] + leak / double(Kb)) + r / double(Kb);
        }
    }
    return A;
}

} // anonymous namespace

//------------------------------------------------------------------------------

std::pair<double, double> CountModel::nb_parameters(size_t copies) const {
    if (copies == 0) {
        return std::make_pair(this->r0, this->eps);
    }
    return std::make_pair(this->r1 * double(copies), this->mu * double(copies));
}

double CountModel::log_probability(size_t count, size_t copies) const {
    std::pair<double, double> nb = this->nb_parameters(copies);
    if (count < MIN_COUNT) {
        double result = -std::numeric_limits<double>::infinity();
        for (size_t k = 0; k < MIN_COUNT; k++) {
            result = log_add(result, nb_log_pmf(double(k), nb.first, nb.second));
        }
        return result;
    }
    if (count >= MAX_COUNT) {
        return nb_log_tail(MAX_COUNT, nb.first, nb.second);
    }
    return nb_log_pmf(double(count), nb.first, nb.second);
}

CountModel CountModel::fit(const std::vector<double>& histogram, double coverage) {
    // Counts above `censor` are only known to be above it. Bin `ys.size()` holds them.
    size_t censor = std::min(size_t(double(FIT_COPIES) * coverage / 2.0), MAX_COUNT - 1);
    std::vector<size_t> ys;
    std::vector<double> weights;
    double above = 0.0;
    for (size_t y = 0; y < histogram.size(); y++) {
        if (y > censor) {
            above += histogram[y];
        } else if (histogram[y] > 0.0) {
            ys.push_back(y);
            weights.push_back(histogram[y]);
        }
    }

    // Negative log-likelihood of x = { log mu, log eps, log r1, log r0, unnormalized log pi }.
    auto nll = [&](const std::vector<double>& x) -> double {
        CountModel model;
        model.mu = std::exp(x[0]); model.eps = std::exp(x[1]); model.r1 = std::exp(x[2]); model.r0 = std::exp(x[3]);
        double top = *std::max_element(x.begin() + 4, x.end()), norm = 0.0;
        for (size_t n = 0; n <= FIT_COPIES; n++) {
            norm += std::exp(x[4 + n] - top);
        }
        double result = 0.0;
        for (size_t b = 0; b <= ys.size(); b++) {
            double weight = (b < ys.size() ? weights[b] : above);
            if (weight == 0.0) {
                continue;
            }
            double total = -std::numeric_limits<double>::infinity();
            for (size_t n = 0; n <= FIT_COPIES; n++) {
                std::pair<double, double> nb = model.nb_parameters(n);
                double lp = (b < ys.size() ? model.log_probability(ys[b], n) : nb_log_tail(censor + 1, nb.first, nb.second));
                total = log_add(total, lp + x[4 + n] - top - std::log(norm));
            }
            result -= weight * total;
        }
        return result;
    };
    std::vector<double> lo = { std::log(coverage / 8.0), std::log(1e-3), std::log(0.5), std::log(1e-3) };
    std::vector<double> hi = { std::log(coverage), std::log(10.0), std::log(1e4), std::log(1e3) };
    std::vector<double> initial_pi = { 0.5, 0.2, 0.2, 0.05, 0.02 };
    for (size_t n = 0; n <= FIT_COPIES; n++) {
        lo.push_back(-30.0);
        hi.push_back(30.0);
    }

    // The likelihood has local optima, so we start from a grid of eps and r1 and keep the best.
    std::vector<double> x;
    double best = std::numeric_limits<double>::infinity();
    for (double eps : { 0.05, 1.0, 3.0 }) {
        for (double r1 : { 5.0, 20.0, 80.0 }) {
            std::vector<double> x0 = { std::log(coverage / 2.0), std::log(eps), std::log(r1), std::log(0.5) };
            for (size_t n = 0; n <= FIT_COPIES; n++) {
                x0.push_back(std::log(n < initial_pi.size() ? initial_pi[n] : 0.01));
            }
            std::vector<double> candidate = minimize_in_box(nll, x0, lo, hi);
            double value = nll(candidate);
            if (value < best) {
                best = value;
                x = candidate;
            }
        }
    }

    CountModel model;
    model.mu = std::exp(x[0]); model.eps = std::exp(x[1]); model.r1 = std::exp(x[2]); model.r0 = std::exp(x[3]);
    double top = *std::max_element(x.begin() + 4, x.end()), norm = 0.0;
    for (size_t n = 0; n <= FIT_COPIES; n++) {
        model.pi.push_back(std::exp(x[4 + n] - top));
        norm += model.pi.back();
    }
    for (double& p : model.pi) {
        p /= norm;
    }
    return model;
}

//------------------------------------------------------------------------------

std::vector<double> count_table(const CountModel& model, const std::vector<double>& histogram) {
    double total = 0.0;
    for (double h : histogram) {
        total += h;
    }
    std::vector<double> table((CountModel::MAX_COUNT + 1) * TABLE_WIDTH);
    for (size_t y = 0; y <= CountModel::MAX_COUNT; y++) {
        double background = std::log(std::max(histogram[y] / total, 1e-300));
        for (size_t m = 0; m < TABLE_WIDTH; m++) {
            double t = model.log_probability(y, m);
            table[y * TABLE_WIDTH + m] = log_add(std::log1p(-CONTAMINATION) + t, std::log(CONTAMINATION) + background);
        }
    }
    return table;
}

//------------------------------------------------------------------------------

namespace {

// Sequences with identical kmer presence and copies, numbered in the order of their
// first sequence.
struct SequenceGroups {
    // For each group, its first sequence.
    std::vector<size_t> representative;

    // For each group, its sequences.
    std::vector<std::vector<size_t>> members;

    size_t size() const { return this->representative.size(); }
};

// Group the usable sequences by kmer presence and copies.
SequenceGroups group_sequences(
    const std::vector<std::vector<std::uint64_t>>& rows,
    const std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>& copies,
    const std::vector<bool>& usable
) {
    SequenceGroups groups;
    std::unordered_map<std::uint64_t, std::vector<size_t>> by_hash;
    for (size_t i = 0; i < rows.size(); i++) {
        if (!usable[i]) {
            continue;
        }
        std::uint64_t h = 0;
        for (std::uint64_t word : rows[i]) {
            h = mix64(h ^ word);
        }
        for (const auto& c : copies[i]) {
            h = mix64(h ^ ((std::uint64_t(c.first) << 8) | c.second));
        }
        std::vector<size_t>& bucket = by_hash[h];
        size_t group = groups.size();
        for (size_t candidate : bucket) {
            size_t rep = groups.representative[candidate];
            if (rows[rep] == rows[i] && copies[rep] == copies[i]) {
                group = candidate;
                break;
            }
        }
        if (group == groups.size()) {
            groups.representative.push_back(i);
            groups.members.emplace_back();
            bucket.push_back(group);
        }
        groups.members[group].push_back(i);
    }
    return groups;
}

// The scores of all pairs of groups with copies clipped to one. The score of a pair
// (g, h) is base + single[g] + single[h] + shared(g, h), where `base` is the
// log-likelihood with no copies, `single[g]` sums the first differences of the
// log-likelihood in copies over the kmers g carries, and shared(g, h) sums the second
// differences over the kmers both carry. A kmer carried by most groups is counted
// through the groups without it, using x_g x_h = 1 - y_g - y_h + y_g y_h for y = 1 - x:
// it adds to `common` and to `missing` of the groups without it.
class ClippedScores {
public:
    ClippedScores(
        const std::vector<std::vector<std::uint64_t>>& rows, const SequenceGroups& groups,
        const std::vector<size_t>& counts, const std::vector<double>& table
    ) : num_groups(groups.size()), single(groups.size(), 0.0), shared(groups.size() * groups.size(), 0.0), missing(groups.size(), 0.0) {
        // The first and second differences of the log-likelihood in copies for each
        // kmer. Ignored kmers contribute nothing.
        size_t J = counts.size(), words = (J + 63) / 64;
        std::vector<double> first(J, 0.0), second(J, 0.0);
        for (size_t j = 0; j < J; j++) {
            if (counts[j] == IGNORED_COUNT) {
                continue;
            }
            const double* t = table.data() + counts[j] * TABLE_WIDTH;
            this->base += t[0];
            first[j] = t[1] - t[0];
            second[j] = t[2] - 2.0 * t[1] + t[0];
        }

        std::vector<size_t> carriers, others;
        for (size_t w = 0; w < words; w++) {
            size_t limit = std::min(size_t(64), J - 64 * w);
            for (size_t b = 0; b < limit; b++) {
                size_t j = 64 * w + b;
                std::uint64_t bit = std::uint64_t(1) << b;
                carriers.clear(); others.clear();
                for (size_t g = 0; g < this->num_groups; g++) {
                    if (rows[groups.representative[g]][w] & bit) {
                        carriers.push_back(g);
                    } else {
                        others.push_back(g);
                    }
                }
                this->add_kmer(carriers, others, first[j], second[j]);
            }
        }
    }

    // Return the score of the pair of groups.
    double operator()(size_t g, size_t h) const {
        size_t a = std::min(g, h), b = std::max(g, h);
        return this->base + this->single[g] + this->single[h] + this->shared[a * this->num_groups + b]
            + this->common - this->missing[g] - this->missing[h];
    }

private:
    size_t num_groups;
    double base = 0.0, common = 0.0;
    std::vector<double> single, shared, missing;

    // Add a kmer with the given differences, carried by the carriers.
    void add_kmer(const std::vector<size_t>& carriers, const std::vector<size_t>& others, double first, double second) {
        for (size_t g : carriers) {
            this->single[g] += first;
        }
        const std::vector<size_t>& listed = (carriers.size() <= others.size() ? carriers : others);
        if (&listed == &others) {
            this->common += second;
            for (size_t g : others) {
                this->missing[g] += second;
            }
        }
        for (size_t x = 0; x < listed.size(); x++) {
            for (size_t y = x; y < listed.size(); y++) {
                this->shared[listed[x] * this->num_groups + listed[y]] += second;
            }
        }
    }
};

// Return the first `CANDIDATES` groups to appear in the best `CANDIDATE_PAIRS` pairs.
// Ties are broken by the order of the pairs (g, h) with g <= h.
std::vector<size_t> choose_candidates(const ClippedScores& clipped, size_t num_groups) {
    std::vector<std::pair<double, size_t>> ranked;
    std::vector<std::pair<size_t, size_t>> pair_of;
    for (size_t g = 0; g < num_groups; g++) {
        for (size_t h = g; h < num_groups; h++) {
            ranked.emplace_back(clipped(g, h), pair_of.size());
            pair_of.emplace_back(g, h);
        }
    }
    size_t top = std::min(CANDIDATE_PAIRS, ranked.size());
    std::partial_sort(ranked.begin(), ranked.begin() + top, ranked.end(),
        [](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {
            return (a.first > b.first || (a.first == b.first && a.second < b.second));
        });
    std::vector<size_t> candidates;
    std::vector<bool> is_candidate(num_groups, false);
    for (size_t k = 0; k < top; k++) {
        for (size_t g : { pair_of[ranked[k].second].first, pair_of[ranked[k].second].second }) {
            if (!is_candidate[g]) {
                is_candidate[g] = true;
                candidates.push_back(g);
            }
        }
    }
    if (candidates.size() > CANDIDATES) {
        candidates.resize(CANDIDATES);
    }
    return candidates;
}

// Correct the K x K candidate pair scores for the kmers that some candidate has more
// than once: replace the term for the clipped copies with the term for the copies.
void add_copy_corrections(
    std::vector<double>& scores, const std::vector<size_t>& representatives,
    const std::vector<std::vector<std::uint64_t>>& rows,
    const std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>& copies,
    const std::vector<size_t>& counts, const std::vector<double>& table
) {
    size_t K = representatives.size();
    std::vector<std::uint32_t> multi;
    for (size_t rep : representatives) {
        for (const auto& c : copies[rep]) {
            if (counts[c.first] != IGNORED_COUNT) {
                multi.push_back(c.first);
            }
        }
    }
    std::sort(multi.begin(), multi.end());
    multi.erase(std::unique(multi.begin(), multi.end()), multi.end());
    if (multi.empty()) {
        return;
    }

    // The copies of those kmers in each candidate.
    std::vector<std::vector<std::uint32_t>> X(K, std::vector<std::uint32_t>(multi.size(), 0));
    for (size_t i = 0; i < K; i++) {
        size_t rep = representatives[i];
        for (size_t c = 0; c < multi.size(); c++) {
            X[i][c] = (rows[rep][multi[c] / 64] >> (multi[c] % 64)) & 1;
        }
        for (const auto& entry : copies[rep]) {
            auto iter = std::lower_bound(multi.begin(), multi.end(), entry.first);
            if (iter != multi.end() && *iter == entry.first) {
                X[i][iter - multi.begin()] = entry.second;
            }
        }
    }

    for (size_t i = 0; i < K; i++) {
        for (size_t k = 0; k < K; k++) {
            double correction = 0.0;
            for (size_t c = 0; c < multi.size(); c++) {
                const double* t = table.data() + counts[multi[c]] * TABLE_WIDTH;
                size_t m = std::min(size_t(X[i][c] + X[k][c]), PAIR_COPIES);
                size_t clipped_m = std::min(X[i][c], std::uint32_t(1)) + std::min(X[k][c], std::uint32_t(1));
                correction += t[m] - t[clipped_m];
            }
            scores[i * K + k] += correction;
        }
    }
}

} // anonymous namespace

Candidates score_pairs(
    const std::vector<std::vector<std::uint64_t>>& rows,
    const std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>& copies,
    const std::vector<bool>& usable,
    const std::vector<size_t>& counts,
    const std::vector<double>& table
) {
    SequenceGroups groups = group_sequences(rows, copies, usable);
    ClippedScores clipped(rows, groups, counts, table);
    std::vector<size_t> candidates = choose_candidates(clipped, groups.size());

    // Score the candidate pairs with the copies.
    Candidates result;
    size_t K = candidates.size();
    result.scores.resize(K * K);
    std::vector<size_t> representatives;
    for (size_t i = 0; i < K; i++) {
        representatives.push_back(groups.representative[candidates[i]]);
        for (size_t k = 0; k < K; k++) {
            result.scores[i * K + k] = clipped(candidates[i], candidates[k]);
        }
    }
    add_copy_corrections(result.scores, representatives, rows, copies, counts, table);

    for (size_t g : candidates) {
        result.members.push_back(groups.members[g]);
    }
    return result;
}

//------------------------------------------------------------------------------

std::vector<std::pair<size_t, size_t>> viterbi(const std::vector<Candidates>& subchains, double tau, double rho) {
    double scale = 1.0 / tau;
    size_t K = subchains.front().size();
    std::vector<double> delta(K * K);
    for (size_t i = 0; i < K * K; i++) {
        delta[i] = subchains.front().scores[i] * scale;
    }

    // Each member of the pair moves independently, so a step is done in two stages.
    // back[i] = (first, second): for the pair (a', b') in subchain i, second[a' * K' + b']
    // is the best b in subchain i - 1, and first[a' * K + b] is the best a given b.
    std::vector<std::pair<std::vector<std::uint16_t>, std::vector<std::uint16_t>>> back(subchains.size());
    for (size_t i = 1; i < subchains.size(); i++) {
        const Candidates& prev = subchains[i - 1];
        const Candidates& next = subchains[i];
        size_t Kp = prev.size(), Kn = next.size();
        double d = 0.5 * (prev.length + next.length);
        std::vector<double> A = transition(prev, next, -std::expm1(-rho * d));
        for (double& x : A) {
            x = std::log(std::max(x, 1e-300));
        }
        // Stage 1: M[a', b] = max_a delta[a, b] + log A[a, a'].
        std::vector<double> M(Kn * Kp);
        std::vector<std::uint16_t>& arg1 = back[i].first;
        arg1.resize(Kn * Kp);
        for (size_t an = 0; an < Kn; an++) {
            for (size_t b = 0; b < Kp; b++) {
                double best = -std::numeric_limits<double>::infinity();
                size_t best_a = 0;
                for (size_t a = 0; a < Kp; a++) {
                    double v = delta[a * Kp + b] + A[a * Kn + an];
                    if (v > best) {
                        best = v; best_a = a;
                    }
                }
                M[an * Kp + b] = best;
                arg1[an * Kp + b] = static_cast<std::uint16_t>(best_a);
            }
        }
        // Stage 2: delta'[a', b'] = max_b M[a', b] + log A[b, b'], plus the scores.
        std::vector<double> updated(Kn * Kn);
        std::vector<std::uint16_t>& arg2 = back[i].second;
        arg2.resize(Kn * Kn);
        for (size_t an = 0; an < Kn; an++) {
            for (size_t bn = 0; bn < Kn; bn++) {
                double best = -std::numeric_limits<double>::infinity();
                size_t best_b = 0;
                for (size_t b = 0; b < Kp; b++) {
                    double v = M[an * Kp + b] + A[b * Kn + bn];
                    if (v > best) {
                        best = v; best_b = b;
                    }
                }
                updated[an * Kn + bn] = best + next.scores[an * Kn + bn] * scale;
                arg2[an * Kn + bn] = static_cast<std::uint16_t>(best_b);
            }
        }
        delta = std::move(updated);
    }

    // Trace the path back.
    size_t last = subchains.back().size();
    size_t best = std::max_element(delta.begin(), delta.end()) - delta.begin();
    size_t a = best / last, b = best % last;
    std::vector<std::pair<size_t, size_t>> path(subchains.size());
    path.back() = { a, b };
    for (size_t i = subchains.size() - 1; i > 0; i--) {
        size_t Kp = subchains[i - 1].size(), Kn = subchains[i].size();
        size_t b_prev = back[i].second[a * Kn + b];
        size_t a_prev = back[i].first[a * Kp + b_prev];
        a = a_prev; b = b_prev;
        path[i - 1] = { a, b };
    }
    return path;
}

std::vector<std::pair<size_t, size_t>> choose_sequences(
    const std::vector<Candidates>& subchains, const std::vector<std::pair<size_t, size_t>>& path
) {
    std::vector<std::pair<size_t, size_t>> result;
    std::uint64_t previous[2] = { 0, 0 };
    bool has_previous[2] = { false, false };
    for (size_t subchain_id = 0; subchain_id < subchains.size(); subchain_id++) {
        const Candidates& subchain = subchains[subchain_id];
        size_t chosen[2] = { 0, 0 };
        for (size_t slot = 0; slot < 2; slot++) {
            size_t group = (slot == 0 ? path[subchain_id].first : path[subchain_id].second);
            const std::vector<size_t>& members = subchain.members[group];
            const std::vector<std::uint64_t>& ids = subchain.haplotypes[group];
            size_t choice = members.size();
            if (has_previous[slot]) {
                for (size_t m = 0; m < members.size(); m++) {
                    if (ids[m] == previous[slot]) {
                        choice = m;
                        break;
                    }
                }
            }
            if (choice == members.size()) {
                bool same_group = (slot == 1 && path[subchain_id].first == path[subchain_id].second);
                for (size_t m = 0; m < members.size(); m++) {
                    if (!(same_group && m == chosen[0])) {
                        choice = m;
                        break;
                    }
                }
                if (choice == members.size()) {
                    choice = 0;
                }
            }
            chosen[slot] = choice;
            previous[slot] = ids[choice];
            has_previous[slot] = true;
        }
        size_t first_group = path[subchain_id].first, second_group = path[subchain_id].second;
        result.emplace_back(subchain.members[first_group][chosen[0]], subchain.members[second_group][chosen[1]]);
    }
    return result;
}

//------------------------------------------------------------------------------

} // namespace pair_model

} // namespace vg
