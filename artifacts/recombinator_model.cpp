#include "recombinator_model.hpp"

#include "../src/statistics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <unordered_map>

namespace vg::pair_model {

//------------------------------------------------------------------------------

// Numerical constants.

constexpr size_t CountModel::MIN_COUNT;
constexpr size_t CountModel::MAX_COUNT;
constexpr size_t CountModel::FIT_COPIES;

//------------------------------------------------------------------------------

namespace {

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
    if (double mode = (r > 1.0 ? (r - 1.0) * m / r : 0.0); mode < double(from)) {
        // The terms decrease from `from` on, so we sum them relative to the first.
        double q = m / (r + m);
        double relative = 1.0;
        double sum = 1.0;
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

//------------------------------------------------------------------------------

// Minimization within a box by BFGS. Unconstrained coordinates z map to the interior
// of the box [lo, hi] by the logistic function x = lo + (hi - lo) / (1 + exp(-z)).

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

// Return a - b.
std::vector<double> difference(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Return the point of the box with unconstrained coordinates z.
std::vector<double> to_box(const std::vector<double>& z, const std::vector<double>& lo, const std::vector<double>& hi) {
    std::vector<double> x(z.size());
    for (size_t i = 0; i < z.size(); i++) {
        x[i] = lo[i] + (hi[i] - lo[i]) / (1.0 + std::exp(-z[i]));
    }
    return x;
}

// Return the unconstrained coordinates of a point in the interior of the box.
std::vector<double> from_box(const std::vector<double>& x, const std::vector<double>& lo, const std::vector<double>& hi) {
    std::vector<double> z(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        double u = (x[i] - lo[i]) / (hi[i] - lo[i]);
        z[i] = std::log(u / (1.0 - u));
    }
    return z;
}

// Return the gradient of `value` at z by central differences.
template<class Function>
std::vector<double> numerical_gradient(const Function& value, const std::vector<double>& z) {
    constexpr double h = 1e-5;
    std::vector<double> g(z.size());
    std::vector<double> t = z;
    for (size_t i = 0; i < z.size(); i++) {
        t[i] = z[i] + h;
        double above = value(t);
        t[i] = z[i] - h;
        double below = value(t);
        t[i] = z[i];
        g[i] = (above - below) / (2.0 * h);
    }
    return g;
}

// A point on a line search and the value of the function there.
struct LinePoint {
    std::vector<double> z;
    double value;
};

// Backtracking line search from z, where the function has value fz, along `step`, on
// which the directional derivative is `slope`. The step is halved until the value
// decreases by at least 1e-4 times the decrease the slope predicts. Return the point,
// or nothing if 80 halvings do not find one.
template<class Function>
std::optional<LinePoint> line_search(
    const Function& value, const std::vector<double>& z, double fz, const std::vector<double>& step, double slope
) {
    LinePoint next { std::vector<double>(z.size()), fz };
    double t = 1.0;
    for (size_t k = 0; k < 80; k++, t *= 0.5) {
        for (size_t i = 0; i < z.size(); i++) {
            next.z[i] = z[i] + t * step[i];
        }
        next.value = value(next.z);
        if (next.value <= fz + 1e-4 * t * slope) {
            return next;
        }
    }
    return std::nullopt;
}

// The BFGS approximation of the inverse Hessian of a function of d variables, which
// turns a gradient into a search direction. It starts as the identity, is scaled at
// the first update to the curvature seen in that step, and learns the curvature from
// each step and the change of the gradient over it.
class InverseHessian {
public:
    explicit InverseHessian(size_t d) : d(d) {
        this->reset();
    }

    // Forget the curvature and become the identity.
    void reset() {
        this->H.assign(this->d, std::vector<double>(this->d, 0.0));
        for (size_t i = 0; i < this->d; i++) {
            this->H[i][i] = 1.0;
        }
        this->scaled = false;
    }

    // Return the search direction for gradient g: -H g.
    std::vector<double> direction(const std::vector<double>& g) const {
        std::vector<double> step = matrix_multiply(this->H, g);
        for (double& x : step) {
            x = -x;
        }
        return step;
    }

    // Learn from step s, over which the gradient changed by y. A step along which the
    // function does not curve upward teaches nothing and is skipped.
    void update(const std::vector<double>& s, const std::vector<double>& y) {
        double sy = dot(s, y);
        if (!(sy > 1e-12 * std::sqrt(dot(s, s) * dot(y, y)))) {
            return;
        }
        if (!this->scaled) {
            double c = sy / dot(y, y);
            for (size_t i = 0; i < this->d; i++) {
                this->H[i][i] = c;
            }
            this->scaled = true;
        }
        double rho = 1.0 / sy;
        std::vector<double> Hy = matrix_multiply(this->H, y);
        double yHy = dot(y, Hy);
        for (size_t i = 0; i < this->d; i++) {
            for (size_t j = 0; j < this->d; j++) {
                this->H[i][j] += -rho * (s[i] * Hy[j] + Hy[i] * s[j]) + (rho * rho * yHy + rho) * s[i] * s[j];
            }
        }
    }

private:
    size_t d;
    std::vector<std::vector<double>> H;
    bool scaled = false;
};

// Take a step from z, where the function has value fz and gradient g: a line search
// along the direction the inverse Hessian gives for g, or along -g if that is not a
// descent direction, in which case the inverse Hessian forgets its curvature. Return
// the new point, or nothing if no step decreases the value.
template<class Function>
std::optional<LinePoint> descent_step(
    const Function& value, const std::vector<double>& z, double fz, const std::vector<double>& g, InverseHessian& hessian
) {
    std::vector<double> step = hessian.direction(g);
    double slope = dot(g, step);
    if (!(slope < 0.0)) {
        hessian.reset();
        for (size_t i = 0; i < step.size(); i++) {
            step[i] = -g[i];
        }
        slope = -dot(g, g);
        if (!(slope < 0.0)) {
            return std::nullopt;
        }
    }
    return line_search(value, z, fz, step, slope);
}

// Minimize `f` within the box [lo, hi], starting from `x0` in the interior of the
// box, by BFGS over the unconstrained coordinates with central-difference gradients.
// Stop after 5000 steps, or after five steps that do not decrease the value. Return
// the minimizer.
template<class Function>
std::vector<double> minimize_in_box(
    const Function& f, const std::vector<double>& x0, const std::vector<double>& lo, const std::vector<double>& hi
) {
    auto value = [&f, &lo, &hi](const std::vector<double>& z) {
        return f(to_box(z, lo, hi));
    };
    std::vector<double> z = from_box(x0, lo, hi);
    double fz = value(z);
    std::vector<double> g = numerical_gradient(value, z);
    InverseHessian hessian(x0.size());
    size_t steps = 0;
    size_t flat = 0;
    while (steps < 5000 && flat < 5) {
        std::optional<LinePoint> next = descent_step(value, z, fz, g, hessian);
        if (!next.has_value()) {
            break;
        }
        steps++;
        std::vector<double> g_next = numerical_gradient(value, next->z);
        hessian.update(difference(next->z, z), difference(g_next, g));
        double decrease = fz - next->value;
        z = std::move(next->z);
        g = std::move(g_next);
        fz = next->value;
        flat = (decrease <= 1e-13 * std::abs(fz) ? flat + 1 : 0);
    }
    return to_box(z, lo, hi);
}

//------------------------------------------------------------------------------

// The histogram the count model is fitted to. Counts above `censor` are only known to
// be above it.
struct CensoredHistogram {
    // Counts above this are censored.
    size_t censor;

    // The counts up to `censor` that some kmers have, in increasing order, and the
    // number of kmers with each.
    std::vector<size_t> counts;
    std::vector<double> kmers;

    // The number of kmers with a count above `censor`.
    double above = 0.0;

    CensoredHistogram(const std::vector<double>& histogram, size_t censor) : censor(censor) {
        for (size_t y = 0; y < histogram.size(); y++) {
            if (y > censor) {
                this->above += histogram[y];
            } else if (histogram[y] > 0.0) {
                this->counts.push_back(y);
                this->kmers.push_back(histogram[y]);
            }
        }
    }
};

// The fit works on the parameter vector x = { log mu, log eps, log r1, log r0,
// log w_0, ..., log w_FIT_COPIES }, where the mixture weight of n copies is
// w_n / (w_0 + ... + w_FIT_COPIES). The log-weights start at this index.
constexpr size_t LOG_WEIGHTS = 4;

// Return the count model with the parameters in x, without the mixture weights.
CountModel model_without_weights(const std::vector<double>& x) {
    return CountModel { .mu = std::exp(x[0]), .eps = std::exp(x[1]), .r1 = std::exp(x[2]), .r0 = std::exp(x[3]) };
}

// The largest log-weight in x, and the sum of the weights divided by its weight.
struct WeightScale {
    double top = 0.0;
    double norm = 0.0;

    explicit WeightScale(const std::vector<double>& x) : top(*std::max_element(x.begin() + LOG_WEIGHTS, x.end())) {
        for (size_t n = 0; n <= CountModel::FIT_COPIES; n++) {
            this->norm += std::exp(x[LOG_WEIGHTS + n] - this->top);
        }
    }
};

// Return the negative log-likelihood of the parameters x given the histogram.
double negative_log_likelihood(const std::vector<double>& x, const CensoredHistogram& histogram) {
    CountModel model = model_without_weights(x);
    WeightScale scale(x);
    // The log-probability of a count under the mixture, given its log-probability for
    // each copy number.
    auto mixture = [&x, &scale](const auto& given_copies) {
        double total = -std::numeric_limits<double>::infinity();
        for (size_t n = 0; n <= CountModel::FIT_COPIES; n++) {
            total = log_add(total, given_copies(n) + x[LOG_WEIGHTS + n] - scale.top - std::log(scale.norm));
        }
        return total;
    };

    double result = 0.0;
    for (size_t b = 0; b < histogram.counts.size(); b++) {
        size_t count = histogram.counts[b];
        result -= histogram.kmers[b] * mixture([&model, count](size_t n) {
            return model.log_probability(count, n);
        });
    }
    if (histogram.above != 0.0) {
        result -= histogram.above * mixture([&model, &histogram](size_t n) {
            auto [size, mean] = model.nb_parameters(n);
            return nb_log_tail(histogram.censor + 1, size, mean);
        });
    }
    return result;
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
    auto [size, mean] = this->nb_parameters(copies);
    if (count < MIN_COUNT) {
        double result = -std::numeric_limits<double>::infinity();
        for (size_t k = 0; k < MIN_COUNT; k++) {
            result = log_add(result, nb_log_pmf(double(k), size, mean));
        }
        return result;
    }
    if (count >= MAX_COUNT) {
        return nb_log_tail(MAX_COUNT, size, mean);
    }
    return nb_log_pmf(double(count), size, mean);
}

CountModel CountModel::fit(const std::vector<double>& histogram, double coverage) {
    auto censor = static_cast<size_t>(std::min(double(FIT_COPIES) * coverage / 2.0, double(MAX_COUNT - 1)));
    CensoredHistogram censored(histogram, censor);
    auto nll = [&censored](const std::vector<double>& x) {
        return negative_log_likelihood(x, censored);
    };
    std::vector<double> lo = { std::log(coverage / 8.0), std::log(1e-3), std::log(0.5), std::log(1e-3) };
    std::vector<double> hi = { std::log(coverage), std::log(10.0), std::log(1e4), std::log(1e3) };
    std::vector<double> initial_weights = { 0.5, 0.2, 0.2, 0.05, 0.02 };
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
                x0.push_back(std::log(n < initial_weights.size() ? initial_weights[n] : 0.01));
            }
            std::vector<double> candidate = minimize_in_box(nll, x0, lo, hi);
            if (double value = nll(candidate); value < best) {
                best = value;
                x = candidate;
            }
        }
    }

    CountModel model = model_without_weights(x);
    WeightScale scale(x);
    for (size_t n = 0; n <= FIT_COPIES; n++) {
        model.pi.push_back(std::exp(x[LOG_WEIGHTS + n] - scale.top) / scale.norm);
    }
    return model;
}

//------------------------------------------------------------------------------

CountTable count_table(const CountModel& model, const std::vector<double>& histogram) {
    double total = std::accumulate(histogram.begin(), histogram.end(), 0.0);
    CountTable table(CountModel::MAX_COUNT + 1);
    for (size_t y = 0; y <= CountModel::MAX_COUNT; y++) {
        double background = std::log(std::max(histogram[y] / total, 1e-300));
        for (size_t m = 0; m <= PAIR_COPIES; m++) {
            double t = model.log_probability(y, m);
            table[y][m] = log_add(std::log1p(-CONTAMINATION) + t, std::log(CONTAMINATION) + background);
        }
    }
    return table;
}

//------------------------------------------------------------------------------

namespace {

// Group the usable sequences by kmer presence and copies. Each group lists its sequences
// in increasing order, and the groups are numbered in the order of their first
// sequences. The first sequence of a group stands for it.
std::vector<std::vector<size_t>> group_sequences(
    const std::vector<std::vector<std::uint64_t>>& rows, const std::vector<KmerCopies>& copies, const std::vector<bool>& usable
) {
    std::vector<std::vector<size_t>> groups;
    std::map<std::pair<std::vector<std::uint64_t>, KmerCopies>, size_t> group_of;
    for (size_t i = 0; i < rows.size(); i++) {
        if (!usable[i]) {
            continue;
        }
        auto [iter, inserted] = group_of.try_emplace(std::make_pair(rows[i], copies[i]), groups.size());
        if (inserted) {
            groups.emplace_back();
        }
        groups[iter->second].push_back(i);
    }
    return groups;
}

// The scores of all pairs of groups with copies clipped to one. For a kmer, let x_g be
// 1 if group g carries it and 0 otherwise, and y_g = 1 - x_g. The log-likelihood of the
// kmer's count for the pair (g, h) is t[0] + (x_g + x_h) d1 + x_g x_h d2, where t[m] is
// its log-likelihood with m copies, d1 = t[1] - t[0], and d2 = t[2] - 2 t[1] + t[0].
// The score of a pair sums this over the kmers: `base` sums t[0], `single[g]` sums d1
// over the kmers g carries, and `shared[g, h]` sums d2 over the kmers both carry. A kmer
// carried by most groups is counted through the groups without it instead, using
// x_g x_h = 1 - y_g - y_h + y_g y_h: its d2 is added to `common`, to `missing[g]` for
// each group g without it, and to `shared[g, h]` for each pair of groups without it.
class ClippedScores {
public:
    ClippedScores(
        const std::vector<std::vector<std::uint64_t>>& rows, const std::vector<std::vector<size_t>>& groups,
        const std::vector<size_t>& counts, const CountTable& table
    ) : num_groups(groups.size()), single(groups.size(), 0.0), shared(groups.size() * (groups.size() + 1) / 2, 0.0), missing(groups.size(), 0.0) {
        // The first and second differences of the log-likelihood in copies for each
        // kmer. Ignored kmers contribute nothing.
        size_t J = counts.size();
        std::vector<double> first(J, 0.0);
        std::vector<double> second(J, 0.0);
        for (size_t j = 0; j < J; j++) {
            if (counts[j] == IGNORED_COUNT) {
                continue;
            }
            const auto& t = table[counts[j]];
            this->base += t[0];
            first[j] = t[1] - t[0];
            second[j] = t[2] - 2.0 * t[1] + t[0];
        }

        std::vector<size_t> carriers;
        std::vector<size_t> others;
        for (size_t j = 0; j < J; j++) {
            std::uint64_t bit = std::uint64_t(1) << (j % 64);
            carriers.clear();
            others.clear();
            for (size_t g = 0; g < this->num_groups; g++) {
                std::vector<size_t>& side = ((rows[groups[g].front()][j / 64] & bit) ? carriers : others);
                side.push_back(g);
            }
            this->add_kmer(carriers, others, first[j], second[j]);
        }
    }

    // Return the score of the pair of groups.
    double operator()(size_t g, size_t h) const {
        size_t a = std::min(g, h);
        size_t b = std::max(g, h);
        return this->base + this->single[g] + this->single[h] + this->shared[this->pair_index(a, b)]
            + this->common - this->missing[g] - this->missing[h];
    }

private:
    size_t num_groups;
    double base = 0.0;
    double common = 0.0;
    std::vector<double> single;
    // The pairs (g, h) with g <= h, ordered by g and then by h.
    std::vector<double> shared;
    std::vector<double> missing;

    // Return the index of the pair (g, h) with g <= h in `shared`.
    size_t pair_index(size_t g, size_t h) const {
        return g * (2 * this->num_groups - g + 1) / 2 + (h - g);
    }

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
        // The groups are listed in increasing order.
        for (size_t x = 0; x < listed.size(); x++) {
            for (size_t y = x; y < listed.size(); y++) {
                this->shared[this->pair_index(listed[x], listed[y])] += second;
            }
        }
    }
};

// Return the first `CANDIDATES` groups to appear in the best `CANDIDATE_PAIRS` pairs.
// Ties are broken by the order of the pairs (g, h) with g <= h.
std::vector<size_t> choose_candidates(const ClippedScores& clipped, size_t num_groups) {
    // The best pairs as (score, (g, h)), kept in a heap whose top is the worst of them.
    // The order is total, so the best pairs and their order are the same as with sorting
    // all pairs.
    using ranked_pair = std::pair<double, std::pair<size_t, size_t>>;
    auto better = [](const ranked_pair& a, const ranked_pair& b) {
        return (a.first > b.first || (a.first == b.first && a.second < b.second));
    };
    std::vector<ranked_pair> best;
    best.reserve(CANDIDATE_PAIRS);
    for (size_t g = 0; g < num_groups; g++) {
        for (size_t h = g; h < num_groups; h++) {
            ranked_pair pair(clipped(g, h), { g, h });
            if (best.size() < CANDIDATE_PAIRS) {
                best.push_back(pair);
                std::ranges::push_heap(best, better);
            } else if (better(pair, best.front())) {
                std::ranges::pop_heap(best, better);
                best.back() = pair;
                std::ranges::push_heap(best, better);
            }
        }
    }
    std::ranges::sort(best, better);
    std::vector<size_t> candidates;
    std::vector<bool> is_candidate(num_groups, false);
    for (const auto& [score, pair] : best) {
        auto [g, h] = pair;
        for (size_t group : { g, h }) {
            if (!is_candidate[group]) {
                is_candidate[group] = true;
                candidates.push_back(group);
            }
        }
    }
    if (candidates.size() > CANDIDATES) {
        candidates.resize(CANDIDATES);
    }
    return candidates;
}

// Correct the pair scores of the candidates for the kmers that some candidate has more
// than once: replace the term for the clipped copies with the term for the copies.
void add_copy_corrections(
    Candidates& candidates, const std::vector<std::vector<std::uint64_t>>& rows, const std::vector<KmerCopies>& copies,
    const std::vector<size_t>& counts, const CountTable& table
) {
    size_t K = candidates.size();
    // The kmers, other than ignored ones, that some candidate has more than once, sorted
    // and distinct.
    std::vector<std::uint32_t> multi;
    for (const std::vector<size_t>& members : candidates.members) {
        for (const auto& [kmer, value] : copies[members.front()]) {
            if (counts[kmer] != IGNORED_COUNT) {
                multi.push_back(kmer);
            }
        }
    }
    gbwt::removeDuplicates(multi, false);
    if (multi.empty()) {
        return;
    }

    // The copies of those kmers in each candidate.
    std::vector<std::vector<std::uint32_t>> X(K, std::vector<std::uint32_t>(multi.size(), 0));
    for (size_t i = 0; i < K; i++) {
        size_t rep = candidates.members[i].front();
        for (size_t c = 0; c < multi.size(); c++) {
            X[i][c] = (rows[rep][multi[c] / 64] >> (multi[c] % 64)) & 1;
        }
        for (const auto& [kmer, value] : copies[rep]) {
            auto iter = std::ranges::lower_bound(multi, kmer);
            if (iter != multi.end() && *iter == kmer) {
                X[i][iter - multi.begin()] = value;
            }
        }
    }

    for (size_t i = 0; i < K; i++) {
        for (size_t k = 0; k < K; k++) {
            double correction = 0.0;
            for (size_t c = 0; c < multi.size(); c++) {
                const auto& t = table[counts[multi[c]]];
                size_t m = std::min(size_t(X[i][c] + X[k][c]), PAIR_COPIES);
                size_t clipped_m = std::min(X[i][c], std::uint32_t(1)) + std::min(X[k][c], std::uint32_t(1));
                correction += t[m] - t[clipped_m];
            }
            candidates.scores[i * K + k] += correction;
        }
    }
}

} // anonymous namespace

Candidates score_pairs(
    const std::vector<std::vector<std::uint64_t>>& rows,
    const std::vector<KmerCopies>& copies,
    const std::vector<std::uint64_t>& haplotype_ids,
    const std::vector<bool>& usable,
    const std::vector<size_t>& counts,
    const CountTable& table,
    double length
) {
    std::vector<std::vector<size_t>> groups = group_sequences(rows, copies, usable);
    ClippedScores clipped(rows, groups, counts, table);
    std::vector<size_t> candidates = choose_candidates(clipped, groups.size());

    // Score the candidate pairs with the copies.
    Candidates result;
    size_t K = candidates.size();
    result.scores.resize(K * K);
    for (size_t i = 0; i < K; i++) {
        for (size_t k = 0; k < K; k++) {
            result.scores[i * K + k] = clipped(candidates[i], candidates[k]);
        }
        result.members.push_back(groups[candidates[i]]);
        std::vector<std::uint64_t>& haplotypes = result.haplotypes.emplace_back();
        for (size_t sequence : result.members.back()) {
            haplotypes.push_back(haplotype_ids[sequence]);
        }
    }
    result.length = length;
    add_copy_corrections(result, rows, copies, counts, table);
    return result;
}

//------------------------------------------------------------------------------

namespace {

// Return the log-probabilities (Ka x Kb) that a member of the pair moves from each group
// of `a` to each group of `b`, when it moves to a group chosen uniformly at random with
// probability r and otherwise keeps its haplotype.
std::vector<double> log_transitions(const Candidates& a, const Candidates& b, double r) {
    size_t Ka = a.size();
    size_t Kb = b.size();
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
            if (auto iter = where.find(x); iter != where.end()) {
                A[ga * Kb + iter->second] += share;
            }
        }
        double kept = std::accumulate(A.begin() + ga * Kb, A.begin() + (ga + 1) * Kb, 0.0);
        double leak = 1.0 - kept;
        for (size_t gb = 0; gb < Kb; gb++) {
            A[ga * Kb + gb] = (1.0 - r) * (A[ga * Kb + gb] + leak / double(Kb)) + r / double(Kb);
        }
    }
    for (double& x : A) {
        x = std::log(std::max(x, 1e-300));
    }
    return A;
}

// Return the largest of value(k) for k < K and the first k that has it. If no value is
// larger than -inf, return (-inf, 0).
template<class Value>
std::pair<double, size_t> best_of(size_t K, const Value& value) {
    double best = -std::numeric_limits<double>::infinity();
    size_t best_k = 0;
    for (size_t k = 0; k < K; k++) {
        if (double v = value(k); v > best) {
            best = v;
            best_k = k;
        }
    }
    return std::make_pair(best, best_k);
}

// A step of the Viterbi algorithm moves the two members of the pair one at a time,
// from the Kp groups of the previous subchain to the Kn groups of the next one. The
// log transitions A are as in `log_transitions()`.

// Move one member of each pair to the next subchain. `scores` has Kp rows of
// X = scores.size() / Kp: row p is the group of the member that moves, and column x is the
// group of the member that stays. Return result[x * Kn + n] = max_p scores[p * X + x] +
// A[p * Kn + n], where n is the new group of the moving member, and store the best p at
// from[x * Kn + n]. The member that stays comes first in the result, so a call on
// delta[a, b] moves a and gives M[b, a'], and a call on M moves b and gives delta[a', b'].
std::vector<double> move_member(
    const std::vector<double>& scores, const std::vector<double>& A, size_t Kp, size_t Kn, std::vector<std::uint16_t>& from
) {
    size_t X = scores.size() / Kp;
    std::vector<double> result(X * Kn);
    from.resize(X * Kn);
    for (size_t x = 0; x < X; x++) {
        for (size_t n = 0; n < Kn; n++) {
            auto [best, best_p] = best_of(Kp, [&scores, &A, X, Kn, x, n](size_t p) {
                return scores[p * X + x] + A[p * Kn + n];
            });
            result[x * Kn + n] = best;
            from[x * Kn + n] = static_cast<std::uint16_t>(best_p);
        }
    }
    return result;
}

// The best previous pair of each pair in a subchain: for the pair (a', b'), the best
// previous b is `second[a' * Kn + b']`, and the best previous a given that b is
// `first[b * Kn + a']`.
struct Backpointers {
    std::vector<std::uint16_t> first;
    std::vector<std::uint16_t> second;
};

} // anonymous namespace

std::vector<std::pair<size_t, size_t>> viterbi(const std::vector<Candidates>& subchains, double tau, double rho) {
    if (std::ranges::any_of(subchains, [](const Candidates& subchain) { return subchain.size() == 0; })) {
        throw std::invalid_argument("pair_model::viterbi(): a subchain has no candidates");
    }
    if (subchains.empty()) {
        return {};
    }
    double scale = 1.0 / tau;
    std::vector<double> delta(subchains.front().scores.size());
    for (size_t i = 0; i < delta.size(); i++) {
        delta[i] = subchains.front().scores[i] * scale;
    }

    std::vector<Backpointers> back(subchains.size());
    for (size_t i = 1; i < subchains.size(); i++) {
        const Candidates& prev = subchains[i - 1];
        const Candidates& next = subchains[i];
        double d = 0.5 * (prev.length + next.length);
        std::vector<double> A = log_transitions(prev, next, -std::expm1(-rho * d));
        std::vector<double> M = move_member(delta, A, prev.size(), next.size(), back[i].first);
        delta = move_member(M, A, prev.size(), next.size(), back[i].second);
        for (size_t j = 0; j < delta.size(); j++) {
            delta[j] += next.scores[j] * scale;
        }
    }

    // Trace the path back.
    size_t last = subchains.back().size();
    auto best = static_cast<size_t>(std::ranges::max_element(delta) - delta.begin());
    size_t a = best / last;
    size_t b = best % last;
    std::vector<std::pair<size_t, size_t>> path(subchains.size());
    path.back() = { a, b };
    for (size_t i = subchains.size() - 1; i > 0; i--) {
        size_t Kn = subchains[i].size();
        b = back[i].second[a * Kn + b];
        a = back[i].first[b * Kn + a];
        path[i - 1] = { a, b };
    }
    return path;
}

namespace {

// Return the member of a group to take, given the haplotypes of its members: the
// member of haplotype `previous`, if there is one, and otherwise the first member
// other than `avoid`, or the first member if there is no other.
size_t choose_member(
    const std::vector<std::uint64_t>& haplotypes, std::optional<std::uint64_t> previous, std::optional<size_t> avoid
) {
    if (previous.has_value()) {
        if (auto iter = std::ranges::find(haplotypes, *previous); iter != haplotypes.end()) {
            return iter - haplotypes.begin();
        }
    }
    if (avoid == std::optional<size_t>(0) && haplotypes.size() > 1) {
        return 1;
    }
    return 0;
}

} // anonymous namespace

std::vector<std::pair<size_t, size_t>> choose_sequences(
    const std::vector<Candidates>& subchains, const std::vector<std::pair<size_t, size_t>>& path
) {
    std::vector<std::pair<size_t, size_t>> result;
    std::optional<std::uint64_t> previous_first;
    std::optional<std::uint64_t> previous_second;
    for (size_t subchain_id = 0; subchain_id < subchains.size(); subchain_id++) {
        const Candidates& subchain = subchains[subchain_id];
        auto [first_group, second_group] = path[subchain_id];
        size_t first = choose_member(subchain.haplotypes[first_group], previous_first, std::nullopt);
        std::optional<size_t> avoid;
        if (second_group == first_group) {
            avoid = first;
        }
        size_t second = choose_member(subchain.haplotypes[second_group], previous_second, avoid);
        previous_first = subchain.haplotypes[first_group][first];
        previous_second = subchain.haplotypes[second_group][second];
        result.emplace_back(subchain.members[first_group][first], subchain.members[second_group][second]);
    }
    return result;
}

//------------------------------------------------------------------------------

} // namespace vg::pair_model
