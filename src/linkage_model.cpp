#include <deque>
#include "linkage_model.hpp"

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>

namespace vg {


/// The distance between two adjacent sites of a chain, in bp, once for each strand, since
/// `transition_apply` takes a switch probability per strand; the two are always equal.
/// Clamped at 1: `switch_probability` reads a gap of 0 as 1, and a negative one would wrap.
static inline std::pair<size_t, size_t> site_gap(const LinkageModel::Site& prev,
                                                const LinkageModel::Site& next) {
    if (prev.unpositioned && next.unpositioned) {
        // Both sites are unpositioned. Sites reach the model grouped by parent, chain, ploidy and
        // strand, so the two are in the same chain, and their positions are offsets along the same parent
        // traversal, whose difference is a distance.
        const size_t d = next.position > prev.position ? (size_t)(next.position - prev.position) : 1;
        return {d, d};
    }
    if (prev.unpositioned || next.unpositioned) {
        // One site is unpositioned and the other is not, so their positions cannot be differenced.
        // SIZE_MAX gives rho = 1, which makes the transition uniform: nothing is assumed about
        // linkage between the two.
        const size_t unknown = numeric_limits<size_t>::max();
        return {unknown, unknown};
    }
    const size_t ref = next.position > prev.position ? (size_t)(next.position - prev.position) : 1;
    return {ref, ref};
}

double LinkageModel::switch_probability(size_t gap) const {
    if (gap < 1) {
        gap = 1;
    }
    double rho = params.rho_min
                 + (1.0 - params.rho_min) * (1.0 - exp(-(double)gap / params.scale));
    rho = min(max(rho, 0.0), 1.0);
    if (params.weight == 1.0) {
        return rho;
    }
    if (params.weight <= 0.0) {
        // Uniform transitions: the chain forgets everything and the posterior is the emission.
        return 1.0;
    }
    return min(max(pow(rho, params.weight), 1e-12), 1.0);
}

/// The allele a panel haplotype carries at a site, or -1 if it carries none there. An allele index
/// out of range is treated as none rather than used.
static inline int allele_at(const LinkageModel::Site& site, size_t h, size_t n_hap) {
    if (h >= n_hap || h >= site.haplotype_allele.size()) {
        return -1;
    }
    int allele = site.haplotype_allele[h];
    return (allele >= 0 && (size_t)allele < site.num_alleles) ? allele : -1;
}

/// The VCF allele a compact allele was written as, or -1 where it was not written as one.
static inline int vcf_allele_of(const vector<int8_t>& allele_arena, size_t allele_offset,
                                size_t num_alleles, size_t compact) {
    if (compact >= num_alleles) {
        return -1;
    }
    size_t at = allele_offset + compact;
    return at < allele_arena.size() ? (int)allele_arena[at] : -1;
}

/// The VCF alleles of a settled compact pair, through the site's traversal-to-allele map.
///
/// Where the map gives no allele for one of the pair, because it was not supplied yet or because
/// the record wrote no ALT for that traversal, the called pair's alleles are used instead and
/// `fell_back` is set; nothing determines the order of that pair.
static inline void render_phase_pair(const vector<int8_t>& allele_arena, size_t allele_offset,
                                     size_t num_alleles, size_t c_first, size_t c_second,
                                     size_t called_i, size_t called_j,
                                     int* out_first, int* out_second, bool* fell_back) {
    const int v_first = vcf_allele_of(allele_arena, allele_offset, num_alleles, c_first);
    const int v_second = vcf_allele_of(allele_arena, allele_offset, num_alleles, c_second);
    if (v_first >= 0 && v_second >= 0) {
        *out_first = v_first;
        *out_second = v_second;
        *fell_back = false;
        return;
    }
    *out_first = vcf_allele_of(allele_arena, allele_offset, num_alleles, called_i);
    *out_second = vcf_allele_of(allele_arena, allele_offset, num_alleles, called_j);
    *fell_back = true;
}

/// What a parent's settled traversals imply about one of its child chains, given the chain's
/// crossing mask. The mask is indexed by candidate traversal, so it is tested against traversals,
/// never against compact allele indices.
LinkageCollector::Relation LinkageCollector::relate_to_parent(uint64_t crossing, int ta, int tb) {
    LinkageCollector::Relation r;
    if (ta < 0 || crossing == 0) {
        return r;   // nothing settled, or descent could not compute the mask
    }
    r.known = true;
    const bool first = ta < 64 && ((crossing >> ta) & 1);
    const bool second = tb >= 0 && tb < 64 && ((crossing >> tb) & 1);
    r.copies = (uint8_t)((int)first + (int)second);
    r.carrying_trav = r.copies == 1 ? (first ? ta : tb) : (r.copies == 2 ? -2 : -1);
    return r;
}

static inline int traversal_of(const vector<uint16_t>& trav_arena, size_t trav_offset,
                               size_t num_alleles, size_t compact) {
    if (compact >= num_alleles) {
        return -1;
    }
    size_t at = trav_offset + compact;
    return at < trav_arena.size() ? (int)trav_arena[at] : -1;
}

/// Relative P(reads | genotype implied by the state) for every ordered pair of panel haplotypes,
/// with the wildcard last, scaled so that the site's best genotype is 1.
static void build_emission(const LinkageModel::Site& site, size_t n_hap, double escape,
                           vector<double>& e, vector<double>& per_genotype) {
    size_t m = n_hap + 1;
    size_t n = site.num_alleles;

    // Genotype likelihoods, shifted so the best is exp(0) = 1.
    double best = -numeric_limits<double>::infinity();
    for (double v : site.genotype_ln_likelihood) {
        if (std::isfinite(v)) {
            best = max(best, v);
        }
    }
    per_genotype.assign(site.genotype_ln_likelihood.size(), 0.0);
    if (std::isfinite(best)) {
        for (size_t g = 0; g < site.genotype_ln_likelihood.size(); ++g) {
            double v = site.genotype_ln_likelihood[g];
            per_genotype[g] = std::isfinite(v) ? exp(v - best) : 0.0;
        }
    }

    // Mean over the other strand's alleles, for a state whose partner is unknown.
    vector<double> marginal(n, 0.0);
    double overall = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double acc = 0.0;
        for (size_t j = 0; j < n; ++j) {
            acc += per_genotype[LinkageModel::genotype_index(i, j)];
        }
        marginal[i] = n ? acc / (double)n : 0.0;
        overall += marginal[i];
    }
    overall = n ? overall / (double)n : 0.0;

    e.assign(m * m, 0.0);
    for (size_t a = 0; a < m; ++a) {
        int ai = allele_at(site, a, n_hap);
        for (size_t b = 0; b < m; ++b) {
            int bi = allele_at(site, b, n_hap);
            double value;
            if (ai >= 0 && bi >= 0) {
                value = per_genotype[LinkageModel::genotype_index((size_t)ai, (size_t)bi)];
            } else if (ai >= 0) {
                // Partner unknown, either the wildcard or a haplotype absent from this site.
                value = marginal[(size_t)ai] * escape;
            } else if (bi >= 0) {
                value = marginal[(size_t)bi] * escape;
            } else {
                value = overall * escape * escape;
            }
            e[a * m + b] = value;
        }
    }
}

/// One Li-Stephens step. T = (1-rho) I + (rho/m) 1, so the sum over previous ordered pairs
/// collapses to four terms and the step is O(m^2) rather than O(m^4).
void transition_apply(const vector<double>& in, size_t m,
                      double rho_a, double rho_b, vector<double>& out) {
    // One switch probability per strand, since the strands switch independently.
    double stay_a = 1.0 - rho_a, stay_b = 1.0 - rho_b;
    double jump_a = rho_a / (double)m, jump_b = rho_b / (double)m;
    vector<double> row(m, 0.0), col(m, 0.0);
    double total = 0.0;
    for (size_t a = 0; a < m; ++a) {
        for (size_t b = 0; b < m; ++b) {
            double v = in[a * m + b];
            row[a] += v;
            col[b] += v;
            total += v;
        }
    }
    out.assign(m * m, 0.0);
    for (size_t a = 0; a < m; ++a) {
        for (size_t b = 0; b < m; ++b) {
            // When the two rhos differ, the row and column terms have different coefficients, so
            // they are written separately.
            out[a * m + b] = stay_a * stay_b * in[a * m + b]
                             + stay_a * jump_b * row[a]
                             + jump_a * stay_b * col[b]
                             + jump_a * jump_b * total;
        }
    }
}

namespace {

const double NEG_INF = -numeric_limits<double>::infinity();

/// Best two values along one axis, with both indices, for the leave-one-out maxima the
/// max-product step needs. `stride` walks a row (1) or a column (m).
struct Top2 {
    double best = NEG_INF;
    double second = NEG_INF;
    size_t arg = 0;
    size_t arg2 = 0;
};

Top2 top2_of(const double* v, size_t n, size_t stride) {
    Top2 t;
    for (size_t i = 0; i < n; ++i) {
        double x = v[i * stride];
        if (x > t.best) {
            t.second = t.best;
            t.arg2 = t.arg;
            t.best = x;
            t.arg = i;
        } else if (x > t.second) {
            t.second = x;
            t.arg2 = i;
        }
    }
    return t;
}

/// One Li-Stephens max-product step, with backpointers.
///
/// A state is a pair (a, b): strand 0 copies haplotype a and strand 1 copies b. For strand x,
/// S_x is the log probability of staying on the same haplotype, and J_x the log probability of
/// moving to one particular other haplotype. `transition_apply` sums, and there the factorisation
/// collapses the pairwise loop into four terms. Maximising does *not* separate the same way,
/// because delta(a,b) couples the strands: max over (a,b) of delta(a,b) + f(a) + g(b) is not a
/// pair of independent 1-D maxima. But each strand's transition takes only two values, so the
/// reduction is by cases on which strands stayed:
///
///     delta'(a',b') = ln e(a',b') + max of
///         delta(a',b')                        + S_0 + S_1   both stayed
///         max_{b != b'} delta(a',b)           + S_0 + J_1   strand 0 stayed
///         max_{a != a'} delta(a,b')           + J_0 + S_1   strand 1 stayed
///         max_{a != a', b != b'} delta(a,b)   + J_0 + J_1   both moved
///
/// Every leave-one-out maximum comes from a top-2 along the relevant axis, so this stays O(m^2)
/// like the forward step rather than O(m^4).
///
/// In logs, unlike the forward pass: sum-product needs rescaling per site to avoid underflow,
/// max-product does not, and in logs the stay-or-jump choice is a comparison of sums.
void viterbi_step(const vector<double>& in, size_t m, double rho_a, double rho_b,
                  const vector<double>& emission,
                  vector<double>& out, vector<uint16_t>& back_a, vector<uint16_t>& back_b) {
    // Per strand, as in `transition_apply`. The four candidates below are already the four
    // stay/jump combinations, so each simply takes the coefficient belonging to its own axis; the
    // leave-one-out maxima are per-axis already and do not change at all.
    double stay_a = 1.0 - rho_a + rho_a / (double)m;
    double stay_b = 1.0 - rho_b + rho_b / (double)m;
    double jump_a = rho_a / (double)m;
    double jump_b = rho_b / (double)m;
    double S_0 = stay_a > 0.0 ? log(stay_a) : NEG_INF;
    double S_1 = stay_b > 0.0 ? log(stay_b) : NEG_INF;
    double J_0 = jump_a > 0.0 ? log(jump_a) : NEG_INF;
    double J_1 = jump_b > 0.0 ? log(jump_b) : NEG_INF;

    vector<Top2> rows(m), cols(m);
    for (size_t a = 0; a < m; ++a) {
        rows[a] = top2_of(&in[a * m], m, 1);
    }
    for (size_t b = 0; b < m; ++b) {
        cols[b] = top2_of(&in[b], m, m);
    }

    // rowExcl[a * m + bp] = max over b != bp of in[a][b], with its argument.
    vector<double> rowExcl(m * m);
    vector<uint16_t> rowExclArg(m * m);
    for (size_t a = 0; a < m; ++a) {
        for (size_t bp = 0; bp < m; ++bp) {
            bool hit = (rows[a].arg == bp);
            rowExcl[a * m + bp] = hit ? rows[a].second : rows[a].best;
            rowExclArg[a * m + bp] = (uint16_t)(hit ? rows[a].arg2 : rows[a].arg);
        }
    }

    out.assign(m * m, NEG_INF);
    back_a.assign(m * m, 0);
    back_b.assign(m * m, 0);

    for (size_t bp = 0; bp < m; ++bp) {
        // Both jumped: max over a != a' of rowExcl[a][bp]. One top-2 per arriving b'.
        Top2 both = top2_of(&rowExcl[bp], m, m);
        for (size_t ap = 0; ap < m; ++ap) {
            double e = emission[ap * m + bp];
            if (!(e > 0.0)) {
                // An impossible state stays impossible. This is also how a constraint is applied:
                // the caller zeroes the emission of every state that does not spell the required
                // genotype, so those simply never become reachable.
                continue;
            }
            double best = NEG_INF;
            size_t ba = 0, bb = 0;

            double c1 = in[ap * m + bp];
            if (c1 > NEG_INF) {
                double v = c1 + S_0 + S_1;
                if (v > best) { best = v; ba = ap; bb = bp; }
            }
            double c2 = rowExcl[ap * m + bp];
            if (c2 > NEG_INF) {
                double v = c2 + S_0 + J_1;
                if (v > best) { best = v; ba = ap; bb = rowExclArg[ap * m + bp]; }
            }
            bool chit = (cols[bp].arg == ap);
            double c3 = chit ? cols[bp].second : cols[bp].best;
            if (c3 > NEG_INF) {
                double v = c3 + J_0 + S_1;
                if (v > best) { best = v; ba = chit ? cols[bp].arg2 : cols[bp].arg; bb = bp; }
            }
            bool bhit = (both.arg == ap);
            double c4 = bhit ? both.second : both.best;
            if (c4 > NEG_INF) {
                double v = c4 + J_0 + J_1;
                if (v > best) {
                    best = v;
                    ba = bhit ? both.arg2 : both.arg;
                    bb = rowExclArg[ba * m + bp];
                }
            }
            if (best > NEG_INF) {
                out[ap * m + bp] = best + log(e);
                back_a[ap * m + bp] = (uint16_t)ba;
                back_b[ap * m + bp] = (uint16_t)bb;
            }
        }
    }
}

}   // anonymous namespace

void LinkageModel::window_posteriors(const vector<Site>& sites, size_t from, size_t to,
                                     vector<vector<double>>& out,
                                     const vector<double>* alpha_in,
                                     const vector<double>* beta_in) const {
    size_t n = to - from;
    if (n == 0) {
        return;
    }
    size_t n_hap = 0;
    for (size_t t = from; t < to; ++t) {
        n_hap = max(n_hap, sites[t].haplotype_allele.size());
    }
    size_t m = n_hap + 1;

    vector<vector<double>> emissions(n), per_genotype(n);
    for (size_t t = 0; t < n; ++t) {
        build_emission(sites[from + t], n_hap, params.escape, emissions[t], per_genotype[t]);
    }

    // Forward, rescaled every step. The scale factors are discarded: only the posterior is
    // wanted here, never the chain's total likelihood.
    vector<vector<double>> alpha(n);
    {
        // Uniform unless the caller handed in the message that reaches this window. Uniform says
        // "nothing is known before here" -- true for a whole chain, false for a segment cut from one.
        vector<double> a(m * m, 1.0 / (double)(m * m));
        if (alpha_in != nullptr && alpha_in->size() == m * m) {
            a = *alpha_in;
        }
        for (size_t k = 0; k < m * m; ++k) {
            a[k] *= emissions[0][k];
        }
        double s = 0.0;
        for (double v : a) {
            s += v;
        }
        if (s <= 0.0) {
            s = 1.0;
        }
        for (double& v : a) {
            v /= s;
        }
        alpha[0] = a;
    }
    for (size_t t = 1; t < n; ++t) {
        // The distance from the previous site; see `site_gap`.
        const std::pair<size_t, size_t> gap = site_gap(sites[from + t - 1], sites[from + t]);
        vector<double> moved;
        // One switch probability per strand.
        const double rho_a = switch_probability(gap.first);
        const double rho_b = switch_probability(gap.second);
        transition_apply(alpha[t - 1], m, rho_a, rho_b, moved);
        double s = 0.0;
        for (size_t k = 0; k < m * m; ++k) {
            moved[k] *= emissions[t][k];
            s += moved[k];
        }
        if (s <= 0.0) {
            s = 1.0;
        }
        for (double& v : moved) {
            v /= s;
        }
        alpha[t] = std::move(moved);
    }

    // Backward, combining as we go so only one beta is held at a time.
    vector<double> beta(m * m, 1.0);
    if (beta_in != nullptr && beta_in->size() == m * m) {
        beta = *beta_in;
    }
    for (size_t t = n; t-- > 0;) {
        const Site& site = sites[from + t];
        // A site's own exponent is bounded so that mass * multiplicity^(f - 1), with mass <= 1 and
        // multiplicity <= n_hap^2, stays under e^690 and sums over genotypes cannot overflow --
        // f 122 on 18 haplotypes, 60 on about 400.
        const double freq_prior =
            site.freq_prior >= 0.0
                ? std::min(site.freq_prior,
                           1.0 + 690.0 / std::log(std::max(4.0, (double)n_hap * (double)n_hap)))
                : params.freq_prior;
        vector<double>& post = out[from + t];
        post.assign(site.genotype_ln_likelihood.size(), 0.0);
        vector<size_t> multiplicity(site.genotype_ln_likelihood.size(), 0);

        // Two accumulators. `known` is mass from states where both strands name an allele, which
        // carries the multiplicity that the frequency prior rescales. `wild` is mass from states
        // with an unknown allele, which has no multiplicity.
        vector<double> known(site.genotype_ln_likelihood.size(), 0.0);
        vector<double> wild(site.genotype_ln_likelihood.size(), 0.0);
        size_t n_alleles = site.num_alleles;

        // A state pairing a panel haplotype with the wildcard occurs once for each haplotype that
        // carries the allele, so the number of carriers acts as a frequency prior on the
        // half-wildcard mass, and is rescaled like the multiplicity.
        vector<size_t> carriers(n_alleles, 0);
        for (size_t h = 0; h < n_hap && h < site.haplotype_allele.size(); ++h) {
            int allele = site.haplotype_allele[h];
            if (allele >= 0 && (size_t)allele < n_alleles) {
                carriers[(size_t)allele] += 1;
            }
        }

        for (size_t a = 0; a < m; ++a) {
            int ai = allele_at(site, a, n_hap);
            for (size_t b = 0; b < m; ++b) {
                int bi = allele_at(site, b, n_hap);
                double g = alpha[t][a * m + b] * beta[a * m + b];
                if (g <= 0.0) {
                    continue;
                }
                if (ai >= 0 && bi >= 0) {
                    size_t idx = genotype_index((size_t)ai, (size_t)bi);
                    known[idx] += g;
                    multiplicity[idx] += 1;
                    continue;
                }
                if (n_alleles == 0) {
                    continue;
                }
                // An unknown allele's mass is shared among the genotypes it could complete, in
                // proportion to their likelihoods rather than uniformly.
                if (ai >= 0 || bi >= 0) {
                    size_t k = (size_t)(ai >= 0 ? ai : bi);
                    double norm = 0.0;
                    for (size_t other = 0; other < n_alleles; ++other) {
                        norm += per_genotype[t][genotype_index(k, other)];
                    }
                    if (norm <= 0.0) {
                        continue;
                    }
                    double share = g;
                    if (k < carriers.size() && carriers[k] > 1) {
                        // At freq_prior 1 this changes nothing; above 1 it strengthens the prior.
                        share /= pow((double)carriers[k], 1.0 - freq_prior);
                    }
                    for (size_t other = 0; other < n_alleles; ++other) {
                        size_t idx = genotype_index(k, other);
                        wild[idx] += share * per_genotype[t][idx] / norm;
                    }
                } else {
                    double norm = 0.0;
                    for (size_t i = 0; i < n_alleles; ++i) {
                        for (size_t j = i; j < n_alleles; ++j) {
                            norm += per_genotype[t][genotype_index(i, j)];
                        }
                    }
                    if (norm <= 0.0) {
                        continue;
                    }
                    for (size_t i = 0; i < n_alleles; ++i) {
                        for (size_t j = i; j < n_alleles; ++j) {
                            size_t idx = genotype_index(i, j);
                            wild[idx] += g * per_genotype[t][idx] / norm;
                        }
                    }
                }
            }
        }

        // Rescale the known-known mass by multiplicity^(freq_prior - 1), the allele-frequency
        // prior.
        double total = 0.0;
        for (size_t idx = 0; idx < post.size(); ++idx) {
            double k = known[idx];
            if (multiplicity[idx] > 1) {
                // As for the carriers above.
                k /= pow((double)multiplicity[idx], 1.0 - freq_prior);
            }
            post[idx] = k + wild[idx];
            total += post[idx];
        }
        if (total > 0.0) {
            for (double& v : post) {
                v /= total;
            }
        } else {
            post.clear();
        }

        if (t == 0) {
            break;
        }
        // The distance from the previous site; see `site_gap`.
        const std::pair<size_t, size_t> gap = site_gap(sites[from + t - 1], sites[from + t]);
        vector<double> weighted(m * m, 0.0);
        for (size_t k = 0; k < m * m; ++k) {
            weighted[k] = beta[k] * emissions[t][k];
        }
        vector<double> next;
        const double rho_a = switch_probability(gap.first);
        const double rho_b = switch_probability(gap.second);
        transition_apply(weighted, m, rho_a, rho_b, next);
        double s = 0.0;
        for (double v : next) {
            s += v;
        }
        if (s <= 0.0) {
            s = 1.0;
        }
        for (double& v : next) {
            v /= s;
        }
        beta = std::move(next);
    }
}

/// Overlapping windows whose interiors are pasted together.
///
/// This works for a per-site quantity: a site's posterior hardly depends on where a window edge
/// fell, so windows can be decoded independently and their interiors kept. A path cannot be
/// decoded this way; see `windowed_path`.
///
/// `window(lo, hi, local)` fills `local` for `[lo, hi)`, indexed from 0 of the whole chain.
template <typename T, typename Window>
static void windowed_marginals(size_t n, size_t step, size_t margin, vector<T>& out,
                               Window window) {
    if (n <= step) {
        window(0, n, out);
        return;
    }
    for (size_t start = 0; start < n; start += step) {
        size_t lo = start > margin ? start - margin : 0;
        size_t hi = min(start + step + margin, n);
        vector<T> local(n);
        window(lo, hi, local);
        size_t keep_to = min(start + step, n);
        for (size_t t = start; t < keep_to; ++t) {
            out[t] = std::move(local[t]);
        }
    }
}

/// Overlapping windows chained by a pin.
///
/// Two windows decoded independently could choose different states where they meet, which would
/// put a spurious switch at every window boundary. So each window after the first is pinned: the
/// state at one index inside its leading margin is fixed to what the previous window chose there.
/// The pin is in the margin, where the previous window's choice was made with context on both
/// sides.
///
/// `window(lo, hi, pin_index, pin, local)` fills `local` for `[lo, hi)`, indexed from `lo`.
template <typename T, typename Window>
static void windowed_path(size_t n, size_t step, size_t margin, vector<T>& out, Window window) {
    bool have_pin = false;
    size_t pin_index = 0;
    T pin{};
    for (size_t start = 0; start < n; start += step) {
        size_t lo = start > margin ? start - margin : 0;
        size_t hi = min(start + step + margin, n);
        vector<T> local;
        window(lo, hi, have_pin ? pin_index : (size_t)-1, pin, local);
        size_t keep_to = min(start + step, n);
        for (size_t t = start; t < keep_to && t - lo < local.size(); ++t) {
            out[t] = local[t - lo];
        }
        if (keep_to == n) {
            break;
        }
        pin_index = keep_to - 1;
        pin = out[pin_index];
        have_pin = true;
    }
}

vector<vector<double>> LinkageModel::posteriors(const vector<Site>& sites, size_t ploidy,
                                                const vector<double>* alpha_in) const {
    vector<vector<double>> out(sites.size());
    if (sites.empty()) {
        return out;
    }
    windowed_marginals(sites.size(), max<size_t>(params.window, 1), params.margin, out,
                       [&](size_t lo, size_t hi, vector<vector<double>>& local) {
                           // Only the window that actually begins the chain gets the entering
                           // message; a later window's left end is an interior seam, and the
                           // margin is what carries the chain across it.
                           const vector<double>* enter = lo == 0 ? alpha_in : nullptr;
                           if (ploidy == 1) {
                               window_haploid_posteriors(sites, lo, hi, local, enter);
                           } else {
                               window_posteriors(sites, lo, hi, local, enter);
                           }
                       });
    return out;
}

vector<LinkageModel::Phase> LinkageModel::phasing(const vector<Site>& sites,
                                                  const vector<size_t>& constraint, size_t ploidy,
                                                  const vector<double>* alpha_in) const {
    vector<Phase> out(sites.size());
    if (sites.empty()) {
        return out;
    }
    size_t n = sites.size();
    if (ploidy == 1) {
        // One strand, so the path is over single haplotypes and `second` stays the wildcard. It
        // is decoded into its own vector because the pin here is one haplotype rather than a pair.
        vector<size_t> single(n, WILDCARD);
        windowed_path(n, max<size_t>(params.window, 1), params.margin, single,
                      [&](size_t lo, size_t hi, size_t pin_index, const size_t& pin,
                          vector<size_t>& local) {
                          // Only the window that begins the chain, as in `posteriors`: a later
                          // window's left end is an interior seam and the margin carries the chain
                          // across it.
                          window_haploid_phasing(sites, lo, hi, constraint, pin_index, pin, local,
                                                 lo == 0 ? alpha_in : nullptr);
                      });
        for (size_t t = 0; t < n; ++t) {
            out[t].first = single[t];
        }
        return out;
    }
    size_t n_hap = 0;
    for (const Site& s : sites) {
        n_hap = max(n_hap, s.haplotype_allele.size());
    }
    size_t m = n_hap + 1;
    if (m < 2) {
        return out;
    }

    windowed_path(n, max<size_t>(params.window, 1), params.margin, out,
                  [&](size_t lo, size_t hi, size_t pin_index, const Phase& pin,
                      vector<Phase>& local) {
                      window_phasing(sites, lo, hi, constraint, pin_index, pin, local);
                  });
    return out;
}

void LinkageModel::window_phasing(const vector<Site>& sites, size_t from, size_t to,
                                  const vector<size_t>& constraint,
                                  size_t pin_index, const Phase& pin,
                                  vector<Phase>& out) const {
    size_t n = to - from;
    out.assign(n, Phase{});
    if (n == 0) {
        return;
    }
    size_t n_hap = 0;
    for (size_t t = from; t < to; ++t) {
        n_hap = max(n_hap, sites[t].haplotype_allele.size());
    }
    size_t m = n_hap + 1;

    // Emissions, with the constraint folded in as zeroes. Zeroing rather than masking keeps the
    // step's own "impossible stays impossible" test doing double duty, and means a constrained
    // run and an unconstrained one differ only in this vector.
    vector<vector<double>> emissions(n), per_genotype(n);
    for (size_t t = 0; t < n; ++t) {
        const Site& site = sites[from + t];
        build_emission(site, n_hap, params.escape, emissions[t], per_genotype[t]);
        size_t want = (from + t) < constraint.size() ? constraint[from + t] : NO_CONSTRAINT;
        if (want == NO_CONSTRAINT) {
            continue;
        }
        // Decode the wanted genotype, because a state with one free strand has to be checked
        // against the individual alleles rather than against the pair.
        size_t wj = 0;
        while (genotype_index(0, wj + 1) <= want) {
            ++wj;
        }
        size_t wi = want - (wj * (wj + 1) / 2);

        for (size_t a = 0; a < m; ++a) {
            int ai = a < site.haplotype_allele.size() ? site.haplotype_allele[a] : -1;
            for (size_t b = 0; b < m; ++b) {
                int bi = b < site.haplotype_allele.size() ? site.haplotype_allele[b] : -1;
                // The wildcard, and a haplotype absent from this site, may carry any allele, which
                // keeps the constrained problem feasible where the panel cannot spell the call.
                // Only the free strand is unconstrained: a known haplotype carrying neither wanted
                // allele cannot be rescued by pairing it with the wildcard.
                bool ok;
                if (ai >= 0 && bi >= 0) {
                    ok = genotype_index((size_t)ai, (size_t)bi) == want;
                } else if (ai >= 0) {
                    ok = ((size_t)ai == wi || (size_t)ai == wj);
                } else if (bi >= 0) {
                    ok = ((size_t)bi == wi || (size_t)bi == wj);
                } else {
                    ok = true;
                }
                if (!ok) {
                    emissions[t][a * m + b] = 0.0;
                    continue;
                }
                // Every surviving state implies the same genotype, so all of them take that
                // genotype's likelihood. build_emission gave the free strands an average over
                // alleles, which here would let the path prefer the wildcard wherever the reads
                // disagree with the call. The escape factor stays, one per free strand, so a
                // genotype the panel can spell is still preferred.
                double e = want < per_genotype[t].size() ? per_genotype[t][want] : 0.0;
                if (ai < 0) {
                    e *= params.escape;
                }
                if (bi < 0) {
                    e *= params.escape;
                }
                emissions[t][a * m + b] = e;
            }
        }
    }

    // Per-site pins, for sites an earlier generation already settled. Same operation as the seam
    // pin below -- zero every state but one -- applied wherever the site asks for it. A settled
    // site's phase is already in the VCF, so letting the path re-orient it here would decode this
    // generation's sites in a frame nothing else uses.
    for (size_t t = 0; t < n; ++t) {
        const Site& site = sites[from + t];
        if (!site.pinned) {
            continue;
        }
        size_t pa = site.pin_first == WILDCARD ? n_hap : min(site.pin_first, n_hap);
        size_t pb = site.pin_second == WILDCARD ? n_hap : min(site.pin_second, n_hap);
        if (emissions[t][pa * m + pb] <= 0.0) {
            // The pinned pair cannot spell this site's constrained genotype, so pinning it would
            // leave the site with no reachable state. Leave it only constrained, and count it: in
            // a group whose only pinned site is its parent, a declined pin leaves the group free
            // to swap its strands relative to the parent.
            ++counters.pin_declined;
            continue;
        }
        ++counters.pin_applied;
        for (size_t a = 0; a < m; ++a) {
            for (size_t b = 0; b < m; ++b) {
                if (a != pa || b != pb) {
                    emissions[t][a * m + b] = 0.0;
                }
            }
        }
    }

    // Pin: every state but the pinned one becomes unreachable at that index. As for the per-site
    // pins above, a pin that is impossible under this window's constraint is skipped, leaving the
    // site only constrained. With margin 0 the pin index falls outside [from, to) and is skipped.
    if (pin_index != (size_t)-1 && pin_index >= from && pin_index < to) {
        size_t t = pin_index - from;
        size_t pa = pin.first == WILDCARD ? n_hap : min(pin.first, n_hap);
        size_t pb = pin.second == WILDCARD ? n_hap : min(pin.second, n_hap);
        double keep = emissions[t][pa * m + pb];
        if (keep > 0.0) {
            emissions[t].assign(m * m, 0.0);
            emissions[t][pa * m + pb] = keep;
        }
    }

    vector<double> delta(m * m, NEG_INF);
    for (size_t k = 0; k < m * m; ++k) {
        if (emissions[0][k] > 0.0) {
            delta[k] = -log((double)(m * m)) + log(emissions[0][k]);
        }
    }
    vector<vector<uint16_t>> back_a(n), back_b(n);
    for (size_t t = 1; t < n; ++t) {
        // The distance from the previous site; see `site_gap`.
        const std::pair<size_t, size_t> gap = site_gap(sites[from + t - 1], sites[from + t]);
        vector<double> next;
        const double rho_a = switch_probability(gap.first);
        const double rho_b = switch_probability(gap.second);
        viterbi_step(delta, m, rho_a, rho_b, emissions[t], next, back_a[t], back_b[t]);
        bool any = false;
        for (double v : next) {
            if (v > NEG_INF) { any = true; break; }
        }
        if (!any) {
            // No state survives: the constraint and the panel disagree beyond what the wildcard
            // can absorb. Restart the chain here rather than abandoning the window, so the rest
            // of it is still phased; the discontinuity is visible as a switch on both strands.
            for (size_t k = 0; k < m * m; ++k) {
                next[k] = emissions[t][k] > 0.0 ? log(emissions[t][k]) : NEG_INF;
                back_a[t][k] = (uint16_t)(k / m);
                back_b[t][k] = (uint16_t)(k % m);
            }
        }
        delta = std::move(next);
    }

    size_t best = m * m;
    double best_v = NEG_INF;
    for (size_t k = 0; k < m * m; ++k) {
        if (delta[k] > best_v) {
            best_v = delta[k];
            best = k;
        }
    }
    if (best == m * m) {
        return;
    }
    size_t a = best / m, b = best % m;
    for (size_t t = n; t-- > 0;) {
        out[t].first = (a == n_hap) ? WILDCARD : a;
        out[t].second = (b == n_hap) ? WILDCARD : b;
        if (t > 0) {
            size_t pa = back_a[t][a * m + b];
            size_t pb = back_b[t][a * m + b];
            a = pa;
            b = pb;
        }
    }
}


//------------------------------------------------------------------------------
// Ploidy-1 chains
//
// The states are single panel haplotypes, so a genotype is one allele.

void LinkageModel::haploid_emission(const Site& site, size_t n_hap, vector<double>& e,
                                    vector<double>& per_allele) const {
    size_t m = n_hap + 1;
    size_t n = site.num_alleles;

    // Shift so the best allele is exp(0) = 1, as build_emission does for genotypes: it keeps the
    // numbers in range without changing any ratio.
    double best = -numeric_limits<double>::infinity();
    for (double v : site.genotype_ln_likelihood) {
        if (std::isfinite(v)) {
            best = max(best, v);
        }
    }
    per_allele.assign(n, 0.0);
    for (size_t a = 0; a < n && a < site.genotype_ln_likelihood.size(); ++a) {
        double v = site.genotype_ln_likelihood[a];
        per_allele[a] = std::isfinite(v) && std::isfinite(best) ? exp(v - best) : 0.0;
    }

    double overall = 0.0;
    for (double v : per_allele) {
        overall += v;
    }
    overall = n ? overall / (double)n : 0.0;

    e.assign(m, 0.0);
    for (size_t a = 0; a < m; ++a) {
        int ai = allele_at(site, a, n_hap);
        // The wildcard, and a haplotype absent from this site, carry an unknown allele: average
        // over the alleles and pay the escape penalty, exactly as the diploid emission does.
        e[a] = (ai >= 0) ? per_allele[(size_t)ai] : overall * params.escape;
    }
}

void LinkageModel::window_haploid_posteriors(const vector<Site>& sites, size_t from, size_t to,
                                             vector<vector<double>>& out,
                                             const vector<double>* alpha_in) const {
    size_t n = to - from;
    if (n == 0) {
        return;
    }
    size_t n_hap = 0;
    size_t n_alleles_max = 0;
    for (size_t t = from; t < to; ++t) {
        n_hap = max(n_hap, sites[t].haplotype_allele.size());
        n_alleles_max = max(n_alleles_max, sites[t].num_alleles);
    }
    size_t m = n_hap + 1;

    vector<vector<double>> emissions(n), per_allele(n);
    for (size_t t = 0; t < n; ++t) {
        haploid_emission(sites[from + t], n_hap, emissions[t], per_allele[t]);
    }

    // Forward, rescaled every step; the scale factors are discarded, as in the diploid pass.
    vector<vector<double>> alpha(n);
    {
        // Uniform unless the caller supplied the message entering the chain.
        vector<double> a(m, 1.0 / (double)m);
        if (alpha_in != nullptr && alpha_in->size() == m) {
            a = *alpha_in;
        }
        double sum = 0.0;
        for (size_t k = 0; k < m; ++k) {
            a[k] *= emissions[0][k];
            sum += a[k];
        }
        if (sum <= 0.0) {
            sum = 1.0;
        }
        for (double& v : a) {
            v /= sum;
        }
        alpha[0] = a;
    }
    for (size_t t = 1; t < n; ++t) {
        // The distance from the previous site; see `site_gap`. One strand, so one distance.
        const double rho = switch_probability(site_gap(sites[from + t - 1], sites[from + t]).first);
        double stay = 1.0 - rho;
        double jump = rho / (double)m;
        double total = 0.0;
        for (double v : alpha[t - 1]) {
            total += v;
        }
        vector<double> next(m, 0.0);
        double sum = 0.0;
        for (size_t k = 0; k < m; ++k) {
            next[k] = (stay * alpha[t - 1][k] + jump * total) * emissions[t][k];
            sum += next[k];
        }
        if (sum <= 0.0) {
            sum = 1.0;
        }
        for (double& v : next) {
            v /= sum;
        }
        alpha[t] = std::move(next);
    }

    // Backward, combined as we go so only one beta is held.
    vector<double> beta(m, 1.0);
    for (size_t t = n; t-- > 0;) {
        const Site& site = sites[from + t];
        const double freq_prior = site.freq_prior >= 0.0 ? site.freq_prior : params.freq_prior;
        vector<double>& post = out[from + t];
        post.assign(site.num_alleles, 0.0);

        // Allele multiplicity is the haploid analogue of the diploid pair count, and is divided
        // out the same way so that freq_prior = 0 means what it says.
        vector<size_t> carriers(site.num_alleles, 0);
        for (size_t h = 0; h < n_hap && h < site.haplotype_allele.size(); ++h) {
            int allele = site.haplotype_allele[h];
            if (allele >= 0 && (size_t)allele < site.num_alleles) {
                carriers[(size_t)allele] += 1;
            }
        }
        vector<double> known(site.num_alleles, 0.0), wild(site.num_alleles, 0.0);
        for (size_t a = 0; a < m; ++a) {
            double g = alpha[t][a] * beta[a];
            if (g <= 0.0) {
                continue;
            }
            int ai = allele_at(site, a, n_hap);
            if (ai >= 0) {
                known[(size_t)ai] += g;
                continue;
            }
            // An unknown allele's mass is shared in proportion to the likelihoods, as in the
            // diploid pass.
            double norm = 0.0;
            for (double v : per_allele[t]) {
                norm += v;
            }
            if (norm <= 0.0) {
                continue;
            }
            for (size_t k = 0; k < site.num_alleles; ++k) {
                wild[k] += g * per_allele[t][k] / norm;
            }
        }
        double total = 0.0;
        for (size_t k = 0; k < site.num_alleles; ++k) {
            double v = known[k];
            if (carriers[k] > 1) {
                v /= pow((double)carriers[k], 1.0 - freq_prior);
            }
            post[k] = v + wild[k];
            total += post[k];
        }
        if (total > 0.0) {
            for (double& v : post) {
                v /= total;
            }
        }

        if (t == 0) {
            break;
        }
        // The distance from the previous site; see `site_gap`. One strand, so one distance.
        const double rho = switch_probability(site_gap(sites[from + t - 1], sites[from + t]).first);
        double stay = 1.0 - rho;
        double jump = rho / (double)m;
        vector<double> weighted(m);
        double total_w = 0.0;
        for (size_t k = 0; k < m; ++k) {
            weighted[k] = beta[k] * emissions[t][k];
            total_w += weighted[k];
        }
        vector<double> prev(m, 0.0);
        double sum = 0.0;
        for (size_t k = 0; k < m; ++k) {
            prev[k] = stay * weighted[k] + jump * total_w;
            sum += prev[k];
        }
        if (sum <= 0.0) {
            sum = 1.0;
        }
        for (double& v : prev) {
            v /= sum;
        }
        beta = std::move(prev);
    }
}

void LinkageModel::window_haploid_phasing(const vector<Site>& sites, size_t from, size_t to,
                                          const vector<size_t>& constraint,
                                          size_t pin_index, size_t pin,
                                          vector<size_t>& out,
                                          const vector<double>* alpha_in) const {
    size_t n = to - from;
    out.assign(n, WILDCARD);
    if (n == 0) {
        return;
    }
    size_t n_hap = 0;
    for (size_t t = from; t < to; ++t) {
        n_hap = max(n_hap, sites[t].haplotype_allele.size());
    }
    size_t m = n_hap + 1;

    vector<vector<double>> emissions(n), per_allele(n);
    for (size_t t = 0; t < n; ++t) {
        const Site& site = sites[from + t];
        haploid_emission(site, n_hap, emissions[t], per_allele[t]);
        size_t want = (from + t) < constraint.size() ? constraint[from + t] : NO_CONSTRAINT;
        if (want == NO_CONSTRAINT) {
            continue;
        }
        // Constrain to states carrying the called allele. As in the diploid case, every surviving
        // state then implies the same call, so they all take that allele's likelihood; a free
        // strand keeps the escape penalty so the panel is still preferred where it can explain.
        double e = want < per_allele[t].size() ? per_allele[t][want] : 0.0;
        for (size_t a = 0; a < m; ++a) {
            int ai = allele_at(site, a, n_hap);
            if (ai >= 0 && (size_t)ai != want) {
                emissions[t][a] = 0.0;
            } else if (ai >= 0) {
                emissions[t][a] = e;
            } else {
                emissions[t][a] = e * params.escape;
            }
        }
    }
    // As in the diploid caller: test before zeroing, and skip the pin outright when the pinned
    // state conflicts with this window's constraints, rather than pinning the forbidden state.
    if (pin_index != (size_t)-1 && pin_index >= from && pin_index < to) {
        size_t t = pin_index - from;
        size_t pa = pin == WILDCARD ? n_hap : min(pin, n_hap);
        double keep = emissions[t][pa];
        if (keep > 0.0) {
            emissions[t].assign(m, 0.0);
            emissions[t][pa] = keep;
        }
    }

    const double NINF = -numeric_limits<double>::infinity();
    vector<double> delta(m, NINF);
    // The entering message, where there is one, multiplies the first site's emission. The
    // collector passes a point mass, which therefore chooses the first state. The message is
    // ignored when missing or of the wrong length for these states.
    const bool have_alpha = alpha_in != nullptr && alpha_in->size() == m;
    for (size_t k = 0; k < m; ++k) {
        if (emissions[0][k] <= 0.0) {
            continue;
        }
        const double prior = have_alpha ? (*alpha_in)[k] : 1.0 / (double)m;
        if (prior > 0.0) {
            delta[k] = log(prior) + log(emissions[0][k]);
        }
    }
    // Every state the message allows is also forbidden by the emission, so the message and the
    // reads disagree outright. Fall back to the reads: a decode with no live state returns
    // WILDCARD everywhere, which reports "the panel cannot explain this strand" for a strand the
    // panel explains perfectly well once the parent's claim is dropped.
    if (have_alpha) {
        bool any = false;
        for (size_t k = 0; k < m && !any; ++k) {
            any = delta[k] > NINF;
        }
        if (!any) {
            for (size_t k = 0; k < m; ++k) {
                if (emissions[0][k] > 0.0) {
                    delta[k] = -log((double)m) + log(emissions[0][k]);
                }
            }
        }
    }
    vector<vector<uint16_t>> back(n);
    for (size_t t = 1; t < n; ++t) {
        // The distance from the previous site; see `site_gap`.
        double rho = switch_probability(site_gap(sites[from + t - 1], sites[from + t]).first);
        rho = min(max(rho, 1e-12), 1.0 - 1e-12);
        double S = log(1.0 - rho + rho / (double)m);
        double J = log(rho / (double)m);

        // One strand, so the reduction is the textbook one: stay on the same haplotype, or jump
        // from whichever was best. No leave-one-out maxima are needed -- that complication in the
        // diploid step comes entirely from the two strands sharing one delta.
        size_t arg_best = 0;
        double best_any = NINF;
        for (size_t k = 0; k < m; ++k) {
            if (delta[k] > best_any) {
                best_any = delta[k];
                arg_best = k;
            }
        }
        vector<double> next(m, NINF);
        back[t].assign(m, 0);
        for (size_t k = 0; k < m; ++k) {
            if (!(emissions[t][k] > 0.0)) {
                continue;
            }
            double stay = delta[k] > NINF ? delta[k] + S : NINF;
            double jump = best_any > NINF ? best_any + J : NINF;
            if (stay >= jump && stay > NINF) {
                next[k] = stay + log(emissions[t][k]);
                back[t][k] = (uint16_t)k;
            } else if (jump > NINF) {
                next[k] = jump + log(emissions[t][k]);
                back[t][k] = (uint16_t)arg_best;
            }
        }
        bool any = false;
        for (double v : next) {
            if (v > NINF) { any = true; break; }
        }
        if (!any) {
            // Constraint and panel disagree beyond what the wildcard absorbs: restart here rather
            // than abandoning the window, as the diploid pass does.
            for (size_t k = 0; k < m; ++k) {
                next[k] = emissions[t][k] > 0.0 ? log(emissions[t][k]) : NINF;
                back[t][k] = (uint16_t)k;
            }
        }
        delta = std::move(next);
    }

    size_t best = m;
    double best_v = NINF;
    for (size_t k = 0; k < m; ++k) {
        if (delta[k] > best_v) {
            best_v = delta[k];
            best = k;
        }
    }
    if (best == m) {
        return;
    }
    size_t cur = best;
    for (size_t t = n; t-- > 0;) {
        out[t] = (cur == n_hap) ? WILDCARD : cur;
        if (t > 0) {
            cur = back[t][cur];
        }
    }
}

//------------------------------------------------------------------------------
// LinkageCollector

vector<int> LinkageCollector::compact_allele_space(
        const map<vector<int>, double>& genotype_ln_likelihood,
        const vector<int>& haplotype_traversal,
        int called_trav_i, int called_trav_j) {
    // Only the traversals that can matter: build_emission reads the likelihoods only at pairs of
    // panel-carried alleles, and the constraint needs only the called pair. Leaving out the other
    // candidates keeps the space and its likelihood vector small.
    set<int> needed;
    if (called_trav_i >= 0) {
        needed.insert(called_trav_i);
    }
    if (called_trav_j >= 0) {
        needed.insert(called_trav_j);
    }
    for (int t : haplotype_traversal) {
        if (t >= 0) {
            needed.insert(t);
        }
    }
    // Sorted by candidate index, so the compact numbering is a function of the site alone and not of
    // the order haplotypes happen to be visited in.
    return vector<int>(needed.begin(), needed.end());
}

LinkageCollector::CompactSite LinkageCollector::compact_site(
        const map<vector<int>, double>& genotype_ln_likelihood,
        const vector<int>& haplotype_traversal,
        int called_trav_i, int called_trav_j, size_t ploidy) const {
    CompactSite cs;
    cs.space = compact_allele_space(genotype_ln_likelihood, haplotype_traversal,
                                    called_trav_i, called_trav_j);
    if (cs.space.empty() || cs.space.size() > 127) {
        return cs;
    }
    cs.ci = cs.compact_of(called_trav_i);
    cs.cj = called_trav_j >= 0 ? cs.compact_of(called_trav_j) : cs.ci;
    if (cs.ci < 0 || cs.cj < 0) {
        return cs;   // cannot happen: the called pair is in the space by construction
    }
    const size_t k = cs.space.size();
    cs.site_ploidy = (ploidy == 1 ? 1 : 2);
    const size_t n_gt = cs.site_ploidy == 1 ? k : k * (k + 1) / 2;
    cs.gls.assign(n_gt, -numeric_limits<float>::infinity());
    for (const auto& kv : genotype_ln_likelihood) {
        if (kv.first.size() != cs.site_ploidy) {
            continue;
        }
        if (cs.site_ploidy == 1) {
            int a = cs.compact_of(kv.first[0]);
            if (a >= 0) {
                cs.gls[(size_t)a] = (float)kv.second;
            }
            continue;
        }
        int a = cs.compact_of(kv.first[0]);
        int b = cs.compact_of(kv.first[1]);
        if (a >= 0 && b >= 0) {
            cs.gls[LinkageModel::genotype_index((size_t)a, (size_t)b)] = (float)kv.second;
        }
    }
    cs.ok = true;
    return cs;
}

bool LinkageModel::run_length_site(const vector<string>& alleles, size_t min_run, size_t ref) {
    for (size_t x = 0; x < alleles.size(); ++x) {
        for (size_t y = 0; y < alleles.size(); ++y) {
            if (ref < alleles.size() && x != ref && y != ref) {
                continue;
            }
            // `a` the longer of the pair; one ordering suffices, the other is its mirror.
            const string& a = alleles[x];
            const string& b = alleles[y];
            if (a.size() <= b.size() || a.size() - b.size() > 49) {
                continue;
            }
            size_t pre = 0;
            while (pre < b.size() && a[pre] == b[pre]) {
                ++pre;
            }
            size_t suf = 0;
            while (suf < b.size() - pre && a[a.size() - 1 - suf] == b[b.size() - 1 - suf]) {
                ++suf;
            }
            // What `a` has and `b` lacks, once the shared flanks are gone.
            const size_t lo = pre;
            const size_t hi = a.size() - suf;
            if (pre + suf != b.size()) {
                continue;  // they differ by more than an insertion
            }
            const char base = a[lo];
            bool one_base = true;
            for (size_t i = lo; i < hi && one_base; ++i) {
                one_base = a[i] == base;
            }
            if (!one_base) {
                continue;
            }
            // The run the extra bases sit in, measured in the longer allele.
            size_t start = lo;
            while (start > 0 && a[start - 1] == base) {
                --start;
            }
            size_t end = hi;
            while (end < a.size() && a[end] == base) {
                ++end;
            }
            // Reaching an end of the allele means the run continues into a boundary node made of
            // this base, where the graph cut it: its length is unknown, so it counts as long.
            if (end - start >= min_run || start == 0 || end == a.size()) {
                return true;
            }
        }
    }
    return false;
}

void LinkageCollector::record(const string& contig, size_t position,
                              const map<vector<int>, double>& genotype_ln_likelihood,
                              const vector<int>& haplotype_traversal,
                              int called_trav_i, int called_trav_j,
                              const vector<int>& traversal_to_allele,
                              size_t record_key,
                              const DirectQuality& direct, size_t ploidy,
                              int64_t start_node, int64_t end_node,
                              const SiteContext& ctx) {
    if (genotype_ln_likelihood.empty() || called_trav_i < 0) {
        return;
    }
    // int8 in the panel arena caps a site at 127 alleles; one that exceeded it would lose linkage
    // rather than be mis-linked, which is the safe direction. `ok` covers that and the two other
    // ways a site cannot be described.
    const CompactSite cs = compact_site(genotype_ln_likelihood, haplotype_traversal,
                                        called_trav_i, called_trav_j, ploidy);
    if (!cs.ok) {
        return;
    }
    const vector<int>& space = cs.space;
    const vector<float>& gls = cs.gls;
    const size_t k = space.size();
    const int ci = cs.ci, cj = cs.cj;
    const size_t site_ploidy = cs.site_ploidy;
    auto compact_of = [&](int trav) { return cs.compact_of(trav); };

    lock_guard<std::mutex> guard(mutex);

    uint32_t contig_id;
    {
        auto it = contig_index.find(contig);
        if (it != contig_index.end()) {
            contig_id = it->second;
        } else {
            contig_id = (uint32_t)contig_names.size();
            contig_names.push_back(contig);
            contig_index.emplace(contig, contig_id);
        }
    }

    Entry e;
    e.position = (uint32_t)position;
    e.contig = contig_id;
    e.num_alleles = (uint16_t)k;
    e.called_i = (uint16_t)ci;
    e.called_j = (uint16_t)cj;
    e.record_key = record_key;
    e.explained_share = (float)direct.explained_share;
    e.gq_factor = (float)direct.gq_factor;
    e.achievable_gap = (float)direct.achievable_gap;
    e.start_node = start_node;
    e.end_node = end_node;
    e.ploidy = (uint8_t)site_ploidy;
    e.nested = ctx.nested;
    e.emitted = ctx.emitted;
    e.parent_record_key = ctx.parent_record_key;
    e.parent_crossing = ctx.parent_crossing;
    e.unpositioned = ctx.unpositioned;
    e.chain_key = ctx.chain_key;
    e.generation = (uint8_t)(ctx.generation > 255 ? 255 : ctx.generation);
    e.freq_prior = (float)ctx.freq_prior;

    e.gl_offset = (uint32_t)gl_arena.size();
    for (float v : gls) {
        // float, not double: these are log-likelihood differences fed to an exp(), and the
        // ratios that survive are nowhere near float's precision limit. It halves the arena.
        gl_arena.push_back(v);
    }
    e.hap_offset = (uint32_t)hap_arena.size();
    for (size_t h = 0; h < n_haplotypes; ++h) {
        int trav = h < haplotype_traversal.size() ? haplotype_traversal[h] : -1;
        int a = trav >= 0 ? compact_of(trav) : -1;
        hap_arena.push_back(a >= 0 ? (int8_t)a : (int8_t)-1);
    }
    e.trav_offset = (uint32_t)trav_arena.size();
    e.allele_offset = (uint32_t)allele_arena.size();
    for (size_t i = 0; i < k; ++i) {
        trav_arena.push_back((uint16_t)space[i]);
        int allele = (size_t)space[i] < traversal_to_allele.size() ? traversal_to_allele[space[i]]
                                                                   : -1;
        allele_arena.push_back(allele >= 0 && allele < 127 ? (int8_t)allele : (int8_t)-1);
    }
    const uint32_t at = (uint32_t)entries.size();
    auto last = last_by_key.find(record_key);
    if (last == last_by_key.end()) {
        first_by_key[record_key] = at;
    } else {
        // Counted before the append, while the chain still describes the existing entries. See
        // `num_duplicate_live_keys`.
        if (live_index(record_key) != NO_ENTRY) {
            ++duplicate_live_keys;
        }
        entries[last->second].next_same_key = at;
    }
    last_by_key[record_key] = at;
    if (e.generation >= by_generation.size()) {
        by_generation.resize((size_t)e.generation + 1);
    }
    by_generation[e.generation].push_back(at);
    entries.push_back(e);
}

uint32_t LinkageCollector::live_index(size_t record_key) const {
    auto it = first_by_key.find(record_key);
    for (uint32_t i = it == first_by_key.end() ? NO_ENTRY : it->second;
         i != NO_ENTRY; i = entries[i].next_same_key) {
        if (!entries[i].retracted) {
            return i;
        }
    }
    return NO_ENTRY;
}

size_t LinkageCollector::num_sites_at(size_t generation) const {
    size_t n = 0;
    for (const Entry& e : entries) {
        n += (e.generation == generation && !e.retracted);
    }
    return n;
}

size_t LinkageCollector::max_generation() const {
    size_t g = 0;
    for (const Entry& e : entries) {
        g = max(g, (size_t)e.generation);
    }
    return g;
}

size_t LinkageCollector::bytes() const {
    // Every arena, and the key index, which costs two hash nodes per site.
    return entries.size() * sizeof(Entry)
           + gl_arena.size() * sizeof(float)
           + hap_arena.size() * sizeof(int8_t)
           + trav_arena.size() * sizeof(uint16_t)
           + allele_arena.size() * sizeof(int8_t)
           + (first_by_key.size() + last_by_key.size())
                 * (sizeof(std::pair<const size_t, uint32_t>) + sizeof(void*));
}


bool LinkageCollector::settled_traversals(size_t record_key, int* first, int* second,
                                         size_t* ploidy) const {
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found != NO_ENTRY) {
        const Entry& e = entries[found];
        const int a = traversal_of(trav_arena, e.trav_offset, e.num_alleles, e.final_i);
        const int b = traversal_of(trav_arena, e.trav_offset, e.num_alleles, e.final_j);
        if (a < 0 || b < 0) {
            return false;
        }
        *first = a;
        *second = b;
        *ploidy = e.ploidy;
        return true;
    }
    return false;
}

std::unordered_set<size_t> LinkageCollector::emitted_records() const {
    lock_guard<std::mutex> guard(mutex);
    std::unordered_set<size_t> emitted;
    for (const Entry& e : entries) {
        if (!e.retracted && e.emitted) {
            emitted.insert(e.record_key);
        }
    }
    return emitted;
}

bool LinkageCollector::set_allele_map(size_t record_key,
                                     const vector<int>& traversal_to_allele, bool emitted) {
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found != NO_ENTRY) {
        Entry& e = entries[found];
        // Only ever set, never cleared: under -A a record can be rendered by more than one path, and
        // an assignment would make the result depend on which thread took the mutex last.
        if (emitted && !e.emitted) {
            e.emitted = true;
        }
        // The span already exists, filled with -1 by `record()`. Rewrite it in place rather than
        // appending: the compact space has not changed, only what each of its alleles is called in
        // the VCF, so there is nothing to re-point and no arena growth.
        for (size_t c = 0; c < e.num_alleles; ++c) {
            const size_t at_trav = e.trav_offset + c;
            const size_t at_allele = e.allele_offset + c;
            if (at_trav >= trav_arena.size() || at_allele >= allele_arena.size()) {
                break;
            }
            const size_t trav = trav_arena[at_trav];
            const int allele = trav < traversal_to_allele.size() ? traversal_to_allele[trav] : -1;
            allele_arena[at_allele] = allele >= 0 && allele < 127 ? (int8_t)allele : (int8_t)-1;
        }
        return true;
    }
    return false;
}



bool LinkageCollector::has_entry(size_t record_key) const {
    lock_guard<std::mutex> guard(mutex);
    return live_index(record_key) != NO_ENTRY;
}

bool LinkageCollector::set_position(size_t record_key, size_t position) {
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found == NO_ENTRY) {
        return false;
    }
    entries[found].position = position;
    return true;
}

bool LinkageCollector::rescore(size_t record_key,
                               const map<vector<int>, double>& genotype_ln_likelihood,
                               const vector<int>& haplotype_traversal,
                               int called_trav_i, int called_trav_j) {
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found == NO_ENTRY) {
        return false;
    }
    Entry& e = entries[found];
    const CompactSite cs = compact_site(genotype_ln_likelihood, haplotype_traversal,
                                        called_trav_i, called_trav_j, e.ploidy);
    if (!cs.ok) {
        return false;
    }
    // Whether the compact allele space changed size. It can: the space is the panel-carried
    // traversals plus the called pair, so moving the call onto or off an allele that no panel
    // haplotype carries adds or removes one.
    bool same_space = cs.space.size() == e.num_alleles;
    for (size_t i = 0; same_space && i < cs.space.size(); ++i) {
        same_space = trav_arena[e.trav_offset + i] == (uint16_t)cs.space[i];
    }
    const size_t n_gt = cs.gls.size();
    if (same_space) {
        // In place. `gl_arena` stores no length -- a reader recomputes it from `num_alleles` --
        // so this is only safe because the shape is provably unchanged.
        const size_t expect = e.ploidy == 1 ? (size_t)e.num_alleles
                                            : (size_t)e.num_alleles * ((size_t)e.num_alleles + 1) / 2;
        if (n_gt != expect) {
            return false;
        }
        for (size_t i = 0; i < n_gt; ++i) {
            gl_arena[e.gl_offset + i] = cs.gls[i];
        }
    } else {
        // Appended, and the offsets moved to the new slice, since a slice of another width does
        // not fit in place; the old slice is left unused, and every other entry's offsets stay
        // valid. The fields that place the site in the snarl tree are not changed.
        e.gl_offset = (uint32_t)gl_arena.size();
        for (float v : cs.gls) {
            gl_arena.push_back(v);
        }
        e.hap_offset = (uint32_t)hap_arena.size();
        for (size_t h = 0; h < n_haplotypes; ++h) {
            const int trav = h < haplotype_traversal.size() ? haplotype_traversal[h] : -1;
            const int a = trav >= 0 ? cs.compact_of(trav) : -1;
            hap_arena.push_back(a >= 0 ? (int8_t)a : (int8_t)-1);
        }
        e.trav_offset = (uint32_t)trav_arena.size();
        e.allele_offset = (uint32_t)allele_arena.size();
        for (size_t i = 0; i < cs.space.size(); ++i) {
            trav_arena.push_back((uint16_t)cs.space[i]);
            // No allele map yet; `set_allele_map` supplies it when the record is written.
            allele_arena.push_back((int8_t)-1);
        }
        e.num_alleles = (uint16_t)cs.space.size();
    }
    // The per-site call moves with the likelihoods. `final_*` is left alone: it is the decode's
    // output and `resolve_generation` resets it from `called_*` on its next pass, which is the
    // pass this rescore exists to feed.
    e.called_i = (uint16_t)cs.ci;
    e.called_j = (uint16_t)cs.cj;
    return true;
}

bool LinkageCollector::retract(size_t record_key) {
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found == NO_ENTRY) {
        return false;
    }
    entries[found].retracted = true;
    // A retracted site has no settled genotype, so it is no longer moved.
    if (live_index(record_key) == NO_ENTRY) {
        moved_quality_by_record.erase(record_key);
    }
    return true;
}

/// Translate the settled compact pair into the traversal on each strand, against which crossing
/// masks are tested, and the VCF allele on each strand, which is what the record writes.
void LinkageCollector::finish_phase_call(PhaseCall& pc, const Entry& e) const {
    const size_t c_first = pc.allele_first, c_second = pc.allele_second;
    pc.trav_first = traversal_of(trav_arena, e.trav_offset, e.num_alleles, c_first);
    pc.trav_second = traversal_of(trav_arena, e.trav_offset, e.num_alleles, c_second);
    int v_first = -1, v_second = -1;
    bool fell_back = false;
    render_phase_pair(allele_arena, e.allele_offset, e.num_alleles, c_first, c_second,
                      e.called_i, e.called_j, &v_first, &v_second, &fell_back);
    pc.allele_first = v_first >= 0 ? (size_t)v_first : LinkageModel::WILDCARD;
    pc.allele_second = v_second >= 0 ? (size_t)v_second : LinkageModel::WILDCARD;
    if (fell_back) {
        pc.order_arbitrary = pc.order_arbitrary || (v_first != v_second);
    }
}


/// Which of its parent's two strands a nested chain at ploidy 1 sits on, or -1.
///
/// `carrying` is the parent candidate traversal that crosses the chain, from `relate_to_parent`.
/// A nested parent at ploidy 1 is on one strand, so everything inside it is on that strand: its
/// own `nested_strand` is the answer. A diploid parent names the strand by which of its two
/// settled traversals carries the chain. A haploid top-level parent has only one strand, so its
/// children get -1 and are not written as `a|.`.
static inline int nested_strand_of(int carrying, size_t parent_ploidy, int parent_trav_first,
                                   int parent_trav_second, int parent_nested_strand) {
    if (carrying < 0) {
        return -1;
    }
    if (parent_ploidy == 1 && parent_nested_strand >= 0) {
        return parent_trav_first == carrying ? parent_nested_strand : -1;
    }
    if (parent_ploidy == 2) {
        if (parent_trav_first == carrying) {
            return 0;
        }
        if (parent_trav_second == carrying) {
            return 1;
        }
    }
    return -1;
}

size_t LinkageCollector::resolve_generation(
        size_t generation, bool last, vector<PhaseCall>* phasing_out) {
    size_t moved = 0;
    if (!model.active() || entries.empty()) {
        return moved;
    }

    // For each site of an earlier generation, by record key: the phase it settled on, so that a
    // clamped site can be pinned to it, and what `nested_strand_of` needs to place a child.
    struct PinnedPhase {
        size_t first;
        size_t second;
        int trav_first;
        int trav_second;
        size_t ploidy;
        int nested_strand;
        bool order_arbitrary;
    };
    // Rebuilt from all of `phasing_out` on each call. It is not read incrementally, because the
    // last call sorts `phasing_out` in place.
    unordered_map<size_t, PinnedPhase> pinned_phase;
    if (phasing_out != nullptr && generation > 0) {
        pinned_phase.reserve(phasing_out->size() * 2);
        for (const PhaseCall& pc : *phasing_out) {
            pinned_phase[pc.record_key] = PinnedPhase{pc.hap_first, pc.hap_second,
                                                      pc.trav_first, pc.trav_second,
                                                      pc.ploidy, (int)pc.nested_strand,
                                                      pc.order_arbitrary};
        }
    }
    // Where a nested haploid chain sits, by record key, so the PhaseCall the chain loop emits can
    // name it. Derived once where the parent's settled pair is in hand, which is the only place
    // both facts are available.
    struct NestedPlacement {
        int strand = -1;
        bool order_arbitrary = false;   // the parent's: this strand IS that coin flip, one level down
        // Whether the site may name a panel haplotype. Under a haploid parent there is one strand
        // and the child sits on it. Under a diploid parent, `nested_strand_of` returns -1 when
        // both settled traversals carry the chain or when the settled pair could not be read, and
        // then no haplotype is named. (A parent carrying no copy cannot occur here, since the
        // barrier retracts such a subtree first.) The two cases are counted separately.
        bool nameable = true;
    };
    unordered_map<size_t, NestedPlacement> unified_strand;

    // This generation's entries, in append order.
    static const vector<uint32_t> no_entries;
    const vector<uint32_t>& this_generation =
        generation < by_generation.size() ? by_generation[generation] : no_entries;

    // Default every site this pass considers to its own per-site call, so that whatever a chain or a
    // sweep fails to reach still has a coherent genotype for a later generation to clamp. Overwritten
    // below wherever something is actually settled.
    for (uint32_t idx : this_generation) {
        Entry& e = entries[idx];
        e.final_i = e.called_i;
        e.final_j = e.ploidy == 1 ? e.called_i : e.called_j;
    }

    // Group by contig, then sort by reference position. Node-ID order is close to reference order
    // but not the same, and the transition probabilities depend on the distances. Only at the top
    // level: below it every live site is nested and is grouped with its parent instead.
    vector<vector<size_t>> by_contig;
    if (generation == 0) {
        by_contig.resize(contig_names.size());
        for (uint32_t i : this_generation) {
            if (entries[i].retracted) {
                continue;   // the settled parent does not carry the chain, so there is no site here
            }
            by_contig[entries[i].contig].push_back(i);
        }
    }

    // A top-level linkage chain is a maximal run of one ploidy on one contig, since strands do not
    // correspond across a ploidy change such as a --ploidy-bed boundary.
    //
    // For each chain, the message to condition it on (empty for none) and its phase set (SIZE_MAX
    // to take it from the chain's first site). `deltas` owns the messages that `chain_context`
    // points to, so it must outlive the decode below; it is a deque so that adding a message does
    // not move the others.
    std::deque<vector<double>> deltas;
    vector<const vector<double>*> chain_context;
    vector<size_t> chain_phase_set;

    vector<vector<size_t>> chains;
    for (auto& contig_indices : by_contig) {
        if (contig_indices.empty()) {
            continue;
        }
        // Position, then the site's own key. Sites arrive in whatever order the threads finished,
        // and two records can share a position, so sorting on position alone leaves their relative
        // order down to scheduling -- which would make the output depend on --threads. The key is
        // derived from the snarl ID, so it is a property of the site rather than of the run.
        sort(contig_indices.begin(), contig_indices.end(), [&](size_t a, size_t b) {
            if (entries[a].position != entries[b].position) {
                return entries[a].position < entries[b].position;
            }
            return entries[a].record_key < entries[b].record_key;
        });
        // Sorted first, so that a run is contiguous in reference order. Nested sites are kept out of
        // these runs and decoded afterwards with their parents, since a nested ploidy-1 site
        // between diploid neighbours would otherwise split the run. A change of the contig's own
        // ploidy still splits it. Every entry here is top level, since `by_contig` is built only
        // at generation 0.
        size_t run_start = 0;
        while (run_start < contig_indices.size()) {
            size_t run_end = run_start + 1;
            while (run_end < contig_indices.size()
                   && entries[contig_indices[run_end]].ploidy
                          == entries[contig_indices[run_start]].ploidy) {
                ++run_end;
            }
            chains.emplace_back(contig_indices.begin() + run_start,
                                contig_indices.begin() + run_end);
            chain_context.push_back(nullptr);
            chain_phase_set.push_back(numeric_limits<size_t>::max());
            run_start = run_end;
        }
    }

    // Sites and groups decoded per parent, for reporting.
    size_t grouped_sites = 0, grouped_groups = 0;
    // After generation 0, each nested chain is decoded on its own, conditioned on its parent's
    // settled state, so the cost grows with the number of children rather than with the length of
    // the contig. Sites within a chain are linked; two chains under the same parent have no
    // transitions between them.
    if (generation > 0 && !pinned_phase.empty()) {
        unordered_map<size_t, size_t> index_of_key;
        index_of_key.reserve(entries.size() * 2);
        for (size_t i = 0; i < entries.size(); ++i) {
            if (!entries[i].retracted) {
                index_of_key[entries[i].record_key] = i;
            }
        }
        // Sites grouped by (parent, chain, ploidy, carrying traversal); see `group_key`.
        map<tuple<size_t, size_t, size_t, int>, vector<size_t>> by_parent;
        // What the parent's settled pair implies about a child, from `relate_to_parent`.
        auto relate = [&](const Entry& child, const Entry& parent) {
            const int ta = traversal_of(trav_arena, parent.trav_offset, parent.num_alleles,
                                        parent.final_i);
            const int tb = parent.ploidy == 2
                               ? traversal_of(trav_arena, parent.trav_offset, parent.num_alleles,
                                              parent.final_j)
                               : -1;
            return LinkageCollector::relate_to_parent(child.parent_crossing, ta, tb);
        };

        // (parent, chain), where the chain half is its boundary pair from the graph. A snarl the
        // decomposition puts in no chain becomes its own group rather than being pooled with every
        // other such snarl under this parent.
        auto group_key = [&](const Entry& e, const Entry& parent) {
            // The chain's boundary pair, from the graph.
            const size_t chain = e.chain_key != 0
                                     ? e.chain_key
                                     : numeric_limits<size_t>::max() - e.record_key;
            // Ploidy is part of the key. A group is decoded at one ploidy, taken from its first
            // member, and a site decoded at the wrong ploidy would index past its likelihood
            // vector. Children in one chain of one parent can differ in ploidy, since a child's
            // ploidy is the number of the parent's settled alleles that cross it.
            //
            // So is, at ploidy 1, the parent's settled traversal that carries the site, which
            // names the strand the site is on. A group is placed on one strand, and two sites of
            // one chain can be carried by different strands of the parent.
            const int carrying = e.ploidy == 1 ? relate(e, parent).carrying_trav : -1;
            return make_tuple(e.parent_record_key, chain, e.ploidy, carrying);
        };
        // Group every live site of this generation with its parent. A site that cannot be grouped
        // is decoded alone. Groups are sorted afterwards on (position, record key), a total order,
        // so membership does not depend on the order in which sites arrived.
        vector<vector<size_t>> ungrouped;
        for (uint32_t idx : this_generation) {
            const Entry& e = entries[idx];
            if (e.retracted) {
                continue;
            }
            auto par = index_of_key.find(e.parent_record_key);
            if (e.parent_record_key == 0) {
                ++model.counters.grp_no_parent;
            } else if (par == index_of_key.end()) {
                ++model.counters.grp_no_entry;
            } else {
                by_parent[group_key(e, entries[par->second])].push_back(idx);
                continue;
            }
            // Decoded alone rather than dropped, so that it is still settled and phased.
            ++model.counters.grp_vetoed;
            ungrouped.push_back(vector<size_t>{idx});
        }
        if (!by_parent.empty()) {
            // The phase set a group belongs to is its parent's, never the group's own first site:
            // a group is a unit of decoding, a phase block is a unit of meaning, and taking one
            // from the other is what fragments the output.
            unordered_map<size_t, size_t> phase_set_of;
            for (const PhaseCall& pc : (phasing_out != nullptr ? *phasing_out
                                                              : vector<PhaseCall>())) {
                phase_set_of[pc.record_key] = pc.phase_set;
            }
            vector<vector<size_t>> groups;
            vector<const vector<double>*> gctx;
            vector<size_t> gps;
            for (auto& kv : by_parent) {
                const size_t pidx = index_of_key[std::get<0>(kv.first)];
                // The parent joins its group when their ploidies match, as for a haploid child of
                // a haploid parent. A haploid child of a diploid parent is decoded without it, and
                // the parent reaches the group through the entering message instead. Every member
                // of the group has this ploidy, since it is part of the key.
                const size_t group_ploidy =
                    kv.second.empty() ? entries[pidx].ploidy : std::get<2>(kv.first);
                const bool parent_in_group = entries[pidx].ploidy == group_ploidy;
                vector<size_t> group;
                if (parent_in_group) {
                    group.push_back(pidx);
                }
                // A total order on (position, record key). A comparator that skips a key when
                // either side lacks it is not a strict weak ordering, and std::sort with one is
                // undefined.
                sort(kv.second.begin(), kv.second.end(), [&](size_t a, size_t c) {
                    const Entry& ea = entries[a];
                    const Entry& ec = entries[c];
                    if (ea.position != ec.position) {
                        return ea.position < ec.position;
                    }
                    return ea.record_key < ec.record_key;
                });
                if (pinned_phase.count(entries[pidx].record_key) != 0) {
                    ++model.counters.group_parent_pinned;
                } else {
                    // No PhaseCall for the parent, so nothing ties this group's orientation to it.
                    ++model.counters.group_parent_unpinned;
                }
                group.insert(group.end(), kv.second.begin(), kv.second.end());
                grouped_sites += group.size();
                ++grouped_groups;
                groups.push_back(std::move(group));
                // The entering message: a point mass at the parent's settled state, owned by
                // `deltas`. The parent's WILDCARD, `(size_t)-1` outside the model, is state
                // `n_haplotypes` inside it, so it is translated before indexing.
                const size_t m = n_haplotypes + 1;
                auto state_of = [&](size_t h) {
                    return h == LinkageModel::WILDCARD ? n_haplotypes : h;
                };
                auto pin = pinned_phase.find(std::get<0>(kv.first));
                if (pin == pinned_phase.end()) {
                    gctx.push_back(nullptr);
                } else if (group_ploidy == 1) {
                    // One haplotype, not a pair: the group sits on one of the parent's strands, the
                    // one `nested_strand_of` names. Every member is carried by the same parent
                    // traversal, which is part of the group key, so the first member stands for
                    // all. A haploid parent records its haplotype in `hap_first`, and
                    // `hap_second` is meaningless at ploidy 1.
                    const Entry& child = entries[kv.second.front()];
                    const int carrying = relate(child, entries[pidx]).carrying_trav;
                    const int strand = nested_strand_of(carrying, pin->second.ploidy,
                                                        pin->second.trav_first,
                                                        pin->second.trav_second,
                                                        pin->second.nested_strand);
                    const size_t hap = pin->second.ploidy == 1
                                           ? pin->second.first
                                           : (strand == 1 ? pin->second.second : pin->second.first);
                    // Give no message where the parent's haplotype does not pass through the child,
                    // since it names no allele there. A haploid parent has one strand, so it needs
                    // no `nested_strand`.
                    const size_t st = state_of(hap);
                    const bool have_hap = pin->second.ploidy == 1 || strand >= 0;
                    bool traversed = have_hap && hap != LinkageModel::WILDCARD
                                     && st < n_haplotypes
                                     && (int)hap_arena[child.hap_offset + hap] >= 0;
                    if (traversed) {
                        deltas.emplace_back(m, 0.0);
                        deltas.back()[st] = 1.0;
                        gctx.push_back(&deltas.back());
                    } else {
                        gctx.push_back(nullptr);
                    }
                    for (size_t idx : kv.second) {
                        unified_strand[entries[idx].record_key] = NestedPlacement{
                            strand, strand >= 0 ? pin->second.order_arbitrary : false, have_hap};
                    }
                    if (strand >= 0) {
                        model.counters.nest_strand += kv.second.size();
                    } else if (have_hap) {
                        model.counters.nest_one_hap += kv.second.size();
                    } else if (carrying == -2) {
                        model.counters.nest_both += kv.second.size();
                    } else {
                        model.counters.nest_unreadable += kv.second.size();
                    }
                } else if (state_of(pin->second.first) >= m || state_of(pin->second.second) >= m) {
                    gctx.push_back(nullptr);
                } else {
                    deltas.emplace_back(m * m, 0.0);
                    deltas.back()[state_of(pin->second.first) * m
                                  + state_of(pin->second.second)] = 1.0;
                    gctx.push_back(&deltas.back());
                }
                auto ps = phase_set_of.find(std::get<0>(kv.first));
                gps.push_back(ps != phase_set_of.end() ? ps->second
                                                       : numeric_limits<size_t>::max());
            }
            // Chains that could not be grouped are kept, so that their sites are still decoded.
            for (vector<size_t>& kc : ungrouped) {
                groups.push_back(std::move(kc));
                gctx.push_back(nullptr);
                gps.push_back(numeric_limits<size_t>::max());
            }
            chains.swap(groups);
            chain_context.swap(gctx);
            chain_phase_set.swap(gps);
        }
    }

    for (size_t chain_i = 0; chain_i < chains.size(); ++chain_i) {
        vector<size_t>& indices = chains[chain_i];
        if (indices.empty()) {
            continue;
        }
        // Count the chain's sites of this generation, and its pinned sites of earlier ones. A chain
        // with none of this generation's sites has nothing to decide, since all its sites are
        // clamped, and is skipped.
        size_t live_here = 0, pinned_here = 0;
        for (size_t idx : indices) {
            if (entries[idx].generation == generation) {
                ++live_here;
            } else if (pinned_phase.count(entries[idx].record_key) != 0) {
                ++pinned_here;
            }
        }
        if (live_here == 0) {
            continue;
        }
        if (generation > 0) {
#pragma omp critical (cerr)
            std::cerr << "[vg call] linkage generation " << generation << ": chain decodes "
                      << indices.size() << " sites for " << live_here << " of its own, "
                      << pinned_here << " pinned" << std::endl;
        }
        // A one-site chain has nothing to link to, so the model cannot change its genotype, but it
        // is still phased, so that it appears in `phasing_out` and the mosaic.

        // Every site of a chain has the same ploidy, and this checks it, since a site decoded at
        // the wrong ploidy would index past its likelihood vector.
        size_t chain_ploidy = entries[indices.front()].ploidy;
        for (size_t idx : indices) {
            assert(entries[idx].ploidy == chain_ploidy
                   && "a decode chain must be homogeneous in ploidy");
        }

        vector<LinkageModel::Site> sites;
        sites.reserve(indices.size());
        for (size_t k = 0; k < indices.size(); ++k) {
            const size_t idx = indices[k];
            const Entry& e = entries[idx];
            LinkageModel::Site s;
            s.position = e.position;
            s.unpositioned = e.unpositioned;
            s.num_alleles = e.num_alleles;
            s.ploidy = e.ploidy;
            s.freq_prior = e.freq_prior;
            size_t n_gt = e.ploidy == 1
                              ? (size_t)e.num_alleles
                              : (size_t)e.num_alleles * ((size_t)e.num_alleles + 1) / 2;
            s.genotype_ln_likelihood.reserve(n_gt);
            for (size_t g = 0; g < n_gt; ++g) {
                s.genotype_ln_likelihood.push_back((double)gl_arena[e.gl_offset + g]);
            }
            s.haplotype_allele.reserve(n_haplotypes);
            for (size_t h = 0; h < n_haplotypes; ++h) {
                s.haplotype_allele.push_back((int)hap_arena[e.hap_offset + h]);
            }
            if (e.generation < generation) {
                // Clamped. A delta emission at the settled genotype, so the site still carries
                // transition context for its neighbours -- which is why it is in the chain at all --
                // while being unable to move. build_emission maps a non-finite entry to zero mass,
                // so this needs nothing from the model.
                size_t settled = e.ploidy == 1
                                     ? (size_t)e.final_i
                                     : LinkageModel::genotype_index(e.final_i, e.final_j);
                for (size_t g = 0; g < s.genotype_ln_likelihood.size(); ++g) {
                    s.genotype_ln_likelihood[g] = (g == settled)
                                                      ? 0.0
                                                      : -numeric_limits<double>::infinity();
                }
                auto pin = pinned_phase.find(e.record_key);
                if (pin != pinned_phase.end()) {
                    s.pinned = true;
                    s.pin_first = pin->second.first;
                    s.pin_second = pin->second.second;
                }
            }
            sites.push_back(std::move(s));
        }

        vector<vector<double>> posteriors;
        // The entering message, a point mass built from the parent's PhaseCall: over ordered pairs
        // for a diploid group, over single haplotypes for a haploid one. Each decode ignores a
        // message of the wrong size for its states, so it is passed without regard to ploidy.
        const vector<double>* ctx =
            chain_i < chain_context.size() ? chain_context[chain_i] : nullptr;
        if (chain_ploidy == 1) {
            posteriors = model.posteriors(sites, 1, ctx);
        } else if (ctx != nullptr) {
            // A child group, conditioned on its parent's settled state.
            model.segment_posteriors(sites, 0, sites.size(), ctx, nullptr, posteriors);
        } else {
            posteriors = model.posteriors(sites, 2);
        }

        // The genotype each site ends up with, whether or not linkage moved it. This is what the
        // phasing is constrained to: phasing the pre-linkage calls would describe a genotype set
        // that never reaches the VCF.
        vector<size_t> final_genotype(indices.size(), LinkageModel::NO_CONSTRAINT);

        for (size_t t = 0; t < indices.size(); ++t) {
            const Entry& e = entries[indices[t]];
            const vector<double>& post = posteriors[t];
            size_t before = e.ploidy == 1 ? (size_t)e.called_i
                                          : LinkageModel::genotype_index(e.called_i, e.called_j);
            if (e.generation < generation) {
                // Clamped: it was settled, reported and emitted at its own generation. Its genotype
                // still has to reach `final_genotype`, because that is what the phasing below is
                // constrained to, but it must not produce a second Change.
                final_genotype[t] = e.ploidy == 1
                                        ? (size_t)e.final_i
                                        : LinkageModel::genotype_index(e.final_i, e.final_j);
                continue;
            }
            if (post.empty()) {
                final_genotype[t] = before;
                moved_quality_by_record.erase(e.record_key);
                continue;
            }
            size_t best = 0;
            for (size_t g = 1; g < post.size(); ++g) {
                if (post[g] > post[best]) {
                    best = g;
                }
            }
            final_genotype[t] = best;
            // Decode the genotype index back to its allele pair. At ploidy 1 the index *is* the
            // allele, and both slots carry it so the change applies through the same path.
            size_t i = best, j = best;
            if (e.ploidy != 1) {
                j = 0;
                while (LinkageModel::genotype_index(0, j) <= best) {
                    ++j;
                }
                --j;
                i = best - (j * (j + 1) / 2);
            }
            // Settled. A later generation clamps the site here instead of reconsidering it, which is
            // what makes a parent's genotype final before any of its children is called.
            entries[indices[t]].final_i = (uint16_t)i;
            entries[indices[t]].final_j = (uint16_t)j;
            if (best != before) {
                // Counted. The record is built later from `final_i`/`final_j`. The posterior is
                // kept, since the record's GQ is computed from it.
                MovedQuality& q = moved_quality_by_record[e.record_key];
                q.posterior = post[best];
                q.direct.explained_share = e.explained_share;
                q.direct.gq_factor = e.gq_factor;
                q.direct.achievable_gap = e.achievable_gap;
                ++moved;
            } else {
                // Unmoved in this resolution, whatever an earlier barrier run concluded.
                moved_quality_by_record.erase(e.record_key);
            }
        }

        if (phasing_out == nullptr) {
            continue;
        }
        // At ploidy 1, `phasing` gives one haplotype per site. It gets the same message as the
        // posteriors above.
        vector<LinkageModel::Phase> phase =
            model.phasing(sites, final_genotype, chain_ploidy, ctx);
        // One phase set per chain, named by its first site's position. The windows are pinned to
        // each other, so the path is continuous across the whole chain.
        size_t phase_set = sites.empty() ? 0 : sites.front().position;
        if (chain_i < chain_phase_set.size()
            && chain_phase_set[chain_i] != numeric_limits<size_t>::max()) {
            phase_set = chain_phase_set[chain_i];
        }
        for (size_t t = 0; t < indices.size() && t < phase.size(); ++t) {
            const Entry& e = entries[indices[t]];
            if (e.generation < generation) {
                continue;   // its PhaseCall was emitted, and pinned above, at its own generation
            }
            const LinkageModel::Phase& ph = phase[t];
            // Read the ordered allele pair off the haplotypes the path chose. Where a strand is
            // on the wildcard the panel does not name its allele, so fall back to the genotype's
            // own order -- the phase is then unsupported at that strand rather than wrong.
            size_t want = final_genotype[t];
            size_t i = want, j = want;
            if (e.ploidy != 1) {
                j = 0;
                while (LinkageModel::genotype_index(0, j) <= want) {
                    ++j;
                }
                --j;
                i = want - (j * (j + 1) / 2);
            }

            int a = -1, b = -1;
            if (ph.first != LinkageModel::WILDCARD
                && ph.first < sites[t].haplotype_allele.size()) {
                a = sites[t].haplotype_allele[ph.first];
            }
            if (ph.second != LinkageModel::WILDCARD
                && ph.second < sites[t].haplotype_allele.size()) {
                b = sites[t].haplotype_allele[ph.second];
            }
            PhaseCall pc;
            pc.ploidy = e.ploidy;
            pc.generation = e.generation;
            // The strand of its parent that a nested ploidy-1 chain sits on, found earlier where the
            // parent's settled pair was at hand.
            int nested_slot = -1;
            bool nameable = true;
            {
                auto us = unified_strand.find(e.record_key);
                if (us != unified_strand.end()) {
                    pc.nested_strand = (int8_t)us->second.strand;
                    nested_slot = us->second.strand;
                    nameable = us->second.nameable;
                    // The strand came off the parent's pair. If the panel did not determine that
                    // pair's order, this strand is that coin flip one level down, and reporting it
                    // as determined overstates what the panel said.
                    pc.order_arbitrary = pc.order_arbitrary || us->second.order_arbitrary;
                }
            }
            pc.record_key = e.record_key;
            pc.contig = contig_names[e.contig];
            pc.position = e.position;
            // Fill the slot that `nested_strand` names, since the mosaic reads that slot and treats
            // the other as empty.
            if (!nameable) {
                pc.hap_first = LinkageModel::WILDCARD;
                pc.hap_second = LinkageModel::WILDCARD;
            } else if (nested_slot == 1) {
                pc.hap_first = LinkageModel::WILDCARD;
                pc.hap_second = ph.first;
            } else {
                pc.hap_first = ph.first;
                pc.hap_second = ph.second;
            }
            pc.start_node = e.start_node;
            pc.end_node = e.end_node;
            pc.phase_set = phase_set;
            if (e.ploidy == 1) {
                // One strand: the called allele sits on it, and the haplotype is whatever the
                // path chose. There is no second slot to fill.
                pc.allele_first = want;
                pc.allele_second = want;
            } else if (a >= 0 && b >= 0
                       && LinkageModel::genotype_index((size_t)a, (size_t)b) == want) {
                pc.allele_first = (size_t)a;
                pc.allele_second = (size_t)b;
            } else if (a >= 0 && ((size_t)a == i || (size_t)a == j)) {
                pc.allele_first = (size_t)a;
                pc.allele_second = ((size_t)a == i) ? j : i;
            } else if (b >= 0 && ((size_t)b == i || (size_t)b == j)) {
                pc.allele_second = (size_t)b;
                pc.allele_first = ((size_t)b == i) ? j : i;
            } else {
                // Neither panel haplotype spells either called allele, so nothing orders the pair.
                // Sorted order is written, which pairs allele_first with hap_first only by accident.
                pc.allele_first = i;
                pc.allele_second = j;
                pc.order_arbitrary = (i != j);
            }
            finish_phase_call(pc, e);
            phasing_out->push_back(pc);
        }
    }

    if (generation > 0 && (model.counters.pin_applied.load() + model.counters.pin_declined.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": pins -- "
                  << model.counters.pin_applied.load() << " applied, " << model.counters.pin_declined.load()
                  << " REFUSED (the pinned pair cannot spell the constrained genotype, so that"
                  << " site's orientation is free); groups whose parent was pinnable: "
                  << model.counters.group_parent_pinned.load() << ", not pinnable: "
                  << model.counters.group_parent_unpinned.load() << std::endl;
    }

    // Report only when something was declined.
    if (generation > 0
        && (model.counters.grp_no_parent.load() + model.counters.grp_no_entry.load() + model.counters.grp_vetoed.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": grouping declines so far -- "
                  << model.counters.grp_no_parent.load() << " sites with no parent key, "
                  << model.counters.grp_no_entry.load() << " whose parent has no live entry; "
                  << model.counters.grp_vetoed.load() << " chains kept ungrouped in total" << std::endl;
    }

    if (grouped_groups > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": decoded " << grouped_sites
                  << " sites in " << grouped_groups << " per-parent groups instead of the contig"
                  << " chain" << std::endl;
    }

    if (generation > 0
        && (model.counters.nest_strand.load() + model.counters.nest_one_hap.load() + model.counters.nest_both.load()
            + model.counters.nest_unreadable.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] nested strands: " << model.counters.nest_strand.load()
                  << " on one of a diploid parent's two strands, " << model.counters.nest_one_hap.load()
                  << " on a haploid parent's single haplotype (no strand to choose), "
                  << model.counters.nest_both.load() << " carried on both parent strands, "
                  << model.counters.nest_unreadable.load()
                  << " whose parent's settled pair could not be read -- the last two name no"
                  << " haplotype" << std::endl;
    }

    // Sort `phasing_out` into reference order. Nested sites were appended after all the chains, so
    // it is not yet in order, and the mosaic writer reads it as one ordered pass. The key matches
    // the chain sort: position, then the site's own key, so that the order does not depend on the
    // threads.
    if (phasing_out != nullptr && last) {
        // Stable: a site re-rendered at the barrier has two PhaseCalls with the same record key, and
        // the readers keep the last one for a key, which must be the later generation's.
        // `phasing_out` is appended one generation at a time, so a stable sort keeps it last.
        std::stable_sort(phasing_out->begin(), phasing_out->end(),
                  [](const PhaseCall& a, const PhaseCall& b) {
                      if (a.contig != b.contig) {
                          return a.contig < b.contig;
                      }
                      if (a.position != b.position) {
                          return a.position < b.position;
                      }
                      return a.record_key < b.record_key;
                  });
    }
    return moved;
}

}
