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

/// The allele a panel haplotype carries at a site, or -1 for "nothing to say".
///
/// Out-of-range indices are absent rather than trusted. A caller that passes an allele index from
/// the wrong numbering -- traversal order instead of VCF allele order, say -- would otherwise index
/// past the genotype vector and corrupt the heap, which is exactly what happened the first time
/// this was wired up.
static inline int allele_at(const LinkageModel::Site& site, size_t h, size_t n_hap) {
    if (h >= n_hap || h >= site.haplotype_allele.size()) {
        return -1;
    }
    int allele = site.haplotype_allele[h];
    return (allele >= 0 && (size_t)allele < site.num_alleles) ? allele : -1;
}

/// The VCF allele a compact allele was emitted as, or -1 where none was.
///
/// The collector works in the genotyper's traversal space; the VCF numbering is a rendering of it.
/// -1 is reachable now that the two differ: a traversal the panel carries can reach the model
/// without the record carrying an ALT for it, in which case there is nothing to write and the change
/// is skipped rather than guessed at.
static inline int vcf_allele_of(const vector<int8_t>& allele_arena, size_t allele_offset,
                                size_t num_alleles, size_t compact) {
    if (compact >= num_alleles) {
        return -1;
    }
    size_t at = allele_offset + compact;
    return at < allele_arena.size() ? (int)allele_arena[at] : -1;
}

/// The VCF alleles a phased GT may name for a settled compact pair.
///
/// Phasing re-orders a genotype; it never re-decides one. That is the property that made it safe to
/// turn on by default, and it is checked again in apply_phasing, which declines any phased GT that
/// is not a permutation of the line's own. So when the settled pair cannot be rendered -- the model
/// chose a traversal this record carries no ALT for, and the genotype patch was correctly skipped --
/// the phase must describe the pair the line still carries, which is the called one. Emitting the
/// settled pair anyway got the patch declined and cost the record its phase and its phase set
/// outright: 3 records on the small fixture, 1,627 on chr20.
///
/// The called pair is renderable by construction: it came from the genotype the record was emitted
/// with. Order across a fallback pair is not determined by anything, so callers flag it.
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

/// The candidate traversal a compact allele stands for. This is the genome fact: what a crossing
/// mask is tested against, and what a per-haplotype path through the snarl is made of.
static inline int traversal_of(const vector<uint16_t>& trav_arena, size_t trav_offset,
                               size_t num_alleles, size_t compact) {
    if (compact >= num_alleles) {
        return -1;
    }
    size_t at = trav_offset + compact;
    return at < trav_arena.size() ? (int)trav_arena[at] : -1;
}

/// Relative P(reads | genotype implied by state) for every ordered pair of panel haplotypes,
/// with the wildcard last. Row-normalised only in the sense that the site's best genotype is 1,
/// which keeps the numbers in range without changing any ratio.
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
static void transition_apply(const vector<double>& in, size_t m, double rho,
                             vector<double>& out) {
    double stay = 1.0 - rho;
    double jump = rho / (double)m;
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
            out[a * m + b] = stay * stay * in[a * m + b]
                             + stay * jump * (row[a] + col[b])
                             + jump * jump * total;
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
/// `transition_apply` sums, and there the factorisation collapses the pairwise loop into four
/// terms. Maximising does *not* separate the same way, because delta(a,b) couples the strands:
/// max over (a,b) of delta(a,b) + f(a) + g(b) is not a pair of independent 1-D maxima. But
/// T(x->y) takes only two values, so the reduction is by cases on which strands stayed:
///
///     delta'(a',b') = ln e(a',b') + max of
///         delta(a',b')                        + 2S      both stayed
///         max_{b != b'} delta(a',b)           + S + J   strand 1 stayed
///         max_{a != a'} delta(a,b')           + J + S   strand 2 stayed
///         max_{a != a', b != b'} delta(a,b)   + 2J      both jumped
///
/// Every leave-one-out maximum comes from a top-2 along the relevant axis, so this stays O(m^2)
/// like the forward step rather than O(m^4). Checked against a literal O(m^4) implementation over
/// random emissions including infeasible states.
///
/// In logs, unlike the forward pass: sum-product needs rescaling per site to avoid underflow,
/// max-product does not, and in logs the stay-or-jump choice is a comparison of sums.
void viterbi_step(const vector<double>& in, size_t m, double rho, const vector<double>& emission,
                  vector<double>& out, vector<uint16_t>& back_a, vector<uint16_t>& back_b) {
    double stay = 1.0 - rho + rho / (double)m;
    double jump = rho / (double)m;
    double S = stay > 0.0 ? log(stay) : NEG_INF;
    double J = jump > 0.0 ? log(jump) : NEG_INF;

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
                double v = c1 + S + S;
                if (v > best) { best = v; ba = ap; bb = bp; }
            }
            double c2 = rowExcl[ap * m + bp];
            if (c2 > NEG_INF) {
                double v = c2 + S + J;
                if (v > best) { best = v; ba = ap; bb = rowExclArg[ap * m + bp]; }
            }
            bool chit = (cols[bp].arg == ap);
            double c3 = chit ? cols[bp].second : cols[bp].best;
            if (c3 > NEG_INF) {
                double v = c3 + J + S;
                if (v > best) { best = v; ba = chit ? cols[bp].arg2 : cols[bp].arg; bb = bp; }
            }
            bool bhit = (both.arg == ap);
            double c4 = bhit ? both.second : both.best;
            if (c4 > NEG_INF) {
                double v = c4 + J + J;
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
                                     vector<vector<double>>& out) const {
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
        vector<double> a(m * m, 1.0 / (double)(m * m));
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
        size_t gap = sites[from + t].position > sites[from + t - 1].position
                         ? sites[from + t].position - sites[from + t - 1].position : 1;
        vector<double> moved;
        transition_apply(alpha[t - 1], m, switch_probability(gap), moved);
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
    for (size_t t = n; t-- > 0;) {
        const Site& site = sites[from + t];
        vector<double>& post = out[from + t];
        post.assign(site.genotype_ln_likelihood.size(), 0.0);
        vector<size_t> multiplicity(site.genotype_ln_likelihood.size(), 0);

        // Two accumulators, kept apart on purpose. `known` is mass from states where both
        // strands name an allele, and it carries the multiplicity that becomes the
        // allele-frequency prior. `wild` is mass from states with a latent allele, which has no
        // multiplicity and must not be divided by one -- conflating the two made a flat site
        // asymmetric between genotypes whose only difference was how many haplotypes spelled them.
        vector<double> known(site.genotype_ln_likelihood.size(), 0.0);
        vector<double> wild(site.genotype_ln_likelihood.size(), 0.0);
        size_t n_alleles = site.num_alleles;

        // The frequency prior leaks through a second channel, which a test caught: a state
        // pairing a panel haplotype with the wildcard also occurs once per carrying haplotype, so
        // the count of carriers weights the half-wildcard mass exactly as multiplicity weights the
        // known-known mass. Neutralising only the first left a flat site tilted toward the common
        // allele. Both have to be divided out for freq_prior = 0 to mean what it says.
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
                // A latent allele is marginalised **conditional on the reads**, not spread
                // uniformly. Spreading it uniformly ignored the emission entirely, and at a site
                // whose reads were decisive for an allele the panel does not carry it handed the
                // win to a genotype the reads had ruled out.
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
                        // No `freq_prior < 1` guard. At exactly 1 the exponent is 0 and this is a
                        // no-op, so the guard cost nothing and looked free -- but it also made
                        // every value above 1 behave as 1, silently. Above 1 the exponent goes
                        // negative and the prior is amplified past the state space's own
                        // multiplicity, which is a real setting and was unreachable.
                        share /= pow((double)carriers[k], 1.0 - params.freq_prior);
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

        // Divide out the allele-frequency prior the state space implies, on the known-known mass
        // only. See the note on Params::freq_prior: over a panel selected against these reads it
        // is the same evidence twice.
        double total = 0.0;
        for (size_t idx = 0; idx < post.size(); ++idx) {
            double k = known[idx];
            if (multiplicity[idx] > 1) {
                // See the note on the carriers guard above: `freq_prior < 1` here silently
                // clamped every larger value to 1.
                k /= pow((double)multiplicity[idx], 1.0 - params.freq_prior);
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
        size_t gap = sites[from + t].position > sites[from + t - 1].position
                         ? sites[from + t].position - sites[from + t - 1].position : 1;
        vector<double> weighted(m * m, 0.0);
        for (size_t k = 0; k < m * m; ++k) {
            weighted[k] = beta[k] * emissions[t][k];
        }
        vector<double> next;
        transition_apply(weighted, m, switch_probability(gap), next);
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

vector<vector<double>> LinkageModel::posteriors(const vector<Site>& sites) const {
    vector<vector<double>> out(sites.size());
    if (sites.empty()) {
        return out;
    }
    size_t n = sites.size();
    size_t step = max<size_t>(params.window, 1);
    size_t margin = params.margin;
    if (n <= step) {
        window_posteriors(sites, 0, n, out);
        return out;
    }
    // Overlapping windows, keeping only interiors, so no posterior is read from a position that
    // can see an artificial window edge.
    for (size_t start = 0; start < n; start += step) {
        size_t lo = start > margin ? start - margin : 0;
        size_t hi = min(start + step + margin, n);
        vector<vector<double>> local(n);
        window_posteriors(sites, lo, hi, local);
        size_t keep_to = min(start + step, n);
        for (size_t t = start; t < keep_to; ++t) {
            out[t] = std::move(local[t]);
        }
    }
    return out;
}

vector<LinkageModel::Phase> LinkageModel::phasing(const vector<Site>& sites,
                                                  const vector<size_t>& constraint) const {
    vector<Phase> out(sites.size());
    if (sites.empty()) {
        return out;
    }
    size_t n = sites.size();
    size_t n_hap = 0;
    for (const Site& s : sites) {
        n_hap = max(n_hap, s.haplotype_allele.size());
    }
    size_t m = n_hap + 1;
    if (m < 2) {
        return out;
    }

    size_t step = max<size_t>(params.window, 1);
    size_t margin = params.margin;

    // Windowing needs more care here than it does for the posteriors. A marginal posterior is a
    // per-site quantity, so windows can be decoded independently and their interiors pasted
    // together. A path is not: decode two windows independently and the state at the seam is
    // chosen twice, by two runs that never saw each other, so the join manufactures a switch at
    // every window boundary -- an artifact that would read as a biological result and would scale
    // with the site count rather than with the genome.
    //
    // So each window after the first is *pinned*: the state at one index inside its leading
    // margin is fixed to what the previous window already decided there. The pin sits in the
    // margin rather than at the boundary precisely because a margin-deep state was decided with
    // the full window's context around it, while a boundary state was not.
    bool have_pin = false;
    size_t pin_index = 0;
    Phase pin{};

    for (size_t start = 0; start < n; start += step) {
        size_t lo = start > margin ? start - margin : 0;
        size_t hi = min(start + step + margin, n);
        vector<Phase> local;
        window_phasing(sites, lo, hi, constraint,
                       have_pin ? pin_index : (size_t)-1, pin, local);
        size_t keep_to = min(start + step, n);
        for (size_t t = start; t < keep_to && t - lo < local.size(); ++t) {
            out[t] = local[t - lo];
        }
        if (keep_to == n) {
            break;
        }
        // Pin the next window at the last index this one is authoritative for.
        pin_index = keep_to - 1;
        pin = out[pin_index];
        have_pin = true;
    }
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
                // The wildcard, and a haplotype absent from this site, may carry any allele --
                // which is what keeps the constrained problem feasible where the panel cannot
                // spell a call. But "free" means the *free* strand is unconstrained, not that the
                // state is: a known haplotype carrying neither of the wanted alleles cannot be
                // rescued by pairing it with a wildcard. Getting that wrong silently emits a
                // phasing that contradicts the very genotype it was constrained to.
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
                // Every surviving state implies the *same* genotype, so they must all carry that
                // genotype's likelihood. build_emission gave the free strands a marginal over
                // alleles, which is right when the genotype is unknown and wrong here: the
                // constraint determines what the free strand carries. Left as a marginal it
                // reliably exceeds the constrained genotype's own likelihood wherever the reads
                // disagree with the call, so the path would route through the wildcard precisely
                // at the sites the constraint exists to pin down.
                //
                // The escape factor stays, one per free strand, so a genotype the panel can spell
                // is still preferred over the wildcard explaining it.
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
            // leave the site with no reachable state and the chain with no path. Leave it merely
            // constrained. Not expected -- the pin comes from a path that was constrained to the
            // same genotype -- but a silent dead chain is a bad way to find out otherwise.
            continue;
        }
        for (size_t a = 0; a < m; ++a) {
            for (size_t b = 0; b < m; ++b) {
                if (a != pa || b != pb) {
                    emissions[t][a * m + b] = 0.0;
                }
            }
        }
    }

    // Pin: everything except the pinned state becomes unreachable at that index.
    //
    // Tested before anything is zeroed, exactly as the per-site pin above does: if the pinned pair
    // is impossible under this window's constraints, the pin is skipped and the site left merely
    // constrained. The earlier shape zeroed first and then "fell back" by setting the forbidden
    // pair to 1.0 -- the opposite of leaving the site free, routing every surviving path through
    // the one state the constraint had ruled out and biasing the kept region's leading sites
    // toward the stale orientation. (With margin == 0 the pin index falls outside [from, to) and
    // the pin is skipped entirely; vg call always runs with a margin, but the Params permit 0.)
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
        size_t gap = sites[from + t].position > sites[from + t - 1].position
                         ? sites[from + t].position - sites[from + t - 1].position : 1;
        vector<double> next;
        viterbi_step(delta, m, switch_probability(gap), emissions[t], next, back_a[t], back_b[t]);
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
// Haploid chains
//
// chrX outside the pseudoautosomal regions, and all of chrY, are haploid in a male sample. The
// diploid model cannot express them -- its state is a pair -- so before this they were dropped
// from the linkage pass entirely, silently, taking their mosaic with them. That is about 5% of a
// genome, and the part where "which haplotype explains this stretch" is the entire answer.

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
                                             vector<vector<double>>& out) const {
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
        vector<double> a(m, 1.0 / (double)m);
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
        size_t gap = sites[from + t].position > sites[from + t - 1].position
                         ? sites[from + t].position - sites[from + t - 1].position : 1;
        double rho = switch_probability(gap);
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
            // A latent allele is marginalised conditional on the reads, not spread uniformly --
            // the same correction the diploid pass needed.
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
                v /= pow((double)carriers[k], 1.0 - params.freq_prior);
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
        size_t gap = sites[from + t].position > sites[from + t - 1].position
                         ? sites[from + t].position - sites[from + t - 1].position : 1;
        double rho = switch_probability(gap);
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

vector<vector<double>> LinkageModel::haploid_posteriors(const vector<Site>& sites) const {
    vector<vector<double>> out(sites.size());
    if (sites.empty()) {
        return out;
    }
    size_t n = sites.size();
    size_t step = max<size_t>(params.window, 1);
    size_t margin = params.margin;
    if (n <= step) {
        window_haploid_posteriors(sites, 0, n, out);
        return out;
    }
    for (size_t start = 0; start < n; start += step) {
        size_t lo = start > margin ? start - margin : 0;
        size_t hi = min(start + step + margin, n);
        vector<vector<double>> local(n);
        window_haploid_posteriors(sites, lo, hi, local);
        size_t keep_to = min(start + step, n);
        for (size_t t = start; t < keep_to; ++t) {
            out[t] = std::move(local[t]);
        }
    }
    return out;
}

void LinkageModel::window_haploid_phasing(const vector<Site>& sites, size_t from, size_t to,
                                          const vector<size_t>& constraint,
                                          size_t pin_index, size_t pin,
                                          vector<size_t>& out) const {
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
    for (size_t k = 0; k < m; ++k) {
        if (emissions[0][k] > 0.0) {
            delta[k] = -log((double)m) + log(emissions[0][k]);
        }
    }
    vector<vector<uint16_t>> back(n);
    for (size_t t = 1; t < n; ++t) {
        size_t gap = sites[from + t].position > sites[from + t - 1].position
                         ? sites[from + t].position - sites[from + t - 1].position : 1;
        double rho = switch_probability(gap);
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

vector<size_t> LinkageModel::haploid_phasing(const vector<Site>& sites,
                                             const vector<size_t>& constraint) const {
    vector<size_t> out(sites.size(), WILDCARD);
    if (sites.empty()) {
        return out;
    }
    size_t n = sites.size();
    size_t step = max<size_t>(params.window, 1);
    size_t margin = params.margin;
    bool have_pin = false;
    size_t pin_index = 0, pin = WILDCARD;

    for (size_t start = 0; start < n; start += step) {
        size_t lo = start > margin ? start - margin : 0;
        size_t hi = min(start + step + margin, n);
        vector<size_t> local;
        window_haploid_phasing(sites, lo, hi, constraint,
                               have_pin ? pin_index : (size_t)-1, pin, local);
        size_t keep_to = min(start + step, n);
        for (size_t t = start; t < keep_to && t - lo < local.size(); ++t) {
            out[t] = local[t - lo];
        }
        if (keep_to == n) {
            break;
        }
        // Pinned a margin deep into what this window was authoritative for, for the same reason
        // the diploid pass does it: a boundary state was decided without full context.
        pin_index = keep_to - 1;
        pin = out[pin_index];
        have_pin = true;
    }
    return out;
}

//------------------------------------------------------------------------------
// LinkageCollector

vector<int> LinkageCollector::compact_allele_space(
        const map<vector<int>, double>& genotype_ln_likelihood,
        const vector<int>& haplotype_traversal,
        int called_trav_i, int called_trav_j) {
    // The reachable set, and nothing more. `build_emission` only ever indexes the likelihood vector
    // at a pair of *panel-carried* alleles, and the genotype constraint needs only the called pair,
    // so a traversal that is neither cannot be spelled by any panel pair and cannot affect the
    // result. Dropping the rest is what keeps this space the same order of size as the emitted-allele
    // space it replaces, instead of the full candidate list -- a 34-haplotype panel offers ~35
    // candidates, whose triangular likelihood vector would be 630 entries against today's handful.
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

void LinkageCollector::record(const string& contig, size_t position,
                              const map<vector<int>, double>& genotype_ln_likelihood,
                              const vector<int>& haplotype_traversal,
                              int called_trav_i, int called_trav_j,
                              const vector<int>& traversal_to_allele,
                              size_t record_key,
                              double explained_share, size_t ploidy,
                              int64_t start_node, int64_t end_node,
                              bool nested, size_t parent_record_key, int parent_trav,
                              uint64_t parent_crossing, size_t generation,
                              bool emitted) {
    if (genotype_ln_likelihood.empty() || called_trav_i < 0) {
        return;
    }
    vector<int> space = compact_allele_space(genotype_ln_likelihood, haplotype_traversal,
                                             called_trav_i, called_trav_j);
    if (space.empty() || space.size() > 127) {
        // int8 in the panel arena caps a site at 127 alleles; one that exceeded it would lose
        // linkage rather than be mis-linked, which is the safe direction.
        return;
    }
    const size_t k = space.size();
    // candidate traversal -> compact allele
    map<int, int> compact;
    for (size_t i = 0; i < k; ++i) {
        compact[space[i]] = (int)i;
    }
    auto compact_of = [&](int trav) -> int {
        auto it = compact.find(trav);
        return it == compact.end() ? -1 : it->second;
    };

    const int ci = compact_of(called_trav_i);
    const int cj = called_trav_j >= 0 ? compact_of(called_trav_j) : ci;
    if (ci < 0 || cj < 0) {
        return;   // cannot happen: the called pair is in the space by construction
    }

    const size_t site_ploidy = (ploidy == 1 ? 1 : 2);
    const size_t n_gt = site_ploidy == 1 ? k : k * (k + 1) / 2;
    vector<float> gls(n_gt, -numeric_limits<float>::infinity());
    for (const auto& kv : genotype_ln_likelihood) {
        if (kv.first.size() != site_ploidy) {
            continue;
        }
        if (site_ploidy == 1) {
            int a = compact_of(kv.first[0]);
            if (a >= 0) {
                gls[(size_t)a] = (float)kv.second;
            }
            continue;
        }
        int a = compact_of(kv.first[0]);
        int b = compact_of(kv.first[1]);
        if (a >= 0 && b >= 0) {
            gls[LinkageModel::genotype_index((size_t)a, (size_t)b)] = (float)kv.second;
        }
    }

    lock_guard<std::mutex> guard(mutex);

    uint32_t contig_id = 0;
    bool found = false;
    for (size_t i = 0; i < contig_names.size(); ++i) {
        if (contig_names[i] == contig) {
            contig_id = (uint32_t)i;
            found = true;
            break;
        }
    }
    if (!found) {
        contig_id = (uint32_t)contig_names.size();
        contig_names.push_back(contig);
    }

    Entry e;
    e.position = (uint32_t)position;
    e.contig = contig_id;
    e.num_alleles = (uint16_t)k;
    e.called_i = (uint16_t)ci;
    e.called_j = (uint16_t)cj;
    e.record_key = record_key;
    e.explained_share = (float)explained_share;
    e.start_node = start_node;
    e.end_node = end_node;
    e.ploidy = (uint8_t)site_ploidy;
    e.nested = nested;
    e.emitted = emitted;
    e.parent_record_key = parent_record_key;
    e.parent_crossing = parent_crossing;
    e.parent_trav = (int16_t)parent_trav;
    e.generation = (uint8_t)(generation > 255 ? 255 : generation);

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
    entries.push_back(e);
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
    return entries.size() * sizeof(Entry)
           + gl_arena.size() * sizeof(float)
           + hap_arena.size() * sizeof(int8_t);
}

bool LinkageCollector::respecify(size_t record_key,
                                 const string& contig, size_t position,
                                 const map<vector<int>, double>& genotype_ln_likelihood,
                                 const vector<int>& haplotype_traversal,
                                 int called_trav_i, int called_trav_j,
                                 const vector<int>& traversal_to_allele,
                                 size_t ploidy, bool nested, int parent_trav,
                                 uint64_t parent_crossing, bool emitted) {
    if (genotype_ln_likelihood.empty() || called_trav_i < 0) {
        return false;
    }
    // The same compact space `record` builds, from the same inputs, so a revised site is described
    // exactly as a freshly recorded one would be.
    vector<int> space = compact_allele_space(genotype_ln_likelihood, haplotype_traversal,
                                            called_trav_i, called_trav_j);
    if (space.empty() || space.size() > 127) {
        return false;
    }
    const size_t k = space.size();
    map<int, int> compact;
    for (size_t i = 0; i < k; ++i) {
        compact[space[i]] = (int)i;
    }
    auto compact_of = [&](int trav) -> int {
        auto it = compact.find(trav);
        return it == compact.end() ? -1 : it->second;
    };
    const int ci = compact_of(called_trav_i);
    const int cj = called_trav_j >= 0 ? compact_of(called_trav_j) : ci;
    if (ci < 0 || cj < 0) {
        return false;
    }
    const size_t site_ploidy = (ploidy == 1 ? 1 : 2);
    const size_t n_gt = site_ploidy == 1 ? k : k * (k + 1) / 2;
    vector<float> gls(n_gt, -numeric_limits<float>::infinity());
    for (const auto& kv : genotype_ln_likelihood) {
        if (kv.first.size() != site_ploidy) {
            continue;
        }
        if (site_ploidy == 1) {
            int a = compact_of(kv.first[0]);
            if (a >= 0) {
                gls[(size_t)a] = (float)kv.second;
            }
            continue;
        }
        int a = compact_of(kv.first[0]);
        int b = compact_of(kv.first[1]);
        if (a >= 0 && b >= 0) {
            gls[LinkageModel::genotype_index((size_t)a, (size_t)b)] = (float)kv.second;
        }
    }

    lock_guard<std::mutex> guard(mutex);
    for (Entry& e : entries) {
        if (e.record_key != record_key || e.retracted) {
            continue;
        }
        // Every vector for the new ploidy is a different length, so none can be written over the old
        // one in place. Appended and the offsets repointed: the arenas only ever grow, and the
        // abandoned spans are a few entries per revised site.
        e.num_alleles = (uint16_t)k;
        e.called_i = (uint16_t)ci;
        e.called_j = (uint16_t)cj;
        e.ploidy = (uint8_t)site_ploidy;
        e.nested = nested;
        // Whether a line exists is a property of the revision, not of the original record. A chain
        // the sweep wrote nothing for is recorded unemitted and may become a real record here; not
        // updating this left it invisible to both the genotype patch and the phase patch.
        e.emitted = emitted;
        // Re-key onto the line that now exists. Re-emitting at a different ploidy changes the
        // emitted allele set, and POS is advanced by the prefix every allele shares, so the
        // replacement is filed under a position the sweep's entry does not name. The patch indices
        // are keyed on (contig, POS), so leaving this alone did not *decline* the patches -- it
        // meant they were never looked up at all.
        {
            uint32_t contig_id = 0;
            bool found_contig = false;
            for (size_t i = 0; i < contig_names.size(); ++i) {
                if (contig_names[i] == contig) {
                    contig_id = (uint32_t)i;
                    found_contig = true;
                    break;
                }
            }
            if (!found_contig) {
                contig_id = (uint32_t)contig_names.size();
                contig_names.push_back(contig);
            }
            e.contig = contig_id;
            e.position = (uint32_t)position;
        }
        e.parent_trav = (int16_t)parent_trav;
        e.parent_crossing = parent_crossing;
        e.gl_offset = (uint32_t)gl_arena.size();
        for (float v : gls) {
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
            int allele = (size_t)space[i] < traversal_to_allele.size()
                             ? traversal_to_allele[space[i]] : -1;
            allele_arena.push_back(allele >= 0 && allele < 127 ? (int8_t)allele : (int8_t)-1);
        }
        return true;
    }
    return false;
}

bool LinkageCollector::set_parent_trav(size_t record_key, int parent_trav) {
    lock_guard<std::mutex> guard(mutex);
    for (Entry& e : entries) {
        if (e.record_key == record_key && !e.retracted) {
            e.parent_trav = (int16_t)parent_trav;
            return true;
        }
    }
    return false;
}

bool LinkageCollector::has_entry(size_t record_key) const {
    lock_guard<std::mutex> guard(mutex);
    for (const Entry& e : entries) {
        if (e.record_key == record_key && !e.retracted) {
            return true;
        }
    }
    return false;
}

bool LinkageCollector::retract(size_t record_key) {
    lock_guard<std::mutex> guard(mutex);
    for (Entry& e : entries) {
        if (e.record_key == record_key && !e.retracted) {
            e.retracted = true;
            return true;
        }
    }
    return false;
}

vector<LinkageCollector::Change> LinkageCollector::resolve_generation(
        size_t generation, bool last, vector<PhaseCall>* phasing_out) {
    // Nested sites, collected across every contig and phased after the diploid chains, so their
    // strand can be read off the parent they hang off rather than guessed.
    vector<size_t> deferred_nested;
    vector<Change> changes;
    // Genotypes the model settled on that the emitted record carries no ALT for. Reachable only
    // because the collector's allele space is the genotyper's rather than the VCF's, so it is
    // counted and reported rather than assumed away.
    size_t unrenderable = 0;
    // Records whose phase names the called pair because the settled pair had no ALT. The phase is
    // still real -- the block and the strand order come from the panel either way -- but the alleles
    // it names are the line's rather than the model's, so the count belongs in the report.
    size_t phase_fallback = 0;
    if (!model.active() || entries.empty()) {
        return changes;
    }

    // The phase an earlier generation settled on, by record key, so a clamped site can be pinned to
    // the pair the VCF already carries. Read before the chain loop appends this generation's own.
    struct PinnedPhase { size_t first; size_t second; };
    unordered_map<size_t, PinnedPhase> pinned_phase;
    if (phasing_out != nullptr && generation > 0) {
        pinned_phase.reserve(phasing_out->size() * 2);
        for (const PhaseCall& pc : *phasing_out) {
            pinned_phase[pc.record_key] = PinnedPhase{pc.hap_first, pc.hap_second};
        }
    }

    // Default every site this pass considers to its own per-site call, so that whatever a chain or a
    // sweep fails to reach still has a coherent genotype for a later generation to clamp. Overwritten
    // below wherever something is actually settled.
    for (Entry& e : entries) {
        if (e.generation != generation) {
            continue;
        }
        e.final_i = e.called_i;
        e.final_j = e.ploidy == 1 ? e.called_i : e.called_j;
    }

    // Group by contig, then sort by reference position. Node-ID order is close to reference order
    // in a reference-first graph but is not guaranteed to be it, and the transition probabilities
    // are computed from the gaps -- so trusting the arrival order would silently feed the model
    // the wrong distances.
    vector<vector<size_t>> by_contig(contig_names.size());
    for (size_t i = 0; i < entries.size(); ++i) {
        if (entries[i].generation > generation) {
            continue;   // not called yet: this generation's descent has not reached it
        }
        if (entries[i].retracted) {
            continue;   // the settled parent does not carry the chain, so there is no site here
        }
        by_contig[entries[i].contig].push_back(i);
    }

    // A chain is a maximal run of one ploidy on one contig, not a whole contig.
    //
    // It used to be a whole contig, on the reasoning that a contig is called at one ploidy and
    // chrX's pseudoautosomal split is two separate runs. --ploidy-bed makes that false: a single
    // run can now carry diploid pseudoautosomal regions and a haploid interior. Splitting is the
    // correct answer rather than a workaround -- the transition model moves probability between
    // adjacent sites through pairs of panel haplotypes, and at a ploidy boundary there is no
    // correspondence to carry across, exactly as there is none between two contigs. It also
    // reproduces what the two-pass splice did before, since each side becomes its own chain with
    // its own phase set.
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
        // Sorted first, so a run is contiguous in reference order rather than in arrival order.
        // Nested sites are held back from the runs entirely and phased afterwards, against the
        // parent they hang off. Letting them into the runs cuts a chain at every one of them,
        // because a nested ploidy-1 site between diploid neighbours is a ploidy change: chr20's
        // autosomal phasing went from 22 blocks to 9,460 and block N50 from 248 Mb to 1.08 Mb, and
        // the switch rate only looked flat because short blocks make switch error cheap.
        //
        // A *regional* ploidy change still cuts, which is the case the rule exists for: across
        // chrX's pseudoautosomal boundary there is no haplotype correspondence to carry, so joining
        // those runs would be wrong rather than merely fragmented.
        vector<size_t> nested_here;
        vector<size_t> chainable;
        chainable.reserve(contig_indices.size());
        for (size_t idx : contig_indices) {
            if (entries[idx].nested) {
                // Only this generation's. An earlier generation's nested site is already placed on a
                // strand and already resolved in its own per-strand chain, and it is held out of the
                // diploid runs by construction -- so there is nothing for it to contribute here, and
                // re-placing it would rewrite a PhaseCall the VCF has already been told about.
                if (entries[idx].generation == generation) {
                    nested_here.push_back(idx);
                }
            } else {
                chainable.push_back(idx);
            }
        }
        deferred_nested.insert(deferred_nested.end(), nested_here.begin(), nested_here.end());

        size_t run_start = 0;
        while (run_start < chainable.size()) {
            size_t run_end = run_start + 1;
            while (run_end < chainable.size()
                   && entries[chainable[run_end]].ploidy
                          == entries[chainable[run_start]].ploidy) {
                ++run_end;
            }
            chains.emplace_back(chainable.begin() + run_start,
                                chainable.begin() + run_end);
            run_start = run_end;
        }
    }

    for (auto& indices : chains) {
        if (indices.empty()) {
            continue;
        }
        // A one-site chain has nothing to link to, so linkage cannot move its genotype -- but it
        // still has to be *phased*, or it never reaches phasing_out and the mosaic stops accounting
        // for every emitted record. Skipping it outright was invisible while chains were maximal
        // runs of one ploidy along a contig, which made singletons vanishingly rare. Propagated
        // ploidy creates them in quantity: an isolated ploidy-1 site between diploid neighbours is
        // a chain of one, and 258 of them went missing from the chr20 mosaic.

        // Every entry in a chain shares a ploidy by construction above.
        size_t chain_ploidy = entries[indices.front()].ploidy;

        vector<LinkageModel::Site> sites;
        sites.reserve(indices.size());
        for (size_t idx : indices) {
            const Entry& e = entries[idx];
            LinkageModel::Site s;
            s.position = e.position;
            s.num_alleles = e.num_alleles;
            s.ploidy = e.ploidy;
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

        vector<vector<double>> posteriors = chain_ploidy == 1 ? model.haploid_posteriors(sites)
                                                             : model.posteriors(sites);

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
            if (best == before) {
                continue;
            }
            if (!e.emitted) {
                continue;   // phased, and it moved, but there is no line to rewrite
            }
            // Rendered into VCF numbering here and nowhere else. All four alleles have to survive
            // the map: a genotype the model settled on can name a traversal the record carries no
            // ALT for, and the patch cannot add one, so such a change is dropped and counted rather
            // than written against the wrong allele.
            const int r_called_i = vcf_allele_of(allele_arena, e.allele_offset, e.num_alleles, e.called_i);
            const int r_called_j = vcf_allele_of(allele_arena, e.allele_offset, e.num_alleles, e.called_j);
            const int r_i = vcf_allele_of(allele_arena, e.allele_offset, e.num_alleles, i);
            const int r_j = vcf_allele_of(allele_arena, e.allele_offset, e.num_alleles, j);
            if (r_called_i < 0 || r_called_j < 0 || r_i < 0 || r_j < 0) {
                ++unrenderable;
                continue;
            }
            Change c;
            c.record_key = e.record_key;
            c.contig = contig_names[e.contig];
            c.position = e.position;
            c.called_i = (size_t)r_called_i;
            c.called_j = (size_t)r_called_j;
            c.allele_i = (size_t)r_i;
            c.allele_j = (size_t)r_j;
            c.posterior = post[best];
            c.explained_share = (double)e.explained_share;
            changes.push_back(c);
        }

        if (phasing_out == nullptr) {
            continue;
        }
        // Haploid chains take the single-strand path: one haplotype per site, no phase to infer,
        // but the mosaic is exactly as meaningful -- it is the whole answer on chrY.
        vector<LinkageModel::Phase> phase;
        if (chain_ploidy == 1) {
            vector<size_t> single = model.haploid_phasing(sites, final_genotype);
            phase.resize(single.size());
            for (size_t t = 0; t < single.size(); ++t) {
                phase[t].first = single[t];
                phase[t].second = LinkageModel::WILDCARD;
            }
        } else {
            phase = model.phasing(sites, final_genotype);
        }
        // One phase set per contig. The windows are pinned to each other, so the path is
        // continuous across the whole chain and the block is the chain -- far longer than a
        // read-based phaser produces, which is a claim the switch-error benchmark is there to
        // test rather than something to assert quietly.
        size_t phase_set = sites.empty() ? 0 : sites.front().position;
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
            pc.emitted = e.emitted;
            pc.record_key = e.record_key;
            pc.contig = contig_names[e.contig];
            pc.position = e.position;
            pc.hap_first = ph.first;
            pc.hap_second = ph.second;
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
            // The compact pair is the collector's own numbering. Split it here into the two things that
            // consume it: the traversal on each strand, which is the genome fact a crossing mask and a
            // haplotype path are expressed in, and the VCF allele on each strand, which is the only form
            // apply_phasing can write. Leaving allele_* compact made the phased-GT guard reject the pair and
            // silently drop the record's phasing -- 1,528 extra strandless records on chr20.
            {
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
                    ++phase_fallback;
                    pc.order_arbitrary = pc.order_arbitrary || (v_first != v_second);
                }
            }
            phasing_out->push_back(pc);
        }
    }

    // Nested sites, phased against the parent they hang off rather than in a chain of their own.
    //
    // A nested site has ploidy 1 exactly because one called parent allele crosses the child chain
    // and the other does not, so which strand it belongs to is *determined* once the parent is
    // phased -- it is not a fit to be estimated. `parent_trav` records which of the parent's two
    // genotype slots did the crossing, and the parent's own PhaseCall says which strand that slot
    // landed on, so the nested allele goes there and the other strand carries the wildcard.
    //
    // Done as a separate pass because it needs the parent already phased, and the parent may be
    // anywhere in any chain.
    if (phasing_out != nullptr && !deferred_nested.empty()) {
        // By value, not by pointer: this loop appends to `phasing_out`, and any pointer into a
        // vector is invalid the moment it reallocates. Only three fields are needed.
        struct ParentPhase {
            size_t phase_set;
            size_t hap_first;
            size_t hap_second;
            // The phased pair, in the two numberings that matter. `allele_*` is compact -- the
            // collector's own space -- and `trav_*` is the candidate traversal each strand is on.
            //
            // The crossing mask is in traversal terms, so the mask must be tested against `trav_*`.
            // Testing it against the compact index is only right when every allele at the site is
            // panel-carried, which is why the unit tests caught it: a site whose panel names one
            // allele has a one-element compact space, so compact 0 is traversal 1 and the bit tested
            // was the wrong one.
            size_t allele_first;
            size_t allele_second;
            int trav_first;
            int trav_second;
            size_t ploidy;
        };
        unordered_map<size_t, ParentPhase> by_key;
        by_key.reserve(phasing_out->size() * 2);
        for (const PhaseCall& pc : *phasing_out) {
            by_key[pc.record_key] = ParentPhase{pc.phase_set, pc.hap_first, pc.hap_second,
                                                pc.allele_first, pc.allele_second,
                                                pc.trav_first, pc.trav_second, pc.ploidy};
        }
        // Nested sites whose parent settled on a pair not containing the traversal they hang off.
        // The barrier retracts these where it can see them, so a nonzero count here is the residue
        // whose parent moved after it looked -- not a coherence disagreement, which the single
        // derivation above makes unrepresentable.
        size_t unplaced_no_strand = 0;
        // Nested sites the settled parent carries on both its traversals: the locus is diploid
        // there and the record names one allele. The barrier re-renders these at ploidy 2 wherever
        // an answer at that ploidy was kept, so this is the residue where none was.
        size_t carried_on_both = 0;
        // Split by whether a VCF line exists, because since collapsed sites started being recorded
        // these counters see line-less entries too, and a figure that mixes records with entries
        // that were never records cannot be read as a defect count.
        size_t carried_on_both_emitted = 0, unplaced_no_strand_emitted = 0;
        // Nested sites that never found a phased parent, so no strand could be derived for them at
        // all. Counted because they are the population that comes out with a bare haploid GT, and
        // the report below used to describe only the sites whose parent *was* found.
        size_t unplaced = 0;
        // Which strand each nested site was placed on, so step three can group them. Filled as the
        // sweeps below resolve parents.
        map<pair<uint32_t, uint8_t>, vector<size_t>> by_strand;
        /// record_key -> what step three settled on, for sites whose genotype it moved. All three
        /// spaces are kept because all three are needed and they are not interchangeable: the
        /// traversal is what a haplotype path walks, the VCF allele is what a GT can name (-1 when
        /// the record carries no ALT for it), and the compact index is what the arenas are keyed by.
        struct Regenotyped {
            size_t compact = 0;
            int vcf_allele = -1;
            int traversal = -1;
            /// The line's own allele, for when the settled one has no ALT -- see render_phase_pair.
            int called_allele = -1;
        };
        unordered_map<size_t, Regenotyped> nested_regenotyped;

        // Resolve shallowest-first so a nested site whose parent is itself nested finds its parent
        // already placed. Depth is bounded by the snarl hierarchy, so a bounded number of sweeps
        // settles it; anything still unresolved after that keeps the wildcard rather than guessing.
        vector<size_t> pending = deferred_nested;
        for (int sweep = 0; sweep < 8 && !pending.empty(); ++sweep) {
            vector<size_t> still;
            size_t placed_this_sweep = 0;
            for (size_t idx : pending) {
                const Entry& e = entries[idx];
                auto found = by_key.find(e.parent_record_key);
                PhaseCall pc;
                pc.ploidy = e.ploidy;
                pc.emitted = e.emitted;
                pc.record_key = e.record_key;
                pc.contig = contig_names[e.contig];
                pc.position = e.position;
                pc.start_node = e.start_node;
                pc.end_node = e.end_node;
                pc.allele_first = e.called_i;
                pc.allele_second = e.called_i;
                if (found == by_key.end()) {
                    still.push_back(idx);
                    continue;
                }
                const ParentPhase& parent = found->second;
                // Which of the parent's strands carries this chain. One lookup, not a derivation:
                // the barrier already decided *which traversal* of the parent's settled pair carries
                // it, and this asks only which haplotype that traversal ended up on.
                //
                // That is the whole of stage 3. Ploidy used to come from `child_ploidy` over the
                // parent's pre-linkage genotype, the strand from a recorded slot index, and their
                // agreement from a third test of the crossing mask against the settled pair -- three
                // readings of one fact, which is exactly why they could disagree and why there were
                // three FILTERs to say so. Here the chain is carried by this traversal, so it has
                // one copy and sits on that traversal's strand, and the two cannot come apart.
                //
                // The traversal rather than the slot index matters: the Viterbi decides which of the
                // parent's traversals lands on which haplotype, and it may decide differently as
                // later generations enlarge the set, so an index can go stale where an identity
                // cannot.
                // -2 is the parent carrying the chain on *both* its settled traversals: the locus
                // is diploid there and this record names one allele where it has two. It keeps no
                // single strand, but it is not strandless either -- both haplotypes carry it -- and
                // the mosaic must say so rather than mark it unexplained.
                const bool on_both = (e.parent_trav == -2);
                int strand = -1;
                if (e.parent_trav >= 0) {
                    if (parent.trav_first == e.parent_trav) {
                        strand = 0;
                    } else if (parent.ploidy == 2 && parent.trav_second == e.parent_trav) {
                        strand = 1;
                    }
                }
                // Inherit the parent's phase set, so the nested site sits in the parent's block
                // instead of starting one of its own -- which is the whole point of the exercise.
                pc.phase_set = parent.phase_set;
                if (on_both) {
                    // Both haplotypes carry it, so both name one. A wildcard here would read as
                    // "the panel cannot explain this strand", which is the opposite of the truth:
                    // the parent explains both. The record's own GT still names one allele, because
                    // it was genotyped at ploidy 1 and the second copy's allele was never scored --
                    // the barrier re-renders such a chain at ploidy 2 whenever an answer at that
                    // ploidy was kept, and this is the residue where none was.
                    ++carried_on_both;
                    carried_on_both_emitted += e.emitted;
                    pc.hap_first = parent.hap_first;
                    pc.hap_second = parent.hap_second;
                } else if (strand < 0) {
                    // No settled parent traversal carries it, so there is no haplotype to name:
                    // claiming one would write a variant into the emitted genome that the parent
                    // record does not carry.
                    ++unplaced_no_strand;
                    unplaced_no_strand_emitted += e.emitted;
                    pc.hap_first = LinkageModel::WILDCARD;
                    pc.hap_second = LinkageModel::WILDCARD;
                } else {
                    pc.hap_first = strand == 1 ? LinkageModel::WILDCARD : parent.hap_first;
                    pc.hap_second = strand == 1 ? parent.hap_second : LinkageModel::WILDCARD;
                    // Recorded rather than left to be inferred from which haplotype is the wildcard:
                    // the parent's own strand can be a wildcard too, and then both sides are, which
                    // is indistinguishable from having no strand at all.
                    pc.nested_strand = (int8_t)strand;
                }
                    by_key[pc.record_key] = ParentPhase{
                    pc.phase_set, pc.hap_first, pc.hap_second,
                    pc.allele_first, pc.allele_second,
                    pc.trav_first, pc.trav_second, pc.ploidy};
                // The compact pair is the collector's own numbering. Split it here into the two things that
                // consume it: the traversal on each strand, which is the genome fact a crossing mask and a
                // haplotype path are expressed in, and the VCF allele on each strand, which is the only form
                // apply_phasing can write. Leaving allele_* compact made the phased-GT guard reject the pair and
                // silently drop the record's phasing -- 1,528 extra strandless records on chr20.
                {
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
                        ++phase_fallback;
                    }
                }
                phasing_out->push_back(pc);
                if (strand >= 0) {
                    // Only sites with one strand are grouped: step three chains within a single
                    // haplotype, and a site with no strand belongs to no such chain.
                    by_strand[make_pair(e.contig, (uint8_t)strand)].push_back(idx);
                }
                ++placed_this_sweep;
            }
            if (placed_this_sweep == 0) {
                pending = still;
                break;
            }
            pending = still;
        }
        // Anything left has no reachable phased parent -- most often because the parent collapsed to
        // the reference symbolically and so emitted no record to phase. Giving each of those its own
        // phase set is what fragments the output: it turned chr20's single block into 285.
        //
        // Instead, attach it to the block that already covers it. A nested site lies inside its
        // parent's span, which lies inside the diploid chain, so the nearest preceding phased site on
        // the same contig is in the right block by construction. The strand stays wildcard on both
        // sides: without the parent's phase there is nothing to place it on, and asserting one would
        // be a guess dressed as a call.
        if (!pending.empty()) {
            unplaced = pending.size();
            map<string, vector<pair<size_t, size_t>>> block_by_contig;   // contig -> (pos, phase_set)
            for (const PhaseCall& pc : *phasing_out) {
                block_by_contig[pc.contig].emplace_back(pc.position, pc.phase_set);
            }
            for (auto& kv : block_by_contig) {
                sort(kv.second.begin(), kv.second.end());
            }
            for (size_t idx : pending) {
                const Entry& e = entries[idx];
                PhaseCall pc;
                pc.ploidy = e.ploidy;
                pc.emitted = e.emitted;
                pc.record_key = e.record_key;
                pc.contig = contig_names[e.contig];
                pc.position = e.position;
                pc.start_node = e.start_node;
                pc.end_node = e.end_node;
                pc.allele_first = e.called_i;
                pc.allele_second = e.called_i;
                pc.hap_first = LinkageModel::WILDCARD;
                pc.hap_second = LinkageModel::WILDCARD;
                pc.phase_set = e.position;
                auto found = block_by_contig.find(pc.contig);
                if (found != block_by_contig.end() && !found->second.empty()) {
                    const vector<pair<size_t, size_t>>& blocks = found->second;
                    auto it = upper_bound(blocks.begin(), blocks.end(),
                                          make_pair(pc.position, numeric_limits<size_t>::max()));
                    if (it != blocks.begin()) {
                        pc.phase_set = prev(it)->second;
                    } else if (!blocks.empty()) {
                        pc.phase_set = blocks.front().second;
                    }
                }
                // The compact pair is the collector's own numbering. Split it here into the two things that
                // consume it: the traversal on each strand, which is the genome fact a crossing mask and a
                // haplotype path are expressed in, and the VCF allele on each strand, which is the only form
                // apply_phasing can write. Leaving allele_* compact made the phased-GT guard reject the pair and
                // silently drop the record's phasing -- 1,528 extra strandless records on chr20.
                {
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
                        ++phase_fallback;
                    }
                }
                phasing_out->push_back(pc);
            }
        }

        // Step three: linkage genotype correction for the nested sites, in per-strand haploid chains.
        //
        // Holding them out of the diploid runs is what restored the phase blocks, but it also cost
        // them linkage entirely -- genotype changes fell 8,829 to 8,441 on chr20, so they were being
        // corrected only per-site. Resolving them here gives that back without reintroducing
        // fragmentation, because the phase set was already fixed in step one.
        //
        // Grouped *by strand*, not merely by contig. Nested sites hanging off opposite parent strands
        // are not on the same haplotype, so chaining them together would link sequences that never
        // co-occur -- a worse error than not linking at all. Each group is ploidy-uniform by
        // construction, so the existing haploid model applies unchanged.
        for (auto& group : by_strand) {
            vector<size_t>& idxs = group.second;
            if (idxs.size() < 2) {
                continue;   // nothing to link against
            }
            sort(idxs.begin(), idxs.end(), [&](size_t a, size_t b) {
                if (entries[a].position != entries[b].position) {
                    return entries[a].position < entries[b].position;
                }
                return entries[a].record_key < entries[b].record_key;
            });
            vector<LinkageModel::Site> sites;
            sites.reserve(idxs.size());
            for (size_t idx : idxs) {
                const Entry& e = entries[idx];
                LinkageModel::Site s;
                s.position = e.position;
                s.num_alleles = e.num_alleles;
                s.ploidy = 1;
                size_t n_gt = (size_t)e.num_alleles;
                s.genotype_ln_likelihood.reserve(n_gt);
                for (size_t g = 0; g < n_gt; ++g) {
                    s.genotype_ln_likelihood.push_back((double)gl_arena[e.gl_offset + g]);
                }
                s.haplotype_allele.reserve(n_haplotypes);
                for (size_t h = 0; h < n_haplotypes; ++h) {
                    s.haplotype_allele.push_back((int)hap_arena[e.hap_offset + h]);
                }
                sites.push_back(std::move(s));
            }
            vector<vector<double>> posteriors = model.haploid_posteriors(sites);
            for (size_t k = 0; k < idxs.size() && k < posteriors.size(); ++k) {
                const Entry& e = entries[idxs[k]];
                const vector<double>& post = posteriors[k];
                if (post.empty()) {
                    continue;
                }
                size_t best = 0;
                for (size_t g = 1; g < post.size(); ++g) {
                    if (post[g] > post[best]) {
                        best = g;
                    }
                }
                if (best == (size_t)e.called_i) {
                    continue;
                }
                // Rendered out of compact traversal space here. `best` indexes the entry's own
                // allele list, which is *not* the record's ALT list: writing it straight through --
                // as this block did -- put allele numbers like 5 on a record carrying one ALT, which
                // is not a parseable VCF. It was the one place stage 1 left the two numberings
                // touching, and both the genotype patch and the phase patch below read it.
                const int r_called = vcf_allele_of(allele_arena, e.allele_offset, e.num_alleles,
                                                   e.called_i);
                const int r_best = vcf_allele_of(allele_arena, e.allele_offset, e.num_alleles, best);
                // Recorded whether or not it can be written: it is what the model chose, and the
                // mosaic and every child's strand are built from it, not from the VCF.
                Regenotyped rg;
                rg.compact = best;
                rg.vcf_allele = r_best;
                rg.traversal = traversal_of(trav_arena, e.trav_offset, e.num_alleles, best);
                rg.called_allele = r_called;
                nested_regenotyped[e.record_key] = rg;
                entries[idxs[k]].final_i = (uint16_t)best;
                entries[idxs[k]].final_j = (uint16_t)best;
                if (!e.emitted) {
                    continue;   // it moved, but there is no line to rewrite
                }
                if (r_called < 0 || r_best < 0) {
                    ++unrenderable;
                    continue;
                }
                Change c;
                c.record_key = e.record_key;
                c.contig = contig_names[e.contig];
                c.position = e.position;
                c.called_i = (size_t)r_called;
                c.called_j = (size_t)r_called;
                c.allele_i = (size_t)r_best;
                c.allele_j = (size_t)r_best;
                c.posterior = post[best];
                c.explained_share = (double)e.explained_share;
                changes.push_back(c);
            }
        }

        // The PhaseCalls for these sites were built before step three ran, so they still describe
        // the pre-linkage allele. apply_phasing refuses a phase that disagrees with the genotype
        // the record ends up carrying -- correctly, since phasing the wrong allele pair would be
        // worse than leaving it unphased -- so a stale PhaseCall silently costs the site its PS.
        // That is what turned 64 unphased records into 466: 402 nested sites whose genotype step
        // three had moved out from under their own phase.
        //
        // Patched here rather than by reordering the passes, because strand assignment needs the
        // parent's phase while genotype resolution needs the strand, so neither can simply go first.
        // How the nested population fared against the parent genotypes linkage actually chose.
        //
        // This used to be a coarser count -- nested sites whose parent's genotype moved at all --
        // which was an upper bound gathered to decide whether a two-pass rebuild was worth its
        // cost. It is superseded: the crossing mask says exactly which sites were placed on the
        // wrong strand and which have the wrong ploidy, so the bound is no longer the best
        // available number.
        if (unrenderable > 0) {
#pragma omp critical (cerr)
            std::cerr << "[vg call] linkage: " << unrenderable
                      << " settled genotypes name a traversal the record carries no ALT for, so the"
                      << " genotype was left as called" << std::endl;
        }
        if (phase_fallback > 0) {
#pragma omp critical (cerr)
            std::cerr << "[vg call] phasing: " << phase_fallback
                      << " records are phased on the alleles the line carries rather than the ones"
                      << " the model settled on, for the same reason" << std::endl;
        }
        if (!deferred_nested.empty()) {
#pragma omp critical (cerr)
            std::cerr << "[vg call] nested strands: " << deferred_nested.size()
                      << " sites, "
                      << (deferred_nested.size() - carried_on_both - unplaced_no_strand
                          - unplaced)
                      << " placed on one strand; " << carried_on_both
                      << " carried on both parent strands (" << carried_on_both_emitted
                      << " with a line), " << unplaced_no_strand << " on neither ("
                      << unplaced_no_strand_emitted << " with a line), " << unplaced
                      << " with no phased parent -- the last two get no strand" << std::endl;
        }

        if (!nested_regenotyped.empty()) {
            for (PhaseCall& pc : *phasing_out) {
                auto found = nested_regenotyped.find(pc.record_key);
                if (found != nested_regenotyped.end()) {
                    // The traversal always; the VCF allele only where the record has one. Writing
                    // the compact index here is what put allele numbers past the end of the ALT list
                    // into phased GTs, and the strand carried the same wrong number with it.
                    pc.trav_first = found->second.traversal;
                    pc.trav_second = found->second.traversal;
                    // The settled allele where the record carries it, else the one the line already
                    // has: a phase that names an allele the record lacks is declined outright, and
                    // the record loses its strand as well as its genotype.
                    const int v = found->second.vcf_allele >= 0
                                      ? found->second.vcf_allele
                                      : found->second.called_allele;
                    if (found->second.vcf_allele < 0 && v >= 0) {
                        ++phase_fallback;
                    }
                    pc.allele_first = v >= 0 ? (size_t)v : LinkageModel::WILDCARD;
                    pc.allele_second = pc.allele_first;
                }
            }
        }
    }


    // Reference order, and only now.
    //
    // The nested sites are appended after every chain has been phased, because placing one needs its
    // parent's phase and the parent can be anywhere in any chain. That leaves `phasing_out` in two
    // runs -- chain sites in order, then nested sites in order -- and write_mosaic reads it as one
    // ordered sweep, closing a segment only where the haplotype changes. Out of order, a nested site
    // 451 kb into chr20 landed in the same run as one 65 Mb along: chr20's mosaic carried five
    // segments spanning tens of megabases, one of them claiming 284 sites between ref_start 451,374
    // and ref_end 65,512,343, with start_node and end_node anchors to match. The site totals still
    // added up, which is what the harness checks, so the file looked complete while being wrong
    // about where those sites were.
    //
    // Sorted here rather than in write_mosaic so every consumer gets one guarantee instead of each
    // having to know. The key matches the chain sort: position, then the site's own key, so two
    // records at one position keep an order that is a property of the sites and not of the threads.
    if (phasing_out != nullptr && last) {
        std::sort(phasing_out->begin(), phasing_out->end(),
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
    return changes;
}

}
