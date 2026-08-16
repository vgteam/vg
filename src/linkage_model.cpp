#include "linkage_model.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

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

    // Pin: everything except the pinned state becomes unreachable at that index.
    if (pin_index != (size_t)-1 && pin_index >= from && pin_index < to) {
        size_t t = pin_index - from;
        size_t pa = pin.first == WILDCARD ? n_hap : min(pin.first, n_hap);
        size_t pb = pin.second == WILDCARD ? n_hap : min(pin.second, n_hap);
        double keep = emissions[t][pa * m + pb];
        emissions[t].assign(m * m, 0.0);
        // If the pinned state is itself impossible here the pin would empty the window, so fall
        // back to leaving the site free rather than returning nothing.
        emissions[t][pa * m + pb] = keep > 0.0 ? keep : 1.0;
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
    if (pin_index != (size_t)-1 && pin_index >= from && pin_index < to) {
        size_t t = pin_index - from;
        size_t pa = pin == WILDCARD ? n_hap : min(pin, n_hap);
        double keep = emissions[t][pa];
        emissions[t].assign(m, 0.0);
        emissions[t][pa] = keep > 0.0 ? keep : 1.0;
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

void LinkageCollector::record(const string& contig, size_t position, size_t num_alleles,
                              const vector<double>& genotype_ln_likelihood,
                              const vector<int>& haplotype_allele,
                              size_t called_i, size_t called_j, size_t record_key,
                              double explained_share, size_t ploidy,
                              int64_t start_node, int64_t end_node) {
    if (num_alleles == 0 || genotype_ln_likelihood.empty()) {
        return;
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
    e.num_alleles = (uint16_t)num_alleles;
    e.called_i = (uint16_t)called_i;
    e.called_j = (uint16_t)called_j;
    e.record_key = record_key;
    e.explained_share = (float)explained_share;
    e.start_node = start_node;
    e.end_node = end_node;
    e.ploidy = (uint8_t)(ploidy == 1 ? 1 : 2);

    e.gl_offset = (uint32_t)gl_arena.size();
    for (double v : genotype_ln_likelihood) {
        // float, not double: these are log-likelihood differences fed to an exp(), and the
        // ratios that survive are nowhere near float's precision limit. It halves the arena.
        gl_arena.push_back((float)v);
    }
    e.hap_offset = (uint32_t)hap_arena.size();
    for (size_t h = 0; h < n_haplotypes; ++h) {
        int allele = h < haplotype_allele.size() ? haplotype_allele[h] : -1;
        // int8 caps alleles per site at 127, which no snarl this caller emits approaches; a site
        // that did would lose linkage rather than be mis-linked, since -1 means "absent".
        hap_arena.push_back(allele >= 0 && allele < 127 && (size_t)allele < num_alleles
                                ? (int8_t)allele : (int8_t)-1);
    }
    entries.push_back(e);
}

size_t LinkageCollector::bytes() const {
    return entries.size() * sizeof(Entry)
           + gl_arena.size() * sizeof(float)
           + hap_arena.size() * sizeof(int8_t);
}

vector<LinkageCollector::Change> LinkageCollector::resolve(vector<PhaseCall>* phasing_out) const {
    vector<Change> changes;
    if (!model.active() || entries.empty()) {
        return changes;
    }

    // Group by contig, then sort by reference position. Node-ID order is close to reference order
    // in a reference-first graph but is not guaranteed to be it, and the transition probabilities
    // are computed from the gaps -- so trusting the arrival order would silently feed the model
    // the wrong distances.
    vector<vector<size_t>> by_contig(contig_names.size());
    for (size_t i = 0; i < entries.size(); ++i) {
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
            run_start = run_end;
        }
    }

    for (auto& indices : chains) {
        if (indices.size() < 2) {
            continue;
        }

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
            if (best == before) {
                continue;
            }
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
            Change c;
            c.record_key = e.record_key;
            c.contig = contig_names[e.contig];
            c.position = e.position;
            c.called_i = e.called_i;
            c.called_j = e.called_j;
            c.allele_i = i;
            c.allele_j = j;
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
                pc.allele_first = i;
                pc.allele_second = j;
            }
            phasing_out->push_back(pc);
        }
    }
    return changes;
}

}
