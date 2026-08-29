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

/// Per-site pins offered to window_phasing, and those it refused because the pinned haplotype pair
/// cannot spell the genotype the site is constrained to. A refused pin on a group's PARENT frees the
/// whole group's orientation while its haploid siblings stay tied to the parent's phase.
static std::atomic<size_t> g_pin_applied(0), g_pin_declined(0);
/// V1: mean -log P(neither strand switches) per step, in units of 1e-9, and the step count.
static std::atomic<size_t> g_margin_neglog(0), g_margin_steps(0);
/// Groups whose parent was never offered a pin at all -- no PhaseCall to pin it to.
static std::atomic<size_t> g_group_parent_unpinned(0), g_group_parent_pinned(0);

/// The distance between two adjacent sites in a chain, in bp, PER STRAND.
///
/// Prefers the explicit `gap_to_previous` the caller measured along the haplotype's own traversal, and
/// falls back to the reference position difference where there is none. Clamped at 1, as the position
/// form was: `switch_probability` reads a gap of 0 as 1 anyway, and a negative one would wrap.
///
/// One slot set and the other not is not an error and not a reference fallback: it is a step measured
/// along one haplotype and unmeasurable along the other -- a parent traversal that does not enter the
/// chain, or a cross-parent step. Using the known distance for both strands is what a single value
/// did, and it is a better estimate than the reference difference wherever the frame exists at all.
static inline std::pair<size_t, size_t> site_gap(const LinkageModel::Site& prev,
                                                const LinkageModel::Site& next) {
    const int64_t a = next.gap_to_previous[0], b = next.gap_to_previous[1];
    auto clamp = [](int64_t g) { return (size_t)(g > 1 ? g : 1); };
    if (a >= 0 && b >= 0) {
        return {clamp(a), clamp(b)};
    }
    if (a >= 0) {
        return {clamp(a), clamp(a)};
    }
    if (b >= 0) {
        return {clamp(b), clamp(b)};
    }
    if (prev.unpositioned || next.unpositioned) {
        // EITHER, not both. A mixed pair is the common case -- a parent's children interleave chains
        // the reference crosses with chains it does not -- and it falls through to the arithmetic
        // below with a 0 standing in for a coordinate: a positioned prev gives `0 > pos` false and
        // so a gap of 1, which is the strongest possible link, and an unpositioned prev gives the
        // next site's whole offset measured from nothing. One coordinate does not make a distance.
        //
        // The clamp below is the wrong answer for that, not a neutral one: 1 bp is
        // `switch_probability`'s minimum rate, which is the strongest possible link -- so two sites
        // about which nothing is known would be asserted to be perfectly linked. SIZE_MAX gives
        // rho = 1.0, at which `transition_apply`'s T = (1-rho)I + (rho/m)11' is exactly uniform:
        // the chain forgets, which is what "unknown" means here.
        return {numeric_limits<size_t>::max(), numeric_limits<size_t>::max()};
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
/// turn on by default, and it is checked again where the record is phased, which refuses any phased GT that
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
void transition_apply(const vector<double>& in, size_t m,
                      double rho_a, double rho_b, vector<double>& out) {
    // One switch probability per strand. The two haplotypes of a diploid sample recombine
    // independently, so the distance each has travelled since the previous site is its own -- which
    // is what the haplotype-frame work needs and what a single scalar cannot express.
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
            // The old grouping `stay * jump * (row[a] + col[b])` cannot survive: once the two
            // coefficients differ the row and column terms carry different factors and must be
            // written separately. That re-association is why this cannot be byte-identical even
            // when both rhos are equal.
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
    double S_a = stay_a > 0.0 ? log(stay_a) : NEG_INF;
    double S_b = stay_b > 0.0 ? log(stay_b) : NEG_INF;
    double J_a = jump_a > 0.0 ? log(jump_a) : NEG_INF;
    double J_b = jump_b > 0.0 ? log(jump_b) : NEG_INF;

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
                double v = c1 + S_a + S_b;
                if (v > best) { best = v; ba = ap; bb = bp; }
            }
            double c2 = rowExcl[ap * m + bp];
            if (c2 > NEG_INF) {
                double v = c2 + S_a + J_b;
                if (v > best) { best = v; ba = ap; bb = rowExclArg[ap * m + bp]; }
            }
            bool chit = (cols[bp].arg == ap);
            double c3 = chit ? cols[bp].second : cols[bp].best;
            if (c3 > NEG_INF) {
                double v = c3 + J_a + S_b;
                if (v > best) { best = v; ba = chit ? cols[bp].arg2 : cols[bp].arg; bb = bp; }
            }
            bool bhit = (both.arg == ap);
            double c4 = bhit ? both.second : both.best;
            if (c4 > NEG_INF) {
                double v = c4 + J_a + J_b;
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
                                     const vector<double>* beta_in,
                                     vector<vector<double>>* alpha_out,
                                     vector<vector<double>>* beta_out,
                                     const vector<char>* harvest_mask) const {
    auto wanted = [&](size_t i) {
        return harvest_mask == nullptr || (i < harvest_mask->size() && (*harvest_mask)[i]);
    };
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
        if (alpha_out != nullptr) {
            alpha_out->assign(sites.size(), vector<double>());
            if (wanted(from)) {
                (*alpha_out)[from] = alpha_in != nullptr && alpha_in->size() == m * m
                                         ? *alpha_in
                                         : vector<double>(m * m, 1.0 / (double)(m * m));
            }
        }
    }
    for (size_t t = 1; t < n; ++t) {
        // An explicit per-step distance where the caller supplied one -- below the top level that is
        // measured along the parent's settled traversal, not along the reference. Falls back to the
        // reference position difference where it is unset, which is the top-level case and the
        // cross-parent case that still needs the parent's reference span.
        const std::pair<size_t, size_t> gap = site_gap(sites[from + t - 1], sites[from + t]);
        vector<double> moved;
        // One switch probability per strand. Equal wherever the two distances are, which is every
        // step that falls back to the reference difference or knows only one frame.
        const double rho_a = switch_probability(gap.first);
        const double rho_b = switch_probability(gap.second);
        transition_apply(alpha[t - 1], m, rho_a, rho_b, moved);
        if (alpha_out != nullptr && wanted(from + t)) {
            // The message ENTERING t, before its own emission -- what `alpha_in` consumes, and
            // symmetric with `beta_out`, which is already the entering message from the right.
            (*alpha_out)[from + t] = moved;
        }
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
    if (beta_out != nullptr) {
        beta_out->assign(sites.size(), vector<double>());
    }
    vector<double> beta(m * m, 1.0);
    if (beta_in != nullptr && beta_in->size() == m * m) {
        beta = *beta_in;
    }
    for (size_t t = n; t-- > 0;) {
        if (beta_out != nullptr && wanted(from + t)) {
            (*beta_out)[from + t] = beta;
        }
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
        // An explicit per-step distance where the caller supplied one -- below the top level that is
        // measured along the parent's settled traversal, not along the reference. Falls back to the
        // reference position difference where it is unset, which is the top-level case and the
        // cross-parent case that still needs the parent's reference span.
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
/// Correct for a *per-site* quantity and only for that: a marginal posterior at one site does not
/// depend on where a window edge fell, so windows can be decoded independently and their middles
/// kept. A path cannot be done this way -- see `windowed_path`.
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

/// Overlapping windows chained by a seam pin.
///
/// Windowing needs more care here than it does for the marginals. Decode two windows
/// independently and the state at the seam is chosen twice, by two runs that never saw each other,
/// so the join manufactures a switch at every window boundary -- an artifact that would read as a
/// biological result and would scale with the site count rather than with the genome.
///
/// So each window after the first is pinned: the state at one index inside its leading margin is
/// fixed to what the previous window already decided there. The pin sits in the margin rather than
/// at the boundary precisely because a margin-deep state was decided with the full window's
/// context around it, while a boundary state was not.
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

vector<vector<double>> LinkageModel::posteriors_with_context(
        const vector<Site>& sites, const vector<char>& harvest_mask,
        vector<vector<double>>& context_out) const {
    vector<vector<double>> out(sites.size());
    context_out.assign(sites.size(), vector<double>());
    if (sites.empty()) {
        return out;
    }
    // The same windowed sweep `posteriors` does, so the posteriors it returns are unchanged; the
    // context is harvested from each window's KEPT range only, for the same reason the posteriors
    // are -- a value from a window's margin was decoded without the full context around it.
    const size_t step = max<size_t>(params.window, 1);
    const size_t margin = params.margin;
    const size_t n = sites.size();
    if (n <= step) {
        vector<vector<double>> a, b;
        window_posteriors(sites, 0, n, out, nullptr, nullptr, &a, &b, &harvest_mask);
        for (size_t i = 0; i < n; ++i) {
            if (i < a.size() && !a[i].empty() && i < b.size() && !b[i].empty()) {
                context_out[i] = a[i];
                for (size_t k = 0; k < context_out[i].size() && k < b[i].size(); ++k) {
                    context_out[i][k] *= b[i][k];
                }
            }
        }
        return out;
    }
    for (size_t start = 0; start < n; start += step) {
        const size_t lo = start > margin ? start - margin : 0;
        const size_t hi = min(start + step + margin, n);
        vector<vector<double>> local(n), a, b;
        window_posteriors(sites, lo, hi, local, nullptr, nullptr, &a, &b, &harvest_mask);
        const size_t keep_to = min(start + step, n);
        for (size_t t = start; t < keep_to; ++t) {
            out[t] = std::move(local[t]);
            if (t < a.size() && !a[t].empty() && t < b.size() && !b[t].empty()) {
                context_out[t] = a[t];
                for (size_t k = 0; k < context_out[t].size() && k < b[t].size(); ++k) {
                    context_out[t][k] *= b[t][k];
                }
                double s = 0.0;
                for (double v : context_out[t]) {
                    s += v;
                }
                if (s > 0.0) {
                    for (double& v : context_out[t]) {
                        v /= s;
                    }
                }
            }
        }
    }
    return out;
}

vector<vector<double>> LinkageModel::posteriors(const vector<Site>& sites) const {
    vector<vector<double>> out(sites.size());
    if (sites.empty()) {
        return out;
    }
    windowed_marginals(sites.size(), max<size_t>(params.window, 1), params.margin, out,
                       [&](size_t lo, size_t hi, vector<vector<double>>& local) {
                           window_posteriors(sites, lo, hi, local);
                       });
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
            //
            // COUNTED, because "not expected" was never measured. In a per-parent group the parent
            // is the ONLY pinned member, so declining its pin leaves the group's Viterbi orientation
            // free -- while the group's haploid siblings are placed by parent-traversal identity
            // against the parent's own phase. Two frames with nothing relating them, which is the
            // incoherence the strand bookkeeping cannot see because each half is self-consistent.
            ++g_pin_declined;
            continue;
        }
        ++g_pin_applied;
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
        // An explicit per-step distance where the caller supplied one -- below the top level that is
        // measured along the parent's settled traversal, not along the reference. Falls back to the
        // reference position difference where it is unset, which is the top-level case and the
        // cross-parent case that still needs the parent's reference span.
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
        // An explicit per-step distance where the caller supplied one -- below the top level that is
        // measured along the parent's settled traversal, not along the reference. Falls back to the
        // reference position difference where it is unset, which is the top-level case and the
        // cross-parent case that still needs the parent's reference span.
        // One strand, so one distance: `site_gap` hands back the same value in both slots wherever
        // only one is known, which is every haploid step.
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
        // An explicit per-step distance where the caller supplied one -- below the top level that is
        // measured along the parent's settled traversal, not along the reference. Falls back to the
        // reference position difference where it is unset, which is the top-level case and the
        // cross-parent case that still needs the parent's reference span.
        // One strand, so one distance: `site_gap` hands back the same value in both slots wherever
        // only one is known, which is every haploid step.
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

vector<vector<double>> LinkageModel::haploid_posteriors(const vector<Site>& sites) const {
    vector<vector<double>> out(sites.size());
    if (sites.empty()) {
        return out;
    }
    windowed_marginals(sites.size(), max<size_t>(params.window, 1), params.margin, out,
                       [&](size_t lo, size_t hi, vector<vector<double>>& local) {
                           window_haploid_posteriors(sites, lo, hi, local);
                       });
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
        // An explicit per-step distance where the caller supplied one -- below the top level that is
        // measured along the parent's settled traversal, not along the reference. Falls back to the
        // reference position difference where it is unset, which is the top-level case and the
        // cross-parent case that still needs the parent's reference span.
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

vector<size_t> LinkageModel::haploid_phasing(const vector<Site>& sites,
                                             const vector<size_t>& constraint) const {
    vector<size_t> out(sites.size(), WILDCARD);
    if (sites.empty()) {
        return out;
    }
    // The unpinned value is WILDCARD here, not a default-constructed T, so the first window is
    // given the same "no pin" the diploid pass gives it -- window_haploid_phasing reads the index,
    // not the value, to decide whether there is one.
    windowed_path(sites.size(), max<size_t>(params.window, 1), params.margin, out,
                  [&](size_t lo, size_t hi, size_t pin_index, const size_t& pin,
                      vector<size_t>& local) {
                      window_haploid_phasing(sites, lo, hi, constraint, pin_index, pin, local);
                  });
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
                              bool emitted, int align_rank, int chain_index,
                              bool chain_backward, bool unpositioned) {
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
    // Passed in, not set afterwards. The barrier computes a chain's alignment rank before the block
    // that can create this entry, so a `set_align_rank` call for a chain the sweep never recorded
    // finds no entry and is discarded -- silently, since its return was ignored. That left the 520
    // chains on chr20 that are reachable only under the settled parent with no rank at all, which is
    // the population the ordering exists for. `parent_trav` is an argument for the same reason.
    e.align_rank = (int32_t)align_rank;
    // A pure graph fact, so it is known at descent and needs none of the barrier's plumbing.
    e.chain_index = (int32_t)chain_index;
    e.chain_backward = chain_backward;
    e.unpositioned = unpositioned;
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
    const uint32_t at = (uint32_t)entries.size();
    auto last = last_by_key.find(record_key);
    if (last == last_by_key.end()) {
        first_by_key[record_key] = at;
    } else {
        entries[last->second].next_same_key = at;
    }
    last_by_key[record_key] = at;
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
    // Every arena, and the key index beside them. The index is two hash nodes a site, which is
    // real -- 21 MB on chr20 -- and leaving it out would make the reported figure drift from the
    // measured one by more than half the total it prints.
    return entries.size() * sizeof(Entry)
           + gl_arena.size() * sizeof(float)
           + hap_arena.size() * sizeof(int8_t)
           + trav_arena.size() * sizeof(uint16_t)
           + allele_arena.size() * sizeof(int8_t)
           + (first_by_key.size() + last_by_key.size())
                 * (sizeof(std::pair<const size_t, uint32_t>) + sizeof(void*));
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
    // exactly as a freshly recorded one would be -- now by construction rather than by two copies
    // agreeing.
    const CompactSite cs = compact_site(genotype_ln_likelihood, haplotype_traversal,
                                        called_trav_i, called_trav_j, ploidy);
    if (!cs.ok) {
        return false;
    }
    const vector<int>& space = cs.space;
    const vector<float>& gls = cs.gls;
    const size_t k = space.size();
    const int ci = cs.ci, cj = cs.cj;
    const size_t site_ploidy = cs.site_ploidy;
    auto compact_of = [&](int trav) { return cs.compact_of(trav); };

    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found != NO_ENTRY) {
        Entry& e = entries[found];
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
        e.emitted = emitted;
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

bool LinkageCollector::set_frame(size_t record_key, int slot, int offset) {
    if (slot < 0 || slot > 1) {
        return false;
    }
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found != NO_ENTRY) {
        Entry& e = entries[found];
        e.frame_offset[slot] = (int32_t)offset;
        return true;
    }
    return false;
}

bool LinkageCollector::set_align_rank(size_t record_key, int rank, bool chain_backward) {
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found == NO_ENTRY) {
        return false;
    }
    entries[found].align_rank = (int32_t)rank;
    entries[found].chain_backward = chain_backward;
    return true;
}

bool LinkageCollector::set_parent_trav(size_t record_key, int parent_trav) {
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found == NO_ENTRY) {
        return false;
    }
    entries[found].parent_trav = (int16_t)parent_trav;
    return true;
}

bool LinkageCollector::has_entry(size_t record_key) const {
    lock_guard<std::mutex> guard(mutex);
    return live_index(record_key) != NO_ENTRY;
}

bool LinkageCollector::retract(size_t record_key) {
    lock_guard<std::mutex> guard(mutex);
    const uint32_t found = live_index(record_key);
    if (found == NO_ENTRY) {
        return false;
    }
    entries[found].retracted = true;
    return true;
}

/// Split the settled compact pair into the two things that consume it: the traversal on each
/// strand, which is the genome fact a crossing mask and a haplotype path are expressed in, and the
/// VCF allele on each strand, which is the only form the renderer can write.
///
/// Leaving `allele_*` compact made the phased-GT guard reject the pair and silently drop the
/// record's phasing -- 1,528 extra strandless records on chr20. It was spelled out three times,
/// in the diploid chain, the nested-strand pass and the haploid pass, with that warning copied
/// alongside each.
///
/// The three copies had drifted: only the diploid one also widened `order_arbitrary` on a
/// fallback. The other two are ploidy-1 populations, where the pair is one allele twice and the
/// test cannot fire, so folding it in changes nothing -- asserted by byte-identity rather than
/// argued.
void LinkageCollector::finish_phase_call(PhaseCall& pc, const Entry& e, size_t& fallbacks) const {
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
        ++fallbacks;
        pc.order_arbitrary = pc.order_arbitrary || (v_first != v_second);
    }
}

static std::atomic<size_t> g_grp_no_parent(0), g_grp_no_context(0);
static std::atomic<size_t> g_grp_no_entry(0), g_grp_ploidy(0), g_grp_vetoed(0);

size_t LinkageCollector::resolve_generation(
        size_t generation, bool last, vector<PhaseCall>* phasing_out) {
    // Nested sites, collected across every contig and phased after the diploid chains, so their
    // strand can be read off the parent they hang off rather than guessed.
    vector<size_t> deferred_nested;
    size_t moved = 0;
    // Genotypes the model settled on that the emitted record carries no ALT for. Reachable only
    // because the collector's allele space is the genotyper's rather than the VCF's, so it is
    // counted and reported rather than assumed away.
    // Records whose phase names the called pair because the settled pair had no ALT. The phase is
    // still real -- the block and the strand order come from the panel either way -- but the alleles
    // it names are the line's rather than the model's, so the count belongs in the report.
    size_t phase_fallback = 0;
    if (!model.active() || entries.empty()) {
        return moved;
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
    // Which sites will be asked for a context message: the parents of everything not yet called.
    // Built once, because a child's parent may sit in any chain on the contig.
    unordered_set<size_t> wanted_parents;
    for (const Entry& e : entries) {
        if (!e.retracted && e.generation > generation && e.parent_record_key != 0) {
            wanted_parents.insert(e.parent_record_key);
        }
    }

    // Per chain: the message to condition it on (empty = none), and the phase set it belongs to
    // (SIZE_MAX = take it from the chain's own first site, as a contig chain does).
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
            chain_context.push_back(nullptr);
            chain_phase_set.push_back(numeric_limits<size_t>::max());
            run_start = run_end;
        }
    }

    // A generation's own sites hang off a bounded set of parents, and a parent's context message is
    // everything the contig chain knows about it. Conditioning each parent's children on that
    // message and decoding them alone is the tree factorisation: nothing is inserted into anything,
    // so nothing downstream goes stale, and the cost is the number of children rather than the
    // length of the contig.
    //
    // Only taken when EVERY live site in the generation is covered -- a parent with no stored
    // context would otherwise have its children decoded with no context at all, which is the
    // uniform boundary that discards everything.
    size_t grouped_sites = 0, grouped_groups = 0;
    // Stage B, per generation: how the diploid nested steps were measured.
    size_t fg_pair = 0, fg_differ = 0, fg_single = 0, fg_one_slot = 0, fg_no_frame = 0;
    size_t fg_cross_parent = 0, fg_mirrored = 0, fg_equals_ref = 0, fg_parent_step = 0;
    size_t fg_one_both_claimed = 0, fg_one_single_claimed = 0, fg_one_undetermined = 0;
    size_t unpos_diploid_sites = 0, unpos_steps = 0, unpos_framed = 0, unpos_uniform = 0;
    // Stage C, per generation: how often the alignment order and the reference order disagree.
    size_t ao_pairs = 0, ao_reordered = 0, ao_unranked = 0, ao_group_unranked = 0;
    // VG_LINKAGE_NO_GROUPING forces the contig-chain path, so one binary can produce both arms of
    // a comparison and the only difference between them is this.
    static const bool no_grouping = getenv("VG_LINKAGE_NO_GROUPING") != nullptr;
    // How a nested sibling's distance from the previous one is measured. One binary, three arms:
    //
    //   0  reference position difference -- what every diploid chain has always used
    //   1  the parent's settled traversals, collapsed to ONE scalar (their mean)
    //   2  the parent's settled traversals, one distance PER STRAND
    //
    // 1 exists to separate the two things stage B changes at once. Spacing along a traversal at all
    // has been measured before -- twice, in the haploid derivation -- and cost about 0.0005 of
    // JointIndel on chr20 as a single scalar. Arm 1 reproduces that form here, so if 2 loses it is
    // possible to say whether the loss is traversal spacing or the per-strand split, which a single
    // on/off switch could not.
    static const int frame_gap_mode = []() {
        const char* v = getenv("VG_LINKAGE_FRAME_GAPS");
        return v != nullptr ? atoi(v) : 0;
    }();
    // Order a parent's children by the alignment of its two settled traversals instead of by
    // reference position. On by default: where both orders exist they are predicted to agree, and
    // where only the alignment exists it is the only order there is.
    //
    // VG_LINKAGE_NO_ALIGN_ORDER forces reference order, so one binary produces both arms.
    static const bool no_align_order = getenv("VG_LINKAGE_NO_ALIGN_ORDER") != nullptr;
    // VG_LINKAGE_PER_CHAIN_GROUPS: decode each nested CHAIN on its own, conditioned on the parent,
    // instead of decoding all of a parent's children together.
    //
    // This is the whole of the "is linkage BETWEEN nested chains worth anything?" question. Within a
    // chain the snarls stay linked, so recombination inside a chain is still modelled; between two
    // chains under one parent there is no transition at all, and each is answered from the parent's
    // context alone. `align_rank` is the chain identity -- every snarl of one chain collapses to one
    // symbol and so shares a column -- so it is exactly the discriminator this needs.
    //
    // If the two arms score the same, inter-chain ordering and inter-chain distance are both dead
    // weight: order and spacing only matter where a transition crosses them.
    static const bool per_chain_groups = getenv("VG_LINKAGE_PER_CHAIN_GROUPS") != nullptr;
    if (generation > 0 && !parent_context.empty() && !no_grouping) {
        unordered_map<size_t, size_t> index_of_key;
        index_of_key.reserve(entries.size() * 2);
        for (size_t i = 0; i < entries.size(); ++i) {
            if (!entries[i].retracted) {
                index_of_key[entries[i].record_key] = i;
            }
        }
        // Keyed by (parent, chain-within-parent). The chain half is 0 unless per_chain_groups, so
        // the default is exactly the previous one-group-per-parent behaviour.
        map<pair<size_t, size_t>, vector<size_t>> by_parent;
        auto group_key = [&](const Entry& e) {
            if (!per_chain_groups) {
                return make_pair(e.parent_record_key, (size_t)0);
            }
            // An unranked chain cannot be pooled with anything, so it becomes its own group rather
            // than being lumped with every other unranked chain under this parent.
            const size_t chain = e.align_rank >= 0 ? (size_t)e.align_rank
                                                   : numeric_limits<size_t>::max() - e.record_key;
            return make_pair(e.parent_record_key, chain);
        };
        // PER CHAIN, not per generation.
        //
        // This was one global flag: a single live site whose parent could not be found sent the
        // WHOLE generation back to the contig chain. It had already been caught once that way, when
        // 51 sites vetoed 21,447, and the off-reference population trips it again -- grouping engaged
        // at generations 1 and 4 on chr20 and not at 2 or 3, and the cost is paid by every other site
        // in the generation. On the contig chain adjacent sites have different parents, so their step
        // has no common frame; and where such a step touches a site with no reference position there
        // is nothing to fall back to, so it goes UNIFORM and severs the chain. That is about 1,100
        // severed links a generation, against 156 where grouping did engage.
        //
        // A chain is an independent decode, so the veto belongs at that scope. A chain whose live
        // sites can all be grouped is replaced by its groups; one that cannot is kept exactly as it
        // was. Nothing is lost either way, which is what the global flag was protecting: a site left
        // out of every group would never be decoded, never settled and never phased.
        vector<vector<size_t>> kept_chains;
        size_t vetoed_chains = 0, grouped_chains = 0;
        for (const vector<size_t>& indices : chains) {
            bool chain_covered = true;
            for (size_t idx : indices) {
                const Entry& e = entries[idx];
                if (e.generation != generation) {
                    continue;
                }
                auto par = index_of_key.find(e.parent_record_key);
                if (e.parent_record_key == 0) {
                    ++g_grp_no_parent;
                    chain_covered = false;
                } else if (par == index_of_key.end()) {
                    ++g_grp_no_entry;
                    chain_covered = false;
                } else if (entries[par->second].ploidy != e.ploidy) {
                    ++g_grp_ploidy;
                    chain_covered = false;
                }
            }
            if (!chain_covered) {
                kept_chains.push_back(indices);
                ++vetoed_chains;
                ++g_grp_vetoed;
                continue;
            }
            ++grouped_chains;
            for (size_t idx : indices) {
                const Entry& e = entries[idx];
                if (e.generation != generation) {
                    continue;
                }
                auto ctx = parent_context.find(e.parent_record_key);
                {
                    // A parent with no stored context still forms a group. Its own emission is in
                    // the group, so its children link to it and to each other; what is missing is
                    // the contig's upstream context, not the parent. Counted, because it is a
                    // weaker decode than the rest and the population should not be invisible.
                    //
                    // These arise from children recorded at the barrier AFTER their generation's
                    // decode -- the chains reachable only under a settled parent -- whose parents
                    // were therefore never masked for harvesting.
                    if (ctx == parent_context.end()) {
                        ++g_grp_no_context;
                    }
                    by_parent[group_key(e)].push_back(idx);
                    continue;
                }
            }
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
                const size_t pidx = index_of_key[kv.first.first];
                vector<size_t> group;
                group.push_back(pidx);
                // Stage C: the alignment of the parent's two settled traversals is the order,
                // where both children have a rank in it. Reference position is the fallback, and
                // the two are compared rather than one silently replacing the other -- a change
                // that reorders nothing is worth being able to say so about.
                for (size_t a = 1; a < kv.second.size(); ++a) {
                    const Entry& p1 = entries[kv.second[a - 1]];
                    const Entry& p2 = entries[kv.second[a]];
                    if (p1.align_rank >= 0 && p2.align_rank >= 0) {
                        ++ao_pairs;
                        if ((p1.align_rank < p2.align_rank) != (p1.position < p2.position)
                            && p1.position != p2.position && p1.align_rank != p2.align_rank) {
                            ++ao_reordered;
                        }
                    } else {
                        ++ao_unranked;
                    }
                }
                // A TOTAL key per entry, not a chain of skippable comparisons.
                //
                // "Compare on this key only if BOTH operands have it, else fall through" is not a
                // strict weak ordering, and std::sort on one is undefined behaviour rather than
                // merely a different order. With absence spelled as -1 and one absent operand
                // skipping the key, A(rank -, pos 2), B(rank 1, pos 3), C(rank 2, pos 1) gives
                // A<B, B<C and C<A: a brute force over small triples finds 204 such cycles here and
                // 81 in the frame-and-position form that predates the alignment rank, so this was
                // latent before the rank was added and the rank widened it.
                //
                // An absent rank therefore sorts LAST, uniformly, rather than deferring to position.
                // It is a convention and not a claim -- a chain the alignment could not place has no
                // defined position among the ones it could -- and its virtue is being consistent.
                // Sorting absence FIRST is the alternative and is known bad: an entry at the head of
                // its group hands its neighbour a spurious multi-megabase gap, which is the shape of
                // the 448-site regression recorded at the strand derivation.
                // The choice of key is made ONCE for the group, not per pair. Every child here
                // hangs off the same parent, so either the alignment placed all of them and their
                // ranks are mutually comparable, or it did not and reference position is the only
                // order that covers the whole group.
                //
                // Deciding per pair is what breaks: an absent rank cannot be compared, so skipping
                // the key for that pair alone mixes two orders inside one group and is intransitive
                // -- 204 cycles over small triples, against 81 for the frame-and-position form that
                // predates the rank. Substituting a sentinel is transitive but displaces the
                // unplaced chain to one end of its group, which distorts its neighbours' gaps; that
                // is the shape of the 448-site regression recorded at the strand derivation, and it
                // moved 241 chr20 lines when tried here.
                //
                // All-or-nothing per group is transitive, moves nobody, and degrades to exactly the
                // previous behaviour wherever the alignment is incomplete.
                bool all_ranked = !no_align_order;
                for (size_t idx : kv.second) {
                    if (entries[idx].align_rank < 0) {
                        all_ranked = false;
                        break;
                    }
                }
                if (!all_ranked) {
                    ++ao_group_unranked;
                }
                bool all_indexed = true;
                for (size_t idx : kv.second) {
                    if (entries[idx].chain_index < 0) {
                        all_indexed = false;
                        break;
                    }
                }
                sort(kv.second.begin(), kv.second.end(), [&](size_t a, size_t c) {
                    const Entry& ea = entries[a];
                    const Entry& ec = entries[c];
                    if (all_ranked && ea.align_rank != ec.align_rank) {
                        return ea.align_rank < ec.align_rank;
                    }
                    // Equal ranks means the same chain, and then the chain's own order decides --
                    // the alignment cannot, having collapsed the chain to one symbol. Reversed where
                    // the parent crosses the chain from its end boundary; both operands are in that
                    // one chain, so they share the flag and either one answers.
                    if (all_ranked && all_indexed && ea.align_rank == ec.align_rank
                        && ea.chain_index != ec.chain_index) {
                        return ea.chain_backward ? ea.chain_index > ec.chain_index
                                                 : ea.chain_index < ec.chain_index;
                    }
                    if (ea.position != ec.position) {
                        return ea.position < ec.position;
                    }
                    return ea.record_key < ec.record_key;
                });
                if (pinned_phase.count(entries[pidx].record_key) != 0) {
                    ++g_group_parent_pinned;
                } else {
                    // No PhaseCall for the parent, so nothing ties this group's orientation to it.
                    ++g_group_parent_unpinned;
                }
                group.insert(group.end(), kv.second.begin(), kv.second.end());
                grouped_sites += group.size();
                ++grouped_groups;
                groups.push_back(std::move(group));
                auto found_ctx = parent_context.find(kv.first.first);
                gctx.push_back(found_ctx != parent_context.end() ? &found_ctx->second : nullptr);
                auto ps = phase_set_of.find(kv.first.first);
                gps.push_back(ps != phase_set_of.end() ? ps->second
                                                       : numeric_limits<size_t>::max());
            }
            // The chains that could not be grouped, carried through unchanged. Without this they
            // would be replaced by groups that do not contain them and their sites would never be
            // decoded at all -- which is the failure the global flag existed to prevent, and the
            // reason relaxing it needs this line and not just a looser test.
            for (vector<size_t>& kc : kept_chains) {
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
        // Nothing here belongs to this generation, so the decode has nothing to produce. Every
        // site is clamped to a delta at its settled genotype, the result loop skips all of them
        // and so does the phasing loop -- the chain is present only to carry transition context
        // for a site of this generation, and there is none.
        //
        // This is not a rare case, it is the usual one. Only nested sites arrive after generation
        // 0, and nested sites are held out of these runs by construction, so every generation
        // after the first re-decoded the whole contig and threw the answer away: five passes over
        // chr20's 192,045 top-level sites, 5.9 s each, 30 s of a 328 s run.
        // How much of this decode is live, and how much of the rest is pinned.
        //
        // Reported because it sizes the one optimisation left in this loop. A per-site pin zeroes
        // every state but one, so it severs the path: a site of this generation depends only on
        // the sites between its bracketing pins, and the chain could be cut there and decoded in
        // pieces. On chr20 that is 209,244 sites decoded for 17,199 live ones at generation 1 and
        // 211,952 for 22 at generation 4 -- about 20 s of a 187 s run.
        //
        // Not done, and the reason is in the two numbers rather than in the idea: it needs the
        // phase set to keep coming from the whole chain rather than from a piece of it, and it
        // rests on every pin being accepted, where `window_phasing` may decline one whose pair
        // cannot spell the site's constrained genotype. Both are checkable; neither is free.
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
        for (size_t k = 0; k < indices.size(); ++k) {
            const size_t idx = indices[k];
            const Entry& e = entries[idx];
            LinkageModel::Site s;
            s.position = e.position;
            s.unpositioned = e.unpositioned;
            s.num_alleles = e.num_alleles;
            s.ploidy = e.ploidy;
            // Stage B: the step from the previous site measured along the parent's settled
            // traversals rather than along the reference, one distance per strand.
            //
            // Only between two children of the SAME parent. The first entry of a per-parent group is
            // the parent itself, whose frames are measured along *its* parent's traversal and so are
            // in a different frame entirely; and a cross-parent step needs the distance from the
            // earlier parent's end to the later parent's start, which is not stored. Both fall back
            // to the reference difference, and both are counted.
            //
            // Slot 0 is strand a by construction, not by convention: `set_frame` writes slot 0 along
            // `PhaseCall::trav_first`, and that traversal sits on `hap_first`, which is the first of
            // the ordered pair the HMM's state is. Nothing re-derives the correspondence.
            //
            // Diploid chains only. A haploid chain has one strand, so the per-strand form has nothing
            // to say there, and the nested per-strand chains already measure their own frame gap.
            if (e.unpositioned && e.ploidy == 2) {
                ++unpos_diploid_sites;
            }
            if (frame_gap_mode > 0 && k > 0 && chain_ploidy == 2 && e.parent_record_key != 0) {
                // Does this step touch a site with no reference position? That is the population
                // whose distance CANNOT fall back, so it is the one where an unframed step goes
                // uniform and severs the chain rather than merely losing precision.
                if (e.unpositioned || entries[indices[k - 1]].unpositioned) {
                    ++unpos_steps;
                }
                const Entry& prev = entries[indices[k - 1]];
                if (prev.record_key == e.parent_record_key) {
                    // The previous site IS this one's parent -- the first entry of a per-parent
                    // group, and the only step in a group that is not same-parent.
                    //
                    // `frame_offset` is measured from the BEGINNING of the parent's settled
                    // traversal, so the child's own offset IS the distance from the parent to it.
                    // There is nothing to compose: the frame_end + frame_total + anchor recipe joins
                    // two DIFFERENT parents' traversals, which cannot arise inside one group. Getting
                    // this wrong left one step per group -- 1,591 of them at generation 1 on chr20,
                    // exactly the group count -- on the reference difference.
                    int64_t g[2] = {-1, -1};
                    for (int slot = 0; slot < 2; ++slot) {
                        if (e.frame_offset[slot] >= 0) {
                            g[slot] = (int64_t)e.frame_offset[slot];
                        }
                    }
                    if (g[0] >= 0 && g[1] >= 0) {
                        if (g[0] != g[1]) {
                            ++fg_differ;
                        }
                        if (frame_gap_mode == 1) {
                            s.gap_to_previous[0] = (g[0] + g[1]) / 2;
                            ++fg_single;
                        } else {
                            s.gap_to_previous[0] = g[0];
                            s.gap_to_previous[1] = g[1];
                            ++fg_pair;
                        }
                        ++fg_parent_step;
                        if (e.unpositioned || prev.unpositioned) { ++unpos_framed; }
                    } else if (g[0] >= 0 || g[1] >= 0) {
                        s.gap_to_previous[g[0] >= 0 ? 0 : 1] = g[0] >= 0 ? g[0] : g[1];
                        ++fg_one_slot;
                        ++fg_parent_step;
                        if (e.unpositioned || prev.unpositioned) { ++unpos_framed; }
                    } else {
                        ++fg_no_frame;
                        if (e.unpositioned || prev.unpositioned) { ++unpos_uniform; }
                    }
                } else if (prev.parent_record_key != e.parent_record_key) {
                    ++fg_cross_parent;
                } else {
                    int64_t g[2] = {-1, -1};
                    bool any_mirrored = false;
                    for (int slot = 0; slot < 2; ++slot) {
                        if (prev.frame_offset[slot] >= 0 && e.frame_offset[slot] >= 0) {
                            // START to START, because that is what the reference fallback measures:
                            // `site_gap` differences two `position` values, and a snarl's position
                            // is where it begins. Measuring the frame end-to-start instead -- which
                            // is how `Entry` describes a same-parent gap -- silently changes the
                            // convention as well as the metric, and for two adjacent siblings the
                            // two answers are 0 and the whole span of the earlier one.
                            //
                            // It also made the gap unmeasurable for most steps. A span is inclusive
                            // of its closing boundary node and adjacent siblings share that node, so
                            // end-to-start comes out negative for them: 10,103 of chr20's 18,235
                            // same-parent steps read as "no frame" and were exempted.
                            //
                            // Absolute, because offsets run in TRAVERSAL order while the sites are
                            // in reference (or alignment) order, and those need not agree in
                            // direction. A step measured against the traversal's own direction is
                            // counted, not refused.
                            const int64_t d =
                                (int64_t)e.frame_offset[slot] - (int64_t)prev.frame_offset[slot];
                            g[slot] = d >= 0 ? d : -d;
                            if (d < 0) {
                                any_mirrored = true;
                            }
                        }
                    }
                    if (any_mirrored) {
                        ++fg_mirrored;
                    }
                    // How often the traversal distance and the reference distance are the SAME
                    // number. Two adjacent siblings share a boundary node, so wherever the parent's
                    // traversal follows the reference through them the two measures coincide
                    // exactly, and this arm can only move the steps where they do not. Measured,
                    // rather than inferred from how a span is taken.
                    if (g[0] >= 0) {
                        const int64_t refd = (int64_t)e.position - (int64_t)prev.position;
                        if (g[0] == (refd >= 0 ? refd : -refd)) {
                            ++fg_equals_ref;
                        }
                    }
                    if (g[0] >= 0 && g[1] >= 0) {
                        // Counted in both modes, so mode 1 can say how many steps it flattened.
                        if (g[0] != g[1]) {
                            ++fg_differ;
                        }
                        if (frame_gap_mode == 1) {
                            // One scalar for both strands. Written to slot 0 alone -- `site_gap`
                            // uses a single known slot for both, so this needs no second copy.
                            s.gap_to_previous[0] = (g[0] + g[1]) / 2;
                            ++fg_single;
                        } else {
                            s.gap_to_previous[0] = g[0];
                            s.gap_to_previous[1] = g[1];
                            ++fg_pair;
                        }
                        if (e.unpositioned || prev.unpositioned) { ++unpos_framed; }
                    } else if (g[0] >= 0 || g[1] >= 0) {
                        // One traversal enters this chain and the other does not, or one step came
                        // out negative. The known distance stands for both strands, which is what a
                        // single value did.
                        //
                        // In a DIPLOID chain this should not happen at all: a chain only one settled
                        // traversal crosses has one copy, so ploidy 1, so it belongs to the
                        // per-strand chains. It was 0 here until off-reference chains were admitted
                        // and is 1,409 with them, so the split below asks which half is wrong -- did
                        // descent claim both traversals carry it (then a frame failed to measure),
                        // or did it claim one (then the ploidy is wrong upstream)?
                        s.gap_to_previous[g[0] >= 0 ? 0 : 1] = g[0] >= 0 ? g[0] : g[1];
                        ++fg_one_slot;
                        if (e.parent_trav == -2) {
                            ++fg_one_both_claimed;
                        } else if (e.parent_trav >= 0) {
                            ++fg_one_single_claimed;
                        } else {
                            ++fg_one_undetermined;
                        }
                    } else {
                        ++fg_no_frame;
                        if (e.unpositioned || prev.unpositioned) { ++unpos_uniform; }
                    }
                }
            }
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

        // V1: how much pair correlation survives a margin of a given length.
        //
        // Per step the probability that neither strand switches is (1-rho_a)(1-rho_b), so over N
        // steps it is exp(-N * mean(-log((1-rho_a)(1-rho_b)))). Accumulating the mean per step
        // answers both questions at once: what the shipped 250-site margin retains, and how many
        // sites would be needed to reach any target. Top-level chains only -- that is where the
        // windowing runs.
        if (generation == 0) {
            for (size_t k = 1; k < sites.size(); ++k) {
                const std::pair<size_t, size_t> gap = site_gap(sites[k - 1], sites[k]);
                const double ra = model.switch_probability(gap.first);
                const double rb = model.switch_probability(gap.second);
                const double stay = max(1e-300, (1.0 - ra) * (1.0 - rb));
                g_margin_neglog.fetch_add((size_t)(-log(stay) * 1e9));
                ++g_margin_steps;
            }
        }

        vector<vector<double>> posteriors;
        if (chain_ploidy == 1) {
            posteriors = model.haploid_posteriors(sites);
        } else if (chain_i < chain_context.size() && chain_context[chain_i] != nullptr) {
            // A child group, conditioned on its parent's context message. One decode over the
            // group, not over the contig.
            //
            // Harvests too. Without this a grouped generation stores nothing, so the NEXT
            // generation finds no context for its parents and falls back to the contig chain --
            // grouping at one depth silently disabling it at the next, which is what the coverage
            // counters showed.
            vector<char> mask(indices.size(), 0);
            size_t masked = 0;
            for (size_t k = 0; k < indices.size(); ++k) {
                if (wanted_parents.count(entries[indices[k]].record_key) != 0) {
                    mask[k] = 1;
                    ++masked;
                }
            }
            vector<vector<double>> ga, gb;
            model.segment_posteriors(sites, 0, sites.size(), chain_context[chain_i], nullptr,
                                     posteriors, masked ? &ga : nullptr, masked ? &gb : nullptr,
                                     masked ? &mask : nullptr);
            for (size_t k = 0; masked && k < indices.size(); ++k) {
                if (k >= ga.size() || ga[k].empty() || k >= gb.size() || gb[k].empty()) {
                    continue;
                }
                vector<double> ctx = ga[k];
                double sum = 0.0;
                for (size_t j = 0; j < ctx.size() && j < gb[k].size(); ++j) {
                    ctx[j] *= gb[k][j];
                    sum += ctx[j];
                }
                if (sum > 0.0) {
                    for (double& v : ctx) {
                        v /= sum;
                    }
                    parent_context[entries[indices[k]].record_key] = std::move(ctx);
                }
            }
        } else if (wanted_parents.empty()) {
            posteriors = model.posteriors(sites);
        } else {
            // Same windowed sweep, harvesting the context message at the parents of sites this
            // generation has not reached yet. Sparse: the mask is what keeps this from costing
            // 9.6 kB a site in both directions.
            vector<char> mask(indices.size(), 0);
            size_t masked = 0;
            for (size_t k = 0; k < indices.size(); ++k) {
                if (wanted_parents.count(entries[indices[k]].record_key) != 0) {
                    mask[k] = 1;
                    ++masked;
                }
            }
            if (masked == 0) {
                posteriors = model.posteriors(sites);
            } else {
                vector<vector<double>> context;
                posteriors = model.posteriors_with_context(sites, mask, context);
                for (size_t k = 0; k < indices.size() && k < context.size(); ++k) {
                    if (!context[k].empty()) {
                        parent_context[entries[indices[k]].record_key] = std::move(context[k]);
                    }
                }
            }
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
                // Counted, not emitted. The settled pair is in `final_i`/`final_j` above and the
                // record is built from it afterwards, so there is no line to rewrite and no VCF
                // numbering to render into -- which is why the four vcf_allele_of lookups and the
                // `unrenderable` counter they fed are gone: they existed only to drop a patch that
                // named an allele the already-written line had no ALT for.
                //
                // The posterior is kept, though. The quality of a moved genotype is the phred
                // complement of this posterior, discounted and capped, and the posterior exists
                // nowhere else -- so it is handed out rather than baked into a patch.
                moved_quality_by_record[e.record_key] =
                    std::make_pair(post[best], (double)e.explained_share);
                ++moved;
            }
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
            finish_phase_call(pc, e, phase_fallback);
            phasing_out->push_back(pc);
        }
    }

    if (!parent_context.empty() && generation == 0) {
        size_t bytes = 0;
        for (const auto& kv : parent_context) {
            bytes += kv.second.capacity() * sizeof(double) + 48;
        }
        #pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": context messages held for "
                  << parent_context.size() << " parents of not-yet-called sites, "
                  << (bytes / (1024.0 * 1024.0)) << " MB" << std::endl;
    }

    if (ao_pairs + ao_unranked > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": " << ao_pairs
                  << " same-parent adjacent pairs ranked by the settled traversals' alignment, "
                  << ao_reordered << " of them in a different order than the reference; "
                  << ao_unranked << " pairs unranked; " << ao_group_unranked
                  << " groups fell back to reference order entire, for want of a rank on every"
                  << " member" << (no_align_order ? " (order forced to the reference)" : "")
                  << std::endl;
    }

    if (unpos_diploid_sites > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": " << unpos_diploid_sites
                  << " decoded diploid sites have no reference position; " << unpos_steps
                  << " steps touch one, of which " << unpos_framed
                  << " took a traversal distance and " << unpos_uniform
                  << " had none and went uniform" << std::endl;
    }

    if (frame_gap_mode > 0 && (fg_pair + fg_single + fg_one_slot + fg_no_frame
                               + fg_cross_parent) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": frame gaps mode "
                  << frame_gap_mode << " -- " << (frame_gap_mode == 1 ? fg_single : fg_pair)
                  << " same-parent steps spaced along the settled traversals ("
                  << fg_differ << " where the two strands' distances differ"
                  << (frame_gap_mode == 1 ? ", flattened to their mean" : "") << "), "
                  << fg_one_slot << " from one traversal only, " << fg_no_frame
                  << " with no frame at all, " << fg_cross_parent
                  << " cross-parent -- the last two on the reference difference; "
                  << fg_mirrored << " measured against the traversal's own direction, "
                  << fg_equals_ref << " exactly equal to the reference difference; "
                  << fg_parent_step << " parent-to-first-child steps taken from the child's own"
                  << " offset along the parent's traversal; of the one-slot steps, "
                  << fg_one_both_claimed << " where descent claimed BOTH traversals carry the chain, "
                  << fg_one_single_claimed << " where it claimed one, " << fg_one_undetermined
                  << " undetermined" << std::endl;
    }

    if (generation == 0 && g_margin_steps.load() > 0) {
        const double mean = (g_margin_neglog.load() / 1e9) / (double)g_margin_steps.load();
#pragma omp critical (cerr)
        std::cerr << "[vg call] margin: " << g_margin_steps.load() << " top-level steps, mean"
                  << " -log P(no switch) = " << mean << " per step; correlation retained across "
                  << params.margin << " sites = " << exp(-(double)params.margin * mean)
                  << "; sites needed for 0.05 = " << (mean > 0 ? (-log(0.05) / mean) : -1.0)
                  << ", for 0.01 = " << (mean > 0 ? (-log(0.01) / mean) : -1.0) << std::endl;
    }

    if (generation > 0 && (g_pin_applied.load() + g_pin_declined.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": pins -- "
                  << g_pin_applied.load() << " applied, " << g_pin_declined.load()
                  << " REFUSED (the pinned pair cannot spell the constrained genotype, so that"
                  << " site's orientation is free); groups whose parent was pinnable: "
                  << g_group_parent_pinned.load() << ", not pinnable: "
                  << g_group_parent_unpinned.load() << std::endl;
    }

    if (generation > 0 && (grouped_groups > 0 || g_grp_vetoed.load() > 0)) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": grouping declines so far -- "
                  << g_grp_no_parent.load() << " sites with no parent key, "
                  << g_grp_no_entry.load() << " whose parent has no live entry, "
                  << g_grp_ploidy.load() << " on a parent/child ploidy difference; "
                  << g_grp_vetoed.load() << " chains kept ungrouped in total" << std::endl;
    }

    if (grouped_groups > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": decoded " << grouped_sites
                  << " sites in " << grouped_groups << " per-parent groups instead of the contig"
                  << " chain; " << g_grp_no_context.load()
                  << " groups had no stored parent context (children recorded after their"
                  << " generation's decode), " << g_grp_ploidy.load()
                  << " declined on a parent/child ploidy difference" << std::endl;
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
            size_t allele_first = 0;
            size_t allele_second = 0;
            int trav_first = -1;
            int trav_second = -1;
            size_t ploidy = 2;
            /// Whether the parent's own strand order was determined by the panel or fell through to
            /// sorting. A child inherits its strand from that order, so an arbitrary parent order
            /// makes the child's strand arbitrary too -- and saying so is the difference between a
            /// determined strand and a coin flip presented as one.
            bool order_arbitrary = false;
            /// The strand the parent itself sits on, when the parent is a nested haploid site. A
            /// nested parent occupies ONE haplotype, so everything inside it is on that haplotype --
            /// which the identity match below cannot discover, because such a parent has
            /// trav_first == trav_second and ploidy 1, so the match can only ever return strand 0 or
            /// no strand, never strand 1.
            int8_t nested_strand = -1;
            /// Whether this entry came from the accumulated seed at the top of the pass or from an
            /// in-loop insert during placement. The seed is re-read every pass before placement
            /// runs, so a stale in-loop value should never be what a child reads -- checked rather
            /// than argued.
            bool from_seed = true;
        };
        // record_key -> entry index, for the nested sites placed below.
        unordered_map<size_t, size_t> nested_entry_of;
        // Strands whose inherited haplotype turned out not to carry the settled allele. Split by
        // whether the site has a VCF line, because only those reach the mosaic -- and without the
        // split the figure cannot be reconciled against anything, which is the same defect this
        // stage is fixing in the mosaic's own wildcard column.
        size_t hap_contradicted = 0, hap_contradicted_emitted = 0;
        unordered_map<size_t, ParentPhase> by_key;
        by_key.reserve(phasing_out->size() * 2);
        for (const PhaseCall& pc : *phasing_out) {
            by_key[pc.record_key] = ParentPhase{pc.phase_set, pc.hap_first, pc.hap_second,
                                                pc.allele_first, pc.allele_second,
                                                pc.trav_first, pc.trav_second, pc.ploidy,
                                                pc.order_arbitrary, pc.nested_strand,
                                                /*from_seed*/ true};
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
        // Keyed by (PHASE SET, strand), not (contig, strand). A phase set is exactly the unit
        // within which a phase is comparable, and a strand only means anything inside one: keying on
        // the contig put every same-depth site under unrelated parents into one chain, and on a
        // contig with several chains it merged two unrelated chains' strand 0 into a single haploid
        // run. Linking sites whose haplotypes have no correspondence is worse than not linking them.
        map<pair<size_t, uint8_t>, vector<size_t>> by_strand;
        // For the gate: how many sites would have been pooled across a phase-set boundary under the
        // old key. Zero is the claim; a count proves the change was not cosmetic.
        map<pair<uint32_t, uint8_t>, size_t> old_key_groups;
        size_t crossed_phase_set = 0;
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
                if (parent.ploidy == 1 && parent.nested_strand >= 0) {
                    // The parent is itself a nested haploid site, so it occupies ONE haplotype and
                    // everything inside it is on that haplotype. Its strand is the answer, and the
                    // identity match below cannot find it: such a parent has
                    // trav_first == trav_second and ploidy 1, so the match can only ever return
                    // strand 0. On chr20 that put all 448 depth->=2 sites under a strand-1 parent on
                    // the wrong haplotype -- while the mosaic, reading the parent's wildcard
                    // hap_first, simultaneously reported them as belonging to no haplotype.
                    if (e.parent_trav >= 0 && parent.trav_first == e.parent_trav) {
                        strand = parent.nested_strand;
                    }
                } else if (e.parent_trav >= 0) {
                    // Both assignments require a DIPLOID parent, not just the second one. A strand
                    // is a claim that the locus has two haplotypes and this allele sits on one of
                    // them, and that claim is only true where the parent has two copies.
                    //
                    // Without the guard on `strand = 0`, a haploid top-level site -- chrX outside the
                    // pseudoautosomal regions, or chrY -- gives its children strand 0, because the
                    // first branch above needs `parent.nested_strand >= 0` and a top-level parent has
                    // none, so this one matches `trav_first` and assigns unconditionally. The child is
                    // then rendered "a|." , which says the locus is diploid and the other strand
                    // carries nothing. On a haploid contig there is no other strand to be empty: the
                    // right rendering is a bare `a`.
                    //
                    // Measured: chrX's haploid interior carried 8,056 "1|." records against the
                    // pre-refactor arm's 1,965, and chrX was the one contig whose F1 fell
                    // (0.94939 -> 0.93643) while all 22 autosomes rose. Autosomes are unaffected --
                    // their parents are ploidy 2, and a nested haploid parent reaches the branch
                    // above instead.
                    if (parent.ploidy == 2) {
                        if (parent.trav_first == e.parent_trav) {
                            strand = 0;
                        } else if (parent.trav_second == e.parent_trav) {
                            strand = 1;
                        }
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
                    // The strand came off the parent's pair. If the panel did not determine that
                    // pair's order, this strand is that coin flip one level down, and it was
                    // previously reported as determined.
                    pc.order_arbitrary = pc.order_arbitrary || parent.order_arbitrary;
                }
                // Which entry this nested PhaseCall came from, so the haplotype inherited above can
                // be checked against the allele the site actually settles on -- which is not known
                // until the per-strand pass below has run.
                nested_entry_of[pc.record_key] = idx;
                finish_phase_call(pc, e, phase_fallback);
                // Registered only now. This used to sit above the block that assigns
                // `pc.trav_first`/`trav_second`, so the entry it left for a deeper site carried the
                // default -1. Harmless today only by accident -- the seed at the top of every pass
                // re-reads the accumulated phasing before placement runs, and a parent and child
                // never resolve in the same pass -- but it made the shallowest-first sweep loop,
                // which exists precisely to serve depth >= 2, dead code that happened to agree with
                // the right answer. Measured: 0 of chr20's 2,409 depth >= 2 sites read an in-loop
                // entry.
                by_key[pc.record_key] = ParentPhase{
                    pc.phase_set, pc.hap_first, pc.hap_second,
                    pc.allele_first, pc.allele_second,
                    pc.trav_first, pc.trav_second, pc.ploidy,
                    pc.order_arbitrary, pc.nested_strand, /*from_seed*/ false};
                phasing_out->push_back(pc);
                if (strand >= 0) {
                    // Only sites with one strand are grouped: step three chains within a single
                    // haplotype, and a site with no strand belongs to no such chain.
                    by_strand[make_pair(pc.phase_set, (uint8_t)strand)].push_back(idx);
                    // What the old contig-keyed grouping would have done, for the gate below.
                    auto ok = make_pair(e.contig, (uint8_t)strand);
                    if (old_key_groups.count(ok) == 0) {
                        old_key_groups[ok] = pc.phase_set;
                    } else if (old_key_groups[ok] != pc.phase_set) {
                        ++crossed_phase_set;
                    }
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
                finish_phase_call(pc, e, phase_fallback);
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
        // Stage 15': the CROSS-PARENT measurement, inert. Stage 14's 0.606% / 91.05% were measured
        // between SIBLINGS -- children of one parent -- while a (phase_set, strand) group spans a
        // whole contig on one haplotype, where adjacent sites usually have different parents. So the
        // neutrality prediction for consuming these frames rested on a statistic about a different
        // population, and this measures the right one.
        //
        // Order key compared here is (parent's position, offset along the parent's settled traversal),
        // which is the snarl-tree key truncated to one level. Distance is formed pairwise -- tail of
        // the earlier site's parent traversal, plus the anchor gap, plus the later site's offset --
        // never as a difference of two absolute coordinates, which is the error that sank the first
        // draft of this design.
        unordered_map<size_t, uint32_t> position_of;
        position_of.reserve(entries.size() * 2);
        for (const Entry& pe : entries) {
            if (!pe.retracted) {
                position_of[pe.record_key] = pe.position;
            }
        }

        size_t singleton_groups = 0;
        size_t total_frame_steps = 0, total_ref_steps = 0;
        // Stage C, in the per-strand haploid chains: alignment order against the key it precedes.
        size_t hao_pairs = 0, hao_reordered = 0, hao_unranked = 0, hao_group_unranked = 0;
        size_t hao_same_chain_indexed = 0, hao_same_chain_unindexed = 0;
        for (auto& group : by_strand) {
            vector<size_t>& idxs = group.second;
            // Singletons are NOT skipped. A one-site group has nothing to link against, so linkage
            // cannot move it on the strength of a neighbour -- but `freq_prior` defaults to 5 and
            // acts on a chain of one, so the posterior still differs from the raw likelihood and the
            // site still has to be phased. Skipping it left the genotype at its per-site value and
            // kept the site out of the phasing entirely. The diploid path fixed exactly this defect
            // for its own singletons, where it cost 258 chr20 sites their place in the mosaic.
            if (idxs.size() == 1) {
                ++singleton_groups;
            }
            // Stage 15': ordered by a TUPLE -- (parent's anchor, offset along the parent's settled
            // traversal, record key) -- not by an arithmetic coordinate. A tuple comparison preserves
            // subtree containment by construction, so no ordering claim is being made that could be
            // violated; composing the two into one number would add a reference anchor to a
            // haplotype-walk length and can invert the order of sites under different parents.
            //
            // Measured on chr20: 0 of 5,540 same-parent adjacent pairs reorder under this key, so the
            // order is not what this change moves -- the distances are.
            const int slot = (int)group.first.second <= 1 ? (int)group.first.second : 0;
            // Same total-key discipline as the diploid groups, and for the same reason: this
            // comparator had TWO skippable keys, the alignment rank and the frame offset, so it was
            // the more exposed of the pair.
            //
            // Same parent: the alignment of its two settled traversals first -- that is the one
            // order defined for a chain only one haplotype crosses, which is what these per-strand
            // chains are made of. Then the offset along the settled traversal, then the reference
            // position. An absent rank or an absent offset sorts last within its key.
            auto anchor_of = [&](const Entry& e) {
                auto it = position_of.find(e.parent_record_key);
                return it != position_of.end() ? it->second : e.position;
            };
            // This bucket spans parents, so the all-or-nothing decision is per PARENT: a key is
            // usable for a pair only if it is usable for every sibling either of them could be
            // compared against, and the anchor is what separates one parent's children from
            // another's. Both keys are decided this way -- the frame offset has exactly the same
            // problem as the rank and predates it.
            unordered_map<size_t, unsigned char> parent_keys;   // bit 0 all ranked, bit 1 all framed
            {
                unordered_map<size_t, unsigned char> missing;
                for (size_t idx : idxs) {
                    const Entry& e = entries[idx];
                    unsigned char& m = missing[e.parent_record_key];
                    if (e.align_rank < 0) { m |= 1; }
                    if (e.frame_offset[slot] < 0) { m |= 2; }
                    if (e.chain_index < 0) { m |= 4; }
                }
                for (const auto& kv : missing) {
                    unsigned char usable = 0;
                    if (!no_align_order && (kv.second & 1) == 0) { usable |= 1; }
                    if ((kv.second & 2) == 0) { usable |= 2; }
                    if (!no_align_order && (kv.second & 4) == 0) { usable |= 4; }
                    parent_keys[kv.first] = usable;
                    if ((usable & 1) == 0) { ++hao_group_unranked; }
                }
            }
            sort(idxs.begin(), idxs.end(), [&](size_t a, size_t b) {
                const Entry& ea = entries[a];
                const Entry& eb = entries[b];
                const uint32_t aa = anchor_of(ea), ab = anchor_of(eb);
                if (aa != ab) {
                    return aa < ab;
                }
                // Same anchor, so the same parent, so one lookup answers for both.
                auto uk = parent_keys.find(ea.parent_record_key);
                const unsigned char usable = uk != parent_keys.end() ? uk->second : 0;
                if ((usable & 1) && ea.align_rank != eb.align_rank) {
                    return ea.align_rank < eb.align_rank;
                }
                // Same chain: its own order, not the frame offset. The offset does separate two
                // snarls of one chain, but as a bp walk along one traversal rather than as anything
                // intrinsic; the chain index is the exact answer the decomposition already holds.
                if ((usable & 1) && (usable & 4) && ea.align_rank == eb.align_rank
                    && ea.chain_index != eb.chain_index) {
                    return ea.chain_backward ? ea.chain_index > eb.chain_index
                                             : ea.chain_index < eb.chain_index;
                }
                if ((usable & 2) && ea.frame_offset[slot] != eb.frame_offset[slot]) {
                    return ea.frame_offset[slot] < eb.frame_offset[slot];
                }
                if (ea.position != eb.position) {
                    return ea.position < eb.position;
                }
                return ea.record_key < eb.record_key;
            });
            // Post-sort, because before it `idxs` is in no meaningful order: how often the
            // alignment rank -- which now leads the tuple -- put a same-parent pair in a different
            // order than the frame offset and reference position it precedes. The diploid groups
            // count the same thing pre-sort; both should read zero wherever both keys exist.
            for (size_t k = 1; k < idxs.size(); ++k) {
                const Entry& e1 = entries[idxs[k - 1]];
                const Entry& e2 = entries[idxs[k]];
                if (e1.parent_record_key != e2.parent_record_key || e1.parent_record_key == 0) {
                    continue;
                }
                if (e1.align_rank >= 0 && e2.align_rank >= 0 && e1.align_rank == e2.align_rank) {
                    // Same chain. Separated by the chain index now, where reference position was
                    // the only thing that could before.
                    if (e1.chain_index >= 0 && e2.chain_index >= 0
                        && e1.chain_index != e2.chain_index) {
                        ++hao_same_chain_indexed;
                    } else {
                        ++hao_same_chain_unindexed;
                    }
                    continue;
                }
                if (e1.align_rank < 0 || e2.align_rank < 0) {
                    ++hao_unranked;
                    continue;
                }
                ++hao_pairs;
                const bool fallback_agrees =
                    (e1.frame_offset[slot] >= 0 && e2.frame_offset[slot] >= 0
                     && e1.frame_offset[slot] != e2.frame_offset[slot])
                        ? e1.frame_offset[slot] < e2.frame_offset[slot]
                        : (e1.position != e2.position ? e1.position < e2.position : true);
                if (!fallback_agrees) {
                    ++hao_reordered;
                }
            }
            vector<LinkageModel::Site> sites;
            sites.reserve(idxs.size());
            size_t frame_steps = 0, ref_steps = 0;
            for (size_t k = 0; k < idxs.size(); ++k) {
                const size_t idx = idxs[k];
                const Entry& e = entries[idx];
                LinkageModel::Site s;
                s.position = e.position;
                s.unpositioned = e.unpositioned;
                // Stage 15': the step from the previous site, measured along the haplotype's own
                // traversal where both sites hang off the same parent and both have a frame. Passed
                // as an explicit distance rather than folded into `position`, so it cannot affect the
                // order that was already fixed above.
                //
                // Cross-parent steps still fall back to the reference difference: forming them in the
                // frame needs the distance from the earlier parent's END to the later parent's START,
                // and only parent start positions are stored. That is 12% of adjacent pairs on chr20
                // and is what remains before the reference dependency is fully gone.
                // Used only where the REFERENCE distance is unavailable, not in preference to it.
                //
                // Measured, twice, by two unrelated derivations: spacing nested chain steps along a
                // traversal instead of along the reference costs about 0.0005 of JointIndel on chr20
                // (0.92390 -> 0.92338 here; 0.92390 -> 0.92329 for stage 15(b), which derived the
                // distance from per-site called alleles instead). Same magnitude, same direction, so
                // it is not the labelling -- 15' fixed that -- and not the derivation. It is the
                // change itself: a traversal distance is a worse predictor for this transition model
                // than the reference distance, wherever a reference distance exists.
                //
                // It exists everywhere today, so this is inert now. Under a covering reference a
                // nested site has no reference position at all, and then the frame is not an
                // improvement to choose but the only measure there is -- which is what this machinery
                // is for.
                if (k > 0) {
                    const Entry& prev = entries[idxs[k - 1]];
                    // `unpositioned`, not `position > 0`. Position is an anchor for a chain the
                    // reference does not cross -- its parent's start -- so it is non-zero and this
                    // test would pass, the frame would never be consulted, and the model would
                    // difference two identical anchors. Absence is a fact to be carried, not
                    // inferred from a coordinate that happens to be zero.
                    const bool have_reference = (!prev.unpositioned && !e.unpositioned
                                                 && prev.position > 0 && e.position > 0);
                    if (!have_reference
                        && prev.parent_record_key == e.parent_record_key
                        && prev.frame_offset[slot] >= 0 && e.frame_offset[slot] >= 0) {
                        // Slot 0 only: this is a per-strand haploid chain, so there is one distance
                        // and `site_gap` uses it for both of the pair's strands.
                        //
                        // Start to start, matching the reference difference this stands in for. The
                        // diploid groups measure the same way, for the same reason.
                        s.gap_to_previous[0] =
                            (int64_t)e.frame_offset[slot] - (int64_t)prev.frame_offset[slot];
                        if (s.gap_to_previous[0] < 0) {
                            // The tuple sort put them in this order, so a negative step means the
                            // frames disagree with it -- do not feed the model a wrapped distance.
                            s.gap_to_previous[0] = -1;
                            ++ref_steps;
                        } else {
                            ++frame_steps;
                        }
                    } else {
                        ++ref_steps;
                    }
                }
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
            total_frame_steps += frame_steps;
            total_ref_steps += ref_steps;
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
                // As above: the answer is recorded in `nested_regenotyped` and in
                // `final_i`/`final_j`, and the record is built from it later. No line exists to
                // patch -- but the posterior is kept, for the quality of the moved genotype.
                moved_quality_by_record[e.record_key] =
                    std::make_pair(post[best], (double)e.explained_share);
                ++moved;
            }
        }

        // The PhaseCalls for these sites were built before step three ran, so they still describe
        // the pre-linkage allele. The renderer refuses a phase that disagrees with the genotype
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
        if (total_frame_steps > 0 || total_ref_steps > 0) {
#pragma omp critical (cerr)
            std::cerr << "[vg call] frames: " << total_frame_steps
                      << " chain steps spaced along the settled parent traversal, "
                      << total_ref_steps << " still by reference difference (cross-parent, or a"
                      << " frame the sort disagreed with)" << std::endl;
        }
        if (hao_pairs + hao_unranked > 0) {
#pragma omp critical (cerr)
            std::cerr << "[vg call] per-strand chains: " << hao_pairs
                      << " same-parent adjacent pairs ordered by the settled traversals'"
                      << " alignment, " << hao_reordered << " of them in a different order than the"
                      << " frame offset and reference position it precedes; " << hao_unranked
                      << " pairs unranked; " << hao_group_unranked
                      << " parents fell back to reference order entire; "
                      << hao_same_chain_indexed << " pairs are two snarls of ONE chain and are"
                      << " ordered by its index, " << hao_same_chain_unindexed
                      << " of those had no index" << std::endl;
        }
        // Stage 15': the cross-parent numbers, which are the ones the decision to consume these
        // frames has to rest on.
        // Stage 15(a): what the pooling fix changed, printed so it is not taken on faith.
        if (crossed_phase_set > 0 || singleton_groups > 0) {
#pragma omp critical (cerr)
            std::cerr << "[vg call] nested strands: " << by_strand.size()
                      << " groups keyed by (phase set, strand); " << crossed_phase_set
                      << " sites the old contig-keyed grouping would have pooled across a phase-set"
                      << " boundary, " << singleton_groups
                      << " groups of one now linked and phased instead of skipped" << std::endl;
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

        // The inherited haplotype is a claim about the panel, and until now nothing checked it.
        //
        // A nested site takes its haplotype from whichever of its parent's strands carries the
        // chain. That says the parent's allele sits on panel haplotype h; it does not say h carries
        // the allele *this* site settled on, and the per-strand pass above may just have moved that
        // allele. So the mosaic could name a haplotype demonstrably carrying something else, which
        // is worse than naming none: a consumer walking the haplotype would read a different
        // sequence than the record states.
        //
        // Checked against the panel matrix the site was genotyped with, and dropped to the wildcard
        // where it fails. That deliberately claims less than before, so the panel-unexplained count
        // rises by exactly this figure.
        for (PhaseCall& pc : *phasing_out) {
            auto at = nested_entry_of.find(pc.record_key);
            if (at == nested_entry_of.end()) {
                continue;
            }
            const Entry& e = entries[at->second];
            // The allele the site settled on, in compact space: the per-strand pass's answer where
            // it moved the genotype, otherwise the call it came in with.
            size_t settled = e.called_i;
            auto moved = nested_regenotyped.find(pc.record_key);
            if (moved != nested_regenotyped.end()) {
                settled = moved->second.compact;
            }
            auto carries = [&](size_t hap) {
                if (hap == LinkageModel::WILDCARD || hap >= n_haplotypes) {
                    return true;   // nothing claimed, nothing to contradict
                }
                const size_t at_hap = e.hap_offset + hap;
                if (at_hap >= hap_arena.size()) {
                    return true;
                }
                const int carried = (int)hap_arena[at_hap];
                // A haplotype absent from this site carries no allele here and cannot contradict
                // the call; only a haplotype that names a *different* allele does.
                return carried < 0 || (size_t)carried == settled;
            };
            if (!carries(pc.hap_first)) {
                pc.hap_first = LinkageModel::WILDCARD;
                ++hap_contradicted;
                hap_contradicted_emitted += pc.emitted;
            }
            if (!carries(pc.hap_second)) {
                pc.hap_second = LinkageModel::WILDCARD;
                ++hap_contradicted;
                hap_contradicted_emitted += pc.emitted;
            }
        }
        if (hap_contradicted > 0) {
#pragma omp critical (cerr)
            std::cerr << "[vg call] phasing: " << hap_contradicted
                      << " nested strands named a panel haplotype that does not carry the allele the"
                      << " site settled on, so the haplotype was dropped ("
                      << hap_contradicted_emitted << " on a site with a line, which is the number"
                      << " the mosaic can show)" << std::endl;
        }
    }


    // `unrenderable` is gone rather than reported as zero. It counted settled genotypes that named
    // a traversal the already-written line had no ALT for, which a patch could not add -- 1,465 on
    // chr20 at its peak, 23% of the contig's false positives. With the record built from the settled
    // genotype the allele list is chosen from that genotype, so the population it counted cannot
    // exist and a counter that can only read zero is worse than no counter.
    // `phase_fallback` is no longer reported, because it no longer describes anything the output
    // shows. It counted PhaseCalls whose settled pair had no ALT in the record's allele map, so the
    // phase had to be rendered against the line's own alleles instead. Phasing is now applied while
    // the record is built, from the traversal pair through the map that record was built with, so
    // `PhaseCall::allele_first/allele_second` are written and never read. Printing a count of a
    // fallback with no consequence -- 192,045 of chr20's 219,600 sites, every top-level one -- is
    // noise that reads like a defect. The dead conversion itself is removed with
    // `render_phase_pair`; it needs a rename rather than a deletion, because `allele_*` carries the
    // compact pair before that block overwrites it with the VCF one.

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
    return moved;
}

}
