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

/// Per-site pins offered to window_phasing, and those it refused because the pinned haplotype pair
/// cannot spell the genotype the site is constrained to. A refused pin on a group's PARENT frees the
/// whole group's orientation while its haploid siblings stay tied to the parent's phase.
static std::atomic<size_t> g_pin_applied(0), g_pin_declined(0);
/// Groups whose parent was never offered a pin at all -- no PhaseCall to pin it to.
static std::atomic<size_t> g_group_parent_unpinned(0), g_group_parent_pinned(0);

/// The distance between two adjacent sites in a chain, in bp.
///
/// A pair, because `transition_apply` takes a switch probability per strand -- the two haplotypes of
/// a diploid sample recombine independently. Both halves are the same number today: the explicit
/// per-strand distance this used to prefer came from the haplotype frame, and the frame was measured
/// to change no byte of the output on three contigs, so nothing sets one any more.
///
/// Clamped at 1: `switch_probability` reads a gap of 0 as 1 anyway, and a negative one would wrap.
static inline std::pair<size_t, size_t> site_gap(const LinkageModel::Site& prev,
                                                const LinkageModel::Site& next) {
    if (prev.unpositioned || next.unpositioned) {
        // EITHER, not both. A mixed pair is the common case -- a parent's children interleave chains
        // the reference crosses with chains it does not -- and differencing an anchor against a real
        // coordinate is not a distance. SIZE_MAX gives rho = 1.0, at which `transition_apply`'s
        // T = (1-rho)I + (rho/m)11' is exactly uniform: the chain forgets, which is what "unknown"
        // means here. A clamp to 1 would be the opposite claim -- perfect linkage between two sites
        // about which nothing is known.
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
/// What a settled parent implies about one of its children, from the parent's two settled
/// traversals and the child's crossing mask. Takes the traversals rather than the parent, so it can
/// be called from wherever the settled pair is in hand -- today only the grouping, which reads them
/// off the parent's Entry.
///
/// The mask is in TRAVERSAL terms, so it must be tested against the traversals and never against a
/// compact allele index -- the two agree only when every allele at the parent is panel-carried.
static inline LinkageCollector::Relation relate_to_parent(uint64_t crossing, int ta, int tb) {
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
        // One strand, so the path is over single haplotypes and `second` stays the wildcard
        // throughout -- there is no second strand for it to describe. Decoded into its own vector
        // rather than into `out` because the seam pin is a haplotype here and a PAIR at ploidy 2,
        // and `windowed_path` pins whatever type it is carrying.
        //
        // No `m < 2` guard, deliberately: the ploidy-2 path needs one because a pair state space
        // of fewer than two states cannot spell a heterozygote, and one strand has no such floor.
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
    // The entering message replaces the uniform start where there is one. It multiplies the first
    // site's emission rather than replacing it, but the collector supplies a POINT MASS, so on this
    // window that is a choice of state and not a nudge -- see `posteriors` for how much the
    // chain's own transitions soften it downstream. Uniform where the caller supplies none, or
    // where the message is the wrong length for this state space.
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
    e.explained_share = (float)explained_share;
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
/// record's phasing -- 1,528 extra strandless records on chr20. It was spelled out three times, in
/// the diploid chain, the nested-strand pass and the haploid pass, with that warning copied
/// alongside each, and the three had drifted: only the diploid one also widened `order_arbitrary`
/// on a fallback. Folding them in changed nothing, because the other two were ploidy-1 populations
/// where the pair is one allele twice and the test cannot fire -- asserted by byte-identity rather
/// than argued. Two of the three call sites have since gone with the per-strand pass; this is the
/// one door that is left, which is why there is nothing to drift against any more.
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

static std::atomic<size_t> g_grp_no_parent(0);
static std::atomic<size_t> g_grp_no_entry(0), g_grp_vetoed(0);
// Where nested HAPLOID chains ended up. Reported because it is how the population is gated: all
// 44,139 "no strand" sites across chr20, chr6, chr17 and chrX were chrX's and none were autosomal,
// so a bug confined to one of these buckets is invisible to any autosome-only check.
static std::atomic<size_t> g_nest_strand(0), g_nest_one_hap(0);
// The two ways a nested haploid chain under a DIPLOID parent ends up on no strand. Counted apart
// because they are different facts with different right answers, and both are empty on every contig
// measured -- so if either ever fires, which one it is decides what to do about it.
static std::atomic<size_t> g_nest_both(0), g_nest_unreadable(0);

/// Which of its parent's two strands a nested haploid chain sits on, or -1.
///
/// `carrying` is the parent candidate traversal that crosses the chain, from `relate_to_parent`.
/// Both guards are load-bearing.
///
/// A NESTED HAPLOID parent occupies ONE haplotype, so everything inside it is on that haplotype.
/// Its own `nested_strand` is the answer, and the identity match cannot find it: such a parent has
/// trav_first == trav_second and ploidy 1, so the match can only ever return strand 0. On chr20 that
/// put all 448 depth->=2 sites under a strand-1 parent on the wrong haplotype.
///
/// A DIPLOID parent names the strand by which of its two settled traversals carries the chain. The
/// `ploidy == 2` guard covers BOTH assignments, not just the second: a strand is a claim that the
/// locus has two haplotypes and this allele sits on one of them, which holds only where the parent
/// has two copies. Without it a haploid top-level parent -- chrX outside the pseudoautosomal
/// regions, or chrY -- gives its children strand 0 and they render "a|.", asserting a second strand
/// that carries nothing. chrX's haploid interior carried 8,056 such records against 1,965, and chrX
/// was the one contig whose F1 fell (0.94939 -> 0.93643) while all 22 autosomes rose.
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
    // Genotypes the model settled on that the emitted record carries no ALT for. Reachable only
    // because the collector's allele space is the genotyper's rather than the VCF's, so it is
    // counted and reported rather than assumed away.
    // Records whose phase names the called pair because the settled pair had no ALT. The phase is
    // still real -- the block and the strand order come from the panel either way -- but the alleles
    // it names are the line's rather than the model's, so the count belongs in the report.
    if (!model.active() || entries.empty()) {
        return moved;
    }

    // The phase an earlier generation settled on, by record key, so a clamped site can be pinned to
    // the pair the VCF already carries. Read before the chain loop appends this generation's own.
    // Everything a child needs of its settled parent: the haplotype pair, and -- for a nested
    // HAPLOID child, which sits on one of the two -- what `nested_strand_of` needs to say which.
    struct PinnedPhase {
        size_t first;
        size_t second;
        int trav_first;
        int trav_second;
        size_t ploidy;
        int nested_strand;
        bool order_arbitrary;
    };
    // Rebuilt from the whole accumulated phasing on every call, which is the O(generations * sites)
    // term this function has left: 49.6 ms of chr20's 50.3 ms preamble, against 0.7 ms for the
    // entry scans now that `by_generation` indexes them.
    //
    // Left alone deliberately. An append cursor into `phasing_out` would make it incremental, but
    // a `last` call SORTS that vector in place and the caller may call again afterwards when a
    // deeper chain grows the generation bound -- so the cursor would be silently invalidated by a
    // reorder, to save 50 ms of a 170 s run. It is also the only thing standing between this and a
    // per-subtree scoping; that is recorded at the generation loop, which is where the decision to
    // stay level-order lives.
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
        // Whether naming a panel haplotype for this site is legitimate at all.
        //
        // Under a HAPLOID parent there is no strand to CHOOSE because there is only one -- the child
        // sits on it and the haplotype is nameable. That is chrX's ordinary case.
        //
        // Under a DIPLOID parent, `nested_strand_of` returns -1 for three situations that
        // `relate_to_parent` spells differently, and only two of them can arrive here:
        //
        //   copies == 0, neither settled traversal carries the chain. NOT reachable: the barrier
        //     retracts the whole subtree (`drop_subtree`) before this generation is grouped, so a
        //     live entry always has a parent that carries it. An earlier version of this comment
        //     claimed this was the case being guarded, and it was the one case that cannot occur.
        //   copies == 2, BOTH carry it. Reachable only when the chain kept a ploidy of 1 that its
        //     parent's settled pair contradicts -- the barrier could not re-render it, so the
        //     revision to ploidy 2 was skipped (387 such chains on chr20, 1,812 on chrX, none of
        //     them this).
        //   the settled pair could not be read at all.
        //
        // All of them name nothing, which for copies == 2 is CONSERVATIVE rather than right: both
        // haplotypes carry the chain, so both could be named, and the per-strand pass did exactly
        // that. Not reproduced here because the population is empty -- 0 on chr20, chr6, chr17 and
        // chrX under the old pass's own counter, and 0 on chr20 and chrX under this one -- and
        // building a rendering for a case nobody has observed is how the last three of these got
        // their behaviour wrong. The counters are split so that a nonzero says which.
        bool nameable = true;
    };
    unordered_map<size_t, NestedPlacement> unified_strand;

    // This generation's entries, in append order -- the same indices the full scans selected, in
    // the same order, so substituting the index for a scan cannot reorder anything.
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
    // in a reference-first graph but is not guaranteed to be it, and the transition probabilities
    // are computed from the gaps -- so trusting the arrival order would silently feed the model
    // the wrong distances.
    //
    // Only at the TOP LEVEL. Below it there are no contig runs to build: every live site is nested,
    // and a nested site is grouped with its parent instead.
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
    // Per chain: the message to condition it on (empty = none), and the phase set it belongs to
    // (SIZE_MAX = take it from the chain's own first site, as a contig chain does).
    //
    // `deltas` owns what `chain_context` points at, so it is declared HERE and not beside the loop
    // that fills it: the pointers are read by the decode below, which is outside that loop's scope.
    // Declared inside it, the deque was destroyed while `chain_context` still pointed into it -- a
    // use-after-free that reads correctly almost every time, because freed heap usually still holds
    // the bytes, and would have been a genuinely intermittent wrong answer. A deque and not a
    // vector, so that pushing another cannot move the ones already pointed at.
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
        // Every entry here is top level -- `by_contig` is built at generation 0 only, and a
        // top-level snarl has no parent to be nested in -- so the whole contig is chainable.
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
    // Each nested CHAIN is decoded on its own, conditioned on its parent -- never pooled with its
    // parent's other chains.
    //
    // This is the whole of the "is linkage BETWEEN nested chains worth anything?" question. Within a
    // chain the snarls stay linked, so recombination inside a chain is still modelled; between two
    // chains under one parent there is no transition at all, and each is answered from the parent's
    // context alone.
    //
    // MEASURED, on three contigs, by decoding each chain alone against pooling a parent's children:
    // chr20 every small-variant class identical to five decimals with matching TP/FP/FN and SV
    // -3.3e-4; chr6 +1.8e-6 ALL and SV identical; chr17 -1.3e-5 ALL and SV +8.4e-4. Every movement
    // is one call and the signs disagree across contigs, and chr17 carries 273 multi-crossing
    // children where the other two carry none -- so the case with the best chance of mattering did
    // not. Inter-chain linkage buys nothing.
    //
    // That is what licenses deleting inter-chain ORDER and inter-chain DISTANCE: both are only
    // observable where a transition crosses between chains, and none does.
    if (generation > 0 && !pinned_phase.empty()) {
        unordered_map<size_t, size_t> index_of_key;
        index_of_key.reserve(entries.size() * 2);
        for (size_t i = 0; i < entries.size(); ++i) {
            if (!entries[i].retracted) {
                index_of_key[entries[i].record_key] = i;
            }
        }
        // Keyed by (parent, chain-within-parent). The chain half is 0 unless per_chain_groups, so
        // the default is exactly the previous one-group-per-parent behaviour.
        map<tuple<size_t, size_t, size_t>, vector<size_t>> by_parent;
        // What the settled parent implies about a child, computed HERE from the parent's settled
        // pair and the child's crossing mask rather than read from fields a separate pass wrote.
        // This replaced `Entry::parent_trav`, which cached the same fact and could go stale.
        auto relate = [&](const Entry& child, const Entry& parent) {
            const int ta = traversal_of(trav_arena, parent.trav_offset, parent.num_alleles,
                                        parent.final_i);
            const int tb = parent.ploidy == 2
                               ? traversal_of(trav_arena, parent.trav_offset, parent.num_alleles,
                                              parent.final_j)
                               : -1;
            return relate_to_parent(child.parent_crossing, ta, tb);
        };

        // (parent, chain), where the chain half is its boundary pair from the graph. A snarl the
        // decomposition puts in no chain becomes its own group rather than being pooled with every
        // other such snarl under this parent.
        auto group_key = [&](const Entry& e) {
            // The chain's own boundary pair, from the graph. Measured against the alignment-rank
            // key it replaces, on three contigs: 6,481 groups partitioned identically and 37
            // differing -- chains the alignment could not identify, which it isolated as singletons
            // and which the graph groups with their siblings. Small-variant F1 unchanged or better,
            // SV unchanged on chr6 and chr17.
            const size_t chain = e.chain_key != 0
                                     ? e.chain_key
                                     : numeric_limits<size_t>::max() - e.record_key;
            // PLOIDY IS PART OF THE KEY, and it has to be. A group is decoded as one chain at one
            // ploidy -- `chain_ploidy` below reads it off `indices.front()` and applies it to every
            // member -- so a group holding both ploidies decodes half of them in the wrong state
            // space. That is not a wrong answer, it is a heap overflow: a ploidy-1 site's
            // `genotype_ln_likelihood` has `num_alleles` entries, the diploid marginals index the
            // parallel `known`/`wild` accumulators at `genotype_index(ai, bi)`, and that runs to
            // n(n+1)/2 - 1.
            //
            // Siblings under one parent in one chain need not share a ploidy: ploidy here is how
            // many of the parent's settled alleles cross the child, so one child crossed by both
            // haplotypes sits beside one crossed by a single haplotype. On the default path the
            // population happened to be homogeneous and this never fired. Admitting the chains the
            // reference does not cross mixes it -- their copies split 249/572/699 across 0/1/2 on
            // chr20 -- and the corruption is then detected by malloc, deterministically, an
            // allocation or two later and nowhere near the write.
            return make_tuple(e.parent_record_key, chain, e.ploidy);
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
        // Every live site of this generation, straight from the entries, grouped with its parent.
        //
        // This used to walk the contig runs and veto a WHOLE run if any one of its sites could not
        // be grouped, sending that run to the contig chain. Two things were wrong with it. The scope
        // was wrong -- a site that cannot be grouped is better decoded alone than chained to
        // unrelated sites under different parents, which is what the contig chain does and what the
        // comment there says costs about 1,100 severed links a generation. And it has never fired:
        // 0 sites with no parent key, 0 whose parent has no live entry, 0 on a parent/child ploidy
        // difference, on chr20, chr6 and chr17.
        //
        // Group membership is order-independent because every group is sorted afterwards, on
        // (position, record key), which is total.
        vector<vector<size_t>> ungrouped;
        for (uint32_t idx : this_generation) {
            const Entry& e = entries[idx];
            if (e.retracted) {
                continue;
            }
            auto par = index_of_key.find(e.parent_record_key);
            if (e.parent_record_key == 0) {
                ++g_grp_no_parent;
            } else if (par == index_of_key.end()) {
                ++g_grp_no_entry;
            } else {
                by_parent[group_key(e)].push_back(idx);
                continue;
            }
            // Decoded alone rather than dropped. What the veto was protecting is that a site left
            // out of every group is never decoded, never settled and never phased.
            ++g_grp_vetoed;
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
                // The parent joins the group IFF its ploidy matches the children's.
                //
                // That is the pre-existing rule -- a chain is a maximal run of one ploidy -- applied
                // uniformly instead of as a special case, and it makes chrX ordinary. On chrX's
                // haploid interior a nested chain is a haploid child of a HAPLOID parent, so the
                // ploidies match and the parent is prepended exactly as it is for a diploid parent
                // and diploid children. Only the mixed case -- a haploid child of a diploid parent,
                // which is the autosomal nested-haploid population -- leaves the parent out, and
                // then the parent reaches the group as its entering message instead, which is all it
                // was contributing by being in it.
                // Every member shares this now, because the key says so.
                const size_t group_ploidy =
                    kv.second.empty() ? entries[pidx].ploidy : std::get<2>(kv.first);
                const bool parent_in_group = entries[pidx].ploidy == group_ploidy;
                vector<size_t> group;
                if (parent_in_group) {
                    group.push_back(pidx);
                }
                // A TOTAL key, not a chain of skippable comparisons.
                //
                // "Compare on this key only if BOTH operands have it, else fall through" is not a
                // strict weak ordering, and std::sort on one is undefined behaviour rather than
                // merely a different order: with absence spelled as -1, A(-, 2), B(1, 3), C(2, 1)
                // gives A<B, B<C and C<A. A brute force over small triples found 204 such cycles
                // here. Both partial keys are gone now -- measured to order nothing anything reads
                // -- and what is left is total for every entry, which is why there is no longer a
                // per-group decision about which key to use.
                sort(kv.second.begin(), kv.second.end(), [&](size_t a, size_t c) {
                    const Entry& ea = entries[a];
                    const Entry& ec = entries[c];
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
                // A DELTA at the pair the parent SETTLED on, built where it is consumed.
                //
                // This was the parent's posterior over ordered pairs, harvested from its own decode
                // with an alpha/beta pass and a sparse mask, and held in a map at m*m doubles per
                // parent -- 16.4 MB on chr20 at generation 0, 28.1 MB on chr6. Measured against it
                // on both contigs: chr6 byte-identical, chr20 identical on every scored metric
                // including SV. A decided parent has no residual uncertainty to pass down, and that
                // turns out not to be an approximation.
                //
                // It also covers MORE than the harvest did. A parent reached only under a settled
                // genotype was never masked for harvesting, so its children were decoded with no
                // context at all -- 51 groups on chr20, 64 on chr6. Those parents have a PhaseCall,
                // so they get a message now, and the "no stored context" population stops existing.
                //
                // Built rather than stored: one m*m buffer per group that has a settled parent,
                // against one per parent held for a whole generation. `deltas` owns them and is a
                // deque so that pushing another cannot move what is already pointed at.
                //
                // WILDCARD is `(size_t)-1` outside the model and state `n_haplotypes` inside it. A
                // parent whose strand the panel cannot explain has still settled -- on the wildcard
                // -- and that is the pair the delta belongs on; left untranslated it wraps and can
                // land on a pair nothing chose.
                const size_t m = n_haplotypes + 1;
                auto state_of = [&](size_t h) {
                    return h == LinkageModel::WILDCARD ? n_haplotypes : h;
                };
                auto pin = pinned_phase.find(std::get<0>(kv.first));
                if (pin == pinned_phase.end()) {
                    gctx.push_back(nullptr);
                } else if (group_ploidy == 1) {
                    // ONE haplotype, not a pair: the chain sits on one of the parent's two strands,
                    // and which one is `nested_strand_of`. A haploid parent records its single
                    // haplotype in `hap_first` and `hap_second` is meaningless at ploidy 1, so the
                    // pair is indexed only where there is a pair to index.
                    const Entry& child = entries[kv.second.front()];
                    const int carrying = relate(child, entries[pidx]).carrying_trav;
                    const int strand = nested_strand_of(carrying, pin->second.ploidy,
                                                        pin->second.trav_first,
                                                        pin->second.trav_second,
                                                        pin->second.nested_strand);
                    const size_t hap = pin->second.ploidy == 1
                                           ? pin->second.first
                                           : (strand == 1 ? pin->second.second : pin->second.first);
                    // DECLINE where the haplotype does not traverse the child: a named haplotype
                    // carrying no allele here is not evidence about which haplotype the sample is
                    // on, and all the mass would land on the latent-allele escape.
                    // A strand is only needed to CHOOSE between two haplotypes. A haploid parent
                    // has one, so it needs no strand -- and a top-level haploid parent has none to
                    // give: `nested_strand` is a nested site's field, and chrX's interior sites are
                    // not nested. Requiring strand >= 0 there declined the whole contig.
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
                        g_nest_strand += kv.second.size();
                    } else if (have_hap) {
                        g_nest_one_hap += kv.second.size();
                    } else if (carrying == -2) {
                        g_nest_both += kv.second.size();
                    } else {
                        g_nest_unreadable += kv.second.size();
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
            // The chains that could not be grouped, carried through unchanged. Without this they
            // would be replaced by groups that do not contain them and their sites would never be
            // decoded at all -- which is the failure the global flag existed to prevent, and the
            // reason relaxing it needs this line and not just a looser test.
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

        // Every entry in a chain shares a ploidy by construction above -- ASSERTED, not assumed.
        // When it was only assumed, a mixed group decoded its minority members in the wrong state
        // space and overflowed the `known`/`wild` accumulators, which malloc reported as heap
        // corruption an allocation later and several frames away. Cheap: one pass over a group.
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
        // ONE context for both ploidies. A group's message is a delta built from its parent's
        // PhaseCall: over ordered PAIRS for a diploid group, over SINGLE haplotypes for a haploid
        // one. Each decode ignores a message of the wrong size for its state space, so nothing here
        // has to dispatch on ploidy to build or pass it -- and `posteriors` and `phasing` are told
        // the ploidy and pick the state space themselves.
        const vector<double>* ctx =
            chain_i < chain_context.size() ? chain_context[chain_i] : nullptr;
        if (chain_ploidy == 1) {
            posteriors = model.posteriors(sites, 1, ctx);
        } else if (ctx != nullptr) {
            // A child group, conditioned on its parent's settled pair. One decode over the group,
            // not over the contig.
            //
            // It used to harvest as it went, so that the next generation would find a message at
            // this group's own parents. Nothing harvests now: the message is a delta built from the
            // PhaseCall wherever it is needed, so there is nothing to carry forward and no way for
            // grouping at one depth to silently disable it at the next.
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
        // Haploid chains take the single-strand path inside `phasing`: one haplotype per site, no
        // phase to infer, but the mosaic is exactly as meaningful -- it is the whole answer on chrY.
        // The SAME message the posteriors above were given, which is the whole reason there is one
        // entry point rather than two.
        vector<LinkageModel::Phase> phase =
            model.phasing(sites, final_genotype, chain_ploidy, ctx);
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
            pc.depth = e.generation;
            // The strand a nested haploid chain sits on. A nested-specific fact -- the chain
            // occupies ONE of its parent's two haplotypes -- and nothing in a chain decode can
            // derive it, so it is carried from where the parent's settled pair was in hand.
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
            pc.emitted = e.emitted;
            pc.record_key = e.record_key;
            pc.contig = contig_names[e.contig];
            pc.position = e.position;
            // WHICH SLOT, not just which haplotype. A haploid decode returns one haplotype and
            // the render used to drop it in `hap_first` unconditionally -- but the mosaic reads the
            // slot NAMED BY `nested_strand` and treats the other as the empty one, so a chain on
            // the parent's second strand came out carried on the strand it is not on and
            // unexplained on the strand it is. About half of the ~6,800 strand-placed nested sites
            // per contig, and invisible to F1 because it is a mosaic-only field.
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

    if (generation > 0 && (g_pin_applied.load() + g_pin_declined.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": pins -- "
                  << g_pin_applied.load() << " applied, " << g_pin_declined.load()
                  << " REFUSED (the pinned pair cannot spell the constrained genotype, so that"
                  << " site's orientation is free); groups whose parent was pinnable: "
                  << g_group_parent_pinned.load() << ", not pinnable: "
                  << g_group_parent_unpinned.load() << std::endl;
    }

    // Only when something actually declined. These have never fired on any contig measured, so
    // printing them every generation is a line that says "zero" forever and trains the reader to
    // skip it -- which is the opposite of what a counter kept as an alarm is for.
    if (generation > 0
        && (g_grp_no_parent.load() + g_grp_no_entry.load() + g_grp_vetoed.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": grouping declines so far -- "
                  << g_grp_no_parent.load() << " sites with no parent key, "
                  << g_grp_no_entry.load() << " whose parent has no live entry; "
                  << g_grp_vetoed.load() << " chains kept ungrouped in total" << std::endl;
    }

    if (grouped_groups > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage generation " << generation << ": decoded " << grouped_sites
                  << " sites in " << grouped_groups << " per-parent groups instead of the contig"
                  << " chain" << std::endl;
    }

    if (generation > 0
        && (g_nest_strand.load() + g_nest_one_hap.load() + g_nest_both.load()
            + g_nest_unreadable.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] nested strands: " << g_nest_strand.load()
                  << " on one of a diploid parent's two strands, " << g_nest_one_hap.load()
                  << " on a haploid parent's single haplotype (no strand to choose), "
                  << g_nest_both.load() << " carried on both parent strands, "
                  << g_nest_unreadable.load()
                  << " whose parent's settled pair could not be read -- the last two name no"
                  << " haplotype" << std::endl;
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
