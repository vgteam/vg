#include "linkage_model.hpp"

#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>

namespace vg {


/// Return the distance in bp from a site at `prev` to the next site of its chain, at `next`, or
/// the largest `size_t` where no distance is known.
///
/// Two positions measured the same way can be differenced. They are both reference positions, or
/// both offsets along the same parent allele for unpositioned sites. A positioned parent and the
/// unpositioned first child of its group can be differenced too, since the child's position is
/// its parent's start plus its offset along the parent's allele.
///
/// An unknown distance makes rho 1, so the transition is uniform. A known distance is at least 1,
/// since `switch_probability` reads a gap of 0 as 1 and a negative one would wrap.
static inline size_t position_gap(size_t prev, bool prev_unpositioned, bool prev_is_parent,
                                  size_t next, bool next_unpositioned) {
    if (prev_unpositioned != next_unpositioned && !(prev_is_parent && next_unpositioned)) {
        return numeric_limits<size_t>::max();
    }
    return next > prev ? next - prev : 1;
}

/// Return the `position_gap` from site `prev` to the next site of its chain, once for each strand.
///
/// `transition_apply` takes a switch probability per strand, though the two gaps are always
/// equal.
static inline std::pair<size_t, size_t> site_gap(const LinkageModel::Site& prev,
                                                const LinkageModel::Site& next) {
    const size_t d = position_gap(prev.position, prev.unpositioned, prev.group_parent,
                                  next.position, next.unpositioned);
    return {d, d};
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
        // Transitions are uniform, so each site's posterior is its emission.
        return 1.0;
    }
    return min(max(pow(rho, params.weight), 1e-12), 1.0);
}

/// Return the allele panel haplotype `h` carries at a site, or -1 if it carries none there or
/// names an allele out of range.
static inline int allele_at(const LinkageModel::Site& site, size_t h, size_t n_hap) {
    if (h >= n_hap || h >= site.haplotype_allele.size()) {
        return -1;
    }
    int allele = site.haplotype_allele[h];
    return (allele >= 0 && (size_t)allele < site.num_alleles) ? allele : -1;
}

/// One strand's distribution over a site's `n` alleles, as a uniform part plus a point mass.
///
/// P(allele i) = base + keep * [i == allele]. Under mutation (`Params::mutation`), a strand that copies a haplotype carrying `allele` keeps it
/// with probability 1 - mutation, and takes each other allele with probability mutation / (n - 1).
/// So `base` is mutation / (n - 1), and `keep` the rest. A strand whose haplotype does not pass
/// through the site has no allele, and takes each allele with probability escape / n.
///
/// Writing the distribution as a uniform part plus a point mass makes a pair's emission a sum of
/// four terms, each a likelihood summed over a row or over the whole table, so it costs O(1) per
/// pair once the row sums are known (see `MutationTerms`).
struct StrandEmission {
    /// The probability of each allele from the uniform part.
    double base = 0.0;
    /// The extra probability of `allele`.
    double keep = 0.0;
    /// The allele of the haplotype the strand copies, or -1 for none.
    int allele = -1;

    /// Return the probability that the strand carries allele `i`.
    double p(size_t i) const {
        return base + (allele >= 0 && (size_t)allele == i ? keep : 0.0);
    }
};

/// Return the distribution of a strand whose haplotype carries `allele`, or -1 for none, at a site
/// with `n` alleles.
static StrandEmission strand_emission(int allele, size_t n, double mutation, double escape) {
    StrandEmission s;
    if (n == 0) {
        return s;
    }
    if (allele < 0) {
        s.base = escape / (double)n;
        return s;
    }
    s.allele = allele;
    if (n == 1) {
        // No other allele to mutate to.
        s.keep = 1.0;
        return s;
    }
    s.base = mutation / (double)(n - 1);
    s.keep = 1.0 - mutation - s.base;
    return s;
}

/// The per-class terms of a diploid site's emission under mutation.
///
/// Strands are classed by their haplotype's allele, with class 0 for a haplotype that does not
/// pass through the site and class x + 1 for allele x. A pair's emission depends only on its two
/// classes.
struct MutationTerms {
    /// The site's number of alleles.
    size_t n = 0;
    /// The allele distribution of a strand of each class.
    vector<StrandEmission> strand;
    /// At x, the sum over j of the likelihood of genotype (x, j). The sum is over ordered pairs,
    /// since each strand picks its allele separately.
    vector<double> row;
    /// The sum of `row`.
    double total = 0.0;
    /// P(reads | the pair of classes (cx, cy)), at cx * (n + 1) + cy.
    vector<double> table;

    /// Return the emission of a pair of strands whose haplotypes carry alleles `ai` and `bi`, or
    /// -1 for none.
    double emission(int ai, int bi) const {
        return table[(size_t)(ai + 1) * (n + 1) + (size_t)(bi + 1)];
    }
};

/// Return the terms for one diploid site with `n` alleles, from its genotype likelihoods
/// `per_genotype`, scaled so that the best is 1.
static MutationTerms mutation_terms(size_t n, const vector<double>& per_genotype,
                                    double mutation, double escape) {
    MutationTerms mt;
    mt.n = n;
    mt.strand.resize(n + 1);
    mt.strand[0] = strand_emission(-1, n, mutation, escape);
    for (size_t x = 0; x < n; ++x) {
        mt.strand[x + 1] = strand_emission((int)x, n, mutation, escape);
    }
    mt.row.assign(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            mt.row[i] += per_genotype[LinkageModel::genotype_index(i, j)];
        }
        mt.total += mt.row[i];
    }
    // Summing P(i | x) P(j | y) L(i, j) over ordered (i, j), with
    // P(i | x) = base_x + keep_x [i == x], gives four terms. They are base_x base_y total,
    // base_x keep_y row[y], keep_x base_y row[x], and keep_x keep_y L(x, y). A strand with no
    // allele has keep 0.
    mt.table.assign((n + 1) * (n + 1), 0.0);
    for (size_t cx = 0; cx <= n; ++cx) {
        const StrandEmission& sx = mt.strand[cx];
        for (size_t cy = 0; cy <= n; ++cy) {
            const StrandEmission& sy = mt.strand[cy];
            double v = sx.base * sy.base * mt.total;
            if (sy.allele >= 0) {
                v += sx.base * sy.keep * mt.row[(size_t)sy.allele];
            }
            if (sx.allele >= 0) {
                v += sx.keep * sy.base * mt.row[(size_t)sx.allele];
            }
            if (sx.allele >= 0 && sy.allele >= 0) {
                v += sx.keep * sy.keep
                     * per_genotype[LinkageModel::genotype_index((size_t)sx.allele,
                                                                 (size_t)sy.allele)];
            }
            mt.table[cx * (n + 1) + cy] = v;
        }
    }
    return mt;
}

/// Split a diploid genotype index into its alleles, with `i` <= `j`, inverting
/// `LinkageModel::genotype_index`.
static void genotype_alleles(size_t genotype, size_t& i, size_t& j) {
    j = 0;
    while (LinkageModel::genotype_index(0, j + 1) <= genotype) {
        ++j;
    }
    i = genotype - (j * (j + 1) / 2);
}

/// Return the genotype index of an allele pair at a ploidy: the allele itself at ploidy 1, and
/// `LinkageModel::genotype_index` at ploidy 2.
static size_t genotype_of(size_t ploidy, size_t i, size_t j) {
    return ploidy == 1 ? i : LinkageModel::genotype_index(i, j);
}

/// Rescale a vector of non-negative messages to sum to 1, leaving an all-zero one as it is.
static void rescale(vector<double>& v) {
    double sum = 0.0;
    for (double x : v) {
        sum += x;
    }
    if (sum <= 0.0) {
        sum = 1.0;
    }
    for (double& x : v) {
        x /= sum;
    }
}

/// The mask is indexed by candidate traversal, so it is tested against traversals, never against
/// compact allele indices.
LinkageCollector::Relation LinkageCollector::relate_to_parent(uint64_t crossing, int ta, int tb) {
    LinkageCollector::Relation r;
    if (ta < 0 || crossing == 0) {
        // Nothing was chosen, or descent could not compute the mask.
        return r;
    }
    r.known = true;
    const bool first = ta < 64 && ((crossing >> ta) & 1);
    const bool second = tb >= 0 && tb < 64 && ((crossing >> tb) & 1);
    r.copies = (uint8_t)((int)first + (int)second);
    if (r.copies == 2) {
        r.carrying_trav = -2;
    } else if (r.copies == 1) {
        r.carrying_trav = first ? ta : tb;
    }
    return r;
}

/// Fill `e` with the relative P(reads | state) of every ordered pair of panel haplotypes, at
/// a * n_hap + b.
///
/// Each strand carries an allele drawn from its `StrandEmission`, and the emission sums the
/// genotype likelihood over them. `per_genotype` gets the site's genotype likelihoods, scaled so
/// that the best is 1.
static void build_emission(const LinkageModel::Site& site, size_t n_hap, double escape,
                           double mutation, vector<double>& e, vector<double>& per_genotype) {
    const size_t m = n_hap;

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

    const MutationTerms mt = mutation_terms(site.num_alleles, per_genotype, mutation, escape);
    e.assign(m * m, 0.0);
    for (size_t a = 0; a < m; ++a) {
        int ai = allele_at(site, a, n_hap);
        for (size_t b = 0; b < m; ++b) {
            e[a * m + b] = mt.emission(ai, allele_at(site, b, n_hap));
        }
    }
}

/// Each strand's transition is T = (1-rho) I + (rho/m) 1. So the sum over previous ordered pairs
/// collapses to four terms, and the step is O(m^2) rather than O(m^4).
void transition_apply(const vector<double>& in, size_t m,
                      double rho_a, double rho_b, vector<double>& out) {
    // One switch probability per strand, since the strands switch independently.
    const double stay_a = 1.0 - rho_a;
    const double stay_b = 1.0 - rho_b;
    const double jump_a = rho_a / (double)m;
    const double jump_b = rho_b / (double)m;
    vector row(m, 0.0);
    vector col(m, 0.0);
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

/// The best two values along one axis of a table, with their indices.
///
/// The leave-one-out maximum of the axis at index i is `second` where i is `arg`, and `best`
/// elsewhere.
struct Top2 {
    double best = NEG_INF;
    double second = NEG_INF;
    size_t arg = 0;
    size_t arg2 = 0;
};

/// Return the best two of the `n` values at `v`, `stride` apart. A stride of 1 walks a row of a
/// table, and a stride of the row length walks a column.
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

/// One Li-Stephens max-product step over ordered pairs, from the previous site's log scores.
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
/// The scores are logs, unlike the forward pass: sum-product needs rescaling per site to avoid
/// underflow, max-product does not, and in logs the stay-or-jump choice is a comparison of sums.
class MaxProductStep {
public:
    /// Prepare a step from the previous site's scores `in`, with `m` haplotypes per strand and a
    /// switch probability per strand. Each strand stays on its haplotype with probability
    /// 1 - rho + rho / m, and moves to each other one with probability rho / m.
    MaxProductStep(const vector<double>& in, size_t m, double rho_a, double rho_b)
        : in(in), m(m) {
        const double stay_a = 1.0 - rho_a + rho_a / (double)m;
        const double stay_b = 1.0 - rho_b + rho_b / (double)m;
        const double jump_a = rho_a / (double)m;
        const double jump_b = rho_b / (double)m;
        S_0 = stay_a > 0.0 ? log(stay_a) : NEG_INF;
        S_1 = stay_b > 0.0 ? log(stay_b) : NEG_INF;
        J_0 = jump_a > 0.0 ? log(jump_a) : NEG_INF;
        J_1 = jump_b > 0.0 ? log(jump_b) : NEG_INF;

        rows.resize(m);
        cols.resize(m);
        for (size_t a = 0; a < m; ++a) {
            rows[a] = top2_of(&in[a * m], m, 1);
        }
        for (size_t b = 0; b < m; ++b) {
            cols[b] = top2_of(&in[b], m, m);
        }
        row_excl.resize(m * m);
        row_excl_arg.resize(m * m);
        for (size_t a = 0; a < m; ++a) {
            for (size_t bp = 0; bp < m; ++bp) {
                const bool hit = (rows[a].arg == bp);
                row_excl[a * m + bp] = hit ? rows[a].second : rows[a].best;
                row_excl_arg[a * m + bp] = (uint16_t)(hit ? rows[a].arg2 : rows[a].arg);
            }
        }
    }

    /// Fill `out` with the next site's scores under its `emission`, and the backpointers with
    /// each state's best predecessor.
    void apply(const vector<double>& emission, vector<double>& out, vector<uint16_t>& back_a,
               vector<uint16_t>& back_b) const {
        out.assign(m * m, NEG_INF);
        back_a.assign(m * m, 0);
        back_b.assign(m * m, 0);
        for (size_t bp = 0; bp < m; ++bp) {
            // The both-moved maxima for every a' arriving at b' come from one top-2.
            const Top2 both = top2_of(&row_excl[bp], m, m);
            for (size_t ap = 0; ap < m; ++ap) {
                const double e = emission[ap * m + bp];
                if (!(e > 0.0)) {
                    // An impossible state stays impossible. This is also how a constraint is
                    // applied: the caller zeroes the emission of every state that does not spell
                    // the required genotype.
                    continue;
                }
                size_t ba = 0;
                size_t bb = 0;
                const double best = best_predecessor(ap, bp, both, ba, bb);
                if (best > NEG_INF) {
                    out[ap * m + bp] = best + log(e);
                    back_a[ap * m + bp] = (uint16_t)ba;
                    back_b[ap * m + bp] = (uint16_t)bb;
                }
            }
        }
    }

private:
    /// Return the best score of a predecessor of state (ap, bp) plus its transition, and set
    /// (ba, bb) to that predecessor. Earlier cases win ties.
    double best_predecessor(size_t ap, size_t bp, const Top2& both, size_t& ba,
                            size_t& bb) const {
        double best = NEG_INF;
        // Add the transition's two strand terms one at a time, in strand order.
        auto consider = [&](double score, double strand_0, double strand_1, size_t a, size_t b) {
            if (score > NEG_INF) {
                const double v = score + strand_0 + strand_1;
                if (v > best) {
                    best = v;
                    ba = a;
                    bb = b;
                }
            }
        };
        // Both strands stayed.
        consider(in[ap * m + bp], S_0, S_1, ap, bp);
        // Strand 0 stayed.
        consider(row_excl[ap * m + bp], S_0, J_1, ap, row_excl_arg[ap * m + bp]);
        // Strand 1 stayed.
        const bool col_hit = (cols[bp].arg == ap);
        consider(col_hit ? cols[bp].second : cols[bp].best, J_0, S_1,
                 col_hit ? cols[bp].arg2 : cols[bp].arg, bp);
        // Both strands moved.
        const bool both_hit = (both.arg == ap);
        const size_t moved_a = both_hit ? both.arg2 : both.arg;
        consider(both_hit ? both.second : both.best, J_0, J_1, moved_a,
                 row_excl_arg[moved_a * m + bp]);
        return best;
    }

    /// The previous site's scores, over ordered pairs at a * m + b.
    const vector<double>& in;
    /// Haplotypes per strand.
    size_t m;
    /// Log probabilities of staying (S) and of moving to one particular other haplotype (J), for
    /// strand 0 and strand 1.
    double S_0;
    double S_1;
    double J_0;
    double J_1;
    /// The best two scores in each row a and each column b of `in`.
    vector<Top2> rows;
    vector<Top2> cols;
    /// At a * m + b', the maximum over b != b' of `in`(a, b).
    vector<double> row_excl;
    /// At a * m + b', the b that gives `row_excl` there.
    vector<uint16_t> row_excl_arg;
};

}   // anonymous namespace

/// Return the most haplotypes any of sites [from, to) has an allele entry for, which is the number
/// of states per strand that a window over them decodes with.
static size_t panel_width(const vector<LinkageModel::Site>& sites, size_t from, size_t to) {
    size_t n_hap = 0;
    for (size_t t = from; t < to; ++t) {
        n_hap = max(n_hap, sites[t].haplotype_allele.size());
    }
    return n_hap;
}

/// Count how many of the first `n_hap` panel haplotypes carry each allele of a site.
static vector<size_t> allele_carriers(const LinkageModel::Site& site, size_t n_hap) {
    vector<size_t> carriers(site.num_alleles, 0);
    for (size_t h = 0; h < n_hap && h < site.haplotype_allele.size(); ++h) {
        const int allele = site.haplotype_allele[h];
        if (allele >= 0 && (size_t)allele < site.num_alleles) {
            carriers[(size_t)allele] += 1;
        }
    }
    return carriers;
}

/// A ledger of one diploid site's posterior mass per genotype, from which the site's posterior is
/// read once every state's mass is added.
///
/// Each state's mass g = alpha * beta is shared among the genotypes its strands could carry, in
/// proportion to P(i | a) P(j | b) L(i, j). The four terms of that product (see `mutation_terms`)
/// are kept apart, because the frequency prior treats them differently:
/// - Both strands keep their haplotypes' alleles. The mass goes to the genotype the pair spells,
///   and counts towards its multiplicity.
/// - One strand keeps allele k and the other takes any allele. The mass is spread over (k, j) by
///   likelihood, and rescaled by the number of k's carriers.
/// - Both take any allele. The mass is spread over every genotype by likelihood, unrescaled.
///
/// The uniform part of a strand's distribution includes its own allele, so a little of the mass
/// of a pair that spells the genotype is counted in the second or third term rather than the first.
/// That is at most mutation / (n - 1) of it.
class DiploidSiteMass {
public:
    /// Start an empty ledger for a site with genotype likelihoods `per_genotype`, scaled so that
    /// the best is 1, decoded over `n_hap` haplotypes per strand.
    DiploidSiteMass(const LinkageModel::Site& site, const vector<double>& per_genotype,
                    size_t n_hap, double mutation, double escape)
        : site(site), per_genotype(per_genotype), n_hap(n_hap),
          mt(mutation_terms(site.num_alleles, per_genotype, mutation, escape)),
          known(site.genotype_ln_likelihood.size(), 0.0),
          multiplicity(site.genotype_ln_likelihood.size(), 0),
          kept(site.num_alleles, 0.0) {}

    /// Add the mass alpha * beta of every ordered pair of haplotypes.
    void add_states(const vector<double>& alpha, const vector<double>& beta) {
        const size_t n = site.num_alleles;
        if (n == 0) {
            return;
        }
        const size_t m = n_hap;
        for (size_t a = 0; a < m; ++a) {
            const int ai = allele_at(site, a, n_hap);
            for (size_t b = 0; b < m; ++b) {
                const double g = alpha[a * m + b] * beta[a * m + b];
                if (g > 0.0) {
                    add_state(ai, allele_at(site, b, n_hap), g);
                }
            }
        }
    }

    /// Write the site's posterior over genotypes to `post`, with the frequency exponent
    /// `freq_prior` applied, or leave `post` empty where the site has no mass.
    void posterior(double freq_prior, vector<double>& post) const {
        const vector<double> wild = spread_wild_mass(freq_prior);
        post.assign(site.genotype_ln_likelihood.size(), 0.0);
        // Rescale the known-known mass by multiplicity^(freq_prior - 1), the allele-frequency
        // prior.
        double total = 0.0;
        for (size_t idx = 0; idx < post.size(); ++idx) {
            double k = known[idx];
            if (multiplicity[idx] > 1) {
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
    }

private:
    /// Add mass `g` of one state whose strands' haplotypes carry alleles `ai` and `bi`, or -1 for
    /// none.
    void add_state(int ai, int bi, double g) {
        const StrandEmission& sa = mt.strand[(size_t)(ai + 1)];
        const StrandEmission& sb = mt.strand[(size_t)(bi + 1)];
        const double e = mt.emission(ai, bi);
        if (!(e > 0.0)) {
            return;
        }
        // g / e is the mass per unit of P(i | a) P(j | b) L(i, j).
        const double w = g / e;
        if (ai >= 0 && bi >= 0) {
            const size_t idx = LinkageModel::genotype_index((size_t)ai, (size_t)bi);
            known[idx] += w * sa.keep * sb.keep * per_genotype[idx];
            multiplicity[idx] += 1;
            kept[(size_t)ai] += w * sa.keep * sb.base;
            kept[(size_t)bi] += w * sa.base * sb.keep;
        } else if (ai >= 0) {
            kept[(size_t)ai] += w * sa.keep * sb.base;
        } else if (bi >= 0) {
            kept[(size_t)bi] += w * sa.base * sb.keep;
        }
        free_mass += w * sa.base * sb.base;
    }

    /// Return the mass of each genotype in which at least one strand took any allele, spread by
    /// likelihood.
    vector<double> spread_wild_mass(double freq_prior) const {
        const size_t n = site.num_alleles;
        vector wild(site.genotype_ln_likelihood.size(), 0.0);
        if (n == 0) {
            return wild;
        }
        // Mass where one strand keeps its haplotype's allele k and the other takes any allele
        // arises once for each haplotype that carries k. So the number of k's carriers acts as a
        // frequency prior on it, and the mass is rescaled by carriers^(freq_prior - 1).
        const vector<size_t> carriers = allele_carriers(site, n_hap);
        for (size_t k = 0; k < n; ++k) {
            if (kept[k] <= 0.0) {
                continue;
            }
            double share = kept[k];
            if (carriers[k] > 1) {
                // At freq_prior 1 this changes nothing; above 1 it strengthens the prior.
                share /= pow((double)carriers[k], 1.0 - freq_prior);
            }
            for (size_t j = 0; j < n; ++j) {
                const size_t idx = LinkageModel::genotype_index(k, j);
                wild[idx] += share * per_genotype[idx];
            }
        }
        if (free_mass > 0.0) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t i = 0; i <= j; ++i) {
                    // Over ordered pairs, so a heterozygote is reached twice.
                    const size_t idx = LinkageModel::genotype_index(i, j);
                    wild[idx] += free_mass * per_genotype[idx] * (i == j ? 1.0 : 2.0);
                }
            }
        }
        return wild;
    }

    const LinkageModel::Site& site;
    /// The site's genotype likelihoods, scaled so the best is 1.
    const vector<double>& per_genotype;
    size_t n_hap;
    /// The site's emission terms, which give each strand's allele distribution.
    MutationTerms mt;
    /// Per genotype, the mass from states whose two strands keep their haplotypes' alleles.
    vector<double> known;
    /// Per genotype, how many ordered pairs of haplotypes spell it.
    vector<size_t> multiplicity;
    /// Per allele k, the mass from states where one strand keeps k and the other takes any allele.
    vector<double> kept;
    /// The mass from states where both strands take any allele.
    double free_mass = 0.0;
};

/// Return the forward messages over ordered pairs for the sites from `from` on, one per entry of
/// `emissions`, each rescaled to sum to 1.
///
/// The scale factors are discarded, since only posteriors are wanted, never the chain's total
/// likelihood. The first site starts from `alpha_in` where it has m * m entries, and from the
/// uniform distribution otherwise.
static vector<vector<double>> pair_forward(const LinkageModel& model,
                                           const vector<LinkageModel::Site>& sites, size_t from,
                                           const vector<vector<double>>& emissions, size_t m,
                                           const vector<double>* alpha_in) {
    const size_t n = emissions.size();
    vector<vector<double>> alpha(n);
    vector a(m * m, 1.0 / (double)(m * m));
    if (alpha_in != nullptr && alpha_in->size() == m * m) {
        a = *alpha_in;
    }
    for (size_t k = 0; k < m * m; ++k) {
        a[k] *= emissions[0][k];
    }
    rescale(a);
    alpha[0] = a;
    for (size_t t = 1; t < n; ++t) {
        const auto [gap_a, gap_b] = site_gap(sites[from + t - 1], sites[from + t]);
        vector<double> moved;
        const double rho_a = model.switch_probability(gap_a);
        const double rho_b = model.switch_probability(gap_b);
        transition_apply(alpha[t - 1], m, rho_a, rho_b, moved);
        for (size_t k = 0; k < m * m; ++k) {
            moved[k] *= emissions[t][k];
        }
        rescale(moved);
        alpha[t] = std::move(moved);
    }
    return alpha;
}

/// Return the backward message over ordered pairs at the site before `site`, from the message
/// `beta` and the `emission` at `site`, rescaled to sum to 1.
static vector<double> pair_backward_step(const LinkageModel& model,
                                         const vector<LinkageModel::Site>& sites, size_t site,
                                         const vector<double>& beta,
                                         const vector<double>& emission, size_t m) {
    const auto [gap_a, gap_b] = site_gap(sites[site - 1], sites[site]);
    vector weighted(m * m, 0.0);
    for (size_t k = 0; k < m * m; ++k) {
        weighted[k] = beta[k] * emission[k];
    }
    vector<double> next;
    const double rho_a = model.switch_probability(gap_a);
    const double rho_b = model.switch_probability(gap_b);
    transition_apply(weighted, m, rho_a, rho_b, next);
    rescale(next);
    return next;
}

void LinkageModel::window_posteriors(const vector<Site>& sites, size_t from, size_t to,
                                     vector<vector<double>>& out,
                                     const vector<double>* alpha_in,
                                     const vector<double>* beta_in,
                                     size_t out_base) const {
    const size_t n = to - from;
    if (n == 0) {
        return;
    }
    const size_t n_hap = panel_width(sites, from, to);
    const size_t m = n_hap;

    vector<vector<double>> emissions(n);
    vector<vector<double>> per_genotype(n);
    for (size_t t = 0; t < n; ++t) {
        build_emission(sites[from + t], n_hap, params.escape, params.mutation, emissions[t],
                       per_genotype[t]);
    }

    const vector<vector<double>> alpha = pair_forward(*this, sites, from, emissions, m, alpha_in);

    // Backward, combining as we go so only one beta is held at a time.
    vector beta(m * m, 1.0);
    if (beta_in != nullptr && beta_in->size() == m * m) {
        beta = *beta_in;
    }
    for (size_t t = n; t-- > 0;) {
        const Site& site = sites[from + t];
        // A site's own exponent is bounded so that mass * multiplicity^(f - 1), with mass <= 1 and
        // multiplicity <= n_hap^2, stays under e^690 and sums over genotypes cannot overflow.
        // That bound is f 122 on 18 haplotypes, and 60 on about 400.
        const double freq_prior =
            site.freq_prior >= 0.0
                ? std::min(site.freq_prior,
                           1.0 + 690.0 / std::log(std::max(4.0, (double)n_hap * (double)n_hap)))
                : params.freq_prior;
        DiploidSiteMass mass(site, per_genotype[t], n_hap, params.mutation, params.escape);
        mass.add_states(alpha[t], beta);
        mass.posterior(freq_prior, out[from + t - out_base]);

        if (t == 0) {
            break;
        }
        beta = pair_backward_step(*this, sites, from + t, beta, emissions[t], m);
    }
}

/// Fill `out` with a per-site quantity over `n` sites, decoded in overlapping windows whose
/// interiors are pasted together.
///
/// A site's posterior hardly depends on where a window edge fell, so windows can be decoded
/// independently and their interiors kept. `window(lo, hi, local)` fills `local` for `[lo, hi)`,
/// indexed from `lo`.
///
/// Each window is a task, decoded into a buffer of its own by whichever thread of the team is
/// free, and keeps only its interior, which no other window writes. Outside a parallel region the
/// tasks run one after another on the calling thread.
template <typename T, typename Window>
static void windowed_marginals(size_t n, size_t step, size_t margin, vector<T>& out,
                               Window window) {
    if (n <= step) {
        window(0, n, out);
        return;
    }
    const size_t windows = (n + step - 1) / step;
#pragma omp taskloop default(shared) grainsize(1)
    for (size_t w = 0; w < windows; ++w) {
        const size_t start = w * step;
        size_t lo = start > margin ? start - margin : 0;
        size_t hi = min(start + step + margin, n);
        vector<T> local(hi - lo);
        window(lo, hi, local);
        size_t keep_to = min(start + step, n);
        for (size_t t = start; t < keep_to; ++t) {
            out[t] = std::move(local[t - lo]);
        }
    }
}

/// Fill `out` with a path over `n` sites, decoded in overlapping windows chained by a pin.
///
/// Two windows decoded independently could choose different states where they meet, which would
/// put a spurious switch at every window boundary. So in each window after the first, the state
/// at one index inside its leading margin is pinned to what the previous window chose there. The
/// previous window chose that state with context on both sides.
///
/// `window(lo, hi, pin_index, pin, local)` fills `local` for `[lo, hi)`, indexed from `lo`.
/// `pin_index` is `(size_t)-1` for the first window.
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
                           // Only the window that begins the chain gets the entering message. A
                           // later window's left end is inside the chain, and its margin carries
                           // the chain across.
                           const vector<double>* enter = lo == 0 ? alpha_in : nullptr;
                           if (ploidy == 1) {
                               window_haploid_posteriors(sites, lo, hi, local, enter, lo);
                           } else {
                               window_posteriors(sites, lo, hi, local, enter, nullptr, lo);
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
    const size_t n = sites.size();
    if (ploidy == 1) {
        // One strand, so the path is over single haplotypes and `second` names none. Decode it
        // into its own vector, since the pin here is one haplotype rather than a pair.
        vector single(n, NO_HAPLOTYPE);
        windowed_path(n, max<size_t>(params.window, 1), params.margin, single,
                      [&](size_t lo, size_t hi, size_t pin_index, const size_t& pin,
                          vector<size_t>& local) {
                          // Only the window that begins the chain gets the entering message, as
                          // in `posteriors`.
                          local = window_haploid_phasing(sites, lo, hi, constraint, pin_index,
                                                         pin, lo == 0 ? alpha_in : nullptr);
                      });
        for (size_t t = 0; t < n; ++t) {
            out[t].first = single[t];
        }
        return out;
    }
    if (panel_width(sites, 0, n) == 0) {
        return out;
    }

    windowed_path(n, max<size_t>(params.window, 1), params.margin, out,
                  [&](size_t lo, size_t hi, size_t pin_index, const Phase& pin,
                      vector<Phase>& local) {
                      window_phasing(sites, lo, hi, constraint, pin_index, pin, local);
                  });
    return out;
}

/// Make every state of a site's emission impossible except `state`.
///
/// Where `state` is out of range, or is itself impossible, the emission is left alone, since
/// pinning it would leave the site with no state to be in.
static void pin_to_state(vector<double>& emission, size_t state) {
    const double keep = state < emission.size() ? emission[state] : 0.0;
    if (keep > 0.0) {
        const size_t states = emission.size();
        emission.assign(states, 0.0);
        emission[state] = keep;
    }
}

/// Replace a diploid site's emission over ordered pairs with the probability that each pair
/// spells the genotype `want`.
///
/// Every state can carry the wanted genotype, by keeping or mutating each strand's allele. So each
/// state takes the probability that its strands carry the genotype's alleles, in either order,
/// times the genotype's likelihood in `per_genotype`.
static void constrain_pair_emission(const LinkageModel::Site& site, size_t want, size_t n_hap,
                                    const vector<double>& per_genotype,
                                    const LinkageModel::Params& params,
                                    vector<double>& emission) {
    const size_t m = n_hap;
    // Each strand's probability is of one allele, so split the genotype into its alleles.
    size_t wi = 0;
    size_t wj = 0;
    genotype_alleles(want, wi, wj);

    const size_t n_alleles = site.num_alleles;
    vector<StrandEmission> strand(n_alleles + 1);
    for (size_t c = 0; c <= n_alleles; ++c) {
        strand[c] = strand_emission((int)c - 1, n_alleles, params.mutation, params.escape);
    }
    const double l = want < per_genotype.size() ? per_genotype[want] : 0.0;
    for (size_t a = 0; a < m; ++a) {
        const StrandEmission& sa = strand[(size_t)(allele_at(site, a, n_hap) + 1)];
        for (size_t b = 0; b < m; ++b) {
            const StrandEmission& sb = strand[(size_t)(allele_at(site, b, n_hap) + 1)];
            double p = sa.p(wi) * sb.p(wj);
            if (wi != wj) {
                p += sa.p(wj) * sb.p(wi);
            }
            emission[a * m + b] = p * l;
        }
    }
}

/// Pin each site that asks for it to its pinned pair, by zeroing the emission of every other state.
///
/// The pinned sites are those an earlier level already chose. A chosen site's phase is already in
/// the VCF, so letting the path re-orient it would decode this level's sites in a frame nothing
/// else uses. `emissions[t]` belongs to site `from + t`.
static void apply_site_pins(const vector<LinkageModel::Site>& sites, size_t from, size_t m,
                            vector<vector<double>>& emissions, LinkageCounters& counters) {
    for (size_t t = 0; t < emissions.size(); ++t) {
        const LinkageModel::Site& site = sites[from + t];
        if (!site.pinned) {
            continue;
        }
        const size_t pa = site.pin_first;
        const size_t pb = site.pin_second;
        if (pa >= m || pb >= m || emissions[t][pa * m + pb] <= 0.0) {
            // The pin names no haplotype on a strand, or the pinned pair cannot carry this site's
            // constrained genotype. Leave the site only constrained, and count it. In a group
            // whose only pinned site is its parent, a declined pin leaves the group free to swap
            // its strands relative to the parent.
            ++counters.pin_declined;
            continue;
        }
        ++counters.pin_applied;
        pin_to_state(emissions[t], pa * m + pb);
    }
}

/// The most haplotypes per strand that `pair_viterbi` can decode, since its backpointers are
/// 16 bits.
static constexpr size_t MAX_PAIR_PANEL = (size_t)numeric_limits<uint16_t>::max() + 1;

/// Return the most probable path of ordered pairs through a window, by max-product decoding.
///
/// `emissions[t]` is the emission of site `from + t` over ordered pairs, with `m` haplotypes per
/// strand. Every site gets `NO_HAPLOTYPE` on both strands where no state is possible at the
/// window's end.
static vector<LinkageModel::Phase> pair_viterbi(const LinkageModel& model,
                                                const vector<LinkageModel::Site>& sites,
                                                size_t from,
                                                const vector<vector<double>>& emissions, size_t m) {
    const size_t n = emissions.size();
    vector<LinkageModel::Phase> out(n);
    vector delta(m * m, NEG_INF);
    for (size_t k = 0; k < m * m; ++k) {
        if (emissions[0][k] > 0.0) {
            delta[k] = -log((double)(m * m)) + log(emissions[0][k]);
        }
    }
    // The backpointers name a haplotype in 16 bits. `window_phasing` refuses a wider panel.
    vector<vector<uint16_t>> back_a(n);
    vector<vector<uint16_t>> back_b(n);
    for (size_t t = 1; t < n; ++t) {
        const auto [gap_a, gap_b] = site_gap(sites[from + t - 1], sites[from + t]);
        vector<double> next;
        const double rho_a = model.switch_probability(gap_a);
        const double rho_b = model.switch_probability(gap_b);
        MaxProductStep(delta, m, rho_a, rho_b).apply(emissions[t], next, back_a[t], back_b[t]);
        if (std::ranges::none_of(next, [](double v) { return v > NEG_INF; })) {
            // No state survives, which happens only where the chosen genotype has no likelihood.
            // Restart the path here, so the rest of the window is still phased. The restart
            // shows as a switch on both strands.
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
        return out;
    }
    // Trace the backpointers from the best final state.
    size_t a = best / m;
    size_t b = best % m;
    for (size_t t = n; t-- > 0;) {
        out[t].first = a;
        out[t].second = b;
        if (t > 0) {
            const size_t pa = back_a[t][a * m + b];
            const size_t pb = back_b[t][a * m + b];
            a = pa;
            b = pb;
        }
    }
    return out;
}

void LinkageModel::window_phasing(const vector<Site>& sites, size_t from, size_t to,
                                  const vector<size_t>& constraint,
                                  size_t pin_index, const Phase& pin,
                                  vector<Phase>& out) const {
    const size_t n = to - from;
    out.assign(n, Phase{});
    if (n == 0) {
        return;
    }
    const size_t n_hap = panel_width(sites, from, to);
    const size_t m = n_hap;
    if (m > MAX_PAIR_PANEL) {
        // A wider panel would need at least 2^32 states per site, and its haplotype indices
        // would not fit in the backpointers of `pair_viterbi`.
        throw runtime_error("vg call: the diploid linkage model can phase against at most "
                            + to_string(MAX_PAIR_PANEL) + " haplotypes, but the panel has "
                            + to_string(m));
    }

    // Emissions, with the constraint folded in.
    vector<vector<double>> emissions(n);
    vector<vector<double>> per_genotype(n);
    for (size_t t = 0; t < n; ++t) {
        const Site& site = sites[from + t];
        build_emission(site, n_hap, params.escape, params.mutation, emissions[t], per_genotype[t]);
        const size_t want = (from + t) < constraint.size() ? constraint[from + t] : NO_CONSTRAINT;
        if (want != NO_CONSTRAINT) {
            constrain_pair_emission(site, want, n_hap, per_genotype[t], params, emissions[t]);
        }
    }

    apply_site_pins(sites, from, m, emissions, counter_set);

    // Pin the seam with the previous window, unless the pin is impossible under this window's
    // constraint. With margin 0 the pin index falls outside [from, to) and is skipped.
    if (pin_index != (size_t)-1 && pin_index >= from && pin_index < to) {
        const size_t pa = pin.first;
        const size_t pb = pin.second;
        pin_to_state(emissions[pin_index - from], (pa < m && pb < m) ? pa * m + pb : m * m);
    }

    out = pair_viterbi(*this, sites, from, emissions, m);
}


//------------------------------------------------------------------------------
// Ploidy-1 chains
//
// The states are single panel haplotypes, so a genotype is one allele.

void LinkageModel::haploid_emission(const Site& site, size_t n_hap, vector<double>& e,
                                    vector<double>& per_allele) const {
    size_t m = n_hap;
    size_t n = site.num_alleles;

    // Shift so the best allele is exp(0) = 1, as build_emission does for genotypes. This keeps
    // the numbers in range without changing any ratio.
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

    double sum = 0.0;
    for (double v : per_allele) {
        sum += v;
    }

    // Sum over the allele the strand carries, as the diploid emission does. That gives
    // base * (sum of the likelihoods) + keep * (the likelihood of the haplotype's own allele).
    e.assign(m, 0.0);
    for (size_t a = 0; a < m; ++a) {
        const int ai = allele_at(site, a, n_hap);
        const StrandEmission s = strand_emission(ai, n, params.mutation, params.escape);
        e[a] = s.base * sum + (ai >= 0 ? s.keep * per_allele[(size_t)ai] : 0.0);
    }
}

/// Return the forward messages over single haplotypes for the sites from `from` on, one per entry
/// of `emissions`, each rescaled to sum to 1.
///
/// The first site starts from `alpha_in` where it has `m` entries, and from the uniform
/// distribution otherwise. Where a given message and the first site's reads disagree outright,
/// so that every state the message allows has zero emission, the first site starts from its reads
/// alone, as `single_viterbi` does.
static vector<vector<double>> single_forward(const LinkageModel& model,
                                             const vector<LinkageModel::Site>& sites, size_t from,
                                             const vector<vector<double>>& emissions, size_t m,
                                             const vector<double>* alpha_in) {
    const size_t n = emissions.size();
    vector<vector<double>> alpha(n);
    vector a(m, 1.0 / (double)m);
    if (alpha_in != nullptr && alpha_in->size() == m) {
        a = *alpha_in;
    }
    double sum = 0.0;
    for (size_t k = 0; k < m; ++k) {
        a[k] *= emissions[0][k];
        sum += a[k];
    }
    if (sum <= 0.0 && alpha_in != nullptr) {
        for (size_t k = 0; k < m; ++k) {
            a[k] = emissions[0][k] / (double)m;
        }
    }
    rescale(a);
    alpha[0] = a;
    for (size_t t = 1; t < n; ++t) {
        // One strand, so one distance.
        const double rho =
            model.switch_probability(site_gap(sites[from + t - 1], sites[from + t]).first);
        const double stay = 1.0 - rho;
        const double jump = rho / (double)m;
        double total = 0.0;
        for (double v : alpha[t - 1]) {
            total += v;
        }
        vector next(m, 0.0);
        for (size_t k = 0; k < m; ++k) {
            next[k] = (stay * alpha[t - 1][k] + jump * total) * emissions[t][k];
        }
        rescale(next);
        alpha[t] = std::move(next);
    }
    return alpha;
}

/// Return the backward message over single haplotypes at the site before `site`, from the message
/// `beta` and the `emission` at `site`, rescaled to sum to 1.
static vector<double> single_backward_step(const LinkageModel& model,
                                           const vector<LinkageModel::Site>& sites, size_t site,
                                           const vector<double>& beta,
                                           const vector<double>& emission, size_t m) {
    const double rho = model.switch_probability(site_gap(sites[site - 1], sites[site]).first);
    const double stay = 1.0 - rho;
    const double jump = rho / (double)m;
    vector<double> weighted(m);
    double total_w = 0.0;
    for (size_t k = 0; k < m; ++k) {
        weighted[k] = beta[k] * emission[k];
        total_w += weighted[k];
    }
    vector prev(m, 0.0);
    for (size_t k = 0; k < m; ++k) {
        prev[k] = stay * weighted[k] + jump * total_w;
    }
    rescale(prev);
    return prev;
}

/// Write a haploid site's posterior over alleles to `post`, from the mass `state_mass[a]` =
/// alpha * beta of each haplotype, or leave `post` empty where no state explains the site.
///
/// A state's mass is split by how its strand got its allele, as in `DiploidSiteMass`. The share
/// where the strand keeps its haplotype's allele goes to that allele, multiplied by the allele's
/// number of carriers to the power `freq_prior` - 1, so that `freq_prior` 0 means what it says.
/// The uniform share is spread over the alleles by likelihood.
static void single_site_posterior(const LinkageModel::Site& site, const vector<double>& emission,
                                  const vector<double>& per_allele,
                                  const vector<double>& state_mass, size_t n_hap,
                                  const LinkageModel::Params& params, vector<double>& post) {
    const size_t m = n_hap;
    const double freq_prior = site.freq_prior >= 0.0 ? site.freq_prior : params.freq_prior;
    post.assign(site.num_alleles, 0.0);
    const vector<size_t> carriers = allele_carriers(site, n_hap);
    vector known(site.num_alleles, 0.0);
    vector wild(site.num_alleles, 0.0);
    for (size_t a = 0; a < m; ++a) {
        const double g = state_mass[a];
        if (g <= 0.0 || !(emission[a] > 0.0)) {
            continue;
        }
        const int ai = allele_at(site, a, n_hap);
        const StrandEmission st =
            strand_emission(ai, site.num_alleles, params.mutation, params.escape);
        const double w = g / emission[a];
        if (ai >= 0) {
            known[(size_t)ai] += w * st.keep * per_allele[(size_t)ai];
        }
        for (size_t k = 0; k < site.num_alleles; ++k) {
            wild[k] += w * st.base * per_allele[k];
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
    } else {
        // The caller keeps the site's own call, as in the diploid pass.
        post.clear();
    }
}

void LinkageModel::window_haploid_posteriors(const vector<Site>& sites, size_t from, size_t to,
                                             vector<vector<double>>& out,
                                             const vector<double>* alpha_in,
                                             size_t out_base) const {
    const size_t n = to - from;
    if (n == 0) {
        return;
    }
    const size_t n_hap = panel_width(sites, from, to);
    const size_t m = n_hap;

    vector<vector<double>> emissions(n);
    vector<vector<double>> per_allele(n);
    for (size_t t = 0; t < n; ++t) {
        haploid_emission(sites[from + t], n_hap, emissions[t], per_allele[t]);
    }

    const vector<vector<double>> alpha =
        single_forward(*this, sites, from, emissions, m, alpha_in);

    // Backward, combined as we go so only one beta is held.
    vector beta(m, 1.0);
    for (size_t t = n; t-- > 0;) {
        vector<double> state_mass(m);
        for (size_t a = 0; a < m; ++a) {
            state_mass[a] = alpha[t][a] * beta[a];
        }
        single_site_posterior(sites[from + t], emissions[t], per_allele[t], state_mass, n_hap,
                              params, out[from + t - out_base]);
        if (t == 0) {
            break;
        }
        beta = single_backward_step(*this, sites, from + t, beta, emissions[t], m);
    }
}

/// Return the log scores of a haploid max-product decode at its first site, from the site's
/// `emission` and the entering message `alpha_in`.
///
/// The message multiplies the emission, as a prior over the first state, where it has `m`
/// entries. Where every state the message allows has zero emission, the message and the reads
/// disagree outright. Then the scores come from the reads alone. A decode with no live state
/// would name no haplotype anywhere, which would say the panel cannot explain a strand that it
/// explains once the parent's claim is dropped.
static vector<double> single_viterbi_start(const vector<double>& emission, size_t m,
                                           const vector<double>* alpha_in) {
    vector delta(m, NEG_INF);
    const bool have_alpha = alpha_in != nullptr && alpha_in->size() == m;
    for (size_t k = 0; k < m; ++k) {
        if (emission[k] <= 0.0) {
            continue;
        }
        const double prior = have_alpha ? (*alpha_in)[k] : 1.0 / (double)m;
        if (prior > 0.0) {
            delta[k] = log(prior) + log(emission[k]);
        }
    }
    if (have_alpha && std::ranges::none_of(delta, [](double v) { return v > NEG_INF; })) {
        for (size_t k = 0; k < m; ++k) {
            if (emission[k] > 0.0) {
                delta[k] = -log((double)m) + log(emission[k]);
            }
        }
    }
    return delta;
}

/// Return the log scores of a haploid max-product decode at the next site, from the scores
/// `delta` at the previous one, a switch probability `rho` and the next site's `emission`. Set
/// `back[k]` to the best predecessor of each state k.
///
/// With one strand, a state either stays on the same haplotype or jumps from whichever was best.
/// The diploid step needs leave-one-out maxima only because its two strands share one delta.
static vector<double> single_max_product_step(const vector<double>& delta, double rho,
                                              const vector<double>& emission,
                                              vector<uint32_t>& back) {
    const size_t m = delta.size();
    rho = min(max(rho, 1e-12), 1.0 - 1e-12);
    const double S = log(1.0 - rho + rho / (double)m);
    const double J = log(rho / (double)m);

    size_t arg_best = 0;
    double best_any = NEG_INF;
    for (size_t k = 0; k < m; ++k) {
        if (delta[k] > best_any) {
            best_any = delta[k];
            arg_best = k;
        }
    }
    vector next(m, NEG_INF);
    back.assign(m, 0);
    for (size_t k = 0; k < m; ++k) {
        if (!(emission[k] > 0.0)) {
            continue;
        }
        const double stay = delta[k] > NEG_INF ? delta[k] + S : NEG_INF;
        const double jump = best_any > NEG_INF ? best_any + J : NEG_INF;
        if (stay >= jump && stay > NEG_INF) {
            next[k] = stay + log(emission[k]);
            back[k] = (uint32_t)k;
        } else if (jump > NEG_INF) {
            next[k] = jump + log(emission[k]);
            back[k] = (uint32_t)arg_best;
        }
    }
    if (std::ranges::none_of(next, [](double v) { return v > NEG_INF; })) {
        // No state survives. Restart here, as in the diploid pass.
        for (size_t k = 0; k < m; ++k) {
            next[k] = emission[k] > 0.0 ? log(emission[k]) : NEG_INF;
            back[k] = (uint32_t)k;
        }
    }
    return next;
}

/// Return the most probable path of single haplotypes through a window, by max-product decoding,
/// starting from the entering message `alpha_in` as `single_viterbi_start` describes.
///
/// `emissions[t]` is the emission of site `from + t` over `m` haplotypes. Every site gets
/// `NO_HAPLOTYPE` where no state is possible at the window's end.
static vector<size_t> single_viterbi(const LinkageModel& model,
                                     const vector<LinkageModel::Site>& sites, size_t from,
                                     const vector<vector<double>>& emissions, size_t m,
                                     const vector<double>* alpha_in) {
    const size_t n = emissions.size();
    vector out(n, LinkageModel::NO_HAPLOTYPE);
    vector<double> delta = single_viterbi_start(emissions[0], m, alpha_in);
    vector<vector<uint32_t>> back(n);
    for (size_t t = 1; t < n; ++t) {
        const double rho =
            model.switch_probability(site_gap(sites[from + t - 1], sites[from + t]).first);
        delta = single_max_product_step(delta, rho, emissions[t], back[t]);
    }

    size_t best = m;
    double best_v = NEG_INF;
    for (size_t k = 0; k < m; ++k) {
        if (delta[k] > best_v) {
            best_v = delta[k];
            best = k;
        }
    }
    if (best == m) {
        return out;
    }
    // Trace the backpointers from the best final state.
    size_t cur = best;
    for (size_t t = n; t-- > 0;) {
        out[t] = cur;
        if (t > 0) {
            cur = back[t][cur];
        }
    }
    return out;
}

vector<size_t> LinkageModel::window_haploid_phasing(const vector<Site>& sites, size_t from,
                                                    size_t to, const vector<size_t>& constraint,
                                                    size_t pin_index, size_t pin,
                                                    const vector<double>* alpha_in) const {
    const size_t n = to - from;
    if (n == 0) {
        return {};
    }
    const size_t n_hap = panel_width(sites, from, to);
    const size_t m = n_hap;

    vector<vector<double>> emissions(n);
    vector<vector<double>> per_allele(n);
    for (size_t t = 0; t < n; ++t) {
        const Site& site = sites[from + t];
        haploid_emission(site, n_hap, emissions[t], per_allele[t]);
        const size_t want = (from + t) < constraint.size() ? constraint[from + t] : NO_CONSTRAINT;
        if (want == NO_CONSTRAINT) {
            continue;
        }
        // As in the diploid case, any state can carry the called allele, by keeping or mutating
        // its own. So each takes the probability that it does times the allele's likelihood.
        const double e = want < per_allele[t].size() ? per_allele[t][want] : 0.0;
        for (size_t a = 0; a < m; ++a) {
            const int ai = allele_at(site, a, n_hap);
            emissions[t][a] =
                strand_emission(ai, site.num_alleles, params.mutation, params.escape).p(want) * e;
        }
    }
    if (pin_index != (size_t)-1 && pin_index >= from && pin_index < to) {
        pin_to_state(emissions[pin_index - from], pin < m ? pin : m);
    }

    return single_viterbi(*this, sites, from, emissions, m, alpha_in);
}

//------------------------------------------------------------------------------
// EntryTable

size_t EntryTable::Entry::called_genotype() const {
    return genotype_of(ploidy, called_i, called_j);
}

size_t EntryTable::Entry::chosen_genotype() const {
    return genotype_of(ploidy, final_i, final_j);
}

uint32_t EntryTable::contig_id(const string& contig) {
    const auto [contig_at, new_contig] =
        contig_index.try_emplace(contig, (uint32_t)contig_names.size());
    if (new_contig) {
        contig_names.push_back(contig);
    }
    return contig_at->second;
}

void EntryTable::add(Entry e, const CompactSite& cs, const vector<int>& haplotype_traversal,
                     const vector<int>& traversal_to_allele) {
    append_site_data(e, cs, haplotype_traversal, traversal_to_allele);

    const auto at = (uint32_t)entries.size();
    const size_t record_key = e.record_key;
    if (auto last = last_by_key.find(record_key); last == last_by_key.end()) {
        first_by_key[record_key] = at;
    } else {
        // Count a duplicate before the append, while the key's chain still describes the existing
        // entries. See `num_duplicate_live_keys`.
        if (live_index(record_key) != NO_ENTRY) {
            ++duplicate_live_keys;
        }
        entries[last->second].next_same_key = at;
    }
    last_by_key[record_key] = at;
    if (e.level >= by_level.size()) {
        by_level.resize((size_t)e.level + 1);
    }
    by_level[e.level].push_back(at);
    entries.push_back(e);
}

void EntryTable::append_site_data(Entry& e, const CompactSite& cs,
                                  const vector<int>& haplotype_traversal,
                                  const vector<int>& traversal_to_allele) {
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
    for (int trav : cs.space) {
        trav_arena.push_back((uint16_t)trav);
        const int allele =
            (size_t)trav < traversal_to_allele.size() ? traversal_to_allele[trav] : -1;
        allele_arena.push_back(allele >= 0 && allele < 127 ? (int8_t)allele : (int8_t)-1);
    }
    e.num_alleles = (uint16_t)cs.space.size();
}

bool EntryTable::replace_site_data(Entry& e, const CompactSite& cs,
                                   const vector<int>& haplotype_traversal) {
    bool same_space = cs.space.size() == e.num_alleles;
    for (size_t i = 0; same_space && i < cs.space.size(); ++i) {
        same_space = trav_arena[e.trav_offset + i] == (uint16_t)cs.space[i];
    }
    if (same_space) {
        // Overwrite the likelihoods in place. `gl_arena` stores no length, since a reader
        // recomputes it from `num_alleles`, so this needs the shape to be unchanged.
        if (const size_t expect = e.ploidy == 1
                                      ? (size_t)e.num_alleles
                                      : (size_t)e.num_alleles * ((size_t)e.num_alleles + 1) / 2;
            cs.gls.size() != expect) {
            return false;
        }
        std::ranges::copy(cs.gls, gl_arena.begin() + e.gl_offset);
    } else {
        // A slice of another width does not fit in place, so append new slices and leave the old
        // ones unused. Every other entry's offsets stay valid. No allele map is known yet, and
        // `set_vcf_alleles` supplies it when the record is written.
        append_site_data(e, cs, haplotype_traversal, {});
    }
    return true;
}

void EntryTable::set_vcf_alleles(const Entry& e, const vector<int>& traversal_to_allele) {
    // Rewrite the span `add` filled in place. The compact space has not changed, only what each
    // of its alleles is called in the VCF.
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
}

bool EntryTable::retract(size_t record_key) {
    const uint32_t found = live_index(record_key);
    if (found == NO_ENTRY) {
        return false;
    }
    entries[found].retracted = true;
    return true;
}

uint32_t EntryTable::live_index(size_t record_key) const {
    auto it = first_by_key.find(record_key);
    for (uint32_t i = it == first_by_key.end() ? NO_ENTRY : it->second;
         i != NO_ENTRY; i = entries[i].next_same_key) {
        if (!entries[i].retracted) {
            return i;
        }
    }
    return NO_ENTRY;
}

const vector<uint32_t>& EntryTable::at_level(size_t level) const {
    static const vector<uint32_t> no_entries;
    return level < by_level.size() ? by_level[level] : no_entries;
}

int EntryTable::traversal_at(const Entry& e, size_t compact) const {
    if (compact >= e.num_alleles) {
        return -1;
    }
    const size_t at = e.trav_offset + compact;
    return at < trav_arena.size() ? (int)trav_arena[at] : -1;
}

int EntryTable::vcf_allele_at(const Entry& e, size_t compact) const {
    if (compact >= e.num_alleles) {
        return -1;
    }
    const size_t at = e.allele_offset + compact;
    return at < allele_arena.size() ? (int)allele_arena[at] : -1;
}

size_t EntryTable::bytes() const {
    // Every arena, and the key index, which costs two hash nodes per site.
    return entries.size() * sizeof(Entry)
           + gl_arena.size() * sizeof(float)
           + hap_arena.size() * sizeof(int8_t)
           + trav_arena.size() * sizeof(uint16_t)
           + allele_arena.size() * sizeof(int8_t)
           + (first_by_key.size() + last_by_key.size())
                 * (sizeof(std::pair<const size_t, uint32_t>) + sizeof(void*));
}

//------------------------------------------------------------------------------
// LinkageCollector

vector<int> LinkageCollector::compact_allele_space(const vector<int>& haplotype_traversal,
                                                   int called_trav_i, int called_trav_j) {
    // Keep only the traversals that can matter, so that the space and its likelihood vector stay
    // small. build_emission reads the likelihoods only at pairs of panel-carried alleles, and the
    // constraint needs only the called pair.
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
    // Sort by candidate index, so the compact numbering is a function of the site alone and not
    // of the order haplotypes happen to be visited in.
    return vector<int>(needed.begin(), needed.end());
}

CompactSite LinkageCollector::compact_site(
        const map<vector<int>, double>& genotype_ln_likelihood,
        const vector<int>& haplotype_traversal,
        int called_trav_i, int called_trav_j, size_t ploidy) const {
    CompactSite cs;
    cs.space = compact_allele_space(haplotype_traversal, called_trav_i, called_trav_j);
    if (cs.space.empty() || cs.space.size() > 127) {
        return cs;
    }
    cs.ci = cs.compact_of(called_trav_i);
    cs.cj = called_trav_j >= 0 ? cs.compact_of(called_trav_j) : cs.ci;
    if (cs.ci < 0 || cs.cj < 0) {
        // The called pair is in the space by construction, so this cannot happen.
        return cs;
    }
    const size_t k = cs.space.size();
    cs.site_ploidy = (ploidy == 1 ? 1 : 2);
    const size_t n_gt = cs.site_ploidy == 1 ? k : k * (k + 1) / 2;
    cs.gls.assign(n_gt, -numeric_limits<float>::infinity());
    for (const auto& [genotype, ln_likelihood] : genotype_ln_likelihood) {
        if (genotype.size() != cs.site_ploidy) {
            continue;
        }
        if (cs.site_ploidy == 1) {
            if (const int a = cs.compact_of(genotype[0]); a >= 0) {
                cs.gls[(size_t)a] = (float)ln_likelihood;
            }
            continue;
        }
        const int a = cs.compact_of(genotype[0]);
        const int b = cs.compact_of(genotype[1]);
        if (a >= 0 && b >= 0) {
            cs.gls[LinkageModel::genotype_index((size_t)a, (size_t)b)] = (float)ln_likelihood;
        }
    }
    cs.ok = true;
    return cs;
}

/// Return whether the allele `a` is the shorter allele `b` with 1-49 extra copies of one base
/// inserted into a run of that base, where the run in `a` is at least `min_run` long or reaches
/// an end of `a`.
///
/// A run that reaches an end of the allele continues into a boundary node made of this base,
/// where the graph cut it. Its length is unknown, so it counts as long.
static bool homopolymer_run_change(const string& a, const string& b, size_t min_run) {
    if (a.size() <= b.size() || a.size() - b.size() > 49) {
        return false;
    }
    size_t pre = 0;
    while (pre < b.size() && a[pre] == b[pre]) {
        ++pre;
    }
    size_t suf = 0;
    while (suf < b.size() - pre && a[a.size() - 1 - suf] == b[b.size() - 1 - suf]) {
        ++suf;
    }
    if (pre + suf != b.size()) {
        // The alleles differ by more than an insertion.
        return false;
    }
    // What `a` has and `b` lacks, once the shared flanks are gone.
    const size_t lo = pre;
    const size_t hi = a.size() - suf;
    const char base = a[lo];
    if (!std::all_of(a.begin() + lo, a.begin() + hi, [base](char c) { return c == base; })) {
        return false;
    }
    // Measure the run the extra bases sit in, in the longer allele.
    size_t start = lo;
    while (start > 0 && a[start - 1] == base) {
        --start;
    }
    size_t end = hi;
    while (end < a.size() && a[end] == base) {
        ++end;
    }
    return end - start >= min_run || start == 0 || end == a.size();
}

bool LinkageModel::run_length_site(const vector<string>& alleles, size_t min_run, size_t ref) {
    for (size_t x = 0; x < alleles.size(); ++x) {
        for (size_t y = 0; y < alleles.size(); ++y) {
            if (ref < alleles.size() && x != ref && y != ref) {
                continue;
            }
            // Test with `x` as the longer allele. The other ordering is tested when y comes
            // round as x.
            if (homopolymer_run_change(alleles[x], alleles[y], min_run)) {
                return true;
            }
        }
    }
    return false;
}

void LinkageCollector::record(const SiteLocus& locus,
                              const map<vector<int>, double>& genotype_ln_likelihood,
                              const vector<int>& haplotype_traversal, const DirectCall& call,
                              const vector<int>& traversal_to_allele, size_t record_key,
                              const SiteContext& ctx) {
    if (genotype_ln_likelihood.empty() || call.called_i < 0) {
        return;
    }
    // A site that cannot be described, such as one with more alleles than the int8 panel arena
    // can name, loses its linkage rather than being mis-linked.
    const CompactSite cs = compact_site(genotype_ln_likelihood, haplotype_traversal,
                                        call.called_i, call.called_j, call.ploidy);
    if (!cs.ok) {
        return;
    }

    std::scoped_lock guard(mutex);

    Entry e;
    e.position = (uint32_t)locus.position;
    e.contig = entries.contig_id(locus.contig);
    e.called_i = (uint16_t)cs.ci;
    e.called_j = (uint16_t)cs.cj;
    e.record_key = record_key;
    e.explained_share = (float)call.quality.explained_share;
    e.gq_factor = (float)call.quality.gq_factor;
    e.achievable_gap = (float)call.quality.achievable_gap;
    e.start_node = ctx.start_node;
    e.end_node = ctx.end_node;
    e.ploidy = (uint8_t)cs.site_ploidy;
    e.nested = ctx.nested;
    e.emitted = ctx.emitted;
    e.parent_record_key = ctx.parent_record_key;
    e.parent_crossing = ctx.parent_crossing;
    e.unpositioned = ctx.unpositioned;
    e.chain_key = ctx.chain_key;
    e.level = (uint8_t)(ctx.level > 255 ? 255 : ctx.level);
    e.freq_prior = (float)ctx.freq_prior;
    entries.add(e, cs, haplotype_traversal, traversal_to_allele);
}

size_t LinkageCollector::num_sites_at(size_t level) const {
    size_t n = 0;
    for (const Entry& e : entries) {
        if (e.level == level && !e.retracted) {
            ++n;
        }
    }
    return n;
}

size_t LinkageCollector::max_level() const {
    size_t g = 0;
    for (const Entry& e : entries) {
        g = max(g, (size_t)e.level);
    }
    return g;
}

bool LinkageCollector::chosen_traversals(size_t record_key, int* first, int* second,
                                         size_t* ploidy) const {
    std::scoped_lock guard(mutex);
    if (const uint32_t found = entries.live_index(record_key); found != EntryTable::NO_ENTRY) {
        const Entry& e = entries[found];
        const int a = entries.traversal_at(e, e.final_i);
        const int b = entries.traversal_at(e, e.final_j);
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

void LinkageCollector::chosen_traversals_for(const vector<size_t>& record_keys,
                                             vector<std::array<int, 3>>& out,
                                             vector<char>& found) const {
    std::scoped_lock guard(mutex);
    out.assign(record_keys.size(), {-1, -1, 0});
    found.assign(record_keys.size(), 0);
    // Each lookup only reads, and each writes its own slot.
#pragma omp parallel for schedule(dynamic, 4096)
    for (size_t i = 0; i < record_keys.size(); ++i) {
        const uint32_t idx = entries.live_index(record_keys[i]);
        if (idx == EntryTable::NO_ENTRY) {
            continue;
        }
        const Entry& e = entries[idx];
        const int a = entries.traversal_at(e, e.final_i);
        const int b = entries.traversal_at(e, e.final_j);
        if (a < 0 || b < 0) {
            continue;
        }
        out[i] = {a, b, (int)e.ploidy};
        found[i] = 1;
    }
}

std::unordered_set<size_t> LinkageCollector::emitted_records() const {
    std::scoped_lock guard(mutex);
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
    std::scoped_lock guard(mutex);
    if (const uint32_t found = entries.live_index(record_key); found != EntryTable::NO_ENTRY) {
        Entry& e = entries[found];
        // Only ever set `emitted`, never clear it. When every snarl is genotyped independently,
        // nested ones included, a record can be rendered by more than one path, and an assignment
        // would make the result depend on which thread took the mutex last.
        if (emitted && !e.emitted) {
            e.emitted = true;
        }
        entries.set_vcf_alleles(e, traversal_to_allele);
        return true;
    }
    return false;
}

bool LinkageCollector::has_entry(size_t record_key) const {
    std::scoped_lock guard(mutex);
    return entries.live_index(record_key) != EntryTable::NO_ENTRY;
}

bool LinkageCollector::set_position(size_t record_key, size_t position) {
    std::scoped_lock guard(mutex);
    const uint32_t found = entries.live_index(record_key);
    if (found == EntryTable::NO_ENTRY) {
        return false;
    }
    entries[found].position = (uint32_t)position;
    return true;
}

bool LinkageCollector::rescore(size_t record_key,
                               const map<vector<int>, double>& genotype_ln_likelihood,
                               const vector<int>& haplotype_traversal,
                               int called_trav_i, int called_trav_j) {
    std::scoped_lock guard(mutex);
    const uint32_t found = entries.live_index(record_key);
    if (found == EntryTable::NO_ENTRY) {
        return false;
    }
    Entry& e = entries[found];
    const CompactSite cs = compact_site(genotype_ln_likelihood, haplotype_traversal,
                                        called_trav_i, called_trav_j, e.ploidy);
    if (!cs.ok) {
        return false;
    }
    if (!entries.replace_site_data(e, cs, haplotype_traversal)) {
        return false;
    }
    // The called pair is re-compacted with the new space. `final_*` is left alone, since
    // `resolve_level` resets it from `called_*` before it decodes.
    e.called_i = (uint16_t)cs.ci;
    e.called_j = (uint16_t)cs.cj;
    return true;
}

bool LinkageCollector::retract(size_t record_key) {
    std::scoped_lock guard(mutex);
    if (!entries.retract(record_key)) {
        return false;
    }
    // A key with no live entry has no chosen genotype, so it is no longer moved.
    if (entries.live_index(record_key) == EntryTable::NO_ENTRY) {
        moved_quality_by_record.erase(record_key);
    }
    return true;
}

void LinkageCollector::finish_phase_call(PhaseCall& pc, const Entry& e, size_t c_first,
                                         size_t c_second) const {
    pc.trav_first = entries.traversal_at(e, c_first);
    pc.trav_second = entries.traversal_at(e, c_second);
    int v_first = entries.vcf_allele_at(e, c_first);
    int v_second = entries.vcf_allele_at(e, c_second);
    if (v_first < 0 || v_second < 0) {
        // The allele map names no VCF allele for one of the pair, because it was not supplied yet
        // or because the record wrote no ALT for that traversal. Write the called pair instead,
        // whose order nothing determined.
        v_first = entries.vcf_allele_at(e, e.called_i);
        v_second = entries.vcf_allele_at(e, e.called_j);
        pc.order_arbitrary = pc.order_arbitrary || (v_first != v_second);
    }
    pc.allele_first = v_first >= 0 ? (size_t)v_first : numeric_limits<size_t>::max();
    pc.allele_second = v_second >= 0 ? (size_t)v_second : numeric_limits<size_t>::max();
}


int LinkageCollector::PinnedPhase::strand_carrying(int carrying) const {
    if (carrying < 0) {
        return -1;
    }
    if (ploidy == 1 && nested_strand >= 0) {
        // A nested parent at ploidy 1 is on one strand, so everything inside it is on that
        // strand.
        return trav_first == carrying ? nested_strand : -1;
    }
    if (ploidy == 2) {
        if (trav_first == carrying) {
            return 0;
        }
        if (trav_second == carrying) {
            return 1;
        }
    }
    // A haploid top-level parent has only one strand, so its children are not written as `a|.`.
    return -1;
}

LinkageCollector::Relation LinkageCollector::relation_to(const Entry& child,
                                                         const Entry& parent) const {
    const int ta = entries.traversal_at(parent, parent.final_i);
    const int tb = parent.ploidy == 2 ? entries.traversal_at(parent, parent.final_j) : -1;
    return relate_to_parent(child.parent_crossing, ta, tb);
}

bool LinkageCollector::before_in_position(size_t a, size_t b) const {
    if (entries[a].position != entries[b].position) {
        return entries[a].position < entries[b].position;
    }
    return entries[a].record_key < entries[b].record_key;
}

unordered_map<size_t, LinkageCollector::PinnedPhase>
LinkageCollector::pinned_phases(const unordered_set<size_t>& parent_keys,
                                const vector<PhaseCall>& phasing) {
    unordered_map<size_t, PinnedPhase> pinned_phase;
    pinned_phase.reserve(parent_keys.size() * 2);
    for (const PhaseCall& pc : phasing) {
        if (parent_keys.contains(pc.record_key)) {
            pinned_phase[pc.record_key] = PinnedPhase{pc.hap_first, pc.hap_second,
                                                      pc.trav_first, pc.trav_second,
                                                      pc.ploidy, (int)pc.nested_strand,
                                                      pc.order_arbitrary, pc.phase_set};
        }
    }
    return pinned_phase;
}

vector<LinkageCollector::DecodeChain>
LinkageCollector::top_level_chains(const vector<uint32_t>& this_level) {
    // Group by contig, then sort by reference position. Node-ID order is close to reference order
    // but not the same, and the transition probabilities depend on the distances.
    vector<vector<size_t>> by_contig(entries.num_contigs());
    for (uint32_t i : this_level) {
        if (entries[i].retracted) {
            continue;
        }
        by_contig[entries[i].contig].push_back(i);
    }
    // A top-level linkage chain is a maximal run of one ploidy on one contig, since strands do not
    // correspond across a ploidy change such as the boundary of a region with its own ploidy.
    vector<DecodeChain> chains;
    for (vector<size_t>& contig_indices : by_contig) {
        std::ranges::sort(contig_indices,
                          [this](size_t a, size_t b) { return before_in_position(a, b); });
        size_t run_start = 0;
        while (run_start < contig_indices.size()) {
            size_t run_end = run_start + 1;
            while (run_end < contig_indices.size()
                   && entries[contig_indices[run_end]].ploidy
                          == entries[contig_indices[run_start]].ploidy) {
                ++run_end;
            }
            DecodeChain chain;
            chain.sites.assign(contig_indices.begin() + run_start,
                               contig_indices.begin() + run_end);
            chains.push_back(std::move(chain));
            run_start = run_end;
        }
    }
    return chains;
}

vector<LinkageCollector::DecodeChain> LinkageCollector::parent_groups(
        const vector<uint32_t>& this_level, const unordered_set<size_t>& parent_keys,
        const unordered_map<size_t, PinnedPhase>& pinned_phase,
        unordered_map<size_t, NestedPlacement>& placements, size_t& grouped_sites,
        size_t& grouped_groups) {
    // The live entry of each parent key; a key's last live entry wins.
    unordered_map<size_t, size_t> index_of_key;
    index_of_key.reserve(parent_keys.size() * 2);
    for (size_t i = 0; i < entries.size(); ++i) {
        if (!entries[i].retracted && parent_keys.contains(entries[i].record_key)) {
            index_of_key[entries[i].record_key] = i;
        }
    }

    // A group is keyed by (parent record key, chain, ploidy, carrying traversal).
    //
    // The chain is its boundary pair from the graph. A snarl the decomposition puts in no chain
    // becomes its own group rather than being pooled with every other such snarl under this
    // parent.
    //
    // A group is decoded at one ploidy, and a site decoded at the wrong ploidy would index past
    // its likelihood vector. Children in one chain of one parent can differ in ploidy, since a
    // child's ploidy is the number of the parent's chosen alleles that cross it.
    //
    // At ploidy 1, the parent's chosen traversal that carries the site names the strand the site
    // is on. A group is placed on one strand, and two sites of one chain can be carried by
    // different strands of the parent.
    using GroupKey = std::tuple<size_t, size_t, size_t, int>;
    auto group_key = [&](const Entry& e, const Entry& parent) {
        const size_t chain = e.chain_key != 0
                                 ? e.chain_key
                                 : numeric_limits<size_t>::max() - e.record_key;
        const int carrying = e.ploidy == 1 ? relation_to(e, parent).carrying_trav : -1;
        return GroupKey(e.parent_record_key, chain, e.ploidy, carrying);
    };
    // Group every live site of this level with its parent. A site that cannot be grouped is
    // decoded alone, so that it is still chosen and phased.
    map<GroupKey, vector<size_t>> by_parent;
    vector<DecodeChain> ungrouped;
    for (uint32_t idx : this_level) {
        const Entry& e = entries[idx];
        if (e.retracted) {
            continue;
        }
        auto par = index_of_key.find(e.parent_record_key);
        if (e.parent_record_key == 0) {
            ++model.counters().grp_no_parent;
        } else if (par == index_of_key.end()) {
            ++model.counters().grp_no_entry;
        } else {
            by_parent[group_key(e, entries[par->second])].push_back(idx);
            continue;
        }
        ++model.counters().grp_vetoed;
        DecodeChain alone;
        alone.sites.push_back(idx);
        ungrouped.push_back(std::move(alone));
    }

    vector<DecodeChain> chains;
    for (auto& [key, members] : by_parent) {
        const size_t parent = index_of_key[std::get<0>(key)];
        chains.push_back(parent_group(parent, std::get<2>(key), members, pinned_phase,
                                      placements));
        grouped_sites += chains.back().sites.size();
        ++grouped_groups;
    }
    for (DecodeChain& alone : ungrouped) {
        chains.push_back(std::move(alone));
    }
    return chains;
}

LinkageCollector::DecodeChain LinkageCollector::parent_group(
        size_t parent, size_t group_ploidy, vector<size_t>& members,
        const unordered_map<size_t, PinnedPhase>& pinned_phase,
        unordered_map<size_t, NestedPlacement>& placements) {
    DecodeChain group;
    // The parent joins its group when their ploidies match, as for a haploid child of a haploid
    // parent. A haploid child of a diploid parent is decoded without it, and the parent reaches
    // the group through the entering message instead.
    const bool parent_in_group = entries[parent].ploidy == group_ploidy;
    if (parent_in_group) {
        group.sites.push_back(parent);
    }
    std::ranges::sort(members, [this](size_t a, size_t b) { return before_in_position(a, b); });
    group.sites.insert(group.sites.end(), members.begin(), members.end());

    const size_t parent_key = entries[parent].record_key;
    auto pin = pinned_phase.find(parent_key);
    if (pin == pinned_phase.end()) {
        // No PhaseCall for the parent, so nothing ties this group's orientation to it.
        ++model.counters().group_parent_unpinned;
        return group;
    }
    ++model.counters().group_parent_pinned;
    // The phase set a group belongs to is its parent's. A group is a unit of decoding and a phase
    // block a unit of meaning, so taking the phase set from the group's own first site would
    // fragment the output.
    group.phase_set = pin->second.phase_set;

    const size_t m = entries.panel_size();
    if (group_ploidy == 1) {
        group.entering = haploid_group_message(parent, parent_in_group, members, pin->second,
                                               placements);
    } else if (pin->second.first < m && pin->second.second < m) {
        // A point mass at the parent's pair, which a parent strand naming no haplotype cannot
        // give.
        group.entering.assign(m * m, 0.0);
        group.entering[pin->second.first * m + pin->second.second] = 1.0;
    }
    return group;
}

vector<double> LinkageCollector::haploid_group_message(
        size_t parent, bool parent_in_group, const vector<size_t>& members,
        const PinnedPhase& pin, unordered_map<size_t, NestedPlacement>& placements) {
    // The group sits on the parent strand that `PinnedPhase::strand_carrying` names. Every member is carried
    // by the same parent traversal, which is part of the group key, so the first member stands
    // for all.
    const Entry& child = entries[members.front()];
    const int carrying = relation_to(child, entries[parent]).carrying_trav;
    const int strand = pin.strand_carrying(carrying);
    // A haploid parent records its haplotype in the slot its own `nested_strand` names:
    // `hap_second` on strand 1 and `hap_first` otherwise.
    const int slot = pin.ploidy == 1 ? pin.nested_strand : strand;
    const size_t hap = slot == 1 ? pin.second : pin.first;
    // A haploid parent has one haplotype, so the child needs no strand of its own to find it.
    const bool have_hap = pin.ploidy == 1 || strand >= 0;

    for (size_t idx : members) {
        placements[entries[idx].record_key] =
            NestedPlacement{strand, strand >= 0 ? pin.order_arbitrary : false, have_hap};
    }
    if (strand >= 0) {
        model.counters().nest_strand += members.size();
    } else if (have_hap) {
        model.counters().nest_one_hap += members.size();
    } else if (carrying == -2) {
        model.counters().nest_both += members.size();
    } else {
        model.counters().nest_unreadable += members.size();
    }

    // Give no message where the parent's haplotype does not pass through the child, since it
    // names no allele there.
    const size_t m = entries.panel_size();
    if (!(have_hap && hap < m && (int)entries.haplotype_allele(child, hap) >= 0)) {
        return {};
    }
    // A group that starts at its parent starts at the parent's haplotype, as a point mass. A group
    // without its parent starts at its first child, so the haplotype is carried through one
    // transition to it, which leaves the child's reads free to overrule it.
    const double rho = parent_in_group
                           ? 0.0
                           : model.switch_probability(position_gap(
                                 entries[parent].position, entries[parent].unpositioned, true,
                                 child.position, child.unpositioned));
    vector message(m, rho / (double)m);
    message[hap] += 1.0 - rho;
    return message;
}

vector<LinkageModel::Site> LinkageCollector::model_sites(
        const vector<size_t>& indices, size_t level,
        const unordered_map<size_t, PinnedPhase>& pinned_phase) const {
    vector<LinkageModel::Site> sites;
    sites.reserve(indices.size());
    for (size_t k = 0; k < indices.size(); ++k) {
        const Entry& e = entries[indices[k]];
        LinkageModel::Site s;
        s.position = e.position;
        s.unpositioned = e.unpositioned;
        s.num_alleles = e.num_alleles;
        s.ploidy = e.ploidy;
        s.freq_prior = e.freq_prior;
        const size_t n_gt = e.ploidy == 1
                                ? (size_t)e.num_alleles
                                : (size_t)e.num_alleles * ((size_t)e.num_alleles + 1) / 2;
        s.genotype_ln_likelihood.reserve(n_gt);
        for (size_t g = 0; g < n_gt; ++g) {
            s.genotype_ln_likelihood.push_back((double)entries.likelihood_at(e, g));
        }
        s.haplotype_allele.reserve(entries.panel_size());
        for (size_t h = 0; h < entries.panel_size(); ++h) {
            s.haplotype_allele.push_back((int)entries.haplotype_allele(e, h));
        }
        if (e.level < level) {
            // Clamp the site with a point-mass emission at its chosen genotype. It still carries
            // transition context for its neighbours, which is why it is in the chain at all, but
            // cannot move. build_emission maps a non-finite entry to zero mass.
            const size_t chosen = e.chosen_genotype();
            for (size_t g = 0; g < s.genotype_ln_likelihood.size(); ++g) {
                s.genotype_ln_likelihood[g] =
                    (g == chosen) ? 0.0 : -numeric_limits<double>::infinity();
            }
            // The one clamped site of a group is its parent, held first.
            s.group_parent = (k == 0);
            if (auto pin = pinned_phase.find(e.record_key); pin != pinned_phase.end()) {
                s.pinned = true;
                s.pin_first = pin->second.first;
                s.pin_second = pin->second.second;
            }
        }
        sites.push_back(std::move(s));
    }
    return sites;
}

void LinkageCollector::choose_genotypes(const vector<size_t>& indices, size_t level,
                                        const vector<vector<double>>& posteriors,
                                        ChainDecode& d) const {
    d.final_genotype.assign(indices.size(), LinkageModel::NO_CONSTRAINT);
    d.best.assign(indices.size(), ChainDecode::NO_POSTERIOR);
    d.best_posterior.assign(indices.size(), 0.0);
    for (size_t t = 0; t < indices.size(); ++t) {
        const Entry& e = entries[indices[t]];
        const vector<double>& post = posteriors[t];
        if (e.level < level) {
            // Clamped. The site was chosen, reported and emitted at its own level, but the
            // phasing is constrained to its genotype too.
            d.final_genotype[t] = e.chosen_genotype();
            continue;
        }
        if (post.empty()) {
            d.final_genotype[t] = e.called_genotype();
            continue;
        }
        size_t best = 0;
        for (size_t g = 1; g < post.size(); ++g) {
            if (post[g] > post[best]) {
                best = g;
            }
        }
        d.final_genotype[t] = best;
        d.best[t] = best;
        d.best_posterior[t] = post[best];
    }
}

LinkageCollector::ChainDecode LinkageCollector::decode_chain(
        const DecodeChain& chain, size_t level,
        const unordered_map<size_t, PinnedPhase>& pinned_phase, bool phase) const {
    ChainDecode d;
    const vector<size_t>& indices = chain.sites;
    if (indices.empty()) {
        return d;
    }
    for (size_t idx : indices) {
        if (entries[idx].level == level) {
            ++d.live_here;
        } else if (pinned_phase.contains(entries[idx].record_key)) {
            ++d.pinned_here;
        }
    }
    if (d.live_here == 0) {
        // Every site is clamped, so there is nothing to decide.
        return d;
    }
    // A one-site chain has nothing to link to, so the model cannot change its genotype. It is
    // still phased, so that it appears in `phasing_out` and the mosaic.

    const size_t chain_ploidy = entries[indices.front()].ploidy;
    for (size_t idx : indices) {
        assert(entries[idx].ploidy == chain_ploidy
               && "a decode chain must be homogeneous in ploidy");
    }

    const vector<LinkageModel::Site> sites = model_sites(indices, level, pinned_phase);

    // The entering message is over ordered pairs for a diploid group and over single haplotypes
    // for a haploid one. Each decode ignores a message of the wrong size for its states.
    const vector<double>* entering = chain.entering.empty() ? nullptr : &chain.entering;
    vector<vector<double>> posteriors;
    if (chain_ploidy == 1) {
        posteriors = model.posteriors(sites, 1, entering);
    } else if (entering != nullptr) {
        model.segment_posteriors(sites, 0, sites.size(), entering, nullptr, posteriors);
    } else {
        posteriors = model.posteriors(sites, 2);
    }
    choose_genotypes(indices, level, posteriors, d);
    vector<vector<double>>().swap(posteriors);

    if (!phase) {
        return d;
    }
    // Phase the genotypes each site ends up with. Phasing the pre-linkage calls would describe a
    // genotype set that never reaches the VCF.
    d.phase = model.phasing(sites, d.final_genotype, chain_ploidy, entering);
    // One phase set per chain, named by its first site's position unless the chain has its
    // parent's. The windows are pinned to each other, so the path is continuous across the chain.
    d.phase_set = chain.phase_set != numeric_limits<size_t>::max() ? chain.phase_set
                                                                    : sites.front().position;
    // Read the alleles the path's haplotypes carry now, so that `sites` need not be kept until
    // the results are applied.
    d.allele_first.assign(d.phase.size(), -1);
    d.allele_second.assign(d.phase.size(), -1);
    for (size_t t = 0; t < indices.size() && t < d.phase.size(); ++t) {
        const LinkageModel::Phase& ph = d.phase[t];
        if (ph.first != LinkageModel::NO_HAPLOTYPE
            && ph.first < sites[t].haplotype_allele.size()) {
            d.allele_first[t] = sites[t].haplotype_allele[ph.first];
        }
        if (ph.second != LinkageModel::NO_HAPLOTYPE
            && ph.second < sites[t].haplotype_allele.size()) {
            d.allele_second[t] = sites[t].haplotype_allele[ph.second];
        }
    }
    return d;
}

vector<LinkageCollector::ChainDecode> LinkageCollector::decode_chains(
        const vector<DecodeChain>& chains, size_t level,
        const unordered_map<size_t, PinnedPhase>& pinned_phase, bool phase) const {
    // Each chain reads only what was fixed before decoding, so the chains are decoded in
    // parallel. Start the longest first, so that the longest decode is not left until last.
    vector<ChainDecode> decoded(chains.size());
    vector<size_t> decode_order(chains.size());
    std::iota(decode_order.begin(), decode_order.end(), 0);
    std::ranges::stable_sort(decode_order, [&](size_t a, size_t b) {
        return chains[a].sites.size() > chains[b].sites.size();
    });
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t order_i = 0; order_i < decode_order.size(); ++order_i) {
        const size_t chain_i = decode_order[order_i];
        decoded[chain_i] = decode_chain(chains[chain_i], level, pinned_phase, phase);
    }
    return decoded;
}

size_t LinkageCollector::apply_chosen_genotypes(const vector<size_t>& indices,
                                                const ChainDecode& d, size_t level) {
    size_t moved = 0;
    for (size_t t = 0; t < indices.size(); ++t) {
        Entry& e = entries[indices[t]];
        if (e.level < level) {
            // Clamped; it was applied at its own level.
            continue;
        }
        if (d.best[t] == ChainDecode::NO_POSTERIOR) {
            moved_quality_by_record.erase(e.record_key);
            continue;
        }
        const size_t before = e.called_genotype();
        const size_t best = d.best[t];
        // At ploidy 1 the genotype index is the allele, and both slots carry it.
        size_t i = best;
        size_t j = best;
        if (e.ploidy != 1) {
            genotype_alleles(best, i, j);
        }
        // A later level clamps the site here instead of reconsidering it, which makes a parent's
        // genotype final before any of its children is called.
        e.final_i = (uint16_t)i;
        e.final_j = (uint16_t)j;
        if (best != before) {
            // Keep the posterior, since the record's GQ is computed from it.
            MovedQuality& q = moved_quality_by_record[e.record_key];
            q.posterior = d.best_posterior[t];
            q.direct.explained_share = e.explained_share;
            q.direct.gq_factor = e.gq_factor;
            q.direct.achievable_gap = e.achievable_gap;
            ++moved;
        } else {
            // The site is unmoved in this resolution, whatever an earlier linkage pass concluded.
            moved_quality_by_record.erase(e.record_key);
        }
    }
    return moved;
}

/// Put a diploid site's chosen genotype `want` on strands, from the alleles `a` and `b` that the
/// haplotypes the path chose carry there, or -1 for none. Set `pc.allele_first` and
/// `pc.allele_second` to compact alleles.
///
/// Where neither haplotype's allele orders the pair, the genotype's own order is written and
/// `pc.order_arbitrary` is set. The phase is then unsupported rather than wrong.
static void order_on_strands(LinkageCollector::PhaseCall& pc, size_t want, int a, int b) {
    size_t i = 0;
    size_t j = 0;
    genotype_alleles(want, i, j);
    if (a >= 0 && b >= 0 && LinkageModel::genotype_index((size_t)a, (size_t)b) == want) {
        pc.allele_first = (size_t)a;
        pc.allele_second = (size_t)b;
    } else if (a >= 0 && a == b && i != j) {
        // Both haplotypes carry the same allele, and a mutation gave one strand another. Either
        // strand could have mutated, so nothing orders the pair.
        pc.allele_first = i;
        pc.allele_second = j;
        pc.order_arbitrary = true;
    } else if (a >= 0 && ((size_t)a == i || (size_t)a == j)) {
        pc.allele_first = (size_t)a;
        pc.allele_second = ((size_t)a == i) ? j : i;
    } else if (b >= 0 && ((size_t)b == i || (size_t)b == j)) {
        pc.allele_second = (size_t)b;
        pc.allele_first = ((size_t)b == i) ? j : i;
    } else {
        // Neither panel haplotype spells either called allele. Sorted order is written, which
        // pairs allele_first with hap_first only by accident.
        pc.allele_first = i;
        pc.allele_second = j;
        pc.order_arbitrary = (i != j);
    }
}

LinkageCollector::PhaseCall LinkageCollector::phase_call_for(
        const Entry& e, const ChainDecode& d, size_t t,
        const unordered_map<size_t, NestedPlacement>& placements) const {
    const LinkageModel::Phase& ph = d.phase[t];
    const size_t want = d.final_genotype[t];
    const int a = d.allele_first[t];
    const int b = d.allele_second[t];

    PhaseCall pc;
    pc.ploidy = e.ploidy;
    pc.level = e.level;
    pc.record_key = e.record_key;
    pc.contig = entries.contig_name(e.contig);
    pc.position = e.position;
    pc.start_node = e.start_node;
    pc.end_node = e.end_node;
    pc.phase_set = d.phase_set;

    // The strand of its parent that a nested ploidy-1 site sits on.
    auto found = placements.find(e.record_key);
    const NestedPlacement placement = found != placements.end() ? found->second
                                                                : NestedPlacement{};
    pc.nested_strand = (int8_t)placement.strand;
    pc.order_arbitrary = placement.order_arbitrary;
    // Fill the slot that `nested_strand` names, since the mosaic reads that slot and treats the
    // other as empty.
    if (!placement.nameable) {
        pc.hap_first = LinkageModel::NO_HAPLOTYPE;
        pc.hap_second = LinkageModel::NO_HAPLOTYPE;
    } else if (placement.strand == 1) {
        pc.hap_first = LinkageModel::NO_HAPLOTYPE;
        pc.hap_second = ph.first;
    } else {
        pc.hap_first = ph.first;
        pc.hap_second = ph.second;
    }

    if (e.ploidy == 1) {
        // The one strand carries the chosen allele.
        pc.allele_first = want;
        pc.allele_second = want;
    } else {
        order_on_strands(pc, want, a, b);
    }
    // A strand whose haplotype carries another allele than the strand's own is a mutation.
    model.counters().mutated_strands += (size_t)(a >= 0 && (size_t)a != pc.allele_first)
                                        + (size_t)(e.ploidy != 1 && b >= 0
                                                   && (size_t)b != pc.allele_second);
    finish_phase_call(pc, e, pc.allele_first, pc.allele_second);
    return pc;
}

void LinkageCollector::append_phase_calls(
        const vector<size_t>& indices, const ChainDecode& d, size_t level,
        const unordered_map<size_t, NestedPlacement>& placements,
        vector<PhaseCall>& phasing_out) const {
    for (size_t t = 0; t < indices.size() && t < d.phase.size(); ++t) {
        const Entry& e = entries[indices[t]];
        if (e.level < level) {
            // Its PhaseCall was emitted at its own level.
            continue;
        }
        phasing_out.push_back(phase_call_for(e, d, t, placements));
    }
}

unordered_set<size_t> LinkageCollector::parent_keys_of(const vector<uint32_t>& this_level) const {
    unordered_set<size_t> parent_keys;
    for (uint32_t idx : this_level) {
        if (!entries[idx].retracted) {
            parent_keys.insert(entries[idx].parent_record_key);
        }
    }
    return parent_keys;
}

void LinkageCollector::report_level(size_t level, size_t grouped_sites,
                                    size_t grouped_groups) const {
    const LinkageCounters& c = model.counters();
    if (level > 0 && (c.pin_applied.load() + c.pin_declined.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage level " << level << ": pins -- "
                  << c.pin_applied.load() << " applied, " << c.pin_declined.load()
                  << " REFUSED (the pinned pair cannot spell the constrained genotype, so that"
                  << " site's orientation is free); groups whose parent was pinnable: "
                  << c.group_parent_pinned.load() << ", not pinnable: "
                  << c.group_parent_unpinned.load() << std::endl;
    }

    // Report grouping only when something was declined.
    if (level > 0
        && (c.grp_no_parent.load() + c.grp_no_entry.load() + c.grp_vetoed.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage level " << level << ": grouping declines so far -- "
                  << c.grp_no_parent.load() << " sites with no parent key, "
                  << c.grp_no_entry.load() << " whose parent has no live entry; "
                  << c.grp_vetoed.load() << " chains kept ungrouped in total" << std::endl;
    }

    if (grouped_groups > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] linkage level " << level << ": decoded " << grouped_sites
                  << " sites in " << grouped_groups << " per-parent groups instead of the contig"
                  << " chain" << std::endl;
    }

    if (level > 0
        && (c.nest_strand.load() + c.nest_one_hap.load() + c.nest_both.load()
            + c.nest_unreadable.load()) > 0) {
#pragma omp critical (cerr)
        std::cerr << "[vg call] nested strands: " << c.nest_strand.load()
                  << " on one of a diploid parent's two strands, " << c.nest_one_hap.load()
                  << " on a haploid parent's single haplotype (no strand to choose), "
                  << c.nest_both.load() << " carried on both parent strands, "
                  << c.nest_unreadable.load()
                  << " whose parent's chosen pair could not be read -- the last two name no"
                  << " haplotype" << std::endl;
    }
}

size_t LinkageCollector::resolve_level(size_t level, bool last, vector<PhaseCall>* phasing_out) {
    if (!model.active() || entries.empty()) {
        return 0;
    }

    // This level's entries, in append order.
    const vector<uint32_t>& this_level = entries.at_level(level);

    // The lookups below are built over the parents' keys only, so that their cost follows the
    // level rather than the genome.
    const unordered_set<size_t> parent_keys =
        level > 0 ? parent_keys_of(this_level) : unordered_set<size_t>();
    // Below the top level, sites are grouped by parent only when an earlier level was phased.
    const bool have_pins = phasing_out != nullptr && level > 0 && !phasing_out->empty();
    // The phase of each parent, rebuilt from `phasing_out` on each call, since the last call sorts
    // `phasing_out` in place.
    unordered_map<size_t, PinnedPhase> pinned_phase;
    if (have_pins) {
        pinned_phase = pinned_phases(parent_keys, *phasing_out);
    }

    // Default every site of this level to its own call, so that a site no chain decodes still has
    // a coherent genotype for a later level to clamp.
    for (uint32_t idx : this_level) {
        Entry& e = entries[idx];
        e.final_i = e.called_i;
        e.final_j = e.ploidy == 1 ? e.called_i : e.called_j;
    }

    // At the top level, decode each contig's linkage chains. Below it, decode each nested chain
    // on its own, conditioned on its parent's chosen state, so the cost grows with the number of
    // children rather than with the length of the contig.
    unordered_map<size_t, NestedPlacement> placements;
    size_t grouped_sites = 0;
    size_t grouped_groups = 0;
    vector<DecodeChain> chains;
    if (level == 0) {
        chains = top_level_chains(this_level);
    } else if (have_pins) {
        chains = parent_groups(this_level, parent_keys, pinned_phase, placements, grouped_sites,
                               grouped_groups);
    }

    vector<ChainDecode> decoded =
        decode_chains(chains, level, pinned_phase, phasing_out != nullptr);

    // Apply the results one chain at a time, in chain order, so the entries,
    // `moved_quality_by_record` and `phasing_out` do not depend on the threads. Write the
    // per-chain lines to the log in one piece, since a level can have tens of thousands of chains.
    size_t moved = 0;
    string chain_lines;
    for (size_t chain_i = 0; chain_i < chains.size(); ++chain_i) {
        const vector<size_t>& indices = chains[chain_i].sites;
        ChainDecode& d = decoded[chain_i];
        if (indices.empty() || d.live_here == 0) {
            continue;
        }
        if (level > 0) {
            chain_lines += "[vg call] linkage level " + std::to_string(level) + ": chain decodes " +
                std::to_string(indices.size()) + " sites for " + std::to_string(d.live_here) +
                " of its own, " + std::to_string(d.pinned_here) + " pinned\n";
        }
        moved += apply_chosen_genotypes(indices, d, level);
        if (phasing_out != nullptr) {
            append_phase_calls(indices, d, level, placements, *phasing_out);
        }
        // Free the applied decode.
        d = ChainDecode();
    }
    if (!chain_lines.empty()) {
#pragma omp critical (cerr)
        std::cerr << chain_lines << std::flush;
    }
    report_level(level, grouped_sites, grouped_groups);

    if (phasing_out != nullptr && last) {
        // Sort `phasing_out` into reference order for the mosaic writer, since nested sites were
        // appended after their parents' chains. Sort stably, because a site re-rendered at the
        // linkage pass has two PhaseCalls with the same record key. Readers keep the last one for
        // a key, which must be the later level's, and `phasing_out` is appended one level at a
        // time.
        std::ranges::stable_sort(*phasing_out, [](const PhaseCall& a, const PhaseCall& b) {
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
