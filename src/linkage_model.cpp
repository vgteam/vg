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

//------------------------------------------------------------------------------
// LinkageCollector

void LinkageCollector::record(const string& contig, size_t position, size_t num_alleles,
                              const vector<double>& genotype_ln_likelihood,
                              const vector<int>& haplotype_allele,
                              size_t called_i, size_t called_j, size_t record_key) {
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

vector<LinkageCollector::Change> LinkageCollector::resolve() const {
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

    for (auto& indices : by_contig) {
        if (indices.size() < 2) {
            continue;
        }
        // Position, then the site's own key. Sites arrive in whatever order the threads finished,
        // and two records can share a position, so sorting on position alone leaves their relative
        // order down to scheduling -- which would make the output depend on --threads. The key is
        // derived from the snarl ID, so it is a property of the site rather than of the run.
        sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
            if (entries[a].position != entries[b].position) {
                return entries[a].position < entries[b].position;
            }
            return entries[a].record_key < entries[b].record_key;
        });

        vector<LinkageModel::Site> sites;
        sites.reserve(indices.size());
        for (size_t idx : indices) {
            const Entry& e = entries[idx];
            LinkageModel::Site s;
            s.position = e.position;
            s.num_alleles = e.num_alleles;
            size_t n_gt = (size_t)e.num_alleles * ((size_t)e.num_alleles + 1) / 2;
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

        vector<vector<double>> posteriors = model.posteriors(sites);
        for (size_t t = 0; t < indices.size(); ++t) {
            const vector<double>& post = posteriors[t];
            if (post.empty()) {
                continue;
            }
            size_t best = 0;
            for (size_t g = 1; g < post.size(); ++g) {
                if (post[g] > post[best]) {
                    best = g;
                }
            }
            const Entry& e = entries[indices[t]];
            if (best == LinkageModel::genotype_index(e.called_i, e.called_j)) {
                continue;
            }
            // Decode the genotype index back to its allele pair.
            size_t j = 0;
            while (LinkageModel::genotype_index(0, j) <= best) {
                ++j;
            }
            --j;
            size_t i = best - (j * (j + 1) / 2);
            Change c;
            c.record_key = e.record_key;
            c.contig = contig_names[e.contig];
            c.position = e.position;
            c.called_i = e.called_i;
            c.called_j = e.called_j;
            c.allele_i = i;
            c.allele_j = j;
            c.posterior = post[best];
            changes.push_back(c);
        }
    }
    return changes;
}

}
