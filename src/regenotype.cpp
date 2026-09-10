#include "regenotype.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "anchor.hpp"

namespace vg {

using std::exp;
using std::fabs;
using std::log;
using std::max;
using std::min;

/// The mixture, from a pair of unnormalised slot weights and the two alleles' responsibilities.
///
/// One definition for both sides of the delta, and that is what makes `tau = 0` exact rather than
/// approximately exact: there, `(s0, s1)` IS `(w_0, w_1)`, so the corrected and uncorrected
/// expressions are the same arithmetic on the same values and cancel bit for bit. Written as a
/// normalised ratio rather than assuming the weights sum to one, because `site_slot_weights`
/// normalises in floating point and `w_0 + w_1` is not exactly 1.
static inline double mixture_of(double s0, double s1, double rel0, double rel1) {
    const double z = s0 + s1;
    if (!(z > 0.0)) {
        return 0.0;
    }
    return (s0 * rel0 + s1 * rel1) / z;
}

/// The tilt a read's strand log-odds applies, as one factor and which side it lands on.
///
/// `Lambda` reaches into the hundreds before tempering, so `exp(tau * Lambda)` overflows for a
/// perfectly ordinary read. Scaling by `e^-max(x, 0)` keeps every intermediate finite and gives
/// the same pair of normalised weights -- and then both branches want `exp(-|x|)`, which is the
/// whole reason this can be hoisted.
struct ReadTilt {
    double factor = 1.0;   ///< e^-|x|
    bool positive = true;  ///< whether x >= 0, so which slot the factor multiplies
};

/// Depends on the READ and not on the candidate genotype, so it is computed once per read rather
/// than twice per (read, genotype). At three genotypes a site that is one `exp` where there were
/// six -- five saved of six, against ~98 M of them on chr20 -- and it is the same value computed
/// once instead of six times, so every sum is bit for bit what it was.
static inline ReadTilt read_tilt(double x) {
    ReadTilt t;
    // `-0.0 >= 0.0` is true and `fabs(-0.0)` is `0.0`, so the zero case takes the same branch and
    // the same factor either way round. That matters: at `tau = 0` every read is this case, and
    // the byte-identity gate rests on it.
    t.positive = x >= 0.0;
    t.factor = exp(-fabs(x));
    return t;
}

static inline void tilted_weights(double w0, double w1, const ReadTilt& t, double& s0, double& s1) {
    if (t.positive) {
        s0 = w0;
        s1 = w1 * t.factor;
    } else {
        s0 = w0 * t.factor;
        s1 = w1;
    }
}

void merge_counters(const RegenotypeCounters& from, RegenotypeCounters& into) {
    into.reads_with_lambda += from.reads_with_lambda;
    into.reads_multi_block += from.reads_multi_block;
    into.reads_singleton += from.reads_singleton;
    into.sites_considered += from.sites_considered;
    into.sites_corrected += from.sites_corrected;
    into.sites_would_move += from.sites_would_move;
    into.moved_hom_to_het += from.moved_hom_to_het;
    into.moved_het_to_hom += from.moved_het_to_hom;
    into.moved_het_to_het += from.moved_het_to_het;
    into.order_reversed += from.order_reversed;
    // `fitted_temper` and the calibration vectors are deliberately not merged: they describe one
    // fit, done once before any of this is parallel, and adding them up would say nothing.
}

double site_read_log_odds(double q0, double p) {
    // The read reports the true slot with probability p and a coin flip otherwise, so neither side
    // can reach 0 while p < 1 and one mismapped read cannot carry an unbounded opinion.
    const double half = 0.5 * (1.0 - p);
    const double a = p * q0 + half;
    const double b = p * (1.0 - q0) + half;
    if (!(a > 0.0) || !(b > 0.0)) {
        return 0.0;
    }
    return log(a / b);
}

void site_own_log_odds(const PhaseSite& site, bool flipped, unordered_map<uint64_t, double>& out) {
    out.clear();
    out.reserve(site.read_key.size() * 2);
    const double sign = flipped ? -1.0 : 1.0;
    for (size_t i = 0; i < site.read_key.size(); ++i) {
        const double l = sign * site_read_log_odds((double)site.q0[i], (double)site.p[i]);
        // Deduplicated by key rather than summed: paired mates share a read name and so a
        // read_key, and a fragment lies on one haplotype, so it is one observation. Keeping the
        // first matches what accumulate_lambda adds, which is what makes the subtraction exact.
        out.emplace(site.read_key[i], l);
    }
}

void accumulate_lambda(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                       LambdaTable& out, RegenotypeCounters& counters) {
    unordered_map<uint64_t, double> own;
    for (const PhaseSite& site : sites) {
        site_own_log_odds(site, flipped.count(site.record_key) != 0, own);
        for (const auto& kv : own) {
            ReadLambda& rl = out[kv.first];
            if (rl.sites == 0) {
                rl.phase_set = site.phase_set;
            } else if (rl.phase_set != site.phase_set) {
                // Two blocks label their strands independently, so there is no sum to take.
                rl.multi_block = true;
            }
            rl.lambda += kv.second;
            ++rl.sites;
        }
    }
    for (const auto& kv : out) {
        ++counters.reads_with_lambda;
        counters.reads_multi_block += kv.second.multi_block ? 1 : 0;
        counters.reads_singleton += kv.second.sites <= 1 ? 1 : 0;
    }
}

double fit_temper(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                  const LambdaTable& lambda, const RegenotypeParams& params,
                  RegenotypeCounters& counters) {
    // Every (read, site) pair where the rest of the read has an opinion and this site has an
    // observation to check it against. `loo` is the prediction, `own` the observation.
    struct Obs { double abs_loo; bool agree; };
    vector<Obs> obs;
    unordered_map<uint64_t, double> own;
    for (const PhaseSite& site : sites) {
        site_own_log_odds(site, flipped.count(site.record_key) != 0, own);
        for (const auto& kv : own) {
            auto found = lambda.find(kv.first);
            if (found == lambda.end() || found->second.multi_block) {
                continue;
            }
            const double loo = found->second.lambda - kv.second;
            if (!(fabs(loo) > 1e-9) || !(fabs(kv.second) > 1e-9)) {
                continue;
            }
            obs.push_back({fabs(loo), (loo > 0.0) == (kv.second > 0.0)});
        }
    }
    if (obs.size() < params.fit_min_per_bin) {
        counters.fitted_temper = 0.0;
        return 0.0;
    }
    // Equal-count bins over |Lambda|, so a long tail does not get one bin to itself.
    std::sort(obs.begin(), obs.end(), [](const Obs& a, const Obs& b) {
        if (a.abs_loo != b.abs_loo) return a.abs_loo < b.abs_loo;
        return (int)a.agree < (int)b.agree;
    });
    const size_t n_bins = max<size_t>(1, params.fit_bins);
    const size_t per_bin = max<size_t>(params.fit_min_per_bin, obs.size() / n_bins);
    struct Bin { double mean_abs; double observed; size_t n; };
    vector<Bin> bins;
    for (size_t start = 0; start < obs.size(); start += per_bin) {
        const size_t stop = min(obs.size(), start + per_bin);
        double sum = 0.0;
        size_t agree = 0;
        for (size_t i = start; i < stop; ++i) {
            sum += obs[i].abs_loo;
            agree += obs[i].agree ? 1 : 0;
        }
        const size_t n = stop - start;
        bins.push_back({sum / (double)n, (double)agree / (double)n, n});
    }
    // tau such that sigmoid(tau * |Lambda|) matches the observed agreement, weighted by bin size.
    // A grid then a refinement: the objective is smooth and one-dimensional, and a closed form
    // would have to assume the link is exactly logistic, which is the thing being tested.
    auto cost = [&](double tau) {
        double c = 0.0;
        for (const Bin& b : bins) {
            const double x = tau * b.mean_abs;
            const double predicted = 1.0 / (1.0 + exp(-x));
            const double d = predicted - b.observed;
            c += (double)b.n * d * d;
        }
        return c;
    };
    double best = 0.0, best_cost = cost(0.0);
    for (int i = 1; i <= 200; ++i) {
        const double tau = i * 0.01;
        const double c = cost(tau);
        if (c < best_cost) {
            best_cost = c;
            best = tau;
        }
    }
    for (double step = 0.005; step > 1e-4; step *= 0.5) {
        for (int s = -1; s <= 1; s += 2) {
            const double tau = best + s * step;
            if (tau < 0.0) {
                continue;
            }
            const double c = cost(tau);
            if (c < best_cost) {
                best_cost = c;
                best = tau;
            }
        }
    }
    counters.fit_abs_lambda.clear();
    counters.fit_observed.clear();
    counters.fit_predicted.clear();
    counters.fit_count.clear();
    for (const Bin& b : bins) {
        counters.fit_abs_lambda.push_back(b.mean_abs);
        counters.fit_observed.push_back(b.observed);
        counters.fit_predicted.push_back(1.0 / (1.0 + exp(-best * b.mean_abs)));
        counters.fit_count.push_back(b.n);
    }
    counters.fitted_temper = best;
    return best;
}

bool phase_aware_correction(const PhaseReadEvidence& ev, const LambdaTable& lambda,
                            const unordered_map<uint64_t, double>& own, double temper,
                            const RegenotypeParams& params, map<vector<int>, double>& gl,
                            RegenotypeCounters& counters) {
    if (gl.empty() || ev.n_alleles == 0 || ev.num_reads() == 0) {
        return false;
    }
    ++counters.sites_considered;

    // Per read, the leave-one-out log-odds and the tilt it implies -- once, rather than per
    // candidate genotype. The tilt is the expensive half and the genotype loop below does not
    // change it.
    vector<double> loo(ev.num_reads(), 0.0);
    vector<ReadTilt> tilt(ev.num_reads());
    bool any_opinion = false;
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        auto found = lambda.find(ev.read_key[r]);
        if (found == lambda.end() || found->second.multi_block) {
            continue;
        }
        double v = found->second.lambda;
        auto mine = own.find(ev.read_key[r]);
        if (mine != own.end()) {
            v -= mine->second;
        }
        if (params.shuffle) {
            // Sign randomised from the read key, keeping |Lambda|: the peakedness distribution is
            // preserved exactly and only the phase content is destroyed. Deterministic, so the
            // arm is reproducible.
            v = (ev.read_key[r] & 1ULL) ? -fabs(v) : fabs(v);
        }
        loo[r] = v;
        if (fabs(v) > 1e-9) {
            any_opinion = true;
        }
    }
    if (!any_opinion) {
        // Every read here spans nothing else in its block. Provably inert; skip the arithmetic
        // rather than compute a column of zeroes.
        //
        // Before the tilts, not after: on short reads most reads span one site, so this is the
        // common case there, and filling a tilt column first would spend an `exp` a read to reach
        // the same return. On ONT it is rare enough not to show up in a timing.
        return false;
    }
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        tilt[r] = read_tilt(temper * loo[r]);
    }

    // The sweep's own argmax, to say afterwards whether the correction moved it.
    const vector<int>* before = nullptr;
    double before_ll = -std::numeric_limits<double>::infinity();
    for (const auto& kv : gl) {
        if (kv.second > before_ll) {
            before_ll = kv.second;
            before = &kv.first;
        }
    }

    bool reversed_called = false;
    size_t corrected = 0;
    for (auto& kv : gl) {
        const vector<int>& g = kv.first;
        if (g.size() != 2) {
            // Ploidy 1 has one slot, so `pi` is 1 and the correction is identically zero.
            continue;
        }
        const int a = g[0], b = g[1];
        if (a == b) {
            // A homozygote's mixture collapses to rel(r, a) whatever the weights are, so the
            // correction is exactly zero. Skipped rather than computed: the algebra says zero and
            // evaluating it would only add rounding to a quantity that must not move.
            continue;
        }
        if (a < 0 || b < 0 || (size_t)a >= ev.n_alleles || (size_t)b >= ev.n_alleles) {
            continue;
        }
        // Weights follow the ALLELE, not the slot, so the uncorrected mixture is invariant under
        // reordering the pair. Pinning them to the slot instead would make a read that spans
        // nothing else order-dependent, and the max below would then pick up a difference with no
        // phase in it -- worth up to 0.37 ln over three reads at an SV-sized length ratio.
        const vector<double> w = site_slot_weights(ev.allele_length, ev.n_alleles,
                                                   ev.mean_read_length, ev.length_weighted,
                                                   vector<int>{a, b});
        if (w.size() != 2) {
            continue;
        }
        double s_fwd = 0.0, s_rev = 0.0, s_w = 0.0;
        for (size_t r = 0; r < ev.num_reads(); ++r) {
            const double e = (double)ev.mismap[r];
            const double ra = (double)ev.rel_at(r, (size_t)a);
            const double rb = (double)ev.rel_at(r, (size_t)b);

            // All three mixtures are evaluated with `ra` first and `rb` second, so the only thing
            // that differs between them is which weight each allele gets. That is not a style
            // choice: at `tau = 0` the three weight pairs are all `(w[0], w[1])`, and identical
            // operands in an identical expression give an identical result, so the correction is
            // exactly zero rather than zero to within whatever the optimiser did to two separate
            // accumulator chains. Written the obvious way -- `(s0, s1, rb, ra)` for the reverse --
            // the two disagreed in the last bits at 957 of chr20's 170,060 sites.
            //
            // The reverse order puts allele `b` on strand 0, so its slot weights come back
            // swapped and are then passed swapped, which restores `ra` to the first position.
            double s0, s1;
            tilted_weights(w[0], w[1], tilt[r], s0, s1);
            s_fwd += log((1.0 - e) * mixture_of(s0, s1, ra, rb) + e);
            tilted_weights(w[1], w[0], tilt[r], s0, s1);
            s_rev += log((1.0 - e) * mixture_of(s1, s0, ra, rb) + e);
            s_w += log((1.0 - e) * mixture_of(w[0], w[1], ra, rb) + e);
        }
        const double correction = max(s_fwd, s_rev) - s_w;
        // Only for the genotype the site is actually called at. At a multi-allelic site most
        // pairs are ones no haplotype carries, and which order "fits" such a pair better is
        // arbitrary -- counting those made this read 87,440 of chr20's 163,396 corrected sites,
        // which is not a statement about the chain, it is a statement about noise.
        if (before != nullptr && g == *before) {
            reversed_called = s_rev > s_fwd;
        }
        kv.second += correction;
        ++corrected;
    }
    if (corrected == 0) {
        return false;
    }
    ++counters.sites_corrected;
    counters.order_reversed += reversed_called ? 1 : 0;

    const vector<int>* after = nullptr;
    double after_ll = -std::numeric_limits<double>::infinity();
    for (const auto& kv : gl) {
        if (kv.second > after_ll) {
            after_ll = kv.second;
            after = &kv.first;
        }
    }
    if (before == nullptr || after == nullptr || *before == *after) {
        return false;
    }
    ++counters.sites_would_move;
    auto is_hom = [](const vector<int>& g) {
        for (size_t i = 1; i < g.size(); ++i) {
            if (g[i] != g[0]) {
                return false;
            }
        }
        return true;
    };
    const bool hom_before = is_hom(*before), hom_after = is_hom(*after);
    if (hom_before && !hom_after) {
        ++counters.moved_hom_to_het;
    } else if (!hom_before && hom_after) {
        ++counters.moved_het_to_hom;
    } else if (!hom_before && !hom_after) {
        ++counters.moved_het_to_het;
    }
    return true;
}

}
