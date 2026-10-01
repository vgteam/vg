#include "regenotype.hpp"

#include <algorithm>
#include <cmath>
#include <limits>


namespace vg {

using std::exp;
using std::fabs;
using std::log;
using std::max;
using std::min;

/// The mixture, from a pair of unnormalised slot weights and the two alleles' relative
/// likelihoods.
///
/// Both sides of the correction use this function, so at tau = 0, where `(s0, s1)` is
/// `(w_0, w_1)`, the two sides are the same arithmetic on the same values and cancel exactly.
/// The ratio is normalised here because the weights need not sum to exactly 1.
static inline double mixture_of(double s0, double s1, double rel0, double rel1) {
    const double z = s0 + s1;
    if (!(z > 0.0)) {
        return 0.0;
    }
    return (s0 * rel0 + s1 * rel1) / z;
}

/// The weighting a read's strand log-odds apply, as one factor and the side it applies to.
///
/// Lambda can reach the hundreds, where `exp(tau * Lambda)` would overflow. Scaling by
/// `e^-max(x, 0)` keeps every value finite and gives the same normalised weights, and both signs
/// then need only `exp(-|x|)`.
struct ReadTilt {
    double factor = 1.0;   ///< e^-|x|
    bool positive = true;  ///< whether x >= 0, so which slot the factor multiplies
};

/// Depends on the read and not on the candidate genotype, so it is computed once per read.
static inline ReadTilt read_tilt(double x) {
    ReadTilt t;
    // `-0.0 >= 0.0` is true and `fabs(-0.0)` is `0.0`, so a zero takes the same branch and factor
    // whatever its sign. At tau = 0 every read is this case.
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
    into.reads_multi_phase_set += from.reads_multi_phase_set;
    into.reads_singleton += from.reads_singleton;
    into.sites_considered += from.sites_considered;
    into.sites_corrected += from.sites_corrected;
    into.sites_would_move += from.sites_would_move;
    into.moved_hom_to_het += from.moved_hom_to_het;
    into.moved_het_to_hom += from.moved_het_to_hom;
    into.moved_het_to_het += from.moved_het_to_het;
    into.order_reversed += from.order_reversed;
    into.haploid_sites += from.haploid_sites;
    into.haploid_would_move += from.haploid_would_move;
    // `fitted_temper` and the calibration vectors describe the one fit done before the parallel
    // region, so they are not added up.
}

double calibrated_log_odds(double lambda, double temper, double ceiling) {
    const double x = temper * lambda;
    if (!(x != 0.0)) {
    // Exactly 0 at tau = 0, whatever the sign of Lambda (including -0.0).
        return 0.0;
    }
    const double c = ceiling >= 1.0 ? 1.0 : (ceiling <= 0.0 ? 0.0 : ceiling);
    // sigmoid without the overflow: at |x| in the hundreds `exp(-x)` is inf on one side.
    const double sig = x >= 0.0 ? 1.0 / (1.0 + exp(-x)) : exp(x) / (1.0 + exp(x));
    double p = c * sig + 0.5 * (1.0 - c);
    // Kept away from 0 and 1: with c == 1 and tau * Lambda around 60, `p` rounds to 1 and the logit
    // would be infinite.
    const double eps = 1e-12;
    p = min(1.0 - eps, max(eps, p));
    return log(p / (1.0 - p));
}

double site_read_log_odds(double q0, double c) {
    // The read reports the true slot with probability c and a coin flip otherwise, so neither side
    // reaches 0 while c < 1, and a mismapped read contributes little.
    const double half = 0.5 * (1.0 - c);
    const double a = c * q0 + half;
    const double b = c * (1.0 - q0) + half;
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
        const double l = sign * site_read_log_odds((double)site.q0[i], (double)site.c[i]);
        // The site's rows hold one per read key (see `merge_mates`). Were a key repeated, the
        // first row would be kept, as `accumulate_lambda` adds, so the subtraction still cancels.
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
                // Two phase sets label their strands independently, so there is no sum to take.
                rl.multi_phase_set = true;
            }
            rl.lambda += kv.second;
            ++rl.sites;
        }
    }
    for (const auto& kv : out) {
        ++counters.reads_with_lambda;
        counters.reads_multi_phase_set += kv.second.multi_phase_set ? 1 : 0;
        counters.reads_singleton += kv.second.sites <= 1 ? 1 : 0;
    }
}

void fit_calibration(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                     const LambdaTable& lambda, const RegenotypeParams& params,
                     double& temper, double& ceiling, RegenotypeCounters& counters) {
    // Every (read, site) pair where the read's other sites give a strand and its allele at this
    // site can check it. `loo` is the prediction (leaving this site out), `own` the observation.
    struct Obs { double abs_loo; bool agree; };
    vector<Obs> obs;
    unordered_map<uint64_t, double> own;
    for (const PhaseSite& site : sites) {
        site_own_log_odds(site, flipped.count(site.record_key) != 0, own);
        for (const auto& kv : own) {
            auto found = lambda.find(kv.first);
            if (found == lambda.end() || !read_strand_usable(found->second, site.phase_set)) {
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
        counters.fitted_ceiling = 1.0;
        temper = 0.0;
        ceiling = 1.0;
        return;
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
    // The squared difference, weighted by bin size, between the predicted and observed agreement
    // in each bin, for a given temper and ceiling. The ceiling is not fitted: `--regeno-ceiling`
    // sets it, and the temper is fitted against it.
    auto cost = [&](double tau, double ceil) {
        double acc = 0.0;
        for (const Bin& b : bins) {
            const double x = tau * b.mean_abs;
            const double sig = x >= 0.0 ? 1.0 / (1.0 + exp(-x)) : exp(x) / (1.0 + exp(x));
            const double predicted = ceil * sig + 0.5 * (1.0 - ceil);
            const double d = predicted - b.observed;
            acc += (double)b.n * d * d;
        }
        return acc;
    };
    // A grid over the temper alone, in steps of 0.01, against the fixed ceiling.
    const double best_ceiling = min(1.0, params.ceiling);
    double best = 0.0, best_cost = cost(0.0, best_ceiling);
    for (int i = 1; i <= 200; ++i) {
        const double tau = i * 0.01;
        if (cost(tau, best_ceiling) < best_cost) {
            best_cost = cost(tau, best_ceiling);
            best = tau;
        }
    }
    counters.fit_abs_lambda.clear();
    counters.fit_observed.clear();
    counters.fit_predicted.clear();
    counters.fit_count.clear();
    for (const Bin& b : bins) {
        counters.fit_abs_lambda.push_back(b.mean_abs);
        counters.fit_observed.push_back(b.observed);
        const double sx = best * b.mean_abs;
        const double sg = sx >= 0.0 ? 1.0 / (1.0 + exp(-sx)) : exp(sx) / (1.0 + exp(sx));
        counters.fit_predicted.push_back(best_ceiling * sg + 0.5 * (1.0 - best_ceiling));
        counters.fit_count.push_back(b.n);
    }
    counters.fitted_temper = best;
    counters.fitted_ceiling = best_ceiling;
    temper = best;
    ceiling = best_ceiling;
}

bool read_strand_usable(const ReadLambda& read, size_t phase_set) {
    if (read.multi_phase_set) {
        return false;
    }
    return phase_set == NO_PHASE_SET || read.phase_set == phase_set;
}

/// A read's strand log-odds leaving this site out, shared by both corrections. 0 for a read whose
/// strand is not usable in `phase_set`.
static bool read_loo(const PhaseReadEvidence& ev, const LambdaTable& lambda, size_t phase_set,
                     const unordered_map<uint64_t, double>& own, const RegenotypeParams& params,
                     vector<double>& loo) {
    loo.assign(ev.num_reads(), 0.0);
    bool any = false;
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        auto found = lambda.find(ev.read_key[r]);
        if (found == lambda.end() || !read_strand_usable(found->second, phase_set)) {
            continue;
        }
        double v = found->second.lambda;
        auto mine = own.find(ev.read_key[r]);
        if (mine != own.end()) {
            v -= mine->second;
        }
        if (params.shuffle) {
            // Sign taken from the read key, keeping |Lambda|, so that the phase information is
            // removed deterministically.
            v = (ev.read_key[r] & 1ULL) ? -fabs(v) : fabs(v);
        }
        loo[r] = v;
        if (fabs(v) > 1e-9) {
            any = true;
        }
    }
    return any;
}

bool haploid_inclusion_correction(const PhaseReadEvidence& ev, const LambdaTable& lambda,
                                  size_t phase_set,
                                  const unordered_map<uint64_t, double>& own, double temper,
                                  double ceiling, int strand_sign,
                                  const RegenotypeParams& params, map<vector<int>, double>& gl,
                                  RegenotypeCounters& counters) {
    if (gl.empty() || ev.n_alleles == 0 || ev.num_reads() == 0 || strand_sign == 0) {
        return false;
    }
    vector<double> loo;
    if (!read_loo(ev, lambda, phase_set, own, params, loo)) {
        return false;
    }
    ++counters.haploid_sites;

    // Per read: how much of it belongs to this chain's strand. Capped at 1, so a read the phase
    // places here is whole and one placed on the other strand is discounted by its odds.
    vector<double> incl(ev.num_reads(), 1.0);
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        const double x = strand_sign * calibrated_log_odds(loo[r], temper, ceiling);
        incl[r] = x >= 0.0 ? 1.0 : exp(x);
    }

    const vector<int>* before = nullptr;
    double before_ll = -std::numeric_limits<double>::infinity();
    for (const auto& kv : gl) {
        if (kv.second > before_ll) {
            before_ll = kv.second;
            before = &kv.first;
        }
    }

    size_t corrected = 0;
    for (auto& kv : gl) {
        const vector<int>& g = kv.first;
        if (g.size() != 1) {
            // This is the ploidy-1 space. A diploid entry here is not ours to touch.
            continue;
        }
        const int a = g[0];
        if (a < 0 || (size_t)a >= ev.n_alleles) {
            continue;
        }
        double s_incl = 0.0, s_base = 0.0;
        for (size_t r = 0; r < ev.num_reads(); ++r) {
            const double e = (double)ev.mismap[r];
            const double ra = (double)ev.rel_at(r, (size_t)a);
            // `(1 - e) * incl * rel + e + (1 - e) * (1 - incl)`, written so that at incl == 1 it is
            // the same expression as the baseline term below, and the two cancel exactly.
            s_incl += log((1.0 - e) * (incl[r] * ra + 1.0 - incl[r]) + e);
            s_base += log((1.0 - e) * (1.0 * ra + 1.0 - 1.0) + e);
        }
        kv.second += s_incl - s_base;
        ++corrected;
    }
    if (corrected == 0) {
        return false;
    }
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
    ++counters.haploid_would_move;
    return true;
}

bool phase_aware_correction(const PhaseReadEvidence& ev, const LambdaTable& lambda,
                            size_t phase_set,
                            const unordered_map<uint64_t, double>& own, double temper,
                            double ceiling, const RegenotypeParams& params,
                            map<vector<int>, double>& gl, RegenotypeCounters& counters) {
    if (gl.empty() || ev.n_alleles == 0 || ev.num_reads() == 0) {
        return false;
    }
    ++counters.sites_considered;

    // For each read, the log-odds leaving this site out and the weighting they imply, computed
    // once rather than once per candidate genotype.
    vector<double> loo;
    const bool any_opinion = read_loo(ev, lambda, phase_set, own, params, loo);
    vector<ReadTilt> tilt(ev.num_reads());
    if (!any_opinion) {
        // No read here spans another site of its phase set, so the correction is zero. Checked
        // before the weightings are computed, since on short reads this is the common case.
        return false;
    }
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        tilt[r] = read_tilt(calibrated_log_odds(loo[r], temper, ceiling));
    }

    // The sweep's best genotype, to report whether the correction changes it.
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
            // A homozygote's mixture is rel(r, a) whatever the weights are, so the correction is
            // zero, and it is skipped rather than computed.
            continue;
        }
        if (a < 0 || b < 0 || (size_t)a >= ev.n_alleles || (size_t)b >= ev.n_alleles) {
            continue;
        }
        // Weights follow the allele, not the slot, so the uncorrected mixture does not change when
        // the pair is reordered; otherwise the larger of the two orders would favour one for a
        // reason unrelated to phase.
        const vector<double> w = allele_length_weights(ev.allele_length, ev.n_alleles,
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

            // All three mixtures are evaluated with `ra` first and `rb` second, so they differ only
            // in the weight each allele gets. At tau = 0 the three weight pairs are the same, so
            // the same expression gives the same result and the correction is exactly zero.
            //
            // The reverse order puts allele `b` on strand 0, so its slot weights come back
            // swapped and are passed swapped, which keeps `ra` first.
            double s0, s1;
            tilted_weights(w[0], w[1], tilt[r], s0, s1);
            s_fwd += log((1.0 - e) * mixture_of(s0, s1, ra, rb) + e);
            tilted_weights(w[1], w[0], tilt[r], s0, s1);
            s_rev += log((1.0 - e) * mixture_of(s1, s0, ra, rb) + e);
            s_w += log((1.0 - e) * mixture_of(w[0], w[1], ra, rb) + e);
        }
        const double correction = max(s_fwd, s_rev) - s_w;
        // Only at the called genotype. At a site with several alleles most pairs are carried by no
        // haplotype, and which order fits such a pair better means nothing.
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
