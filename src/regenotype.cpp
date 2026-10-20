#include "regenotype.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <omp.h>

#include "parallel_sort.hpp"


namespace vg {

using std::exp;
using std::fabs;
using std::log;
using std::max;
using std::min;

/// Return the mixture of the two alleles' relative likelihoods under a pair of unnormalised slot
/// weights.
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

/// Return the weighting that strand log-odds `x` apply.
///
/// It depends on the read and not on the candidate genotype, so it is computed once per read.
static inline ReadTilt read_tilt(double x) {
    ReadTilt t;
    // `-0.0 >= 0.0` is true and `fabs(-0.0)` is `0.0`, so a zero takes the same branch and factor
    // whatever its sign. At tau = 0 every read is this case.
    t.positive = x >= 0.0;
    t.factor = exp(-fabs(x));
    return t;
}

/// Apply a read's weighting to the slot weights `w0` and `w1`, giving unnormalised weights `s0` and
/// `s1`.
static inline void tilted_weights(double w0, double w1, const ReadTilt& t, double& s0, double& s1) {
    if (t.positive) {
        s0 = w0;
        s1 = w1 * t.factor;
    } else {
        s0 = w0 * t.factor;
        s1 = w1;
    }
}

/// Return the logistic function of `x`, computed so that it does not overflow.
///
/// At |x| in the hundreds, `exp(-x)` is infinite for negative x.
static inline double stable_sigmoid(double x) {
    return x >= 0.0 ? 1.0 / (1.0 + exp(-x)) : exp(x) / (1.0 + exp(x));
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
    // `calibration` and the calibration table describe the one fit done before the parallel
    // region, so they are not added up.
}

double StrandCalibration::log_odds(double lambda) const {
    const double x = temper * lambda;
    if (!(x != 0.0)) {
        // Return exactly 0 at tau = 0, whatever the sign of Lambda (including -0.0).
        return 0.0;
    }
    double c = ceiling;
    if (c >= 1.0) {
        c = 1.0;
    } else if (c <= 0.0) {
        c = 0.0;
    }
    const double sig = stable_sigmoid(x);
    double p = c * sig + 0.5 * (1.0 - c);
    // Keep `p` away from 0 and 1. With c == 1 and tau * Lambda around 60, `p` rounds to 1 and the
    // logit would be infinite.
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
        // The site's rows hold one per read key (see `merge_mates`). If a key were repeated, we
        // would keep its first row, which is also the one `accumulate_lambda` adds, so the
        // subtraction would still cancel.
        out.emplace(site.read_key[i], l);
    }
}

void accumulate_lambda(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                       LambdaTable& out, RegenotypeCounters& counters) {
    unordered_map<uint64_t, double> own;
    for (const PhaseSite& site : sites) {
        site_own_log_odds(site, flipped.contains(site.record_key), own);
        for (const auto& [read, log_odds] : own) {
            ReadLambda& rl = out[read];
            if (rl.sites == 0) {
                rl.phase_set = site.phase_set;
            } else if (rl.phase_set != site.phase_set) {
                // Two phase sets label their strands independently, so there is no sum to take.
                rl.multi_phase_set = true;
            }
            rl.lambda += log_odds;
            ++rl.sites;
        }
    }
    for (const auto& [read, rl] : out) {
        ++counters.reads_with_lambda;
        counters.reads_multi_phase_set += rl.multi_phase_set ? 1 : 0;
        counters.reads_singleton += rl.sites <= 1 ? 1 : 0;
    }
}

namespace {

/// One read at one site, where the read's other sites give it a strand and its allele here can
/// check that strand.
struct StrandCheck {
    /// |Lambda| without this site, which predicts the strand.
    double abs_loo;
    /// Whether the allele here points to the strand Lambda points to.
    bool agree;
};

/// Return every strand check from `sites`, sorted by `abs_loo` and then `agree`.
///
/// A whole genome has over 10^8 of them, so they are collected and sorted on several threads.
/// Sorting on every field makes the result independent of the order they were collected in.
vector<StrandCheck> collect_strand_checks(const vector<PhaseSite>& sites,
                                          const unordered_set<size_t>& flipped,
                                          const LambdaTable& lambda) {
    vector<vector<StrandCheck>> collected(max(1, omp_get_max_threads()));
#pragma omp parallel
    {
        vector<StrandCheck>& mine = collected[omp_get_thread_num()];
        unordered_map<uint64_t, double> own;
#pragma omp for schedule(dynamic, 1024)
        for (const PhaseSite& site : sites) {
            site_own_log_odds(site, flipped.contains(site.record_key), own);
            for (const auto& [read, own_log_odds] : own) {
                auto found = lambda.find(read);
                if (found == lambda.end() || !read_strand_usable(found->second, site.phase_set)) {
                    continue;
                }
                const double loo = found->second.lambda - own_log_odds;
                if (!(fabs(loo) > 1e-9) || !(fabs(own_log_odds) > 1e-9)) {
                    continue;
                }
                mine.push_back({fabs(loo), (loo > 0.0) == (own_log_odds > 0.0)});
            }
        }
    }
    size_t total = 0;
    for (const vector<StrandCheck>& part : collected) {
        total += part.size();
    }
    vector<StrandCheck> checks;
    checks.reserve(total);
    for (vector<StrandCheck>& part : collected) {
        checks.insert(checks.end(), part.begin(), part.end());
        vector<StrandCheck>().swap(part);
    }
    parallel_sort(checks.begin(), checks.end(), [](const StrandCheck& a, const StrandCheck& b) {
        if (a.abs_loo != b.abs_loo) {
            return a.abs_loo < b.abs_loo;
        }
        return (int)a.agree < (int)b.agree;
    });
    return checks;
}

/// A range of |Lambda| in the calibration table.
struct CalibrationBin {
    /// Mean |Lambda| of the bin's strand checks.
    double mean_abs;
    /// The fraction of the bin's strand checks that agree.
    double observed;
    /// The number of strand checks in the bin.
    size_t n;
};

/// Group sorted strand checks into bins of equal count, so a long tail does not get one bin to
/// itself.
vector<CalibrationBin> bin_strand_checks(const vector<StrandCheck>& checks,
                                         const RegenotypeParams& params) {
    const size_t n_bins = max<size_t>(1, params.fit_bins);
    const size_t per_bin = max<size_t>(params.fit_min_per_bin, checks.size() / n_bins);
    vector<CalibrationBin> bins;
    for (size_t start = 0; start < checks.size(); start += per_bin) {
        const size_t stop = min(checks.size(), start + per_bin);
        double sum = 0.0;
        size_t agree = 0;
        for (size_t i = start; i < stop; ++i) {
            sum += checks[i].abs_loo;
            agree += checks[i].agree ? 1 : 0;
        }
        const size_t n = stop - start;
        bins.emplace_back(sum / (double)n, (double)agree / (double)n, n);
    }
    return bins;
}

/// Return the agreement a bin's mean |Lambda| predicts at temper `tau` and ceiling `ceil`.
double predicted_agreement(const CalibrationBin& bin, double tau, double ceil) {
    const double sig = stable_sigmoid(tau * bin.mean_abs);
    return ceil * sig + 0.5 * (1.0 - ceil);
}

/// Return the squared difference between predicted and observed agreement, weighted by bin size
/// and summed over the bins.
double calibration_cost(const vector<CalibrationBin>& bins, double tau, double ceil) {
    double acc = 0.0;
    for (const CalibrationBin& b : bins) {
        const double d = predicted_agreement(b, tau, ceil) - b.observed;
        acc += (double)b.n * d * d;
    }
    return acc;
}

}

StrandCalibration fit_calibration(const vector<PhaseSite>& sites,
                                  const unordered_set<size_t>& flipped, const LambdaTable& lambda,
                                  const RegenotypeParams& params, RegenotypeCounters& counters) {
    const vector<StrandCheck> checks = collect_strand_checks(sites, flipped, lambda);
    if (checks.size() < params.fit_min_per_bin) {
        counters.calibration = StrandCalibration{0.0, 1.0};
        return counters.calibration;
    }
    const vector<CalibrationBin> bins = bin_strand_checks(checks, params);
    // Search a grid over the temper alone, in steps of 0.01. `RegenotypeParams::ceiling` sets the
    // ceiling.
    const double best_ceiling = min(1.0, params.ceiling);
    double best = 0.0;
    double best_cost = calibration_cost(bins, 0.0, best_ceiling);
    for (int i = 1; i <= 200; ++i) {
        const double tau = i * 0.01;
        if (calibration_cost(bins, tau, best_ceiling) < best_cost) {
            best_cost = calibration_cost(bins, tau, best_ceiling);
            best = tau;
        }
    }
    counters.fit_abs_lambda.clear();
    counters.fit_observed.clear();
    counters.fit_predicted.clear();
    counters.fit_count.clear();
    for (const CalibrationBin& b : bins) {
        counters.fit_abs_lambda.push_back(b.mean_abs);
        counters.fit_observed.push_back(b.observed);
        counters.fit_predicted.push_back(predicted_agreement(b, best, best_ceiling));
        counters.fit_count.push_back(b.n);
    }
    counters.calibration = StrandCalibration{best, best_ceiling};
    return counters.calibration;
}

bool read_strand_usable(const ReadLambda& read, size_t phase_set) {
    if (read.multi_phase_set) {
        return false;
    }
    return phase_set == NO_PHASE_SET || read.phase_set == phase_set;
}

/// Fill `loo` with each read's strand log-odds leaving this site out, and return whether any is
/// nonzero.
///
/// A read whose strand is not usable in `phase_set` gets 0.
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
        if (auto mine = own.find(ev.read_key[r]); mine != own.end()) {
            v -= mine->second;
        }
        if (params.shuffle) {
            // Take the sign from the read key, keeping |Lambda|, so that the phase information is
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

/// Return the genotype with the highest likelihood in `gl`, the first of any tie, or null if
/// none is above negative infinity.
static const vector<int>* best_genotype(const map<vector<int>, double>& gl) {
    const vector<int>* best = nullptr;
    double best_ll = -std::numeric_limits<double>::infinity();
    for (const auto& [genotype, ll] : gl) {
        if (ll > best_ll) {
            best_ll = ll;
            best = &genotype;
        }
    }
    return best;
}

bool haploid_inclusion_correction(const PhaseReadEvidence& ev, const LambdaTable& lambda,
                                  size_t phase_set,
                                  const unordered_map<uint64_t, double>& own,
                                  const StrandCalibration& calibration, int strand_sign,
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
    vector incl(ev.num_reads(), 1.0);
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        const double x = strand_sign * calibration.log_odds(loo[r]);
        incl[r] = x >= 0.0 ? 1.0 : exp(x);
    }

    const vector<int>* before = best_genotype(gl);

    size_t corrected = 0;
    for (auto& [g, ll] : gl) {
        if (g.size() != 1) {
            // The inclusion correction applies only at ploidy 1.
            continue;
        }
        const int a = g[0];
        if (a < 0 || (size_t)a >= ev.n_alleles) {
            continue;
        }
        double s_incl = 0.0;
        double s_base = 0.0;
        for (size_t r = 0; r < ev.num_reads(); ++r) {
            const double e = ev.mismap[r];
            const double ra = ev.rel_at(r, (size_t)a);
            // `(1 - e) * incl * rel + e + (1 - e) * (1 - incl)`, written so that at incl == 1 it is
            // the same expression as the baseline term below, and the two cancel exactly.
            s_incl += log((1.0 - e) * (incl[r] * ra + 1.0 - incl[r]) + e);
            s_base += log((1.0 - e) * (1.0 * ra + 1.0 - 1.0) + e);
        }
        ll += s_incl - s_base;
        ++corrected;
    }
    if (corrected == 0) {
        return false;
    }
    if (const vector<int>* after = best_genotype(gl);
        before == nullptr || after == nullptr || *before == *after) {
        return false;
    }
    ++counters.haploid_would_move;
    return true;
}

namespace {

/// The read terms of one heterozygous genotype, summed over a site's reads, under three
/// weightings.
struct PairReadTerms {
    /// With strand weights for allele `a` on strand 0.
    double forward = 0.0;
    /// With strand weights for allele `b` on strand 0.
    double reverse = 0.0;
    /// With the allele-length weights alone.
    double unphased = 0.0;
};

/// Sum the read terms of heterozygous genotype {a, b} under each weighting.
///
/// `w` holds the allele-length weights of `a` and `b`, and `tilt` each read's weighting.
PairReadTerms pair_read_terms(const PhaseReadEvidence& ev, const vector<ReadTilt>& tilt, int a,
                              int b, const vector<double>& w) {
    PairReadTerms terms;
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        const double e = ev.mismap[r];
        const double ra = ev.rel_at(r, (size_t)a);
        const double rb = ev.rel_at(r, (size_t)b);

        // All three mixtures take `ra` first and `rb` second, so they differ only in the weight
        // each allele gets. At tau = 0 the three weight pairs are the same, so the same
        // expression gives the same result and the correction is exactly zero.
        //
        // The reverse order puts allele `b` on strand 0, so its slot weights come back swapped
        // and are passed swapped, which keeps `ra` first.
        double s0;
        double s1;
        tilted_weights(w[0], w[1], tilt[r], s0, s1);
        terms.forward += log((1.0 - e) * mixture_of(s0, s1, ra, rb) + e);
        tilted_weights(w[1], w[0], tilt[r], s0, s1);
        terms.reverse += log((1.0 - e) * mixture_of(s1, s0, ra, rb) + e);
        terms.unphased += log((1.0 - e) * mixture_of(w[0], w[1], ra, rb) + e);
    }
    return terms;
}

/// Return whether every allele of `g` is the same.
bool is_homozygous(const vector<int>& g) {
    for (size_t i = 1; i < g.size(); ++i) {
        if (g[i] != g[0]) {
            return false;
        }
    }
    return true;
}

/// Count, in `counters`, the kind of move from best genotype `before` to `after`.
void count_move(const vector<int>& before, const vector<int>& after,
                RegenotypeCounters& counters) {
    const bool hom_before = is_homozygous(before);
    const bool hom_after = is_homozygous(after);
    if (hom_before && !hom_after) {
        ++counters.moved_hom_to_het;
    } else if (!hom_before && hom_after) {
        ++counters.moved_het_to_hom;
    } else if (!hom_before && !hom_after) {
        ++counters.moved_het_to_het;
    }
}

}

bool phase_aware_correction(const PhaseReadEvidence& ev, const LambdaTable& lambda,
                            size_t phase_set, int strand0_allele,
                            const unordered_map<uint64_t, double>& own,
                            const StrandCalibration& calibration, const RegenotypeParams& params,
                            map<vector<int>, double>& gl, RegenotypeCounters& counters) {
    if (gl.empty() || ev.n_alleles == 0 || ev.num_reads() == 0) {
        return false;
    }
    ++counters.sites_considered;

    // Compute each read's log-odds leaving this site out, and the weighting they imply, once
    // rather than once per candidate genotype.
    vector<double> loo;
    if (!read_loo(ev, lambda, phase_set, own, params, loo)) {
        // No read here spans another site of its phase set, so the correction is zero.
        return false;
    }
    vector<ReadTilt> tilt(ev.num_reads());
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        tilt[r] = read_tilt(calibration.log_odds(loo[r]));
    }

    const vector<int>* before = best_genotype(gl);

    bool reversed_called = false;
    size_t corrected = 0;
    for (auto& [g, ll] : gl) {
        if (g.size() != 2) {
            // Ploidy 1 has one slot, so `pi` is 1 and the correction is identically zero.
            continue;
        }
        const int a = g[0];
        const int b = g[1];
        if (a == b) {
            // A homozygote's mixture is rel(r, a) whatever the weights are, so the correction is
            // zero.
            continue;
        }
        if (a < 0 || b < 0 || (size_t)a >= ev.n_alleles || (size_t)b >= ev.n_alleles) {
            continue;
        }
        // Weights follow the allele, not the slot, so the unphased mixture does not change when
        // the pair is reordered. Otherwise taking the larger of the two orders would favour one
        // for a reason unrelated to phase.
        const vector<double> w = allele_length_weights(ev.allele_length, ev.n_alleles,
                                                       ev.mean_read_length, ev.length_weighted,
                                                       vector<int>{a, b});
        if (w.size() != 2) {
            continue;
        }
        const PairReadTerms terms = pair_read_terms(ev, tilt, a, b, w);
        const double correction = max(terms.forward, terms.reverse) - terms.unphased;
        // Compare orders only at the called genotype. At a site with several alleles most pairs
        // are carried by no haplotype, and which order fits such a pair better means nothing.
        if (before != nullptr && g == *before) {
            if (strand0_allele == a) {
                reversed_called = terms.reverse > terms.forward;
            } else if (strand0_allele == b) {
                reversed_called = terms.forward > terms.reverse;
            }
        }
        ll += correction;
        ++corrected;
    }
    if (corrected == 0) {
        return false;
    }
    ++counters.sites_corrected;
    counters.order_reversed += reversed_called ? 1 : 0;

    const vector<int>* after = best_genotype(gl);
    if (before == nullptr || after == nullptr || *before == *after) {
        return false;
    }
    ++counters.sites_would_move;
    count_move(*before, *after, counters);
    return true;
}

}
