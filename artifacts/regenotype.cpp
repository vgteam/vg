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
using std::pair;

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

/// Return whether strand log-odds `x` say anything about a read's strand.
static inline bool has_strand(double x) {
    return fabs(x) > 1e-9;
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
    into.sites_split_runs += from.sites_split_runs;
    into.haploid_split_runs += from.haploid_split_runs;
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
            if (site.run == NO_RUN || rl.many_runs) {
                continue;
            }
            const auto part = std::ranges::find_if(rl.run_parts, [&site](const auto& p) {
                return p.run == NO_RUN || p.run == site.run;
            });
            if (part == rl.run_parts.end()) {
                rl.many_runs = true;
                continue;
            }
            part->run = site.run;
            ++part->sites;
            part->lambda += log_odds;
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
                if (!has_strand(loo) || !has_strand(own_log_odds)) {
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

RunCorrelationTable::RunCorrelationTable(const vector<PhaseSite>& sites) {
    for (const PhaseSite& site : sites) {
        if (site.run != NO_RUN) {
            correlation.try_emplace(RunKey {site.phase_set, site.run}, site.run_correlation);
        }
    }
}

double RunCorrelationTable::agreement(const RunKey& a, const RunKey& b) const {
    if (a.phase_set != b.phase_set || a.run == NO_RUN || b.run == NO_RUN) {
        // Nothing orients a run against a run of another phase set, or a site in no run.
        return 0.5;
    }
    double product = 1.0;
    const size_t last = max(a.run, b.run);
    for (size_t run = min(a.run, b.run) + 1; run <= last; ++run) {
        const auto found = correlation.find(RunKey {a.phase_set, run});
        product *= found == correlation.end() ? 0.0 : found->second;
        if (product <= 1e-12) {
            // Every correlation is in [0, 1], so the product can only shrink, and stopping here
            // moves the agreement by at most 5e-13.
            break;
        }
    }
    return (1.0 + product) / 2.0;
}

ReadStrand read_strand(const ReadLambda& read, size_t site_run, bool has_own, double own_log_odds) {
    ReadStrand strand;
    strand.run.phase_set = read.phase_set;
    if (read.many_runs) {
        // The read's sites lie in more than two runs, and `run_parts` holds only two of them, so it
        // has no strand.
        return strand;
    }
    const size_t own_run = has_own ? site_run : NO_RUN;
    size_t kept = 0;
    for (const ReadLambda::RunPart& part : read.run_parts) {
        const bool own_here = own_run != NO_RUN && part.run == own_run;
        if (part.sites - (own_here ? 1 : 0) == 0) {
            continue;
        }
        ++kept;
        const double log_odds = part.lambda - (own_here ? own_log_odds : 0.0);
        const size_t kept_run = strand.run.run;
        if (kept_run == NO_RUN || part.run == site_run
            || (kept_run != site_run && fabs(log_odds) > fabs(strand.log_odds))) {
            strand.log_odds = log_odds;
            strand.run.run = part.run;
        }
    }
    if (kept < 2) {
        // The read's other sites lie in one run, or none, so its whole sum is in one frame.
        strand.log_odds = read.lambda - (has_own ? own_log_odds : 0.0);
    }
    return strand;
}

double reoriented_log_odds(double y, double agreement) {
    // logit(p sigma(y) + (1 - p) sigma(-y)), written with e^-|y| so that it cannot overflow, and so
    // that at p = 1/2 the two sides are the same arithmetic and the result is exactly 0.
    const double z = exp(-fabs(y));
    const double toward = agreement + (1.0 - agreement) * z;
    const double away = (1.0 - agreement) + agreement * z;
    const double shrunk = log(toward / away);
    return y >= 0.0 ? shrunk : -shrunk;
}

bool read_strand_usable(const ReadLambda& read, size_t phase_set) {
    if (read.multi_phase_set) {
        return false;
    }
    return phase_set == NO_PHASE_SET || read.phase_set == phase_set;
}

namespace {

/// A table of each read's strand at one site, leaving the site out.
struct SiteStrandTable {
    /// Each read's strand, with log-odds 0 for a read whose strand is not usable at the site.
    vector<ReadStrand> of_read;
    /// Whether any read has a strand.
    bool any = false;
};

/// Return each read's strand at a site of run `site` (see `read_strand`), leaving out the site's
/// own contribution, `own`.
SiteStrandTable site_strands(const PhaseReadEvidence& ev, const StrandModel& strands,
                             const RunKey& site, const unordered_map<uint64_t, double>& own) {
    SiteStrandTable reads;
    reads.of_read.assign(ev.num_reads(), ReadStrand {});
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        const auto found = strands.lambda.find(ev.read_key[r]);
        if (found == strands.lambda.end() || !read_strand_usable(found->second, site.phase_set)) {
            continue;
        }
        const auto mine = own.find(ev.read_key[r]);
        const bool has_own = mine != own.end();
        ReadStrand& strand = reads.of_read[r];
        strand = read_strand(found->second, site.run, has_own, has_own ? mine->second : 0.0);
        if (strands.params.shuffle) {
            // Take the sign from the read key, keeping |Lambda|, so that the phase information is
            // removed deterministically.
            const double size = fabs(strand.log_odds);
            strand.log_odds = (ev.read_key[r] & 1ULL) ? -size : size;
        }
        reads.any = reads.any || has_strand(strand.log_odds);
    }
    return reads;
}

/// Return the run whose frame a site's heterozygotes are scored in: the site's own run if a read
/// with a strand comes from it, and otherwise the run that most such reads come from.
RunKey frame_run(const SiteStrandTable& reads, const RunKey& site) {
    map<RunKey, size_t> count;
    for (const ReadStrand& strand : reads.of_read) {
        if (has_strand(strand.log_odds)) {
            ++count[strand.run];
        }
    }
    if (count.empty() || count.contains(site)) {
        return site;
    }
    return std::ranges::max_element(count, [](const auto& x, const auto& y) {
               return x.second < y.second;
           })->first;
}

/// A division of a site's reads into groups by the run their strands come from.
struct ReadGrouping {
    /// Each read's group. Group 0, the frame group, holds the reads of the frame's run and the
    /// reads with no strand.
    vector<size_t> of_read;
    /// For each group, the probability that the panel orients its run against the frame's as
    /// reads would.
    vector<double> agreement {1.0};
};

/// Group a site's reads by the run their strands come from, in the frame of run `frame`.
ReadGrouping group_reads(const SiteStrandTable& reads, const RunKey& frame,
                          const RunCorrelationTable& correlations) {
    ReadGrouping grouping;
    grouping.of_read.assign(reads.of_read.size(), 0);
    vector<RunKey> group_run {frame};
    for (size_t r = 0; r < reads.of_read.size(); ++r) {
        const ReadStrand& strand = reads.of_read[r];
        if (!has_strand(strand.log_odds) || strand.run == frame) {
            continue;
        }
        const auto known = std::ranges::find(group_run, strand.run);
        grouping.of_read[r] = known - group_run.begin();
        if (known == group_run.end()) {
            group_run.push_back(strand.run);
            grouping.agreement.push_back(correlations.agreement(frame, strand.run));
        }
    }
    return grouping;
}

/// Return ln(p e^same + (1 - p) e^opposite), the term of a group of reads whose strands are
/// oriented as the frame group's are with probability `p`, at least 1/2, and the other way
/// otherwise.
///
/// Equal terms, as at temper 0, give `same` exactly, since p + (1 - p) is exactly 1 for p in
/// [1/2, 1].
double oriented_sum(double same, double opposite, double p) {
    const double larger = same >= opposite ? p : 1.0 - p;
    const double smaller = same >= opposite ? 1.0 - p : p;
    if (larger == 0.0) {
        // At p = 1 only `same` is possible, and `exp` might underflow.
        return min(same, opposite) + log(smaller);
    }
    return max(same, opposite) + log(larger + smaller * exp(-fabs(same - opposite)));
}

/// Return the frame group's term in `same` plus each other group's `oriented_sum`, where `same` and
/// `opposite` hold each group's term with its strands oriented as the frame group's and the other
/// way.
double grouped_sum(const vector<double>& same, const vector<double>& opposite,
                   const ReadGrouping& grouping) {
    double total = same[0];
    for (size_t k = 1; k < grouping.agreement.size(); ++k) {
        total += oriented_sum(same[k], opposite[k], grouping.agreement[k]);
    }
    return total;
}

/// Return the sum of a term over the groups, in group order.
double over_groups(const vector<double>& term) {
    double total = term[0];
    for (size_t k = 1; k < term.size(); ++k) {
        total += term[k];
    }
    return total;
}

/// Return the genotype with the highest likelihood in `gl`, the first of any tie, or null if
/// none is above negative infinity.
const vector<int>* best_genotype(const map<vector<int>, double>& gl) {
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

/// Return the correction to allele `a` of a nested chain at ploidy 1, given how much of each read
/// belongs to the chain's strand, `incl`, and how much would if its group were oriented the other
/// way, `incl_reversed`.
double haploid_allele_correction(const PhaseReadEvidence& ev, int a, const vector<double>& incl,
                                 const vector<double>& incl_reversed,
                                 const ReadGrouping& grouping) {
    const size_t groups = grouping.agreement.size();
    vector kept(groups, 0.0);
    vector reversed(groups, 0.0);
    vector base(groups, 0.0);
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        const double e = ev.mismap[r];
        const double ra = ev.rel_at(r, (size_t)a);
        const size_t k = grouping.of_read[r];
        // `(1 - e) * incl * rel + e + (1 - e) * (1 - incl)`, written so that at incl == 1 it is the
        // same expression as the baseline term below, and the two cancel exactly.
        kept[k] += log((1.0 - e) * (incl[r] * ra + 1.0 - incl[r]) + e);
        if (k > 0) {
            const double other = incl_reversed[r];
            reversed[k] += log((1.0 - e) * (other * ra + 1.0 - other) + e);
        }
        base[k] += log((1.0 - e) * (1.0 * ra + 1.0 - 1.0) + e);
    }
    return grouped_sum(kept, reversed, grouping) - over_groups(base);
}

/// The read terms of one heterozygous genotype, summed over each group of a site's reads, under
/// three weightings.
struct PairReadTerms {
    /// With strand weights for allele `a` on strand 0.
    vector<double> forward;
    /// With strand weights for allele `b` on strand 0.
    vector<double> reverse;
    /// With the allele-length weights alone.
    vector<double> unphased;
};

/// Sum the read terms of heterozygous genotype {a, b} under each weighting.
///
/// `w` holds the allele-length weights of `a` and `b`, and `tilt` each read's weighting.
PairReadTerms pair_read_terms(const PhaseReadEvidence& ev, const vector<ReadTilt>& tilt, int a,
                              int b, const vector<double>& w, const ReadGrouping& grouping) {
    PairReadTerms terms;
    const size_t groups = grouping.agreement.size();
    terms.forward.assign(groups, 0.0);
    terms.reverse.assign(groups, 0.0);
    terms.unphased.assign(groups, 0.0);
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        const double e = ev.mismap[r];
        const double ra = ev.rel_at(r, (size_t)a);
        const double rb = ev.rel_at(r, (size_t)b);
        const size_t k = grouping.of_read[r];

        // All three mixtures take `ra` first and `rb` second, so they differ only in the weight
        // each allele gets. At tau = 0 the three weight pairs are the same, so the same
        // expression gives the same result and the correction is exactly zero.
        //
        // The reverse order puts allele `b` on strand 0, so its slot weights come back swapped
        // and are passed swapped, which keeps `ra` first.
        double s0;
        double s1;
        tilted_weights(w[0], w[1], tilt[r], s0, s1);
        terms.forward[k] += log((1.0 - e) * mixture_of(s0, s1, ra, rb) + e);
        tilted_weights(w[1], w[0], tilt[r], s0, s1);
        terms.reverse[k] += log((1.0 - e) * mixture_of(s1, s0, ra, rb) + e);
        terms.unphased[k] += log((1.0 - e) * mixture_of(w[0], w[1], ra, rb) + e);
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


/// Return whether the reads prefer the other order of a site's two alleles `a` and `b` to the one
/// that puts `strand0_allele` on strand 0, where `forward` scores `a` on strand 0 and `reverse`
/// scores `b` there. A site whose strand-0 allele is neither prefers no other order.
bool prefers_other_order(int strand0_allele, int a, int b, double forward, double reverse) {
    if (strand0_allele == a) {
        return reverse > forward;
    }
    if (strand0_allele == b) {
        return forward > reverse;
    }
    return false;
}
}

bool haploid_inclusion_correction(const PhaseReadEvidence& ev, const StrandModel& strands,
                                  const RunKey& chain, int strand_sign,
                                  const unordered_map<uint64_t, double>& own,
                                  map<vector<int>, double>& gl, RegenotypeCounters& counters) {
    if (gl.empty() || ev.n_alleles == 0 || ev.num_reads() == 0 || strand_sign == 0) {
        return false;
    }
    const SiteStrandTable reads = site_strands(ev, strands, chain, own);
    if (!reads.any) {
        return false;
    }
    ++counters.haploid_sites;
    // The chain's strand is named in its run, so that run's reads keep `strand_sign`, whether or
    // not any has a strand.
    const bool run_known = chain.phase_set != NO_PHASE_SET && chain.run != NO_RUN;
    const RunKey frame = run_known ? chain : frame_run(reads, chain);
    const ReadGrouping grouping = group_reads(reads, frame, strands.correlations);
    counters.haploid_split_runs += grouping.agreement.size() > 1 ? 1 : 0;

    // Per read: how much of it belongs to this chain's strand, and how much would if its group
    // were oriented the other way. Capped at 1, so a read the phase places here is whole and one
    // placed on the other strand is discounted by its odds.
    vector incl(ev.num_reads(), 1.0);
    vector incl_reversed(ev.num_reads(), 1.0);
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        const double x = strand_sign * strands.calibration.log_odds(reads.of_read[r].log_odds);
        incl[r] = x >= 0.0 ? 1.0 : exp(x);
        incl_reversed[r] = x <= 0.0 ? 1.0 : exp(-x);
    }

    const vector<int>* before = best_genotype(gl);
    size_t corrected = 0;
    for (auto& [g, ll] : gl) {
        // The inclusion correction applies only at ploidy 1.
        if (g.size() == 1 && g[0] >= 0 && (size_t)g[0] < ev.n_alleles) {
            ll += haploid_allele_correction(ev, g[0], incl, incl_reversed, grouping);
            ++corrected;
        }
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

bool phase_aware_correction(const PhaseReadEvidence& ev, const StrandModel& strands,
                            const RunKey& site, int strand0_allele,
                            const unordered_map<uint64_t, double>& own,
                            map<vector<int>, double>& gl, RegenotypeCounters& counters) {
    if (gl.empty() || ev.n_alleles == 0 || ev.num_reads() == 0) {
        return false;
    }
    ++counters.sites_considered;

    // Compute each read's log-odds leaving this site out, and the weighting they imply, once
    // rather than once per candidate genotype.
    const SiteStrandTable reads = site_strands(ev, strands, site, own);
    if (!reads.any) {
        // No read here spans another site of its phase set, so the correction is zero.
        return false;
    }
    vector<ReadTilt> tilt(ev.num_reads());
    for (size_t r = 0; r < ev.num_reads(); ++r) {
        tilt[r] = read_tilt(strands.calibration.log_odds(reads.of_read[r].log_odds));
    }
    const RunKey frame = frame_run(reads, site);
    const ReadGrouping grouping = group_reads(reads, frame, strands.correlations);
    // The site's order is in its own run's frame, so only reads scored in that frame say whether
    // they prefer the other order.
    const bool frame_is_own = site.run == NO_RUN || frame == site;

    const vector<int>* before = best_genotype(gl);
    bool reversed_called = false;
    size_t corrected = 0;
    for (auto& [g, ll] : gl) {
        // Ploidy 1 has one slot, so `pi` is 1 and the correction is identically zero. So is a
        // homozygote's, whose mixture is rel(r, a) whatever the weights are.
        if (g.size() != 2 || g[0] == g[1] || g[0] < 0 || g[1] < 0
            || (size_t)g[0] >= ev.n_alleles || (size_t)g[1] >= ev.n_alleles) {
            continue;
        }
        const int a = g[0];
        const int b = g[1];
        // Weights follow the allele, not the slot, so the unphased mixture does not change when
        // the pair is reordered. Otherwise taking the larger of the two orders would favour one
        // for a reason unrelated to phase.
        const vector<double> w = allele_length_weights(ev.allele_length, ev.n_alleles,
                                                       ev.mean_read_length, ev.length_weighted,
                                                       vector<int>{a, b});
        if (w.size() != 2) {
            continue;
        }
        const PairReadTerms terms = pair_read_terms(ev, tilt, a, b, w, grouping);
        // Under each assignment of the frame group, every other group takes the same assignment
        // with the probability that read phasing orients its run as reads would, and the other
        // otherwise.
        const double forward = grouped_sum(terms.forward, terms.reverse, grouping);
        const double reverse = grouped_sum(terms.reverse, terms.forward, grouping);
        // Compare orders only at the called genotype. At a site with several alleles most pairs
        // are carried by no haplotype, and which order fits such a pair better means nothing.
        if (before != nullptr && g == *before && frame_is_own) {
            reversed_called = prefers_other_order(strand0_allele, a, b, forward, reverse);
        }
        ll += max(forward, reverse) - over_groups(terms.unphased);
        ++corrected;
    }
    if (corrected == 0) {
        return false;
    }
    ++counters.sites_corrected;
    counters.sites_split_runs += grouping.agreement.size() > 1 ? 1 : 0;
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
