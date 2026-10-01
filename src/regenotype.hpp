#ifndef VG_REGENOTYPE_HPP_INCLUDED
#define VG_REGENOTYPE_HPP_INCLUDED

/** \file regenotype.hpp
 *
 * Re-genotyping from the phase (--regenotype): once read phasing has ordered the heterozygous
 * sites, a read that spans other heterozygous sites shows which strand it came from, and that is
 * used to correct each site's genotype likelihoods.
 *
 * One site cannot tell which strand a read came from, so the site likelihood weights every read's
 * haplotypes by the same mixture weights. Here each read at a heterozygous genotype gets its own
 * weights instead, from the allele-length weights `v` (`allele_length_weights`):
 *
 *     pi_r^0 = v_0 e^{y_r} / (v_0 e^{y_r} + v_1),   pi_r^1 = 1 - pi_r^0
 *
 * where `y_r` is the read's tempered strand log-odds (see `calibrated_log_odds`), computed from
 * the other sites the read spans in its phase set. At temper 0, and for a read that spans no other
 * site, `pi_r = v` and the read's correction is zero. A homozygous genotype's read term does not
 * depend on the weights, so its likelihood does not change.
 *
 * The correction is the read term with `pi` minus the read term with `v`, and it is added to the
 * likelihood the sweep computed. Everything else in the sweep's likelihood, such as the depth
 * term, is left as it was, so only `PhaseReadEvidence` has to be kept. The sweep's own read term
 * used the mixture weights, which equal `v` when the two alleles have the same length. The method
 * is described in doc/read-likelihood-genotyping.md, under "Re-genotyping from the phase".
 */

#include <cstddef>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "read_phasing.hpp"

namespace vg {

using std::map;
using std::size_t;
using std::uint64_t;
using std::unordered_map;
using std::unordered_set;
using std::vector;

/**
 * One read's accumulated evidence about which strand of its phase set it came from.
 *
 * Keyed by read alone rather than by (read, phase set): a read normally lies in one phase set. A
 * read found in two has no comparable strand between them, and is marked unusable.
 */
struct ReadLambda {
    /// Natural-log odds of strand 0 against strand 1, with each site's contribution signed by
    /// its settled order.
    double lambda = 0.0;
    /// How many sites contributed.
    size_t sites = 0;
    size_t phase_set = 0;
    /// Seen in more than one phase set, so `lambda` mixes two unrelated strand labellings. Such a
    /// read is given the site's own weights, as if it spanned no other site.
    bool multi_phase_set = false;
};

using LambdaTable = unordered_map<uint64_t, ReadLambda>;

struct RegenotypeParams {
    /// The temper tau (--regeno-temper); 0 leaves the likelihoods unchanged, and a negative value
    /// means fit it (see `fit_calibration`). Reads are not independent, so the summed strand
    /// log-odds overstate how sure the strand is, and tau scales them down.
    double temper = -1.0;
    /// The ceiling c on a read's strand probability (--regeno-ceiling), in (0, 1]:
    /// `P(on strand 0) = c * sigmoid(tau * Lambda) + (1 - c) / 2`. A value below 1 keeps either
    /// strand's probability further from 1; 1 turns the ceiling off. At tau = 0 the probability
    /// is 1/2 whatever c is.
    double ceiling = 1.0;
    /// At a nested chain at ploidy 1, make reads less informative when the phase places them on
    /// the parent's other strand (on unless --no-regeno-haploid). The chain has one strand, so
    /// the mixture correction does nothing there. See `haploid_inclusion_correction`.
    bool haploid_include = true;
    /// Bins for the calibration fit, over |Lambda|.
    size_t fit_bins = 12;
    /// The fewest observations the fit needs, and the fewest per bin: each bin holds the larger of
    /// this and a `fit_bins`-th of the observations, except the last, which holds the remainder.
    size_t fit_min_per_bin = 200;
    /// For testing (--regeno-shuffle): randomise the sign of each read's Lambda, keeping |Lambda|,
    /// which removes the phase information and keeps the rest. The temper is fitted before the
    /// signs are randomised, so the control tests a correction of the same strength with no phase
    /// information in it.
    bool shuffle = false;
};

struct RegenotypeCounters {
    size_t reads_with_lambda = 0;
    size_t reads_multi_phase_set = 0;
    /// Reads whose Lambda, leaving out the site itself, is zero because they span no other site
    /// in the phase set. The correction does nothing for these.
    size_t reads_singleton = 0;
    size_t sites_considered = 0;
    size_t sites_corrected = 0;
    /// Sites whose corrected argmax names a different unordered genotype than the sweep's.
    size_t sites_would_move = 0;
    /// Of those, how many are homozygous-to-heterozygous and the reverse.
    size_t moved_hom_to_het = 0;
    size_t moved_het_to_hom = 0;
    size_t moved_het_to_het = 0;
    /// Sites where the reads prefer the reversed order at the called genotype, so they disagree
    /// with the chain about this site's phase. Counted only at the called pair, since at a site
    /// with several alleles most candidate pairs are carried by no haplotype.
    size_t order_reversed = 0;
    double fitted_temper = 0.0;
    double fitted_ceiling = 1.0;
    /// Nested chains at ploidy 1 that the inclusion weight reached, and how many it moved.
    size_t haploid_sites = 0;
    size_t haploid_would_move = 0;
    /// Calibration table: per |Lambda| bin, predicted and observed agreement with the chain.
    vector<double> fit_abs_lambda, fit_predicted, fit_observed;
    vector<size_t> fit_count;
};

/// Add one thread's counters into another's. The calibration table is filled once, before the
/// parallel region, so it is not merged.
void merge_counters(const RegenotypeCounters& from, RegenotypeCounters& into);

/// One site's contribution to a read's strand log-odds, from its `q0` and `c` (see `PhaseSite`).
/// It uses the same `c * x + (1 - c) / 2` form as `phase_link`, so a probably mismapped read
/// contributes little.
double site_read_log_odds(double q0, double c);

/// Accumulate each read's Lambda over every site, into `out`.
///
/// `flipped` is what `read_phase_flips` returned: a site in it has had its pair swapped, so its
/// contribution enters with the opposite sign. Each site contributes once per read key, so paired
/// mates count once (see `merge_mates`).
void accumulate_lambda(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                       LambdaTable& out, RegenotypeCounters& counters);

/// Fit the temper from the run's own data, against the given ceiling.
///
/// Reads are grouped into bins by |Lambda|, and in each bin we measure how often the strand that
/// Lambda points to agrees with the one the read's own allele points to. The temper is chosen so
/// that the predicted agreement matches. The chain was built from these same reads, so the fit is
/// somewhat optimistic; `--regeno-temper` sets the temper directly instead.
void fit_calibration(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                     const LambdaTable& lambda, const RegenotypeParams& params,
                     double& temper, double& ceiling, RegenotypeCounters& counters);

/// A read's tempered strand log-odds, `logit(c * sigmoid(tau * Lambda) + (1 - c) / 2)`. It is 0
/// at tau = 0 for every ceiling. The probability is kept away from 0 and 1, where the logit is
/// infinite.
double calibrated_log_odds(double lambda, double temper, double ceiling);

/// The correction at a nested chain at ploidy 1, where the mixture correction is zero.
///
/// The chain sits on one of its parent's two strands, and reads from the other strand can reach
/// it only as mismapped reads. `strand_sign` is +1 when the chain is on strand 0 of its phase set
/// and -1 on strand 1, so that a read's log-odds are for or against this chain:
///
///     incl  = min(1, exp(x))          x = strand_sign * calibrated_log_odds(...)
///     e_r'  = e_r + (1 - e_r)(1 - incl)
///     term  = (1 - e_r) * incl * rel(r, a) + e_r'
///
/// `incl` is capped at 1 rather than being the read's probability of this strand, which would
/// be 1/2 at tau = 0 and change the result; capped, the term is unchanged at tau = 0. A read that
/// points to the other strand counts less for or against the allele.
///
/// Returns true if the corrected best genotype differs from the called one.
bool haploid_inclusion_correction(const PhaseReadEvidence& evidence, const LambdaTable& lambda,
                                  const unordered_map<uint64_t, double>& own, double temper,
                                  double ceiling, int strand_sign,
                                  const RegenotypeParams& params, map<vector<int>, double>& gl,
                                  RegenotypeCounters& counters);

/// Add the phase-aware correction to one site's genotype likelihoods, in place.
///
/// `gl` is keyed by the sorted allele multiset, as `ReadLikelihoodCallInfo::genotype_lls` is.
/// Both assignments of a heterozygote's alleles to the strands are scored and the larger
/// correction kept, rather than summed over, so the phase information is not averaged away.
///
/// `own` is this site's own per-read contribution to Lambda, signed as `accumulate_lambda`
/// signed it, which is subtracted so that a site does not confirm its own genotype. Empty for a
/// site with no `PhaseSite`, such as a homozygote.
///
/// Returns true if the corrected best genotype differs from the called one.
bool phase_aware_correction(const PhaseReadEvidence& evidence, const LambdaTable& lambda,
                            const unordered_map<uint64_t, double>& own, double temper,
                            double ceiling, const RegenotypeParams& params,
                            map<vector<int>, double>& gl, RegenotypeCounters& counters);

/// This site's per-read contribution to Lambda, signed as `accumulate_lambda` signed it, for
/// subtracting from each read's total.
void site_own_log_odds(const PhaseSite& site, bool flipped, unordered_map<uint64_t, double>& out);

}

#endif
