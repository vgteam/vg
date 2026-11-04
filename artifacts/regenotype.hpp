#ifndef VG_REGENOTYPE_HPP_INCLUDED
#define VG_REGENOTYPE_HPP_INCLUDED

/** \file regenotype.hpp
 *
 * Re-genotyping from the phase: correct each site's genotype likelihoods with the strand each read
 * came from.
 *
 * Once read phasing has ordered the heterozygous sites, a read that spans other heterozygous sites
 * shows which strand it came from. One site alone cannot tell, so the site likelihood weights every
 * read's haplotypes by the same mixture weights. Here each read at a heterozygous genotype gets its
 * own weights instead, from the allele-length weights `v` (`allele_length_weights`):
 *
 *     pi_r^0 = v_0 e^{y_r} / (v_0 e^{y_r} + v_1),   pi_r^1 = 1 - pi_r^0
 *
 * where `y_r` is the read's tempered strand log-odds (see `StrandCalibration`), computed from
 * the other sites the read spans in its phase set. At temper 0, and for a read that spans no other
 * site, `pi_r = v` and the read's correction is zero. A homozygous genotype's read term does not
 * depend on the weights, so its likelihood does not change.
 *
 * The correction is the read term with `pi` minus the read term with `v`, and it is added to the
 * likelihood the direct pass computed. The rest of that likelihood, such as the depth term, stays
 * as it was, so only `PhaseReadEvidence` has to be kept. The direct pass's own read term used the
 * mixture weights, which equal `v` when the two alleles have the same length. The method is
 * described in artifacts/doc/read-likelihood-read-phasing.md, under "Re-genotyping from the phase".
 */

#include <cstddef>
#include <cstdint>
#include <limits>
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
 * A read normally lies in one phase set, so this is kept per read. Each phase set labels its
 * strands independently, so a read's strand is usable only at a site of the phase set it was found
 * in (see `read_strand_usable`).
 */
struct ReadLambda {
    /// Lambda, the natural-log odds of strand 0 against strand 1.
    ///
    /// Each site's contribution is signed by its chosen order.
    double lambda = 0.0;
    /// How many sites contributed.
    size_t sites = 0;
    /// The phase set the read was found in, as an ID unique across contigs.
    size_t phase_set = 0;
    /// Whether the read was seen in more than one phase set.
    ///
    /// Then `lambda` mixes unrelated strand labellings, and the read's strand is usable nowhere.
    bool multi_phase_set = false;
};

/// Each read's `ReadLambda`, by read key.
using LambdaTable = unordered_map<uint64_t, ReadLambda>;

/// The phase set of a site that has none, because the linkage model did not phase it.
constexpr size_t NO_PHASE_SET = std::numeric_limits<size_t>::max();

/// Return whether a read's strand can be used at a site of phase set `phase_set`.
///
/// It can when the read was found in one phase set, and that is `phase_set`. At a site with no
/// phase set (`NO_PHASE_SET`), any read found in one phase set is usable. An unusable read is given
/// the site's own weights, as if it spanned no other site.
bool read_strand_usable(const ReadLambda& read, size_t phase_set);

/// The settings of re-genotyping from the phase.
struct RegenotypeParams {
    /// The temper tau, or a negative value to fit it (see `fit_calibration`).
    ///
    /// Reads are not independent, so the summed strand log-odds overstate how sure the strand is,
    /// and tau scales them down. At 0 the likelihoods do not change.
    double temper = -1.0;
    /// The ceiling c on a read's strand probability, in (0, 1].
    ///
    /// `P(on strand 0) = c * sigmoid(tau * Lambda) + (1 - c) / 2`. A value below 1 keeps either
    /// strand's probability further from 1, and 1 turns the ceiling off.
    double ceiling = 1.0;
    /// Whether to correct nested chains at ploidy 1 (see `haploid_inclusion_correction`).
    ///
    /// The chain has one strand, so the mixture correction does nothing there.
    bool haploid_include = true;
    /// The number of bins for the calibration fit, over |Lambda|.
    size_t fit_bins = 12;
    /// The fewest observations the fit needs, and the fewest per bin.
    ///
    /// Each bin holds the larger of this and a `fit_bins`-th of the observations, except the last,
    /// which holds the remainder.
    size_t fit_min_per_bin = 200;
    /// Whether to randomise the sign of each read's Lambda, for testing.
    ///
    /// This removes the phase information and keeps |Lambda|. The temper is fitted before the
    /// signs are randomised, so the control tests a correction of the same strength with no phase
    /// information in it.
    bool shuffle = false;
};

/**
 * A StrandCalibration is the fitted curve that turns a read's summed strand log-odds into the
 * log-odds it is believed at.
 *
 * The curve is `logit(c * sigmoid(tau * Lambda) + (1 - c) / 2)`, for temper tau and ceiling c.
 * Reads are not independent, so the summed strand log-odds overstate how sure the strand is. The
 * temper scales them down, and the ceiling keeps either strand's probability away from 1.
 */
struct StrandCalibration {
    /// The temper tau. At 0 every read's log-odds are 0, so the strand of no read is known.
    double temper = 0.0;
    /// The ceiling c on a read's strand probability, in (0, 1]. At 1 there is no ceiling.
    double ceiling = 1.0;

    /// Return the log-odds a read with summed strand log-odds `lambda` is believed at.
    ///
    /// This is exactly 0 at temper 0, for every ceiling. The probability is kept away from 0 and 1,
    /// where the logit is infinite.
    double log_odds(double lambda) const;
};

/// A report of what re-genotyping did.
struct RegenotypeCounters {
    /// Reads with a Lambda.
    size_t reads_with_lambda = 0;
    /// Reads seen in more than one phase set.
    size_t reads_multi_phase_set = 0;
    /// Reads seen at only one site, which the correction does nothing for.
    size_t reads_singleton = 0;
    /// Sites with reads and genotype likelihoods.
    size_t sites_considered = 0;
    /// Sites whose likelihoods were corrected.
    size_t sites_corrected = 0;
    /// Sites whose corrected argmax names a different unordered genotype than the direct pass's.
    size_t sites_would_move = 0;
    /// Of those, the moves from a homozygous to a heterozygous genotype.
    size_t moved_hom_to_het = 0;
    /// Of those, the moves from a heterozygous to a homozygous genotype.
    size_t moved_het_to_hom = 0;
    /// Of those, the moves from one heterozygous genotype to another.
    size_t moved_het_to_het = 0;
    /// Sites where the reads prefer the reversed order at the called genotype, so they disagree
    /// with the chain about this site's phase. Counted only at the called pair, since at a site
    /// with several alleles most candidate pairs are carried by no haplotype.
    size_t order_reversed = 0;
    /// The calibration the round used, which `fit_calibration` fills when it fits one.
    StrandCalibration calibration;
    /// Nested chains at ploidy 1 where some read's strand log-odds are nonzero.
    size_t haploid_sites = 0;
    /// Of those, the chains whose best genotype the correction changed.
    size_t haploid_would_move = 0;
    /// The mean |Lambda| of each bin of the calibration fit.
    vector<double> fit_abs_lambda;
    /// The agreement with the phase chain the fit predicts for each bin.
    vector<double> fit_predicted;
    /// The agreement with the phase chain observed in each bin.
    vector<double> fit_observed;
    /// The number of observations in each bin.
    vector<size_t> fit_count;
};

/// Add one thread's counters into another's.
///
/// The fit's results are filled once, before the threads start, so they are not added.
void merge_counters(const RegenotypeCounters& from, RegenotypeCounters& into);

/// Return one site's contribution to a read's strand log-odds, from its `q0` and `c` (see
/// `PhaseSite`).
///
/// A probably mismapped read contributes little.
double site_read_log_odds(double q0, double c);

/// Accumulate each read's Lambda over every site, into `out`.
///
/// `flipped` holds the record keys of the sites read phasing flipped, whose contributions enter
/// with the opposite sign. Each site contributes once per read key, so paired mates count once.
void accumulate_lambda(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                       LambdaTable& out, RegenotypeCounters& counters);

/// Fit the temper from the run's own data, against the ceiling in `params`.
///
/// Reads are grouped into bins by |Lambda|, and in each bin we measure how often the strand that
/// Lambda points to agrees with the one the read's own allele points to. The temper is chosen so
/// that the predicted agreement matches. The phase chain was built from these same reads, so the
/// fit is somewhat optimistic.
///
/// With too few observations, the temper is 0 and the ceiling 1.
///
/// The fitted calibration is returned and also stored in `counters.calibration`, next to the
/// calibration table of the fit.
StrandCalibration fit_calibration(const vector<PhaseSite>& sites,
                                  const unordered_set<size_t>& flipped, const LambdaTable& lambda,
                                  const RegenotypeParams& params, RegenotypeCounters& counters);

/// Add the correction at a nested chain at ploidy 1 to its genotype likelihoods, in place.
///
/// The chain sits on one of its parent's two strands, and reads from the other strand can reach
/// it only as mismapped reads. `strand_sign` is +1 when the chain is on strand 0 of its phase set
/// `phase_set` and -1 on strand 1, so that a read's log-odds are for or against this chain. A read
/// whose strand is not usable in `phase_set` (see `read_strand_usable`) is left unchanged:
///
///     incl  = min(1, exp(x))          x = strand_sign * calibration.log_odds(...)
///     e_r'  = e_r + (1 - e_r)(1 - incl)
///     term  = (1 - e_r) * incl * rel(r, a) + e_r'
///
/// We cap `incl` at 1 so that the term is unchanged at tau = 0, where the read's probability of
/// this strand would be 1/2. A read that points to the other strand counts less for or against the
/// allele.
///
/// Returns true if the corrected best genotype differs from the called one.
bool haploid_inclusion_correction(const PhaseReadEvidence& evidence, const LambdaTable& lambda,
                                  size_t phase_set,
                                  const unordered_map<uint64_t, double>& own,
                                  const StrandCalibration& calibration, int strand_sign,
                                  const RegenotypeParams& params, map<vector<int>, double>& gl,
                                  RegenotypeCounters& counters);

/// Add the phase-aware correction to one site's genotype likelihoods, in place.
///
/// `gl` is keyed by the sorted allele multiset, as `ReadLikelihoodCallInfo::genotype_lls` is.
/// Both assignments of a heterozygote's alleles to the strands are scored, and the larger
/// correction is kept, so the phase information is not averaged away.
///
/// `own` is this site's own per-read contribution to Lambda, negated where read phasing flipped
/// the site. It is subtracted so that a site does not confirm its own genotype. Empty for a site
/// with no `PhaseSite`, such as a homozygote.
///
/// `phase_set` is the site's phase set, or `NO_PHASE_SET`; a read whose strand is not usable there
/// (see `read_strand_usable`) gets the site's own weights.
///
/// `strand0_allele` is the allele the chain puts on strand 0 at this site, or -1 where the chain
/// gives it no order. `order_reversed` counts the site when the reads prefer the other order.
///
/// Returns true if the corrected best genotype differs from the called one.
bool phase_aware_correction(const PhaseReadEvidence& evidence, const LambdaTable& lambda,
                            size_t phase_set, int strand0_allele,
                            const unordered_map<uint64_t, double>& own,
                            const StrandCalibration& calibration, const RegenotypeParams& params,
                            map<vector<int>, double>& gl, RegenotypeCounters& counters);

/// Fill `out` with a site's contribution to each read's Lambda, for subtracting from each read's
/// total. Negate the contribution if `flipped`, which says that read phasing flipped the site.
void site_own_log_odds(const PhaseSite& site, bool flipped, unordered_map<uint64_t, double>& out);

}

#endif
