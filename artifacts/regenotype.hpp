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
 * the other sites the read spans in its phase set. Where those sites lie in more than one run of
 * read phasing, `y_r` comes from one run's sites (see `read_strand`), and a site's reads are scored
 * in groups by run, which the run correlations relate (see `phase_aware_correction`). At temper
 * 0, and for a read that spans no other site, `pi_r = v` and the read's correction is zero. A
 * homozygous genotype's read term does not depend on the weights, so its likelihood does not
 * change.
 *
 * The correction is the read term with `pi` minus the read term with `v`, and it is added to the
 * likelihood the direct pass computed. The rest of that likelihood, such as the depth term, stays
 * as it was, so only `PhaseReadEvidence` has to be kept. The direct pass's own read term used the
 * mixture weights, which equal `v` when the two alleles have the same length. The method is
 * described in artifacts/doc/read-likelihood-read-phasing.md, under "Re-genotyping from the phase".
 */

#include <array>
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
using std::pair;
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
    /// The part of `lambda` from the sites of one run (see `PhaseSite::run`).
    struct RunPart {
        /// The run, or `NO_RUN` for a part not yet used.
        size_t run = NO_RUN;
        /// How many of the read's sites lie in the run.
        size_t sites = 0;
        /// Their sum of the read's log-odds.
        double lambda = 0.0;
    };
    /// The parts of the first two runs the read's sites lie in. Sites with no run are not counted.
    ///
    /// Only the panel relates one run's strands to another's, so the read's strand is comparable
    /// only with those of the run its other sites lie in (see `read_strand`).
    std::array<RunPart, 2> run_parts;
    /// Whether the read's sites lie in more than two runs.
    bool many_runs = false;
};

/// Each read's `ReadLambda`, by read key.
using LambdaTable = unordered_map<uint64_t, ReadLambda>;

/// The phase set of a site that has none, because the linkage model did not phase it.
constexpr size_t NO_PHASE_SET = std::numeric_limits<size_t>::max();

/// A run of read phasing: its phase set, and its number there (see `PhaseSite::run`).
struct RunKey {
    size_t phase_set = NO_PHASE_SET;
    size_t run = NO_RUN;
    auto operator<=>(const RunKey&) const = default;
};

/// A read's strand log-odds at one site, and the run they come from.
struct ReadStrand {
    /// The log-odds, leaving the site's own contribution out.
    double log_odds = 0.0;
    /// The run they come from, with `NO_RUN` for none.
    RunKey run;
};

/// Return a read's strand log-odds at a site of run `site_run`, leaving out the site's own
/// contribution, `own_log_odds`, which `read` holds if `has_own`.
///
/// Only the panel relates one run's strands to another's, so the parts of a read's log-odds from
/// different runs do not add up. A read whose other sites lie in two runs keeps the part from
/// `site_run` if it is one of them, and otherwise the larger part. A read whose sites, this one
/// included, lie in more than two runs has no strand.
ReadStrand read_strand(const ReadLambda& read, size_t site_run, bool has_own, double own_log_odds);

/// Return log-odds `y` of a read's strand in its own run's frame, as seen from a run that the
/// panel orients against its own as reads would with probability `agreement`, at least 1/2.
///
/// They shrink towards 0 as `agreement` falls to 1/2, where they are exactly 0.
double reoriented_log_odds(double y, double agreement);

/**
 * A table of each run's correlation with the run before it, which gives how likely read phasing
 * orients one run against another as reads would.
 *
 * Reads fix the orientation of a run's sites against each other. Read phasing orients each run
 * against the one before it rightly with a probability p, and gives the run's correlation as
 * 2p - 1 (see `PhaseSite::run_correlation`). The table treats the errors at different breaks as
 * independent, so between two runs the correlations multiply. They are not independent where one
 * site the panel orders is the nearest to two breaks: its flip then counts twice, which understates
 * the agreement.
 */
class RunCorrelationTable {
public:
    /// Orient no run against another.
    RunCorrelationTable() = default;

    /// Take the correlation of each run of `sites`.
    explicit RunCorrelationTable(const vector<PhaseSite>& sites);

    /// Return the probability that read phasing orients run `b` against run `a` as reads would:
    /// (1 + the product of the correlations of the runs after the earlier one, up to the later)
    /// / 2. It is 1/2 for runs of different phase sets, and where either is in no run.
    double agreement(const RunKey& a, const RunKey& b) const;

private:
    /// The correlation of each run with the run before it.
    map<RunKey, double> correlation;
};

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
    /// Corrected sites whose reads were scored in more than one group, because some of their
    /// strands come from a run other than the frame group's (see `phase_aware_correction`).
    size_t sites_split_runs = 0;
    /// Nested chains at ploidy 1 whose reads were scored in more than one group, because some of
    /// their strands come from a run other than the frame group's (see
    /// `haploid_inclusion_correction`).
    size_t haploid_split_runs = 0;
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

/**
 * What one round of re-genotyping knows about the strand each read came from, which every site it
 * corrects shares.
 */
struct StrandModel {
    /// Each read's summed strand log-odds.
    const LambdaTable& lambda;
    /// How far each read's log-odds are believed.
    StrandCalibration calibration;
    /// How likely read phasing orients each run against another as reads would.
    const RunCorrelationTable& correlations;
    const RegenotypeParams& params;
};

/// Add the correction at a nested chain at ploidy 1 to its genotype likelihoods, in place.
///
/// The chain sits on one of its parent's two strands, and reads from the other strand can reach
/// it only as mismapped reads. `strand_sign` is +1 when the chain is on strand 0 of its phase set
/// and -1 on strand 1, so that a read's log-odds are for or against this chain. A read whose
/// strand is not usable in the chain's phase set (see `read_strand_usable`) is left unchanged:
///
///     incl  = min(1, exp(x))          x = strand_sign * calibration.log_odds(...)
///     e_r'  = e_r + (1 - e_r)(1 - incl)
///     term  = (1 - e_r) * incl * rel(r, a) + e_r'
///
/// We cap `incl` at 1 so that the term is unchanged at tau = 0, where the read's probability of
/// this strand would be 1/2. A read that points to the other strand counts less for or against the
/// allele.
///
/// `chain` is the chain's phase set and the run its strand is named in, the run of the parent's
/// order, or `NO_RUN` if that is not known. The reads are grouped by the run their strands come
/// from. The frame group, which keeps `strand_sign`, is that of the chain's run, or, where that is
/// not known, of the run most of the reads' strands come from, and holds the reads with no strand.
/// Each other group takes `strand_sign` with the probability `strands.correlations` gives for its
/// run against the frame group's, and the opposite sign otherwise.
///
/// Returns true if the corrected best genotype differs from the called one.
bool haploid_inclusion_correction(const PhaseReadEvidence& evidence, const StrandModel& strands,
                                  const RunKey& chain, int strand_sign,
                                  const unordered_map<uint64_t, double>& own,
                                  map<vector<int>, double>& gl, RegenotypeCounters& counters);

/// Add the phase-aware correction to one site's genotype likelihoods, in place.
///
/// `gl` is keyed by the sorted allele multiset, as `ReadLikelihoodCallInfo::genotype_lls` is.
/// Both assignments of a heterozygote's alleles to the strands are scored, and the larger
/// correction is kept, so the phase information is not averaged away.
///
/// Reads fix the strands of one run's sites relative to each other but not to another run's (see
/// `PhaseSite::run`), so the reads are scored in groups by the run their strands come from (see
/// `read_strand`). Reads with no strand
/// join the frame group: the group of the site's own run if a read with a strand comes from it,
/// and otherwise the largest group.
///
/// Each assignment of the frame group is scored with the other groups: each takes the same
/// assignment with the probability `strands.correlations` gives for its run against the frame
/// group's, and the other otherwise. The larger of the two totals is kept. A site whose reads with
/// a strand come from one run is scored as one group, as if there were no runs.
///
/// `site` is the site's phase set, or `NO_PHASE_SET`, and its own run. A read whose strand is not
/// usable in the phase set (see `read_strand_usable`) gets the site's own weights.
///
/// `own` is this site's own per-read contribution to Lambda, negated where read phasing flipped
/// the site. It is subtracted so that a site does not confirm its own genotype. Empty for a site
/// with no `PhaseSite`, such as a homozygote.
///
/// `strand0_allele` is the allele the chain puts on strand 0 at this site, or -1 where the chain
/// gives it no order. `order_reversed` counts the site when the reads prefer the other order, if
/// the site has no run or the frame group is its own run.
///
/// Returns true if the corrected best genotype differs from the called one.
bool phase_aware_correction(const PhaseReadEvidence& evidence, const StrandModel& strands,
                            const RunKey& site, int strand0_allele,
                            const unordered_map<uint64_t, double>& own,
                            map<vector<int>, double>& gl, RegenotypeCounters& counters);

/// Fill `out` with a site's contribution to each read's Lambda, for subtracting from each read's
/// total. Negate the contribution if `flipped`, which says that read phasing flipped the site.
void site_own_log_odds(const PhaseSite& site, bool flipped, unordered_map<uint64_t, double>& out);

}

#endif
