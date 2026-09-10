#ifndef VG_REGENOTYPE_HPP_INCLUDED
#define VG_REGENOTYPE_HPP_INCLUDED

/** \file regenotype.hpp
 *
 * Spend the reads' phase on the genotype rather than only on the order of an already-settled pair.
 *
 * `AlleleReadLikelihoods::genotype_likelihood` marginalises each read over the haplotypes of a
 * candidate genotype with a weight that belongs to the SITE:
 *
 *     mixture = sum_k w_k * rel(r, g_k);   total += log((1 - e_r) * mixture + e_r)
 *
 * `w_k` cannot depend on the read, because nothing at the site knows which haplotype the read came
 * from. The phasing chain does know, from the other sites the read touches -- a median of 40 of
 * them on ONT. So the weight becomes per-read:
 *
 *     pi_r^0 = w_0 e^{tau Lambda_r} / (w_0 e^{tau Lambda_r} + w_1),   pi_r^1 = 1 - pi_r^0
 *
 * where `Lambda_r` is read r's log-odds of lying on strand 0, summed over the OTHER sites it
 * touches. Three properties make this safe, and all three are asserted in the unit tests:
 *
 *   - At `tau = 0`, `pi_r == w` identically, so the correction is exactly zero and the caller's
 *     output is unchanged bit for bit. That is a gate, not an aspiration.
 *   - A read that touches no other site in its block has `Lambda_r == 0` after the leave-one-out,
 *     so `pi_r == w` for that read and its term cancels. Sites whose reads all span nothing else
 *     are genotyped by exactly the function that genotypes them today -- which is also why this
 *     degenerates harmlessly on short reads.
 *   - A homozygous genotype's mixture collapses to `rel(r, g_0)` whatever the weights are, so
 *     homozygous likelihoods never move. The correction can only re-rank hets against each other
 *     and hets against homs.
 *
 * The correction is a DELTA against the likelihood the sweep already computed, never a
 * re-derivation. `AlleleReadLikelihoods::mixture_weights` weights an allele by the sequence no
 * other member of the genotype carries, while `site_slot_weights` -- which is what is retained --
 * weights by full spelled length; the two disagree at every indel. Applying the same weight
 * function to both sides of the delta leaves only the `pi/w` tilt, and leaves the sweep's own
 * weighting untouched inside `ll(g)`. It also cancels the depth term, the `-ln n!` normaliser and
 * the un-retained row divisor, so nothing has to be retained beyond `PhaseReadEvidence`.
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
 * One read's accumulated evidence about which strand of its phase block it came from.
 *
 * Keyed by read alone rather than by (read, block): a read's sites are contiguous, so it lies in
 * one block, and the hom case needs a lookup that has no block to offer. A read that does turn up
 * in two blocks has no comparable strand between them and is marked unusable rather than averaged.
 */
struct ReadLambda {
    /// log-odds of strand 0 against strand 1, natural log, signed by each site's settled
    /// orientation. `tau` absorbs the base, so this differs from `phase_link`'s log10 by a
    /// constant that the temper fit swallows.
    double lambda = 0.0;
    /// How many sites contributed. A read at one site has nothing left after leave-one-out.
    size_t sites = 0;
    size_t phase_set = 0;
    /// Seen in more than one phase block, so `lambda` mixes two incomparable strand labellings.
    /// Such a read is given the site's own weights, exactly as if it spanned nothing.
    bool multi_block = false;
};

using LambdaTable = unordered_map<uint64_t, ReadLambda>;

struct RegenotypeParams {
    /// The temper. 0 reproduces the uncorrected caller exactly; negative means "fit it".
    ///
    /// `Lambda` sums a per-site log-odds over ~40 sites each worth up to a log10 unit, while the
    /// measured read-only switch error is 10^-2.4. Reads are not independent, because mapping
    /// error is correlated, so an untempered `pi` asserts a certainty the data does not support.
    double temper = -1.0;
    /// Bins for the calibration fit, over |Lambda|.
    size_t fit_bins = 12;
    /// A bin needs this many reads before it can move the fit.
    size_t fit_min_per_bin = 200;
    /// DEBUG. Randomise the sign of each read's Lambda, keeping |Lambda| -- so the peakedness
    /// distribution is preserved exactly and only the phase information is destroyed. A
    /// permutation across reads would not do that: it preserves the bimodal strand structure and
    /// can re-correlate by accident.
    ///
    /// The temper is fitted on the UNSHUFFLED evidence and then applied to the shuffled signs,
    /// and that is the point rather than an oversight. The control asks what a tilt of this
    /// strength buys with no phase in it, so the strength has to be held fixed; fitting on
    /// shuffled data would find no signal, return a temper near zero, and make the arm inert for
    /// a reason that has nothing to do with what it is controlling for.
    bool shuffle = false;
};

struct RegenotypeCounters {
    size_t reads_with_lambda = 0;
    size_t reads_multi_block = 0;
    /// Reads whose leave-one-out `Lambda` is zero because they touch no other site in the block.
    /// The correction is provably inert for these, and on short reads it is most of them.
    size_t reads_singleton = 0;
    size_t sites_considered = 0;
    size_t sites_corrected = 0;
    /// Sites whose corrected argmax names a different unordered genotype than the sweep's.
    size_t sites_would_move = 0;
    /// Of those, how many are homozygous-to-heterozygous and the reverse.
    size_t moved_hom_to_het = 0;
    size_t moved_het_to_hom = 0;
    size_t moved_het_to_het = 0;
    /// Sites where the reads prefer the reversed order AT THE CALLED GENOTYPE -- they disagree
    /// with the chain about this site's phase. A diagnostic on the chain, not an output, and
    /// counted only at the called pair: at a multi-allelic site most candidate pairs are ones no
    /// haplotype carries, and which order fits such a pair better is arbitrary.
    size_t order_reversed = 0;
    double fitted_temper = 0.0;
    /// Calibration table: per |Lambda| bin, predicted and observed agreement with the chain.
    vector<double> fit_abs_lambda, fit_predicted, fit_observed;
    vector<size_t> fit_count;
};

/// Fold one thread's counters into another's. Sums only: the calibration table is filled once,
/// before the parallel region, and must not be duplicated by the merge.
void merge_counters(const RegenotypeCounters& from, RegenotypeCounters& into);

/// One site's contribution to a read's strand log-odds, from the two numbers `PhaseSite` keeps.
///
/// `p` is the probability the read came from one of the two settled haplotypes at all, so the
/// mismapping escape survives into `Lambda` and a mismapped read cannot reach a confident strand.
/// Same `p * x + (1 - p) / 2` shape `phase_link` uses, and for the same reason.
double site_read_log_odds(double q0, double p);

/// Accumulate `Lambda` over every site, into `out`.
///
/// `flipped` is what `read_phase_flips` returned: a site in it has had its slot order swapped, so
/// its contribution enters with the opposite sign. Without that the accumulation would describe
/// the panel's frame rather than the settled one.
///
/// A read appearing twice at one site contributes once. Paired mates share a name and therefore a
/// `read_key`, and the read source deduplicates on name plus first mapping position, so both mates
/// are separate rows under one key -- 108,536 such names on chr20. `phase_link`'s sorted merge
/// pairs them 1:1 and is unharmed by this; a running sum is not.
void accumulate_lambda(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                       LambdaTable& out, RegenotypeCounters& counters);

/// Fit the temper from the run's own data, with no truth.
///
/// Bin reads by |Lambda| and measure, in each bin, how often a read's implied strand agrees with
/// the chain the other reads settled on; choose `tau` so the predicted agreement matches. The
/// chain was built from these same reads, so this is not fully out of sample -- it aggregates over
/// a whole block and one read's contribution to it is small, and the leave-one-out removes the
/// site-level part, but it will run optimistic. `--regeno-temper` overrides it so the sensitivity
/// to that optimism can be measured.
double fit_temper(const vector<PhaseSite>& sites, const unordered_set<size_t>& flipped,
                  const LambdaTable& lambda, const RegenotypeParams& params,
                  RegenotypeCounters& counters);

/// Add the phase-aware correction to one site's genotype likelihoods, in place.
///
/// `gl` is keyed by the sorted allele multiset, exactly as `ReadLikelihoodCallInfo::genotype_lls`
/// is, and stays that way: the two orders are scored and the better one taken, so the vector keeps
/// its layout and no consumer learns that an ordered quantity was computed. Taking the max rather
/// than the log-sum is deliberate -- marginalising the phase away would discard the information
/// being added -- and it also keeps `derive`'s runner-up from becoming a het's own mirror, which
/// would silently turn GQ into phase confidence.
///
/// `own` is this site's own per-read contribution to `Lambda`, signed as `accumulate_lambda`
/// signed it, for the leave-one-out subtraction. Empty for a site that contributed none: a
/// homozygote has no `PhaseSite`, and for it the leave-one-out is automatic.
///
/// Returns true if the corrected argmax names a different unordered genotype.
bool phase_aware_correction(const PhaseReadEvidence& evidence, const LambdaTable& lambda,
                            const unordered_map<uint64_t, double>& own, double temper,
                            const RegenotypeParams& params, map<vector<int>, double>& gl,
                            RegenotypeCounters& counters);

/// This site's per-read contribution to `Lambda`, for the leave-one-out. Signed the way
/// `accumulate_lambda` signed it, so the two cancel exactly.
void site_own_log_odds(const PhaseSite& site, bool flipped, unordered_map<uint64_t, double>& out);

}

#endif
