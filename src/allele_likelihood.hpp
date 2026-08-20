#ifndef VG_ALLELE_LIKELIHOOD_HPP_INCLUDED
#define VG_ALLELE_LIKELIHOOD_HPP_INCLUDED

/** \file allele_likelihood.hpp
 *
 * A reads x alleles matrix of P(read | allele) for one site, and the machinery
 * to build it from the reads' existing alignments in the graph.
 *
 * This is the evidence layer for read-level genotyping: where the support-based
 * callers ask "how many reads cover this allele", this asks "how well does each
 * individual read fit each allele", which is what a proper P(reads | genotype)
 * needs.
 *
 * The model this feeds -- the objective, every term, every parameter and what the
 * reported qualities mean -- is written up in doc/read-likelihood-genotyping.md.
 * Read that first; the comments here assume it and cover only what is local.
 */

#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <vg/vg.pb.h>

#include "alignment_scorer.hpp"
#include "handle.hpp"
#include "site_read_source.hpp"
#include "snarls.hpp"

namespace vg {

using namespace std;

/**
 * Per-read relative likelihoods over the alleles at one site.
 *
 * Rows are **normalised by their own maximum**, so every entry is in [0,1] and
 * each row's maximum is exactly 1. Each read's likelihoods are therefore
 * expressed relative to that read's best available explanation at this site.
 *
 * Normalising this way matters for more than tidiness. What the scoring produces
 * is ln P(read | allele) only up to a per-read constant (the score to likelihood
 * conversion is just `log_base * score`, with an unknown and uncomputable
 * normaliser per read). That constant cancels when comparing genotypes, but only
 * as long as nothing else enters the per-read expression. The moment a
 * mismapping term is added it stops cancelling, unless every component sits on
 * the same scale. Dividing each row by its own maximum enforces that by
 * construction: the background becomes dimensionless, exactly 1, and there is no
 * scale left to mix.
 *
 * It also removes every numerical hazard. "No valid placement" becomes 0 rather
 * than -inf, underflow gives the correct answer instead of needing care, and the
 * per-read term is provably finite, so no logsumexp is required anywhere. As a
 * side benefit it makes rows scored under different scorers safely comparable,
 * which is what lets a read without base qualities sit in the same matrix as one
 * with them.
 *
 * Do NOT normalise the other way. Making the alleles sum to 1 within a row is a
 * plausible-sounding mistake: that quantity is a posterior over alleles under a
 * uniform prior, not P(read | allele). It double counts the prior and destroys
 * the absolute-fit information the mismapping term depends on, making a read
 * that fits one allele badly and the rest worse indistinguishable from a read
 * that fits one allele perfectly.
 *
 * These values are valid for argmax, for GQ (a likelihood difference), and for a
 * posterior over genotypes. They are NOT calibrated absolute probabilities.
 */
class AlleleReadLikelihoods {
public:

    AlleleReadLikelihoods() = default;

    size_t num_reads() const { return n_reads; }
    size_t num_alleles() const { return n_alleles; }

    /// Likelihood of read r under allele a, relative to read r's best allele at
    /// this site. In [0,1]. Exactly 0 means this allele cannot place the read at
    /// all, which is strong evidence *against* the allele, not missing data.
    double rel(size_t r, size_t a) const;

    /// The read's mismapping probability e_r, derived from MAPQ and already
    /// clamped. Guaranteed strictly inside (0,1) so the per-read term is finite.
    double mismap_prob(size_t r) const;

    /// ln of the row's divisor: the read's absolute best fit at this site.
    /// The model does not use this. It is the only surviving record of absolute
    /// fit, kept for diagnostics and as a realignment trigger.
    double best_ln_likelihood(size_t r) const;

    /// Name of read r, for the likelihood dump. May be empty.
    const string& read_name(size_t r) const;

    /// How many reads were dropped for placing on no allele at all. A rising
    /// count means the read source is over-fetching, or the reads and the graph
    /// do not match.
    size_t num_unplaceable() const { return unplaceable; }

    /// Weight the mixture by how many reads each haplotype is *expected* to
    /// contribute at this site, instead of a flat 1/|G|.
    ///
    /// The flat weight asserts that each haplotype of a genotype produced half
    /// the reads. Over an interval where one haplotype carries a deletion that is
    /// simply false -- the deleted haplotype produces no reads there at all -- and
    /// the error grows with the length imbalance. It is why a read lying inside a
    /// heterozygous deletion argues for the homozygous long genotype by ln 2, and
    /// why large heterozygous deletions are lost (and large heterozygous
    /// insertions mis-genotyped, the same failure mirrored).
    ///
    /// The number of read start positions that yield a read overlapping this site
    /// from a haplotype presenting an allele of length L is L + R - 1, for read
    /// length R. So
    ///
    ///     w_h = (L_h + R - 1) / sum_{h' in G} (L_h' + R - 1)
    ///
    /// Three properties this has and a plain maximum does not. The weights sum to
    /// 1, so adding an allele still costs and a clean homozygote still beats the
    /// heterozygote. Equal-length alleles give exactly 1/2, so SNVs and balanced
    /// indels are unchanged bit for bit. And it is symmetric in the direction of
    /// the imbalance, so it addresses insertions and deletions with one rule.
    ///
    /// `allele_lengths` is indexed by allele, `mean_read_length` is R. Passing an
    /// empty vector, or a zero R, falls back to the flat 1/|G|.
    void set_length_weights(vector<size_t> allele_lengths, double mean_read_length) {
        this->allele_lengths = std::move(allele_lengths);
        this->mean_read_length = mean_read_length;
    }

    /// Sharpen the weights by counting only sequence *unique* to each allele.
    ///
    /// Whole traversal length over-counts the shorter allele, because a traversal
    /// includes the site's shared sequence -- at one measured 2648 bp deletion the
    /// traversals are 296 and 2945 bp, so the raw ratio is 6.9 where the reads
    /// actually split about 14.6. Reads landing in shared sequence fit every allele
    /// equally, contribute the same factor to every genotype, and cancel; only
    /// reads overlapping sequence unique to one allele can move a genotype
    /// comparison at all. So the quantity the weight wants is unique content:
    ///
    ///     w_h = (U_h + R - 1) / sum_{h' in G} (U_h' + R - 1)
    ///
    /// where U_h is the sequence h visits and the other allele of the genotype does
    /// not. Balanced alleles still give exactly 1/2 -- a SNV has U = 1 on both
    /// sides -- so the no-op-on-small-variants property is unaffected.
    ///
    /// `unique_lengths[a][b]` is the sequence in allele a but not allele b. For a
    /// diploid genotype that is exact. Above diploid the minimum over the other
    /// members is used, which over-states uniqueness when three or more alleles
    /// overlap partially; `vg call` genotypes diploid, so this is a bound rather
    /// than a live approximation.
    void set_unique_lengths(vector<vector<size_t>> unique_lengths) {
        this->unique_lengths = std::move(unique_lengths);
    }

    /// Add ln P(N | G) to the genotype likelihood: does this genotype predict the
    /// number of reads actually seen?
    ///
    /// The model otherwise scores P(reads | G) conditioned on the reads it was
    /// handed, and never asks whether that many reads should be there. A complete
    /// generative model factorises as P(N | G) * P(reads | N, G), so the missing
    /// piece is additive:
    ///
    ///     ln P(data | G) = w * ln Poisson(N ; lambda_G) + sum_r ln[...]
    ///     lambda_G       = rate * sum_{h in G} (T_h + R - 1)
    ///
    /// `T_h` is the allele's **interior** traversal length -- the sequence over which
    /// a read can become a row of this matrix. It is neither the whole traversal, whose
    /// two boundary nodes recruit no rows because a read inside one of them cannot
    /// discriminate, nor the unique content the mixture weights use. All three are
    /// different quantities and conflating any two is an easy, invisible error: a
    /// mixture weight asks which allele a read could *distinguish*, while lambda asks
    /// how much sequence can put a read in front of the question at all.
    ///
    /// `rate` is reads per position per haplotype, estimated locally. It must be
    /// measured through the same fetch and placement path that produced `N`, since
    /// `N` is rows in this matrix -- neither coverage nor what a `vg pack` index
    /// reports -- and a rate in different units carries a silent scale error.
    ///
    /// `effective_count` decides what `N` means. The read term already believes each
    /// read only to the extent of `1 - e_r`; counting that same read as a whole read
    /// of depth asserts something the read term explicitly declines to. Under
    /// `effective_count` the observation is `N_eff = sum_r (1 - e_r)`, the expected
    /// number of reads genuinely from this locus, and `rate` must be measured the
    /// same way or the two carry a constant scale factor between them.
    ///
    /// Off unless `weight` is positive.
    /// Rescale the per-haplotype depth rate.
    ///
    /// The rate is the local read rate divided by the ploidy the site is being genotyped at, so
    /// scoring the same matrix at a different ploidy needs it scaled by the ratio -- 0.5 to go from
    /// haploid to diploid. Exact, not an approximation: nothing else in the matrix depends on
    /// ploidy, since `rel(r,a)` is a per-read-per-allele fit and lambda_G already sums over the
    /// genotype's haplotypes. Scoring the other ploidy without this leaves lambda wrong by the
    /// ploidy ratio, in the one term where that error is invisible in the output.
    void scale_depth_rate(double factor) {
        this->depth_rate *= factor;
    }

    void set_depth_context(vector<size_t> traversal_lengths, double rate,
                           double read_length, double weight,
                           bool effective_count = true) {
        this->traversal_lengths = std::move(traversal_lengths);
        this->depth_rate = rate;
        this->depth_read_length = read_length;
        this->depth_weight = weight;
        this->depth_effective = effective_count;
    }

    bool uses_depth_term() const {
        return depth_weight > 0.0 && depth_rate > 0.0 && !traversal_lengths.empty();
    }

    /// Expected reads at this site under G, from the same geometry the term uses.
    double expected_reads(const vector<int>& genotype) const;

    /// The read count the depth term compares against: `sum_r (1 - e_r)` when the
    /// context was set with `effective_count`, and the plain row count otherwise.
    /// Fractional by construction, which is why the Poisson uses lgamma.
    double observed_reads() const;

    /// This allele's interior traversal length, as lambda uses it. Zero if the depth
    /// context was never set or the index is out of range, so a difference between two
    /// of these is only meaningful when uses_depth_term() or a rate was supplied.
    size_t traversal_length(size_t allele) const {
        return allele < traversal_lengths.size() ? traversal_lengths[allele] : 0;
    }

    /// Observed over expected, for the genotype given. 1.0 is a site whose read
    /// count is exactly what the call predicts; 7.0 is a collapsed repeat. Emitted
    /// as `DR` whether or not the term is switched on, so it can be measured as a
    /// ranking signal before it is trusted as a likelihood.
    double depth_ratio(const vector<int>& genotype) const;

    /// Mean read length in this site's own matrix, for the depth term's geometry.
    double mean_read_length_estimate() const { return mean_read_length; }

    /// Supply the mean read length on its own, for the depth term's lambda = rate * (L + R - 1).
    /// set_length_weights also sets it, but only runs under the length-weighted mixture -- and the
    /// depth term stays armed under --flat-mixture, where leaving R at 0 made lambda ~150x too
    /// small for an SNV at 150 bp reads and erased the missing-reads signal for indels.
    void set_mean_read_length(double mean_read_length) {
        this->mean_read_length = mean_read_length;
    }

    /// True when set_length_weights supplied usable data.
    bool uses_length_weights() const {
        return !allele_lengths.empty() && mean_read_length > 0.0;
    }

    /**
     * ln P(reads | G), where G is a multiset of allele indices of size ploidy.
     *
     * **This is the objective the caller maximises.** Both terms are on by default,
     * so a description of only the first describes a configuration nobody runs.
     *
     *   ln P(reads | G) =  sum_r ln [ (1 - e_r) * sum_{h in G} w_h * rel(r,h) + e_r ]
     *                    + w_d * ln Poisson( N_eff ; lambda_G )
     *
     * where r runs over reads overlapping the site, h over the haplotypes of G
     * (twice over the same allele for a homozygote), rel(r,h) is the read's fit to
     * allele h relative to its own best allele here, e_r its mismapping probability,
     * w_h the mixture weight, and the second term a Poisson on read count.
     *
     * **Every symbol, every parameter and the reasoning behind each is specified in
     * doc/read-likelihood-genotyping.md**, which is the canonical description of the
     * model. Kept there rather than here because the model spans this header, the
     * linkage layer and the VCF writer, and no one of them can hold all of it. What
     * follows is only what a reader *modifying this function* needs.
     *
     * The bracket lies in [e_r, 1], so its log is finite and bounded below by
     * ln(e_r): no single read can penalise a genotype without limit. That bound is
     * the whole reason for the floor, and it is why min_mismap_prob reads as
     * "P(this read's evidence here is unreliable)" rather than as mismapping alone.
     * Anything that removes the bound -- an e_r of zero, a mixture allowed outside
     * [0,1] -- reintroduces the unbounded single-read veto that the floor exists to
     * stop, and it will present as spurious heterozygosity rather than as an error.
     *
     * rel(r,h) = 0 means h cannot place the read at all. That is evidence *against*
     * h, not missing data, and must not be skipped or imputed.
     *
     * The read term has no opinion about reads that are *absent*, which is exactly
     * the evidence a homozygous deletion presents; the depth term supplies it, at a
     * deliberately small weight because it is a much cruder statistic and is
     * dominated by the read term wherever reads are informative.
     *
     * **Not in this expression:** the linkage layer. `--linkage-weight` re-decides
     * genotypes *after* calling, from forward-backward posteriors over pairs of
     * GBWT panel haplotypes, using this quantity as the per-site emission. See
     * linkage_model.hpp. Nothing here changes when it is on.
     *
     * Reads are assumed independent. Mates overlapping the same site are not, and
     * the product therefore accumulates confidence like R rather than sqrt(R), so
     * derived GQ will be over-confident and increasingly so with depth. Treat
     * GQ/GL as useful for ranking, not as calibrated probabilities.
     */
    double genotype_likelihood(const vector<int>& genotype) const;

    /**
     * The largest likelihood gap this site could produce between these two genotypes:
     * what `genotype_likelihood(called) - genotype_likelihood(runner_up)` would come
     * to if every read fitted `called` perfectly and nothing else.
     *
     * This exists because the raw gap -- which is what GQ reports -- is not comparable
     * between sites, and it moves on two axes: depth and **ploidy**. Depth is the
     * obvious one and the easy one. Ploidy is neither: at ploidy 1 the runner-up is a
     * different allele outright, so every read discriminates fully, while at ploidy 2
     * a heterozygote's runner-up differs on one strand only, so a read discriminates
     * about half as much and only half the reads discriminate at all. Per read, with
     * the mismap floor at its 0.02 default, that is 3.91 nats haploid against 1.28 for
     * a diploid het and 0.67 for a diploid hom -- so a hemizygous call carries some six
     * times the GQ of a diploid homozygote on identical evidence. Measured on HG002,
     * chrX hemizygous calls run a median GQ of 247 where chr7 diploid homs at the same
     * depth run 46.
     *
     * Dividing the observed gap by this one gives a fraction in [0,1] meaning the same
     * thing at any depth and any ploidy. Measured over a coverage titration of chr20
     * (diploid, 5-30x) and chrX (haploid, 2.5-14.6x), the mean spread of observed
     * precision at a claimed score falls from 0.348 for raw GQ to 0.266.
     *
     * Dividing by depth alone does **not** work, and looks like it does if only one
     * ploidy is checked: GQ/DP halves the spread on the diploid series (0.101 to 0.050)
     * and makes the pooled figure worse than doing nothing (0.348 to 0.496). It removes
     * the smaller axis, leaves the larger one, and compresses the range so that what
     * remains does more damage.
     *
     * Each haplotype of `called` is taken to contribute its mixture weight's share of
     * the reads, and such a read to have rel 1 for its own allele and 0 for every other.
     * Both terms of the difference are kept: reads from an allele the runner-up also
     * carries make the gap *smaller*, and dropping them overstates a heterozygote's
     * achievable gap by about a quarter.
     *
     * **The reads' own e_r is deliberately not used**, which is the easy thing to get
     * wrong here -- the first version of this used it. An ideal read is well-fitting
     * *and* well-mapped, so the denominator is built at the mismap floor. Using each
     * read's own e_r puts the site's unreliability into both sides of the ratio, where
     * it cancels: a window of MAPQ-0 reads offers -ln(0.7) = 0.36 per read against
     * -ln(0.02) = 3.91, so a badly mapped site would be scored against a denominator
     * small enough to make a weak call look strong. Measured over the titration, that
     * version scored 0.427 against 0.347 for raw GQ -- worse than not normalising at
     * all -- where the floor version scores 0.260.
     *
     * Deliberately excludes the depth term. That term judges the read *count* against
     * what a genotype predicts and is already scale-free; folding it in here would put
     * a quantity that does not vanish under a perfect pileup into a denominator that is
     * supposed to represent one.
     *
     * Returns 0 when there are no reads or the two genotypes are equal, so callers must
     * guard the division rather than assume a positive denominator.
     */
    double achievable_gap(const vector<int>& called, const vector<int>& runner_up) const;

    /// Every genotype of the given ploidy over num_alleles alleles, as sorted
    /// non-decreasing index multisets in VCF genotype-ordering (colex) order, so
    /// a genotype's position in the returned vector is its GL field index.
    static vector<vector<int>> enumerate_genotypes(size_t num_alleles, int ploidy);

    /// genotype_likelihood over every genotype of this ploidy, in VCF GL order.
    vector<pair<vector<int>, double>> score_genotypes(int ploidy) const;

    /// Write the matrix as TSV for debugging. One row per read.
    void dump(ostream& out, const string& site_name) const;

    /// Populate the matrix. Only for AlleleReadLikelihoodsBuilder.
    void set_contents(size_t n_reads, size_t n_alleles, vector<double>&& matrix,
                      vector<double>&& mismap, vector<double>&& best_ln,
                      vector<string>&& names, size_t unplaceable);

    /// The --mismap-min floor these reads were clamped to. Only achievable_gap uses it,
    /// and only as the reliability an *ideal* read would have; see there for why the
    /// reads' own e_r is the wrong thing to put in that denominator.
    void set_mismap_floor(double floor) { this->mismap_floor = floor; }

private:
    /// Expected share of this site's reads per haplotype of the genotype: flat 1/|G|
    /// unless lengths were supplied, in which case each haplotype is weighted by the
    /// sequence no other member of the genotype carries. Shared by
    /// genotype_likelihood and achievable_gap so the two cannot drift apart -- a gap
    /// measured against different weights than the likelihood it normalises would be
    /// meaningless.
    vector<double> mixture_weights(const vector<int>& genotype) const;

    /// Row major, n_reads * n_alleles, every entry in [0,1], row max exactly 1.
    vector<double> matrix;
    vector<double> read_mismap_prob;
    /// --mismap-min, the reliability an ideal read would have. achievable_gap only.
    double mismap_floor = 0.02;
    /// sum_r (1 - e_r), filled in by set_contents.
    double effective_read_total = 0.0;
    vector<size_t> allele_lengths;
    vector<vector<size_t>> unique_lengths;
    double mean_read_length = 0.0;
    vector<size_t> traversal_lengths;
    double depth_rate = 0.0;
    bool depth_effective = true;
    double depth_read_length = 0.0;
    double depth_weight = 0.0;
    vector<double> read_best_ln;
    vector<string> read_names;
    size_t n_reads = 0;
    size_t n_alleles = 0;
    size_t unplaceable = 0;
};

/**
 * Accumulates raw per-read scores and produces a normalised AlleleReadLikelihoods.
 *
 * Reads are added one at a time with their raw (unnormalised) ln-likelihoods
 * against every allele. A read whose every entry is -inf placed on nothing at
 * all and is dropped, because normalising it would divide by zero; those are
 * counted rather than silently discarded.
 */
class AlleleReadLikelihoodsBuilder {
public:
    /// Mismapping probabilities are clamped into [min_mismap, max_mismap].
    ///
    /// The upper clamp is load-bearing rather than hygiene: many mappers use
    /// MAPQ 0 to mean "multi-mapping" rather than P(wrong) = 1, and an e_r of 1
    /// would collapse the read's term to ln(1) = 0 for every genotype, silently
    /// contributing nothing at all. That may even be the desired behaviour for a
    /// MAPQ 0 read, but it should be a deliberate clamp rather than an artefact
    /// of the phred conversion. It also keeps the per-read term's log finite.
    AlleleReadLikelihoodsBuilder(size_t num_alleles, double min_mismap = 0.01,
                                 double max_mismap = 0.1);

    /// Add a read. raw_ln_likelihood must have one entry per allele and may
    /// contain -inf for alleles that cannot place the read.
    /// read_length feeds the mean R used by the length-weighted mixture. Zero
    /// means "unknown"; if every read is unknown the mixture stays flat.
    void add_read(const vector<double>& raw_ln_likelihood, double mismap_prob,
                  const string& name = "", size_t read_length = 0);

    size_t num_reads_added() const { return rows.size(); }
    size_t num_unplaceable() const { return unplaceable; }


    /// Allele lengths for the length-weighted mixture, indexed by allele. The
    /// read length is accumulated from the reads themselves, so only this is
    /// needed from the caller. Carried through build().
    void set_allele_lengths(vector<size_t> lengths) {
        allele_lengths = std::move(lengths);
    }

    /// See AlleleReadLikelihoods::set_unique_lengths. Carried through build().
    void set_unique_lengths(vector<vector<size_t>> lengths) {
        unique_lengths = std::move(lengths);
    }

    /// Normalise every row by its own maximum and produce the matrix.
    AlleleReadLikelihoods build();

private:
    size_t n_alleles;
    double min_mismap;
    double max_mismap;
    vector<size_t> allele_lengths;
    vector<vector<size_t>> unique_lengths;
    double read_length_total = 0.0;
    size_t read_length_count = 0;
    vector<vector<double>> rows;
    vector<double> mismap_probs;
    vector<string> names;
    vector<double> best_lns;
    size_t unplaceable = 0;
};

/**
 * Tuning for GraphAlignedAlleleLikelihoodCalculator.
 *
 * At namespace scope rather than nested in the calculator so that its default
 * member initializers can be used in a defaulted argument.
 */
struct AlleleLikelihoodParams {
    /// Clamps on the MAPQ-derived mismapping probability. See the builder.
    ///
    /// The *floor* is the more consequential of the two on real data, and its
    /// meaning is broader than mismapping. e_r bounds how much one read can veto
    /// an allele: a read fitting allele A perfectly and B not at all costs B
    /// exactly ln(e_r). MAPQ answers "is this read in the right place", not "is
    /// its path through this site right", so a MAPQ 60 read with a locally
    /// misaligned indel still gets e_r ~ 1e-6 and a -13.8 nat veto it has not
    /// earned. Raising the floor caps that veto: read it as
    /// P(this read's evidence at this site is unreliable), of which mismapping is
    /// only one cause and local misalignment is another.
    ///
    /// Default raised from 1e-8 to 0.01, then to 0.05. 1e-8 asserted that a MAPQ 60
    /// read is locally misaligned once in 10^8 sites, which nothing supports; it let a
    /// single read veto an allele by -13.8 nats, and a few misaligned reads then forced
    /// spurious heterozygous calls. On HG002 chr20 against the GIAB draft benchmark,
    /// 0.01 moved 1,493 genotypes -- 94% of them het to hom -- and improved every class.
    ///
    /// 0.02 was then measured once max_mismap_prob had been corrected to 0.5, because
    /// the two clamps interact and a floor chosen under the old cap is not evidence
    /// about the new one. Raising the cap absorbed most of what the floor used to cost
    /// SNVs -- under the old cap the same move cost about four times as much SNV
    /// precision, which is why an earlier pass rejected it.
    ///
    /// Swept at cap 0.5 on chr20. **Small-variant GT F1 and SV F1 peak in different
    /// places, so this is a choice and not an optimum:**
    ///
    ///   floor         0.01     0.02     0.03     0.05     0.10
    ///   4-hap  small 0.9479   0.9490   0.9495   0.9492   0.9449
    ///   4-hap  SV    0.5164   0.5145   0.5121   0.5074   0.5033
    ///   34-hap small 0.9520   0.9547   0.9563   0.9571   0.9535
    ///   34-hap SV    0.4906   0.4912   0.4860   0.4789   0.4619
    ///
    /// 0.02 is chosen because on the 34-haplotype graph it is the only point that
    /// improves **both** metrics against 0.01 -- everything above it buys small-variant
    /// accuracy by selling SV accuracy, at a rate that gets worse the higher it goes.
    /// It is also the best point on an equal-weight sum of all four numbers.
    ///
    /// That choice is not forced. Weighted by record count the answer is 0.05: there
    /// are 94,691 small-variant truth records on chr20 against 1,680 SVs, so its
    /// +0.0024 small-variant F1 is worth ~227 records while its -0.0123 SV F1 costs
    /// ~22. Anyone who cares only about small variants should set 0.05 explicitly.
    ///
    /// The knob itself is indel-shaped: from 0.01 to 0.02 on the 34-haplotype graph,
    /// insertion GT F1 goes 0.8662 -> 0.8764 and deletion 0.8743 -> 0.8883, against SNV
    /// precision 0.9937 -> 0.9934. Above 0.10 everything degrades on both graphs.
    double min_mismap_prob = 0.02;

    /// The *cap* is what stops the model believing MAPQ when MAPQ is low, and how much
    /// it matters is a property of the graph rather than a constant.
    ///
    /// Default raised from 0.1 to 0.5 on measurement. 0.1 was set when the evaluation
    /// graph had 4 haplotypes, where it was correctly measured as inert -- only 6.3% of
    /// reads sat at MAPQ <= 9, so the cap almost never bound. On a 34-haplotype graph
    /// that population is a quarter of the evidence at the sites that go wrong. Extra
    /// haplotypes do not mainly produce *unmappable* reads -- MAPQ 0 gets rarer -- they
    /// produce **two-way ties** between near-identical placements: reads that fit the
    /// graph better than before (identity 0.921 -> 0.965 at those sites) but cannot be
    /// placed. There, MAPQ 1 alone is 23.3% of reads; MAPQ 1 means p(wrong) = 0.79,
    /// which 0.1 understates 7.9x, and across all reads where the cap binds it discards
    /// 8.1x of the mapper's stated doubt. The model then hears confident support for a
    /// second allele and calls a heterozygote -- 95% of the spurious calls were het.
    ///
    /// At 0.5, on HG002 chr20 against the GIAB draft benchmark: false-positive SNVs
    /// 1,597 -> 443 on the 34-haplotype graph (94% of the excess over the 4-haplotype
    /// one), SNV precision 0.9776 -> 0.9937, overall GT F1 0.9460 -> 0.9520. On the
    /// 4-haplotype graph it is neutral -- 375 -> 376 false SNVs, F1 0.9482 -> 0.9479 --
    /// so it costs nothing where it does not apply. The effect saturates above 0.5 (0.9
    /// gives 0.9517), because beyond that the cap only reaches MAPQ 0-3, so the exact
    /// value is not delicate.
    ///
    /// The cap cannot simply be removed, and the reason is structural rather than
    /// empirical: at e_r = 1 the per-read term becomes log((1-1)*mixture + 1) = 0 for
    /// every genotype, so the read contributes nothing to any of them and silently
    /// disappears. Many mappers also use MAPQ 0 to mean "multi-mapping" rather than
    /// literally P(wrong) = 1. So the cap has to stay strictly below 1; 0.5 says
    /// believe the mapper, but never let one read count for less than half a read.
    ///
    /// The two clamps are not interchangeable. The cap governs *placement* ambiguity
    /// and shows up in SNVs; the floor governs how hard one read may veto an allele and
    /// shows up in indels -- and raising the floor makes SNV precision slightly worse.
    /// Tuning either against an aggregate F1 hides what the other is doing.
    double max_mismap_prob = 0.7;

    /// Turn the mismapping term off entirely, so its contribution can be
    /// measured rather than assumed. With it off, e_r is pinned to the minimum,
    /// which is as close to "trust every read fully" as the model can get while
    /// keeping the log finite.
    bool use_mismap_term = true;
    /// Weight the mixture by each haplotype's expected read contribution at the
    /// site rather than a flat 1/|G|. See AlleleReadLikelihoods::set_length_weights.
    ///
    /// **On by default.** The flat weight asserts each haplotype of a genotype
    /// produced half the site's reads, which is false wherever the alleles differ
    /// in length, and badly so for large events: it lost 94% of heterozygous
    /// deletions above 1 kb and mis-genotyped two thirds of heterozygous insertions
    /// above 1 kb. Correcting it costs nothing measurable elsewhere -- equal-length
    /// alleles give exactly 1/2, so SNV genotype F1 is unchanged to four decimal
    /// places on both graphs tested -- and about 2% of runtime.
    /// The weight counts sequence *unique* to each allele among the genotype's members,
    /// not whole traversal length. That choice was measured against the alternative and
    /// the alternative lost: a traversal includes the site's shared sequence, so at one
    /// 2648 bp deletion the traversals are 296 and 2945 bp -- a ratio of 6.9 where the
    /// reads actually split about 14.6. Reads in shared sequence fit every allele
    /// equally and cancel, so only unique content can move a genotype. There was briefly
    /// a `--length-weight-whole-traversal` flag to select the older form; it is gone,
    /// since the comparison is settled and a knob inside a knob is not worth the surface.
    bool length_weighted_mixture = true;
    /// Weight on ln P(N | G). Zero disables the depth term entirely; the `DR`
    /// diagnostic is still computed, so the observable can be measured before the
    /// model is allowed to act on it.
    ///
    /// **On by default at 0.1.** The read term is conditioned on the reads it was
    /// handed and never asks whether that many reads should be there, which is why
    /// collapsed-repeat pile-ups survived it and why the Poisson caller led on large
    /// heterozygous deletions. Measured across the full five-arm matrix on two
    /// chromosomes and two graphs: structural-variant F1 rises on 4 of 4 datasets for
    /// haplotype enumeration (+0.0067 to +0.0108) and 3 of 4 for support enumeration,
    /// with the fourth flat; small-variant genotype F1 is unchanged to four decimal
    /// places on every read arm; and both Poisson arms are byte-identical, as they
    /// must be. Heterozygous deletion recall above 1 kb roughly doubles.
    ///
    /// 0.1 rather than a heavier weight, from a 27-point grid searched on chr20 and
    /// validated on chr6. The weight is unimodal and turns over sharply above 0.25 --
    /// false positives climb far faster than true ones -- and the optimum is
    /// graph-dependent in a consistent way: 4-haplotype graphs prefer 0.25 and
    /// 34-haplotype graphs prefer 0.1, two for two. 0.1 ties on mean F1, wins on
    /// precision, and wins on the richer graph, which is the direction pangenomes are
    /// going. It gives back some heterozygous deletion recall against 0.25; raise it if
    /// that class is the objective.
    double depth_weight = 0.1;
    /// Count reads toward depth in proportion to `1 - e_r` -- the probability the
    /// read came from this locus at all -- rather than one apiece.
    ///
    /// **On by default.** A read the mapper places with MAPQ 0 is evidence that
    /// *something* is here, not that a read is here: the read term already believes
    /// it only to the extent of `1 - e_r`, and counting it as a whole read of depth
    /// asserts precisely what that term declines to. The same weighting is applied
    /// when the local rate is measured, so a site whose mapping quality matches its
    /// neighbourhood's is unaffected -- the correction is relative, and it moves
    /// only where a site is more or less ambiguously mapped than the sequence
    /// around it.
    bool depth_effective_reads = true;
    /// Ploidy to assume when the caller does not supply one. Only a fallback: the
    /// site's own ploidy is passed to `compute` and used in preference, because
    /// `vg call` varies ploidy by contig (-d, --ploidy-regex) and a haploid region
    /// scored against a diploid rate gets a lambda that is wrong by the ploidy
    /// ratio, in a term that is otherwise one of the better-behaved parts of the
    /// model.
    int depth_ploidy = 2;
};

/**
 * Produces the reads x alleles matrix for a site.
 */
class AlleleLikelihoodCalculator {
public:
    virtual ~AlleleLikelihoodCalculator() = default;

    /// Build the matrix for one site. traversals are the candidate alleles, in
    /// the order the caller will genotype them. `ploidy` is the site's ploidy, which
    /// the depth term needs to turn a local read count into a per-haplotype rate;
    /// `vg call` varies it by contig through -d and --ploidy-regex.
    virtual AlleleReadLikelihoods compute(const Snarl& snarl,
                                          const vector<SnarlTraversal>& traversals,
                                          int ploidy) = 0;
};

/**
 * Scores each read against each allele from the read's *existing* alignment in
 * the graph, by walking the read's path against the allele's path and accounting
 * the differences. No dynamic programming.
 *
 * This is not merely a speed shortcut. In a variation graph the snarl
 * decomposition already is a multiple alignment: two traversals of a snarl share
 * its boundary nodes by construction, so they are aligned to each other through
 * the topology, and "how many differences between this read and this allele" is
 * already well defined structurally. Re-deriving it with DP could find a higher
 * scoring alignment corresponding to no path in the graph, which is not the
 * quantity we want -- our alleles *are* paths, and we need P(read | this path).
 * Counting differences also tells us exactly which read bases mismatch, so each
 * mismatch is charged its own base quality rather than a length average.
 *
 * The cost it carries is that it inherits the mapper's placement: if a read was
 * placed wrongly we score against the wrong path and cannot recover. That is
 * what an optional realigning calculator would be for.
 *
 * ## The scoring window
 *
 * The window is a property of the read and the site, never of the allele. It is
 * the span of read bases the read's own alignment places inside the site, and
 * **every allele is scored over that same span**. Read bases an allele cannot
 * place are charged as insertions rather than omitted.
 *
 * This is an invariant, not a preference. Scoring allele *a* over 150 read bases
 * and allele *b* over only the 100 it can place makes their difference contain
 * 50 x match-reward: a fabricated likelihood ratio of arbitrary magnitude that
 * merely happens to point the right way. P(read | allele) may differ between
 * alleles only in how well the *same* observed bases are explained.
 */
class GraphAlignedAlleleLikelihoodCalculator : public AlleleLikelihoodCalculator {
public:

    /// Defined at namespace scope; aliased here for readability.
    using Params = AlleleLikelihoodParams;

    /**
     * Two scorers are needed, not one.
     *
     * A quality-adjusted scorer charges each mismatch its own base quality,
     * which is the whole point of scoring per read. But it is inapplicable to a
     * read with no base qualities, which happens with GAF that has no quality
     * column. Rather than fabricate qualities or silently mis-score, pick per
     * read: qual_scorer when the read has qualities, plain_scorer when it does
     * not. Their log bases differ, but that is harmless here because every row
     * is normalised by its own maximum before use, so a row never mixes scorers
     * and cross-row comparisons are dimensionless.
     */
    GraphAlignedAlleleLikelihoodCalculator(const PathHandleGraph& graph,
                                           SnarlManager& snarl_manager,
                                           const SiteReadSource& read_source,
                                           const EditAlignmentScorer& qual_scorer,
                                           const EditAlignmentScorer& plain_scorer,
                                           const Params& params = Params());

    AlleleReadLikelihoods compute(const Snarl& snarl,
                                  const vector<SnarlTraversal>& traversals,
                                  int ploidy) override;

protected:

    /// One node visit on an allele, with its sequence materialised.
    struct AlleleStep {
        nid_t node_id;
        bool backward;
        string sequence;
    };

    /// One node visit by a read inside the site, tied back to its mapping.
    struct ReadStep {
        nid_t node_id;
        bool backward;
        /// Offset into the read sequence at which this mapping's bases start.
        size_t read_offset;
        /// How many read bases this mapping consumes.
        size_t read_length;
        const Mapping* mapping;
    };

    /// Materialise an allele's node visits and sequences. Per allele, not per
    /// (read, allele), so cheap enough to do once per site.
    vector<AlleleStep> get_allele_steps(const SnarlTraversal& traversal) const;

    /// Extract the read's visits inside the site, in read order. Returns false
    /// if the read has no informative overlap: a read touching only the site's
    /// boundary nodes would contribute an identical constant to every allele,
    /// since every allele shares those nodes by construction.
    ///
    /// That is NOT the same as a read failing to place on *some* allele, which
    /// is informative and must be kept. Conflating the two would systematically
    /// destroy deletion and structural variant genotyping.
    bool get_read_steps(const Alignment& aln, const unordered_set<nid_t>& site_nodes,
                        const unordered_set<nid_t>& boundary_nodes,
                        vector<ReadStep>& steps_out) const;

    /// True if the read traverses the site against the direction the alleles read
    /// in, so it must be reverse-complemented before being compared to them.
    /// Decided by vote over shared nodes rather than from a single step.
    bool read_is_reverse_of_alleles(const vector<ReadStep>& read_steps,
                                   const unordered_map<nid_t, bool>& allele_orientations) const;

    /// Score the read's own edits on a node the read and the allele share.
    int32_t score_shared_node(const Alignment& aln, const ReadStep& step,
                              const EditAlignmentScorer& read_scorer) const;

    /// Score read bases against allele bases of the same length, base by base,
    /// charging each mismatch its own quality.
    int32_t score_substitution(const Alignment& aln, size_t read_offset, size_t length,
                               const string& allele_bases, size_t allele_offset,
                               const EditAlignmentScorer& read_scorer) const;

    /// Score one read against one allele over the read's window.
    int32_t score_read_against_allele(const Alignment& aln, const vector<ReadStep>& read_steps,
                                      const vector<AlleleStep>& allele_steps,
                                      const EditAlignmentScorer& read_scorer,
                                      bool& placed_out) const;

    /// Reads per position per haplotype, measured over the read source's own fetch
    /// window around this site.
    ///
    /// The window is the right denominator and the snarl is not. A snarl's shared
    /// sequence is its boundary nodes, which are too short to contain a read --
    /// measured on ten large deletions, the number of reads fitting every allele
    /// equally was zero at nine of them. A window is thousands of node IDs wide, and
    /// crucially it is already fetched and cached to answer the site's own query, so
    /// asking for it a second time is a cache hit rather than new I/O.
    ///
    /// Returns 0 if the source has no window (an in-memory source answers exactly),
    /// which switches the depth term off rather than inventing a rate.
    /// Reads per base pair over the window containing these ranges, weighted by
    /// `1 - e_r` when `depth_effective_reads` is set.
    ///
    /// Deliberately **not** divided by ploidy. The window statistic is a property of
    /// the data; ploidy is a property of how the site is being genotyped. Folding
    /// the second into the first put a genotyping assumption inside a cache keyed
    /// only on node range, so a rate computed under one ploidy could be reused under
    /// another. The division happens at the point of use instead.
    double local_read_rate(const vector<pair<nid_t, nid_t>>& site_ranges) const;

    const PathHandleGraph& graph;
    SnarlManager& snarl_manager;
    const SiteReadSource& read_source;
    mutable unordered_map<size_t, double> window_rate;
    mutable std::mutex window_bp_mutex;
    const EditAlignmentScorer& qual_scorer;
    const EditAlignmentScorer& plain_scorer;
    Params params;
};

}

#endif
