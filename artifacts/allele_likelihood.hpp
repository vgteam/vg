#ifndef VG_ALLELE_LIKELIHOOD_HPP_INCLUDED
#define VG_ALLELE_LIKELIHOOD_HPP_INCLUDED

/** \file allele_likelihood.hpp
 *
 * The read-by-allele likelihood matrix of one site, and the code that builds it
 * from the reads' existing alignments to the graph.
 *
 * Each entry says how well one read fits one candidate allele. The read-likelihood
 * genotyper combines these entries into a likelihood for each genotype. The model
 * is described in artifacts/doc/read-likelihood-direct-genotyping.md.
 */

#include <functional>
#include <iostream>
#include <limits>
#include <atomic>
#include <map>
#include <mutex>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <vg/vg.pb.h>

#include "../src/alignment_scorer.hpp"
#include "../src/handle.hpp"
#include "materialized_site.hpp"
#include "site_read_source.hpp"
#include "site_values.hpp"
#include "../src/snarls.hpp"

#include "anchor.hpp"
#include "read_phasing.hpp"

namespace vg {

using namespace std;

/**
 * A record of the reads seen at a site and of how many reads each of its genotypes predicts.
 *
 * A genotype G predicts lambda_G = rate * sum_{h in G} (T_h + R - 1) reads, where T_h is the
 * length of allele h without the site's two boundary nodes and R is the read length. A copy kept
 * after the matrix is gone can still give DR for a genotype other than the call.
 */
struct DepthContext {
    /// T_h for each allele, in matrix column order.
    ///
    /// An allele past the end, or a negative marker, has T_h = 0.
    vector<size_t> lengths;
    /// The read starts per base per haplotype near the site.
    double rate = 0.0;
    /// The read length R.
    double read_length = 0.0;
    /// The read count N seen at the site, which need not be a whole number.
    double observed = 0.0;

    /// Return lambda_G, the number of reads the genotype predicts.
    double expected_reads(const vector<int>& genotype) const;

    /// Return N / lambda_G, or -1 where lambda_G is 0.
    double ratio(const vector<int>& genotype) const;
};

/**
 * A matrix of the relative likelihood of each read of one site under each of its alleles.
 *
 * Each row belongs to a read and is divided by its own maximum, so every entry is in [0, 1] and
 * the read's best allele scores 1. The alignment score gives ln P(read | allele) only up to a
 * constant that differs between reads, and dividing by the row maximum removes that constant. It
 * also puts every row on the scale of the mismapping term, whose background likelihood is then 1,
 * and lets rows from scorers with different log bases share one matrix.
 *
 * The values can be used to compare genotypes, to choose one, to compute GQ, or to form a
 * posterior. They are not calibrated absolute probabilities.
 */
class AlleleReadLikelihoods {
public:

    AlleleReadLikelihoods() = default;

    /// Get the number of reads, which is the number of rows.
    size_t num_reads() const { return n_reads; }
    /// Get the number of alleles, which is the number of columns.
    size_t num_alleles() const { return n_alleles; }

    /**
     * Get the likelihood of read r under allele a, relative to the best allele of read r.
     *
     * It is 0 where allele a fits the read so much worse than the best that the ratio
     * underflows. It is also 0 where allele a has no node visits, as the empty traversal of a
     * star allele has, and cannot place the read.
     */
    double rel(size_t r, size_t a) const;

    /// Get the mismapping probability e_r of read r, from its MAPQ, clamped inside (0, 1).
    double mismap_prob(size_t r) const;

    /**
     * Get the best ln-likelihood of read r at this site, in nats, which its row was divided by.
     *
     * The genotype likelihood does not use it. It is written as the BL output field.
     */
    double best_ln_likelihood(size_t r) const;

    /// Get the number of reads dropped because they placed on no allele.
    size_t num_unplaceable() const { return unplaceable; }

    /**
     * Weight the mixture by length in place of a flat 1/|G|, and set R.
     *
     * Each haplotype of a genotype is weighted by the number of start positions from which a read
     * of length R overlaps a stretch of length X:
     *
     *     w_h = (X_h + R - 1) / sum_{h' in G} (X_h' + R - 1)
     *
     * X_h is the unique length U_h from `set_unique_lengths`, which the calculator gives every
     * matrix. Without unique lengths it is the full length L_h of the allele, from `lengths`.
     *
     * R should be the mean length of reads in the neighbourhood rather than of the reads at this
     * site, since a long read overlaps more sites and so is over-represented at each. An empty
     * `lengths` or a zero R gives the flat 1/|G|.
     *
     * @param lengths The length L_h of each allele.
     * @param read_length The mean read length R.
     */
    void set_length_weights(vector<size_t> lengths, double read_length) {
        this->allele_lengths = std::move(lengths);
        this->mean_read_length = read_length;
    }

    /**
     * Weight the mixture by the sequence unique to each allele (see `set_length_weights`).
     *
     * Reads that lie in sequence both alleles share fit both equally and cannot change a genotype
     * comparison. So the weight uses U_h, the total length of the nodes that the allele of
     * haplotype h visits and the other allele of the genotype does not, in either orientation.
     *
     * @param lengths U_h for each allele against each other allele, as `lengths[a][b]`.
     */
    void set_unique_lengths(vector<vector<size_t>> lengths) {
        this->unique_lengths = std::move(lengths);
    }

    /**
     * Set up the depth term, which asks whether genotype G predicts the number of reads seen, N.
     *
     * The term is w * ln Poisson(N ; lambda_G) - ln Z_w(lambda_G), where Z_w(lambda) normalises
     * Poisson(n ; lambda)^w over n >= 0 (see `genotype_likelihood`), and
     *
     *     lambda_G = rate * sum_{h in G} (T_h + R - 1)
     *
     * The term is off unless `weight` is positive.
     *
     * @param lengths T_h for each allele, its length without the site's two boundary nodes, since
     * a read inside a boundary node is not a row of this matrix.
     * @param rate Read starts per base per haplotype near the site.
     * @param read_length R.
     * @param effective_count Whether N counts each read as 1 - e_r. `rate` must count reads the
     * same way.
     */
    void set_depth_context(vector<size_t> lengths, double rate, double read_length,
                           double weight, bool effective_count = true) {
        this->depth.lengths = std::move(lengths);
        this->depth.rate = rate;
        this->depth.read_length = read_length;
        this->depth_weight = weight;
        this->depth_effective = effective_count;
        count_observed();
    }

    /// Return whether the depth term is on.
    bool uses_depth_term() const {
        return depth_weight > 0.0 && depth.rate > 0.0 && !depth.lengths.empty();
    }

    /// Return the expected number of reads at this site under the genotype, lambda_G.
    double expected_reads(const vector<int>& genotype) const {
        return depth.expected_reads(genotype);
    }

    /// Return whether `set_depth_context` has supplied traversal lengths.
    bool has_depth_context() const { return !depth.lengths.empty(); }

    /// Get the depth context, with N as `observed_reads` gives it.
    const DepthContext& depth_context() const { return depth; }

    /**
     * Return the read count N that the depth term compares with lambda_G.
     *
     * N is sum_r (1 - e_r) with `effective_count`, and the number of rows otherwise. Since N need
     * not be a whole number, the depth term treats it with a continuous analogue of the Poisson
     * distribution over read counts n >= 0, using lgamma in place of ln n!.
     */
    double observed_reads() const { return depth.observed; }

    /// Get the length T_h of an allele without the site's boundary nodes, or 0 if the depth
    /// context was not set or the index is out of range.
    size_t traversal_length(size_t allele) const {
        return allele < depth.lengths.size() ? depth.lengths[allele] : 0;
    }

    /**
     * Return the observed over the expected read count, N / lambda_G, for the given genotype.
     *
     * It is -1 where lambda_G is 0, as it is when no read begins in the rate window. It is written
     * as the DR output field whether or not the depth term is on.
     */
    double depth_ratio(const vector<int>& genotype) const { return depth.ratio(genotype); }

    /// Get the mean read length R of the mixture weights, which is also copied into the
    /// evidence kept for read phasing and anchors.
    double mean_read_length_estimate() const { return mean_read_length; }

    /// Set R for the mixture weights, and for the evidence copied from this matrix, without
    /// turning on length weighting.
    ///
    /// The depth term takes its own R from `set_depth_context`.
    void set_mean_read_length(double read_length) {
        this->mean_read_length = read_length;
    }

    /// Return whether `set_length_weights` supplied usable data.
    bool uses_length_weights() const {
        return !allele_lengths.empty() && mean_read_length > 0.0;
    }

    /**
     * Return ln P(reads | G), where G is a multiset of allele indices of size ploidy.
     *
     *   ln P(reads | G) =  sum_r ln [ (1 - e_r) * sum_{h in G} w_h * rel(r,h) + e_r ]
     *                    + w_d * ln Poisson( N ; lambda_G ) - ln Z_{w_d}(lambda_G)
     *
     * r runs over the site's reads and h over the haplotypes of G, so a homozygote counts its
     * allele twice. e_r is the read's mismapping probability and w_h the mixture weight. The rest
     * is the depth term (see `set_depth_context`). It is the log of the density
     * Poisson(n ; lambda)^w_d / Z_{w_d}(lambda) over n >= 0 at the observed count, with Z_{w_d}
     * taken from its normal approximation, ((1 - w_d) / 2) ln(2 pi lambda) - (1/2) ln w_d.
     *
     * Since rel(r,h) lies in [0, 1], the term of each read lies between ln(e_r) and 0, so the floor
     * on e_r limits how much one read can count against a genotype. Reads are treated as
     * independent, so GL and GQ grow over-confident with depth.
     */
    double genotype_likelihood(const vector<int>& genotype) const;

    /**
     * Return the largest difference the read term could give between the called genotype and the
     * runner-up at this site, which is the denominator of GQN.
     *
     * It is the difference for ideal reads. Each haplotype of `called` contributes the share of
     * the reads its mixture weight gives it. Each read fits the allele of its own haplotype with
     * rel 1 and every other allele with 0, and has e_r at the mismapping floor. Dividing the
     * observed difference by this one gives a value in [0, 1] that does not depend on depth or on
     * ploidy.
     *
     * The depth term is left out, since a perfect set of reads does not maximise it.
     *
     * @return The difference, or 0 when there are no reads or the two genotypes are equal, so
     * callers must check before dividing.
     */
    double achievable_gap(const vector<int>& called, const vector<int>& runner_up) const;

    /**
     * Return every genotype of the given ploidy over `num_alleles` alleles.
     *
     * Each genotype is a non-decreasing list of allele indices. They are in VCF genotype order,
     * which is colexicographic, so the position of a genotype in the result is its GL field index.
     */
    static vector<vector<int>> enumerate_genotypes(size_t num_alleles, int ploidy);

    /// Return `genotype_likelihood` for every genotype of this ploidy, in VCF GL order.
    vector<pair<vector<int>, double>> score_genotypes(int ploidy) const;

    /// Write the matrix as TSV for debugging, with one line per read.
    void dump(ostream& out, const string& site_name) const;

    /**
     * Keep the anchor evidence of the site, and fill its rows from this matrix.
     *
     * The calculator attaches it only when anchors are being written.
     *
     * @pre `evidence` has one read for each row, in row order.
     */
    void attach_anchor_evidence(unique_ptr<AnchorSiteEvidence> evidence);

    /**
     * Keep the read-phasing evidence of the site, and fill its rows from this matrix.
     *
     * The calculator attaches it only when read phasing is on and anchors are not being written,
     * since the anchor evidence holds the same rows.
     *
     * @pre `evidence` has one read key for each row, in row order.
     */
    void attach_phase_evidence(unique_ptr<PhaseReadEvidence> evidence);

    /// Give up the anchor evidence, or null if none was attached.
    unique_ptr<AnchorSiteEvidence> take_anchor_evidence() { return std::move(anchor_evidence); }

    /// Give up the read-phasing evidence, or null if none was attached.
    unique_ptr<PhaseReadEvidence> take_phase_evidence() { return std::move(phase_evidence); }

    /// Fill the matrix, for `AlleleReadLikelihoodsBuilder`.
    void set_contents(size_t reads, size_t alleles, vector<double>&& rows,
                      vector<double>&& mismap, vector<double>&& best_ln,
                      vector<string>&& names, size_t dropped);

    /// Set the floor on e_r, which `achievable_gap` gives its ideal reads.
    void set_mismap_floor(double floor) { this->mismap_floor = floor; }

private:
    /**
     * Return the mixture weight of each haplotype of the genotype (see `set_length_weights`).
     *
     * The weights are a flat 1/|G| unless lengths were supplied. `genotype_likelihood` and
     * `achievable_gap` must use the same weights.
     */
    vector<double> mixture_weights(const vector<int>& genotype) const;

    /// Return the length X_h that weights haplotype `i` of the genotype in the mixture.
    double weighted_length(const vector<int>& genotype, size_t i) const;

    /// Set N in `depth` from the reads and from whether they count as 1 - e_r.
    void count_observed();

    /// The anchor evidence of the site, or null.
    unique_ptr<AnchorSiteEvidence> anchor_evidence;
    /// The read-phasing evidence of the site, or null.
    unique_ptr<PhaseReadEvidence> phase_evidence;

    /// The relative likelihoods, row major, n_reads * n_alleles, with a maximum of 1 in each row.
    vector<double> matrix;
    /// The mismapping probability e_r of each read.
    vector<double> read_mismap_prob;
    /// The floor on e_r, used by `achievable_gap` only.
    double mismap_floor = 0.02;
    /// The expected number of the reads that came from this locus, sum_r (1 - e_r).
    double effective_read_total = 0.0;
    vector<size_t> allele_lengths;
    vector<vector<size_t>> unique_lengths;
    double mean_read_length = 0.0;
    /// The depth context.
    ///
    /// Its N is `effective_read_total` under `depth_effective`, and `n_reads` otherwise.
    DepthContext depth;
    bool depth_effective = true;
    double depth_weight = 0.0;
    vector<double> read_best_ln;
    vector<string> read_names;
    size_t n_reads = 0;
    size_t n_alleles = 0;
    size_t unplaceable = 0;
};

/**
 * A builder that collects the raw scores of each read and makes an `AlleleReadLikelihoods`.
 *
 * Reads are added one at a time, each with its raw ln-likelihood under every allele. A read whose
 * entries are all -inf placed on no allele and has no row maximum to divide by, so it is dropped
 * and counted.
 */
class AlleleReadLikelihoodsBuilder {
public:
    /**
     * Make a builder for a site with `num_alleles` alleles.
     *
     * Mismapping probabilities are clamped into [min_mismap, max_mismap]. The upper clamp matters
     * because many mappers give MAPQ 0 to a read with several equally good placements. Its
     * unclamped e_r of 1 would make the read's term 0 under every genotype, so the read would
     * count for nothing.
     */
    explicit AlleleReadLikelihoodsBuilder(size_t num_alleles, double min_mismap = 0.02,
                                          double max_mismap = 0.95);

    /**
     * Add a read.
     *
     * @param raw_ln_likelihood The ln-likelihood of the read under each allele, or -inf for an
     * allele that cannot place the read.
     * @param read_length The length of the read, for the mean R of the length-weighted mixture,
     * or 0 if it is unknown. If every read's length is unknown, the mixture stays flat.
     * @param start Where the read's alignment begins. With the name it identifies the alignment,
     * since paired mates share a name, and `build` orders the rows by the two.
     * @return False when the read placed on no allele and was dropped, so that a caller keeping
     * anything alongside the rows can stay in step with them.
     */
    bool add_read(const vector<double>& raw_ln_likelihood, double mismap_prob,
                  const string& name = "", size_t read_length = 0,
                  const Position* start = nullptr);

    /// Set the length of each allele for the length-weighted mixture.
    ///
    /// The mean read length comes from the reads themselves.
    void set_allele_lengths(vector<size_t> lengths) {
        allele_lengths = std::move(lengths);
    }

    /// Set the unique lengths for the mixture (see `AlleleReadLikelihoods::set_unique_lengths`).
    void set_unique_lengths(vector<vector<size_t>> lengths) {
        unique_lengths = std::move(lengths);
    }

    /**
     * Divide every row by its own maximum and make the matrix.
     *
     * The rows are put in a canonical order, by read name, then by where the alignment begins,
     * then by the row's own values. Then the matrix, and every sum over its reads, does not depend
     * on the order the reads were added in, which can change with the fetch window of the read
     * source.
     */
    AlleleReadLikelihoods build();

    /**
     * Get, for each row of the built matrix, the index of its read among the reads `add_read` kept.
     *
     * A caller that kept something per read alongside the rows reorders it by this.
     */
    const vector<size_t>& row_order() const {
        return order;
    }

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
    /// Where the alignment of each kept read begins, as node, offset and strand.
    vector<tuple<nid_t, int64_t, bool>> starts;
    vector<size_t> order;
    size_t unplaceable = 0;
};

/**
 * The tuning of a `GraphAlignedAlleleLikelihoodCalculator`.
 *
 * It is at namespace scope rather than nested in the calculator, so that its default member
 * initializers can be used in a defaulted argument.
 */
struct AlleleLikelihoodParams {
    /**
     * The nats added to a read's log-likelihood for each gap in which the read has bases the
     * allele lacks.
     *
     * A positive value makes extra read bases count against an allele less than missing ones. It
     * is applied after the integer alignment score is converted to nats, so it can be fractional.
     */
    double insertion_gap_nats = 0.0;

    /**
     * The rate at which a read shows an insertion that the sequence it was read from lacks, per
     * base of the repeat tract it shows it in.
     *
     * An allele whose step sequences differ from its walk by one indel, as an off-panel indel's
     * do, is scored with this rate and `indel_deletion_rate`. A read covering the indel's tract
     * and a base on each side counts for or against the allele by these rates, and any other read
     * for neither (see `observe_indel`). The rates were measured on short reads. On long reads the
     * mismapping floor binds before them outside long homopolymers.
     */
    double indel_insertion_rate = 2.8e-5;
    /// The rate at which a read shows a deletion that the sequence it was read from lacks, per
    /// base of the repeat tract it shows it in (see `indel_insertion_rate`).
    double indel_deletion_rate = 1.3e-5;

    /// Choose each read's pairing with an allele by optimal pairing rather than greedy pairing
    /// (see `GraphAlignedAlleleLikelihoodCalculator`).
    bool optimal_pairing = false;

    /**
     * The floor on the mismapping probability e_r.
     *
     * A read that fits allele A perfectly and allele B not at all lowers the likelihood of B by
     * at most -ln(e_r), so the floor limits how strongly one read can count against an allele.
     * MAPQ measures whether the read is at the right locus, and not whether its path through this
     * site is right, so the floor also stands for a well-mapped read that is misaligned locally.
     */
    double min_mismap_prob = 0.02;

    /**
     * The ceiling on e_r.
     *
     * It applies to reads with MAPQ 0 or close to it, and decides how much such a read still
     * counts. It must stay below 1, since at e_r = 1 the read's term is 0 under every genotype.
     */
    double max_mismap_prob = 0.95;

    /// Use the mismapping term.
    ///
    /// When false, every e_r is the floor, the closest the model can come to trusting every read
    /// fully while keeping the log finite.
    bool use_mismap_term = true;

    /// Weight each haplotype of a genotype by its share of the reads that can tell the genotype's
    /// alleles apart, from the sequence unique to its allele, rather than by a flat 1/|G| (see
    /// `AlleleReadLikelihoods::set_length_weights`).
    bool length_weighted_mixture = true;

    /// The weight of the depth term, ln P(N | G), or 0 to turn the term off.
    ///
    /// DR is computed either way.
    double depth_weight = 0.1;

    /// Count each read towards depth as 1 - e_r, the probability that it came from this locus,
    /// rather than as 1. The local rate is counted the same way.
    bool depth_effective_reads = true;

    /// The ploidy used for the depth rate when `compute` gets a region ploidy of 0 or less.
    ///
    /// If this is 0 or less too, the rate window is not measured.
    int depth_ploidy = 2;

    /// Collect the anchor evidence of each read, resolving its position while its alignment is
    /// available.
    bool collect_anchors = false;
    /// Where to count the position resolutions of this calculator, or null to not count them.
    ///
    /// It is not owned (see `AnchorParams::counters`).
    AnchorCounters* anchor_counters = nullptr;

    /// Keep each read's row of relative likelihoods for read phasing.
    ///
    /// Phasing needs only the rows, not the positions that `collect_anchors` also resolves.
    bool collect_read_phasing = false;
};

/**
 * A calculator that makes the matrix of reads by alleles for a site.
 */
class AlleleLikelihoodCalculator {
public:
    virtual ~AlleleLikelihoodCalculator() = default;

    /**
     * Make the matrix for one site.
     *
     * Nothing in the matrix depends on the ploidy the site is then genotyped at. An allele with
     * step sequences is scored against them rather than against the nodes of its walk.
     *
     * @param site The site, whose alleles are in the order the caller will genotype them.
     * @param fallback_nodes Nodes to place the site's depth-rate window by, most preferred
     * first, when neither bound of the site is on a reference path.
     * @param region_ploidy The number of the sample's haplotypes in the region around the site.
     * The depth term divides the local read rate by it to get a rate per haplotype. It is the
     * site's own ploidy at a top-level site. At a nested site that only some of its parent's
     * alleles cross, it is more, since the reads counted near the site come from all of them.
     */
    virtual AlleleReadLikelihoods compute(const MaterializedSite& site,
                                          const vector<nid_t>& fallback_nodes,
                                          int region_ploidy) = 0;
};

/**
 * A calculator that scores each read against each allele from the read's existing alignment to
 * the graph, rather than by aligning the read to each allele again.
 *
 * ## Pairings
 *
 * A read's placement and an allele are both sequences of node visits. A read is scored against
 * an allele through a *pairing* of the two sequences, in order. Each read visit inside the site
 * is either paired with an allele visit or left unpaired. The pairing's score is the sum of these
 * parts:
 *
 *   - A read visit paired with the same allele visit scores the read's own edits in that node.
 *   - A read visit that the allele never makes, paired with an allele visit that the read never
 *     makes, is a substitution. The read's bases in its node are compared one by one with the
 *     allele node's sequence, from the first base of each, plus a gap for the difference in
 *     length.
 *   - A run of read visits left unpaired is an insertion, and a run of allele visits skipped
 *     between two pairs is a deletion, each scored as one gap. A read visit with no bases, where
 *     the read deletes its whole node, adds nothing when left unpaired and does not break a run.
 *
 * Allele visits before the read's first same-visit pair, or after its last pair, lie outside the
 * read and score nothing. Before that first same-visit pair, each unpaired read visit is a gap of
 * its own, since there is no pair yet for an insertion to extend from.
 *
 * Every read base inside the site is scored under every allele. So all alleles are scored over
 * the same read bases, the read's *scoring window*, and differ only in how well they explain
 * them.
 *
 * ## Two ways to choose the pairing
 *
 *   - *Greedy pairing* (`score_by_greedy_pairing`, the default) makes one pass along the read's
 *     visits. Each is paired with the next occurrence of the same visit in the allele after the
 *     last pair, and a pair is never revised.
 *   - *Optimal pairing* (`score_by_optimal_pairing`) finds the highest-scoring pairing that the
 *     rules above allow, by dynamic programming over the two sequences. At large sites the search
 *     is restricted to a band.
 *
 * Both score the pairing they find by the rules above. Inside a node that the read and the allele
 * share, both use the mapper's edits and neither aligns bases again.
 */
class GraphAlignedAlleleLikelihoodCalculator : public AlleleLikelihoodCalculator {
public:

    using Params = AlleleLikelihoodParams;

    /**
     * Make a calculator that scores the reads of `read_source` on `graph`.
     *
     * The two scorers have different log bases. That does not matter, because each row is divided
     * by its own maximum, and a row is scored by one scorer only.
     *
     * @param qual_scorer Scores reads that have base qualities, charging each mismatch at its
     * base quality.
     * @param plain_scorer Scores reads without base qualities, such as reads from a GAF file with
     * no quality column.
     */
    GraphAlignedAlleleLikelihoodCalculator(const PathHandleGraph& graph,
                                           const SiteReadSource& read_source,
                                           const EditAlignmentScorer& qual_scorer,
                                           const EditAlignmentScorer& plain_scorer,
                                           const Params& params = Params());

    AlleleReadLikelihoods compute(const MaterializedSite& site,
                                  const vector<nid_t>& fallback_nodes,
                                  int region_ploidy) override;

    /**
     * Place rate windows on reference paths.
     *
     * Until this is called, windows are blocks of node IDs (see `local_read_stats`).
     *
     * @param position_graph The graph the calculator was built on, or a view of it with the same
     * nodes.
     * @param reference_paths The reference paths of `position_graph` to place windows on.
     */
    void set_rate_reference(const PathPositionHandleGraph* position_graph,
                            const vector<path_handle_t>& reference_paths);

    /// Get the number of sites whose rate window fell back to a block of node IDs.
    size_t rate_id_fallbacks() const {
        return id_fallbacks.load();
    }

    /// The width, in bp of reference, of a rate-window bucket. A window is three buckets.
    static const int64_t RATE_BUCKET = 16384;

    /// The width, in node IDs, of the fallback rate window.
    static const nid_t RATE_ID_WINDOW = 4096;

protected:

    /// One visit of an allele to a node, with its sequence.
    struct AlleleStep {
        nid_t node_id;
        bool backward;
        string sequence;
        /// Whether `sequence` is the allele's own step sequence and differs from the node's (see
        /// `StepSequences`).
        ///
        /// A read's visit to the node is then scored against it, and not against the node through
        /// the read's alignment (see `score_shared_node`).
        bool own_sequence = false;
    };

    /// One visit of a read to a node inside the site, tied back to its mapping.
    struct ReadStep {
        nid_t node_id;
        bool backward;
        /// The offset in the read sequence at which this mapping's bases start.
        size_t read_offset;
        /// The number of read bases this mapping consumes.
        size_t read_length;
        const Mapping* mapping;
    };

    /// A record of how the alleles of a site inside nodes change the one walk they all take.
    struct WalkChanges;

    /// The two pins of a site, where the anchors place each read in the site's bounds.
    struct SitePins;

    /// A table of the alleles of one site, as node visits with sequences, ready for scoring reads.
    struct SiteAlleles;

    /// A dynamic-programming grid that finds the optimal pairing of the visits of one read with
    /// those of one allele.
    class PairingGrid;

    /// Return the node visits of an allele with their sequences, which are its step sequences
    /// `own` when not null, or else its nodes'.
    vector<AlleleStep> get_allele_steps(const Traversal& walk,
                                        const StepSequences* own = nullptr) const;

    /// Return the number of bases spelled by the steps [from, to) of an allele.
    static size_t spelled_length(const vector<AlleleStep>& steps, size_t from, size_t to);

    /**
     * Find the visits of a read inside the site, in read order.
     *
     * A read that fails to place on some allele is kept, since that is evidence against the
     * allele.
     *
     * @return False if the read cannot tell the alleles apart, because it lies within one
     * boundary node, which every allele shares.
     */
    bool get_read_steps(const SiteRead& read, const unordered_set<nid_t>& site_nodes,
                        const unordered_set<nid_t>& boundary_nodes,
                        vector<ReadStep>& steps_out) const;

    /**
     * Return whether the read traverses the site against the direction the alleles read in.
     *
     * Such a read must be reverse-complemented before being compared to the alleles. We decide by
     * a vote over the nodes in `allele_orientations`, the nodes every allele visits in the same
     * orientation. A tie leaves the read as it is.
     */
    bool read_is_reverse_of_alleles(const vector<ReadStep>& read_steps,
                                    const unordered_map<nid_t, bool>& allele_orientations) const;

    /**
     * Return the alignment of a read in the direction the alleles read in, or null if the read
     * cannot tell the alleles apart.
     *
     * @param read_steps Gets the read's visits inside the site, in that direction.
     * @param flipped Holds the flipped alignment of a read that needs flipping.
     */
    const Alignment* orient_read(const SiteRead& read, const unordered_set<nid_t>& site_nodes,
                                 const unordered_set<nid_t>& boundary_nodes,
                                 const SiteAlleles& alleles, vector<ReadStep>& read_steps,
                                 Alignment& flipped) const;

    /**
     * Score the read's own edits on a node the read and the allele share.
     *
     * Where `allele_step` has its own sequence, the read's bases are compared with that sequence
     * instead of with the node's. Each aligned base of the read is compared with the base of the
     * step sequence at the same offset in the node, and each gap of the alignment is charged as
     * it is. The offsets are the node's, so a step sequence that changes the node's length is
     * compared correctly only up to the change.
     *
     * @param nat_adjust Gets real-valued corrections added to it, which the integer score cannot
     * carry. The caller adds it after converting the score to nats.
     */
    int32_t score_shared_node(const Alignment& aln, const ReadStep& step,
                              const EditAlignmentScorer& read_scorer, double& nat_adjust,
                              const AlleleStep* allele_step = nullptr) const;

    /// Score `length` of the read's own bases, from `read_offset`, against as many allele bases,
    /// from `allele_offset`, base by base, charging each mismatch at its own base quality.
    int32_t score_substitution(const Alignment& aln, size_t read_offset, size_t length,
                               const string& allele_bases, size_t allele_offset,
                               const EditAlignmentScorer& read_scorer) const;

    /**
     * Score pairing a read visit with a different allele visit, as a substitution.
     *
     * @param nat_adjust Gets `insertion_gap_nats` added to it if the read has more bases.
     */
    int32_t score_node_substitution(const Alignment& aln, const ReadStep& read_step,
                                    const AlleleStep& allele_step,
                                    const EditAlignmentScorer& read_scorer,
                                    double& nat_adjust) const;

    /// The parts of scoring a read against an allele that depend on the read alone.
    struct ReadScratch {
        /// The score of the read's own edits on the node of each read step.
        vector<int32_t> own;
        /// The real-valued nats of each read step that its score in `own` cannot carry.
        vector<double> own_nats;
        /// The read's visits as step keys, sorted for binary search.
        vector<int64_t> keys;
    };

    /// Fill `scratch` for one read.
    void prepare_read_scratch(const Alignment& aln, const vector<ReadStep>& read_steps,
                              const EditAlignmentScorer& read_scorer,
                              ReadScratch& scratch) const;

    /// Return the allele's visits as step keys, sorted for binary search.
    static vector<int64_t> sorted_allele_keys(const vector<AlleleStep>& allele_steps);

    /// The positions at which each step key occurs in an allele, in ascending order.
    using AlleleStepPositions = unordered_map<int64_t, vector<uint32_t>>;

    /// Return where each step key occurs in the allele, for greedy pairing.
    static AlleleStepPositions index_allele_steps(const vector<AlleleStep>& allele_steps);

    /**
     * Score one read against one allele by greedy pairing.
     *
     * @pre The allele has at least one step.
     * @param nat_adjust Gets real-valued corrections added to it (see `score_shared_node`).
     */
    int32_t score_by_greedy_pairing(const Alignment& aln,
                                    const vector<ReadStep>& read_steps,
                                    const vector<AlleleStep>& allele_steps,
                                    const AlleleStepPositions& allele_positions,
                                    const EditAlignmentScorer& read_scorer,
                                    double& nat_adjust) const;

    /**
     * Score one read against one allele by optimal pairing.
     *
     * @pre The allele has at least one step.
     * @param nat_adjust Gets real-valued corrections added to it (see `score_shared_node`).
     */
    int32_t score_by_optimal_pairing(const Alignment& aln, const vector<ReadStep>& read_steps,
                                     const vector<AlleleStep>& allele_steps,
                                     const ReadScratch& scratch,
                                     const vector<int64_t>& allele_keys,
                                     const EditAlignmentScorer& read_scorer,
                                     double& nat_adjust) const;

    /// Return the mismapping probability e_r of a read, before clamping.
    double mismap_probability(const Alignment& aln) const;

    /// Fill `row` with the ln-likelihood of a read under each allele, or -inf for an allele that
    /// cannot place it.
    void score_row(const Alignment& aln, const vector<ReadStep>& read_steps,
                   const EditAlignmentScorer& read_scorer, const SiteAlleles& alleles,
                   ReadScratch& scratch, vector<double>& row) const;

    /// The statistics of the reads that begin in a rate window.
    struct WindowReadStats {
        /// The number of reads that begin on a reference node in the window, per base of
        /// reference in the window.
        ///
        /// Each read counts as 1 - e_r when `depth_effective_reads` is set. The rate is per base,
        /// not per haplotype. It is 0 if no read begins in the window, which turns the depth term
        /// off.
        double start_rate = 0.0;
        /// The mean length of the reads that begin in the window.
        double mean_read_length = 0.0;
    };

    /**
     * Return the statistics of the reads in the rate window of a site.
     *
     * The window is placed on reference coordinates, so that renumbering the graph's nodes does
     * not change it. The reference path is cut into buckets of RATE_BUCKET bp. A site falls in the
     * bucket holding the reference position of its start bound, or of its end bound, or failing
     * both of the nearest ancestor with a bound on a reference path. Its window is that bucket and
     * one bucket on each side. Only reference nodes count, in the numerator and the denominator
     * alike. An off-reference node has no position to place it in a window, and counting its
     * length without its reads, or its reads without its length, would bias the rate.
     *
     * Without reference positions (see `set_rate_reference`), or for a site with no reference
     * bound in its ancestry, the window is the block of RATE_ID_WINDOW consecutive node IDs
     * holding the site's lowest node ID. That does depend on the numbering, and such sites are
     * counted in `rate_id_fallbacks`.
     *
     * The depth term's lambda = rate * (L + R - 1) counts the reads whose start position places
     * them over an interval of length L, so the rate must count read starts, and not reads that
     * overlap the window. R should be the mean length of all reads, and not of the reads at a
     * site, because a long read overlaps more sites. Counting only the reads that begin in the
     * window gives both.
     *
     * The window is fixed, rather than taken from the read source's fetch window, so that DR and
     * the depth term do not depend on how the reads were supplied. It is computed once and shared
     * by the sites in it.
     */
    WindowReadStats local_read_stats(const MaterializedSite& site,
                                     const vector<nid_t>& fallback_nodes,
                                     const vector<pair<nid_t, nid_t>>& site_ranges) const;

    /// Return the statistics of the reads in the node-ID fallback window of a site (see
    /// `local_read_stats`).
    WindowReadStats id_window_read_stats(const vector<pair<nid_t, nid_t>>& site_ranges) const;

    /// The counts of the reads that begin on a set of nodes, and the nodes' total length.
    ///
    /// Windows add up the counts of their buckets before dividing.
    struct StartCounts {
        /// The number of reads, each counted as 1 - e_r when `depth_effective_reads` is set, and
        /// as 1 otherwise.
        double reads = 0.0;
        /// The total length of the reads.
        double length_total = 0.0;
        /// The number of reads whose length is in `length_total`.
        size_t length_count = 0;
        /// The total length of the nodes.
        size_t bp = 0;
        /// Add another set of counts to these.
        void add(const StartCounts& other) {
            reads += other.reads;
            length_total += other.length_total;
            length_count += other.length_count;
            bp += other.bp;
        }
        /// Return the statistics of a window with these counts.
        WindowReadStats stats() const;
    };

    /// Count the reads that begin on the given nodes, whose total length is `bp`.
    StartCounts count_starts(const vector<nid_t>& nodes, size_t bp) const;

    /// Return the counts for one reference bucket, computing them only once.
    StartCounts bucket_counts(size_t path_index, int64_t bucket) const;

    /**
     * Find the reference path index and position that place `site` in a rate window.
     *
     * The position is that of one of the site's bounds on a reference path, or else of the first
     * of `fallback_nodes` on one.
     *
     * @return False if there is no such position.
     */
    bool rate_position(const MaterializedSite& site, const vector<nid_t>& fallback_nodes,
                       size_t& path_index, int64_t& position) const;

    const PathHandleGraph& graph;

private:
    const SiteReadSource& read_source;
    /// The statistics of each node-ID fallback window, by window index.
    mutable unordered_map<size_t, WindowReadStats> window_rate;
    /// The counts of each reference bucket, by reference path index and bucket.
    ///
    /// Each bucket is counted once, although three windows use it.
    mutable unordered_map<pair<size_t, int64_t>, StartCounts> ref_bucket_counts;
    /// The statistics of the rate window centred on each reference bucket, by reference path
    /// index and bucket.
    mutable unordered_map<pair<size_t, int64_t>, WindowReadStats> ref_window_rate;
    /// Guards `window_rate`, `ref_bucket_counts` and `ref_window_rate`.
    mutable std::mutex window_bp_mutex;
    const PathPositionHandleGraph* rate_graph = nullptr;
    vector<path_handle_t> rate_paths;
    mutable std::atomic<size_t> id_fallbacks{0};
    const EditAlignmentScorer& qual_scorer;
    const EditAlignmentScorer& plain_scorer;
    Params params;
};

}

#endif
