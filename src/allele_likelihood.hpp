#ifndef VG_ALLELE_LIKELIHOOD_HPP_INCLUDED
#define VG_ALLELE_LIKELIHOOD_HPP_INCLUDED

/** \file allele_likelihood.hpp
 *
 * The read-by-allele likelihood matrix of one site, and the code that builds it
 * from the reads' existing alignments to the graph.
 *
 * Each entry says how well one read fits one candidate allele. The read-likelihood
 * genotyper combines these entries into a likelihood for each genotype. The model
 * is described in doc/read-likelihood-genotyping.md.
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

#include "alignment_scorer.hpp"
#include "handle.hpp"
#include "site_read_source.hpp"
#include "snarls.hpp"

#include "anchor.hpp"
#include "read_phasing.hpp"

namespace vg {

using namespace std;

/**
 * Per-read relative likelihoods over the alleles of one site.
 *
 * Each row belongs to a read and is divided by its own maximum, so every entry is
 * in [0, 1] and the read's best allele scores 1. The alignment score gives
 * ln P(read | allele) only up to a constant that differs between reads; dividing
 * by the row maximum removes that constant. It also puts every row on the scale of
 * the mismapping term, whose background likelihood is then 1, and lets rows from
 * scorers with different log bases share one matrix.
 *
 * The values can be used to compare genotypes (to choose one, to compute GQ, or to
 * form a posterior), but they are not calibrated absolute probabilities.
 */
class AlleleReadLikelihoods {
public:

    AlleleReadLikelihoods() = default;

    size_t num_reads() const { return n_reads; }
    size_t num_alleles() const { return n_alleles; }

    /// Likelihood of read r under allele a, relative to read r's best allele at this
    /// site, in [0, 1]. It is 0 where allele a fits the read so much worse than the best
    /// that the ratio underflows, or where allele a has no node visits at all, as the
    /// empty traversal of a star allele has, and cannot place the read.
    double rel(size_t r, size_t a) const;

    /// The read's mismapping probability e_r, derived from its MAPQ and clamped to lie
    /// strictly inside (0, 1).
    double mismap_prob(size_t r) const;

    /// ln of the row's divisor, the read's best log-likelihood score at this site, in
    /// nats. The genotype likelihood does not use it; it feeds the BL output field.
    double best_ln_likelihood(size_t r) const;

    /// How many reads were dropped because they placed on no allele at all.
    size_t num_unplaceable() const { return unplaceable; }

    /// Turn on length-weighted mixture weights in place of a flat 1/|G|, and set R.
    ///
    /// Each haplotype of a genotype is weighted by the number of start positions from
    /// which a read of length R overlaps a stretch of length X:
    ///
    ///     w_h = (X_h + R - 1) / sum_{h' in G} (X_h' + R - 1)
    ///
    /// X_h is the unique length U_h given by `set_unique_lengths`, which is how the
    /// calculator builds every matrix. Without unique lengths it is the allele's full
    /// length L_h, from `allele_lengths`, which is indexed by allele.
    ///
    /// `mean_read_length` is R, which should be the mean length of reads in the
    /// neighbourhood rather than of the reads at this site, since a long read overlaps
    /// more sites and so is over-represented at each. `compute` replaces it with such a
    /// neighbourhood mean after `build`. An empty vector or a zero R gives the flat 1/|G|.
    void set_length_weights(vector<size_t> allele_lengths, double mean_read_length) {
        this->allele_lengths = std::move(allele_lengths);
        this->mean_read_length = mean_read_length;
    }

    /// Weight the mixture by the sequence unique to each allele (see
    /// `set_length_weights`).
    ///
    /// Reads that lie in sequence both alleles share fit both equally and cannot
    /// change a genotype comparison, so the weight uses U_h: the total length of the
    /// nodes that haplotype h's allele visits and the genotype's other allele does not,
    /// in either orientation. `unique_lengths[a][b]` is that length for allele a against
    /// allele b.
    void set_unique_lengths(vector<vector<size_t>> unique_lengths) {
        this->unique_lengths = std::move(unique_lengths);
    }

    /// Set up the depth term, w * ln Poisson(N ; lambda_G), which asks whether genotype
    /// G predicts the number of reads seen, N. Here
    ///
    ///     lambda_G = rate * sum_{h in G} (T_h + R - 1)
    ///
    /// `traversal_lengths` gives T_h for each allele: its length without the site's two
    /// boundary nodes, since a read inside a boundary node is not a row of this matrix.
    /// `rate` is read starts per base per haplotype near the site, and `read_length` is
    /// R. With `effective_count`, N counts each read as 1 - e_r, and `rate` must count
    /// reads the same way. The term is off unless `weight` is positive.
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

    /// Expected number of reads at this site under the genotype, lambda_G.
    double expected_reads(const vector<int>& genotype) const;

    /// lambda_G = rate * sum_{h in G} (T_h + R - 1) from the depth context's parts, for a caller
    /// that keeps them after the matrix is gone. `expected_reads` is this with the matrix's own.
    /// An allele outside `traversal_lengths`, such as a negative marker, has T_h = 0.
    static double expected_reads_from(const vector<size_t>& traversal_lengths, double rate,
                                      double read_length, const vector<int>& genotype);

    /// True once set_depth_context has supplied traversal lengths.
    bool has_depth_context() const { return !traversal_lengths.empty(); }

    /// The read-start rate per base per haplotype that set_depth_context was given.
    double depth_rate_per_haplotype() const { return depth_rate; }

    /// The mean read length R that set_depth_context was given.
    double depth_read_length_used() const { return depth_read_length; }

    /// The read count N that the depth term compares with lambda_G: sum_r (1 - e_r)
    /// with `effective_count`, and the number of rows otherwise. Since N need not be a
    /// whole number, the depth term treats it with a continuous analogue of the Poisson
    /// distribution over read counts n >= 0, using lgamma in place of ln n!. That
    /// density is used only up to its normalising constant, which depends on lambda_G
    /// and is left out.
    double observed_reads() const;

    /// This allele's length without the site's boundary nodes, T_h. 0 if the depth
    /// context was not set or the index is out of range.
    size_t traversal_length(size_t allele) const {
        return allele < traversal_lengths.size() ? traversal_lengths[allele] : 0;
    }

    /// Observed over expected read count, N / lambda_G, for the given genotype, or -1
    /// where lambda_G is 0, as it is when no read begins in the rate window. It is
    /// written as the DR output field whether or not the depth term is on.
    double depth_ratio(const vector<int>& genotype) const;

    /// The mean read length R used by the mixture weights and copied into the evidence
    /// kept for read phasing and anchors.
    double mean_read_length_estimate() const { return mean_read_length; }

    /// Set R for the mixture weights, and for the evidence copied from this matrix,
    /// without turning on length weighting. The depth term takes its own R from
    /// set_depth_context. `compute` passes both the rate window's mean read length.
    void set_mean_read_length(double mean_read_length) {
        this->mean_read_length = mean_read_length;
    }

    /// True when set_length_weights supplied usable data.
    bool uses_length_weights() const {
        return !allele_lengths.empty() && mean_read_length > 0.0;
    }

    /**
     * ln P(reads | G), where G is a multiset of allele indices of size ploidy:
     *
     *   ln P(reads | G) =  sum_r ln [ (1 - e_r) * sum_{h in G} w_h * rel(r,h) + e_r ]
     *                    + w_d * ln Poisson( N ; lambda_G )
     *
     * r runs over the site's reads and h over the haplotypes of G, so a homozygote
     * counts its allele twice. e_r is the read's mismapping probability and w_h the
     * mixture weight. The second term is the depth term (see set_depth_context).
     *
     * Since rel(r,h) lies in [0, 1], each read's term lies between ln(e_r) and 0, so
     * the floor on e_r limits how much one read can count against a genotype.
     * rel(r,h) = 0 is the strongest evidence against allele h that one read can give.
     * Reads are treated as independent, so GL and GQ grow over-confident with depth.
     */
    double genotype_likelihood(const vector<int>& genotype) const;

    /**
     * The largest difference the read term could give between the called genotype and
     * the runner-up at this site, the denominator of GQN.
     *
     * It is the difference that would result if each haplotype of `called`
     * contributed its mixture weight's share of the reads, each read fitted its own
     * haplotype's allele with rel 1 and every other allele with 0, and every read had
     * e_r at the mismap floor. Dividing the observed difference by this one gives a
     * value in [0, 1] that does not depend on depth or on ploidy.
     *
     * The depth term is left out, since a perfect set of reads does not maximise it.
     * Returns 0 when there are no reads or the two genotypes are equal, so callers
     * must check before dividing.
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

    /// The site's anchor evidence, filled only when anchors are being written. The
    /// caller moves it into the CallInfo it keeps for the site.
    unique_ptr<AnchorSiteEvidence> anchor_evidence;

    /// The site's read-phasing evidence, filled only when read phasing is on and
    /// `anchor_evidence` is not filled, since that holds the same rows.
    unique_ptr<PhaseReadEvidence> phase_evidence;

    /// Populate the matrix. Only for AlleleReadLikelihoodsBuilder.
    void set_contents(size_t n_reads, size_t n_alleles, vector<double>&& matrix,
                      vector<double>&& mismap, vector<double>&& best_ln,
                      vector<string>&& names, size_t unplaceable);

    /// The --mismap-min floor on e_r, which achievable_gap uses for its ideal reads.
    void set_mismap_floor(double floor) { this->mismap_floor = floor; }

private:
    /// The mixture weight of each haplotype of the genotype (see `set_length_weights`):
    /// flat 1/|G| unless lengths were supplied. Shared by genotype_likelihood and
    /// achievable_gap, which must use the same weights.
    vector<double> mixture_weights(const vector<int>& genotype) const;

    /// Row major, n_reads * n_alleles, every entry in [0,1], row max exactly 1.
    vector<double> matrix;
    vector<double> read_mismap_prob;
    /// --mismap-min, used by achievable_gap only.
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
 * Reads are added one at a time, each with its raw ln-likelihood against every
 * allele. A read whose entries are all -inf placed on no allele and has no row
 * maximum to divide by, so it is dropped and counted.
 */
class AlleleReadLikelihoodsBuilder {
public:
    /// Mismapping probabilities are clamped into [min_mismap, max_mismap].
    ///
    /// The upper clamp matters because many mappers give MAPQ 0 to a read with several
    /// equally good placements. Its unclamped e_r of 1 would make the read's term 0
    /// under every genotype, so the read would count for nothing.
    AlleleReadLikelihoodsBuilder(size_t num_alleles, double min_mismap = 0.02,
                                 double max_mismap = 0.95);

    /// Add a read. raw_ln_likelihood must have one entry per allele and may
    /// contain -inf for alleles that cannot place the read.
    /// read_length feeds the mean R used by the length-weighted mixture. Zero
    /// means "unknown"; if every read is unknown the mixture stays flat.
    /// `start` is where the read's alignment begins, which with its name identifies the
    /// alignment, since paired mates share a name; build() orders the rows by the two.
    /// Returns false when the read placed on no allele at all and was dropped, so a caller
    /// accumulating anything alongside the rows can stay in step with them.
    bool add_read(const vector<double>& raw_ln_likelihood, double mismap_prob,
                  const string& name = "", size_t read_length = 0,
                  const Position* start = nullptr);



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
    ///
    /// The rows are put in a canonical order, by read name, then by where the alignment begins,
    /// then by the row's own values, so that the matrix, and every sum over its reads, does not
    /// depend on the order the reads were added in. That order is the read source's, which can
    /// change with the fetch window (--read-window).
    AlleleReadLikelihoods build();

    /// After build(): for each row of the matrix, the index among the reads add_read kept, in
    /// the order it kept them. A caller that kept something per read alongside the rows
    /// reorders it by this.
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
    /// Per kept read, where its alignment begins: node, offset and strand.
    vector<tuple<nid_t, int64_t, bool>> starts;
    vector<size_t> order;
    size_t unplaceable = 0;
};

/**
 * Tuning for GraphAlignedAlleleLikelihoodCalculator.
 *
 * At namespace scope rather than nested in the calculator so that its default
 * member initializers can be used in a defaulted argument.
 */
struct AlleleLikelihoodParams {
    /// Nats added to a read's log-likelihood for each gap in which the read has bases
    /// the allele lacks (--insertion-nats). A positive value makes extra read bases
    /// count against an allele less than missing ones. It is applied after the integer
    /// alignment score is converted to nats, so it can take fractional values.
    double insertion_gap_nats = 0.0;

    /// Choose each read's pairing with an allele by optimal pairing rather than greedy
    /// pairing (--optimal-pairing); see GraphAlignedAlleleLikelihoodCalculator.
    bool optimal_pairing = false;

    /// The floor on the mismapping probability e_r (--mismap-min).
    ///
    /// A read that fits allele A perfectly and allele B not at all lowers B's
    /// likelihood by at most -ln(e_r), so the floor limits how strongly one read can
    /// count against an allele. MAPQ measures whether the read is at the right locus,
    /// not whether its path through this site is right, so the floor also stands for
    /// a well-mapped read that is misaligned locally.
    double min_mismap_prob = 0.02;

    /// The ceiling on e_r (--mismap-max). It applies to reads with MAPQ 0 or close to
    /// it, and decides how much such a read still counts. It must stay below 1, since
    /// at e_r = 1 the read's term is 0 under every genotype.
    double max_mismap_prob = 0.95;

    /// Use the mismapping term. When false, every e_r is set to the floor, the closest
    /// the model can come to trusting every read fully while keeping the log finite.
    bool use_mismap_term = true;

    /// Weight each haplotype of a genotype by its share of the reads that can tell the
    /// genotype's alleles apart, from the sequence unique to its allele, rather than by a
    /// flat 1/|G| (--flat-mixture turns it off). See
    /// AlleleReadLikelihoods::set_length_weights.
    bool length_weighted_mixture = true;

    /// Weight of the depth term, ln P(N | G) (--depth-term). Zero turns the term off;
    /// DR is computed either way.
    double depth_weight = 0.1;

    /// Count each read toward depth as 1 - e_r, the probability that it came from this
    /// locus, rather than as 1. The local rate is counted the same way.
    bool depth_effective_reads = true;

    /// Ploidy used for the depth rate when `compute` is called with a region ploidy of 0
    /// or less; at 0 or less here too, the rate window is not measured.
    int depth_ploidy = 2;

    /// Collect per-read anchor evidence while the reads are in memory (--anchors-out).
    /// Each read's position is resolved here, while its alignment is available.
    bool collect_anchors = false;
    /// Where to count the position resolutions this calculator performs. Not owned;
    /// see AnchorParams::counters. Null means do not count.
    AnchorCounters* anchor_counters = nullptr;

    /// Keep each read's row of relative likelihoods for read-backed phasing
    /// (--read-phasing). Phasing needs only the rows, not the positions that
    /// `collect_anchors` also resolves.
    bool collect_read_phasing = false;
};

/**
 * Produces the reads x alleles matrix for a site.
 */
class AlleleLikelihoodCalculator {
public:
    virtual ~AlleleLikelihoodCalculator() = default;

    /// Build the matrix for one site. traversals are the candidate alleles, in
    /// the order the caller will genotype them. `region_ploidy` is the number of the
    /// sample's haplotypes in the region around the site, which the depth term divides the
    /// local read rate by to get a per-haplotype rate. It is the site's own ploidy at a
    /// top-level site, and more than it at a nested site that only some of its parent's
    /// alleles cross, since the reads counted near the site come from all of them.
    /// Nothing in the matrix depends on the ploidy the site is then genotyped at.
    virtual AlleleReadLikelihoods compute(const Snarl& snarl,
                                          const vector<SnarlTraversal>& traversals,
                                          int region_ploidy) = 0;
};

/**
 * Scores each read against each allele from the read's existing alignment to the graph, rather
 * than by aligning the read to each allele again.
 *
 * ## Pairings
 *
 * A read's placement and an allele are both sequences of node visits, a visit being a node in one
 * orientation. A read is scored against an allele through a *pairing* of the two sequences, in
 * order: each read visit inside the site is either paired with an allele visit or left unpaired.
 * The pairing's score is the sum of these parts:
 *
 *   - a read visit paired with the same allele visit scores the read's own edits in that node;
 *   - a read visit that the allele never makes, paired with an allele visit that the read never
 *     makes, is a substitution: the read's bases in its node are compared one by one with the
 *     allele node's sequence, from the first base of each, plus a gap for the difference in
 *     length;
 *   - a run of read visits left unpaired is an insertion, and a run of allele visits skipped
 *     between two pairs is a deletion, each scored as one gap.
 *
 * Allele visits before the read's first same-visit pair, or after its last pair, lie outside the
 * read and score nothing. Before that first same-visit pair, each unpaired read visit is a gap of
 * its own, since there is no pair yet for an insertion to extend from.
 *
 * Every read base inside the site is scored under every allele, so all alleles are scored over the
 * same read bases, the read's *scoring window*, and differ only in how well they explain them.
 *
 * ## Two ways to choose the pairing
 *
 *   - *Greedy pairing* (`score_by_greedy_pairing`, the default) makes one pass along the read's
 *     visits. Each is paired with the same visit's next occurrence in the allele after the last
 *     pair, and a pair is never revised. The pairing it finds is scored by the rules above, but it
 *     can pair different visits as a substitution even where one of them is shared.
 *   - *Optimal pairing* (`score_by_optimal_pairing`, --optimal-pairing) finds the highest-scoring pairing
 *     that the rules above allow, by dynamic programming over the two sequences. A run of unpaired
 *     read visits is one gap. At large sites the search is restricted to a band.
 *
 * Inside a node that the read and the allele share, both use the mapper's edits and neither aligns
 * bases again.
 */
class GraphAlignedAlleleLikelihoodCalculator : public AlleleLikelihoodCalculator {
public:

    /// Defined at namespace scope; aliased here for readability.
    using Params = AlleleLikelihoodParams;

    /**
     * `qual_scorer` charges each mismatch at its base quality, and scores reads that
     * have base qualities. `plain_scorer` scores reads without them, such as reads
     * from a GAF file with no quality column. The two have different log bases, which
     * does not matter because each row is divided by its own maximum and a row is
     * scored by one scorer only.
     */
    GraphAlignedAlleleLikelihoodCalculator(const PathHandleGraph& graph,
                                           SnarlManager& snarl_manager,
                                           const SiteReadSource& read_source,
                                           const EditAlignmentScorer& qual_scorer,
                                           const EditAlignmentScorer& plain_scorer,
                                           const Params& params = Params());

    AlleleReadLikelihoods compute(const Snarl& snarl,
                                  const vector<SnarlTraversal>& traversals,
                                  int region_ploidy) override;

    /// Place rate windows on these reference paths of `position_graph`, which must be the
    /// graph the calculator was built on, or a view of it with the same nodes. Until this is
    /// called, windows are blocks of node IDs (see local_read_stats).
    void set_rate_reference(const PathPositionHandleGraph* position_graph,
                            const vector<path_handle_t>& reference_paths);

    /// Sites whose rate window fell back to a block of node IDs.
    size_t rate_id_fallbacks() const {
        return id_fallbacks.load();
    }

    /// The width, in bp of reference, of a rate-window bucket. A window is three buckets.
    static const int64_t RATE_BUCKET = 16384;

    /// The width, in node IDs, of the fallback rate window.
    static const nid_t RATE_ID_WINDOW = 4096;

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

    /// Extract the read's visits inside the site, in read order. Returns false if the
    /// read cannot tell the alleles apart because it lies within one boundary node,
    /// which every allele shares. A read that fails to place on some allele is kept,
    /// since that is evidence against the allele.
    bool get_read_steps(const SiteRead& read, const unordered_set<nid_t>& site_nodes,
                        const unordered_set<nid_t>& boundary_nodes,
                        vector<ReadStep>& steps_out) const;

    /// True if the read traverses the site against the direction the alleles read
    /// in, so it must be reverse-complemented before being compared to them.
    /// Decided by a vote over the nodes that every allele visits in the same
    /// orientation; a tie leaves the read as it is.
    bool read_is_reverse_of_alleles(const vector<ReadStep>& read_steps,
                                   const unordered_map<nid_t, bool>& allele_orientations) const;

    /// Score the read's own edits on a node the read and the allele share.
    /// `nat_adjust` accumulates real-valued corrections that cannot be expressed in the
    /// integer score; the caller adds it after converting the score to nats.
    int32_t score_shared_node(const Alignment& aln, const ReadStep& step,
                              const EditAlignmentScorer& read_scorer, double& nat_adjust) const;

    /// Score `length` of the read's own bases, from `read_offset`, against as many
    /// allele bases, from `allele_offset`, base by base, charging each mismatch at
    /// its own base quality.
    int32_t score_substitution(const Alignment& aln, size_t read_offset, size_t length,
                               const string& allele_bases, size_t allele_offset,
                               const EditAlignmentScorer& read_scorer) const;

    /// The parts of scoring a read against an allele that depend on the read alone,
    /// computed once per read rather than once per allele.
    struct ReadScratch {
        /// Per read step: the score of the read's own edits on that node, and the
        /// real-valued nats the integer score cannot carry.
        vector<int32_t> own;
        vector<double> own_nats;
        /// The read's (node, orientation) keys, sorted for binary search.
        vector<int64_t> keys;
    };

    /// Fill `scratch` for one read. Once per read, not once per (read, allele).
    void prepare_read_scratch(const Alignment& aln, const vector<ReadStep>& read_steps,
                              const EditAlignmentScorer& read_scorer,
                              ReadScratch& scratch) const;

    /// The allele's (node, orientation) keys, sorted for binary search. Computed once
    /// per allele per site.
    static vector<int64_t> sorted_allele_keys(const vector<AlleleStep>& allele_steps);

    /// The positions at which each (node, orientation) occurs in the allele, in
    /// ascending order. Used by greedy pairing and computed once per allele.
    using AlleleStepPositions = unordered_map<int64_t, vector<uint32_t>>;
    static AlleleStepPositions index_allele_steps(const vector<AlleleStep>& allele_steps);

    /// Score one read against one allele by greedy pairing.
    int32_t score_by_greedy_pairing(const Alignment& aln,
                                    const vector<ReadStep>& read_steps,
                                    const vector<AlleleStep>& allele_steps,
                                    const AlleleStepPositions& allele_positions,
                                    const EditAlignmentScorer& read_scorer,
                                    bool& placed_out, double& nat_adjust) const;

    /// Score one read against one allele by optimal pairing.
    int32_t score_by_optimal_pairing(const Alignment& aln, const vector<ReadStep>& read_steps,
                                     const vector<AlleleStep>& allele_steps,
                                     const ReadScratch& scratch,
                                     const vector<int64_t>& allele_keys,
                                     const EditAlignmentScorer& read_scorer,
                                     bool& placed_out, double& nat_adjust) const;

    /// Read statistics for a site's rate window: the number of reads whose alignment begins
    /// on a reference node in the window, per base of reference in the window, and the mean
    /// length of those reads. Computed once per window and shared by the sites in it.
    ///
    /// The window is placed on reference coordinates, so that renumbering the graph's nodes
    /// does not change it. The reference path is cut into buckets of RATE_BUCKET bp; a site
    /// falls in the bucket holding the reference position of its start boundary, or of its
    /// end boundary, or failing both of the nearest ancestor snarl with a boundary on a
    /// reference path. Its window is that bucket and one bucket on each side. Only reference
    /// nodes count, in the numerator and the denominator alike: an off-reference node has no
    /// position to place it in a window, and counting its length without its reads, or its
    /// reads without its length, would bias the rate.
    ///
    /// Without reference positions (see set_rate_reference), or for a site with no reference
    /// boundary in its ancestry, the window falls back to the block of RATE_ID_WINDOW
    /// consecutive node IDs holding the site's lowest node ID, which does depend on the
    /// numbering; such sites are counted in rate_id_fallbacks.
    ///
    /// The depth term's lambda = rate * (L + R - 1) counts the reads whose start
    /// position places them over an interval of length L, so the rate must count read
    /// starts, not reads that overlap the window. R should be the mean length of all
    /// reads, not of the reads that reach a site, because a long read reaches more
    /// sites. Counting only the reads that begin in the window gives both. Where no
    /// read begins in the window, the caller falls back to the site's own reads for R.
    ///
    /// Each read is weighted by 1 - e_r when `depth_effective_reads` is set. The rate
    /// is per base, not per haplotype: each site divides it by the region's ploidy. It is 0
    /// if no read begins in the window, which turns that site's depth term off.
    struct WindowReadStats {
        double start_rate = 0.0;
        double mean_read_length = 0.0;
    };
    WindowReadStats local_read_stats(const Snarl& snarl,
                                     const vector<pair<nid_t, nid_t>>& site_ranges) const;

    /// The node-ID fallback of local_read_stats.
    WindowReadStats id_window_read_stats(const vector<pair<nid_t, nid_t>>& site_ranges) const;

    /// The reads that begin on a set of nodes, as counted for a rate window, and the nodes'
    /// total length. Windows add them up before dividing.
    struct StartCounts {
        double reads = 0.0;
        double length_total = 0.0;
        size_t length_count = 0;
        size_t bp = 0;
        void add(const StartCounts& other) {
            reads += other.reads;
            length_total += other.length_total;
            length_count += other.length_count;
            bp += other.bp;
        }
        WindowReadStats stats() const;
    };

    /// Count the reads that begin on the given nodes, whose total length is `bp`.
    StartCounts count_starts(const vector<nid_t>& nodes, size_t bp) const;

    /// The counts for one reference bucket, computed once and kept.
    StartCounts bucket_counts(size_t path_index, int64_t bucket) const;

    /// The reference path index and position that place `snarl` in a rate window, if any.
    bool rate_position(const Snarl& snarl, size_t& path_index, int64_t& position) const;

    const PathHandleGraph& graph;
    SnarlManager& snarl_manager;
    const SiteReadSource& read_source;
    mutable unordered_map<size_t, WindowReadStats> window_rate;
    /// Both keyed by reference path index and bucket: a bucket's own counts, and the rate
    /// window centred on it. Each bucket's reads are fetched once, although three windows
    /// use them.
    mutable unordered_map<pair<size_t, int64_t>, StartCounts> ref_bucket_counts;
    mutable unordered_map<pair<size_t, int64_t>, WindowReadStats> ref_window_rate;
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
