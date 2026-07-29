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
 */

#include <functional>
#include <iostream>
#include <limits>
#include <string>
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

    /**
     * ln P(reads | G), where G is a multiset of allele indices of size ploidy.
     *
     *   ln P(reads | G) = sum_r ln [ (1 - e_r) * sum_{h in G} (1/|G|) * rel(r,h)
     *                                + e_r * 1 ]
     *
     * The inner sum marginalises over which haplotype produced the read; for a
     * homozygous genotype its weights collapse to 1, so hom and het need no
     * special casing and allele balance comes out of the 1/|G| weighting rather
     * than a tuned bias parameter.
     *
     * The `+ e_r` term is "this read came from somewhere else entirely, where it
     * fits about as well as its best explanation here". Because rel is in [0,1],
     * the bracket lies in [e_r, 1], so its log is always finite and bounded below
     * by ln e_r: no single read can penalise a genotype without limit.
     *
     * e_r is genotype-independent, so it cannot favour one genotype over
     * another; a wrong e_r costs power rather than creating bias.
     *
     * Reads are assumed independent. Mates overlapping the same site are not,
     * and the product therefore accumulates confidence like R rather than
     * sqrt(R), so derived GQ will be over-confident and increasingly so with
     * depth. Treat GQ/GL as useful for ranking, not as calibrated probabilities.
     */
    double genotype_likelihood(const vector<int>& genotype) const;

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

private:
    /// Row major, n_reads * n_alleles, every entry in [0,1], row max exactly 1.
    vector<double> matrix;
    vector<double> read_mismap_prob;
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
    AlleleReadLikelihoodsBuilder(size_t num_alleles, double min_mismap = 1e-8,
                                 double max_mismap = 0.1);

    /// Add a read. raw_ln_likelihood must have one entry per allele and may
    /// contain -inf for alleles that cannot place the read.
    void add_read(const vector<double>& raw_ln_likelihood, double mismap_prob,
                  const string& name = "");

    size_t num_reads_added() const { return rows.size(); }
    size_t num_unplaceable() const { return unplaceable; }

    /// Normalise every row by its own maximum and produce the matrix.
    AlleleReadLikelihoods build();

private:
    size_t n_alleles;
    double min_mismap;
    double max_mismap;
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
    double min_mismap_prob = 1e-8;
    double max_mismap_prob = 0.1;
    /// Turn the mismapping term off entirely, so its contribution can be
    /// measured rather than assumed. With it off, e_r is pinned to the minimum,
    /// which is as close to "trust every read fully" as the model can get while
    /// keeping the log finite.
    bool use_mismap_term = true;
};

/**
 * Produces the reads x alleles matrix for a site.
 */
class AlleleLikelihoodCalculator {
public:
    virtual ~AlleleLikelihoodCalculator() = default;

    /// Build the matrix for one site. traversals are the candidate alleles, in
    /// the order the caller will genotype them.
    virtual AlleleReadLikelihoods compute(const Snarl& snarl,
                                          const vector<SnarlTraversal>& traversals) = 0;
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
                                  const vector<SnarlTraversal>& traversals);

    /// Reads seen at the last site, before dropping uninformative ones.
    size_t get_last_site_read_count() const { return last_site_reads; }

    /// Reads dropped at the last site for having no informative overlap.
    size_t get_last_site_uninformative_count() const { return last_site_uninformative; }

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

    const PathHandleGraph& graph;
    SnarlManager& snarl_manager;
    const SiteReadSource& read_source;
    const EditAlignmentScorer& qual_scorer;
    const EditAlignmentScorer& plain_scorer;
    Params params;
    size_t last_site_reads = 0;
    size_t last_site_uninformative = 0;
};

}

#endif
