#ifndef VG_READ_LIKELIHOOD_CALLER_HPP_INCLUDED
#define VG_READ_LIKELIHOOD_CALLER_HPP_INCLUDED

/** \file read_likelihood_caller.hpp
 *
 * A SnarlCaller that genotypes from an explicit P(reads | genotype) model rather
 * than from aggregate read depth.
 *
 * The model, its parameters, and the meaning of every field this writes to the VCF are
 * specified in doc/read-likelihood-genotyping.md. The quality fields in particular are
 * easy to misread -- GQ here is scaled and is a ranking score, not a calibrated
 * posterior -- and that page is where the distinction is spelled out.
 */

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "allele_likelihood.hpp"
#include "snarl_caller.hpp"

namespace vg {

using namespace std;

/**
 * Genotypes a site by building a reads x alleles likelihood matrix and scoring
 * every allele combination as a proper likelihood.
 *
 * ## Why this subclasses SupportBasedSnarlCaller
 *
 * Not for the genotyping, which uses no support at all, but because FlowCaller,
 * NestedFlowCaller and LegacyCaller all dynamic_cast their caller to
 * SupportBasedSnarlCaller to reach get_support_finder(), and dereference the
 * result without a null check. That support finder supplies the node and edge
 * weights FlowTraversalFinder uses to *enumerate* alleles. So the pack file
 * stays required -- for enumeration, not for genotyping.
 *
 * ## Relationship to the VCF layer
 *
 * VCFOutputCaller::emit_variant deduplicates alleles by sequence string and
 * drops uncalled ones before calling update_vcf_info, so by then the traversal
 * indices no longer match the matrix. Rather than change that interface, this
 * caller keeps a copy of the traversals it scored in its CallInfo and maps the
 * deduplicated ones back by structural comparison. That is exact regardless of
 * how the allele strings were flattened, needs no change to graph_caller.cpp,
 * and degrades cleanly: an allele it cannot match (the star allele's empty
 * placeholder traversal, say) simply gets no GL entry.
 */
class ReadLikelihoodSnarlCaller : public SupportBasedSnarlCaller {
public:

    ReadLikelihoodSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                              TraversalSupportFinder& support_finder,
                              AlleleLikelihoodCalculator& likelihood_calculator);

    virtual ~ReadLikelihoodSnarlCaller();

    struct ReadLikelihoodCallInfo : public SnarlCaller::CallInfo {
        virtual ~ReadLikelihoodCallInfo() = default;

        /// Phred-scaled gap between the best and second-best genotype.
        double gq = 0;
        /// ln posterior of the called genotype under a uniform prior.
        double posterior = 0;
        /// Reads seen at the site, and how many were informative enough to keep.
        size_t n_reads = 0;
        size_t n_informative = 0;
        /// Reads dropped for placing on no allele at all.
        size_t n_unplaceable = 0;
        /// The ploidy this site was genotyped at.
        int ploidy = 2;

        /// ln P(reads | G) for every genotype scored, keyed by the sorted allele
        /// index multiset so the VCF layer can look up by remapped indices.
        map<vector<int>, double> genotype_lls;
        /// For a site genotyped at ploidy 1, the genotype it would take at ploidy 2, as traversal
        /// indices; empty unless that measurement is switched on.
        ///
        /// Nested descent gives a child ploidy 1 when one parent allele crosses the chain, and
        /// linkage can then move the parent so that both do. This says what the site's own reads
        /// would call there, which is the thing a re-genotype would have to produce -- computed on
        /// the matrix that is already built, so it costs a few passes over memory and no re-reading.
        vector<int> alt_ploidy_best;

        /// The same site derived at the other ploidy, from the same reads-by-alleles matrix.
        ///
        /// Nested calling needs both answers from one visit to the reads: a chain's ploidy comes from
        /// its parent's *settled* genotype, which is not known while the reads are resident, and going
        /// back to find out is what made post-linkage descent cost half again as much read I/O. This
        /// carries everything a record needs -- `genotype_lls`, `gq`, `gq_fraction`, `explained_share`,
        /// `depth_ratio`, `posterior` -- so the record can be built later at whichever ploidy the
        /// barrier settles on, with no re-reading and no re-scoring.
        ///
        /// Null unless the alternate ploidy was requested (see `set_want_alt_ploidy`). Owned, and it is the only heavyweight member
        /// here: it shares the ploidy-independent halves (`scored_traversals`, `allele_support`) by
        /// copy rather than recomputing them.
        unique_ptr<ReadLikelihoodCallInfo> alt_ploidy_info;

        /// The traversals that were scored, in matrix column order. Kept so the
        /// deduplicated traversals handed to update_vcf_info can be mapped back.
        vector<SnarlTraversal> scored_traversals;

        /// Reads whose best-fitting allele is each scored allele, in matrix column
        /// order. A read fitting several alleles equally splits its vote between them
        /// rather than going to the lowest index: at multi-allelic sites many reads are
        /// genuinely undecided, and awarding those to one column would manufacture
        /// allele balance that the reads do not support.
        vector<double> allele_support;

        /// Mean over reads of the best raw score any allele gave them -- the row
        /// divisor, which the genotype likelihood divides out and never uses again. It
        /// measures whether reads fit *anything* here, where GQ measures only how far
        /// apart the top two genotypes are, so the two are close to independent.
        double mean_best_ln = 0;

        /// Fraction of reads whose best-fitting allele is one of the *called* alleles,
        /// derived from allele_support. 1.0 when the called genotype accounts for every
        /// read at the site.
        double explained_share = 1.0;

        /// Observed reads over the number the called genotype predicts, from the local
        /// rate and the called alleles' geometry. 1.0 means the read count is exactly
        /// what the call implies; 7.0 is a collapsed repeat. Negative means unavailable.
        ///
        /// Emitted whether or not `--depth-term` is armed, so it can be measured as a
        /// ranking signal before the model is allowed to act on it -- the same order the
        /// explained-share discount was established in.
        double depth_ratio = -1.0;

        /// GQ before the explained-share discount, so the discount stays auditable and a
        /// consumer that wants the raw likelihood-ratio quality can still have it.
        double gq_undiscounted = 0;

        /// The likelihood-ratio gap as a fraction of the gap this site could have
        /// produced -- see AlleleReadLikelihoods::achievable_gap. In [0,1], and unlike GQ
        /// it means the same thing at any depth and any ploidy, which is what makes a
        /// single threshold usable across a 5x diploid and a 15x haploid contig.
        /// Negative when there was no gap to normalise (no reads, or a site offering one
        /// genotype), which is not the same as 0 and must not be filtered as though it
        /// were.
        double gq_fraction = -1.0;

    };

    /// Also score a site genotyped at ploidy 1 at ploidy 2, and report what it would call there in
    /// `ReadLikelihoodCallInfo::alt_ploidy_best`. Off by default.
    ///
    /// Only meaningful for nested sites, whose ploidy comes from a parent genotype that linkage can
    /// afterwards invalidate. Costs a second pass over the matrix that is already built -- the reads
    /// are neither re-fetched nor re-scored, since `rel(r,a)` does not depend on ploidy.

    /// Ask for the other ploidy's answer on the next `genotype` call, on this thread only.
    ///
    /// Only a nested chain can have its ploidy revised: it comes from a parent genotype the linkage
    /// layer may still move. A top-level site takes its ploidy from the contig or a --ploidy-bed
    /// region and can never change, so computing an alternate for one is a second pass over the
    /// matrix and a retained second CallInfo for nothing.
    ///
    /// Getting this wrong is not cheap and does not show up where you would look for it. Generalising
    /// the old `ploidy == 1` gate to every site left all 165,408 of chr20's top-level snarls carrying
    /// an alternate, which put peak memory up by 1.4 GB while the retention counter -- which only
    /// counts nested chains -- still read 72 MB.
    ///
    /// Thread-local because the sweep is parallel over node-ID windows and this is set per call.
    static void set_want_alt_ploidy(bool on) { want_alt_ploidy = on; }

    virtual pair<vector<int>, unique_ptr<CallInfo>> genotype(const Snarl& snarl,
                                                             const vector<SnarlTraversal>& traversals,
                                                             int ref_trav_idx,
                                                             int ploidy,
                                                             const string& ref_path_name,
                                                             pair<size_t, size_t> ref_range);

    virtual void update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const unique_ptr<CallInfo>& call_info,
                                 const string& sample_name,
                                 vcflib::Variant& variant);

    virtual void update_vcf_header(string& header) const;

    /**
     * Never skip an allele.
     *
     * The inherited support-based version prunes any allele whose read support is
     * below a threshold, which exists to stop VCFTraversalFinder's brute-force
     * enumeration exploding at dense multi-allelic sites. Two reasons to drop it
     * here:
     *
     * - Scoring is DP-free, so the tractability argument for pruning is gone. We
     *   can afford to genotype against every enumerated traversal, including the
     *   long tail a support prefilter would have discarded.
     * - When there is no pack file the support finder reports zero everywhere, so
     *   the inherited version would prune *every* allele at *every* site.
     *
     * Only takes effect when support is unavailable, so behaviour with a pack file
     * is unchanged.
     */
    virtual function<bool(const SnarlTraversal&, int iteration)> get_skip_allele_fn() const;

    /// Tell the caller that its support finder reports nothing real, so anything
    /// support-derived must be skipped rather than believed.
    void set_support_available(bool available);

    /// Write the matrix for every site to this stream as TSV. Not owned.
    void set_likelihood_dump(ostream* dump_stream);

    /**
     * Scale GQ by the fraction of reads the called genotype explains.
     *
     * The genotype likelihood compares genotypes using only *which* alleles each read
     * fits, so a read that fits some allele the call does not contain is counted in every
     * genotype's likelihood and cancels out of the comparison. GQ therefore says nothing
     * about whether the called genotype accounts for the reads at all -- only about how
     * far ahead of its nearest rival it is. Sites where a third of the reads prefer an
     * uncalled allele get the same GQ as sites where none do.
     *
     * Discounting by that fraction improved the ranking of calls in every case measured:
     * two chromosomes, a 4- and a 34-haplotype graph, and both a small-variant and a
     * structural benchmark, at both moderate and high recall. The linear form was chosen
     * over stronger ones (share^2, share^4, and a phred cap on the unexplained fraction)
     * which score better on AUC but lose ground at high recall on some of those eight
     * combinations; linear was the only form that never made any of them worse.
     *
     * This does cost something real: a discounted GQ is no longer the phred-scaled
     * posterior odds of the top two genotypes, so it is a quality score rather than a
     * calibrated probability. GQI keeps the undiscounted value for anyone who needs it.
     * Off restores the previous behaviour exactly.
     */
    void set_share_discount(bool discount);

    /**
     * Scale GQ by how far the site's read count is from what the call predicts, at
     * records whose called alleles change length by at least `min_length` bp.
     *
     *     GQ' = GQ * exp(-exponent * |ln DR|)
     *
     * The sibling of set_share_discount, aimed at the model's other structural
     * blindness. That one asks whether the called genotype explains the reads that are
     * here; this asks whether the right *number* of reads is here at all, which
     * `P(reads | G)` cannot see because it is conditioned on the reads it was handed.
     * Both can only lower GQ, and GQI keeps the undiscounted value.
     *
     * **Off by default, and the size gate is why.** Ungated, this gains on structural
     * variants on all four datasets measured and *loses* on small variants on the two
     * 34-haplotype graphs -- the sign reversal that has blocked every previous attempt
     * to use depth here. The gate is not a fitted threshold: at a SNV, lambda's geometry
     * is dominated by the read length, so DR mostly reports local coverage scatter and
     * ranking on it adds noise to a GQ that already ranks well there; at a large event
     * the geometry is dominated by the allele, so DR reports whether the called sequence
     * is present at all. 50 bp is the boundary the two benchmarks already draw.
     *
     * Gated, it reaches 15 of 16 operating points improving or tying -- the standard the
     * explained-share discount met -- but structural-variant AUC still falls on one of
     * the eight cells, at every exponent tried, so it does not quite clear the bar that
     * put the share discount on by default. Hence a flag rather than a default.
     *
     * An `exponent` of 0 disables it and restores the previous behaviour exactly.
     */
    void set_depth_quality(double exponent, size_t min_length = 50);

    /**
     * Mark records whose GQN falls below `threshold` as `FILTER=lowconf`.
     *
     * This is the threshold GQN exists to make possible: one number that means the same
     * thing on a 5x diploid contig and a 15x haploid one. A raw GQ threshold cannot do
     * it. Measured over a coverage titration, requiring GQ >= 10 takes a 5x diploid
     * contig's F1 from 0.888 to 0.669, because at that depth half the true calls sit
     * below it; requiring GQN >= 0.05 costs 0.009 on the same arm and gains 0.017 on a
     * haploid one.
     *
     * At GQN >= 0.05 precision rises on **every** arm measured -- 0.9712 to 0.9759 and
     * 0.9828 to 0.9863 on diploid, 0.8863 to 0.9256 and 0.9389 to 0.9840 on haploid --
     * for one to two points of recall.
     *
     * **Off by default, and there is no good default to be had.** F1 weights precision
     * and recall equally, and under that weighting a gate helps the haploid arms and
     * hurts the diploid ones; no single threshold wins everywhere. What the threshold is
     * worth depends on which error a caller would rather make, which is not something
     * this code can know. 0 disables it.
     *
     * Marks rather than drops, for the reason given in update_vcf_info: records are
     * rewritten by the linkage layer after this runs, and a low-confidence site is
     * exactly the kind that layer exists to fix.
     */
    void set_min_confidence(double threshold);

protected:

    /// True if two traversals visit exactly the same nodes in the same
    /// orientations. Used to map deduplicated VCF alleles back to matrix columns.
    static bool traversals_equal(const SnarlTraversal& a, const SnarlTraversal& b);

    AlleleLikelihoodCalculator& likelihood_calculator;

    /// Optional TSV dump of every site's matrix, for development.
    ostream* dump_stream = nullptr;

    /// False when running without a pack file, so the support finder is a
    /// NullTraversalSupportFinder and reports zero for everything.
    bool support_available = true;

    /// See set_want_alt_ploidy.
    static thread_local bool want_alt_ploidy;


    /// Whether GQ is scaled by the explained-read fraction. See set_share_discount.
    bool share_discount = true;

    /// Exponent on |ln DR| in the depth-implausibility discount, and the minimum called
    /// allele length change that arms it. Zero exponent disables. See set_depth_quality.
    double depth_quality = 0.0;
    size_t depth_quality_min_length = 50;

    /// GQN below which a record is marked lowconf. Zero disables. See set_min_confidence.
    double min_confidence = 0.0;
};

}

#endif
