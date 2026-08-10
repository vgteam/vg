#ifndef VG_READ_LIKELIHOOD_CALLER_HPP_INCLUDED
#define VG_READ_LIKELIHOOD_CALLER_HPP_INCLUDED

/** \file read_likelihood_caller.hpp
 *
 * A SnarlCaller that genotypes from an explicit P(reads | genotype) model rather
 * than from aggregate read depth.
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
    };

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

    /// Whether GQ is scaled by the explained-read fraction. See set_share_discount.
    bool share_discount = true;
};

}

#endif
