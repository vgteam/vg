#ifndef VG_READ_LIKELIHOOD_CALLER_HPP_INCLUDED
#define VG_READ_LIKELIHOOD_CALLER_HPP_INCLUDED

/** \file read_likelihood_caller.hpp
 *
 * A SnarlCaller that genotypes a site from the likelihood of its reads under each
 * genotype, P(reads | genotype), rather than from read depth.
 *
 * The model and the VCF fields this caller writes are described in
 * doc/read-likelihood-genotyping.md.
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
 * Genotypes a site by building its reads x alleles likelihood matrix and scoring
 * every genotype over the candidate alleles.
 *
 * ## Why this subclasses SupportBasedSnarlCaller
 *
 * The genotyping uses no read support. It derives from SupportBasedSnarlCaller
 * only because the graph callers reach their traversal finder's support through
 * that interface.
 *
 * ## Relationship to the VCF layer
 *
 * VCFOutputCaller::emit_variant merges alleles with the same sequence and drops
 * uncalled ones before calling update_vcf_info, so the traversal indices it passes
 * do not match the matrix. This caller keeps the traversals it scored in its
 * CallInfo and matches each remaining traversal to one of them node by node. A
 * record with an allele that matches none, such as the empty traversal of a star
 * allele, is written without GL, since GL needs an entry for every genotype.
 */
class ReadLikelihoodSnarlCaller : public SupportBasedSnarlCaller {
public:

    ReadLikelihoodSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                              TraversalSupportFinder& support_finder,
                              AlleleLikelihoodCalculator& likelihood_calculator);

    virtual ~ReadLikelihoodSnarlCaller();

    struct ReadLikelihoodCallInfo : public SnarlCaller::CallInfo {
        virtual ~ReadLikelihoodCallInfo() = default;

        /// Phred-scaled gap between the best and second-best genotype, times the discounts
        /// (see `discounted_gq`). Written as GQ.
        double gq = 0;
        /// ln posterior, under a uniform prior over genotypes, of the genotype with the
        /// highest likelihood. Written as GP.
        double posterior = 0;
        /// The number of reads in the site's matrix. Written as DP.
        size_t n_reads = 0;
        /// The ploidy this site was genotyped at.
        int ploidy = 2;

        /// ln P(reads | G) for every genotype scored, keyed by the sorted allele
        /// index multiset so the VCF layer can look up by remapped indices.
        map<vector<int>, double> genotype_lls;
        /// The best genotype at the other ploidy (2 for a site genotyped at ploidy 1, and 1
        /// otherwise), as traversal indices. Empty unless requested with set_want_alt_ploidy.
        vector<int> alt_ploidy_best;

        /// The whole call at the other ploidy, computed from the same matrix, so that the
        /// site's record can be built at whichever ploidy the barrier settles on. Null unless
        /// requested with set_want_alt_ploidy. Its ploidy-independent fields, such as
        /// `scored_traversals` and `allele_support`, are copies of this one's.
        unique_ptr<ReadLikelihoodCallInfo> alt_ploidy_info;

        /// The traversals that were scored, in matrix column order. Kept so the
        /// deduplicated traversals handed to update_vcf_info can be mapped back.
        vector<SnarlTraversal> scored_traversals;

        /// Per-read anchor evidence, with --anchors-out; null otherwise. It depends on the
        /// matrix, not on the ploidy, so when the barrier replaces this CallInfo with
        /// `alt_ploidy_info` it must move it across.
        unique_ptr<AnchorSiteEvidence> anchor_evidence;
        /// Per-read phasing evidence, with read phasing and no anchors; null otherwise. Moved
        /// across like `anchor_evidence`.
        unique_ptr<PhaseReadEvidence> phase_evidence;

        /// `genotype_lls` as the sweep computed them, before re-genotyping corrected them.
        /// Saved at the first correction, and each later round of re-genotyping corrects
        /// these rather than the previous round's values. Null until re-genotyping changes
        /// the site.
        unique_ptr<map<vector<int>, double>> uncorrected_lls;

        /// For each scored allele, in matrix column order, the number of reads that fit it
        /// best. A read that fits several alleles equally well splits its count between
        /// them.
        vector<double> allele_support;

        /// Mean over reads of the best log-likelihood score any allele gave them, in nats,
        /// the row divisor. It says whether the reads fit any allele here, where GQ says how
        /// far apart the top two genotypes are. Written as BL.
        double mean_best_ln = 0;

        /// Fraction of reads whose best-fitting allele is one of the called alleles,
        /// derived from allele_support. 1.0 when the called genotype accounts for every
        /// read at the site.
        double explained_share = 1.0;

        /// Observed reads over the number the called genotype predicts, from the local
        /// rate and the lengths of the called alleles; 1.0 when the two agree. Written as DR
        /// whether or not the depth term is on. Negative, and not written, where no read
        /// begins in the rate window.
        double depth_ratio = -1.0;

        /// GQ before any discount, the explained share's or the depth discount. Written as GQI.
        double gq_undiscounted = 0;

        /// The factor the depth discount multiplies GQ by (see set_depth_quality); 1.0 where
        /// it does not apply.
        double depth_discount = 1.0;

        /// The ln-likelihood difference between the called genotype and the runner-up, as a
        /// fraction of the largest difference the read term could give between the two (see
        /// AlleleReadLikelihoods::achievable_gap), held at 1 or less and multiplied by the
        /// explained share. In [0, 1], and comparable across depths and ploidies. Negative
        /// when there was nothing to normalise (no reads, or only one possible genotype),
        /// which is different from 0. Written as GQN, except on records whose genotype the
        /// linkage model changed, which get a GQN of their own.
        double gq_fraction = -1.0;

        /// The divisor of `gq_fraction`: the largest ln-likelihood difference the read term
        /// could give between the called genotype and the runner-up (see
        /// AlleleReadLikelihoods::achievable_gap), in nats. 0 where there is no runner-up.
        double achievable_gap = 0.0;

    };

    /// While on, every `genotype` call on this thread that finds a genotype at a site with more
    /// than one candidate allele also scores the site at the other ploidy, filling
    /// `alt_ploidy_best` and `alt_ploidy_info`. The matrix is reused, since it does not
    /// depend on ploidy. The caller turns it off again after the call.
    ///
    /// Only a site in a nested chain needs this, because its ploidy comes from its parent's
    /// genotype, which the barrier settles later. A top-level site's ploidy is fixed, so asking
    /// for it there only costs memory. Thread-local because the sweep calls sites in parallel.
    static void set_want_alt_ploidy(bool on) { want_alt_ploidy = on; }

    /// Tell the `genotype` calls on this thread how many of the sample's haplotypes the region
    /// around the site holds, or 0 to take the site's own ploidy. The depth term needs it: its
    /// rate is per haplotype, and the reads it is measured from come from every haplotype in the
    /// region. Only a nested site needs it set, since a nested site's ploidy counts only the
    /// parent alleles that cross it, where a top-level site's ploidy is the region's.
    /// Thread-local for the same reason as set_want_alt_ploidy.
    static void set_region_ploidy(int ploidy) { region_ploidy = ploidy; }

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
     * Skip no allele when read support is unavailable.
     *
     * Only the traversal finder of `vg call -v`, which enumerates the alleles of an
     * input VCF, asks for this. The inherited version skips alleles whose read support
     * is below a threshold, to keep that enumeration small at dense sites. Without a
     * pack file the support finder reports zero everywhere, so it would skip every
     * allele. With a pack file the inherited version is used.
     */
    virtual function<bool(const SnarlTraversal&, int iteration)> get_skip_allele_fn() const;

    /// Tell the caller that its support finder reports nothing real, so anything
    /// support-derived must be skipped rather than believed.
    void set_support_available(bool available);

    /// Write the matrix for every site to this stream as TSV. Not owned.
    void set_likelihood_dump(ostream* dump_stream);

    /**
     * Scale GQ by the explained share, the fraction of reads whose best allele is a
     * called allele (on unless --no-share-quality).
     *
     * A read that fits an uncalled allele best fits the called genotype and its
     * runner-up about equally, so it barely changes GQ. The discount lowers GQ when
     * the call leaves reads unexplained. The discounted GQ is a score for ranking
     * calls rather than a posterior; GQI keeps the undiscounted value.
     */
    void set_share_discount(bool discount);

    /**
     * Scale GQ by how far the site's read count is from what the call predicts, at
     * records whose called alleles change length by at least `min_length` bp
     * (--depth-quality):
     *
     *     GQ' = GQ * exp(-exponent * |ln DR|)
     *
     * It can only lower GQ, and GQI keeps the undiscounted value. Small variants are
     * left alone because at a short allele the expected read count depends mostly on
     * read length, so DR there reflects coverage noise more than the call. An
     * `exponent` of 0 turns it off.
     */
    void set_depth_quality(double exponent, size_t min_length = 50);

    /**
     * The explained share of the called genotype `called`: the fraction of the site's reads
     * whose best allele is one of its alleles.
     */
    static double explained_share(const ReadLikelihoodCallInfo& info, const vector<int>& called);

    /**
     * GQ for the called genotype `called`, given `gap`, the phred difference between the two
     * best genotypes: `gap` multiplied by the explained share of `called`, unless
     * --no-share-quality, and by `info.depth_discount`.
     */
    double discounted_gq(const ReadLikelihoodCallInfo& info, const vector<int>& called,
                         double gap) const;

    /**
     * The factor `discounted_gq` multiplies the gap by at the called genotype: `info`'s
     * explained share, unless --no-share-quality, times `info.depth_discount`.
     */
    double gq_factor(const ReadLikelihoodCallInfo& info) const;

    /**
     * Recompute `info.gq` from `info.genotype_lls`, for the genotype with the highest
     * likelihood, after something has changed the likelihoods. GQI, GQN, DR and the depth
     * discount keep the values computed with the matrix.
     */
    void recompute_gq(ReadLikelihoodCallInfo& info) const;

    /**
     * Set FILTER to `lowconf` on records whose GQN is below `threshold`
     * (--min-confidence); 0 turns it off.
     *
     * GQN, unlike GQ, does not depend on depth or ploidy, so one threshold suits
     * every contig. Records are marked rather than dropped, because the linkage model
     * can change a record's genotype and GQN after this runs.
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
    /// See set_region_ploidy.
    static thread_local int region_ploidy;


    /// Whether GQ is scaled by the explained-read fraction. See set_share_discount.
    bool share_discount = true;

    /// Exponent on |ln DR| in the depth discount, and the smallest change in allele length
    /// it applies to. Zero exponent disables. See set_depth_quality.
    double depth_quality = 0.0;
    size_t depth_quality_min_length = 50;

    /// GQN below which a record is marked lowconf. Zero disables. See set_min_confidence.
    double min_confidence = 0.0;
};

}

#endif
