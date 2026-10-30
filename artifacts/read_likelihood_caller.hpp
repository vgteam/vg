#ifndef VG_READ_LIKELIHOOD_CALLER_HPP_INCLUDED
#define VG_READ_LIKELIHOOD_CALLER_HPP_INCLUDED

/** \file read_likelihood_caller.hpp
 *
 * A SnarlCaller that genotypes a site from the likelihood of its reads under each
 * genotype, P(reads | genotype), rather than from read depth.
 *
 * The model is described in artifacts/doc/read-likelihood-direct-genotyping.md, and the VCF fields
 * this caller writes in artifacts/doc/read-likelihood-genotyping.md.
 */

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "allele_likelihood.hpp"
#include "linkage_model.hpp"
#include "../src/snarl_caller.hpp"

namespace vg {

using namespace std;

/// The ploidies one site is genotyped with.
struct Ploidies {
    /// The ploidy to genotype the site at.
    int ploidy = 2;
    /// The number of the sample's haplotypes in the region around the site, or 0 to take
    /// `ploidy`.
    ///
    /// The depth term needs it, since its rate is per haplotype, and its reads come from every
    /// haplotype in the region. Only a nested site needs it set, since a nested site's ploidy
    /// counts only the parent alleles that cross it.
    int region_ploidy = 0;
    /// Whether to also score the site at the other ploidy, 2 for ploidy 1 and 1 otherwise.
    ///
    /// That fills `alt_ploidy_best` and `alt_ploidy_info` where the site has more than one
    /// candidate allele, from the same matrix. Only a site in a nested chain needs it, because its
    /// ploidy comes from its parent's genotype, which the linkage pass chooses later.
    bool also_score_other = false;
};

/**
 * A genotyper that scores every genotype of a site from its reads x alleles likelihood matrix.
 *
 * The genotyping uses no read support. A record can write its alleles in another order than the
 * one they were scored in, so the caller writing it says which scored allele each written allele
 * is (see `write_vcf_fields`). A record with an allele that was not scored, such as a star allele,
 * is written without GL, since GL needs an entry for every genotype.
 */
class ReadLikelihoodSnarlCaller {
public:

    /// Make a genotyper of the sites of `graph`, from the likelihoods `likelihood_calculator`
    /// computes. Neither is owned.
    ReadLikelihoodSnarlCaller(const PathHandleGraph& graph,
                              AlleleLikelihoodCalculator& likelihood_calculator);

    ~ReadLikelihoodSnarlCaller();

    /// The genotyping result of one site: the scores of its genotypes, and the fields the record
    /// is written from.
    struct ReadLikelihoodCallInfo : public SnarlCaller::CallInfo {
        ~ReadLikelihoodCallInfo() override = default;

        /// The phred-scaled gap between the best and second-best genotypes, times the discounts
        /// (see `discounted_gq`). Written as GQ.
        double gq = 0;
        /// The ln posterior of the genotype with the highest likelihood, under a uniform prior
        /// over genotypes. Written as GP.
        double posterior = 0;
        /// The number of reads in the site's matrix. Written as DP.
        size_t n_reads = 0;
        /// The ploidy this site was genotyped at.
        int ploidy = 2;

        /// The ln P(reads | G) of every genotype scored, by the sorted multiset of its allele
        /// indices.
        map<vector<int>, double> genotype_lls;
        /// The best genotype at the other ploidy (2 for a site genotyped at ploidy 1, and 1
        /// otherwise), as traversal indices. Empty unless requested with
        /// `Ploidies::also_score_other`.
        vector<int> alt_ploidy_best;

        /// The whole call at the other ploidy, from the same matrix, or null unless requested
        /// with `Ploidies::also_score_other`.
        ///
        /// The site's record can then be built at whichever ploidy the linkage pass gives it. Its
        /// fields that do not depend on ploidy, such as `scored_traversals` and
        /// `allele_support`, are copies of this one's.
        unique_ptr<ReadLikelihoodCallInfo> alt_ploidy_info;

        /// The walks that were scored, in matrix column order.
        ///
        /// The deduplicated traversals handed to update_vcf_info are mapped back to them.
        vector<Traversal> scored_traversals;

        /// The per-read anchor evidence when anchors are written, or null.
        ///
        /// It depends on the matrix, not on the ploidy, so when the linkage pass replaces this
        /// CallInfo with `alt_ploidy_info` it must move it across.
        unique_ptr<AnchorSiteEvidence> anchor_evidence;
        /// The per-read phasing evidence with read phasing and no anchors, or null.
        ///
        /// It depends on the matrix, not on the ploidy, so when the linkage pass replaces this
        /// CallInfo with `alt_ploidy_info` it must move it across too.
        unique_ptr<PhaseReadEvidence> phase_evidence;

        /// Get the per-read evidence that read phasing and re-genotyping need, or null where the
        /// site kept none, or has no reads or alleles.
        ///
        /// When anchors are written the site keeps only `anchor_evidence`, which is converted
        /// into `scratch`. Then `scratch` must outlive the returned pointer.
        const PhaseReadEvidence* read_phasing_evidence(PhaseReadEvidence& scratch) const;

        /// The `genotype_lls` as the direct pass computed them, or null until re-genotyping
        /// changes the site.
        ///
        /// They are saved at the first correction. Each later round of re-genotyping corrects
        /// these, rather than the previous round's values.
        unique_ptr<map<vector<int>, double>> uncorrected_lls;

        /// The number of reads that fit each scored allele best, in matrix column order.
        ///
        /// A read that fits several alleles equally well splits its count between them.
        vector<double> allele_support;

        /// The mean over reads of the best log-likelihood score any allele gave them, in nats.
        /// Written as BL.
        ///
        /// It says whether the reads fit any allele here, where GQ says how far apart the top two
        /// genotypes are.
        double mean_best_ln = 0;

        /// The fraction of reads whose best-fitting allele is one of the called alleles, from
        /// `allele_support`.
        double explained_share = 1.0;

        /// The observed reads over the number the called genotype predicts, from the local rate
        /// and the lengths of the called alleles, or negative if no read begins in the rate
        /// window.
        ///
        /// DR is written from `depth` for the written genotype where that can be computed, and
        /// from this otherwise.
        double depth_ratio = -1.0;

        /// The depth context of the matrix, which gives DR for a genotype given as matrix
        /// columns.
        ///
        /// Where the matrix had no depth context it is empty, with a rate of 0, so it gives a DR
        /// of -1 for every genotype.
        DepthContext depth;

        /// GQ before the explained share and the depth discount. Written as GQI.
        double gq_undiscounted = 0;

        /// The factor the depth discount multiplies GQ by (see set_depth_quality), or 1.0 where
        /// it does not apply.
        double depth_discount = 1.0;

        /// The ln-likelihood difference between the called genotype and the runner-up, as a
        /// fraction of `achievable_gap`, capped at 1 and multiplied by the explained share.
        ///
        /// It is in [0, 1], and comparable across depths and ploidies. It is negative when there
        /// was nothing to normalise, because there were no reads, or only one possible genotype.
        /// Written as GQN, except on records whose genotype the linkage model changed, which get
        /// a GQN of their own.
        double gq_fraction = -1.0;

        /// The largest ln-likelihood difference the read term could give between the called
        /// genotype and the runner-up (see AlleleReadLikelihoods::achievable_gap), in nats, or 0
        /// where there is no runner-up.
        double achievable_gap = 0.0;
    };

    /**
     * Genotype the site `site` at `ploidies`, choosing among its alleles.
     *
     * @param fallback_nodes Nodes to place the site's depth-rate window by, most preferred first,
     * when neither bound of the site is on a reference path.
     * @return The genotype, as a sorted multiset of indices into the site's walks, or empty if the
     * site cannot be genotyped. Also the score, which is never null.
     */
    pair<vector<int>, unique_ptr<ReadLikelihoodCallInfo>> genotype_at(
        const MaterializedSite& site, const Ploidies& ploidies,
        const vector<nid_t>& fallback_nodes);

    /**
     * Write the record's fields from the scores in `call_info`: DP, BL, DR, AD, GL and the rest.
     *
     * Do nothing when `call_info` is not this genotyper's.
     *
     * @param scored_allele For each written allele, the index of the scored allele it is, or -1
     * for one that was not scored, such as a star allele.
     * @param genotype The genotype, in written alleles.
     */
    void write_vcf_fields(const vector<int>& scored_allele, const vector<int>& genotype,
                          const SnarlCaller::CallInfo* call_info, const string& sample_name,
                          vcflib::Variant& variant) const;

    /// Declare in the header of `writer` the FORMAT and FILTER fields that this genotyper
    /// writes.
    void add_header_lines(SortedVCFWriter& writer) const;

    /// Write the matrix for every site to this stream as TSV. Not owned.
    void set_likelihood_dump(ostream* dump_stream);

    /**
     * Scale GQ by the explained share, the fraction of reads whose best allele is a
     * called allele. On by default.
     *
     * A read that fits an uncalled allele best fits the called genotype and its
     * runner-up about equally, so it barely changes GQ. The discount lowers GQ when
     * the call leaves reads unexplained. The discounted GQ is a score for ranking
     * calls rather than a posterior; GQI keeps the undiscounted value.
     */
    void set_share_discount(bool discount);

    /**
     * Scale GQ by how far the site's read count is from what the call predicts, at
     * records whose called alleles change length by at least `min_length` bp:
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
     * Get the explained share of the called genotype `called`, the fraction of the site's reads
     * whose best allele is one of its alleles.
     */
    static double explained_share(const ReadLikelihoodCallInfo& info, const vector<int>& called);

    /**
     * Get GQ for the called genotype `called`, given `gap`, the phred difference between the two
     * best genotypes.
     *
     * That is `gap` times the explained share of `called`, when the share discount is on, times
     * `info.depth_discount`.
     */
    double discounted_gq(const ReadLikelihoodCallInfo& info, const vector<int>& called,
                         double gap) const;

    /**
     * Get the factor `discounted_gq` multiplies the gap by at the called genotype.
     *
     * That is the explained share of `info`, when the share discount is on, times
     * `info.depth_discount`.
     */
    double gq_factor(const ReadLikelihoodCallInfo& info) const;

    /**
     * Recompute `info.gq` from `info.genotype_lls`, for the genotype with the highest
     * likelihood, after something has changed the likelihoods. GQI, GQN, DR and the depth
     * discount keep the values computed with the matrix.
     */
    void recompute_gq(ReadLikelihoodCallInfo& info) const;

    /**
     * Set FILTER to `lowconf` on records whose GQN is below `threshold`; 0 turns it off.
     *
     * GQN, unlike GQ, does not depend on depth or ploidy, so one threshold suits
     * every contig. Records are marked rather than dropped, because the linkage model
     * can change a record's genotype and GQN after this runs.
     */
    void set_min_confidence(double threshold);

    /**
     * Rewrite the GQ and GQN fields of `sample` and the FILTER column of `variant`, a record
     * written by a ReadLikelihoodSnarlCaller, whose GT holds a genotype the linkage model chose in
     * place of the direct call. Before the rewrite, these three describe the direct call.
     *
     * GQ is recomputed from `moved.posterior`, the linkage model's posterior of the chosen
     * genotype. It is multiplied by the GQ factor of the direct call and capped at the record's
     * GQI where it has one. GQN is recomputed from the record's GL, as the margin of the chosen
     * genotype over the best other, divided by the achievable gap of the direct call, multiplied
     * by its explained share and clamped to [-1, 1]. It is "." where it cannot be computed.
     * `moved.direct` holds the factor, gap and share of the direct call. FILTER becomes `lowconf`
     * when the new GQN is below `lowconf_threshold` and PASS otherwise; where there is no new GQN,
     * or `lowconf_threshold` is 0, a `lowconf` FILTER becomes PASS. Every field is read as the
     * record's line prints it.
     *
     * Returns false, and leaves the record unchanged, if it has no FORMAT fields for `sample`.
     */
    static bool rewrite_quality_for_chosen_genotype(vcflib::Variant& variant, const string& sample,
                                                    const LinkageCollector::MovedQuality& moved,
                                                    double lowconf_threshold);

private:
    /**
     * Call the site at ploidy `ploidy` from its matrix, filling the ploidy-dependent fields of
     * `info`.
     *
     * Break exact ties in favour of the all-reference genotype.
     *
     * @param ref_trav_idx The index of the reference allele, or -1 for none.
     * @param trav_count The number of candidate alleles.
     * @return The best genotype, or empty if no genotype was scored.
     */
    vector<int> call_at_ploidy(const AlleleReadLikelihoods& matrix, int ploidy, int ref_trav_idx,
                               size_t trav_count, ReadLikelihoodCallInfo& info) const;

    /**
     * Get the depth discount of the called genotype `called`, or 1.0 where it does not apply.
     *
     * A call's size is the largest change in length of a called allele against the reference
     * allele. The discount applies to calls of at least `depth_quality_min_length`.
     */
    double depth_discount_of(const AlleleReadLikelihoods& matrix, const vector<int>& called,
                             int ref_trav_idx, double depth_ratio) const;

    /// The graph the sites are in.
    const PathHandleGraph& graph;

    /// Builds the matrix of each site.
    AlleleLikelihoodCalculator& likelihood_calculator;

    /// The stream to dump every site's matrix to as TSV, for development, or null.
    ostream* dump_stream = nullptr;

    /// Whether GQ is scaled by the explained share. See set_share_discount.
    bool share_discount = true;

    /// The exponent on |ln DR| in the depth discount, or 0 to turn it off. See
    /// set_depth_quality.
    double depth_quality = 0.0;
    /// The smallest change in allele length the depth discount applies to.
    size_t depth_quality_min_length = 50;

    /// GQN below which a record is marked lowconf. Zero disables. See set_min_confidence.
    double min_confidence = 0.0;
};

namespace multipass {

/// The genotyping result of one site from the reads.
///
/// It holds the likelihood of every genotype at the site's ploidy, and at the other ploidy where
/// asked for, and everything the later passes and the VCF record need.
using SiteScore = ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo;

}

}

#endif
