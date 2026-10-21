#ifndef VG_SITE_GENOTYPER_HPP_INCLUDED
#define VG_SITE_GENOTYPER_HPP_INCLUDED

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "read_likelihood_caller.hpp"

namespace vg::multipass {

using namespace std;

/// The genotyping result of one site from the reads.
///
/// It holds the likelihood of every genotype at the site's ploidy, and at the other ploidy where
/// asked for, and everything the later passes and the VCF record need.
using SiteScore = ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo;

/**
 * A genotyper of single sites from the reads, which gives each site a `SiteScore`.
 *
 * It is a view of a `ReadLikelihoodSnarlCaller`, which it does not own. It returns the score
 * with its own type, so the passes of the multi-pass caller can read the score without casting.
 */
class SiteGenotyper {
public:
    explicit SiteGenotyper(ReadLikelihoodSnarlCaller& genotyper);

    /**
     * Genotype `site` at `ploidies`, choosing among the alleles `travs`.
     *
     * @param ref_trav_idx The index of the reference allele in `travs`.
     * @param enclosing The bounds of the sites enclosing `site`, innermost first.
     * @param sequences Sequences to spell alleles by, in place of their walks, if given.
     * @return The genotype, as a sorted multiset of indices into `travs`, or empty if the site
     * cannot be genotyped. Also the score, which is never null.
     */
    pair<vector<int>, unique_ptr<SiteScore>> genotype(const SiteBounds& site,
                                                      const vector<Traversal>& travs,
                                                      int ref_trav_idx, const Ploidies& ploidies,
                                                      const vector<SiteBounds>& enclosing,
                                                      const AlleleSequences* sequences
                                                      = nullptr) const;

    /// Get the factor that GQ multiplies the likelihood gap by, at the called genotype.
    ///
    /// The gap is the phred difference between the two best genotypes. The factor is the depth
    /// discount, times the explained share of the called genotype when the share discount is on.
    double gq_factor(const SiteScore& score) const;

    /// Recompute `score.gq` from `score.genotype_lls`, after something has changed them.
    void recompute_gq(SiteScore& score) const;

private:
    /// Does the genotyping.
    ReadLikelihoodSnarlCaller& genotyper;
};

}

#endif
