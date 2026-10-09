#ifndef VG_SITE_GENOTYPER_HPP_INCLUDED
#define VG_SITE_GENOTYPER_HPP_INCLUDED

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "read_likelihood_caller.hpp"

namespace vg {

using namespace std;

/// One site's genotyping result from the reads: the likelihood of every genotype, at the site's
/// ploidy and, where asked for, at the other, and everything the later passes and the record read.
using SiteScore = ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo;

/**
 * Genotypes one site from the reads, with explicit ploidies, and returns its typed `SiteScore`,
 * so that the passes after the direct pass read the score without asking which genotyper made
 * it. It is the read-likelihood genotyper (`ReadLikelihoodSnarlCaller`), which it does not own.
 */
class SiteGenotyper {
public:
    explicit SiteGenotyper(ReadLikelihoodSnarlCaller& genotyper);

    /// Genotype `site`, whose candidate alleles are `travs`, the reference among them at
    /// `ref_trav_idx`, at `ploidies`. The genotype is a sorted multiset of indices into `travs`,
    /// empty where the site cannot be genotyped. The score is never null.
    pair<vector<int>, unique_ptr<SiteScore>> genotype(const Snarl& site,
                                                      const vector<SnarlTraversal>& travs,
                                                      int ref_trav_idx, const Ploidies& ploidies,
                                                      const string& ref_path_name,
                                                      pair<size_t, size_t> ref_range) const;

    /// The factor GQ is the gap times, at the called genotype: the explained share, unless
    /// --no-share-quality, times the depth discount.
    double gq_factor(const SiteScore& score) const;

    /// Recompute `score.gq` from `score.genotype_lls`, after something has changed them.
    void recompute_gq(SiteScore& score) const;

private:
    ReadLikelihoodSnarlCaller& genotyper;
};

}

#endif
