#include "site_genotyper.hpp"

namespace vg {

SiteGenotyper::SiteGenotyper(ReadLikelihoodSnarlCaller& genotyper) : genotyper(genotyper) {
}

pair<vector<int>, unique_ptr<SiteScore>> SiteGenotyper::genotype(
    const Snarl& site, const vector<SnarlTraversal>& travs, int ref_trav_idx,
    const Ploidies& ploidies, const string& ref_path_name, pair<size_t, size_t> ref_range) const {
    return genotyper.genotype_at(site, travs, ref_trav_idx, ploidies, ref_path_name, ref_range);
}

double SiteGenotyper::gq_factor(const SiteScore& score) const {
    return genotyper.gq_factor(score);
}

void SiteGenotyper::recompute_gq(SiteScore& score) const {
    genotyper.recompute_gq(score);
}

}
