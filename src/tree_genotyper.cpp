#include "flow_caller.hpp"

namespace vg {

pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> FlowCaller::genotype_site(
    const Snarl& site, const vector<SnarlTraversal>& travs, int ref_trav_idx,
    const Ploidies& ploidies, const string& ref_path_name, pair<size_t, size_t> ref_range,
    SiteScore*& score) {
    score = nullptr;
    if (site_genotyper == nullptr) {
        // Another genotyper, which takes the ploidy alone.
        return snarl_caller.genotype(site, travs, ref_trav_idx, ploidies.ploidy, ref_path_name,
                                     ref_range);
    }
    auto called =
        site_genotyper->genotype(site, travs, ref_trav_idx, ploidies, ref_path_name, ref_range);
    score = called.second.get();
    return make_pair(std::move(called.first),
                     unique_ptr<SnarlCaller::CallInfo>(std::move(called.second)));
}

// The CallInfo is kept because update_vcf_info reads it when the record is rendered, to map the
// written alleles back to matrix columns, index GL and compute QUAL.
unique_ptr<StagedSite> FlowCaller::stage_render_record(
        const Snarl& snarl, const vector<int>& trav_genotype, int ref_trav_idx,
        unique_ptr<SnarlCaller::CallInfo>& call_info, SiteScore* score,
        const string& ref_path_name, int ref_offset, int ploidy) {
    if (!staged_sites.active()) {
        return nullptr;
    }
    unique_ptr<StagedSite> rec(new StagedSite());
    rec->snarl = snarl;
    rec->ref_path_name = ref_path_name;
    rec->ref_offset = ref_offset;
    rec->ref_trav_idx = ref_trav_idx;
    rec->genotype = trav_genotype;
    rec->ploidy = ploidy;
    rec->record_key = record_key_of(snarl);
    rec->level = 0;
    rec->set_call(std::move(call_info), score);
    // `travs` is not moved here: descent runs after the emit and reads `travs` to find which
    // children the called alleles reach. The caller completes the record after descent.
    return rec;
}

}
