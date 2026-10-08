#include <limits>

#include "vcf_output_caller.hpp"

namespace vg {

// Each read's strand log-odds for the render, used by the anchors. Built here rather than taken
// from re-genotyping, which may not have run and whose table is built before `phase_sites` is
// final.
void VCFOutputCaller::build_render_lambda() {
    render_lambda.clear();
    render_lambda_site.clear();
    render_lambda_phase_set.clear();
    render_lambda_temper = 0.0;
    render_lambda_ceiling = 1.0;
    if (phase_sites.empty()) {
        return;
    }
    RegenotypeCounters scratch;
    accumulate_lambda(phase_sites, phase_flips, render_lambda, scratch);
    for (const PhaseSite& site : phase_sites) {
        render_lambda_site[site.record_key] = &site;
    }
    // The last PhaseCall written winning, as in `build_render_phases`.
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        render_lambda_phase_set[pc.record_key] = phase_set_id(pc.contig, pc.phase_set);
    }
    // The summed strand log-odds overstate how sure the strand is, so they are tempered. Use the
    // temper re-genotyping fitted, where it ran; otherwise fit one here.
    if (regenotype_counters.fitted_temper > 0.0) {
        render_lambda_temper = regenotype_counters.fitted_temper;
        render_lambda_ceiling = regenotype_counters.fitted_ceiling;
    } else {
        double temper = -1.0;
        double ceiling = regenotype_params.ceiling < 0.0 ? 1.0 : regenotype_params.ceiling;
        RegenotypeCounters fit_scratch;
        fit_calibration(phase_sites, phase_flips, render_lambda, regenotype_params, temper, ceiling,
                        fit_scratch);
        if (fit_scratch.fitted_temper > 0.0) {
            render_lambda_temper = fit_scratch.fitted_temper;
            render_lambda_ceiling = fit_scratch.fitted_ceiling;
        }
    }
}

bool VCFOutputCaller::site_own_strand_log_odds(size_t record_key, unordered_map<uint64_t, double>& out) const {
    if (render_lambda.empty() || render_lambda_temper <= 0.0) {
        return false;
    }
    const auto site = render_lambda_site.find(record_key);
    if (site == render_lambda_site.end()) {
        return false;
    }
    site_own_log_odds(*site->second, phase_flips.count(record_key) != 0, out);
    return true;
}

double VCFOutputCaller::read_strand_log_odds(size_t record_key, std::string_view read_name,
                                             const unordered_map<uint64_t, double>* site_own) const {
    if (render_lambda.empty() || render_lambda_temper <= 0.0) {
        return 0.0;
    }
    const uint64_t key = (uint64_t)std::hash<std::string_view>{}(read_name);
    const auto found = render_lambda.find(key);
    const auto ps = render_lambda_phase_set.find(record_key);
    const size_t phase_set = ps != render_lambda_phase_set.end() ? ps->second : NO_PHASE_SET;
    if (found == render_lambda.end()) {
        // The read reached no phased site.
        return 0.0;
    }
    if (!read_strand_usable(found->second, phase_set)) {
        // The read has a strand, but in another phase set, whose strands do not correspond to
        // this site's. NaN rather than 0, because a split homozygous site drops such a read but
        // places one with no strand by a coin (see `build_site_anchors`).
        return std::numeric_limits<double>::quiet_NaN();
    }
    double value = found->second.lambda;
    size_t sites = found->second.sites;
    // Subtract this record's own contribution, so that a site is not judged by its own evidence;
    // if it was the only one, there is nothing left.
    unordered_map<uint64_t, double> built;
    if (site_own == nullptr && site_own_strand_log_odds(record_key, built)) {
        site_own = &built;
    }
    if (site_own != nullptr) {
        const auto mine = site_own->find(key);
        if (mine != site_own->end()) {
            value -= mine->second;
            if (sites > 0) {
                --sites;
            }
        }
    }
    if (sites == 0) {
        return 0.0;
    }
    return calibrated_log_odds(value, render_lambda_temper, render_lambda_ceiling);
}

}
