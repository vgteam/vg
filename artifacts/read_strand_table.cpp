#include <limits>

#include "read_strand_table.hpp"

namespace vg::multipass {

void TemperFit::keep(const RegenotypeCounters& counters) {
    calibration = counters.calibration;
    abs_lambda = counters.fit_abs_lambda;
    observed = counters.fit_observed;
    predicted = counters.fit_predicted;
    count = counters.fit_count;
}

void TemperFit::restore(RegenotypeCounters& counters) const {
    counters.calibration = calibration;
    counters.fit_abs_lambda = abs_lambda;
    counters.fit_observed = observed;
    counters.fit_predicted = predicted;
    counters.fit_count = count;
}

void ReadStrandTable::build_lambda(
    const vector<PhaseCall>& calls,
    const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
    const StrandCalibration& fitted, const RegenotypeParams& params) {
    lambda.clear();
    lambda_site.clear();
    lambda_phase_set.clear();
    calibration = StrandCalibration();
    if (phase_sites.empty()) {
        return;
    }
    RegenotypeCounters scratch;
    accumulate_lambda(phase_sites, phase_flips, lambda, scratch);
    for (const PhaseSite& site : phase_sites) {
        lambda_site[site.record_key] = &site;
    }
    // Let the last PhaseCall for a record win, as `PhaseTable::freeze_for_render` does.
    for (const PhaseCall& pc : calls) {
        lambda_phase_set[pc.record_key] = phase_set_id(pc.contig, pc.phase_set);
    }
    // The summed strand log-odds overstate how sure the strand is, so we temper them. Use the
    // calibration that re-genotyping used, if it ran, or else fit one here.
    if (fitted.temper > 0.0) {
        calibration = fitted;
    } else {
        RegenotypeCounters fit_scratch;
        const StrandCalibration fit =
            fit_calibration(phase_sites, phase_flips, lambda, params, fit_scratch);
        if (fit.temper > 0.0) {
            calibration = fit;
        }
    }
}

bool ReadStrandTable::site_own_strand_log_odds(size_t record_key,
                                               unordered_map<uint64_t, double>& out) const {
    if (lambda.empty() || calibration.temper <= 0.0) {
        return false;
    }
    const auto site = lambda_site.find(record_key);
    if (site == lambda_site.end()) {
        return false;
    }
    site_own_log_odds(*site->second, phase_flips.contains(record_key), out);
    return true;
}

double ReadStrandTable::read_strand_log_odds(size_t record_key, std::string_view read_name,
                                             const unordered_map<uint64_t, double>* site_own) const {
    if (lambda.empty() || calibration.temper <= 0.0) {
        return 0.0;
    }
    const auto key = (uint64_t)std::hash<std::string_view>{}(read_name);
    const auto found = lambda.find(key);
    const auto ps = lambda_phase_set.find(record_key);
    const size_t phase_set = ps != lambda_phase_set.end() ? ps->second : NO_PHASE_SET;
    if (found == lambda.end()) {
        // The read lies on no phased site.
        return 0.0;
    }
    if (!read_strand_usable(found->second, phase_set)) {
        // The strand of the read belongs to another phase set, which labels its strands
        // independently of this site's phase set. Return NaN, because a split homozygous site
        // drops such a read, but places a read with no strand by a coin (see
        // `build_site_anchors`).
        return std::numeric_limits<double>::quiet_NaN();
    }
    double value = found->second.lambda;
    size_t sites = found->second.sites;
    // Subtract the contribution of this site, so that a site is not judged by its own evidence.
    // Keep `built` in this scope, because `site_own` can point to it.
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
    return calibration.log_odds(value);
}

}
