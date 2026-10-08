#include <limits>

#include "read_strand_table.hpp"

namespace vg {

void TemperFit::keep(const RegenotypeCounters& counters) {
    temper = counters.fitted_temper;
    abs_lambda = counters.fit_abs_lambda;
    observed = counters.fit_observed;
    predicted = counters.fit_predicted;
    count = counters.fit_count;
}

void TemperFit::restore(RegenotypeCounters& counters) const {
    counters.fitted_temper = temper;
    counters.fit_abs_lambda = abs_lambda;
    counters.fit_observed = observed;
    counters.fit_predicted = predicted;
    counters.fit_count = count;
}

// Built for the render rather than taken from re-genotyping, which may not have run and whose
// table is built before the sites are final.
void ReadStrandTable::build_lambda(
    const vector<PhaseCall>& calls,
    const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
    double fitted_temper, double fitted_ceiling, const RegenotypeParams& params) {
    lambda.clear();
    lambda_site.clear();
    lambda_phase_set.clear();
    lambda_temper = 0.0;
    lambda_ceiling = 1.0;
    if (phase_sites.empty()) {
        return;
    }
    RegenotypeCounters scratch;
    accumulate_lambda(phase_sites, phase_flips, lambda, scratch);
    for (const PhaseSite& site : phase_sites) {
        lambda_site[site.record_key] = &site;
    }
    // The last PhaseCall written winning, as in `PhaseTable::freeze_for_render`.
    for (const PhaseCall& pc : calls) {
        lambda_phase_set[pc.record_key] = phase_set_id(pc.contig, pc.phase_set);
    }
    // The summed strand log-odds overstate how sure the strand is, so they are tempered. Use the
    // temper re-genotyping fitted, where it ran; otherwise fit one here.
    if (fitted_temper > 0.0) {
        lambda_temper = fitted_temper;
        lambda_ceiling = fitted_ceiling;
    } else {
        double temper = -1.0;
        double ceiling = params.ceiling < 0.0 ? 1.0 : params.ceiling;
        RegenotypeCounters fit_scratch;
        fit_calibration(phase_sites, phase_flips, lambda, params, temper, ceiling, fit_scratch);
        if (fit_scratch.fitted_temper > 0.0) {
            lambda_temper = fit_scratch.fitted_temper;
            lambda_ceiling = fit_scratch.fitted_ceiling;
        }
    }
}

bool ReadStrandTable::site_own_strand_log_odds(size_t record_key,
                                               unordered_map<uint64_t, double>& out) const {
    if (lambda.empty() || lambda_temper <= 0.0) {
        return false;
    }
    const auto site = lambda_site.find(record_key);
    if (site == lambda_site.end()) {
        return false;
    }
    site_own_log_odds(*site->second, phase_flips.count(record_key) != 0, out);
    return true;
}

double ReadStrandTable::read_strand_log_odds(size_t record_key, std::string_view read_name,
                                             const unordered_map<uint64_t, double>* site_own) const {
    if (lambda.empty() || lambda_temper <= 0.0) {
        return 0.0;
    }
    const uint64_t key = (uint64_t)std::hash<std::string_view>{}(read_name);
    const auto found = lambda.find(key);
    const auto ps = lambda_phase_set.find(record_key);
    const size_t phase_set = ps != lambda_phase_set.end() ? ps->second : NO_PHASE_SET;
    if (found == lambda.end()) {
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
    return calibrated_log_odds(value, lambda_temper, lambda_ceiling);
}

}
