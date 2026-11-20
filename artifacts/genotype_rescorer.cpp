#include <algorithm>
#include <limits>

#include <omp.h>

#include "genotype_rescorer.hpp"
#include "read_likelihood_caller.hpp"

namespace vg::multipass {

namespace {

/// What the phase says about one staged site, for correcting its likelihoods.
struct SitePhase {
    /// The run-wide ID of the phase set of the site, or `NO_PHASE_SET`.
    ///
    /// A read's strand is usable only at sites of the phase set it was found in.
    size_t phase_set = NO_PHASE_SET;
    /// The run of read phasing (see `PhaseSite::run`) that the site's order is in, or, for a nested
    /// chain at ploidy 1, that its strand is named in, or `NO_RUN`.
    size_t run = NO_RUN;
    /// The traversal on strand 0 at a diploid site, or -1 where the phase gives no order.
    int strand0_allele = -1;
    /// For a nested chain at ploidy 1, +1 if it is on strand 0 of its phase set and -1 if it is on
    /// strand 1. 0 for other sites, which take the diploid correction.
    int haploid_strand = 0;
    /// The read evidence read phasing used at the site, or null if it used none, as at a
    /// homozygote.
    const PhaseSite* phase_site = nullptr;
};

/**
 * Collect what the phase says about each site, keyed by record key.
 *
 * Where a site has several phase calls, the last one is its phase, as in `PhaseTable::index()`.
 * A nested chain's `nested_strand` is in the frame of its phase set, because read phasing swapped
 * it along with its parent's pair.
 *
 * @param haploid_include Whether to give nested ploidy-1 chains their strand, so that they take
 * the haploid correction.
 */
unordered_map<size_t, SitePhase> collect_site_phases(
    const vector<LinkageCollector::PhaseCall>& calls, const vector<PhaseSite>& phase_sites,
    const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
    bool haploid_include) {
    unordered_map<size_t, SitePhase> phases;
    phases.reserve(calls.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : calls) {
        SitePhase& phase = phases[pc.record_key];
        phase.phase_set = phase_set_id(pc.contig, pc.phase_set);
        phase.strand0_allele = pc.ploidy == 2 ? pc.trav_first : -1;
        if (haploid_include && pc.ploidy == 1 && pc.nested_strand >= 0) {
            phase.haploid_strand = pc.nested_strand == 0 ? 1 : -1;
        }
    }
    for (const PhaseSite& ps : phase_sites) {
        phases[ps.record_key].phase_site = &ps;
        phases[ps.record_key].run = ps.run;
    }
    return phases;
}

/**
 * Give each nested chain at ploidy 1 in `phases` the run its strand is named in.
 *
 * Its strand is one of its parent's two, so it is named in the run of the parent's order. A
 * parent at ploidy 1 has the strand of its own parent, so we go up through such parents to the
 * first site that is not one, and take its run, which is `NO_RUN` if it has no read evidence.
 */
void set_haploid_runs(unordered_map<size_t, SitePhase>& phases, StagedSiteTable& staged) {
    unordered_map<size_t, size_t> parent_of;
    staged.for_each([&parent_of](const StagedSite& rec) {
        if (rec.placement.level > 0) {
            parent_of[rec.record_key] = rec.placement.parent_record_key;
        }
    });
    // The run of the first site enclosing `key` that is not a haploid chain, which names the
    // chain's strand.
    auto enclosing_run = [&parent_of, &phases](size_t key) {
        for (auto parent = parent_of.find(key); parent != parent_of.end();
             parent = parent_of.find(parent->second)) {
            const auto above = phases.find(parent->second);
            if (above == phases.end()) {
                return NO_RUN;
            }
            if (above->second.haploid_strand == 0) {
                return above->second.run;
            }
        }
        return NO_RUN;
    };
    for (auto& [key, phase] : phases) {
        if (phase.haploid_strand != 0) {
            phase.run = enclosing_run(key);
        }
    }
}

/**
 * Put the direct pass's likelihoods back into `score` before it is corrected again.
 *
 * Each round corrects the direct pass's likelihoods, not the previous round's. The first round
 * saves them in `uncorrected_lls`, and later rounds restore them. GQ is restored with them, so
 * that a GQ recomputed in an earlier round does not outlive the correction it came from.
 */
void restore_direct_pass(SiteScore& score, const ReadLikelihoodSnarlCaller& genotyper) {
    if (score.uncorrected_lls == nullptr) {
        score.uncorrected_lls.reset(new map<vector<int>, double>(score.genotype_lls));
    } else {
        score.genotype_lls = *score.uncorrected_lls;
        genotyper.recompute_gq(score);
    }
}

/// The inputs of one round of correction, shared by every site.
struct CorrectionRound {
    /// What the round knows about each read's strand.
    const StrandModel& strands;
    /// What the phase says about each site.
    const unordered_map<size_t, SitePhase>& site_phases;
    /// The record keys of the sites whose strands read phasing swapped.
    const unordered_set<size_t>& phase_flips;
    /// The genotyper that recomputes GQ.
    const ReadLikelihoodSnarlCaller& genotyper;
    /// Whether the corrected likelihoods replace each site's own.
    ///
    /// When false, they are computed and counted, and then thrown away. `genotype_lls` is what
    /// GL and QUAL are written from, so correcting it in place would change them while the
    /// genotypes stood still.
    bool keep;
};

/**
 * A worker that corrects sites' genotype likelihoods for one round, on one thread.
 *
 * It tallies the sites it corrects in its own counters, and the sites whose best genotype
 * changed, for merging once every worker is done.
 */
class SiteCorrector {
public:
    explicit SiteCorrector(const CorrectionRound& round) : round(round) {}

    /// Correct one staged site's likelihoods from its reads' strands, and tally the result.
    void correct(StagedSite& rec);

    /// The counters of the corrections made.
    RegenotypeCounters counters;
    /// The number of sites whose best genotype the correction changed.
    size_t moved = 0;

private:
    /// Correct the likelihoods at the other ploidy too, since the linkage pass can move a chain
    /// from ploidy 1 to 2.
    ///
    /// It takes the diploid correction, which is zero at ploidy 1. Its counts are not kept.
    void correct_other_ploidy(SiteScore& alt, const PhaseReadEvidence& evidence,
                              const SitePhase& phase);

    const CorrectionRound& round;
    /// The contribution of the current site to each read's log-odds.
    unordered_map<uint64_t, double> own;
};

void SiteCorrector::correct(StagedSite& rec) {
    SiteScore* score = rec.score;
    if (score == nullptr) {
        return;
    }
    // `converted` must outlive `pe`.
    PhaseReadEvidence converted;
    const PhaseReadEvidence* pe = score->read_phasing_evidence(converted);
    if (pe == nullptr) {
        return;
    }
    SitePhase phase;
    if (const auto found = round.site_phases.find(rec.record_key);
        found != round.site_phases.end()) {
        phase = found->second;
    }
    // Subtract the site's own term from its reads' log-odds, so that the site does not confirm
    // its own genotype. A site with no `PhaseSite` added nothing to the log-odds.
    own.clear();
    if (phase.phase_site != nullptr) {
        site_own_log_odds(*phase.phase_site, round.phase_flips.contains(rec.record_key), own);
    }

    map<vector<int>, double> scratch;
    if (round.keep) {
        restore_direct_pass(*score, round.genotyper);
    } else {
        scratch = score->genotype_lls;
    }
    map<vector<int>, double>& target = round.keep ? score->genotype_lls : scratch;
    const bool site_moved =
        phase.haploid_strand != 0
            ? haploid_inclusion_correction(*pe, round.strands, RunKey {phase.phase_set, phase.run},
                                           phase.haploid_strand, own, target, counters)
            : phase_aware_correction(*pe, round.strands, RunKey {phase.phase_set, phase.run},
                                     phase.strand0_allele, own, target, counters);
    if (round.keep && site_moved) {
        // Recompute GQ from the corrected likelihoods. GQI and GQN stay, because the achievable
        // gap of GQN assumes the site's own mixture weights, not per-read ones.
        round.genotyper.recompute_gq(*score);
    }
    if (round.keep && score->alt_ploidy_info != nullptr) {
        correct_other_ploidy(*score->alt_ploidy_info, *pe, phase);
    }
    if (site_moved) {
        ++moved;
    }
}

void SiteCorrector::correct_other_ploidy(SiteScore& alt, const PhaseReadEvidence& evidence,
                                         const SitePhase& phase) {
    restore_direct_pass(alt, round.genotyper);
    RegenotypeCounters ignored;
    if (phase_aware_correction(evidence, round.strands, RunKey {phase.phase_set, phase.run},
                               phase.strand0_allele, own, alt.genotype_lls, ignored)) {
        round.genotyper.recompute_gq(alt);
    }
}

}

void GenotypeRescorer::configure(bool enable, const RegenotypeParams& settings, size_t passes) {
    on = enable;
    regenotype_params = settings;
    max_passes = passes;
}

void GenotypeRescorer::set_site_reader(SiteReader site_reader) {
    reader = std::move(site_reader);
}

StrandCalibration GenotypeRescorer::choose_calibration(const vector<PhaseSite>& phase_sites,
                                                       const unordered_set<size_t>& phase_flips,
                                                       const LambdaTable& lambda) {
    const StrandCalibration given {
        regenotype_params.temper,
        regenotype_params.ceiling < 0.0 ? 1.0 : regenotype_params.ceiling};
    if (given.temper < 0.0) {
        // A negative temper asks for a fit. Fit on the first round only, because fitting again
        // each round would feed each round's genotypes into the next fit.
        if (counters.calibration.temper > 0.0) {
            return counters.calibration;
        }
        return fit_calibration(phase_sites, phase_flips, lambda, regenotype_params, counters);
    }
    counters.calibration = given;
    return given;
}

bool GenotypeRescorer::rescore(
    StagedSiteTable& staged, const PhaseTable& phases, const ReadStrandTable& strands,
    TemperFit& fit, const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
    bool show_calibration) {
    const vector<LinkageCollector::PhaseCall>& calls = phases.calls();
    if (!on || calls.empty()) {
        return false;
    }
    const vector<PhaseSite>& phase_sites = strands.sites();
    const unordered_set<size_t>& phase_flips = strands.flips();
    // Reset the counters before `accumulate_lambda` counts reads into them, so that the report
    // describes this round. Bring back the calibration fitted on the first round.
    counters = RegenotypeCounters();
    fit.restore(counters);

    LambdaTable lambda;
    accumulate_lambda(phase_sites, phase_flips, lambda, counters);

    unordered_map<size_t, SitePhase> site_phases = collect_site_phases(
        calls, phase_sites, phase_set_id, regenotype_params.haploid_include);
    set_haploid_runs(site_phases, staged);

    const StrandCalibration calibration = choose_calibration(phase_sites, phase_flips, lambda);
    fit.keep(counters);

    // Only the read-likelihood genotyper makes the scores corrected here.
    const ReadLikelihoodSnarlCaller* genotyper = reader.genotyper;
    size_t moved = 0;
    if (genotyper != nullptr) {
        const StrandModel model {lambda, calibration, strands.run_correlations(),
                                 regenotype_params};
        const CorrectionRound round {model, site_phases, phase_flips, *genotyper,
                                     max_passes >= 2};
        // Stride the records across one worker per queue, and merge the workers in queue order.
        const vector<StagedSite*> all_records = staged.in_order();
        const size_t n_queues = max<size_t>(1, staged.queue_count());
        vector<SiteCorrector> workers(n_queues, SiteCorrector(round));
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t qi = 0; qi < n_queues; ++qi) {
            for (size_t ri = qi; ri < all_records.size(); ri += n_queues) {
                workers[qi].correct(*all_records[ri]);
            }
        }
        for (SiteCorrector& worker : workers) {
            merge_counters(worker.counters, counters);
            moved += worker.moved;
        }
    }

    report(calibration, show_calibration);
    return moved > 0;
}

void GenotypeRescorer::report(const StrandCalibration& calibration, bool show_calibration) const {
    const RegenotypeCounters& c = counters;
    cerr << "[vg call] re-genotyping: temper " << calibration.temper << ", ceiling "
         << calibration.ceiling << ", "
         << c.reads_with_lambda
         << " reads carry a strand log-odds (" << c.reads_singleton
         << " span one site, so are inert; " << c.reads_multi_phase_set << " span two blocks), "
         << c.sites_corrected << " of " << c.sites_considered << " sites corrected, "
         << c.sites_would_move << " would move (" << c.moved_hom_to_het << " hom->het, "
         << c.moved_het_to_hom << " het->hom, " << c.moved_het_to_het << " het->het), "
         << c.order_reversed << " where the reads prefer the other order, "
         << c.sites_split_runs << " whose reads were scored in more than one group by run" << endl;
    if (regenotype_params.haploid_include) {
        cerr << "[vg call] re-genotyping: " << c.haploid_sites
             << " nested haploid chains weighted by whether the reads belong to their strand, "
             << c.haploid_would_move << " would move, " << c.haploid_split_runs
             << " whose reads were scored in more than one group by run" << endl;
    }
    if (show_calibration && !c.fit_count.empty()) {
        cerr << "[vg call] re-genotyping calibration, |Lambda| / observed / predicted / n:";
        for (size_t i = 0; i < c.fit_count.size(); ++i) {
            cerr << "  " << c.fit_abs_lambda[i] << " " << c.fit_observed[i] << " "
                 << c.fit_predicted[i] << " " << c.fit_count[i];
        }
        cerr << endl;
    }
}

}
