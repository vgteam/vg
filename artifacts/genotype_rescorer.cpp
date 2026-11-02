#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>

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
    }
    return phases;
}

/// Find the genotype with the highest likelihood, or null for no genotypes.
///
/// The first of tied genotypes, in map order, wins.
pair<const vector<int>*, double> best_genotype(const map<vector<int>, double>& gl) {
    const vector<int>* best = nullptr;
    double best_ll = -numeric_limits<double>::infinity();
    for (const auto& [genotype, ll] : gl) {
        if (ll > best_ll) {
            best_ll = ll;
            best = &genotype;
        }
    }
    return {best, best_ll};
}

/// Spell a genotype as unphased VCF GT, as `0/1`, or `.` for null.
string spell_genotype(const vector<int>* genotype) {
    if (genotype == nullptr) {
        return ".";
    }
    string out;
    for (size_t i = 0; i < genotype->size(); ++i) {
        if (i) {
            out += "/";
        }
        out += std::to_string((*genotype)[i]);
    }
    return out;
}

/// One row of the ledger file, with the fields it is sorted by.
struct LedgerRow {
    string contig;
    size_t position;
    /// The site's ID, which orders rows at the same position.
    string snarl;
    /// The tab-separated row, without a newline.
    string text;
};

/// Open the ledger file at `path` and write its header, or exit if it cannot be written.
void open_ledger(ofstream& ledger, const string& path, double temper) {
    ledger.open(path);
    if (!ledger) {
        cerr << "error [vg call]: cannot write --regeno-ledger " << path << endl;
        exit(1);
    }
    ledger << "#regeno-ledger-version\t1" << endl;
    ledger << "#temper\t" << temper << endl;
    ledger << "#snarl\tcontig\tposition\tploidy\tcalled\tproposed\tdelta_ln\treads" << endl;
}

/// Write `rows` to the ledger sorted by position, so the file does not depend on threading.
void write_ledger(ofstream& ledger, vector<LedgerRow>& rows) {
    std::ranges::sort(rows, [](const LedgerRow& x, const LedgerRow& y) {
        if (x.contig != y.contig) return x.contig < y.contig;
        if (x.position != y.position) return x.position < y.position;
        return x.snarl < y.snarl;
    });
    for (const LedgerRow& row : rows) {
        ledger << row.text << endl;
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
    /// The summed strand log-odds of each read.
    const LambdaTable& lambda;
    /// What the phase says about each site.
    const unordered_map<size_t, SitePhase>& site_phases;
    /// The record keys of the sites whose strands read phasing swapped.
    const unordered_set<size_t>& phase_flips;
    /// The genotyper that recomputes GQ.
    const ReadLikelihoodSnarlCaller& genotyper;
    const RegenotypeParams& params;
    /// The calibration of every read's strand log-odds.
    StrandCalibration calibration;
    /// Whether the corrected likelihoods replace each site's own.
    ///
    /// When false, they are computed and counted, and then thrown away. `genotype_lls` is what
    /// GL and QUAL are written from, so correcting it in place would change them while the
    /// genotypes stood still.
    bool keep;
    /// Whether to make a ledger row for each site whose best genotype changes.
    bool want_ledger;
};

/**
 * A worker that corrects sites' genotype likelihoods for one round, on one thread.
 *
 * It tallies the sites it corrects in its own counters, the sites whose best genotype changed,
 * and their ledger rows, for merging once every worker is done.
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
    /// A ledger row for each site whose best genotype changed, if the round wants them.
    vector<LedgerRow> ledger_rows;

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

    map<vector<int>, double> before;
    if (round.want_ledger) {
        before = score->genotype_lls;
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
            ? haploid_inclusion_correction(*pe, round.lambda, phase.phase_set, own,
                                           round.calibration, phase.haploid_strand, round.params,
                                           target, counters)
            : phase_aware_correction(*pe, round.lambda, phase.phase_set, phase.strand0_allele, own,
                                     round.calibration, round.params, target, counters);
    if (round.keep && site_moved) {
        // Recompute GQ from the corrected likelihoods. GQI and GQN stay, because the achievable
        // gap of GQN assumes the site's own mixture weights, not per-read ones.
        round.genotyper.recompute_gq(*score);
    }
    if (round.keep && score->alt_ploidy_info != nullptr) {
        correct_other_ploidy(*score->alt_ploidy_info, *pe, phase);
    }
    if (!site_moved) {
        return;
    }
    ++moved;
    if (round.want_ledger) {
        const auto [called, called_ll] = best_genotype(before);
        const auto [proposed, proposed_ll] = best_genotype(target);
        std::ostringstream row;
        row << rec.id << "\t" << rec.reference.path_name << "\t"
            << rec.reference.offset << "\t" << rec.ploidy << "\t" << spell_genotype(called) << "\t"
            << spell_genotype(proposed) << "\t" << (proposed_ll - called_ll) << "\t"
            << pe->num_reads();
        ledger_rows.emplace_back(rec.reference.path_name, (size_t)rec.reference.offset, rec.id,
                                 row.str());
    }
}

void SiteCorrector::correct_other_ploidy(SiteScore& alt, const PhaseReadEvidence& evidence,
                                         const SitePhase& phase) {
    restore_direct_pass(alt, round.genotyper);
    RegenotypeCounters ignored;
    if (phase_aware_correction(evidence, round.lambda, phase.phase_set, phase.strand0_allele, own,
                               round.calibration, round.params, alt.genotype_lls,
                               ignored)) {
        round.genotyper.recompute_gq(alt);
    }
}

}

void GenotypeRescorer::configure(bool enable, const RegenotypeParams& settings, size_t passes,
                                 std::string_view ledger) {
    on = enable;
    regenotype_params = settings;
    max_passes = passes;
    ledger_path = ledger;
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

    const unordered_map<size_t, SitePhase> site_phases = collect_site_phases(
        calls, phase_sites, phase_set_id, regenotype_params.haploid_include);

    const StrandCalibration calibration = choose_calibration(phase_sites, phase_flips, lambda);
    fit.keep(counters);

    ofstream ledger;
    const bool want_ledger = !ledger_path.empty();
    if (want_ledger) {
        open_ledger(ledger, ledger_path, calibration.temper);
    }

    // Only the read-likelihood genotyper makes the scores corrected here.
    const ReadLikelihoodSnarlCaller* genotyper = reader.genotyper;
    size_t moved = 0;
    if (genotyper != nullptr) {
        const CorrectionRound round {lambda, site_phases, phase_flips, *genotyper,
                                     regenotype_params, calibration, max_passes >= 2,
                                     want_ledger};
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
        vector<LedgerRow> rows;
        for (SiteCorrector& worker : workers) {
            merge_counters(worker.counters, counters);
            moved += worker.moved;
            std::ranges::move(worker.ledger_rows, std::back_inserter(rows));
        }
        if (ledger.is_open()) {
            write_ledger(ledger, rows);
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
         << c.order_reversed << " where the reads prefer the other order" << endl;
    if (regenotype_params.haploid_include) {
        cerr << "[vg call] re-genotyping: " << c.haploid_sites
             << " nested haploid chains weighted by whether the reads belong to their strand, "
             << c.haploid_would_move << " would move" << endl;
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
