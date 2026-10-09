#include <fstream>
#include <limits>
#include <sstream>

#include <omp.h>

#include "genotype_rescorer.hpp"
#include "read_likelihood_caller.hpp"

namespace vg {

void GenotypeRescorer::configure(bool on, const RegenotypeParams& params, size_t passes,
                                 const string& ledger) {
    this->on = on;
    this->regenotype_params = params;
    this->max_passes = passes;
    this->ledger_path = ledger;
}

void GenotypeRescorer::set_site_reader(SiteReader reader) {
    this->reader = std::move(reader);
}

bool GenotypeRescorer::rescore(StagedSiteTable& staged, const PhaseTable& phases,
                               const ReadStrandTable& strands, TemperFit& fit,
                               const function<size_t(const string& contig, size_t phase_set)>& phase_set_id,
                               bool show_calibration) {
    const vector<LinkageCollector::PhaseCall>& calls = phases.calls();
    if (!on || calls.empty()) {
        return false;
    }
    const vector<PhaseSite>& phase_sites = strands.sites();
    const unordered_set<size_t>& phase_flips = strands.flips();
    // Reset the counters first, before `accumulate_lambda` fills the read counts, so that the report
    // describes this round. The calibration table and fitted temper are kept: they are set once, on
    // the first round.
    this->counters = RegenotypeCounters();
    fit.restore(this->counters);

    // Lambda over every site read phasing covered, in one pass, into a table keyed by read.
    LambdaTable lambda;
    accumulate_lambda(phase_sites, phase_flips, lambda, this->counters);

    // Each site's phase set, the last PhaseCall written winning, as in
    // `PhaseTable::freeze_for_render`. A read's strand is usable only at sites of the phase set it
    // was found in.
    unordered_map<size_t, size_t> site_phase_set;
    site_phase_set.reserve(calls.size() * 2);
    // And the allele the chain puts on strand 0 at each diploid site, against which the reads'
    // preferred order is reported.
    unordered_map<size_t, int> site_strand0;
    site_strand0.reserve(calls.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : calls) {
        site_phase_set[pc.record_key] = phase_set_id(pc.contig, pc.phase_set);
        site_strand0[pc.record_key] = pc.ploidy == 2 ? pc.trav_first : -1;
    }

    // Which strand of its parent each nested ploidy-1 chain sits on. `nested_strand` was set in the
    // linkage pass and corrected when its parent's pair was swapped, so it is in the same frame as
    // Lambda for reads of the chain's phase set: strand 0 of that phase set.
    unordered_map<size_t, int> haploid_strand;
    if (regenotype_params.haploid_include) {
        for (const LinkageCollector::PhaseCall& pc : calls) {
            if (pc.ploidy == 1 && pc.nested_strand >= 0) {
                haploid_strand[pc.record_key] = pc.nested_strand == 0 ? 1 : -1;
            }
        }
    }

    double temper = regenotype_params.temper;
    double ceiling = regenotype_params.ceiling < 0.0 ? 1.0 : regenotype_params.ceiling;
    if (temper < 0.0) {
        // The temper is fitted once, on the first round, and kept. It describes how reliable the
        // reads' summed strand log-odds are, not which genotypes are called, and fitting it again
        // each round would feed each round's result into the next fit.
        if (this->counters.fitted_temper > 0.0) {
            temper = this->counters.fitted_temper;
            ceiling = this->counters.fitted_ceiling;
        } else {
            fit_calibration(phase_sites, phase_flips, lambda, regenotype_params, temper, ceiling,
                            this->counters);
        }
    } else {
        this->counters.fitted_temper = temper;
        this->counters.fitted_ceiling = ceiling;
    }
    fit.keep(this->counters);

    // Each site's own PhaseSite, so that its term can be subtracted from its reads' log-odds.
    unordered_map<size_t, const PhaseSite*> site_by_key;
    site_by_key.reserve(phase_sites.size() * 2);
    for (const PhaseSite& ps : phase_sites) {
        site_by_key[ps.record_key] = &ps;
    }

    ofstream ledger;
    const bool want_ledger = !ledger_path.empty();
    if (want_ledger) {
        ledger.open(ledger_path);
        if (!ledger) {
            cerr << "error [vg call]: cannot write --regeno-ledger " << ledger_path << endl;
            exit(1);
        }
        ledger << "#regeno-ledger-version\t1" << endl;
        ledger << "#temper\t" << temper << endl;
        ledger << "#snarl\tcontig\tposition\tploidy\tcalled\tproposed\tdelta_ln\treads" << endl;
    }

    // Parallel over one flat list of records, strided across threads; `lambda`, `site_by_key` and
    // `phase_flips` are read only. Counters and ledger rows are kept per thread and merged
    // afterwards.
    const vector<StagedSite*> all_records = staged.in_order();
    const size_t n_queues = max<size_t>(1, staged.queue_count());
    // Only a read-likelihood caller makes the CallInfos corrected below.
    const auto* rl_caller = dynamic_cast<const ReadLikelihoodSnarlCaller*>(reader.caller);
    vector<RegenotypeCounters> thread_counters(n_queues);
    // Ledger rows are sorted before writing, since which thread handles a record depends on
    // scheduling.
    struct LedgerRow { string contig; size_t position; string snarl; string text; };
    vector<vector<LedgerRow>> thread_ledger(n_queues);
    vector<size_t> thread_moved(n_queues, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t qi = 0; qi < n_queues; ++qi) {
        RegenotypeCounters& counters = thread_counters[qi];
        unordered_map<uint64_t, double> own;
        size_t moved = 0;
        for (size_t ri = qi; ri < all_records.size(); ri += n_queues) {
            StagedSite& rec = *all_records[ri];
            auto* info = dynamic_cast<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                rec.call_info.get());
            if (info == nullptr || rl_caller == nullptr) {
                continue;
            }
            // `converted` belongs to this iteration; nothing may point into it afterwards.
            PhaseReadEvidence converted;
            const PhaseReadEvidence* pe = info->read_phasing_evidence(converted);
            if (pe == nullptr) {
                continue;
            }
            // This site's own term, or none: a homozygote has no PhaseSite and contributed nothing to
            // Lambda, so there is nothing to subtract, and it can be corrected into a heterozygote.
            own.clear();
            auto found_site = site_by_key.find(rec.record_key);
            if (found_site != site_by_key.end()) {
                site_own_log_odds(*found_site->second, phase_flips.count(rec.record_key) != 0,
                                  own);
            }

            // The likelihoods before correction, copied only when something reads them.
            map<vector<int>, double> before;
            if (want_ledger) {
                before = info->genotype_lls;
            }
            // At --regeno-passes 1 the correction is computed and reported, and nothing is kept.
            // `genotype_lls` is what GL and QUAL are written from, so correcting it in place would
            // change them while the genotypes stood still.
            const bool keep = max_passes >= 2;
            map<vector<int>, double> scratch;
            if (keep) {
                // Correct the direct pass's likelihoods every round, not the previous round's: the first
                // round saves them, and later rounds restore them before correcting. GQ is
                // restored with them, so that a GQ recomputed in an earlier round does not outlive
                // the correction it was computed from.
                if (info->uncorrected_lls == nullptr) {
                    info->uncorrected_lls.reset(
                        new map<vector<int>, double>(info->genotype_lls));
                } else {
                    info->genotype_lls = *info->uncorrected_lls;
                    rl_caller->recompute_gq(*info);
                }
            } else {
                scratch = info->genotype_lls;
            }
            map<vector<int>, double>& target = keep ? info->genotype_lls : scratch;
            const auto ps = site_phase_set.find(rec.record_key);
            const size_t phase_set = ps != site_phase_set.end() ? ps->second : NO_PHASE_SET;
            const auto s0 = site_strand0.find(rec.record_key);
            const int strand0_allele = s0 != site_strand0.end() ? s0->second : -1;
            const auto hap = haploid_strand.find(rec.record_key);
            const bool site_moved =
                hap != haploid_strand.end()
                    ? haploid_inclusion_correction(*pe, lambda, phase_set, own, temper, ceiling,
                                                   hap->second, regenotype_params, target,
                                                   counters)
                    : phase_aware_correction(*pe, lambda, phase_set, strand0_allele, own, temper,
                                             ceiling, regenotype_params, target, counters);
            if (keep && site_moved) {
                // The correction changed the best genotype, so GQ is recomputed from the corrected
                // likelihoods, as the direct pass computes it. GQI and GQN are not: GQN's achievable gap
                // assumes the site's own mixture weights, not per-read ones.
                rl_caller->recompute_gq(*info);
            }
            // Both ploidies, as in the direct pass, since the linkage pass can move a chain from
            // ploidy 1 to 2. At ploidy 1 the correction is zero, but it is applied the same way.
            if (keep && info->alt_ploidy_info != nullptr) {
                auto& alt = *info->alt_ploidy_info;
                if (alt.uncorrected_lls == nullptr) {
                    alt.uncorrected_lls.reset(new map<vector<int>, double>(alt.genotype_lls));
                } else {
                    alt.genotype_lls = *alt.uncorrected_lls;
                    rl_caller->recompute_gq(alt);
                }
                RegenotypeCounters ignored;
                if (phase_aware_correction(*pe, lambda, phase_set, strand0_allele, own, temper,
                                           ceiling, regenotype_params, alt.genotype_lls,
                                           ignored)) {
                    rl_caller->recompute_gq(alt);
                }
            }
            if (!site_moved) {
                continue;
            }
            ++moved;
            if (want_ledger) {
                auto best_of = [](const map<vector<int>, double>& gl) {
                    const vector<int>* b = nullptr;
                    double v = -numeric_limits<double>::infinity();
                    for (const auto& kv : gl) {
                        if (kv.second > v) { v = kv.second; b = &kv.first; }
                    }
                    return std::make_pair(b, v);
                };
                auto spell = [](const vector<int>* g) {
                    string out;
                    if (g == nullptr) {
                        return string(".");
                    }
                    for (size_t i = 0; i < g->size(); ++i) {
                        out += (i ? "/" : "") + std::to_string((*g)[i]);
                    }
                    return out;
                };
                const auto a = best_of(before);
                const auto b = best_of(target);
                const string snarl_id = reader.name(rec.snarl);
                std::ostringstream row;
                row << snarl_id << "\t" << rec.ref_path_name << "\t"
                    << rec.ref_offset << "\t" << rec.ploidy << "\t" << spell(a.first) << "\t"
                    << spell(b.first) << "\t" << (b.second - a.second) << "\t"
                    << pe->num_reads();
                thread_ledger[qi].push_back(
                    LedgerRow{rec.ref_path_name, (size_t)rec.ref_offset, snarl_id, row.str()});
            }
        }
        thread_moved[qi] = moved;
    }
    size_t moved = 0;
    vector<LedgerRow> rows;
    for (size_t qi = 0; qi < n_queues; ++qi) {
        merge_counters(thread_counters[qi], this->counters);
        moved += thread_moved[qi];
        std::move(thread_ledger[qi].begin(), thread_ledger[qi].end(), std::back_inserter(rows));
    }
    if (ledger.is_open()) {
        // The snarl ID breaks ties between records at the same position.
        std::sort(rows.begin(), rows.end(), [](const LedgerRow& x, const LedgerRow& y) {
            if (x.contig != y.contig) return x.contig < y.contig;
            if (x.position != y.position) return x.position < y.position;
            return x.snarl < y.snarl;
        });
        for (const LedgerRow& row : rows) {
            ledger << row.text << endl;
        }
        ledger.close();
    }

    const RegenotypeCounters& c = this->counters;
    cerr << "[vg call] re-genotyping: temper " << temper << ", ceiling " << ceiling << ", "
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
        // Only under --progress: the calibration table is a diagnostic.
        cerr << "[vg call] re-genotyping calibration, |Lambda| / observed / predicted / n:";
        for (size_t i = 0; i < c.fit_count.size(); ++i) {
            cerr << "  " << c.fit_abs_lambda[i] << " " << c.fit_observed[i] << " "
                 << c.fit_predicted[i] << " " << c.fit_count[i];
        }
        cerr << endl;
    }
    return moved > 0;
}

}
