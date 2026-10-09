#include <cmath>
#include <limits>
#include <set>

#include <omp.h>

#include "flow_caller.hpp"
#include "read_likelihood_caller.hpp"

namespace vg {

// The anchor gqn column for a staged site.
//
// `gq_fraction` was computed in the direct pass for the reads' best genotype, so on a record whose
// genotype the linkage model changed, it describes the abandoned genotype. Such a record gets the
// signed margin of its chosen genotype instead, as the VCF's GQN does (see
// ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype).
//
// Returns the direct pass's value when the model did not change the call; the recomputed signed margin
// when it did; and NaN, written as ".", when it did but the margin cannot be recomputed.
double FlowCaller::anchor_gqn_for(const StagedSite& rec,
                                  const vector<int>& chosen) const {
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
    // The direct pass's value, with its "no gap to normalise" value (-1) turned into NaN, written as ".",
    // so that it stays distinct from the signed range [-1, 1].
    const double direct_value = (info == nullptr || info->gq_fraction < 0.0)
        ? std::numeric_limits<double>::quiet_NaN()
        : info->gq_fraction;
    const double blank = std::numeric_limits<double>::quiet_NaN();
    if (!linker.enabled()) {
        return direct_value;
    }
    const auto& moved = linker.collector()->moved_quality();
    const auto found = moved.find(rec.record_key);
    if (found == moved.end()) {
        return direct_value;   // linkage left the call alone, so the direct pass's value still holds
    }
    // The model changed the call, so the direct pass's value describes the wrong genotype, and any
    // failure below gives NaN rather than falling back to it.
    if (info == nullptr || info->genotype_lls.empty()) {
        return blank;
    }
    // The divisor and share the VCF's GQN uses (see
    // ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype), so that the two agree.
    const LinkageCollector::DirectQuality& direct = found->second.direct;
    if (!(direct.achievable_gap > 0.0)) {
        return blank;   // no scale, and no honest pre-linkage value to fall back on
    }
    const double achievable_phred = 10.0 * direct.achievable_gap / log(10.0);

    // The chosen genotype, not rec.genotype, which is the direct pass's call before the linkage model:
    // the reads prefer that call, so its margin would have the wrong sign.
    vector<int> called = chosen;
    sort(called.begin(), called.end());
    const auto mine = info->genotype_lls.find(called);
    if (mine == info->genotype_lls.end()) {
        return blank;
    }
    // Only genotypes over the written alleles, as the VCF's GL has: the reference traversal and
    // the ones the chosen genotype names. A traversal with no ALT could otherwise beat the call.
    set<int> emitted(called.begin(), called.end());
    if (rec.ref_trav_idx >= 0) {
        emitted.insert(rec.ref_trav_idx);
    }
    double best_other = -numeric_limits<double>::infinity();
    for (const auto& entry : info->genotype_lls) {
        if (entry.first == called) {
            continue;
        }
        bool all_emitted = true;
        for (int a : entry.first) {
            if (emitted.count(a) == 0) {
                all_emitted = false;
                break;
            }
        }
        if (all_emitted) {
            best_other = max(best_other, entry.second);
        }
    }
    if (!std::isfinite(best_other)) {
        return blank;
    }
    // Nats to phred, matching the VCF's GL, which is log10.
    const double margin_phred = 10.0 * (mine->second - best_other) / log(10.0);
    return min(1.0, max(-1.0, margin_phred / achievable_phred * direct.explained_share));
}

void FlowCaller::collect_anchors_for_record(const StagedSite& rec,
                                            const vector<int>& genotype) {
    if (!anchor_collector.is_enabled() || rec.call_info == nullptr) {
        return;
    }
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
    if (info == nullptr || info->anchor_evidence == nullptr) {
        // A genotype derived from a parent rather than scored here, or a run whose caller is not the
        // read-likelihood one. There are no per-read responsibilities to partition on.
        return;
    }
    anchor_collector.collect(*info->anchor_evidence, info->explained_share,
                             phase_table.phase_ordered_genotype(rec.record_key, genotype),
                             phase_table.haploid_slot(rec.record_key, genotype),
                             print_snarl(rec.snarl),
                             anchor_collector.wants_leaf_test() ? snarl_is_leaf(rec.snarl) : true,
                             anchor_gqn_for(rec, genotype), read_strands, rec.record_key);
}

void FlowCaller::hand_off_deferred_records() {
    if (!staged_sites.active()) {
        return;
    }
    const vector<StagedSite>& pending = staged_sites.nested();
    // A chain that gets no line still gets anchors: one whose variation an enclosing block's ALT
    // already spells, and one with no reference path, so no REF or POS, whose anchors are placed by
    // node ID. Each chain's anchors are its own and the anchor writer sorts every anchor before
    // writing, so they are collected on several threads, before the chains are handed over.
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < pending.size(); ++i) {
        const StagedSite& pr = pending[i];
        if (!pr.dropped && (pr.reported_inline || pr.no_reference)) {
            collect_anchors_for_record(pr, linker.chosen_genotype(pr));
        }
    }
    // Hand every surviving chain to the render, so that nested and top-level records are written
    // in one place from their chosen genotypes. A dropped chain is not handed over, since the
    // sample has no copy of it. A chain an enclosing block's ALT spells, or one with no reference
    // path, gets no line, since the render calls emit_variant for every record it is given; its
    // anchors were collected above.
    const StagedSiteTable::HandOff held = staged_sites.hand_off();
    if (held.inline_unrendered > 0 && show_progress) {
        cerr << "[vg call] block emission: " << held.inline_unrendered
             << " chains genotyped, recorded and phased, but left unrendered because an enclosing"
             << " block's ALT already spells them out" << endl;
    }
    if (held.no_ref_unrendered > 0 && show_progress) {
        cerr << "[vg call] off-reference nested: " << held.no_ref_unrendered
             << " chains chosen and left unrendered, having no reference position to write" << endl;
    }
}

}
