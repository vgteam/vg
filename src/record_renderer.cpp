#include <cmath>
#include <limits>
#include <set>

#include <omp.h>

#include "record_renderer.hpp"
#include "read_likelihood_caller.hpp"

namespace vg {

void RecordRenderer::configure(SiteReader reader, function<bool(const Snarl&)> is_leaf) {
    this->reader = std::move(reader);
    this->is_leaf = std::move(is_leaf);
}

// The anchor gqn column for a staged site.
//
// `gq_fraction` was computed in the direct pass for the reads' best genotype, so on a record whose
// genotype the linkage model changed, it describes the abandoned genotype. Such a record gets the
// signed margin of its chosen genotype instead, as the VCF's GQN does (see
// ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype).
//
// Returns the direct pass's value when the model did not change the call; the recomputed signed margin
// when it did; and NaN, written as ".", when it did but the margin cannot be recomputed.
double RecordRenderer::anchor_gqn(const StagedSite& rec, const vector<int>& chosen,
                                  const LinkageCollector* model) {
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
    // The direct pass's value, with its "no gap to normalise" value (-1) turned into NaN, written as ".",
    // so that it stays distinct from the signed range [-1, 1].
    const double direct_value = (info == nullptr || info->gq_fraction < 0.0)
        ? std::numeric_limits<double>::quiet_NaN()
        : info->gq_fraction;
    const double blank = std::numeric_limits<double>::quiet_NaN();
    if (model == nullptr) {
        return direct_value;
    }
    const auto& moved = model->moved_quality();
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

void RecordRenderer::collect_anchors(const StagedSite& rec, const vector<int>& genotype,
                                     const PhaseTable& phases, const ReadStrandTable& strands,
                                     const LinkageCollector* model,
                                     AnchorCollector& anchors) const {
    if (rec.call_info == nullptr) {
        return;
    }
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(rec.call_info.get());
    if (info == nullptr || info->anchor_evidence == nullptr) {
        // A genotype derived from a parent rather than scored here, or a run whose caller is not the
        // read-likelihood one. There are no per-read responsibilities to partition on.
        return;
    }
    anchors.collect(*info->anchor_evidence, info->explained_share,
                    phases.phase_ordered_genotype(rec.record_key, genotype),
                    phases.haploid_slot(rec.record_key, genotype), reader.name(rec.snarl),
                    anchors.wants_leaf_test() ? is_leaf(rec.snarl) : true,
                    anchor_gqn(rec, genotype, model), strands, rec.record_key);
}

void RecordRenderer::hand_off(StagedSiteTable& staged, const PhaseTable& phases,
                              const ReadStrandTable& strands, const GenotypeLinker& linker,
                              AnchorCollector* anchors, bool show_progress) const {
    if (!staged.active()) {
        return;
    }
    const vector<StagedSite>& pending = staged.nested();
    // A chain that gets no line still gets anchors: one whose variation an enclosing block's ALT
    // already spells, and one with no reference path, so no REF or POS, whose anchors are placed by
    // node ID. Each chain's anchors are its own and the anchor writer sorts every anchor before
    // writing, so they are collected on several threads, before the chains are handed over.
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < pending.size(); ++i) {
        const StagedSite& pr = pending[i];
        if (anchors != nullptr && !pr.dropped && (pr.reported_inline || pr.no_reference)) {
            collect_anchors(pr, linker.chosen_genotype(pr), phases, strands, linker.collector(),
                            *anchors);
        }
    }
    // Hand every surviving chain to the render, so that nested and top-level records are written
    // in one place from their chosen genotypes. A dropped chain is not handed over, since the
    // sample has no copy of it. A chain an enclosing block's ALT spells, or one with no reference
    // path, gets no line, since the render writes a line for every record it is given; its anchors
    // were collected above.
    const StagedSiteTable::HandOff held = staged.hand_off();
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

void RecordRenderer::render(StagedSiteTable& staged, const PhaseTable& phases,
                            const ReadStrandTable& strands, const GenotypeLinker& linker,
                            AnchorCollector* anchors, const LineWriter& write_line,
                            bool show_progress) const {
    // Every linkage pass is done, so the records move to the render, once, which also keeps their
    // anchors from being collected twice.
    hand_off(staged, phases, strands, linker, anchors, show_progress);
    if (!staged.active()) {
        return;
    }
    // The records are nested chains as well as top-level sites, and each carries its own nesting
    // in its `StagedSite`.
    const size_t n_threads = staged.queue_count();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t t = 0; t < n_threads; ++t) {
        for (StagedSite& rec : staged.queue(t)) {
            // The chosen pair, not the direct pass's. The ALT list, whether a line is written at
            // all, QUAL, and the arity of AD, GL and GQI are all built from the genotype passed in,
            // so they agree with the call.
            vector<int> genotype = linker.chosen_genotype(rec);
            // Before the line, whose writer passes the CallInfo on to update_vcf_info. The anchors
            // are collected in phase order, while `genotype` itself stays sorted, since the line's
            // ALT list, AD, GL and QUAL are built from its order.
            if (anchors != nullptr) {
                collect_anchors(rec, genotype, phases, strands, linker.collector(), *anchors);
            }
            write_line(rec, genotype);
        }
    }
    if (show_progress) {
        cerr << "[vg call] rendered " << staged.queued_count()
             << " retained records after the direct pass" << endl;
    }
}

}
