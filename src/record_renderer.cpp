#include <cmath>
#include <limits>
#include <set>

#include <omp.h>

#include "record_renderer.hpp"
#include "read_likelihood_caller.hpp"

namespace vg::multipass {

void RecordRenderer::configure(SiteReader site_reader) {
    this->reader = std::move(site_reader);
}

double RecordRenderer::anchor_gqn(const StagedSite& rec, const vector<int>& chosen,
                                  const LinkageCollector* model) {
    const SiteScore* score = rec.score;
    // The value of the direct pass, with -1, which means there was no gap to normalise, turned
    // into NaN, so that it stays apart from the signed range [-1, 1].
    const double direct_value = (score == nullptr || score->gq_fraction < 0.0)
        ? std::numeric_limits<double>::quiet_NaN()
        : score->gq_fraction;
    const double blank = std::numeric_limits<double>::quiet_NaN();
    if (model == nullptr) {
        return direct_value;
    }
    const auto& moved = model->moved_quality();
    const auto found = moved.find(rec.record_key);
    if (found == moved.end()) {
        // The linkage model left the call alone, so the value of the direct pass still holds.
        return direct_value;
    }
    // The model changed the call, so the value of the direct pass describes the wrong genotype.
    // Any failure below gives NaN.
    if (score == nullptr || score->genotype_lls.empty()) {
        return blank;
    }
    // Scale the margin as GQN in the VCF does, so that the two agree.
    const LinkageCollector::DirectQuality& direct = found->second.direct;
    if (!(direct.achievable_gap > 0.0)) {
        return blank;
    }
    const double achievable_phred = 10.0 * direct.achievable_gap / log(10.0);

    // Use the chosen genotype. `rec.genotype` is the call of the direct pass, which the reads
    // prefer, so its margin would have the wrong sign.
    vector<int> called = chosen;
    std::ranges::sort(called);
    const auto mine = score->genotype_lls.find(called);
    if (mine == score->genotype_lls.end()) {
        return blank;
    }
    // Compare only with genotypes over the written alleles, as GL in the VCF does. Those are the
    // reference traversal and the traversals the chosen genotype names.
    set<int> emitted(called.begin(), called.end());
    if (rec.ref_trav_idx >= 0) {
        emitted.insert(rec.ref_trav_idx);
    }
    double best_other = -numeric_limits<double>::infinity();
    for (const auto& [other, ll] : score->genotype_lls) {
        if (other == called) {
            continue;
        }
        bool all_emitted = true;
        for (int a : other) {
            if (!emitted.contains(a)) {
                all_emitted = false;
                break;
            }
        }
        if (all_emitted) {
            best_other = max(best_other, ll);
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
    const SiteScore* score = rec.score;
    if (score == nullptr || score->anchor_evidence == nullptr) {
        // The genotype came from a parent, or the run uses another caller, so there are no
        // per-read responsibilities to collect anchors from.
        return;
    }
    anchors.collect(*score->anchor_evidence, score->explained_share,
                    phases.phase_ordered_genotype(rec.record_key, genotype),
                    phases.haploid_slot(rec.record_key, genotype), rec.id,
                    anchors.wants_leaf_test() ? rec.leaf : true,
                    anchor_gqn(rec, genotype, model), strands, rec.record_key);
}

void RecordRenderer::hand_off(StagedSiteTable& staged, const PhaseTable& phases,
                              const ReadStrandTable& strands, const GenotypeLinker& linker,
                              AnchorCollector* anchors, bool show_progress) const {
    if (!staged.active()) {
        return;
    }
    const vector<StagedSite>& pending = staged.nested();
    // Collect anchors for the chains that get no line. Those are the chains whose variation the ALT
    // of an enclosing block already spells, and the chains with no reference path, whose anchors
    // are placed by node ID. The anchor writer sorts every anchor before writing, so the chains can
    // be done on several threads.
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < pending.size(); ++i) {
        const StagedSite& pr = pending[i];
        if (anchors != nullptr && !pr.dropped && (pr.reported_inline || pr.no_reference)) {
            collect_anchors(pr, linker.chosen_genotype(pr), phases, strands, linker.collector(),
                            *anchors);
        }
    }
    // Hand the other chains to the render queues, so that nested and top-level records are
    // written in one place from their chosen genotypes.
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
                            bool show_progress, const function<void()>& before_lines) const {
    hand_off(staged, phases, strands, linker, anchors, show_progress);
    if (!staged.active()) {
        return;
    }
    if (before_lines) {
        before_lines();
    }
    const size_t n_threads = staged.queue_count();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t t = 0; t < n_threads; ++t) {
        for (const StagedSite& rec : staged.queue(t)) {
            // Write the line from the genotype the linkage model chose, so that the ALT list, QUAL,
            // AD, GL and GQI agree with the call. It stays sorted, since the line is built from its
            // order, and the anchors put it in phase order themselves.
            vector<int> genotype = linker.chosen_genotype(rec);
            // Collect the anchors before writing the line, which hands the site's CallInfo on to
            // update_vcf_info.
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
