#include <chrono>
#include <functional>

#include <omp.h>

#include "genotype_linker.hpp"
#include "read_likelihood_caller.hpp"
#include "vcf_record.hpp"

namespace vg::multipass {

void GenotypeLinker::configure(LinkageCollector* collector, const PanelLookup* panel) {
    model = collector;
    lookup = panel;
}

void GenotypeLinker::set_site_reader(SiteReader site_reader) {
    reader = std::move(site_reader);
}

vector<int> GenotypeLinker::chosen_genotype(const StagedSite& rec) const {
    vector<int> genotype = rec.genotype;
    int chosen_a = -1;
    int chosen_b = -1;
    if (size_t chosen_ploidy = 0;
        model != nullptr
        && model->chosen_traversals(rec.record_key, &chosen_a, &chosen_b, &chosen_ploidy)
        && chosen_ploidy == genotype.size()) {
        genotype.assign(1, chosen_a);
        if (chosen_ploidy > 1) {
            genotype.push_back(chosen_b);
        }
    }
    return genotype;
}

/// Get the quality inputs of a direct call, which the model keeps for rewriting the record if
/// it moves the genotype.
///
/// `genotyper` gives the GQ factor, since the share and depth discounts are its settings.
static LinkageCollector::DirectQuality direct_quality_of(const ReadLikelihoodSnarlCaller* genotyper,
                                                         const SiteScore& score) {
    return LinkageCollector::DirectQuality{
        .explained_share = score.explained_share,
        .gq_factor = genotyper != nullptr ? genotyper->gq_factor(score) : score.explained_share,
        .achievable_gap = score.achievable_gap,
    };
}

bool GenotypeLinker::file(StagedSite& site, bool nested) const {
    if (model == nullptr || site.score == nullptr) {
        return false;
    }
    // File only a genotype of one or two traversals, none of them missing or a star allele.
    const vector<int>& genotype = site.genotype;
    const size_t site_ploidy = genotype.size();
    if (site_ploidy != 1 && site_ploidy != 2) {
        return false;
    }
    if (std::ranges::any_of(genotype, [](int allele) { return allele < 0; })) {
        return false;
    }
    const NestingPlacement& placement = site.placement;
    // The written alleles are chosen when the record is built, and `set_allele_map` supplies the
    // allele map then.
    static const vector<int> no_allele_map;
    model->record(
        site.reference.locus(*reader.graph, *reader.sites, site.net), site.score->genotype_lls,
        site.panel_alleles(*lookup),
        LinkageCollector::DirectCall{
            .called_i = genotype[0],
            .called_j = site_ploidy > 1 ? genotype[1] : genotype[0],
            .ploidy = site_ploidy,
            .quality = direct_quality_of(reader.genotyper, *site.score),
        },
        no_allele_map, site.record_key,
        LinkageCollector::SiteContext{
            .start_node = (int64_t)reader.graph->get_id(reader.sites->get_site_start(site.net)),
            .end_node = (int64_t)reader.graph->get_id(reader.sites->get_site_end(site.net)),
            .nested = nested,
            .parent_record_key = placement.parent_record_key,
            .parent_crossing = placement.parent_crossing,
            .level = placement.level,
            .emitted = false,
            .unpositioned = site.reference.no_reference,
            .chain_key = placement.chain_key,
            .freq_prior =
                freq_prior(site.alleles),
        });
    return true;
}

double GenotypeLinker::freq_prior(const CandidateAlleles& candidates) const {
    const LinkageModel::Params& params = model->model_params();
    if (params.hp_prior <= 0.0) {
        return -1.0;
    }
    vector<string> alleles;
    alleles.reserve(candidates.walks.size());
    for (size_t a = 0; a < candidates.walks.size(); ++a) {
        alleles.push_back(candidates.sequence(*reader.graph, a));
    }
    const size_t ref = candidates.ref_allele >= 0 ? (size_t)candidates.ref_allele : (size_t)-1;
    return LinkageModel::run_length_site(alleles, params.hp_prior_run, ref) ? params.hp_prior : -1.0;
}

void GenotypeLinker::resync(StagedSiteTable& sites) const {
    if (model == nullptr) {
        return;
    }
    const vector<StagedSite*> records = sites.in_order();
    // Fill each record's cache of panel alleles in parallel first. Those GBWT lookups would take
    // most of the time of the serial loop.
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < records.size(); ++i) {
        records[i]->panel_alleles(*lookup);
    }
    size_t rescored = 0;
    size_t refused = 0;
    for (StagedSite* recp : records) {
        StagedSite& rec = *recp;
        const SiteScore* score = rec.score;
        if (score == nullptr || score->genotype_lls.empty()) {
            continue;
        }
        // The corrected best genotype becomes the called pair of the entry.
        const vector<int>* best = nullptr;
        double best_ll = -numeric_limits<double>::infinity();
        for (const auto& [genotype, ll] : score->genotype_lls) {
            if (ll > best_ll) {
                best_ll = ll;
                best = &genotype;
            }
        }
        if (best == nullptr || best->empty()) {
            continue;
        }
        const int called_i = (*best)[0];
        const int called_j = best->size() > 1 ? (*best)[1] : called_i;
        if (model->rescore(rec.record_key, score->genotype_lls, rec.panel_alleles(*lookup),
                           called_i, called_j)) {
            ++rescored;
        } else {
            // The key has no live entry, as for a chain this round has not reinstated, or the
            // corrected likelihoods cannot be compacted.
            ++refused;
        }
    }
    cerr << "[vg call] re-genotyping: " << rescored << " sites re-scored into the layer, "
         << refused << " refused for want of a live entry or a compactable space" << endl;
}

GenotypeLinker::PassCounts GenotypeLinker::link(StagedSiteTable& sites, PhaseTable& phases,
                                                bool keep_phase) {
    PassCounts counts;
    if (model != nullptr) {
        counts.levels = model->max_level();
    }
    // Gather the nested sites on the first pass, and keep them for the passes re-genotyping runs.
    vector<StagedSite>& pending = sites.gather_nested();

    // A later pass derives everything again from the chosen genotypes. Clear `dropped`, since a
    // correction can move a parent onto an allele that crosses a dropped chain. The levels then
    // drop again the chains still not crossed.
    if (passes_run > 0) {
        for (StagedSite& pr : pending) {
            pr.dropped = false;
        }
        phases.calls().clear();
        total_moved = 0;
    }
    ++passes_run;

    unordered_map<size_t, StagedSite*> record_by_key = sites.by_key();
    // The crossing finder of each parent traversal, keyed by address. The addresses stay valid
    // for the pass, and the traversals do not change after the direct pass.
    unordered_map<const Traversal*, ChildCrossingFinder> crossing_finders;

    // Choose each level, then revise the chains of the next level from the parents just chosen.
    // Revising can add a chain at a deeper level, so `counts.levels` can grow.
    size_t level = 0;
    while (level <= counts.levels) {
        // The last level builds the phase calls and the mosaic from everything chosen. If a
        // deeper chain appears, a later level builds them again.
        resolve_level(level, level == counts.levels, keep_phase ? &phases.calls() : nullptr);

        const unordered_map<size_t, const PhaseCall*> chosen =
            parent_phase_calls(pending, phases, level + 1);
        for (size_t i = 0; i < pending.size(); ++i) {
            const StagedSite& pr = pending[i];
            if (pr.placement.level == level + 1 && !pr.dropped) {
                revise_chain(sites, i, record_by_key, chosen, crossing_finders, counts);
            }
        }
        if (model != nullptr) {
            counts.levels = max(counts.levels, model->max_level());
        }
        ++level;
    }

    return counts;
}

unordered_map<size_t, const GenotypeLinker::PhaseCall*> GenotypeLinker::parent_phase_calls(
    const vector<StagedSite>& pending, const PhaseTable& phases, size_t level) {
    unordered_set<size_t> parents;
    for (const StagedSite& pr : pending) {
        if (pr.placement.level == level) {
            parents.insert(pr.placement.parent_record_key);
        }
    }
    unordered_map<size_t, const PhaseCall*> chosen;
    chosen.reserve(parents.size() * 2);
    for (const PhaseCall& pc : phases.calls()) {
        if (parents.contains(pc.record_key)) {
            chosen[pc.record_key] = &pc;
        }
    }
    return chosen;
}

void GenotypeLinker::revise_chain(
    StagedSiteTable& sites, size_t index, const unordered_map<size_t, StagedSite*>& record_by_key,
    const unordered_map<size_t, const PhaseCall*>& chosen,
    unordered_map<const Traversal*, ChildCrossingFinder>& crossing_finders, PassCounts& counts) {
    StagedSite& pr = sites.nested()[index];
    const auto parent_record = record_by_key.find(pr.placement.parent_record_key);
    const StagedSite* parent =
        parent_record != record_by_key.end() ? parent_record->second : nullptr;
    if (model != nullptr && parent != nullptr) {
        // The parent was placed when its own level was revised.
        place_chain(pr, *parent, crossing_finders);
    }
    if (!pr.placement.crossing_known) {
        // The parent has more candidate traversals than a 64-bit crossing mask can hold. Count the
        // chain and leave it as it is, since its mask does not mean that no allele crosses it.
        ++counts.crossing_unknown;
        return;
    }
    if (pr.placement.parent_crossing == 0) {
        // No candidate traversal of the parent crosses the chain, so the sample has no copy.
        ++counts.no_crossing;
        counts.retracted += drop_subtree(sites, index);
        return;
    }
    int chosen_first;
    int chosen_second;
    if (!chosen_pair_of_parent(pr, parent, chosen, chosen_first, chosen_second)) {
        // The chain keeps the ploidy it was called at.
        ++counts.no_chosen;
        return;
    }
    const LinkageCollector::Relation rel = LinkageCollector::relate_to_parent(
        pr.placement.parent_crossing, chosen_first, chosen_second);
    const auto copies = (int)rel.copies;
    if (copies == 0) {
        counts.retracted += drop_subtree(sites, index);
        return;
    }
    if (copies == pr.ploidy && model != nullptr && model->has_entry(pr.record_key)) {
        // The chain is filed at the ploidy its parent implies. A chain with no entry still has to
        // be filed: the direct pass stages a chain that no called parent allele crosses, but does
        // not file it.
        return;
    }
    if (!buildable(pr)) {
        ++counts.unrenderable;
        return;
    }
    if (!adopt_ploidy(pr, copies, counts)) {
        return;
    }
    const bool had_entry = model != nullptr && model->has_entry(pr.record_key);
    if (pr.score != nullptr) {
        // Retract before filing, because `live_index` finds the first live entry for a key.
        if (had_entry) {
            model->retract(pr.record_key);
        }
        file(pr, pr.ploidy == 1);
        if (had_entry && !model->has_entry(pr.record_key)) {
            // `record` files nothing for a site whose compact allele space it cannot describe.
            // The old entry stays retracted, and the record keeps its own call.
            ++counts.unspecifiable;
        }
    }
    if (had_entry) {
        ++counts.revised;
    } else {
        ++counts.gained;
    }
    // The chosen pair of this chain has changed, or the chain is new, so the crossing masks of
    // its children are out of date.
    recompute_child_crossings(sites, pr);
}

void GenotypeLinker::place_chain(
    StagedSite& pr, const StagedSite& parent,
    unordered_map<const Traversal*, ChildCrossingFinder>& crossing_finders) const {
    const size_t offset =
        parent.placement.parent_offset
        + entry_bases_along_genotype(*reader.sites, *reader.graph, parent.alleles.walks,
                                     chosen_genotype(parent), pr.net, crossing_finders);
    if (offset == pr.placement.parent_offset) {
        return;
    }
    if (pr.reference.no_reference) {
        pr.reference.position_from_parent += (int64_t)offset - (int64_t)pr.placement.parent_offset;
        model->set_position(pr.record_key,
                            pr.reference.locus(*reader.graph, *reader.sites, pr.net).position);
    }
    pr.placement.parent_offset = offset;
}

bool GenotypeLinker::chosen_pair_of_parent(const StagedSite& pr, const StagedSite* parent,
                                           const unordered_map<size_t, const PhaseCall*>& chosen,
                                           int& first, int& second) const {
    first = -1;
    second = -1;
    if (const auto found = chosen.find(pr.placement.parent_record_key); found != chosen.end()) {
        const PhaseCall& call = *found->second;
        first = call.trav_first;
        second = call.ploidy == 2 ? call.trav_second : -1;
        return true;
    }
    if (model != nullptr && parent != nullptr && !parent->genotype.empty()) {
        // The model gave the parent no phase call, so the parent is rendered at the genotype
        // `chosen_genotype()` reads.
        const vector<int> parent_genotype = chosen_genotype(*parent);
        first = parent_genotype[0];
        second = parent_genotype.size() > 1 ? parent_genotype[1] : -1;
        if (first < 0) {
            std::swap(first, second);
        }
        return first >= 0;
    }
    return false;
}

bool GenotypeLinker::buildable(const StagedSite& pr) {
    if (pr.alleles.walks.empty()) {
        return false;
    }
    // `SiteEmitter::emit` spells the reference allele without checking its number. A record with no
    // reference path is never emitted.
    if (!pr.reference.no_reference
        && (pr.alleles.ref_allele < 0
            || (size_t)pr.alleles.ref_allele >= pr.alleles.walks.size())) {
        return false;
    }
    if (pr.genotype.empty()) {
        return false;
    }
    return std::ranges::none_of(pr.genotype, [&](int allele) {
        return allele >= 0 && (size_t)allele >= pr.alleles.walks.size();
    });
}

bool GenotypeLinker::adopt_ploidy(StagedSite& pr, int copies, PassCounts& counts) {
    if (copies == pr.ploidy) {
        return true;
    }
    SiteScore* rl = pr.score;
    if (rl == nullptr || rl->alt_ploidy_info == nullptr || rl->alt_ploidy_info->ploidy != copies
        || (int)rl->alt_ploidy_best.size() != copies) {
        // The direct pass scored no answer at that ploidy, because the chain offers too few
        // traversals for a second genotype. The chain keeps a ploidy its chosen parent
        // contradicts.
        ++counts.ploidy_unscored;
        return false;
    }
    if (!std::ranges::all_of(rl->alt_ploidy_best, [&](int allele) {
            return allele >= 0 && (size_t)allele < pr.alleles.walks.size();
        })) {
        return false;
    }
    // Exchange the two answers, so that a chain can follow its parent to either ploidy however
    // often the linkage pass runs.
    vector<int> promoted_genotype = rl->alt_ploidy_best;
    const vector<int> demoted_genotype = pr.genotype;
    // `demoted` now owns what `rl` points to.
    unique_ptr<SiteScore> demoted = pr.take_score();
    unique_ptr<SiteScore> promoted = std::move(rl->alt_ploidy_info);
    // The answer at the other ploidy does not copy the evidence, which does not depend on
    // ploidy, so the evidence moves to whichever answer is in front.
    promoted->anchor_evidence = std::move(rl->anchor_evidence);
    promoted->phase_evidence = std::move(rl->phase_evidence);
    promoted->alt_ploidy_best = demoted_genotype;
    promoted->alt_ploidy_info = std::move(demoted);
    pr.genotype = std::move(promoted_genotype);
    pr.ploidy = copies;
    pr.set_score(std::move(promoted));
    return true;
}

void GenotypeLinker::recompute_child_crossings(StagedSiteTable& sites,
                                               const StagedSite& pr) const {
    const vector<size_t>* kids = sites.children_of(pr.record_key);
    if (kids == nullptr) {
        return;
    }
    vector<ChildCrossingFinder> pr_walks;
    pr_walks.reserve(pr.alleles.walks.size());
    for (const Traversal& t : pr.alleles.walks) {
        pr_walks.emplace_back(*reader.sites, *reader.graph, t);
    }
    vector<StagedSite>& pending = sites.nested();
    for (size_t ci : *kids) {
        StagedSite& child = pending[ci];
        bool known = true;
        child.placement.parent_crossing = crossing_mask(pr_walks, child.net, &known);
        child.placement.crossing_known = known;
    }
}

size_t GenotypeLinker::drop_subtree(StagedSiteTable& sites, size_t root) {
    vector<StagedSite>& pending = sites.nested();
    size_t dropped_here = 0;
    // Use an explicit stack, since the depth of the snarl tree depends on the data.
    vector<size_t> stack{root};
    while (!stack.empty()) {
        size_t idx = stack.back();
        stack.pop_back();
        StagedSite& victim = pending[idx];
        if (victim.dropped) {
            continue;
        }
        victim.dropped = true;
        if (model != nullptr && model->retract(victim.record_key)) {
            ++dropped_here;
        }
        if (const vector<size_t>* kids = sites.children_of(victim.record_key)) {
            for (size_t k : *kids) {
                if (k != idx) {
                    stack.push_back(k);
                }
            }
        }
    }
    return dropped_here;
}

void GenotypeLinker::resolve_level(size_t level, bool last, vector<PhaseCall>* calls) {
    if (model == nullptr) {
        return;
    }
    auto start = std::chrono::steady_clock::now();
    const size_t moved = model->resolve_level(level, last, calls);
    double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    total_seconds += seconds;
    total_moved += moved;
    if (!last) {
        cerr << "[vg call] linkage level " << level << ": "
             << model->num_sites_at(level) << " sites, "
             << moved << " genotypes moved by linkage, " << seconds << " s" << endl;
    }
}

void GenotypeLinker::report() const {
    if (model == nullptr) {
        return;
    }
    cerr << "[vg call] linkage: " << model->num_sites() << " sites, "
         << (model->bytes() / (1024.0 * 1024.0)) << " MB retained, "
         << total_moved << " genotypes moved by linkage, " << total_seconds << " s" << endl;
    if (model->num_duplicate_live_keys() > 0) {
        // `retract` retracts only the first live entry of a key, so it cannot revise these sites.
        cerr << "[vg call] linkage: " << model->num_duplicate_live_keys()
             << " sites recorded onto a key that already had a live entry; the retract path cannot"
             << " address these" << endl;
    }
    if (model->model_params().hp_prior > 0.0) {
        cerr << "[vg call] linkage: " << model->num_site_prior_entries()
             << " live entries decoded at a run-length site's own frequency exponent (--hp-prior)"
             << endl;
    }
    cerr << "[vg call] linkage: " << model->mutated_strands()
         << " phased strands carry an allele other than their panel haplotype's, explained as"
         << " a mutation (--linkage-mutation)" << endl;
}

}
