#include <chrono>
#include <functional>

#include <omp.h>

#include "genotype_linker.hpp"
#include "read_likelihood_caller.hpp"
#include "vcf_record.hpp"

namespace vg {

void GenotypeLinker::configure(LinkageCollector* collector, const PanelLookup* panel) {
    this->model = collector;
    this->lookup = panel;
}

void GenotypeLinker::set_site_reader(SiteReader reader) {
    this->reader = std::move(reader);
}

// Used by both anchor-collection paths, the render and the hand-off, so that records with no VCF
// line (`reported_inline` and `no_reference`) also get anchors for their chosen genotype.
vector<int> GenotypeLinker::chosen_genotype(const StagedSite& rec) const {
    vector<int> genotype = rec.genotype;
    int chosen_a = -1, chosen_b = -1;
    size_t chosen_ploidy = 0;
    if (model != nullptr
        && model->chosen_traversals(rec.record_key, &chosen_a, &chosen_b, &chosen_ploidy)
        && chosen_ploidy == genotype.size()) {
        genotype.assign(1, chosen_a);
        if (chosen_ploidy > 1) {
            genotype.push_back(chosen_b);
        }
    }
    return genotype;
}

/// The locus of a (path name, position) pair: the contig as the VCF names it, and the position,
/// held at 0 or above.
static SiteLocus locus_of(pair<string, int64_t> pos_info) {
    const string locus = PathMetadata::parse_locus_name(pos_info.first);
    if (locus != PathMetadata::NO_LOCUS_NAME) {
        pos_info.first = locus;
    }
    return SiteLocus{pos_info.first, (size_t)max((int64_t)0, pos_info.second)};
}

SiteLocus GenotypeLinker::site_locus(const SiteBounds& site, const string& ref_path_name,
                                     int ref_offset) const {
    // The position before the record's alleles are trimmed, which can move POS.
    // `get_ref_position` names the base path, as in "CHM13#0#chr20".
    return locus_of(get_ref_position(*reader.graph, site.start, site.end, ref_path_name,
                                     ref_offset));
}

SiteLocus GenotypeLinker::off_reference_site_locus(const string& ref_path_name,
                                                   int64_t stand_in_position) {
    // Not `get_ref_position`: `get_ref_interval` asserts on a snarl the path does not pass
    // through.
    return locus_of(make_pair(ref_path_name, stand_in_position));
}

/// The quality inputs of the direct call `info`, which the linkage collector keeps for rewriting
/// the record if the model moves it. `genotyper` decides the GQ factor, since --no-share-quality
/// and --depth-quality are its settings.
static LinkageCollector::DirectQuality direct_quality_of(const SiteGenotyper* genotyper,
                                                         const SiteScore& info) {
    return LinkageCollector::DirectQuality{
        .explained_share = info.explained_share,
        .gq_factor = genotyper != nullptr ? genotyper->gq_factor(info) : info.explained_share,
        .achievable_gap = info.achievable_gap,
    };
}

bool GenotypeLinker::add(const SiteBounds& site, const vector<Traversal>& travs,
                         const vector<int>& trav_genotype, const SiteScore* score,
                         int ref_trav_idx, const string& ref_path_name, int ref_offset,
                         size_t record_key, const NestingPlacement& placement,
                         bool no_reference, int64_t position_from_parent,
                         vector<int>* panel_out) const {
    if (model == nullptr) {
        return false;
    }
    const SiteScore* rl_info = score;
    if (rl_info == nullptr) {
        return false;
    }
    // The same test the emitter uses, in traversal space: a genotype of one or two alleles, none of
    // them a missing or star marker. Haploid chains are included.
    const size_t site_ploidy = trav_genotype.size();
    if (site_ploidy != 1 && site_ploidy != 2) {
        return false;
    }
    for (int allele : trav_genotype) {
        if (allele < 0) {
            return false;
        }
    }
    const SiteLocus locus = no_reference
                                ? off_reference_site_locus(ref_path_name, position_from_parent)
                                : site_locus(site, ref_path_name, ref_offset);
    const int called_i = trav_genotype[0];
    const int called_j = site_ploidy > 1 ? trav_genotype[1] : called_i;
    // No allele map yet: the written alleles are chosen when the record is built, and
    // `set_allele_map` supplies the map then.
    static const vector<int> no_allele_map;
    vector<int> panel = lookup->alleles(travs);
    model->record(
        locus.contig, locus.position,
        rl_info->genotype_lls,
        panel,
        called_i, called_j, no_allele_map,
        record_key,
        direct_quality_of(reader.genotyper, *rl_info), site_ploidy,
        (int64_t)reader.graph->get_id(site.start), (int64_t)reader.graph->get_id(site.end),
        // `nested` only when one copy of the chain is present, as for any other chain; a chain with
        // two copies joins its parent's diploid group.
        LinkageCollector::SiteContext{
            .nested = placement.one_copy,
            .parent_record_key = placement.parent_record_key,
            .parent_crossing = placement.parent_crossing,
            .level = placement.level,
            .emitted = false,
            .unpositioned = no_reference,
            .chain_key = placement.chain_key,
            .freq_prior = freq_prior(travs, ref_trav_idx),
        });
    if (panel_out != nullptr) {
        *panel_out = std::move(panel);
    }
    return true;
}

double GenotypeLinker::freq_prior(const vector<Traversal>& travs, int ref_trav_idx) const {
    const LinkageModel::Params& params = model->model_params();
    if (params.hp_prior <= 0.0) {
        return -1.0;
    }
    vector<string> alleles;
    alleles.reserve(travs.size());
    for (const Traversal& walk : travs) {
        alleles.push_back(reader.spell(walk));
    }
    const size_t ref = ref_trav_idx >= 0 ? (size_t)ref_trav_idx : (size_t)-1;
    return LinkageModel::run_length_site(alleles, params.hp_prior_run, ref) ? params.hp_prior : -1.0;
}

void GenotypeLinker::resync(StagedSiteTable& sites) const {
    if (model == nullptr) {
        return;
    }
    const vector<StagedSite*> records = sites.in_order();
    // The loop below is serial, and most of its time would go to each record's first
    // `PanelLookup::alleles`, a GBWT lookup per allele. Each record's lookup is independent of the
    // others', so fill the caches in parallel first.
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < records.size(); ++i) {
        records[i]->panel_alleles(*lookup);
    }
    size_t rescored = 0, refused = 0;
    for (StagedSite* recp : records) {
        StagedSite& rec = *recp;
        const SiteScore* info = rec.score;
        if (info == nullptr || info->genotype_lls.empty()) {
            continue;
        }
        // The corrected best genotype, which becomes the entry's called pair.
        const vector<int>* best = nullptr;
        double best_ll = -numeric_limits<double>::infinity();
        for (const auto& kv : info->genotype_lls) {
            if (kv.second > best_ll) {
                best_ll = kv.second;
                best = &kv.first;
            }
        }
        if (best == nullptr || best->empty()) {
            continue;
        }
        const int called_i = (*best)[0];
        const int called_j = best->size() > 1 ? (*best)[1] : called_i;
        if (model->rescore(rec.record_key, info->genotype_lls, rec.panel_alleles(*lookup),
                           called_i, called_j)) {
            ++rescored;
        } else {
            // The key has no live entry, as for a chain this round has not reinstated, or the
            // corrected likelihoods cannot be compacted. A changed allele space is handled by
            // `rescore`.
            ++refused;
        }
    }
    cerr << "[vg call] re-genotyping: " << rescored << " sites re-scored into the layer, "
         << refused << " refused for want of a live entry or a compactable space" << endl;
}

GenotypeLinker::PassCounts GenotypeLinker::link(StagedSiteTable& sites, PhaseTable& phases,
                                                bool keep_phase) {
    // Descent already happened during the direct pass. What is left is to choose the chains' genotypes
    // in the order their ploidies depend on: a level's parents before its children. The reads are not
    // used.
    PassCounts counts;
    size_t& levels = counts.levels;
    if (model != nullptr) {
        levels = model->max_level();
    }
    // Gathered from the per-thread lists on the first pass only, and kept, since re-genotyping
    // runs the linkage pass again.
    vector<StagedSite>& pending = sites.gather_nested();

    // On a later pass, everything a pass concludes is derived again from the chosen genotypes.
    // `dropped` is cleared, since a correction can move a parent onto an allele that crosses a
    // dropped chain; the level loop then drops again the chains still not crossed, and records
    // the others afresh.
    if (passes_run > 0) {
        for (StagedSite& pr : pending) {
            pr.dropped = false;
        }
        // Appended to by every resolve, so cleared here.
        phases.calls().clear();
        // Accumulated by every resolve too.
        total_moved = 0;
    }
    ++passes_run;

    // Counters for the report.
    size_t& revise_unrenderable = counts.unrenderable;
    size_t& pass_no_crossing = counts.no_crossing;
    size_t& pass_no_chosen = counts.no_chosen;
    size_t& pass_ploidy_unscored = counts.ploidy_unscored;
    unordered_map<size_t, StagedSite*> record_by_key = sites.by_key();
    // Each parent traversal's child offsets, for placing its chains. Keyed by address, which stays
    // valid for the pass as record_by_key's do, and the traversals do not change after the direct
    // pass.
    unordered_map<const Traversal*, ChildPlacer::ChildOffsets> child_offsets;

    size_t& revised = counts.revised;
    size_t& retracted = counts.retracted;
    size_t& gained = counts.gained;
    size_t& crossing_unknown = counts.crossing_unknown;
    size_t& unspecifiable = counts.unspecifiable;
    // One linkage pass per level, in order. `levels` is read again after each pass, since
    // a pass can add a chain at a deeper level, which must still be chosen.
    for (size_t gen = 0; gen <= levels; ++gen) {
        // The final pass has last=true and builds the phasing map and the mosaic from everything
        // accumulated. If it adds a deeper chain, the bound grows and a later pass rebuilds them.
        resolve_level(gen, gen == levels, keep_phase ? &phases.calls() : nullptr);

        // This level's parents are chosen, so each chain under one can be given the ploidy
        // its parent's chosen genotype implies before the chain's own level resolves. The
        // direct pass kept the answer at both ploidies, so this is a revision, not a new call.
        // Only the next level's parents are looked up, so only they are indexed; a key's last
        // PhaseCall wins, as it would in an index over every PhaseCall.
        unordered_set<size_t> next_parents;
        for (const StagedSite& pr : pending) {
            if (pr.level == gen + 1) {
                next_parents.insert(pr.parent_record_key);
            }
        }
        unordered_map<size_t, const LinkageCollector::PhaseCall*> chosen;
        chosen.reserve(next_parents.size() * 2);
        for (const LinkageCollector::PhaseCall& pc : phases.calls()) {
            if (next_parents.count(pc.record_key) != 0) {
                chosen[pc.record_key] = &pc;
            }
        }
        for (size_t i = 0; i < pending.size(); ++i) {
            StagedSite& pr = pending[i];
            if (pr.level != gen + 1 || pr.dropped) {
                continue;
            }
            auto parent_record = record_by_key.find(pr.parent_record_key);
            if (model != nullptr && parent_record != record_by_key.end()) {
                // Place the chain along the allele chosen for its parent, before this level's
                // linkage pass orders and spaces its sites by position. The parent's own offset
                // was placed in the previous level's iteration.
                const StagedSite& par = *parent_record->second;
                const size_t offset =
                    par.chain_offset
                    + ChildPlacer::offset_along_genotype(*reader.graph, par.travs,
                                                         chosen_genotype(par), pr.bounds,
                                                         child_offsets);
                if (offset != pr.chain_offset) {
                    if (pr.no_reference) {
                        pr.position_from_parent += (int64_t)offset - (int64_t)pr.chain_offset;
                        model->set_position(
                            pr.record_key,
                            off_reference_site_locus(pr.ref_path_name, pr.position_from_parent)
                                .position);
                    }
                    pr.chain_offset = offset;
                }
            }
            if (!pr.crossing_known) {
                // The direct pass could not compute this chain's crossing mask, because its parent has
                // more candidate traversals than a 64-bit mask can hold. Left as it is and
                // counted, rather than read as "no allele crosses".
                ++crossing_unknown;
                continue;
            }
            if (pr.parent_crossing == 0) {
                // No candidate traversal of the parent crosses the chain, so no chosen genotype
                // can carry it: the sample has no copy, as at ploidy 0 below.
                ++pass_no_crossing;
                retracted += drop_subtree(sites, i);
                continue;
            }
            // The chosen pair as traversals, which the crossing mask is indexed by, through
            // `LinkageCollector::relate_to_parent`, which `resolve_level` also uses to set
            // `nested_strand`.
            int chosen_first = -1, chosen_second = -1;
            bool have_pair = false;
            auto found = chosen.find(pr.parent_record_key);
            if (found != chosen.end()) {
                const LinkageCollector::PhaseCall& parent = *found->second;
                chosen_first = parent.trav_first;
                chosen_second = parent.ploidy == 2 ? parent.trav_second : -1;
                have_pair = true;
            } else if (model != nullptr && parent_record != record_by_key.end()
                       && !parent_record->second->genotype.empty()) {
                // The linkage model gave the parent no PhaseCall, so the parent is rendered at its
                // own chosen genotype, which `chosen_genotype` reads.
                const vector<int> parent_genotype = chosen_genotype(*parent_record->second);
                chosen_first = parent_genotype[0];
                chosen_second = parent_genotype.size() > 1 ? parent_genotype[1] : -1;
                if (chosen_first < 0) {
                    std::swap(chosen_first, chosen_second);
                }
                have_pair = chosen_first >= 0;
            }
            if (!have_pair) {
                // Neither a PhaseCall nor a called allele of the parent can be read, and the chain
                // keeps the ploidy it was called at. Counted.
                ++pass_no_chosen;
                continue;
            }
            const LinkageCollector::Relation rel = LinkageCollector::relate_to_parent(
                pr.parent_crossing, chosen_first, chosen_second);
            int copies = (int)rel.copies;

            // How many copies of the chain the sample has, from the parent's chosen pair. Computed
            // again on every linkage pass, since the ploidy at descent came from the direct pass's
            // genotype.
            if (copies == 0) {
                // The sample has no copy of this chain, and everything inside it is missing too, so
                // the whole subtree is dropped, whether or not it had lines.
                retracted += drop_subtree(sites, i);
                continue;
            }
            if (copies == pr.ploidy && model != nullptr && model->has_entry(pr.record_key)) {
                // The chain was called at the ploidy its parent's chosen genotype implies, so
                // nothing needs revising. `has_entry` matters: a chain that no called parent allele
                // reached in the direct pass is staged but not recorded, and falling through records it.
                continue;
            }

            // Leave alone a record that cannot be built: no traversals, or a genotype out of range.
            // emit_variant indexes `called_traversals[ref_trav_idx]` unchecked, so a record with a
            // reference path also needs a valid ref_trav_idx. A record with no reference path skips
            // that check, since it is never emitted.
            if (pr.travs.empty()
                || (!pr.no_reference
                    && (pr.ref_trav_idx < 0 || (size_t)pr.ref_trav_idx >= pr.travs.size()))) {
                ++revise_unrenderable;
                continue;
            }
            bool genotype_in_range = !pr.genotype.empty();
            for (int allele : pr.genotype) {
                if (allele >= 0 && (size_t)allele >= pr.travs.size()) {
                    genotype_in_range = false;
                }
            }
            if (!genotype_in_range) {
                ++revise_unrenderable;
                continue;
            }

            // Build the record at the ploidy the chosen parent implies, from the answers kept in the
            // direct pass; `alt_ploidy_info` holds the other ploidy's.
            SiteScore* rl = pr.score;
            unique_ptr<SiteScore> use_info;
            vector<int> use_genotype;
            if (copies == pr.ploidy) {
                use_genotype = pr.genotype;
            } else if (rl != nullptr && rl->alt_ploidy_info != nullptr
                       && (int)rl->alt_ploidy_info->ploidy == copies
                       && (int)rl->alt_ploidy_best.size() == copies) {
                bool ok = true;
                for (int allele : rl->alt_ploidy_best) {
                    if (allele < 0 || (size_t)allele >= pr.travs.size()) {
                        ok = false;
                    }
                }
                if (!ok) {
                    continue;
                }
                // The two answers are exchanged, not one discarded, so that a chain can follow its
                // parent to either ploidy however often the linkage pass runs.
                use_genotype = rl->alt_ploidy_best;
                const vector<int> demoted_genotype = pr.genotype;
                unique_ptr<SiteScore> demoted = pr.take_score();
                // `rl` still points at it -- `demoted` owns what `pr.call_info` did.
                unique_ptr<SiteScore> promoted = std::move(rl->alt_ploidy_info);
                // The fields that do not depend on ploidy go with whichever answer is in front, since
                // the alternate does not copy them.
                promoted->anchor_evidence = std::move(rl->anchor_evidence);
                promoted->phase_evidence = std::move(rl->phase_evidence);
                // The replaced answer becomes the new alternate, with the genotype it was called at.
                promoted->alt_ploidy_best = demoted_genotype;
                promoted->alt_ploidy_info = std::move(demoted);
                use_info = std::move(promoted);
            } else {
                // No answer at that ploidy: the direct pass computed none, because the chain offers too
                // few traversals for a second genotype. The chain keeps a ploidy its chosen parent
                // contradicts, which is counted.
                ++pass_ploidy_unscored;
                continue;
            }
            // Whether this chain was already in the linkage model: whether it is being revised or
            // added, and whether there is an old entry to retract.
            const bool had_entry = model != nullptr && model->has_entry(pr.record_key);
            // Revise the staged site; the render builds its line once, at the end, from the chosen
            // genotype.
            pr.genotype = use_genotype;
            pr.ploidy = copies;
            if (use_info != nullptr) {
                pr.set_score(std::move(use_info));
            }
            const SiteScore* used = pr.score;
            if (used != nullptr) {
                // In traversal space, as `add` files it, so the linkage pass and the direct pass
                // describe a site the same way. No allele map yet, as in `add`.
                static const vector<int> no_allele_map;
                const vector<int>& trav_to_allele_vec = no_allele_map;
                const SiteLocus locus =
                    pr.no_reference
                        ? off_reference_site_locus(pr.ref_path_name, pr.position_from_parent)
                        : site_locus(pr.bounds, pr.ref_path_name, pr.ref_offset);
                int called_i = use_genotype.empty() ? -1 : use_genotype[0];
                int called_j = use_genotype.size() > 1 ? use_genotype[1] : called_i;
                const vector<int>& panel = pr.panel_alleles(*lookup);
                // Retract the old entry and record the site again, so that there is one way a site
                // enters the linkage model. Retract first: `live_index` returns the first live entry
                // for a key, so the new entry is the live one.
                if (had_entry) {
                    model->retract(pr.record_key);
                }
                model->record(
                    locus.contig, locus.position, used->genotype_lls, panel,
                    called_i, called_j, trav_to_allele_vec,
                    // The explained share the old entry carried, or 1.0 for a chain that had none.
                    // The quality inputs of the direct call recorded here, which the record's
                    // GQI and GL also come from.
                    pr.record_key, direct_quality_of(reader.genotyper, *used),
                    (size_t)copies, reader.graph->get_id(pr.bounds.start),
                    reader.graph->get_id(pr.bounds.end),
                    LinkageCollector::SiteContext{
                        .nested = copies == 1,
                        .parent_record_key = pr.parent_record_key,
                        .parent_crossing = pr.parent_crossing,
                        .level = pr.level,
                        .emitted = false,
                        .unpositioned = pr.no_reference,
                        .chain_key = pr.chain_key,
                        .freq_prior = freq_prior(pr.travs, pr.ref_trav_idx),
                    });
                if (!model->has_entry(pr.record_key)) {
                    // `record` adds nothing for a site whose compact space it cannot describe: no called
                    // traversal, no likelihoods, or more than 127 alleles. The old entry stays
                    // retracted, and the record keeps its per-site call.
                    if (had_entry) {
                        ++unspecifiable;
                    }
                }
            }
            if (had_entry) {
                ++revised;
            } else {
                ++gained;
            }

            // This chain's chosen pair has changed, or the chain is new, so its children's crossing
            // masks are computed again.
            if (const vector<size_t>* kids = sites.children_of(pr.record_key)) {
                // Once for this parent: see TraversalNodeIndex.
                vector<ChildPlacer::TraversalNodeIndex> pr_visits;
                pr_visits.reserve(pr.travs.size());
                for (const Traversal& t : pr.travs) {
                    pr_visits.push_back(ChildPlacer::index_traversal_nodes(*reader.graph, t));
                }
                for (size_t ci : *kids) {
                    StagedSite& child = pending[ci];
                    bool known = true;
                    child.parent_crossing =
                        ChildPlacer::child_crossing_mask(*reader.graph, pr_visits, child.bounds,
                                                         &known);
                    child.crossing_known = known;
                }
            }
        }
        if (model != nullptr) {
            levels = max(levels, model->max_level());
        }
    }

    return counts;
}

// Iterative, over an explicit stack, since the depth depends on the data.
size_t GenotypeLinker::drop_subtree(StagedSiteTable& sites, size_t root) {
    vector<StagedSite>& pending = sites.nested();
    size_t dropped_here = 0;
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
    // Time the pass and report the collector's size.
    auto start = std::chrono::steady_clock::now();
    // The phase calls accumulate across levels, since the model needs the earlier ones: a
    // nested site's strand is read from its parent's PhaseCall, and a clamped site's phase is
    // pinned to its chosen pair.
    const size_t moved = model->resolve_level(level, last, calls);
    double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    total_seconds += seconds;
    // How many sites the model moved off the genotype the reads alone chose.
    total_moved += moved;
    if (!last) {
        // One line per level except the last: its site count, how many of its genotypes the
        // linkage model moved, and the seconds it took.
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
        // Duplicate keys need not change the output, but `retract` cannot handle those sites, since
        // it retracts only the first live entry.
        cerr << "[vg call] linkage: " << model->num_duplicate_live_keys()
             << " sites recorded onto a key that already had a live entry; the retract path cannot"
             << " address these" << endl;
    }
    if (model->model_params().hp_prior > 0.0) {
        cerr << "[vg call] linkage: " << model->num_site_prior_entries()
             << " live entries decoded at a run-length site's own frequency exponent (--hp-prior)"
             << endl;
    }
}

}
