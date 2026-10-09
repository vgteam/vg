#include <chrono>
#include <functional>

#include <omp.h>

#include "flow_caller.hpp"
#include "read_likelihood_caller.hpp"

namespace vg {

/// The genotype the linkage model chose, or the direct pass's own where it chose none. Used by both
/// anchor-collection paths, the render and `hand_off_deferred_records`, so that records with no
/// VCF line (`reported_inline` and `no_reference`) also get anchors for their chosen genotype.
vector<int> FlowCaller::chosen_genotype_for(const StagedSite& rec) const {
    vector<int> genotype = rec.genotype;
    int chosen_a = -1, chosen_b = -1;
    size_t chosen_ploidy = 0;
    if (linkage_collector != nullptr
        && linkage_collector->chosen_traversals(rec.record_key, &chosen_a, &chosen_b,
                                                 &chosen_ploidy)
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
static FlowCaller::SiteLocus locus_of(pair<string, int64_t> pos_info) {
    const string locus = PathMetadata::parse_locus_name(pos_info.first);
    if (locus != PathMetadata::NO_LOCUS_NAME) {
        pos_info.first = locus;
    }
    return FlowCaller::SiteLocus{pos_info.first, (size_t)max((int64_t)0, pos_info.second)};
}

FlowCaller::SiteLocus FlowCaller::site_locus(const Snarl& snarl, const string& ref_path_name,
                                             int ref_offset) const {
    // The position before the record's alleles are trimmed, which can move POS.
    // `get_ref_position` names the base path, as in "CHM13#0#chr20".
    return locus_of(get_ref_position(graph, snarl, ref_path_name, ref_offset));
}

FlowCaller::SiteLocus FlowCaller::off_reference_site_locus(const string& ref_path_name,
                                                           int64_t stand_in_position) const {
    // Not `get_ref_position`: `get_ref_interval` asserts on a snarl the path does not pass
    // through.
    return locus_of(make_pair(ref_path_name, stand_in_position));
}

/// The quality inputs of the direct call `info`, which the linkage collector keeps for rewriting
/// the record if the model moves it. `caller` decides the GQ factor, since --no-share-quality and
/// --depth-quality are its settings.
static LinkageCollector::DirectQuality direct_quality_of(
    const SnarlCaller& caller, const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo& info) {
    const auto* rl_caller = dynamic_cast<const ReadLikelihoodSnarlCaller*>(&caller);
    return LinkageCollector::DirectQuality{
        .explained_share = info.explained_share,
        .gq_factor = rl_caller != nullptr ? rl_caller->gq_factor(info) : info.explained_share,
        .achievable_gap = info.achievable_gap,
    };
}

bool FlowCaller::record_site(const Snarl& snarl, const vector<SnarlTraversal>& travs,
                            const vector<int>& trav_genotype,
                            const unique_ptr<SnarlCaller::CallInfo>& call_info, int ref_trav_idx,
                            const string& ref_path_name, int ref_offset,
                            bool no_reference, int64_t position_from_parent,
                            vector<int>* panel_out) {
    if (linkage_collector == nullptr) {
        return false;
    }
    const auto* rl_info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info.get());
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
                                : site_locus(snarl, ref_path_name, ref_offset);
    const int called_i = trav_genotype[0];
    const int called_j = site_ploidy > 1 ? trav_genotype[1] : called_i;
    // No allele map yet: the written alleles are chosen when the record is built, and
    // `set_allele_map` supplies the map then.
    static const vector<int> no_allele_map;
    vector<int> panel = panel_lookup.alleles(travs);
    linkage_collector->record(
        locus.contig, locus.position,
        rl_info->genotype_lls,
        panel,
        called_i, called_j, no_allele_map,
        record_key_of(snarl),
        direct_quality_of(snarl_caller, *rl_info), site_ploidy,
        (int64_t)snarl.start().node_id(), (int64_t)snarl.end().node_id(),
        // `nested` only when one copy of the chain is present, as for any other chain; a chain with
        // two copies joins its parent's diploid group.
        LinkageCollector::SiteContext{
            .nested = nested_context.one_copy,
            .parent_record_key = nested_context.parent_record_key,
            .parent_crossing = nested_context.parent_crossing,
            .level = current_level,
            .emitted = false,
            .unpositioned = no_reference,
            .chain_key = nested_context.chain_key,
            .freq_prior = site_freq_prior(travs, ref_trav_idx),
        });
    if (panel_out != nullptr) {
        *panel_out = std::move(panel);
    }
    return true;
}

double FlowCaller::site_freq_prior(const vector<SnarlTraversal>& travs, int ref_trav_idx) const {
    const LinkageModel::Params& params = linkage_collector->model_params();
    if (params.hp_prior <= 0.0) {
        return -1.0;
    }
    vector<string> alleles;
    alleles.reserve(travs.size());
    for (const SnarlTraversal& trav : travs) {
        alleles.push_back(trav_string(graph, trav));
    }
    const size_t ref = ref_trav_idx >= 0 ? (size_t)ref_trav_idx : (size_t)-1;
    return LinkageModel::run_length_site(alleles, params.hp_prior_run, ref) ? params.hp_prior : -1.0;
}

void FlowCaller::rerun_linkage_pass() {
    if (linkage_collector == nullptr) {
        return;
    }
    // Give the linkage model the corrected likelihoods, then run the linkage pass again in full, so that
    // every child is reassessed against its parent's new chosen pair, as on the first pass.
    const vector<StagedSite*> records = staged_sites.in_order();
    // The loop below is serial, and most of its time would go to each record's first
    // `panel_lookup.alleles`, a GBWT lookup per allele. Each record's lookup is independent of the
    // others', so fill the caches in parallel first.
#pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < records.size(); ++i) {
        records[i]->panel_alleles(panel_lookup);
    }
    size_t rescored = 0, refused = 0;
    for (StagedSite* recp : records) {
        StagedSite& rec = *recp;
        const auto* info = dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
            rec.call_info.get());
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
        if (linkage_collector->rescore(rec.record_key, info->genotype_lls,
                                       rec.panel_alleles(panel_lookup), called_i, called_j)) {
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
    // The linkage pass again, in full.
    run_linkage_pass();
}

void FlowCaller::run_linkage_pass() {
    if (!staged_sites.active()) {
        return;
    }
    // Descent already happened during the direct pass. What is left is to choose the chains' genotypes
    // in the order their ploidies depend on: a level's parents before its children. The reads are not
    // used.
    size_t levels = 0;
    if (linkage_collector != nullptr) {
        levels = linkage_collector->max_level();
    }
    // Gathered from the per-thread lists on the first pass only, and kept, since re-genotyping
    // runs the linkage pass again.
    vector<StagedSite>& pending = staged_sites.gather_nested();

    // On a later pass, everything a pass concludes is derived again from the chosen genotypes.
    // `dropped` is cleared, since a correction can move a parent onto an allele that crosses a
    // dropped chain; the level loop then drops again the chains still not crossed, and records
    // the others afresh.
    if (linkage_passes_run > 0) {
        for (StagedSite& pr : pending) {
            pr.dropped = false;
        }
        // Appended to by every resolve, so cleared here.
        phase_table.calls().clear();
        // Accumulated by every resolve too.
        linkage_changed = 0;
    }
    ++linkage_passes_run;

    // Counters for the report.
    size_t revise_unrenderable = 0, pass_no_crossing = 0, pass_no_chosen = 0, pass_ploidy_unscored = 0;
    size_t pass_inline_rederived = 0;
    unordered_map<size_t, StagedSite*> record_by_key = staged_sites.by_key();
    // Each parent traversal's child offsets, for placing its chains. Keyed by address, which stays
    // valid for the pass as record_by_key's do, and the traversals do not change after the direct
    // pass.
    unordered_map<const SnarlTraversal*, ChildOffsets> child_offsets;
    // Drop a chain and its whole subtree: the chosen parent does not carry the chain, so the
    // sample has no copy of it or of anything inside it. Returns how many entries were retracted.
    // Iterative, over an explicit stack, since the depth depends on the data.
    std::function<size_t(size_t)> drop_subtree = [&](size_t root) -> size_t {
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
            if (linkage_collector != nullptr && linkage_collector->retract(victim.record_key)) {
                ++dropped_here;
            }
            if (const vector<size_t>* kids = staged_sites.children_of(victim.record_key)) {
                for (size_t k : *kids) {
                    if (k != idx) {
                        stack.push_back(k);
                    }
                }
            }
        }
        return dropped_here;
    };

    size_t revised = 0, retracted = 0, gained = 0, crossing_unknown = 0, unspecifiable = 0;
    // One linkage pass per level, in order. `levels` is read again after each pass, since
    // a pass can add a chain at a deeper level, which must still be chosen.
    for (size_t gen = 0; gen <= levels; ++gen) {
        // The final pass has last=true and builds the phasing map and the mosaic from everything
        // accumulated. If it adds a deeper chain, the bound grows and a later pass rebuilds them.
        resolve_linkage_level(gen, gen == levels);

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
        for (const LinkageCollector::PhaseCall& pc : phase_table.calls()) {
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
            if (linkage_collector != nullptr && parent_record != record_by_key.end()) {
                // Place the chain along the allele chosen for its parent, before this level's
                // linkage pass orders and spaces its sites by position. The parent's own offset
                // was placed in the previous level's iteration.
                const StagedSite& par = *parent_record->second;
                const size_t offset =
                    par.chain_offset
                    + offset_along_genotype(par.travs, chosen_genotype_for(par), pr.snarl,
                                            child_offsets);
                if (offset != pr.chain_offset) {
                    if (pr.no_reference) {
                        pr.position_from_parent += (int64_t)offset - (int64_t)pr.chain_offset;
                        linkage_collector->set_position(
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
                retracted += drop_subtree(i);
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
            } else if (linkage_collector != nullptr && parent_record != record_by_key.end()
                       && !parent_record->second->genotype.empty()) {
                // The linkage model gave the parent no PhaseCall, so the parent is rendered at its
                // own chosen genotype, which `chosen_genotype_for` reads.
                const vector<int> parent_genotype = chosen_genotype_for(*parent_record->second);
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
                retracted += drop_subtree(i);
                continue;
            }
            if (copies == pr.ploidy && linkage_collector != nullptr
                && linkage_collector->has_entry(pr.record_key)) {
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
            ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo* rl =
                dynamic_cast<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(pr.call_info.get());
            unique_ptr<SnarlCaller::CallInfo> use_info;
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
                unique_ptr<SnarlCaller::CallInfo> demoted = std::move(pr.call_info);
                // `rl` still points at it -- `demoted` owns what `pr.call_info` did.
                unique_ptr<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo> promoted(
                    rl->alt_ploidy_info.release());
                // The fields that do not depend on ploidy go with whichever answer is in front, since
                // the alternate does not copy them.
                promoted->anchor_evidence = std::move(rl->anchor_evidence);
                promoted->phase_evidence = std::move(rl->phase_evidence);
                // The replaced answer becomes the new alternate, with the genotype it was called at.
                promoted->alt_ploidy_best = demoted_genotype;
                promoted->alt_ploidy_info.reset(
                    static_cast<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                        demoted.release()));
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
            const bool had_entry = linkage_collector != nullptr
                                   && linkage_collector->has_entry(pr.record_key);
            // Revise the staged site; the render builds its line once, at the end, from the chosen
            // genotype.
            pr.genotype = use_genotype;
            pr.ploidy = copies;
            if (use_info != nullptr) {
                pr.call_info = std::move(use_info);
            }
            const unique_ptr<SnarlCaller::CallInfo>& info = pr.call_info;

            const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo* used =
                dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(info.get());
            if (used != nullptr) {
                // In traversal space, as `record_site` records it, so the linkage pass and the
                // direct pass describe a site the same way. No allele map yet, as in `record_site`.
                static const vector<int> no_allele_map;
                const vector<int>& trav_to_allele_vec = no_allele_map;
                const SiteLocus locus =
                    pr.no_reference
                        ? off_reference_site_locus(pr.ref_path_name, pr.position_from_parent)
                        : site_locus(pr.snarl, pr.ref_path_name, pr.ref_offset);
                int called_i = use_genotype.empty() ? -1 : use_genotype[0];
                int called_j = use_genotype.size() > 1 ? use_genotype[1] : called_i;
                const vector<int>& panel = pr.panel_alleles(panel_lookup);
                // Retract the old entry and record the site again, so that there is one way a site
                // enters the linkage model. Retract first: `live_index` returns the first live entry
                // for a key, so the new entry is the live one.
                if (had_entry) {
                    linkage_collector->retract(pr.record_key);
                }
                linkage_collector->record(
                    locus.contig, locus.position, used->genotype_lls, panel,
                    called_i, called_j, trav_to_allele_vec,
                    // The explained share the old entry carried, or 1.0 for a chain that had none.
                    // The quality inputs of the direct call recorded here, which the record's
                    // GQI and GL also come from.
                    pr.record_key, direct_quality_of(snarl_caller, *used),
                    (size_t)copies, pr.snarl.start().node_id(), pr.snarl.end().node_id(),
                    LinkageCollector::SiteContext{
                        .nested = copies == 1,
                        .parent_record_key = pr.parent_record_key,
                        .parent_crossing = pr.parent_crossing,
                        .level = pr.level,
                        .emitted = false,
                        .unpositioned = pr.no_reference,
                        .chain_key = pr.chain_key,
                        .freq_prior = site_freq_prior(pr.travs, pr.ref_trav_idx),
                    });
                if (!linkage_collector->has_entry(pr.record_key)) {
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
            if (const vector<size_t>* kids = staged_sites.children_of(pr.record_key)) {
                // Once for this parent: see TraversalNodeIndex.
                vector<TraversalNodeIndex> pr_visits;
                pr_visits.reserve(pr.travs.size());
                for (const SnarlTraversal& t : pr.travs) {
                    pr_visits.push_back(index_traversal_nodes(t));
                }
                for (size_t ci : *kids) {
                    StagedSite& child = pending[ci];
                    bool known = true;
                    child.parent_crossing = child_crossing_mask(pr_visits, child.snarl, &known);
                    child.crossing_known = known;
                }
            }
        }
        if (linkage_collector != nullptr) {
            levels = max(levels, linkage_collector->max_level());
        }
    }

    // The exactly-once test, from the chosen genotypes the render builds each parent's blocks
    // from: `reported_inline` holds back a chain's line where an enclosing block's ALT spells it.
    // Without the linkage model every chosen genotype is the direct call the direct pass tested, so
    // there is nothing to redo.
    if (linkage_collector != nullptr && staged_sites.has_children()) {
        // Counted again from here, so that the report gives the chains held back now.
        block_records.restart_inline_count();
        // Parents before their children, so that a chain inherits its parent's final flag.
        staged_sites.for_each_parent_top_down(record_by_key, [&](const StagedSite& parent,
                                                                 const vector<size_t>& children) {
            if (parent.dropped) {
                return;   // its children were dropped with it
            }
            // The parts of the test that do not depend on the child, built once for this parent;
            // see BlockRecordWriter::ChainInlineContext.
            const BlockRecordWriter::ChainInlineContext ctx = block_records.chain_inline_context(
                parent.snarl, parent.travs, chosen_genotype_for(parent), parent.ref_trav_idx);
            for (size_t ci : children) {
                StagedSite& child = pending[ci];
                if (child.dropped) {
                    continue;
                }
                const bool was = child.reported_inline;
                child.reported_inline = parent.reported_inline
                                        || block_records.chain_reported_inline(ctx, child.snarl);
                if (was != child.reported_inline) {
                    ++pass_inline_rederived;
                }
            }
        });
    }
    if (show_progress) {
        // The bytes kept for the staged sites, counted by walking the objects. They are walked in
        // parallel; the totals are sums, so they do not depend on how the walk is split.
        size_t retained_bytes = 0, retained_visits = 0, retained_gls = 0;
        auto measure = [](const StagedSite& rec, size_t& bytes, size_t& visits, size_t& gls) {
            bytes += sizeof(StagedSite) + rec.ref_path_name.capacity()
                     + rec.genotype.capacity() * sizeof(int)
                     + rec.panel_cache.capacity() * sizeof(int);
            bytes += rec.travs.capacity() * sizeof(SnarlTraversal);
            for (const SnarlTraversal& t : rec.travs) {
                visits += (size_t)t.visit_size();
                bytes += (size_t)t.visit_size() * sizeof(Visit);
            }
            const auto* rl = dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
                rec.call_info.get());
            if (rl != nullptr) {
                for (const auto& kv : rl->genotype_lls) {
                    ++gls;
                    bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                }
                if (rl->anchor_evidence != nullptr) {
                    bytes += rl->anchor_evidence->bytes();
                }
                if (rl->phase_evidence != nullptr) {
                    bytes += rl->phase_evidence->bytes();
                }
                // The parts re-genotyping adds.
                auto gl_bytes = [](const map<vector<int>, double>& gl) {
                    size_t n = 0;
                    for (const auto& kv : gl) {
                        n += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                    }
                    return n;
                };
                if (rl->uncorrected_lls != nullptr) {
                    bytes += gl_bytes(*rl->uncorrected_lls);
                }
                bytes += rl->scored_traversals.capacity() * sizeof(SnarlTraversal)
                         + rl->allele_support.capacity() * sizeof(double);
                if (rl->alt_ploidy_info != nullptr) {
                    // The alternate answer is kept too, with all its parts.
                    const auto& alt = *rl->alt_ploidy_info;
                    bytes += alt.scored_traversals.capacity() * sizeof(SnarlTraversal)
                             + alt.allele_support.capacity() * sizeof(double);
                    if (alt.uncorrected_lls != nullptr) {
                        bytes += gl_bytes(*alt.uncorrected_lls);
                    }
                }

                if (rl->alt_ploidy_info != nullptr) {
                    for (const auto& kv : rl->alt_ploidy_info->genotype_lls) {
                        ++gls;
                        bytes += 48 + kv.first.capacity() * sizeof(int) + sizeof(double);
                    }
                }
            }
        };
#pragma omp parallel reduction(+ : retained_bytes, retained_visits, retained_gls)
        {
            for (size_t q = 0; q < staged_sites.queue_count(); ++q) {
                const vector<StagedSite>& queue = staged_sites.queue(q);
#pragma omp for schedule(dynamic, 4096) nowait
                for (size_t r = 0; r < queue.size(); ++r) {
                    measure(queue[r], retained_bytes, retained_visits, retained_gls);
                }
            }
#pragma omp for schedule(dynamic, 4096) nowait
            for (size_t r = 0; r < pending.size(); ++r) {
                measure(pending[r], retained_bytes, retained_visits, retained_gls);
            }
        }
        // The read-phasing evidence. In the report below, the snarls the linkage pass will not revise
        // are the top-level ones and the children RecurseOnFail reaches without a ploidy override.
        for (const PhaseSite& ps : read_strands.sites()) {
            retained_bytes += sizeof(PhaseSite) + ps.read_key.capacity() * sizeof(uint64_t)
                              + ps.q0.capacity() * sizeof(float) + ps.c.capacity() * sizeof(float);
        }
        retained_bytes += read_strands.flips().size() * (sizeof(size_t) + 16);
        cerr << "[vg call] retained for rendering: " << staged_sites.queued_count()
             << " snarls the linkage pass will not revise, plus " << pending.size()
             << " nested chains; " << (retained_bytes / (1024.0 * 1024.0)) << " MB over "
             << retained_visits << " traversal visits and " << retained_gls
             << " genotype likelihoods" << endl;
        cerr << "[vg call] linkage pass exits: " << pass_no_crossing
             << " dropped because no parent candidate crosses them, " << pass_no_chosen
             << " whose parent's chosen pair could not be read, " << revise_unrenderable
             << " unrenderable so left unrevised, " << pass_ploidy_unscored
             << " stranded at a ploidy the direct pass never scored" << endl;
        if (pass_inline_rederived > 0) {
            cerr << "[vg call] linkage pass: " << pass_inline_rederived
                 << " children whose exactly-once suppression changed with their parent's"
                 << " chosen genotype" << endl;
        }
        cerr << "[vg call] single direct pass: " << pending.size() << " nested chains retained over "
             << (levels + 1) << " levels; " << revised << " revised, " << gained
             << " reachable only under the chosen parent, " << retracted << " retracted";
        if (crossing_unknown > 0) {
            cerr << ", " << crossing_unknown << " with a crossing mask the direct pass could not compute";
        }
        if (unspecifiable > 0) {
            cerr << ", " << unspecifiable << " dropped from the layer because the site's "
                 << "compact allele space could not be built";
        }
        cerr << endl;
    }

}

void VCFOutputCaller::resolve_linkage() {
    if (linkage_resolved) {
        return;
    }
    if (linkage_collector == nullptr) {
        resolve_linkage_level(0, true);
        return;
    }
    // Resolve every level, since chain construction skips entries of later levels
    // than the one being resolved. `max_level()` is read again on each pass, since a pass can
    // add a chain at a deeper level.
    for (size_t gen = 0;; ++gen) {
        const size_t deepest = linkage_collector->max_level();
        resolve_linkage_level(gen, gen >= deepest);
        if (gen >= deepest) {
            break;
        }
    }
}

void VCFOutputCaller::resolve_linkage_level(size_t level, bool last) {
    linkage_resolved = true;
    if (linkage_collector == nullptr) {
        return;
    }
    // Time the pass and report the collector's size.
    auto start = std::chrono::steady_clock::now();
    // The phase calls accumulate across levels, since the model needs the earlier ones: a
    // nested site's strand is read from its parent's PhaseCall, and a clamped site's phase is
    // pinned to its chosen pair.
    const size_t moved =
        linkage_collector->resolve_level(level, last,
                                              emit_phasing ? &phase_table.calls() : nullptr);
    double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    linkage_seconds += seconds;
    // How many sites the model moved off the genotype the reads alone chose.
    linkage_changed += moved;
    if (!last) {
        // One line per level except the last: its site count, how many of its genotypes the
        // linkage model moved, and the seconds it took.
        cerr << "[vg call] linkage level " << level << ": "
             << linkage_collector->num_sites_at(level) << " sites, "
             << moved << " genotypes moved by linkage, " << seconds << " s" << endl;
        return;
    }

}

}
