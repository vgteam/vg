#include <algorithm>
#include <limits>

#include "tree_genotyper.hpp"
#include "graph_caller.hpp"
#include "site_values.hpp"
#include "vcf_record.hpp"

//#define debug

namespace vg {

void TreeGenotyper::configure(const Parts& parts, const Options& options) {
    this->parts = parts;
    this->options = options;
}

bool TreeGenotyper::genotype(const Snarl& managed_snarl) {
    // A top-level site has no parent context.
    return genotype_tree(managed_snarl, "", make_pair(0, 0), nullptr, -1, NestingPlacement());
}

bool TreeGenotyper::genotype_tree(const Snarl& managed_snarl, const string& parent_ref_path_name,
                                  pair<size_t, size_t> parent_ref_interval,
                                  const ChildTraversalSets* parent_child_trav_sets,
                                  int ploidy_override, const NestingPlacement& placement) {
    // Staged in the nested branch below and completed after descent, which reads `travs`, since
    // this record then takes ownership of them.
    unique_ptr<StagedSite> pending_this;
    // The panel alleles `GenotypeLinker::add` looked up for this snarl, if it filed the site. The staged
    // record keeps them, so that re-genotyping does not look them up again; the record's
    // traversals are this snarl's `travs`, which do not change after the site is recorded.
    vector<int> site_panel;
    bool site_panel_set = false;
    // The same, for a snarl the linkage pass will not revise.
    unique_ptr<StagedSite> render_this;

    // The reference path, the site oriented along it, and its candidate traversals.
    CandidateFinder::Site site;
    if (!parts.candidates->find(managed_snarl, parent_ref_path_name, parent_ref_interval,
                                parent_child_trav_sets, placement.no_reference, site)) {
        return false;
    }
    Snarl& snarl = site.snarl;
    const string& ref_path_name = site.ref_path_name;
    const tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t>& ref_interval =
        site.ref_interval;
    const bool use_parent_interval = site.use_parent_interval;
    vector<SnarlTraversal>& travs = site.travs;
    int ref_trav_idx = site.ref_trav_idx;
    // The site as plain values: its bounds, what it holds, and the sites around it. `walks` are
    // `travs` as walks, filled once `travs` is complete.
    const PathPositionHandleGraph& graph = *parts.graph;
    const SiteBounds bounds = bounds_of(graph, snarl);
    const SiteChildren children = site_children(*parts.snarl_manager, graph, snarl);
    const vector<SiteBounds> enclosing = enclosing_sites(*parts.snarl_manager, graph, snarl);
    vector<Traversal> walks;

    bool ret_val = true;
    vector<int> trav_genotype;  // Declared outside block so we can pass to children

    // A ploidy from the parent overrides the contig's or the region BED's: it is the number of the
    // parent's called alleles that reach this child. The region's is still the number of the
    // sample's haplotypes here, which the depth term needs.
    const int region_ploidy =
        parts.ploidy_regions->ploidy_at(ref_path_name, get<0>(ref_interval),
                                        ref_offset_of(*parts.ref_offsets, ref_path_name),
                                        ref_ploidy_of(*parts.ref_ploidies, ref_path_name));
    int ploidy = ploidy_override >= 0 ? ploidy_override : region_ploidy;

    // What both the parent-traversal-set branch and the top-level branch do with their genotype:
    // give it to the linkage model and stage it, so that `render_retained_records` writes its
    // record after the linkage pass. `trav_call_info` differs between the branches, so it is a
    // parameter. `snarl` is captured by reference; the candidate finder may have flipped it.
    auto stage = [&](unique_ptr<SnarlCaller::CallInfo>& trav_call_info, SiteScore* trav_score) {
        site_panel_set = parts.linker->add(bounds, walks, trav_genotype, trav_score, ref_trav_idx,
                                           ref_path_name,
                                           ref_offset_of(*parts.ref_offsets, ref_path_name),
                                           parts.record_key_of(snarl), placement, false, 0,
                                           &site_panel);
        render_this = stage_render_record(snarl, trav_genotype, ref_trav_idx, trav_call_info,
                                          trav_score, ref_path_name,
                                          ref_offset_of(*parts.ref_offsets, ref_path_name), ploidy);
    };

    if (parent_child_trav_sets != nullptr && !parent_child_trav_sets->empty()) {
        // Genotype using bounded search over traversal sets from parent
        // Each set contains traversals consistent with one parent allele
        ploidy = parent_child_trav_sets->size();

        // Track which set each traversal index belongs to (for phase consistency)
        // set_membership[i] = which parent allele set traversal i came from, or -1 if from finder
        vector<int> set_membership(travs.size(), -1);

        // Merge traversals from sets into travs, tracking membership

        for (int set_idx = 0; set_idx < ploidy; ++set_idx) {
            const TraversalSet& tset = (*parent_child_trav_sets)[set_idx];

            if (tset.empty()) {
                // Empty set means parent allele doesn't traverse this child (star allele)
                continue;
            }

            // Add traversals from this set to travs (avoiding duplicates)
            // Keep track of indices for this set
            for (const SnarlTraversal& trav : tset) {
                // Check if this traversal already exists in travs
                int match_idx = -1;
                for (int i = 0; i < travs.size() && match_idx < 0; ++i) {
                    if (travs[i] == trav) {
                        match_idx = i;
                    }
                }

                if (match_idx < 0) {
                    // New traversal - add it
                    match_idx = travs.size();
                    travs.push_back(trav);
                    set_membership.push_back(set_idx);
                } else if (set_membership[match_idx] < 0) {
                    // Traversal was from finder, now claim it for this set
                    set_membership[match_idx] = set_idx;
                }
                // Note: if already claimed by another set, that's fine (shared region)
            }
        }
        walks = walks_of(graph, travs);

        // Which parent haplotypes actually traverse this child? A parent allele with
        // an empty traversal set skips the child entirely, and gets a star or
        // missing allele rather than a genotype.
        vector<int> traversing_sets;
        for (int set_idx = 0; set_idx < ploidy; ++set_idx) {
            if (!(*parent_child_trav_sets)[set_idx].empty()) {
                traversing_sets.push_back(set_idx);
            }
        }

        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        SiteScore* trav_score = nullptr;
        int marker = options.star_allele ? STAR_ALLELE_MARKER : MISSING_ALLELE_MARKER;

        if (traversing_sets.empty()) {
            // No parent allele traverses this child at all.
            trav_genotype.assign(ploidy, marker);
        } else {
            // Genotype at the ploidy that passes through the site, not at the parent's ploidy, since
            // a site only one strand reaches is not diploid. genotype() returns a sorted multiset, so
            // the alleles are then placed on the strands that pass through.
            int effective_ploidy = (int)traversing_sets.size();
            vector<int> called_alleles;
            std::tie(called_alleles, trav_call_info) = genotype_site(
                bounds, walks, ref_trav_idx,
                Ploidies{.ploidy = effective_ploidy, .region_ploidy = region_ploidy}, enclosing,
                ref_path_name, make_pair(get<0>(ref_interval), get<1>(ref_interval)), trav_score);

            // Scatter the called alleles back onto the traversing haplotypes,
            // leaving the others as star/missing.
            trav_genotype.assign(ploidy, marker);
            for (size_t j = 0; j < traversing_sets.size() && j < called_alleles.size(); ++j) {
                trav_genotype[traversing_sets[j]] = called_alleles[j];
            }
        }

        // Staged only if the snarl is on a reference path.
        if (!use_parent_interval) {
            stage(trav_call_info, trav_score);
        }

        ret_val = trav_genotype.size() == ploidy;
    } else if (ploidy_override >= 0) {
        // A nested chain, reached by descent, at the ploidy its parent implied. Only a nested chain
        // can have its ploidy revised at the linkage pass, so only it needs the other ploidy's answer.
        walks = walks_of(graph, travs);
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        SiteScore* trav_score = nullptr;
        std::tie(trav_genotype, trav_call_info) = genotype_site(
            bounds, walks, ref_trav_idx,
            Ploidies{.ploidy = ploidy, .region_ploidy = region_ploidy, .also_score_other = true},
            enclosing, ref_path_name, make_pair(get<0>(ref_interval), get<1>(ref_interval)),
            trav_score);

        const bool retain_only = placement.retain_only;
        // Whether this snarl's own boundaries are on no reference path, checked from the graph for
        // each snarl.
        const bool no_ref_position = use_parent_interval;

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);
        if (no_ref_position) {
            // Genotyped and recorded, never written. Checked before retain_only, which does not
            // record.
            site_panel_set = parts.linker->add(
                bounds, walks, trav_genotype, trav_score, ref_trav_idx, ref_path_name,
                ref_offset_of(*parts.ref_offsets, ref_path_name), parts.record_key_of(snarl),
                placement, /*no_reference*/ true,
                // The parent's position, as `get_ref_position` gives it from the interval
                // `use_parent_interval` set, plus the chain's offset along its parent, as
                // `StagedSite::position_from_parent` has it.
                base_path_position(ref_path_name,
                                   get<0>(ref_interval)
                                       + ref_offset_of(*parts.ref_offsets, ref_path_name))
                    + (int64_t)placement.parent_offset,
                &site_panel);
            ++parts.descent_counters->no_ref_recorded;
            {
                int copies = 0;
                for (int a : trav_genotype) {
                    copies += (a >= 0);
                }
                parts.descent_counters->no_ref_copies[copies < 3 ? copies : 2].fetch_add(1);
            }
        } else if (retain_only) {
            // No called parent allele reaches this chain, so nothing about it is written yet. It is
            // genotyped and kept, since the linkage model may still move the parent onto an allele
            // that reaches it.
        } else if (placement.reported_inline) {
            // An enclosing block's ALT already spells this chain, so it gets no line, but it is
            // genotyped and recorded, since its allele pair phases everything inside it. Checked
            // after retain_only, which does not record.
            site_panel_set = parts.linker->add(bounds, walks, trav_genotype, trav_score,
                                               ref_trav_idx, ref_path_name,
                                               ref_offset_of(*parts.ref_offsets, ref_path_name),
                                               parts.record_key_of(snarl), placement, false, 0,
                                               &site_panel);
        } else {
            // Recorded here, and staged below, as at top level: the line is written after the
            // linkage pass, from the chosen genotype. A retained chain, on the path above, is
            // recorded only if the linkage pass later finds that the sample carries it.
            site_panel_set = parts.linker->add(bounds, walks, trav_genotype, trav_score,
                                               ref_trav_idx, ref_path_name,
                                               ref_offset_of(*parts.ref_offsets, ref_path_name),
                                               parts.record_key_of(snarl), placement, false, 0,
                                               &site_panel);
        }

        // Stage the nested site without its traversals: descent below still reads `travs` to find
        // which children the called alleles reach, and they are moved in once descent is done.
        pending_this.reset(new StagedSite());
        pending_this->bounds = bounds;
        fill_tree_fields(snarl, children, *pending_this);
        pending_this->ref_path_name = ref_path_name;
        pending_this->ref_offset = ref_offset_of(*parts.ref_offsets, ref_path_name);
        pending_this->ref_trav_idx = ref_trav_idx;
        pending_this->genotype = trav_genotype;
        pending_this->ploidy = ploidy;
        pending_this->record_key = parts.record_key_of(snarl);
        pending_this->parent_record_key = placement.parent_record_key;
        pending_this->parent_crossing = placement.parent_crossing;
        pending_this->chain_key = placement.chain_key;
        pending_this->no_reference = no_ref_position;
        pending_this->reported_inline = placement.reported_inline;
        pending_this->position_from_parent =
            no_ref_position
                ? base_path_position(ref_path_name,
                                     get<0>(ref_interval)
                                         + ref_offset_of(*parts.ref_offsets, ref_path_name))
                      + (int64_t)placement.parent_offset
                : 0;
        pending_this->chain_offset = placement.parent_offset;
        pending_this->crossing_known = placement.crossing_known;
        pending_this->level = (uint8_t)min(placement.level, (size_t)255);
        pending_this->set_call(std::move(trav_call_info), trav_score);
        ret_val = trav_genotype.size() == ploidy;
    } else {
        // Top-level snarl or no parent context - genotype from scratch using support
        walks = walks_of(graph, travs);
        unique_ptr<SnarlCaller::CallInfo> trav_call_info;
        SiteScore* trav_score = nullptr;
        std::tie(trav_genotype, trav_call_info) = genotype_site(
            bounds, walks, ref_trav_idx, Ploidies{.ploidy = ploidy}, enclosing, ref_path_name,
            make_pair(get<0>(ref_interval), get<1>(ref_interval)), trav_score);

        assert(trav_genotype.empty() || trav_genotype.size() == ploidy);

        stage(trav_call_info, trav_score);

        ret_val = trav_genotype.size() == ploidy;
    }

    // Nested calling: descend into each child the called alleles reach, at the ploidy they reach
    // it with.
    //
    // Descent does not depend on whether a line was written: a parent written as the reference
    // still has children to call. Children are genotyped independently, with no parent traversal
    // sets. Only a successful call descends, since a failed snarl has no genotype to take a child's
    // ploidy from. RecurseOnFail calls the children of a failed top-level snarl as top-level
    // snarls, but nothing does so for a failed nested snarl: its children are not called.
    if (ret_val && options.nested_calling && !trav_genotype.empty() &&
        parent_child_trav_sets == nullptr) {
        // A child no called allele reaches is kept only where the linkage pass can come back to
        // it: with the linkage model.
        parts.child_placer->place(
            snarl, children, parts.record_key_of(snarl), walks, trav_genotype, ref_trav_idx, ploidy,
            placement, options.off_reference, parts.linker->enabled(),
            [&](const ChildPlacer::Placed& child) {
                if (child.placement.level < 16) {
                    ++parts.descent_counters->depth_hist[child.placement.level];
                }
                // The other ploidy's answer is computed as well, so the linkage pass can change it
                // later.
                genotype_tree(*child.snarl, ref_path_name,
                              make_pair(get<0>(ref_interval), get<1>(ref_interval)), nullptr,
                              child.ploidy, child.placement);
            });
    }


    // --top-down: genotype each child against the traversals its parent's called alleles allow.
    if (options.top_down && !trav_genotype.empty()) {
        // Find the managed snarl pointer so we can get its children
        const Snarl* managed_ptr = parts.snarl_manager->into_which_snarl(snarl.start().node_id(), snarl.start().backward());
        if (managed_ptr) {
            for (const Snarl* child : parts.snarl_manager->children_of(managed_ptr)) {
                if (child && !parts.snarl_manager->is_trivial(child, *parts.graph)) {
                    // Build ChildTraversalSets: one set per parent allele
                    // Each set contains all traversals through child consistent with that parent allele
                    ChildTraversalSets child_trav_sets;
                    bool any_real_traversals = false;

                    for (int allele_idx : trav_genotype) {
                        if (allele_idx >= 0 && allele_idx < travs.size()) {
                            // Find all traversals through child consistent with this parent traversal
                            TraversalSet tset = parts.candidates->find_child_traversal_set(travs[allele_idx], *child);
                            if (!tset.empty()) {
                                any_real_traversals = true;
                            }
                            child_trav_sets.push_back(std::move(tset));
                        } else {
                            // Star/missing allele - pass empty set
                            child_trav_sets.push_back(TraversalSet());
                        }
                    }

                    // If no genotyped alleles traverse the child, skip it
                    if (!any_real_traversals) {
                        continue;
                    }

                    // Recursively call child with traversal sets
                    genotype_tree(*child, ref_path_name,
                                  make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                                  &child_trav_sets, -1, placement);
                }
            }
        }
    }

    // Descent above and the --top-down recursion, which builds each child's ChildTraversalSets
    // from `travs`, are done, so the staged site can take the traversals. At most one of these
    // is set.
    if (pending_this != nullptr) {
        pending_this->travs = std::move(walks);
        if (site_panel_set) {
            pending_this->panel_cache = std::move(site_panel);
            pending_this->panel_cached = true;
        }
        parts.staged_sites->add_nested(std::move(*pending_this));
        pending_this.reset();
    } else if (render_this != nullptr) {
        render_this->travs = std::move(walks);
        fill_tree_fields(snarl, children, *render_this);
        if (site_panel_set) {
            render_this->panel_cache = std::move(site_panel);
            render_this->panel_cached = true;
        }
        parts.staged_sites->add_top_level(std::move(*render_this));
        render_this.reset();
    }

    return ret_val;
}

pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> TreeGenotyper::genotype_site(
    const SiteBounds& site, const vector<Traversal>& travs, int ref_trav_idx,
    const Ploidies& ploidies, const vector<SiteBounds>& enclosing, const string& ref_path_name,
    pair<size_t, size_t> ref_range, SiteScore*& score) const {
    auto called = parts.genotyper->genotype(site, travs, ref_trav_idx, ploidies, enclosing,
                                            ref_path_name, ref_range);
    score = called.second.get();
    return make_pair(std::move(called.first),
                     unique_ptr<SnarlCaller::CallInfo>(std::move(called.second)));
}

// The CallInfo is kept because update_vcf_info reads it when the record is rendered, to map the
// written alleles back to matrix columns, index GL and compute QUAL.
unique_ptr<StagedSite> TreeGenotyper::stage_render_record(
        const Snarl& snarl, const vector<int>& trav_genotype, int ref_trav_idx,
        unique_ptr<SnarlCaller::CallInfo>& call_info, SiteScore* score,
        const string& ref_path_name, int ref_offset, int ploidy) const {
    unique_ptr<StagedSite> rec(new StagedSite());
    rec->bounds = bounds_of(*parts.graph, snarl);
    rec->ref_path_name = ref_path_name;
    rec->ref_offset = ref_offset;
    rec->ref_trav_idx = ref_trav_idx;
    rec->genotype = trav_genotype;
    rec->ploidy = ploidy;
    rec->record_key = parts.record_key_of(snarl);
    rec->level = 0;
    rec->set_call(std::move(call_info), score);
    // `travs` is not moved here: descent runs after the emit and reads `travs` to find which
    // children the called alleles reach. The caller completes the record after descent.
    return rec;
}

void TreeGenotyper::fill_tree_fields(const Snarl& snarl, const SiteChildren& children,
                                     StagedSite& site) const {
    site.children = children;
    // The snarl a walk enters by the site's start, as the record steps and the linkage pass
    // looked it up from the site itself: it decides whether the site is a leaf, and its chain.
    const Snarl* entered = parts.snarl_manager->into_which_snarl(snarl.start().node_id(),
                                                                 snarl.start().backward());
    site.leaf = entered == nullptr || parts.snarl_manager->children_of(entered).empty();
    site.in_chain = entered != nullptr;
    if (entered != nullptr) {
        site.chain = chain_of_site(*parts.snarl_manager, *parts.graph, entered);
    }
}

}
