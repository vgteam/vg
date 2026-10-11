#include <algorithm>
#include <limits>

#include "tree_genotyper.hpp"
#include "graph_caller.hpp"
#include "site_values.hpp"
#include "vcf_record.hpp"

//#define debug

namespace vg {
namespace multipass {

/// The outcome of a candidate as it is placed, or why it is not, for the dump.
static EditOutcome unplaced_outcome(const EditCandidate& c, const string& why) {
    EditOutcome outcome;
    outcome.node = c.node;
    outcome.offset = c.offset;
    outcome.alt = c.alt;
    outcome.kind = c.kind;
    outcome.sequence = c.sequence;
    outcome.placed = why;
    return outcome;
}

void TreeGenotyper::configure(const Parts& parts, const Options& options) {
    this->parts = parts;
    this->options = options;
}

bool TreeGenotyper::genotype(const SiteView& site) {
    // A top-level site has no parent context.
    const bool genotyped = genotype_tree(site, "", make_pair(0, 0), nullptr, -1, NestingPlacement());
    if (parts.edit_source != nullptr && site.enclosing.empty()) {
        // The edit sites on the chain nodes this site lies between, which no site holds, while
        // their reads are in memory: its own read query fetched both bounds. A node between two
        // sites is taken by the first to reach it.
        const nid_t first = std::min(parts.graph->get_id(site.bounds.start),
                                     parts.graph->get_id(site.bounds.end));
        const nid_t last = std::max(parts.graph->get_id(site.bounds.start),
                                    parts.graph->get_id(site.bounds.end));
        vector<nid_t> bounds{first};
        if (last != first) {
            bounds.push_back(last);
        }
        for (nid_t bound : bounds) {
            vector<EditCandidate> on_bound;
            parts.edit_source->for_each_edit_candidate(
                {bound}, [&](const EditCandidate& c) { on_bound.push_back(c); });
            genotype_free_edits(std::move(on_bound));
        }
        // Then the nodes no site holds that run along the reference from either bound to the
        // next site, through the trivial sites the walk does not visit, as far as their reads come
        // with the lower bound's.
        vector<nid_t> near;
        for (nid_t bound : {first, last}) {
            free_nodes_along_reference(bound, near);
        }
        for (nid_t node : near) {
            if (node != first && node != last && parts.edit_source->fetched_together(first, node)) {
                vector<EditCandidate> on_node;
                parts.edit_source->for_each_edit_candidate(
                    {node}, [&](const EditCandidate& c) { on_node.push_back(c); });
                genotype_free_edits(std::move(on_node));
            }
        }
    }
    return genotyped;
}

bool TreeGenotyper::genotype_tree(const SiteView& view, const string& parent_ref_path_name,
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
    const PathPositionHandleGraph& graph = *parts.graph;
    if (!parts.candidates->find(snarl_of(graph, view.bounds), parent_ref_path_name,
                                parent_ref_interval, parent_child_trav_sets,
                                placement.no_reference, site)) {
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
    const SiteBounds bounds = bounds_of(graph, snarl);
    // The site's ID and record key, from its bounds as the candidate finder oriented them.
    const string site_id = parts.site_id(bounds);
    const size_t record_key = record_key_of(site_id);
    const ChildPlacer::Nested nested = parts.child_placer->nested_in(view.net, bounds);
    const SiteChildren& children = nested.children;
    const vector<SiteBounds>& enclosing = view.enclosing;
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
    // parameter.
    auto stage = [&](unique_ptr<SnarlCaller::CallInfo>& trav_call_info, SiteScore* trav_score) {
        site_panel_set = parts.linker->add(bounds, walks, trav_genotype, trav_score, ref_trav_idx,
                                           ref_path_name,
                                           ref_offset_of(*parts.ref_offsets, ref_path_name),
                                           record_key, placement, false, 0,
                                           &site_panel);
        render_this = stage_render_record(bounds, site_id, trav_genotype, ref_trav_idx,
                                          trav_call_info, trav_score, ref_path_name,
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
                ref_offset_of(*parts.ref_offsets, ref_path_name), record_key,
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
                                               record_key, placement, false, 0,
                                               &site_panel);
        } else {
            // Recorded here, and staged below, as at top level: the line is written after the
            // linkage pass, from the chosen genotype. A retained chain, on the path above, is
            // recorded only if the linkage pass later finds that the sample carries it.
            site_panel_set = parts.linker->add(bounds, walks, trav_genotype, trav_score,
                                               ref_trav_idx, ref_path_name,
                                               ref_offset_of(*parts.ref_offsets, ref_path_name),
                                               record_key, placement, false, 0,
                                               &site_panel);
        }

        // Stage the nested site without its traversals: descent below still reads `travs` to find
        // which children the called alleles reach, and they are moved in once descent is done.
        pending_this.reset(new StagedSite());
        pending_this->bounds = bounds;
        fill_tree_fields(view, children, *pending_this);
        pending_this->ref_path_name = ref_path_name;
        pending_this->ref_offset = ref_offset_of(*parts.ref_offsets, ref_path_name);
        pending_this->ref_trav_idx = ref_trav_idx;
        pending_this->genotype = trav_genotype;
        pending_this->ploidy = ploidy;
        pending_this->id = site_id;
        pending_this->record_key = record_key;
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
        // The off-panel candidates this site holds itself are children of it too, each in the
        // chain its node is in, placed after the decomposition's children.
        vector<HeldEdit> held;
        if (parts.edit_source != nullptr) {
            held = edits_held_by(view, bounds, walks, ref_trav_idx);
        }
        ChildPlacer::Nested with_edits;
        if (!held.empty()) {
            // Each before the first child that starts after its node, so that the children are
            // visited along the node IDs and each window's reads are still in memory when its
            // edit sites are genotyped. `held` is sorted by node.
            with_edits.children = nested.children;
            size_t e = 0;
            for (const ChildSite& child : nested.sites) {
                const nid_t child_first = std::min(graph.get_id(child.bounds.start),
                                                   graph.get_id(child.bounds.end));
                for (; e < held.size() && held[e].candidate.node < child_first; ++e) {
                    with_edits.sites.push_back(held[e].child);
                }
                with_edits.sites.push_back(child);
            }
            for (; e < held.size(); ++e) {
                with_edits.sites.push_back(held[e].child);
            }
        }
        vector<bool> edit_placed(held.size(), false);
        // A child no called allele reaches is kept only where the linkage pass can come back to
        // it: with the linkage model.
        parts.child_placer->place(
            view, held.empty() ? nested : with_edits, record_key, walks, trav_genotype,
            ref_trav_idx, ploidy, placement, options.off_reference, parts.linker->enabled(),
            [&](const ChildPlacer::Placed& child) {
                if (child.placement.level < 16) {
                    ++parts.descent_counters->depth_hist[child.placement.level];
                }
                if (child.site.bounds.inside_node) {
                    for (size_t i = 0; i < held.size(); ++i) {
                        if (!edit_placed[i] && held[i].child.bounds == child.site.bounds) {
                            edit_placed[i] = true;
                            genotype_edit_site(held[i], child.site.enclosing, ref_path_name,
                                               make_pair(get<0>(ref_interval),
                                                         get<1>(ref_interval)),
                                               child.ploidy, child.placement, false);
                            break;
                        }
                    }
                    return;
                }
                // The other ploidy's answer is computed as well, so the linkage pass can change it
                // later.
                genotype_tree(child.site, ref_path_name,
                              make_pair(get<0>(ref_interval), get<1>(ref_interval)), nullptr,
                              child.ploidy, child.placement);
            });
        for (size_t i = 0; i < held.size(); ++i) {
            if (!edit_placed[i]) {
                // Left out as any child would be: off the reference without off-reference
                // descent, or crossed by no called allele without the linkage model.
                parts.edit_outcomes->add(unplaced_outcome(held[i].candidate, "not_placed"));
            }
        }
    }


    // Top-down calling: genotype each child against the traversals its parent's called alleles
    // allow.
    if (options.top_down && !trav_genotype.empty() && children.known) {
        for (const ChildSite& child : nested.sites) {
            // Build ChildTraversalSets: one set per parent allele
            // Each set contains all traversals through child consistent with that parent allele
            ChildTraversalSets child_trav_sets;
            bool any_real_traversals = false;

            for (int allele_idx : trav_genotype) {
                if (allele_idx >= 0 && allele_idx < travs.size()) {
                    // Find all traversals through child consistent with this parent traversal
                    TraversalSet tset = parts.candidates->find_child_traversal_set(
                        travs[allele_idx], snarl_of(graph, child.bounds));
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
            genotype_tree(ChildPlacer::view_of(view, child), ref_path_name,
                          make_pair(get<0>(ref_interval), get<1>(ref_interval)),
                          &child_trav_sets, -1, placement);
        }
    }

    // Descent above and the top-down recursion, which builds each child's ChildTraversalSets
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
        fill_tree_fields(view, children, *render_this);
        if (site_panel_set) {
            render_this->panel_cache = std::move(site_panel);
            render_this->panel_cached = true;
        }
        parts.staged_sites->add_top_level(std::move(*render_this));
        render_this.reset();
    }

    return ret_val;
}

vector<EditCandidate> TreeGenotyper::admitted(vector<EditCandidate> candidates) const {
    // Sorted by node, offset, kind and ALT, so each base's substitutions are consecutive, as are
    // the indels anchored on it; one of each is kept.
    vector<EditCandidate> chosen;
    for (const EditCandidate& c : candidates) {
        if (!options.edit_bar.admits(c)) {
            continue;
        }
        if (!chosen.empty() && chosen.back().node == c.node && chosen.back().offset == c.offset
            && chosen.back().is_indel() == c.is_indel()) {
            if (c.alt_fragments > chosen.back().alt_fragments) {
                chosen.back() = c;
            }
            continue;
        }
        chosen.push_back(c);
    }
    return chosen;
}

bool TreeGenotyper::edit_site_of(const EditCandidate& c, handle_t read_as,
                                 const NodePlacement& where, EditSite& out) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    if (!c.is_indel()) {
        out = snv_edit_site(graph, read_as, c.offset, c.alt);
        return true;
    }
    if (graph.get_is_reverse(read_as)) {
        return false;
    }
    // Every node of its walk in the same site as its anchor.
    for (nid_t id : c.walk) {
        const NodePlacement other = parts.child_placer->placement_of(id);
        if (other.enclosing.size() != where.enclosing.size()
            || (!other.enclosing.empty() && other.enclosing.front() != where.enclosing.front())) {
            return false;
        }
    }
    // Its record is spelled along the walk and placed where a reference path reads the walk's
    // first node, so a reference path must read the whole walk there, forward and once, or not
    // read it at all.
    size_t on_reference = 0;
    int64_t expected = -1;
    for (nid_t id : c.walk) {
        const handle_t h = graph.get_handle(id);
        size_t steps = 0;
        int64_t position = -1;
        graph.for_each_step_on_handle(h, [&](const step_handle_t& step) {
            if (parts.ref_offsets->count(graph.get_path_name(graph.get_path_handle_of_step(step)))) {
                ++steps;
                position = graph.get_is_reverse(graph.get_handle_of_step(step))
                               ? -2
                               : (int64_t)graph.get_position_of_step(step);
            }
        });
        if (steps > 1 || position == -2 || (expected >= 0 && steps == 1 && position != expected)) {
            return false;
        }
        on_reference += steps;
        expected = steps == 1 ? position + (int64_t)graph.get_length(h) : -1;
    }
    if (on_reference != 0 && on_reference != c.walk.size()) {
        return false;
    }
    out = indel_edit_site(graph, c.walk, c.offset, c.kind == 'D', c.sequence);
    return true;
}

/// The bounds `bounds` read the other way.
static SiteBounds turned(const HandleGraph& graph, const SiteBounds& bounds) {
    SiteBounds out = bounds;
    out.start = graph.flip(bounds.end);
    out.end = graph.flip(bounds.start);
    return out;
}

vector<TreeGenotyper::HeldEdit> TreeGenotyper::edits_held_by(const SiteView& view,
                                                             const SiteBounds& bounds,
                                                             const vector<Traversal>& walks,
                                                             int ref_trav_idx) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    // The nodes the walks visit inside the site. Its own bounds are in its parent's chain.
    vector<nid_t> nodes;
    for (const Traversal& walk : walks) {
        for (const handle_t& h : walk) {
            nodes.push_back(graph.get_id(h));
        }
    }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    vector<EditCandidate> found;
    const nid_t start = graph.get_id(bounds.start);
    const nid_t end = graph.get_id(bounds.end);
    // The site's own read query counted them.
    parts.edit_source->for_each_edit_candidate(nodes, [&](const EditCandidate& c) {
        if (c.node != start && c.node != end) {
            found.push_back(c);
        }
    });
    vector<HeldEdit> held;
    if (found.empty()) {
        return held;
    }
    const SiteBounds view_turned = turned(graph, view.bounds);
    for (const EditCandidate& c : admitted(std::move(found))) {
        // Only the innermost site holding the node holds the edit.
        NodePlacement where = parts.child_placer->placement_of(c.node);
        if (where.enclosing.empty()
            || (where.enclosing.front() != view.bounds && where.enclosing.front() != view_turned)) {
            continue;
        }
        // Read the node as the reference walk does, or else as the first walk that visits it.
        handle_t read_as = graph.get_handle(c.node);
        bool seen = false;
        for (int pass = 0; pass < 2 && !seen; ++pass) {
            for (size_t w = 0; w < walks.size() && !seen; ++w) {
                if ((pass == 0) != ((int)w == ref_trav_idx)) {
                    continue;
                }
                for (const handle_t& h : walks[w]) {
                    if (graph.get_id(h) == c.node) {
                        read_as = h;
                        seen = true;
                        break;
                    }
                }
            }
        }
        HeldEdit edit;
        if (!edit_site_of(c, read_as, where, edit.site)) {
            parts.edit_outcomes->add(unplaced_outcome(c, "walk_leaves_site"));
            continue;
        }
        edit.candidate = c;
        edit.child = ChildSite{view.net, edit.site.bounds, where.chain};
        held.push_back(std::move(edit));
    }
    return held;
}

void TreeGenotyper::free_nodes_along_reference(nid_t bound, vector<nid_t>& out) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    // At most this many nodes each way, which bounds the work per site.
    const size_t limit = 64;
    graph.for_each_step_on_handle(graph.get_handle(bound), [&](const step_handle_t& step) {
        if (!parts.ref_offsets->count(graph.get_path_name(graph.get_path_handle_of_step(step)))) {
            return true;
        }
        for (bool forward : {true, false}) {
            step_handle_t at = step;
            for (size_t k = 0; k < limit; ++k) {
                if (forward ? !graph.has_next_step(at) : !graph.has_previous_step(at)) {
                    break;
                }
                at = forward ? graph.get_next_step(at) : graph.get_previous_step(at);
                const nid_t node = graph.get_id(graph.get_handle_of_step(at));
                if (!parts.child_placer->placement_of(node).enclosing.empty()) {
                    break;   // into a site
                }
                out.push_back(node);
            }
        }
        return false;   // one reference step is enough
    });
}

void TreeGenotyper::genotype_top_level_edits(vector<EditCandidate> candidates) {
    vector<EditCandidate> left;
    for (EditCandidate& c : candidates) {
        if (!free_edit_nodes.count(c.node)) {
            left.push_back(std::move(c));
        }
    }
    // By node, since each node's are genotyped together.
    vector<vector<EditCandidate>> by_node;
    for (EditCandidate& c : admitted(std::move(left))) {
        if (by_node.empty() || by_node.back().front().node != c.node) {
            by_node.emplace_back();
        }
        by_node.back().push_back(std::move(c));
    }
    after_pass_nodes = by_node.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < by_node.size(); ++i) {
        genotype_free_edits(by_node[i]);
    }
}

void TreeGenotyper::genotype_free_edits(vector<EditCandidate> candidates) {
    if (candidates.empty()) {
        return;
    }
    {
        // Once per node: a node starts one site, but the site may be visited again as the child
        // of a failed site under RecurseOnFail.
        std::lock_guard<std::mutex> lock(free_edit_mutex);
        if (!free_edit_nodes.insert(candidates.front().node).second) {
            return;
        }
    }
    const PathPositionHandleGraph& graph = *parts.graph;
    const vector<EditCandidate> chosen = admitted(std::move(candidates));
    for (size_t i = 0; i < chosen.size(); ++i) {
        const EditCandidate& c = chosen[i];
        NodePlacement where = parts.child_placer->placement_of(c.node);
        if (!where.enclosing.empty()) {
            // A site holds it, so its descent placed it, if anything did.
            continue;
        }
        // Read the node as a reference path does.
        handle_t read_as = graph.get_handle(c.node);
        graph.for_each_step_on_handle(read_as, [&](const step_handle_t& step) {
            const string name = graph.get_path_name(graph.get_path_handle_of_step(step));
            if (parts.ref_offsets->count(name)) {
                read_as = graph.get_handle_of_step(step);
                return false;
            }
            return true;
        });
        HeldEdit edit;
        if (!edit_site_of(c, read_as, where, edit.site)) {
            parts.edit_outcomes->add(unplaced_outcome(c, "walk_leaves_site"));
            continue;
        }
        edit.candidate = c;
        edit.child = ChildSite{where.site, edit.site.bounds, where.chain};
        genotype_edit_site(edit, {}, "", make_pair(0, 0), -1, NestingPlacement(), true);
    }
}

void TreeGenotyper::genotype_edit_site(const HeldEdit& edit, const vector<SiteBounds>& enclosing,
                                       const string& parent_ref_path_name,
                                       pair<size_t, size_t> parent_ref_interval,
                                       int ploidy_override, const NestingPlacement& placement,
                                       bool top_level) {
    const PathPositionHandleGraph& graph = *parts.graph;
    const SiteBounds& bounds = edit.site.bounds;
    const vector<Traversal>& walks = edit.site.travs;
    const AlleleSequences& sequences = edit.site.sequences;
    const EditCandidate& c = edit.candidate;
    EditOutcome outcome = unplaced_outcome(c, top_level ? "top" : "nested");
    outcome.site = edit_site_id(c.node, c.offset, c.kind,
                                c.kind == 'S'                      ? string(1, c.alt)
                                : c.kind == 'I' && !c.is_large() ? c.sequence
                                                                   : std::to_string(c.sequence.size()));
    const string& site_id = outcome.site;
    const size_t record_key = record_key_of(site_id);

    // The reference path through the node, preferring the parent's, and the node's interval on
    // it. A node no reference path visits takes its parent's, as a snarl off the reference does.
    string ref_path_name;
    size_t node_start = 0;
    graph.for_each_step_on_handle(bounds.start, [&](const step_handle_t& step) {
        const string name = graph.get_path_name(graph.get_path_handle_of_step(step));
        if (parts.ref_offsets->count(name)
            && (ref_path_name.empty() || name == parent_ref_path_name)) {
            ref_path_name = name;
            node_start = graph.get_position_of_step(step);
        }
        return name != parent_ref_path_name;
    });
    const bool no_ref_position = ref_path_name.empty() || placement.no_reference;
    if (ref_path_name.empty()) {
        if (top_level) {
            outcome.placed = "no_reference";
            parts.edit_outcomes->add(std::move(outcome));
            return;
        }
        ref_path_name = parent_ref_path_name;
    }
    size_t walk_length = 0;
    for (const handle_t& h : walks.front()) {
        walk_length += graph.get_length(h);
    }
    const pair<size_t, size_t> ref_range =
        no_ref_position ? parent_ref_interval : make_pair(node_start, node_start + walk_length);
    const int ref_offset = ref_offset_of(*parts.ref_offsets, ref_path_name);
    const int region_ploidy =
        parts.ploidy_regions->ploidy_at(ref_path_name, ref_range.first, ref_offset,
                                        ref_ploidy_of(*parts.ref_ploidies, ref_path_name));
    const int ploidy = ploidy_override >= 0 ? ploidy_override : region_ploidy;

    Ploidies ploidies{.ploidy = ploidy};
    if (!top_level) {
        ploidies.region_ploidy = region_ploidy;
        ploidies.also_score_other = true;
    }
    auto called = parts.genotyper->genotype(bounds, walks, 0, ploidies, enclosing, ref_path_name,
                                            ref_range, &sequences);
    vector<int>& genotype = called.first;
    SiteScore* score = called.second.get();
    unique_ptr<SnarlCaller::CallInfo> call_info(std::move(called.second));
    outcome.genotype = genotype;
    outcome.gq = score != nullptr ? score->gq : 0.0;
    outcome.depth_ratio = score != nullptr ? score->depth_ratio : -1.0;
    const bool carries_alt = std::count(genotype.begin(), genotype.end(), 1) > 0;
    const double max_dr = options.edit_bar.max_depth_ratio;
    if (score == nullptr || genotype.size() != (size_t)ploidy || !carries_alt
        || outcome.gq < options.edit_bar.min_gq
        || (max_dr > 0.0 && outcome.depth_ratio > max_dr)) {
        // Called without its ALT, not confidently, or from more reads than the called genotype
        // accounts for: nothing to stage.
        parts.edit_outcomes->add(std::move(outcome));
        return;
    }
    outcome.staged = true;

    vector<int> site_panel;
    bool site_panel_set = false;
    int64_t position_from_parent = 0;
    unique_ptr<StagedSite> staged;
    if (top_level) {
        site_panel_set = parts.linker->add(bounds, walks, genotype, score, 0, ref_path_name,
                                           ref_offset, record_key, placement, false, 0,
                                           &site_panel, &sequences);
        staged = stage_render_record(bounds, site_id, genotype, 0, call_info, score,
                                     ref_path_name, ref_offset, ploidy);
    } else {
        // As a nested snarl is: recorded unless no called parent allele reaches it yet, and with
        // a stand-in position where its node is off the reference.
        if (no_ref_position) {
            position_from_parent =
                base_path_position(ref_path_name, parent_ref_interval.first + ref_offset)
                + (int64_t)placement.parent_offset;
            site_panel_set = parts.linker->add(bounds, walks, genotype, score, 0, ref_path_name,
                                               ref_offset, record_key, placement, true,
                                               position_from_parent, &site_panel, &sequences);
        } else if (!placement.retain_only) {
            site_panel_set = parts.linker->add(bounds, walks, genotype, score, 0, ref_path_name,
                                               ref_offset, record_key, placement, false, 0,
                                               &site_panel, &sequences);
        }
        staged.reset(new StagedSite());
        staged->bounds = bounds;
        staged->ref_path_name = ref_path_name;
        staged->ref_offset = ref_offset;
        staged->ref_trav_idx = 0;
        staged->genotype = genotype;
        staged->ploidy = ploidy;
        staged->id = site_id;
        staged->record_key = record_key;
        staged->parent_record_key = placement.parent_record_key;
        staged->parent_crossing = placement.parent_crossing;
        staged->chain_key = placement.chain_key;
        staged->no_reference = no_ref_position;
        staged->reported_inline = placement.reported_inline;
        staged->position_from_parent = position_from_parent;
        staged->chain_offset = placement.parent_offset;
        staged->crossing_known = placement.crossing_known;
        staged->level = (uint8_t)min(placement.level, (size_t)255);
        staged->set_call(std::move(call_info), score);
    }
    // A leaf, in the chain its node is in.
    staged->off_panel = true;
    staged->leaf = true;
    staged->in_chain = true;
    staged->chain = edit.child.chain;
    staged->enclosing = enclosing;
    staged->travs = walks;
    staged->sequences = sequences;
    if (site_panel_set) {
        staged->panel_cache = std::move(site_panel);
        staged->panel_cached = true;
    }
    if (top_level) {
        parts.staged_sites->add_top_level(std::move(*staged));
    } else {
        parts.staged_sites->add_nested(std::move(*staged));
    }
    parts.edit_outcomes->add(std::move(outcome));
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
        const SiteBounds& bounds, const string& site_id, const vector<int>& trav_genotype,
        int ref_trav_idx,
        unique_ptr<SnarlCaller::CallInfo>& call_info, SiteScore* score,
        const string& ref_path_name, int ref_offset, int ploidy) const {
    unique_ptr<StagedSite> rec(new StagedSite());
    rec->bounds = bounds;
    rec->ref_path_name = ref_path_name;
    rec->ref_offset = ref_offset;
    rec->ref_trav_idx = ref_trav_idx;
    rec->genotype = trav_genotype;
    rec->ploidy = ploidy;
    rec->id = site_id;
    rec->record_key = record_key_of(site_id);
    rec->level = 0;
    rec->set_call(std::move(call_info), score);
    // `travs` is not moved here: descent runs after the emit and reads `travs` to find which
    // children the called alleles reach. The caller completes the record after descent.
    return rec;
}

void TreeGenotyper::fill_tree_fields(const SiteView& view, const SiteChildren& children,
                                     StagedSite& site) const {
    site.children = children;
    // A site holds a child site wherever one of its child chains has two nodes.
    site.leaf = children.entries.empty();
    site.in_chain = true;
    site.chain = parts.child_placer->chain_of(view.net);
    site.enclosing = view.enclosing;
}

}
}
