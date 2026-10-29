#include <algorithm>
#include <limits>

#include "tree_genotyper.hpp"
#include "../src/graph_caller.hpp"
#include "site_values.hpp"
#include "vcf_record.hpp"

//#define debug

namespace vg::multipass {

void TreeGenotyper::configure(const Parts& parts, const Options& options) {
    this->parts = parts;
    this->options = options;
    placer.configure(EditSitePlacer::Parts{.graph = parts.graph,
                                           .descent_policy = parts.descent_policy,
                                           .edit_source = parts.edit_source,
                                           .edit_outcomes = parts.edit_outcomes,
                                           .ref_paths = parts.ref_paths,
                                           .sites = parts.sites},
                     options.edit_bar);
}

namespace {

/**
 * A hold on the reads of a read source for one subtree, while the hold is in scope (see
 * `SiteReadSource::hold_reads`).
 *
 * A site too wide for one window of the read cache keeps the reads of its own fetch. The sites
 * below it that are as wide are genotyped from them, rather than by fetching them again. With
 * off-panel detection, so are the edit sites on the nodes of the candidates the call bar admits.
 */
struct SubtreeReadHold {
    /// The source holding the reads, or null for none.
    const SiteReadSource* source;
    SubtreeReadHold(const SubtreeReadHold&) = delete;
    SubtreeReadHold& operator=(const SubtreeReadHold&) = delete;
    SubtreeReadHold(const SiteReadSource* source, const EditCallBar& bar) : source(source) {
        if (source != nullptr) {
            source->hold_reads(bar);
        }
    }
    ~SubtreeReadHold() {
        if (source != nullptr) {
            source->release_reads();
        }
    }
};

/**
 * Get the boundary nodes of the sites in `sites` that enclose the site or edit `net`, innermost
 * site first.
 *
 * They place a site's depth-rate window when the site's own bounds are off the reference. A nested
 * site's bounds often are. Its ancestors' are not, and an ancestor lies within a window's width of
 * the site unless it is very large.
 */
vector<nid_t> fallback_nodes_of(const EditSnarlDecompositionOverlay& sites,
                                const HandleGraph& graph, const net_handle_t& net) {
    vector<nid_t> nodes;
    for (net_handle_t here = enclosing_site(sites, net); !sites.is_root(here);
         here = enclosing_site(sites, here)) {
        const auto [start, end] = oriented_ends(sites, graph, here);
        nodes.push_back(graph.get_id(start));
        nodes.push_back(graph.get_id(end));
    }
    return nodes;
}

/// Write out the site or edit `net` of `sites`, read the way `net` reads, with its alleles
/// `alleles`.
MaterializedSite materialize(const EditSnarlDecompositionOverlay& sites, const net_handle_t& net,
                             const CandidateAlleles& alleles) {
    return MaterializedSite{sites.get_site_start(net), sites.get_site_start_offset(net),
                            sites.get_site_end(net), sites.get_site_end_offset(net), alleles};
}

} // namespace

bool TreeGenotyper::genotype(const net_handle_t& site) {
    // Genotype the sites below, edit sites among them, while the hold has the reads.
    SubtreeReadHold hold(parts.read_source != nullptr ? parts.read_source : parts.edit_source,
                         options.edit_bar);
    const bool genotyped = genotype_tree(site, "", make_pair(0, 0), nullptr, -1, NestingPlacement());
    if (parts.edit_source != nullptr && parts.sites->is_root(enclosing_site(*parts.sites, site))) {
        // Genotype the free edit sites on the boundary nodes of this site, whose reads its own read
        // query fetched. The first site to reach a node between two sites takes it.
        const auto [start, end] = oriented_ends(*parts.sites, *parts.graph, site);
        const nid_t first = std::min(parts.graph->get_id(start), parts.graph->get_id(end));
        const nid_t last = std::max(parts.graph->get_id(start), parts.graph->get_id(end));
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
        // Then genotype those on the nodes that no site holds along the reference from either
        // bound to the next site, as far as their reads come with the reads of the lower bound.
        vector<nid_t> near;
        for (nid_t bound : {first, last}) {
            placer.free_nodes_along_reference(bound, near);
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

/// One site as the direct pass genotypes it.
///
/// It holds the site's candidates, its genotype, and its staged site. The staged site waits here
/// until the sites below have read the candidates, and then takes them.
struct TreeGenotyper::SiteCall {
    /// The site as the walk handed it over.
    net_handle_t view;
    /// The place of the site in the nesting tree.
    const NestingPlacement& placement;
    /// The site oriented along its reference path, and its candidate walks.
    CandidateFinder::WalkSite site;
    /// The site's net handle in the overlay of sites, read the way `site` reads.
    net_handle_t net;
    /// The ID of the site's records (see `StagedSite::id`).
    string id;
    /// The record key of `id`.
    size_t record_key = 0;
    /// The sites nested in the site.
    DescentPolicy::Nested nested;
    /// The candidate traversals as walks, once the candidates are complete, until the staged site
    /// takes them.
    vector<Traversal> walks;
    /// The ploidy of the region the site is in, from the contig or the region BED.
    int region_ploidy = 0;
    /// The ploidy the site is genotyped at.
    int ploidy = 0;
    /// The genotype, as indices into `walks`.
    vector<int> genotype;
    /// The staged site, or null if the site is not staged.
    unique_ptr<StagedSite> staged;
    /// Whether `staged` goes to the nested list rather than a render queue.
    bool staged_nested = false;

    /// Get where the site lies on its reference path, as a range of positions.
    pair<size_t, size_t> ref_range() const {
        return make_pair(get<0>(site.ref_interval), get<1>(site.ref_interval));
    }
};

bool TreeGenotyper::genotype_tree(const net_handle_t& view, const string& parent_ref_path_name,
                                  pair<size_t, size_t> parent_ref_interval,
                                  const vector<WalkSet>* parent_child_walk_sets,
                                  int ploidy_override, const NestingPlacement& placement) {
    SiteCall call{view, placement};
    const auto [view_start, view_end] = oriented_ends(*parts.sites, *parts.graph, view);
    if (!parts.candidates->find(view_start, view_end, parent_ref_path_name,
                                parent_ref_interval, parent_child_walk_sets,
                                placement.no_reference, call.site)) {
        return false;
    }
    call.net = parts.sites->get_site_start(view) == call.site.start ? view
                                                                   : parts.sites->flip(view);
    // Name the site from its bounds as the candidate finder oriented them.
    call.id = parts.site_id(call.site.start, call.site.end);
    call.record_key = record_key_of(call.id);
    call.nested = parts.descent_policy->nested_in(view, call.site.start, call.site.end);

    // A ploidy from the parent overrides the contig's or the region BED's. The region's is still
    // the number of the sample's haplotypes here, which the depth term needs.
    const string& ref_path_name = call.site.ref_path_name;
    call.region_ploidy =
        parts.ploidy_regions->ploidy_at(ref_path_name, get<0>(call.site.ref_interval),
                                        parts.ref_paths->offset_of(ref_path_name),
                                        parts.ref_paths->ploidy_of(ref_path_name));
    call.ploidy = ploidy_override >= 0 ? ploidy_override : call.region_ploidy;
    if (call.ploidy == 0) {
        // The sample has no copy of this region, so there is nothing to genotype and no record.
        // A parent always passes at least 1, so the ploidy came from the region. Report success,
        // so that the walk does not call the children in the site's place.
        return true;
    }

    bool genotyped;
    if (parent_child_walk_sets != nullptr && !parent_child_walk_sets->empty()) {
        genotyped = genotype_in_parent_sets(call, *parent_child_walk_sets);
    } else if (ploidy_override >= 0) {
        genotyped = genotype_nested_chain(call);
    } else {
        genotyped = genotype_top_level_site(call);
    }

    // Descend only from a successful call, since a failed site has no genotype to take a child's
    // ploidy from. RecurseOnFail calls the children of a failed top-level site as top-level sites,
    // but the children of a failed nested site are not called.
    if (genotyped && options.nested_calling && !call.genotype.empty()
        && parent_child_walk_sets == nullptr) {
        descend_nested(call);
    }
    if (options.top_down && !call.genotype.empty() && call.nested.children.known) {
        descend_top_down(call);
    }
    file_staged(call);
    return genotyped;
}

/// Add each walk in the parent's walk sets that `walks` does not already hold.
static void merge_parent_sets(vector<Traversal>& walks, const vector<WalkSet>& sets) {
    for (const WalkSet& set : sets) {
        for (const Traversal& walk : set) {
            if (std::find(walks.begin(), walks.end(), walk) == walks.end()) {
                walks.push_back(walk);
            }
        }
    }
}

bool TreeGenotyper::genotype_in_parent_sets(SiteCall& call, const vector<WalkSet>& sets) {
    // The child has one strand per parent allele.
    call.ploidy = sets.size();
    merge_parent_sets(call.site.walks, sets);
    call.walks = call.site.walks;

    // A parent allele with an empty traversal set skips the child, and gets a star or missing
    // allele rather than a genotype.
    vector<int> traversing_sets;
    for (int set_idx = 0; set_idx < call.ploidy; ++set_idx) {
        if (!sets[set_idx].empty()) {
            traversing_sets.push_back(set_idx);
        }
    }

    unique_ptr<SnarlCaller::CallInfo> trav_call_info;
    SiteScore* trav_score = nullptr;
    const int marker = options.star_allele ? STAR_ALLELE_MARKER : MISSING_ALLELE_MARKER;
    call.genotype.assign(call.ploidy, marker);
    if (!traversing_sets.empty()) {
        // Genotype at the number of strands that cross the site, since a site only one strand
        // crosses is not diploid. Then put the called alleles, a sorted multiset, on those
        // strands.
        const auto effective_ploidy = (int)traversing_sets.size();
        vector<int> called_alleles;
        std::tie(called_alleles, trav_call_info) = genotype_site(
            call.site.start, call.site.end, call.walks, call.site.ref_trav_idx,
            Ploidies{.ploidy = effective_ploidy, .region_ploidy = call.region_ploidy},
            call.net, trav_score);
        for (size_t j = 0; j < traversing_sets.size() && j < called_alleles.size(); ++j) {
            call.genotype[traversing_sets[j]] = called_alleles[j];
        }
    }

    // Stage the site only if it is on a reference path.
    if (!call.site.use_parent_interval) {
        record_and_stage(call, trav_call_info, trav_score);
    }
    return call.genotype.size() == (size_t)call.ploidy;
}

/// Give a staged site its place in the nesting tree, with the level capped at 255.
static void place_in_tree(StagedSite& staged, const NestingPlacement& placement) {
    staged.placement = placement;
    staged.placement.level = min(placement.level, (size_t)255);
}

bool TreeGenotyper::genotype_nested_chain(SiteCall& call) {
    call.walks = call.site.walks;
    unique_ptr<SnarlCaller::CallInfo> trav_call_info;
    SiteScore* trav_score = nullptr;
    // Score the other ploidy too, since the linkage pass can revise the ploidy of a nested chain.
    std::tie(call.genotype, trav_call_info) = genotype_site(
        call.site.start, call.site.end, call.walks, call.site.ref_trav_idx,
        Ploidies{.ploidy = call.ploidy, .region_ploidy = call.region_ploidy,
                 .also_score_other = true},
        call.net, trav_score);
    assert(call.genotype.empty() || call.genotype.size() == call.ploidy);

    const NestingPlacement& placement = call.placement;
    // The boundaries of a site that took its parent's interval are on no reference path.
    const bool no_ref_position = call.site.use_parent_interval;
    const int64_t ref_offset = parts.ref_paths->offset_of(call.site.ref_path_name);
    // Stand in for the missing position with the parent's position, plus the chain's offset along
    // the parent.
    const int64_t position_from_parent =
        no_ref_position ? base_path_position(call.site.ref_path_name,
                                             get<0>(call.site.ref_interval) + ref_offset)
                              + (int64_t)placement.parent_offset
                        : 0;

    call.staged = stage_render_record(
        call.net, call.id, call.genotype, std::move(call.walks),
        call.site.ref_trav_idx,
        trav_call_info, trav_score,
        ReferencePlacement{.path_name = call.site.ref_path_name,
                           .offset = ref_offset,
                           .no_reference = no_ref_position,
                           .position_from_parent = position_from_parent},
        call.ploidy);
    place_in_tree(*call.staged, placement);
    call.staged_nested = true;
    record_nested_chain(call);
    return call.genotype.size() == (size_t)call.ploidy;
}

void TreeGenotyper::record_nested_chain(SiteCall& call) const {
    const NestingPlacement& placement = call.placement;
    if (call.site.use_parent_interval) {
        // The site has no reference position, so it is recorded but never written. This is
        // tested before retain_only, which would not record it.
        parts.linker->file(*call.staged, placement.one_copy);
        ++parts.descent_counters->no_ref_recorded;
        int copies = 0;
        for (int a : call.genotype) {
            copies += a >= 0 ? 1 : 0;
        }
        parts.descent_counters->no_ref_copies[copies < 3 ? copies : 2].fetch_add(1);
    } else if (!placement.retain_only) {
        // Record the site. With `reported_inline`, an enclosing block's ALT already spells the
        // chain, so it gets no line, but its allele pair still phases everything inside it.
        // Otherwise the line is written after the linkage pass, from the chosen genotype.
        parts.linker->file(*call.staged, placement.one_copy);
    }
    // A retain_only chain is crossed by no called parent allele, so it is not recorded yet. It is
    // kept, since the linkage model may still move the parent onto an allele that crosses it.
}

bool TreeGenotyper::genotype_top_level_site(SiteCall& call) {
    call.walks = call.site.walks;
    unique_ptr<SnarlCaller::CallInfo> trav_call_info;
    SiteScore* trav_score = nullptr;
    std::tie(call.genotype, trav_call_info) = genotype_site(
        call.site.start, call.site.end, call.walks, call.site.ref_trav_idx,
        Ploidies{.ploidy = call.ploidy}, call.net, trav_score);
    assert(call.genotype.empty() || call.genotype.size() == call.ploidy);
    record_and_stage(call, trav_call_info, trav_score);
    return call.genotype.size() == (size_t)call.ploidy;
}

void TreeGenotyper::record_and_stage(SiteCall& call,
                                     unique_ptr<SnarlCaller::CallInfo>& trav_call_info,
                                     SiteScore* trav_score) const {
    call.staged = stage_render_record(
        call.net, call.id, call.genotype, std::move(call.walks),
        call.site.ref_trav_idx,
        trav_call_info, trav_score,
        ReferencePlacement{.path_name = call.site.ref_path_name,
                           .offset = (int64_t)parts.ref_paths->offset_of(call.site.ref_path_name)},
        call.ploidy);
    call.staged_nested = false;
    // TODO: is this right? A child of top-down calling takes the placement of its parent, which
    // under nested calling can be a nested chain. The child is filed with that placement, but
    // staged and rendered as a top-level site, with none.
    place_in_tree(*call.staged, call.placement);
    parts.linker->file(*call.staged, call.placement.one_copy);
    call.staged->placement = NestingPlacement();
}

DescentPolicy::Nested TreeGenotyper::with_held_edits(const DescentPolicy::Nested& nested,
                                                   const vector<HeldEdit>& held) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    DescentPolicy::Nested out;
    out.children = nested.children;
    size_t e = 0;
    for (const net_handle_t& child : nested.sites) {
        const nid_t child_first = std::min(graph.get_id(parts.sites->get_site_start(child)),
                                           graph.get_id(parts.sites->get_site_end(child)));
        for (; e < held.size() && held[e].candidate.node < child_first; ++e) {
            out.sites.push_back(held[e].net);
        }
        out.sites.push_back(child);
    }
    for (; e < held.size(); ++e) {
        out.sites.push_back(held[e].net);
    }
    return out;
}

void TreeGenotyper::descend_nested(const SiteCall& call) {
    // The edit sites this site holds itself, rather than through a child site, are children of
    // it too, each in the chain its node is in.
    vector<HeldEdit> held;
    // A site that descends is staged, so its staged site holds its walks.
    const vector<Traversal>& walks = call.staged->alleles.walks;
    if (parts.edit_source != nullptr) {
        held = placer.held_by(call.view, walks, call.site.ref_trav_idx);
    }
    DescentPolicy::Nested with_edits;
    if (!held.empty()) {
        with_edits = with_held_edits(call.nested, held);
    }
    vector<bool> edit_placed(held.size(), false);
    // Keep a child no called allele crosses only where the linkage pass can come back to it,
    // which needs the linkage model.
    const DescentPolicy::Parent parent{.site = call.view,
                                     .nested = held.empty() ? call.nested : with_edits,
                                     .record_key = call.record_key,
                                     .travs = walks,
                                     .genotype = call.genotype,
                                     .ref_trav_idx = call.site.ref_trav_idx,
                                     .ploidy = call.ploidy,
                                     .placement = call.placement};
    parts.descent_policy->place(
        parent, options.off_reference, parts.linker->enabled(),
        [&](const DescentPolicy::Placed& child) {
            genotype_placed_child(call, child, held, edit_placed);
        });
    for (size_t i = 0; i < held.size(); ++i) {
        if (!edit_placed[i]) {
            // The placer left the edit site out, as it would any child: off the reference without
            // off-reference descent, or crossed by no called allele without the linkage model.
            parts.edit_outcomes->add(unplaced_outcome(held[i].candidate, "not_placed"));
        }
    }
}

void TreeGenotyper::genotype_placed_child(const SiteCall& call, const DescentPolicy::Placed& child,
                                          const vector<HeldEdit>& held,
                                          vector<bool>& edit_placed) {
    if (child.placement.level < 16) {
        ++parts.descent_counters->depth_hist[child.placement.level];
    }
    if (parts.sites->is_edit(child.site)) {
        // The child is one of the held edit sites.
        for (size_t i = 0; i < held.size(); ++i) {
            if (held[i].net == child.site) {
                edit_placed[i] = true;
                genotype_edit_site(held[i], call.site.ref_path_name, call.ref_range(),
                                   child.ploidy, child.placement, false);
                return;
            }
        }
        return;
    }
    SubtreeReadHold hold(parts.edit_source, options.edit_bar);
    genotype_tree(child.site, call.site.ref_path_name, call.ref_range(), nullptr, child.ploidy,
                  child.placement);
}

void TreeGenotyper::descend_top_down(const SiteCall& call) {
    const vector<Traversal>& travs = call.site.walks;
    for (const net_handle_t& child : call.nested.sites) {
        const auto [child_start, child_end] = oriented_ends(*parts.sites, *parts.graph, child);
        // Find, for each parent allele, the traversals of the child that agree with it. A star
        // or missing allele gets an empty set.
        vector<WalkSet> child_trav_sets;
        bool any_real_traversals = false;
        for (int allele_idx : call.genotype) {
            if (allele_idx >= 0 && allele_idx < travs.size()) {
                WalkSet tset = parts.candidates->find_child_walk_set(
                    travs[allele_idx], child_start, child_end);
                any_real_traversals = any_real_traversals || !tset.empty();
                child_trav_sets.push_back(std::move(tset));
            } else {
                child_trav_sets.emplace_back();
            }
        }
        if (any_real_traversals) {
            genotype_tree(child, call.site.ref_path_name,
                          call.ref_range(), &child_trav_sets, -1, call.placement);
        }
    }
}

void TreeGenotyper::file_staged(SiteCall& call) {
    if (call.staged == nullptr) {
        return;
    }
    call.staged->children = call.nested.children;
    if (call.staged_nested) {
        parts.staged_sites->add_nested(std::move(*call.staged));
    } else {
        parts.staged_sites->add_top_level(std::move(*call.staged));
    }
    call.staged.reset();
}

/// Get the ALT part of the ID of an edit candidate's site.
///
/// That is the base of a substitution, the inserted sequence of a small insertion, and the length
/// of any other change.
static string alt_label(const EditCandidate& c) {
    if (c.kind == 'S') {
        return string(1, c.alt);
    }
    if (c.kind == 'I' && !c.is_large()) {
        return c.sequence;
    }
    return std::to_string(c.sequence.size());
}

string TreeGenotyper::reference_through(handle_t node, const string& preferred,
                                        size_t& node_start) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    string ref_path_name;
    graph.for_each_step_on_handle(node, [&](const step_handle_t& step) {
        const string name = graph.get_path_name(graph.get_path_handle_of_step(step));
        if (parts.ref_paths->contains(name) && (ref_path_name.empty() || name == preferred)) {
            ref_path_name = name;
            node_start = graph.get_position_of_step(step);
        }
        return name != preferred;
    });
    return ref_path_name;
}

void TreeGenotyper::genotype_top_level_edits(vector<EditCandidate> candidates) {
    vector<vector<EditCandidate>> by_node = placer.group_unclaimed(std::move(candidates));
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < by_node.size(); ++i) {
        genotype_free_edits(by_node[i]);
    }
}

void TreeGenotyper::genotype_free_edits(vector<EditCandidate> candidates) {
    placer.place_free_edits(std::move(candidates), [&](const HeldEdit& edit) {
        genotype_edit_site(edit, "", make_pair(0, 0), -1, NestingPlacement(), true);
    });
}

void TreeGenotyper::genotype_edit_site(const HeldEdit& edit,
                                       const string& parent_ref_path_name,
                                       pair<size_t, size_t> parent_ref_interval,
                                       int ploidy_override, const NestingPlacement& placement,
                                       bool top_level) {
    const PathPositionHandleGraph& graph = *parts.graph;
    const vector<Traversal>& walks = edit.alleles.walks;
    const EditCandidate& c = edit.candidate;
    EditOutcome outcome = unplaced_outcome(c, top_level ? "top" : "nested");
    outcome.site = edit_site_id(c.node, c.offset, c.kind, alt_label(c));
    const string& site_id = outcome.site;

    // A node no reference path visits takes its parent's path, as a site off the reference does.
    size_t node_start = 0;
    string ref_path_name =
        reference_through(parts.sites->get_site_start(edit.net), parent_ref_path_name, node_start);
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
    const int64_t ref_offset = parts.ref_paths->offset_of(ref_path_name);
    const int region_ploidy =
        parts.ploidy_regions->ploidy_at(ref_path_name, ref_range.first, ref_offset,
                                        parts.ref_paths->ploidy_of(ref_path_name));
    const int ploidy = ploidy_override >= 0 ? ploidy_override : region_ploidy;
    if (ploidy == 0) {
        // The sample has no copy of this region, so the edit is not genotyped.
        outcome.placed = "ploidy_zero";
        parts.edit_outcomes->add(std::move(outcome));
        return;
    }

    Ploidies ploidies{.ploidy = ploidy};
    if (!top_level) {
        ploidies.region_ploidy = region_ploidy;
        ploidies.also_score_other = true;
    }
    auto called = parts.genotyper->genotype_at(
        materialize(*parts.sites, edit.net, edit.alleles), ploidies,
        fallback_nodes_of(*parts.sites, graph, edit.net));
    vector<int>& genotype = called.first;
    SiteScore* score = called.second.get();
    unique_ptr<SnarlCaller::CallInfo> call_info(std::move(called.second));
    outcome.genotype = genotype;
    outcome.gq = score != nullptr ? score->gq : 0.0;
    outcome.depth_ratio = score != nullptr ? score->depth_ratio : -1.0;
    const bool carries_alt =
        std::any_of(genotype.begin(), genotype.end(), [](int a) { return a > 0; });
    const double max_dr = options.edit_bar.max_depth_ratio;
    if (score == nullptr || genotype.size() != (size_t)ploidy || !carries_alt
        || outcome.gq < options.edit_bar.min_gq
        || (max_dr > 0.0 && outcome.depth_ratio > max_dr)) {
        // Stage nothing for a call without its ALT, a call without confidence, or a call from
        // more reads than the called genotype accounts for.
        parts.edit_outcomes->add(std::move(outcome));
        return;
    }
    outcome.staged = true;

    // A top-level edit site always has a reference position, so only a nested one stands in at
    // the position of its parent.
    const int64_t position_from_parent =
        no_ref_position ? base_path_position(ref_path_name, parent_ref_interval.first + ref_offset)
                              + (int64_t)placement.parent_offset
                        : 0;
    unique_ptr<StagedSite> staged = stage_render_record(
        edit.net, site_id, genotype, walks, 0, call_info, score,
        ReferencePlacement{.path_name = ref_path_name,
                           .offset = ref_offset,
                           .no_reference = no_ref_position,
                           .position_from_parent = position_from_parent},
        ploidy);
    if (!top_level) {
        place_in_tree(*staged, placement);
    }
    staged->alleles.edits = edit.alleles.edits;
    // A retain_only site is crossed by no called parent allele, so it is not filed yet. A site
    // with no reference position is filed anyway.
    if (top_level || no_ref_position || !placement.retain_only) {
        parts.linker->file(*staged, placement.one_copy);
    }
    if (top_level) {
        parts.staged_sites->add_top_level(std::move(*staged));
    } else {
        parts.staged_sites->add_nested(std::move(*staged));
    }
    parts.edit_outcomes->add(std::move(outcome));
}

pair<vector<int>, unique_ptr<SnarlCaller::CallInfo>> TreeGenotyper::genotype_site(
    const handle_t& start, const handle_t& end, const vector<Traversal>& travs, int ref_trav_idx,
    const Ploidies& ploidies, const net_handle_t& net, SiteScore*& score) const {
    const CandidateAlleles alleles{travs, ref_trav_idx, {}};
    // A snarl's bounds lie on the inward sides of its boundary nodes.
    const MaterializedSite site{start, parts.graph->get_length(start), end, 0, alleles};
    auto called = parts.genotyper->genotype_at(site, ploidies,
                                               fallback_nodes_of(*parts.sites, *parts.graph, net));
    score = called.second.get();
    return make_pair(std::move(called.first),
                     unique_ptr<SnarlCaller::CallInfo>(std::move(called.second)));
}

unique_ptr<StagedSite> TreeGenotyper::stage_render_record(
        const net_handle_t& net, const string& site_id,
        const vector<int>& trav_genotype,
        vector<Traversal> travs, int ref_trav_idx, unique_ptr<SnarlCaller::CallInfo>& call_info,
        SiteScore* score, const ReferencePlacement& reference, int ploidy) const {
    auto rec = std::make_unique<StagedSite>();
    rec->net = net;
    rec->reference = reference;
    rec->alleles.walks = std::move(travs);
    rec->alleles.ref_allele = ref_trav_idx;
    rec->genotype = trav_genotype;
    rec->ploidy = ploidy;
    rec->id = site_id;
    rec->record_key = record_key_of(site_id);
    // Keep the call info, which the render reads to map the written alleles back to matrix
    // columns, index GL and compute QUAL.
    rec->set_call(std::move(call_info), score);
    return rec;
}

}
