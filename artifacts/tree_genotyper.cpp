#include <algorithm>
#include <limits>

#include "tree_genotyper.hpp"
#include "graph_caller.hpp"
#include "site_values.hpp"
#include "vcf_record.hpp"

//#define debug

namespace vg::multipass {

/// Make the outcome of a candidate for the dump, saying where it was placed, or why it was not.
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

} // namespace

bool TreeGenotyper::genotype(const SiteView& site) {
    // Genotype the sites below, edit sites among them, while the hold has the reads.
    SubtreeReadHold hold(parts.read_source != nullptr ? parts.read_source : parts.edit_source,
                         options.edit_bar);
    const bool genotyped = genotype_tree(site, "", make_pair(0, 0), nullptr, -1, NestingPlacement());
    if (parts.edit_source != nullptr && site.enclosing.empty()) {
        // Genotype the free edit sites on the boundary nodes of this site, whose reads its own read
        // query fetched. The first site to reach a node between two sites takes it.
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
        // Then genotype those on the nodes that no site holds along the reference from either
        // bound to the next site, as far as their reads come with the reads of the lower bound.
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

/// One site as the direct pass genotypes it.
///
/// It holds the site's candidates, its genotype, and its staged site. The staged site waits here
/// until the sites below have read the candidates, and then takes them.
struct TreeGenotyper::SiteCall {
    /// The site as the walk handed it over.
    const SiteView& view;
    /// The place of the site in the nesting tree.
    const NestingPlacement& placement;
    /// The site oriented along its reference path, and its candidate walks.
    CandidateFinder::WalkSite site;
    /// The bounds of `site`.
    SiteBounds bounds;
    /// The ID of the site's records (see `StagedSite::id`).
    string id;
    /// The record key of `id`.
    size_t record_key = 0;
    /// The sites nested in the site.
    ChildPlacer::Nested nested;
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

bool TreeGenotyper::genotype_tree(const SiteView& view, const string& parent_ref_path_name,
                                  pair<size_t, size_t> parent_ref_interval,
                                  const ChildWalkSets* parent_child_walk_sets,
                                  int ploidy_override, const NestingPlacement& placement) {
    SiteCall call{view, placement};
    if (!parts.candidates->find(view.bounds.start, view.bounds.end, parent_ref_path_name,
                                parent_ref_interval, parent_child_walk_sets,
                                placement.no_reference, call.site)) {
        return false;
    }
    call.bounds = SiteBounds{call.site.start, call.site.end};
    // Name the site from its bounds as the candidate finder oriented them.
    call.id = parts.site_id(call.bounds);
    call.record_key = record_key_of(call.id);
    call.nested = parts.child_placer->nested_in(view.net, call.bounds);

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
static void merge_parent_sets(vector<Traversal>& walks, const ChildWalkSets& sets) {
    for (const WalkSet& set : sets) {
        for (const Traversal& walk : set) {
            if (std::find(walks.begin(), walks.end(), walk) == walks.end()) {
                walks.push_back(walk);
            }
        }
    }
}

bool TreeGenotyper::genotype_in_parent_sets(SiteCall& call, const ChildWalkSets& sets) {
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
            call.bounds, call.walks, call.site.ref_trav_idx,
            Ploidies{.ploidy = effective_ploidy, .region_ploidy = call.region_ploidy},
            call.view.enclosing, trav_score);
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
        call.bounds, call.walks, call.site.ref_trav_idx,
        Ploidies{.ploidy = call.ploidy, .region_ploidy = call.region_ploidy,
                 .also_score_other = true},
        call.view.enclosing, trav_score);
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
        call.bounds, call.id, call.genotype, std::move(call.walks), call.site.ref_trav_idx,
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
        call.bounds, call.walks, call.site.ref_trav_idx, Ploidies{.ploidy = call.ploidy},
        call.view.enclosing, trav_score);
    assert(call.genotype.empty() || call.genotype.size() == call.ploidy);
    record_and_stage(call, trav_call_info, trav_score);
    return call.genotype.size() == (size_t)call.ploidy;
}

void TreeGenotyper::record_and_stage(SiteCall& call,
                                     unique_ptr<SnarlCaller::CallInfo>& trav_call_info,
                                     SiteScore* trav_score) const {
    call.staged = stage_render_record(
        call.bounds, call.id, call.genotype, std::move(call.walks), call.site.ref_trav_idx,
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

ChildPlacer::Nested TreeGenotyper::with_held_edits(const ChildPlacer::Nested& nested,
                                                   const vector<HeldEdit>& held) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    ChildPlacer::Nested out;
    out.children = nested.children;
    size_t e = 0;
    for (const ChildSite& child : nested.sites) {
        const nid_t child_first = std::min(graph.get_id(child.bounds.start),
                                           graph.get_id(child.bounds.end));
        for (; e < held.size() && held[e].candidate.node < child_first; ++e) {
            out.sites.push_back(held[e].child);
        }
        out.sites.push_back(child);
    }
    for (; e < held.size(); ++e) {
        out.sites.push_back(held[e].child);
    }
    return out;
}

void TreeGenotyper::descend_nested(const SiteCall& call) {
    // The edit sites this site holds itself, rather than through a child site, are children of
    // it too, each in the chain its node is in.
    vector<HeldEdit> held;
    // A site that descends is staged, so its staged site holds its walks.
    const vector<Traversal>& walks = call.staged->travs;
    if (parts.edit_source != nullptr) {
        held = edits_held_by(call.view, call.bounds, walks, call.site.ref_trav_idx);
    }
    ChildPlacer::Nested with_edits;
    if (!held.empty()) {
        with_edits = with_held_edits(call.nested, held);
    }
    vector<bool> edit_placed(held.size(), false);
    // Keep a child no called allele crosses only where the linkage pass can come back to it,
    // which needs the linkage model.
    const ChildPlacer::Parent parent{.site = call.view,
                                     .nested = held.empty() ? call.nested : with_edits,
                                     .record_key = call.record_key,
                                     .travs = walks,
                                     .genotype = call.genotype,
                                     .ref_trav_idx = call.site.ref_trav_idx,
                                     .ploidy = call.ploidy,
                                     .placement = call.placement};
    parts.child_placer->place(
        parent, options.off_reference, parts.linker->enabled(),
        [&](const ChildPlacer::Placed& child) {
            genotype_placed_child(call, child, held, edit_placed);
        });
    for (size_t i = 0; i < held.size(); ++i) {
        if (!edit_placed[i]) {
            // The placer left the edit site out, as it would any child: off the reference without
            // off-reference descent, or crossed by no called allele without the linkage model.
            parts.edit_outcomes->add(unplaced_outcome(held[i].candidate, "not_placed"));
        }
        parts.sites->remove_edit(held[i].child.net);
    }
}

void TreeGenotyper::genotype_placed_child(const SiteCall& call, const ChildPlacer::Placed& child,
                                          const vector<HeldEdit>& held,
                                          vector<bool>& edit_placed) {
    if (child.placement.level < 16) {
        ++parts.descent_counters->depth_hist[child.placement.level];
    }
    if (parts.sites->is_edit(child.site.net)) {
        // The child is one of the held edit sites.
        for (size_t i = 0; i < held.size(); ++i) {
            if (held[i].child.net == child.site.net) {
                edit_placed[i] = true;
                genotype_edit_site(held[i], child.site.enclosing, call.site.ref_path_name,
                                   call.ref_range(), child.ploidy, child.placement, false);
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
    for (const ChildSite& child : call.nested.sites) {
        // Find, for each parent allele, the traversals of the child that agree with it. A star
        // or missing allele gets an empty set.
        ChildWalkSets child_trav_sets;
        bool any_real_traversals = false;
        for (int allele_idx : call.genotype) {
            if (allele_idx >= 0 && allele_idx < travs.size()) {
                WalkSet tset = parts.candidates->find_child_walk_set(
                    travs[allele_idx], child.bounds.start, child.bounds.end);
                any_real_traversals = any_real_traversals || !tset.empty();
                child_trav_sets.push_back(std::move(tset));
            } else {
                child_trav_sets.emplace_back();
            }
        }
        if (any_real_traversals) {
            genotype_tree(ChildPlacer::view_of(call.view, child), call.site.ref_path_name,
                          call.ref_range(), &child_trav_sets, -1, call.placement);
        }
    }
}

void TreeGenotyper::file_staged(SiteCall& call) {
    if (call.staged == nullptr) {
        return;
    }
    fill_tree_fields(call.view, call.nested.children, *call.staged);
    if (call.staged_nested) {
        parts.staged_sites->add_nested(std::move(*call.staged));
    } else {
        parts.staged_sites->add_top_level(std::move(*call.staged));
    }
    call.staged.reset();
}

vector<EditCandidate> TreeGenotyper::admitted(const vector<EditCandidate>& candidates) const {
    // The candidates of each base, substitutions and indels separately, are consecutive. Keep one
    // of each.
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

/// Get the ALTs of an edit candidate, along the forward strand of its node.
static vector<string> alts_of(const EditCandidate& c) {
    switch (c.kind) {
    case 'M':
        return c.alts;
    case 'S':
        return {string(1, c.alt)};
    case 'I':
        return {c.sequence};
    default:
        return {""};
    }
}

/// Count the bases that each ALT of an edit candidate replaces.
static size_t replaced_length(const EditCandidate& c) {
    switch (c.kind) {
    case 'M':
        return c.span;
    case 'S':
        return 1;
    case 'D':
        return c.sequence.size();
    default:
        return 0;
    }
}

bool TreeGenotyper::walk_in_one_site(const EditCandidate& c, const NodePlacement& where) const {
    for (nid_t id : c.walk) {
        const NodePlacement other = parts.child_placer->placement_of(id);
        if (other.enclosing.size() != where.enclosing.size()
            || (!other.enclosing.empty() && other.enclosing.front() != where.enclosing.front())) {
            return false;
        }
    }
    return true;
}

bool TreeGenotyper::reference_reads_whole_walk(const EditCandidate& c) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    size_t on_reference = 0;
    // The position the next node must have on the reference, or -1 if the last node was not on
    // it.
    int64_t expected = -1;
    for (nid_t id : c.walk) {
        const handle_t h = graph.get_handle(id);
        size_t steps = 0;
        // The position of the node on the reference, or -2 if the reference reads it backward.
        int64_t position = -1;
        graph.for_each_step_on_handle(h, [&](const step_handle_t& step) {
            const path_handle_t path = graph.get_path_handle_of_step(step);
            if (parts.ref_paths->contains(graph.get_path_name(path))) {
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
    return on_reference == 0 || on_reference == c.walk.size();
}

net_handle_t TreeGenotyper::add_edit(const EditSite& site) const {
    return parts.sites->add_edit(site.change_start, site.change_start_offset, site.change_end,
                                 site.change_end_offset);
}

bool TreeGenotyper::edit_site_of(const EditCandidate& c, handle_t read_as,
                                 const NodePlacement& where, EditSite& out) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    vector<string> alts = alts_of(c);
    const size_t length = replaced_length(c);
    if (c.kind == 'S' || (c.kind == 'M' && c.walk.size() == 1)) {
        // The edit lies inside one node. Read it as `read_as` reads the node.
        const bool backward = graph.get_is_reverse(read_as);
        for (string& alt : alts) {
            alt = backward ? reverse_complement(alt) : alt;
        }
        out = edit_site(graph, {read_as}, backward ? c.offset + length - 1 : c.offset,
                        backward ? graph.get_length(read_as) - c.offset - length : c.offset, length,
                        alts);
        return true;
    }
    if (graph.get_is_reverse(read_as)) {
        return false;
    }
    // The record is spelled along the walk, and placed where a reference path reads the first
    // node of the walk.
    if (!walk_in_one_site(c, where) || !reference_reads_whole_walk(c)) {
        return false;
    }
    Traversal steps;
    for (nid_t id : c.walk) {
        steps.push_back(graph.get_handle(id));
    }
    out = edit_site(graph, steps, c.offset, c.kind == 'M' ? c.offset : c.offset + 1, length, alts);
    return true;
}

/// Find the handle a site's walks read node `node` by.
///
/// Prefer the reference walk, then the walks in order. If no walk visits the node, read it
/// forward.
static handle_t read_as_in_walks(const HandleGraph& graph, const vector<Traversal>& walks,
                                 int ref_trav_idx, nid_t node) {
    auto visit_in = [&](const Traversal& walk) -> const handle_t* {
        for (const handle_t& h : walk) {
            if (graph.get_id(h) == node) {
                return &h;
            }
        }
        return nullptr;
    };
    if (ref_trav_idx >= 0 && ref_trav_idx < (int)walks.size()) {
        if (const handle_t* h = visit_in(walks[ref_trav_idx])) {
            return *h;
        }
    }
    for (size_t w = 0; w < walks.size(); ++w) {
        if ((int)w == ref_trav_idx) {
            continue;
        }
        if (const handle_t* h = visit_in(walks[w])) {
            return *h;
        }
    }
    return graph.get_handle(node);
}

/// Get the bounds `bounds` read the other way.
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
    // Collect the nodes the walks visit.
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
    // Leave out the candidates on the bounds, which are in the parent's chain. The site's own read
    // query counted the rest.
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
    for (const EditCandidate& c : admitted(found)) {
        // Only the innermost site holding the node holds the edit.
        NodePlacement where = parts.child_placer->placement_of(c.node);
        if (where.enclosing.empty()
            || (where.enclosing.front() != view.bounds && where.enclosing.front() != view_turned)) {
            continue;
        }
        HeldEdit edit;
        if (!edit_site_of(c, read_as_in_walks(graph, walks, ref_trav_idx, c.node), where,
                          edit.site)) {
            parts.edit_outcomes->add(unplaced_outcome(c, "walk_leaves_site"));
            continue;
        }
        edit.candidate = c;
        edit.child = ChildSite{add_edit(edit.site), edit.site.bounds, where.chain};
        held.push_back(std::move(edit));
    }
    return held;
}

void TreeGenotyper::free_nodes_along_reference(nid_t bound, vector<nid_t>& out) const {
    const PathPositionHandleGraph& graph = *parts.graph;
    // Bound the work per site.
    const size_t limit = 64;
    graph.for_each_step_on_handle(graph.get_handle(bound), [&](const step_handle_t& step) {
        if (!parts.ref_paths->contains(graph.get_path_name(graph.get_path_handle_of_step(step)))) {
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
                    // A site holds the node.
                    break;
                }
                out.push_back(node);
            }
        }
        // One reference step is enough.
        return false;
    });
}

void TreeGenotyper::genotype_top_level_edits(vector<EditCandidate> candidates) {
    vector<EditCandidate> left;
    for (EditCandidate& c : candidates) {
        if (!free_edit_nodes.count(c.node)) {
            left.push_back(std::move(c));
        }
    }
    // Group the candidates by node, since each node's are genotyped together.
    vector<vector<EditCandidate>> by_node;
    for (EditCandidate& c : admitted(left)) {
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
        // Genotype each node's once. A site may be visited again as the child of a failed site
        // under RecurseOnFail.
        std::scoped_lock lock(free_edit_mutex);
        if (!free_edit_nodes.insert(candidates.front().node).second) {
            return;
        }
    }
    const PathPositionHandleGraph& graph = *parts.graph;
    const vector<EditCandidate> chosen = admitted(candidates);
    for (size_t i = 0; i < chosen.size(); ++i) {
        const EditCandidate& c = chosen[i];
        NodePlacement where = parts.child_placer->placement_of(c.node);
        if (!where.enclosing.empty()) {
            // A site holds it, so that site's descent placed it, if anything did.
            continue;
        }
        // Read the node as a reference path does.
        handle_t read_as = graph.get_handle(c.node);
        graph.for_each_step_on_handle(read_as, [&](const step_handle_t& step) {
            const string name = graph.get_path_name(graph.get_path_handle_of_step(step));
            if (parts.ref_paths->contains(name)) {
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
        edit.child = ChildSite{add_edit(edit.site), edit.site.bounds, where.chain};
        genotype_edit_site(edit, {}, "", make_pair(0, 0), -1, NestingPlacement(), true);
        parts.sites->remove_edit(edit.child.net);
    }
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
    outcome.site = edit_site_id(c.node, c.offset, c.kind, alt_label(c));
    const string& site_id = outcome.site;

    // A node no reference path visits takes its parent's path, as a site off the reference does.
    size_t node_start = 0;
    string ref_path_name = reference_through(bounds.start, parent_ref_path_name, node_start);
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
    auto called = parts.genotyper->genotype(bounds, walks, 0, ploidies, enclosing, &sequences);
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
        bounds, site_id, genotype, walks, 0, call_info, score,
        ReferencePlacement{.path_name = ref_path_name,
                           .offset = ref_offset,
                           .no_reference = no_ref_position,
                           .position_from_parent = position_from_parent},
        ploidy);
    if (!top_level) {
        place_in_tree(*staged, placement);
    }
    // Put the edit site in the chain its node is in.
    staged->chain = edit.child.chain;
    staged->enclosing = enclosing;
    staged->sequences = sequences;
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
    const SiteBounds& site, const vector<Traversal>& travs, int ref_trav_idx,
    const Ploidies& ploidies, const vector<SiteBounds>& enclosing, SiteScore*& score) const {
    auto called = parts.genotyper->genotype(site, travs, ref_trav_idx, ploidies, enclosing);
    score = called.second.get();
    return make_pair(std::move(called.first),
                     unique_ptr<SnarlCaller::CallInfo>(std::move(called.second)));
}

unique_ptr<StagedSite> TreeGenotyper::stage_render_record(
        const SiteBounds& bounds, const string& site_id, const vector<int>& trav_genotype,
        vector<Traversal> travs, int ref_trav_idx, unique_ptr<SnarlCaller::CallInfo>& call_info,
        SiteScore* score, const ReferencePlacement& reference, int ploidy) const {
    auto rec = std::make_unique<StagedSite>();
    rec->bounds = bounds;
    rec->reference = reference;
    rec->travs = std::move(travs);
    rec->ref_trav_idx = ref_trav_idx;
    rec->genotype = trav_genotype;
    rec->ploidy = ploidy;
    rec->id = site_id;
    rec->record_key = record_key_of(site_id);
    // Keep the call info, which the render reads to map the written alleles back to matrix
    // columns, index GL and compute QUAL.
    rec->set_call(std::move(call_info), score);
    return rec;
}

void TreeGenotyper::fill_tree_fields(const SiteView& view, const SiteChildren& children,
                                     StagedSite& site) const {
    site.children = children;
    site.chain = parts.child_placer->chain_of(view.net);
    site.enclosing = view.enclosing;
}

}
