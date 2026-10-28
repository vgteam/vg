#include <algorithm>

#include "edit_site_placer.hpp"
#include "../src/utility.hpp"

namespace vg::multipass {

/// Make the outcome of a candidate for the dump, saying where it was placed, or why it was not.
EditOutcome unplaced_outcome(const EditCandidate& c, const string& why) {
    EditOutcome outcome;
    outcome.node = c.node;
    outcome.offset = c.offset;
    outcome.alt = c.alt;
    outcome.kind = c.kind;
    outcome.sequence = c.sequence;
    outcome.placed = why;
    return outcome;
}

vector<EditCandidate> EditSitePlacer::admitted(const vector<EditCandidate>& candidates) const {
    // The candidates of each base, substitutions and indels separately, are consecutive. Keep one
    // of each.
    vector<EditCandidate> chosen;
    for (const EditCandidate& c : candidates) {
        if (!bar.admits(c)) {
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

bool EditSitePlacer::walk_in_one_site(const EditCandidate& c, const NodePlacement& where) const {
    const net_handle_t holder = parts.sites->canonical(where.site);
    for (nid_t id : c.walk) {
        if (parts.sites->canonical(parts.descent_policy->placement_of(id).site) != holder) {
            return false;
        }
    }
    return true;
}

bool EditSitePlacer::reference_reads_whole_walk(const EditCandidate& c) const {
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

net_handle_t EditSitePlacer::add_edit(const EditSite& site) const {
    const Traversal& walk = site.alleles.walks.front();
    // The site runs from its first base to the end of its walk.
    return parts.sites->add_edit(walk.front(), site.start_offset, walk.back(),
                                 parts.graph->get_length(walk.back()));
}

bool EditSitePlacer::edit_site_of(const EditCandidate& c, handle_t read_as,
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

vector<HeldEdit> EditSitePlacer::held_by(const net_handle_t& site,
                                         const vector<Traversal>& walks, int ref_trav_idx) const {
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
    const nid_t start = graph.get_id(parts.sites->get_site_start(site));
    const nid_t end = graph.get_id(parts.sites->get_site_end(site));
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
    const net_handle_t holder = parts.sites->canonical(site);
    for (const EditCandidate& c : admitted(found)) {
        // Only the innermost site holding the node holds the edit.
        NodePlacement where = parts.descent_policy->placement_of(c.node);
        if (parts.sites->canonical(where.site) != holder) {
            continue;
        }
        EditSite site;
        if (!edit_site_of(c, read_as_in_walks(graph, walks, ref_trav_idx, c.node), where, site)) {
            parts.edit_outcomes->add(unplaced_outcome(c, "walk_leaves_site"));
            continue;
        }
        HeldEdit edit;
        edit.candidate = c;
        edit.net = add_edit(site);
        edit.alleles = std::move(site.alleles);
        held.push_back(std::move(edit));
    }
    return held;
}

void EditSitePlacer::free_nodes_along_reference(nid_t bound, vector<nid_t>& out) const {
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
                if (!parts.sites->is_root(parts.descent_policy->placement_of(node).site)) {
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

vector<vector<EditCandidate>> EditSitePlacer::group_unclaimed(vector<EditCandidate> candidates) {
    vector<EditCandidate> left;
    for (EditCandidate& c : candidates) {
        if (!free_edit_nodes.count(c.node)) {
            left.push_back(std::move(c));
        }
    }
    // Group the candidates by node, since each node's are placed together.
    vector<vector<EditCandidate>> by_node;
    for (EditCandidate& c : admitted(left)) {
        if (by_node.empty() || by_node.back().front().node != c.node) {
            by_node.emplace_back();
        }
        by_node.back().push_back(std::move(c));
    }
    after_pass_nodes = by_node.size();
    return by_node;
}

void EditSitePlacer::place_free_edits(vector<EditCandidate> candidates,
                                      const function<void(const HeldEdit& edit)>& genotype) {
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
        NodePlacement where = parts.descent_policy->placement_of(c.node);
        if (!parts.sites->is_root(where.site)) {
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
        EditSite site;
        if (!edit_site_of(c, read_as, where, site)) {
            parts.edit_outcomes->add(unplaced_outcome(c, "walk_leaves_site"));
            continue;
        }
        HeldEdit edit;
        edit.candidate = c;
        edit.net = add_edit(site);
        edit.alleles = std::move(site.alleles);
        genotype(edit);
    }
}

}
