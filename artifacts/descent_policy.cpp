#include <algorithm>
#include <limits>

#include "descent_policy.hpp"
#include "block_records.hpp"
#include "site_values.hpp"
#include "symbolic_allele.hpp"

namespace vg::multipass {

DescentPolicy::TraversalNodeIndex DescentPolicy::index_traversal_nodes(const HandleGraph& graph,
                                                                   const Traversal& walk) {
    TraversalNodeIndex visits;
    for (int i = 0; i < (int)walk.size(); ++i) {
        visits[graph.get_id(walk[i])].push_back(i);
    }
    return visits;
}

int DescentPolicy::crossings_of_child(const HandleGraph& graph, const TraversalNodeIndex& visits,
                                   const SiteBounds& child) {
    const nid_t start = graph.get_id(child.start);
    const nid_t end = graph.get_id(child.end);
    if (inside_one_node(child)) {
        // Each visit passes both of the child's boundaries.
        auto found = visits.find(start);
        return found == visits.end() ? 0 : (int)found->second.size();
    }
    // A crossing is a visit to one boundary node followed by a visit to the other. Only visits to
    // the two boundary nodes matter, so walk their positions merged in ascending order. When both
    // boundaries are the same node, one list plays both roles.
    static const vector<int> none;
    auto s = visits.find(start);
    auto e = visits.find(end);
    const vector<int>& sp = (s == visits.end()) ? none : s->second;
    const vector<int>& ep = (e == visits.end() || end == start) ? none : e->second;
    int crossings = 0;
    // The boundary node that would complete the current crossing, or 0 outside a crossing.
    nid_t open = 0;
    size_t i = 0;
    size_t j = 0;
    while (i < sp.size() || j < ep.size()) {
        nid_t node;
        if (j >= ep.size() || (i < sp.size() && sp[i] <= ep[j])) {
            node = start;
            ++i;
        } else {
            node = end;
            ++j;
        }
        if (open == 0 && (node == start || node == end)) {
            open = (node == start) ? end : start;
        } else if (open != 0 && node == open) {
            ++crossings;
            open = 0;
        }
    }
    return crossings;
}

int DescentPolicy::offset_of_child(const HandleGraph& graph, const Traversal& walk,
                                 const SiteBounds& child) {
    const nid_t start = graph.get_id(child.start);
    const nid_t end = graph.get_id(child.end);
    if (inside_one_node(child)) {
        // One visit is a whole crossing.
        for (int i = 0; i < (int)walk.size(); ++i) {
            if (graph.get_id(walk[i]) == start) {
                return i;
            }
        }
        return -1;
    }
    nid_t open = 0;
    int entry = -1;
    for (int i = 0; i < (int)walk.size(); ++i) {
        nid_t node = graph.get_id(walk[i]);
        if (open == 0 && (node == start || node == end)) {
            open = (node == start) ? end : start;
            entry = i;
        } else if (open != 0 && node == open) {
            // Report the entry of the first complete crossing.
            return entry;
        }
    }
    return -1;
}

DescentPolicy::ChildOffsets::ChildOffsets(const HandleGraph& graph, const Traversal& indexed)
    : walk(&indexed) {
    bases_before.assign(indexed.size() + 1, 0);
    for (int i = 0; i < (int)indexed.size(); ++i) {
        visits_of[graph.get_id(indexed[i])].push_back(i);
        bases_before[i + 1] = bases_before[i] + (int64_t)graph.get_length(indexed[i]);
    }
}

int64_t DescentPolicy::ChildOffsets::base_offset(const HandleGraph& graph,
                                               const SiteBounds& child) const {
    // Follow the rule of `offset_of_child`. The entry is the first visit to either boundary node,
    // and it counts only if the other boundary node is visited after it.
    const nid_t start = graph.get_id(child.start);
    const nid_t end = graph.get_id(child.end);
    auto start_visits = visits_of.find(start);
    auto end_visits = visits_of.find(end);
    if (inside_one_node(child)) {
        // The first visit is a whole crossing, and the child starts at its base inside the node,
        // read the way this visit reads the node.
        if (start_visits == visits_of.end()) {
            return -1;
        }
        const int visit = start_visits->second.front();
        SiteBounds as_visited = child;
        as_visited.start = (*walk)[visit];
        return bases_before[visit] + (int64_t)bases_before_site(graph, as_visited);
    }
    const int none = numeric_limits<int>::max();
    const int first_start = start_visits == visits_of.end() ? none : start_visits->second.front();
    const int first_end = end_visits == visits_of.end() ? none : end_visits->second.front();
    const int entry = min(first_start, first_end);
    if (entry == none) {
        return -1;
    }
    if (const auto& closing = first_start <= first_end ? end_visits : start_visits;
        closing == visits_of.end()
        || std::ranges::upper_bound(closing->second, entry) == closing->second.end()) {
        return -1;
    }
    if (child.inside_node) {
        // A child inside nodes starts at its base inside the first.
        SiteBounds as_visited = child;
        as_visited.start = (*walk)[entry];
        return bases_before[entry] + (int64_t)bases_before_site(graph, as_visited);
    }
    return bases_before[entry];
}

size_t DescentPolicy::offset_along_genotype(
    const HandleGraph& graph, const vector<Traversal>& travs, const vector<int>& genotype,
    const SiteBounds& child, unordered_map<const Traversal*, ChildOffsets>& offsets) {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const Traversal* walk = &travs[allele];
        auto found = offsets.find(walk);
        if (found == offsets.end()) {
            found = offsets.try_emplace(walk, graph, *walk).first;
        }
        const int64_t within = found->second.base_offset(graph, child);
        if (within >= 0) {
            return (size_t)within;
        }
    }
    return 0;
}

uint64_t DescentPolicy::child_crossing_mask(const HandleGraph& graph,
                                         const vector<TraversalNodeIndex>& visits,
                                         const SiteBounds& child, bool* known) {
    if (known != nullptr) {
        *known = true;
    }
    if (visits.size() > 64) {
        if (known != nullptr) {
            *known = false;
        }
        return 0;
    }
    uint64_t mask = 0;
    for (size_t i = 0; i < visits.size(); ++i) {
        if (crossings_of_child(graph, visits[i], child) > 0) {
            mask |= (uint64_t)1 << i;
        }
    }
    return mask;
}

int DescentPolicy::child_ploidy(const vector<TraversalNodeIndex>& visits,
                              const vector<int>& genotype, const SiteBounds& child,
                              int cap) const {
    int copies = 0;
    bool capped = false;

    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)visits.size()) {
            // A star or missing allele carries no copy.
            continue;
        }
        int crossings = crossings_of_child(*graph, visits[allele], child);
        if (crossings > 1) {
            capped = true;
            crossings = 1;
        }
        copies += crossings;
    }
    if (capped) {
        ++counters->child_multi_crossing;
    }
    return min(copies, cap);
}

void DescentPolicy::configure(const HandleGraph* new_graph,
                            const SnarlDecomposition* new_decomposition,
                            const BlockRecordWriter* new_blocks, DescentCounters* new_counters) {
    graph = new_graph;
    decomposition = new_decomposition;
    blocks = new_blocks;
    counters = new_counters;
}

DescentPolicy::Nested DescentPolicy::nested_in(const net_handle_t& site,
                                           const SiteBounds& as_called) const {
    Nested nested;
    nested.children = site_children(*decomposition, *graph, site, as_called, &nested.sites);
    return nested;
}

SiteView DescentPolicy::view_of(const SiteView& parent, const ChildSite& child) {
    SiteView view{child.net, child.bounds, {}};
    // The sites enclosing the child are its parent, then those enclosing the parent.
    view.enclosing.reserve(parent.enclosing.size() + 1);
    view.enclosing.push_back(parent.bounds);
    view.enclosing.insert(view.enclosing.end(), parent.enclosing.begin(), parent.enclosing.end());
    return view;
}

ChildChain DescentPolicy::chain_of(const net_handle_t& site) const {
    return chain_of_site(*decomposition, *graph, site);
}

NodePlacement DescentPolicy::placement_of(nid_t id) const {
    return placement_of_node(*decomposition, *graph, id);
}

bool DescentPolicy::reference_crosses(const Parent& parent,
                                    const vector<TraversalNodeIndex>& visits,
                                    const SiteBounds& child) const {
    if (parent.ref_trav_idx < 0 || parent.ref_trav_idx >= (int)parent.travs.size()) {
        // With no reference walk, there is nothing to test.
        return true;
    }
    const vector<int> ref_only(1, parent.ref_trav_idx);
    return child_ploidy(visits, ref_only, child, 1) != 0;
}

void DescentPolicy::place(const Parent& parent, bool off_reference, bool keep_uncrossed,
                        const function<void(const Placed& child)>& visit) const {
    if (!parent.nested.children.known) {
        return;
    }
    const vector<Traversal>& travs = parent.travs;
    const NestingPlacement& placement = parent.placement;
    // Build once for the parent the parts of the test for an inline-reported chain that do not
    // depend on the child.
    const BlockRecordWriter::ChainInlineContext inline_ctx = blocks->chain_inline_context(
        *graph, parent.nested.children, travs, parent.genotype, parent.ref_trav_idx);
    vector<TraversalNodeIndex> trav_visits;
    trav_visits.reserve(travs.size());
    for (const Traversal& t : travs) {
        trav_visits.push_back(index_traversal_nodes(*graph, t));
    }
    // Index the base offsets of each walk the first time a child needs that walk.
    unordered_map<const Traversal*, ChildOffsets> walk_offsets;
    for (const ChildSite& child : parent.nested.sites) {
        const SiteBounds& child_bounds = child.bounds;
        // A chain the reference walk does not cross has no REF or POS for its records, so leave
        // it out unless off-reference descent is on. Everything under a chain off the reference
        // is off it too.
        bool child_off_reference = placement.no_reference;
        if (!reference_crosses(parent, trav_visits, child_bounds)) {
            if (!off_reference) {
                ++counters->skipped_no_ref;
                continue;
            }
            child_off_reference = true;
            ++counters->off_reference;
        }

        // Under block emission, a difference block's ALT already spells a chain that every called
        // strand crosses only inside that block. That holds back the chain's line, but the chain
        // is still genotyped, filed and phased.
        const bool child_reported_inline =
            placement.reported_inline
            || blocks->chain_reported_inline(*graph, inline_ctx, child.chain);

        const int copies = child_ploidy(trav_visits, parent.genotype, child_bounds, parent.ploidy);
        bool retain_only = placement.retain_only;
        if (copies <= 0) {
            ++counters->skipped_no_copy;
            if (!keep_uncrossed) {
                // Without the linkage model nothing moves the parent after the direct pass, so
                // the sample has no copy of this chain.
                continue;
            }
            // Genotype it anyway, while the reads of this window are in memory. Nothing about it
            // is written unless the linkage pass says so.
            retain_only = true;
        }

        Placed out;
        out.site = view_of(parent.site, child);
        out.chain = child.chain;
        // A chain no called allele crosses takes the parent's ploidy, the most copies it can have.
        out.ploidy = copies >= 1 ? copies : parent.ploidy;
        NestingPlacement& next = out.placement;
        next.one_copy = (copies == 1);
        next.parent_record_key = parent.record_key;
        next.retain_only = retain_only;
        next.no_reference = child_off_reference;
        // Only an off-reference chain uses the offset, but compute it for every chain, so that
        // offsets add up down the tree.
        next.parent_offset =
            placement.parent_offset
            + offset_along_genotype(*graph, travs, parent.genotype, child_bounds, walk_offsets);
        next.reported_inline = child_reported_inline;
        next.chain_key = chain_key_of(*graph, out.chain);
        bool crossing_known = true;
        next.parent_crossing =
            child_crossing_mask(*graph, trav_visits, child_bounds, &crossing_known);
        next.crossing_known = crossing_known;
        next.level = placement.level + 1;
        visit(out);
    }
}

}
