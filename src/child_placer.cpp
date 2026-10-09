#include <algorithm>
#include <limits>

#include "child_placer.hpp"
#include "flow_caller.hpp"

namespace vg {

ChildPlacer::TraversalNodeIndex ChildPlacer::index_traversal_nodes(const SnarlTraversal& trav) {
    TraversalNodeIndex visits;
    for (int i = 0; i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        visits[trav.visit(i).node_id()].push_back(i);
    }
    return visits;
}

int ChildPlacer::crossings_of_child(const TraversalNodeIndex& visits, const Snarl& child) {
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    // Count crossings: an entry at one boundary followed by the other. Order matters: testing for
    // the two boundaries separately would count a traversal that touches both on unrelated
    // excursions, as find_child_traversal_set does. Only visits to the two boundary nodes can
    // change the state, so we walk their positions merged in ascending order. When both
    // boundaries are the same node, one list plays both roles.
    static const vector<int> none;
    auto s = visits.find(start);
    auto e = visits.find(end);
    const vector<int>& sp = (s == visits.end()) ? none : s->second;
    const vector<int>& ep = (e == visits.end() || end == start) ? none : e->second;
    int crossings = 0;
    nid_t open = 0;
    size_t i = 0, j = 0;
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

int ChildPlacer::offset_of_child(const SnarlTraversal& trav, const Snarl& child) {
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    nid_t open = 0;
    int entry = -1;
    for (int i = 0; i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        nid_t node = trav.visit(i).node_id();
        if (open == 0 && (node == start || node == end)) {
            open = (node == start) ? end : start;
            entry = i;
        } else if (open != 0 && node == open) {
            return entry;   // first complete crossing, entry side
        }
    }
    return -1;
}

ChildPlacer::ChildOffsets::ChildOffsets(const HandleGraph& graph, const SnarlTraversal& trav) {
    bases_before.assign((size_t)trav.visit_size() + 1, 0);
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        int64_t bases = 0;
        if (!visit.has_snarl()) {
            visits_of[visit.node_id()].push_back(i);
            bases = (int64_t)graph.get_length(graph.get_handle(visit.node_id()));
        }
        bases_before[i + 1] = bases_before[i] + bases;
    }
}

int64_t ChildPlacer::ChildOffsets::base_offset(const Snarl& child) const {
    // `offset_of_child`'s rule: the entry is the first visit to either boundary node, and it counts
    // only if the other boundary node is visited after it.
    const nid_t start = child.start().node_id();
    const nid_t end = child.end().node_id();
    auto start_visits = visits_of.find(start);
    auto end_visits = visits_of.find(end);
    const int none = numeric_limits<int>::max();
    const int first_start = start_visits == visits_of.end() ? none : start_visits->second.front();
    const int first_end = end_visits == visits_of.end() ? none : end_visits->second.front();
    const int entry = min(first_start, first_end);
    if (entry == none) {
        return -1;
    }
    const auto& closing = first_start <= first_end ? end_visits : start_visits;
    if (closing == visits_of.end()
        || std::upper_bound(closing->second.begin(), closing->second.end(), entry)
               == closing->second.end()) {
        return -1;
    }
    return bases_before[entry];
}

size_t ChildPlacer::offset_along_genotype(
    const HandleGraph& graph, const vector<SnarlTraversal>& travs, const vector<int>& genotype,
    const Snarl& child, unordered_map<const SnarlTraversal*, ChildOffsets>& offsets) {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const SnarlTraversal* trav = &travs[allele];
        auto found = offsets.find(trav);
        if (found == offsets.end()) {
            found = offsets.emplace(trav, ChildOffsets(graph, *trav)).first;
        }
        const int64_t within = found->second.base_offset(child);
        if (within >= 0) {
            return (size_t)within;
        }
    }
    return 0;
}

uint64_t ChildPlacer::child_crossing_mask(const vector<TraversalNodeIndex>& visits,
                                         const Snarl& child, bool* known) {
    if (known != nullptr) {
        *known = true;
    }
    // One bit per candidate traversal, not per VCF allele, since the linkage model chooses a site's
    // genotype as a pair of traversals.
    if (visits.size() > 64) {
        // Unknown rather than none: the mask cannot index this site's candidates.
        if (known != nullptr) {
            *known = false;
        }
        return 0;
    }
    uint64_t mask = 0;
    for (size_t i = 0; i < visits.size(); ++i) {
        if (crossings_of_child(visits[i], child) > 0) {
            mask |= (uint64_t)1 << i;
        }
    }
    return mask;
}

int64_t FlowCaller::base_offset_of_child(const SnarlTraversal& trav, const Snarl& child) const {
    const int entry = ChildPlacer::offset_of_child(trav, child);
    if (entry < 0) {
        return -1;
    }
    int64_t bases = 0;
    for (int i = 0; i < entry && i < trav.visit_size(); ++i) {
        if (trav.visit(i).has_snarl()) {
            continue;
        }
        bases += (int64_t)graph.get_length(graph.get_handle(trav.visit(i).node_id()));
    }
    return bases;
}

size_t FlowCaller::offset_along_genotype(const vector<SnarlTraversal>& travs,
                                         const vector<int>& genotype, const Snarl& child) const {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const int64_t within = base_offset_of_child(travs[allele], child);
        if (within >= 0) {
            return (size_t)within;
        }
    }
    return 0;
}

int FlowCaller::child_ploidy(const vector<ChildPlacer::TraversalNodeIndex>& visits,
                             const vector<int>& genotype,
                             const Snarl& child, int cap) const {
    int copies = 0;
    bool capped = false;

    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)visits.size()) {
            continue;   // star or missing: that haplotype contributes no copy here
        }
        int crossings = ChildPlacer::crossings_of_child(visits[allele], child);
        if (crossings > 1) {
            capped = true;
            crossings = 1;   // a cycle or tandem duplication; see the header comment
        }
        copies += crossings;
    }
    if (capped) {
        // Counted and reported once per run.
        ++descent_counters.child_multi_crossing;
    }
    return min(copies, cap);
}

}
