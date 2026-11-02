/** \file
 * Implements ChildCrossingFinder.
 */

#include "child_crossing_finder.hpp"

namespace vg::multipass {

ChildCrossingFinder::ChildCrossingFinder(const EditSnarlDecompositionOverlay& overlay,
                                         const HandleGraph& walk_graph, const Traversal& walked)
    : sites(&overlay), graph(&walk_graph), walk(&walked) {
    bases_before.assign(walked.size() + 1, 0);
    for (int i = 0; i < (int)walked.size(); ++i) {
        visits_of[walk_graph.get_id(walked[i])].push_back(i);
        bases_before[i + 1] = bases_before[i] + (int64_t)walk_graph.get_length(walked[i]);
    }
}

int64_t ChildCrossingFinder::bases_before_child(const net_handle_t& child, int step) const {
    // An edit read forward has its first base at its start.
    const net_handle_t forward = sites->canonical(child);
    const handle_t start = sites->get_site_start(forward);
    const auto before = (int64_t)sites->get_site_start_offset(forward);
    const handle_t visit = (*walk)[step];
    // Read the other way, the bases before the first base are the ones after it.
    return graph->get_is_reverse(visit) == graph->get_is_reverse(start)
               ? before
               : (int64_t)graph->get_length(visit) - 1 - before;
}

vector<ChildCrossing> ChildCrossingFinder::crossings(const net_handle_t& child) const {
    const net_handle_t forward = sites->canonical(child);
    const handle_t start_handle = sites->get_site_start(forward);
    const handle_t end_handle = sites->get_site_end(forward);
    const nid_t start = graph->get_id(start_handle);
    const nid_t end = graph->get_id(end_handle);
    const bool inside_nodes = sites->is_edit(child);
    auto crossing_at = [&](int step) {
        return ChildCrossing{step, bases_before[step]
                                       + (inside_nodes ? bases_before_child(child, step) : 0)};
    };

    vector<ChildCrossing> out;
    static const vector<int> none;
    auto s = visits_of.find(start);
    if (inside_nodes && start_handle == end_handle) {
        // Each visit to the node passes both of the child's boundaries.
        for (int step : s == visits_of.end() ? none : s->second) {
            out.push_back(crossing_at(step));
        }
        return out;
    }
    // A crossing is a visit to one boundary node followed by a visit to the other. Only visits to
    // the two boundary nodes matter, so walk their steps merged in ascending order. When both
    // boundaries are the same node, one list plays both roles.
    auto e = visits_of.find(end);
    const vector<int>& sp = (s == visits_of.end()) ? none : s->second;
    const vector<int>& ep = (e == visits_of.end() || end == start) ? none : e->second;
    // The boundary node that would complete the current crossing, or 0 outside a crossing.
    nid_t open = 0;
    int entry = -1;
    size_t i = 0;
    size_t j = 0;
    while (i < sp.size() || j < ep.size()) {
        nid_t node;
        int step;
        if (j >= ep.size() || (i < sp.size() && sp[i] <= ep[j])) {
            node = start;
            step = sp[i++];
        } else {
            node = end;
            step = ep[j++];
        }
        if (open == 0) {
            open = (node == start) ? end : start;
            entry = step;
        } else if (node == open) {
            out.push_back(crossing_at(entry));
            open = 0;
        }
    }
    return out;
}

uint64_t crossing_mask(const vector<ChildCrossingFinder>& walks, const net_handle_t& child,
                       bool* known) {
    if (known != nullptr) {
        *known = walks.size() <= 64;
    }
    if (walks.size() > 64) {
        return 0;
    }
    uint64_t mask = 0;
    for (size_t i = 0; i < walks.size(); ++i) {
        if (!walks[i].crossings(child).empty()) {
            mask |= (uint64_t)1 << i;
        }
    }
    return mask;
}

size_t entry_bases_along_genotype(const EditSnarlDecompositionOverlay& sites,
                                  const HandleGraph& graph, const vector<Traversal>& travs,
                                  const vector<int>& genotype, const net_handle_t& child,
                                  unordered_map<const Traversal*, ChildCrossingFinder>& finders) {
    for (int allele : genotype) {
        if (allele < 0 || allele >= (int)travs.size()) {
            continue;
        }
        const Traversal* walk = &travs[allele];
        auto found = finders.find(walk);
        if (found == finders.end()) {
            found = finders.try_emplace(walk, sites, graph, *walk).first;
        }
        const vector<ChildCrossing> crossed = found->second.crossings(child);
        if (!crossed.empty()) {
            return (size_t)crossed.front().entry_bases;
        }
    }
    return 0;
}

}
