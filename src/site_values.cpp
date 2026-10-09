#include "site_values.hpp"

#include <algorithm>
#include <map>

namespace vg {

const ChildChain* SiteChildren::entered_by(nid_t id, bool backward) const {
    const pair<nid_t, bool> key(id, backward);
    auto found = std::lower_bound(entries.begin(), entries.end(), key,
                                  [](const pair<pair<nid_t, bool>, size_t>& entry,
                                     const pair<nid_t, bool>& k) { return entry.first < k; });
    if (found == entries.end() || found->first != key) {
        return nullptr;
    }
    return &chains[found->second];
}

SiteBounds bounds_of(const HandleGraph& graph, const Snarl& snarl) {
    return SiteBounds{graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
                      graph.get_handle(snarl.end().node_id(), snarl.end().backward())};
}

Snarl snarl_of(const HandleGraph& graph, const SiteBounds& site) {
    Snarl snarl;
    snarl.mutable_start()->set_node_id(graph.get_id(site.start));
    snarl.mutable_start()->set_backward(graph.get_is_reverse(site.start));
    snarl.mutable_end()->set_node_id(graph.get_id(site.end));
    snarl.mutable_end()->set_backward(graph.get_is_reverse(site.end));
    return snarl;
}

Traversal walk_of(const HandleGraph& graph, const SnarlTraversal& trav) {
    Traversal walk;
    walk.reserve(trav.visit_size());
    for (const Visit& visit : trav.visit()) {
        walk.push_back(graph.get_handle(visit.node_id(), visit.backward()));
    }
    return walk;
}

vector<Traversal> walks_of(const HandleGraph& graph, const vector<SnarlTraversal>& travs) {
    vector<Traversal> walks;
    walks.reserve(travs.size());
    for (const SnarlTraversal& trav : travs) {
        walks.push_back(walk_of(graph, trav));
    }
    return walks;
}

SnarlTraversal snarl_traversal_of(const HandleGraph& graph, const Traversal& walk) {
    SnarlTraversal trav;
    for (const handle_t& handle : walk) {
        Visit* visit = trav.add_visit();
        visit->set_node_id(graph.get_id(handle));
        visit->set_backward(graph.get_is_reverse(handle));
    }
    return trav;
}

vector<SnarlTraversal> snarl_traversals_of(const HandleGraph& graph,
                                           const vector<Traversal>& walks) {
    vector<SnarlTraversal> travs;
    travs.reserve(walks.size());
    for (const Traversal& walk : walks) {
        travs.push_back(snarl_traversal_of(graph, walk));
    }
    return travs;
}

bool same_walk(const HandleGraph& graph, const SnarlTraversal& trav, const Traversal& walk) {
    if ((size_t)trav.visit_size() != walk.size()) {
        return false;
    }
    for (size_t i = 0; i < walk.size(); ++i) {
        const Visit& visit = trav.visit(i);
        if (visit.has_snarl() || visit.node_id() != graph.get_id(walk[i])
            || visit.backward() != graph.get_is_reverse(walk[i])) {
            return false;
        }
    }
    return true;
}

/// The manager's snarl that `site` is, or null when there is none. `out_reversed` is set when it
/// is that snarl only turned round.
static const Snarl* managed_site(const Snarl& site, const SnarlManager& manager,
                                 bool* out_reversed) {
    *out_reversed = false;
    const Snarl* site_ptr = manager.into_which_snarl(site.start().node_id(),
                                                     site.start().backward());
    if (site_ptr == nullptr) {
        return nullptr;
    }
    // The start leads into this snarl, but it is the same snarl only if the other bound agrees,
    // read either way round. The forward pairing compares node IDs; the reversed one compares
    // whole visits, since the turned-round snarl has its bounds reversed and swapped.
    auto same_visit = [](const Visit& a, const Visit& b) {
        return a.node_id() == b.node_id() && a.backward() == b.backward();
    };
    const bool forward = site_ptr->start().node_id() == site.start().node_id() &&
                         site_ptr->end().node_id() == site.end().node_id();
    const bool reversed = same_visit(site_ptr->start(), reverse(site.end())) &&
                          same_visit(site_ptr->end(), reverse(site.start()));
    if (!forward && !reversed) {
        return nullptr;
    }
    *out_reversed = reversed && !forward;
    return site_ptr;
}

ChildChain chain_of_site(const SnarlManager& manager, const HandleGraph& graph,
                         const Snarl* site) {
    const Chain* chain = manager.chain_of(site);
    const bool in_chain = chain != nullptr && !chain->empty();
    const Visit first = in_chain ? get_start_of(*chain) : site->start();
    const Visit last = in_chain ? get_end_of(*chain) : site->end();
    return ChildChain{graph.get_handle(first.node_id(), first.backward()),
                      graph.get_handle(last.node_id(), last.backward())};
}

SiteChildren site_children(const SnarlManager& manager, const HandleGraph& graph,
                           const Snarl& site) {
    SiteChildren out;
    const Snarl* site_ptr = managed_site(site, manager, &out.reversed);
    if (site_ptr == nullptr) {
        return out;
    }
    out.known = true;
    // Each chain once, however many of its sites are children.
    map<pair<nid_t, nid_t>, size_t> chain_index;
    for (const Snarl* child : manager.children_of(site_ptr)) {
        const ChildChain chain = chain_of_site(manager, graph, child);
        auto inserted = chain_index.emplace(make_pair(graph.get_id(chain.start),
                                                      graph.get_id(chain.end)),
                                            out.chains.size());
        if (inserted.second) {
            out.chains.push_back(chain);
        }
        // The two ways in, each kept only where the manager's own lookup gives this child.
        for (const pair<nid_t, bool>& way_in :
             {make_pair(child->start().node_id(), child->start().backward()),
              make_pair(child->end().node_id(), !child->end().backward())}) {
            if (manager.into_which_snarl(way_in.first, way_in.second) == child) {
                out.entries.emplace_back(way_in, inserted.first->second);
            }
        }
    }
    std::sort(out.entries.begin(), out.entries.end());
    return out;
}

vector<SiteBounds> enclosing_sites(const SnarlManager& manager, const HandleGraph& graph,
                                   const Snarl& site) {
    vector<SiteBounds> out;
    // Up from the snarl a walk enters by the site's start, one parent at a time.
    const Snarl* managed = manager.into_which_snarl(site.start().node_id(),
                                                    site.start().backward());
    const Snarl* current = managed == nullptr ? nullptr : manager.parent_of(managed);
    while (current != nullptr) {
        out.push_back(bounds_of(graph, *current));
        managed = manager.into_which_snarl(current->start().node_id(),
                                           current->start().backward());
        current = managed == nullptr ? nullptr : manager.parent_of(managed);
    }
    return out;
}

}
