#include "site_values.hpp"

#include <algorithm>

namespace vg {
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

namespace multipass {

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

SiteBounds oriented_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                           const net_handle_t& net) {
    // The bounds ignore the traversal: the start bound read in, and the end bound read out.
    return SiteBounds{decomposition.get_handle(decomposition.get_bound(net, false, true), &graph),
                      decomposition.get_handle(decomposition.get_bound(net, true, false), &graph)};
}

SiteChildren site_children(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                           const net_handle_t& site, const SiteBounds& as_called,
                           vector<ChildSite>* child_sites) {
    SiteChildren out;
    // The caller's site is this one as it is, or turned round. The forward pairing compares node
    // IDs; the reversed one compares whole handles, since the turned-round site has its bounds
    // reversed and swapped.
    const SiteBounds own = oriented_bounds(decomposition, graph, site);
    const bool forward = graph.get_id(own.start) == graph.get_id(as_called.start)
                         && graph.get_id(own.end) == graph.get_id(as_called.end);
    const bool reversed = own.start == graph.flip(as_called.end)
                          && own.end == graph.flip(as_called.start);
    if (!forward && !reversed) {
        return out;
    }
    out.known = true;
    out.reversed = reversed && !forward;
    decomposition.for_each_child(site, [&](const net_handle_t& chain) {
        if (!decomposition.is_chain(chain)) {
            return true;
        }
        // The chain's nodes in order, read the way the chain is read. A site lies between each
        // two consecutive nodes, and between the last and the first of a cyclic chain; it is
        // entered by the first node read forward or by the second read backward.
        const ChildChain bounds = oriented_bounds(decomposition, graph, chain);
        vector<handle_t> nodes;
        decomposition.for_each_child(chain, [&](const net_handle_t& child) {
            if (decomposition.is_node(child)) {
                nodes.push_back(decomposition.get_handle(child, &graph));
            } else if (child_sites != nullptr && decomposition.is_snarl(child)) {
                child_sites->push_back(
                    ChildSite{child, oriented_bounds(decomposition, graph, child), bounds});
            }
            return true;
        });
        const bool cyclic = nodes.size() > 1
                            && graph.get_id(bounds.start) == graph.get_id(bounds.end);
        const size_t pairs = nodes.empty() ? 0 : (cyclic ? nodes.size() : nodes.size() - 1);
        if (pairs == 0) {
            return true;
        }
        const size_t index = out.chains.size();
        out.chains.push_back(bounds);
        for (size_t k = 0; k < pairs; ++k) {
            const handle_t& first = nodes[k];
            const handle_t& second = nodes[(k + 1) % nodes.size()];
            out.entries.emplace_back(make_pair(graph.get_id(first), graph.get_is_reverse(first)),
                                     index);
            out.entries.emplace_back(make_pair(graph.get_id(second), !graph.get_is_reverse(second)),
                                     index);
        }
        return true;
    });
    std::sort(out.entries.begin(), out.entries.end());
    return out;
}

vector<SiteBounds> enclosing_sites(const SnarlDecomposition& decomposition,
                                   const HandleGraph& graph, const net_handle_t& site) {
    vector<SiteBounds> out;
    for (net_handle_t parent = parent_site(decomposition, site); !decomposition.is_root(parent);
         parent = parent_site(decomposition, parent)) {
        out.push_back(oriented_bounds(decomposition, graph, parent));
    }
    return out;
}

ChildChain chain_of_site(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& site) {
    return oriented_bounds(decomposition, graph, decomposition.get_parent(site));
}

}
}
