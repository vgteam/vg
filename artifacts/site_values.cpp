#include "site_values.hpp"
#include "../src/utility.hpp"

#include <algorithm>

namespace vg {
EditSite edit_site(const HandleGraph& graph, const Traversal& steps, uint32_t offset, size_t at,
                   size_t length, const vector<string>& alts) {
    EditSite site;
    site.start_offset = graph.get_is_reverse(steps.front())
                            ? graph.get_length(steps.front()) - 1 - offset
                            : offset;
    site.alleles.walks.assign(alts.size() + 1, steps);
    site.alleles.ref_allele = 0;
    site.alleles.edits.resize(alts.size() + 1);
    const size_t host = length > 0 ? at : at - 1;
    for (size_t a = 0; a < alts.size(); ++a) {
        size_t start = 0;
        for (size_t k = 0; k < steps.size(); ++k) {
            const size_t end = start + graph.get_length(steps[k]);
            // Cut the replaced bases this step holds, and put the ALT in the step holding the
            // base at `at`, or for an insertion the base before it, where the cut starts.
            const size_t from = std::max(start, at);
            const size_t to = std::min(end, at + length);
            const size_t removed = from < to ? to - from : 0;
            if (host >= start && host < end) {
                site.alleles.edits[a + 1].push_back(AlleleEdit{k, at - start, removed, alts[a]});
            } else if (removed > 0) {
                site.alleles.edits[a + 1].push_back(AlleleEdit{k, from - start, removed, ""});
            }
            start = end;
        }
    }
    return site;
}

string edit_site_id(nid_t node, uint32_t offset, char kind, const string& alt) {
    return "op" + std::to_string(node) + "." + std::to_string(offset) + "." + kind + alt;
}

namespace multipass {

const ChildChain* SiteChildren::entered_by(nid_t id, bool backward) const {
    const pair<nid_t, bool> key(id, backward);
    auto found = std::ranges::lower_bound(entries, key, {},
                                          &pair<pair<nid_t, bool>, size_t>::first);
    if (found == entries.end() || found->first != key) {
        return nullptr;
    }
    return &chains[found->second];
}

pair<handle_t, handle_t> oriented_ends(const SnarlDecomposition& decomposition,
                                       const HandleGraph& graph, const net_handle_t& net) {
    // get_bound() ignores the traversal, so this is the start bound read in and the end bound
    // read out.
    return {decomposition.get_handle(decomposition.get_bound(net, false, true), &graph),
            decomposition.get_handle(decomposition.get_bound(net, true, false), &graph)};
}

/**
 * Add a child chain of a site to the site's `SiteChildren`, if a site lies in the chain.
 *
 * Every two consecutive nodes of the chain, read the way the chain is read, have a site between
 * them, and so do the last and first nodes of a looping chain. A walk enters that site by reading
 * the first of the two nodes forward, or the second backward. Each of those is an entry.
 *
 * If `child_sites` is not null, also add the chain's child sites to it, in chain order.
 */
static void add_child_chain(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                            const net_handle_t& chain, SiteChildren& out,
                            vector<net_handle_t>* child_sites) {
    const ChildChain bounds = chain_ends(decomposition, graph, chain);
    if (child_sites != nullptr) {
        for (const net_handle_t& child : vg::child_sites(decomposition, chain)) {
            child_sites->push_back(child);
        }
    }
    vector<handle_t> nodes;
    decomposition.for_each_child(chain, [&](const net_handle_t& child) {
        if (decomposition.is_node(child)) {
            nodes.push_back(decomposition.get_handle(child, &graph));
        }
        return true;
    });
    const bool looping = nodes.size() > 1
                         && graph.get_id(bounds.start) == graph.get_id(bounds.end);
    size_t pairs = 0;
    if (looping) {
        pairs = nodes.size();
    } else if (!nodes.empty()) {
        pairs = nodes.size() - 1;
    }
    if (pairs == 0) {
        return;
    }
    const size_t index = out.chains.size();
    out.chains.push_back(bounds);
    for (size_t k = 0; k < pairs; ++k) {
        const handle_t& first = nodes[k];
        const handle_t& second = nodes[(k + 1) % nodes.size()];
        out.entries.push_back({{graph.get_id(first), graph.get_is_reverse(first)}, index});
        out.entries.push_back({{graph.get_id(second), !graph.get_is_reverse(second)}, index});
    }
}

SiteChildren site_children(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                           const net_handle_t& site, const handle_t& called_start,
                           const handle_t& called_end, vector<net_handle_t>* child_sites) {
    SiteChildren out;
    // The caller has this site as the decomposition orients it, or turned round. Compare node IDs
    // for the first, and whole handles for the second, since turning round flips the bounds as
    // well as swapping them.
    const auto [own_start, own_end] = oriented_ends(decomposition, graph, site);
    const bool forward = graph.get_id(own_start) == graph.get_id(called_start)
                         && graph.get_id(own_end) == graph.get_id(called_end);
    const bool reversed = own_start == graph.flip(called_end)
                          && own_end == graph.flip(called_start);
    if (!forward && !reversed) {
        return out;
    }
    out.known = true;
    out.reversed = reversed && !forward;
    decomposition.for_each_child(site, [&](const net_handle_t& chain) {
        if (decomposition.is_chain(chain)) {
            add_child_chain(decomposition, graph, chain, out, child_sites);
        }
        return true;
    });
    std::ranges::sort(out.entries);
    return out;
}

ChildChain chain_ends(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                      const net_handle_t& chain) {
    // get_bound() ignores the traversal, so this is the start bound read in and the end bound
    // read out.
    return ChildChain{
        decomposition.get_handle(decomposition.get_bound(chain, false, true), &graph),
        decomposition.get_handle(decomposition.get_bound(chain, true, false), &graph)};
}

ChildChain chain_of_site(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& site) {
    return chain_ends(decomposition, graph, decomposition.get_parent(site));
}

NodePlacement placement_of_node(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                nid_t id) {
    // A node's parent is its chain, which is the node alone when it bounds no site. The chain's
    // parent is a site or the root.
    const net_handle_t chain = decomposition.get_parent(
        decomposition.get_net(graph.get_handle(id, false), &graph));
    NodePlacement out;
    out.site = decomposition.get_parent(chain);
    out.chain = chain_ends(decomposition, graph, chain);
    return out;
}

ChildChain chain_holding(const EditSnarlDecompositionOverlay& sites, const HandleGraph& graph,
                         const net_handle_t& net) {
    if (!sites.is_edit(net)) {
        return chain_of_site(sites, graph, net);
    }
    const handle_t start = sites.get_site_start(sites.canonical(net));
    return chain_ends(sites, graph,
                      sites.get_parent(sites.get_net(graph.forward(start), &graph)));
}

size_t chain_key_of(const HandleGraph& graph, const ChildChain& chain) {
    // Every node is in exactly one chain, so the chain's first node names it alone.
    return (size_t)graph.get_id(chain.start);
}

}
}
