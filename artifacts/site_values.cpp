#include "site_values.hpp"
#include "utility.hpp"

#include <algorithm>

namespace vg {
const StepSequences* own_sequences(const AlleleSequences* sequences, size_t allele) {
    if (sequences == nullptr || allele >= sequences->size() || (*sequences)[allele].empty()) {
        return nullptr;
    }
    return &(*sequences)[allele];
}

size_t bases_before_site(const HandleGraph& graph, const SiteBounds& site) {
    if (!site.inside_node) {
        return 0;
    }
    return graph.get_is_reverse(site.start) ? graph.get_length(site.start) - 1 - site.offset
                                            : site.offset;
}

EditSite edit_site(const HandleGraph& graph, const Traversal& steps, uint32_t offset, size_t at,
                   size_t length, const vector<string>& alts) {
    EditSite site;
    site.bounds.start = steps.front();
    site.bounds.end = steps.back();
    site.bounds.inside_node = true;
    site.bounds.offset = offset;
    site.travs.assign(alts.size() + 1, steps);
    site.sequences.resize(alts.size() + 1);
    const size_t host = length > 0 ? at : at - 1;
    for (size_t a = 0; a < alts.size(); ++a) {
        size_t start = 0;
        for (const handle_t& h : steps) {
            string seq = graph.get_sequence(h);
            const size_t end = start + seq.size();
            const size_t from = std::max(start, at);
            if (const size_t to = std::min(end, at + length); from < to) {
                seq.erase(from - start, to - from);
            }
            if (host >= start && host < end) {
                seq.insert(at - start, alts[a]);
            }
            site.sequences[a + 1].push_back(std::move(seq));
            start = end;
        }
    }
    return site;
}

bool inside_one_node(const SiteBounds& site) {
    return site.inside_node && site.start == site.end;
}

string edit_site_id(nid_t node, uint32_t offset, char kind, const string& alt) {
    return "op" + std::to_string(node) + "." + std::to_string(offset) + "." + kind + alt;
}

string allele_sequence(const HandleGraph& graph, const Traversal& walk, const StepSequences* own) {
    string sequence;
    if (own != nullptr) {
        for (const string& step : *own) {
            sequence += step;
        }
    } else {
        for (const handle_t& handle : walk) {
            sequence += graph.get_sequence(handle);
        }
    }
    return sequence;
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

SiteBounds oriented_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                           const net_handle_t& net) {
    // get_bound() ignores the traversal, so this is the start bound read in and the end bound
    // read out.
    return SiteBounds{decomposition.get_handle(decomposition.get_bound(net, false, true), &graph),
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
                           const net_handle_t& site, const SiteBounds& as_called,
                           vector<net_handle_t>* child_sites) {
    SiteChildren out;
    // The caller has this site as the decomposition orients it, or turned round. Compare node IDs
    // for the first, and whole handles for the second, since turning round flips the bounds as
    // well as swapping them.
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
