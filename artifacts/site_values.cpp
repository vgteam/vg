#include "site_values.hpp"
#include "utility.hpp"

#include <algorithm>

namespace vg {
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
    // Find the node holding base `position` of the walk, and the bases of the walk before that
    // node.
    auto holding = [&](size_t position) {
        size_t before = 0;
        for (const handle_t& h : steps) {
            const size_t length = graph.get_length(h);
            if (position < before + length) {
                return make_pair(h, before);
            }
            before += length;
        }
        return make_pair(steps.back(), before - graph.get_length(steps.back()));
    };
    // The boundary after the last changed base, or, for an insertion, after the base before it.
    const auto [end_node, end_before] = holding(at + length - 1);
    site.change_end = end_node;
    site.change_end_offset = at + length - end_before;
    if (length > 0) {
        const auto [start_node, start_before] = holding(at);
        site.change_start = start_node;
        site.change_start_offset = at - start_before;
    } else {
        site.change_start = site.change_end;
        site.change_start_offset = site.change_end_offset;
    }
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
                            vector<ChildSite>* child_sites) {
    const ChildChain bounds = oriented_bounds(decomposition, graph, chain);
    if (child_sites != nullptr) {
        for (const net_handle_t& child : vg::child_sites(decomposition, chain)) {
            child_sites->push_back(
                ChildSite{child, oriented_bounds(decomposition, graph, child), bounds});
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
                           vector<ChildSite>* child_sites) {
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

NodePlacement placement_of_node(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                nid_t id) {
    // A node's parent is its chain, which is the node alone when it bounds no site. The chain's
    // parent is a site or the root.
    const net_handle_t chain = decomposition.get_parent(
        decomposition.get_net(graph.get_handle(id, false), &graph));
    NodePlacement out;
    out.site = decomposition.get_parent(chain);
    out.chain = oriented_bounds(decomposition, graph, chain);
    if (!decomposition.is_root(out.site)) {
        out.enclosing.push_back(oriented_bounds(decomposition, graph, out.site));
        const vector<SiteBounds> above = enclosing_sites(decomposition, graph, out.site);
        out.enclosing.insert(out.enclosing.end(), above.begin(), above.end());
    }
    return out;
}

size_t chain_key_of(const HandleGraph& graph, const ChildChain& chain) {
    // Every node is in exactly one chain, so the chain's first node names it alone.
    return (size_t)graph.get_id(chain.start);
}

}
}
