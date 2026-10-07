/** \file
 * Implements the free functions for walking the sites of a snarl
 * decomposition.
 */

#include "decomposition_sites.hpp"

namespace vg {

using namespace std;

SiteBounds site_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph, const net_handle_t& net) {
    // The bounds ignore the traversal, so they come out start-to-end.
    SiteBounds bounds {
        decomposition.get_handle(decomposition.get_bound(net, false, true), &graph),
        decomposition.get_handle(decomposition.get_bound(net, true, false), &graph)
    };
    if (decomposition.starts_at(net) == SnarlDecomposition::END) {
        bounds = {graph.flip(bounds.end), graph.flip(bounds.start)};
    }
    return bounds;
}

vector<net_handle_t> top_level_sites(const SnarlDecomposition& decomposition) {
    vector<net_handle_t> sites;
    decomposition.for_each_child(decomposition.get_root(), [&](const net_handle_t& chain) {
        for (const net_handle_t& site : child_sites(decomposition, chain)) {
            sites.push_back(site);
        }
    });
    return sites;
}

vector<net_handle_t> child_chains(const SnarlDecomposition& decomposition, const net_handle_t& site) {
    vector<net_handle_t> chains;
    decomposition.for_each_child(site, [&](const net_handle_t& chain) {
        bool has_site = !decomposition.for_each_child(chain, [&](const net_handle_t& child) {
            return !decomposition.is_snarl(child);
        });
        if (has_site) {
            chains.push_back(chain);
        }
    });
    return chains;
}

vector<net_handle_t> child_sites(const SnarlDecomposition& decomposition, const net_handle_t& chain) {
    vector<net_handle_t> sites;
    decomposition.for_each_child(chain, [&](const net_handle_t& child) {
        if (decomposition.is_snarl(child)) {
            sites.push_back(child);
        }
    });
    return sites;
}

net_handle_t parent_site(const SnarlDecomposition& decomposition, const net_handle_t& site) {
    // A site's parent is its chain, and the chain's parent is a site or the root.
    return decomposition.get_parent(decomposition.get_parent(site));
}

size_t depth(const SnarlDecomposition& decomposition, const net_handle_t& site) {
    size_t enclosing = 0;
    for (net_handle_t parent = parent_site(decomposition, site);
         !decomposition.is_root(parent);
         parent = parent_site(decomposition, parent)) {
        enclosing++;
    }
    return enclosing;
}

bool is_leaf(const SnarlDecomposition& decomposition, const net_handle_t& site) {
    return child_chains(decomposition, site).empty();
}

}
