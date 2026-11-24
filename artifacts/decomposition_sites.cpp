/** \file
 * Implements the free functions for walking the sites of a snarl
 * decomposition.
 */

#include "decomposition_sites.hpp"

namespace vg {

using namespace std;

bool is_trivial_site(const SnarlDecomposition& decomposition, const net_handle_t& snarl) {
    if (!decomposition.for_each_child(snarl, [](const net_handle_t&) { return false; })) {
        // The iteration stopped at a child.
        return false;
    }
    return decomposition.for_each_traversal(snarl, [&](const net_handle_t& traversal) {
        SnarlDecomposition::endpoint_t start = decomposition.starts_at(traversal);
        SnarlDecomposition::endpoint_t end = decomposition.ends_at(traversal);
        return (start == SnarlDecomposition::START && end == SnarlDecomposition::END) ||
            (start == SnarlDecomposition::END && end == SnarlDecomposition::START);
    });
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
            return !decomposition.is_snarl(child) || is_trivial_site(decomposition, child);
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
        if (decomposition.is_snarl(child) && !is_trivial_site(decomposition, child)) {
            sites.push_back(child);
        }
    });
    return sites;
}

net_handle_t parent_site(const SnarlDecomposition& decomposition, const net_handle_t& site) {
    // A site's parent is its chain, and the chain's parent is a site or the root.
    return decomposition.get_parent(decomposition.get_parent(site));
}

net_handle_t enclosing_site(const SnarlDecomposition& decomposition, const net_handle_t& net) {
    net_handle_t here = decomposition.get_parent(net);
    while (!decomposition.is_root(here) && !decomposition.is_snarl(here)) {
        here = decomposition.get_parent(here);
    }
    return here;
}

}
