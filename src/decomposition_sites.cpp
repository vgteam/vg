/** \file
 * Implements the free functions for walking the sites of a snarl
 * decomposition.
 */

#include "decomposition_sites.hpp"

#include <algorithm>
#include <iterator>

namespace vg {

using namespace std;

SiteBounds site_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph, const net_handle_t& net) {
    // Read in through the bound the traversal starts at, and out through the
    // bound it ends at.
    bool starts_at_end = decomposition.starts_at(net) == SnarlDecomposition::END;
    bool ends_at_end = decomposition.ends_at(net) == SnarlDecomposition::END;
    return {
        decomposition.get_handle(decomposition.get_bound(net, starts_at_end, true), &graph),
        decomposition.get_handle(decomposition.get_bound(net, ends_at_end, false), &graph)
    };
}

vector<net_handle_t> top_level_sites(const SnarlDecomposition& decomposition) {
    vector<net_handle_t> sites;
    decomposition.for_each_child(decomposition.get_root(), [&decomposition, &sites](const net_handle_t& chain) {
        std::ranges::copy(child_sites(decomposition, chain), back_inserter(sites));
    });
    return sites;
}

/// Return true if a chain holds at least one site.
static bool has_site(const SnarlDecomposition& decomposition, const net_handle_t& chain) {
    bool found = false;
    decomposition.for_each_child(chain, [&decomposition, &found](const net_handle_t& child) {
        found = decomposition.is_snarl(child);
        // Stop at the first site.
        return !found;
    });
    return found;
}

vector<net_handle_t> child_chains(const SnarlDecomposition& decomposition, const net_handle_t& site) {
    vector<net_handle_t> chains;
    decomposition.for_each_child(site, [&decomposition, &chains](const net_handle_t& chain) {
        if (has_site(decomposition, chain)) {
            chains.push_back(chain);
        }
    });
    return chains;
}

vector<net_handle_t> child_sites(const SnarlDecomposition& decomposition, const net_handle_t& chain) {
    vector<net_handle_t> sites;
    decomposition.for_each_child(chain, [&decomposition, &sites](const net_handle_t& child) {
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
