#ifndef VG_DECOMPOSITION_SITES_HPP_INCLUDED
#define VG_DECOMPOSITION_SITES_HPP_INCLUDED

/** \file
 * Free functions for walking the sites of a snarl decomposition.
 *
 * A site is a snarl other than the root, so every site has a start bound and
 * an end bound. Each site sits in a chain, and its child chains hold its
 * child sites.
 *
 * A SnarlDistanceIndex can group root chains into root snarls, which have no
 * bounds. These functions do not support decompositions with root snarls.
 */

#include <vector>

#include <handlegraph/snarl_decomposition.hpp>

#include "handle.hpp"

namespace vg {

using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

/**
 * The two bounds of a site or chain, oriented the way a Snarl's start and
 * end Visits are: the start handle reads into it, and the end handle reads
 * out of it.
 */
struct SiteBounds {
    handle_t start;
    handle_t end;

    bool operator==(const SiteBounds& other) const = default;
};

/**
 * Get the bounds of a traversal of a site or chain: the start handle reads
 * in through the bound the traversal starts at, and the end handle reads out
 * through the bound it ends at. Must not be called on the root, or on a
 * traversal that starts or ends at an internal tip.
 */
SiteBounds site_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph, const net_handle_t& net);

/**
 * Get the sites in the chains of the root, in chain order, each read the way
 * its chain reads.
 */
std::vector<net_handle_t> top_level_sites(const SnarlDecomposition& decomposition);

/**
 * Get the chains directly inside a site that hold at least one site.
 */
std::vector<net_handle_t> child_chains(const SnarlDecomposition& decomposition, const net_handle_t& site);

/**
 * Get the sites in a chain, in the order and orientation of the given
 * traversal of the chain.
 */
std::vector<net_handle_t> child_sites(const SnarlDecomposition& decomposition, const net_handle_t& chain);

/**
 * Get the site that directly encloses a site, start-to-end, or a root handle
 * if the site is top-level.
 */
net_handle_t parent_site(const SnarlDecomposition& decomposition, const net_handle_t& site);

/**
 * Count the sites enclosing a site. A top-level site has depth 0.
 */
size_t depth(const SnarlDecomposition& decomposition, const net_handle_t& site);

/**
 * Return true if a site has no child sites.
 */
bool is_leaf(const SnarlDecomposition& decomposition, const net_handle_t& site);

}

#endif
