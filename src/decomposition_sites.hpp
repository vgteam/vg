#ifndef VG_DECOMPOSITION_SITES_HPP_INCLUDED
#define VG_DECOMPOSITION_SITES_HPP_INCLUDED

/** \file
 * Free functions for walking the sites of a snarl decomposition.
 *
 * Here a site is any snarl other than the root and the trivial snarls. Each site sits in a chain,
 * and its child chains hold its child sites.
 */

#include <vector>

#include <handlegraph/snarl_decomposition.hpp>

#include "handle.hpp"

namespace vg {

using namespace std;
using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

/**
 * The two bounds of a site or chain, oriented like the start and end of a Snarl.
 *
 * `start` reads into the site and `end` reads out of it. Most sites lie between the inward sides
 * of `start` and `end`, which can be one node, as in a looping chain. A site made from an edit in
 * the reads lies inside nodes instead (`inside_node`). Its boundaries then lie inside the nodes of
 * `start` and `end`, and are not node sides. The decomposition gives only bounds between nodes.
 */
struct SiteBounds {
    handle_t start;
    handle_t end;
    /// Whether the site's boundaries lie inside the nodes of `start` and `end`.
    bool inside_node = false;
    /// For a site inside nodes, the offset of its first base along the forward strand of the node
    /// of `start`, whichever way `start` reads the node. 0 for any other site.
    uint32_t offset = 0;

    bool operator==(const SiteBounds& other) const = default;
};

/**
 * Get the bounds of a site or chain, read the way the net handle traverses it.
 *
 * For an end-to-start traversal, the start is the end bound read inward.
 *
 * @pre `net` is not the root, and does not start or end at a tip.
 */
SiteBounds site_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph, const net_handle_t& net);

/**
 * Return true if a snarl is a trivial snarl, which holds nothing to call.
 *
 * A snarl counts as trivial when it has no children and every realizable traversal of it goes
 * from one bound to the other.
 */
bool is_trivial_site(const SnarlDecomposition& decomposition, const net_handle_t& snarl);

/**
 * Get the sites in the chains of the root, in chain order, each read the way
 * its chain reads.
 */
vector<net_handle_t> top_level_sites(const SnarlDecomposition& decomposition);

/**
 * Get the chains directly inside a site that hold at least one site.
 */
vector<net_handle_t> child_chains(const SnarlDecomposition& decomposition, const net_handle_t& site);

/**
 * Get the sites in a chain, in the order and orientation of the given
 * traversal of the chain.
 */
vector<net_handle_t> child_sites(const SnarlDecomposition& decomposition, const net_handle_t& chain);

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
