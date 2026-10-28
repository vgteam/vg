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

#include "../src/handle.hpp"

namespace vg {

using namespace std;
using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

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
 * Get the innermost site that holds the element `net`, or a root handle if no site does.
 *
 * For a site, this is its parent site. It works for any element, such as an edit, whose parent
 * need not be a chain.
 */
net_handle_t enclosing_site(const SnarlDecomposition& decomposition, const net_handle_t& net);

}

#endif
