#ifndef VG_DECOMPOSITION_SITES_HPP_INCLUDED
#define VG_DECOMPOSITION_SITES_HPP_INCLUDED

/** \file
 * Free functions for walking the sites of a snarl decomposition.
 *
 * A site is a snarl in the decomposition, for which is_snarl() is true. Every
 * site has a start bound and an end bound, sits in a chain, and holds child
 * chains that hold its child sites. The root is not a site, and neither are
 * the root snarls that a SnarlDistanceIndex groups root chains under, since
 * is_root() is true for them.
 */

#include <vector>

#include <handlegraph/snarl_decomposition.hpp>

#include "handle.hpp"

namespace vg {

using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

/**
 * @brief The two bounds of a site or chain, as graph handles that read in
 * the direction of a traversal from one bound to the other.
 */
struct SiteBounds {
    /// @brief The handle reading into the site or chain at its start.
    handle_t start;
    /// @brief The handle reading out of the site or chain at its end.
    handle_t end;

    bool operator==(const SiteBounds& other) const = default;
};

/**
 * @brief Get the bounds of a traversal of a site or chain.
 *
 * @return The handle reading in through the bound where the traversal
 * starts, and the handle reading out through the bound where it ends. For a
 * traversal that starts and ends at the same bound, the two read the same
 * node in opposite orientations.
 *
 * @pre The traversal neither starts nor ends at a tip, and is not of the
 * root.
 */
SiteBounds site_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph, const net_handle_t& net);

/**
 * @brief Get the sites in the root's chains.
 *
 * @post Sites come chain by chain, in the order each chain lists them, each
 * read the way its chain reads.
 */
std::vector<net_handle_t> top_level_sites(const SnarlDecomposition& decomposition);

/**
 * @brief Get the chains directly inside a site that hold at least one site.
 */
std::vector<net_handle_t> child_chains(const SnarlDecomposition& decomposition, const net_handle_t& site);

/**
 * @brief Get the sites in a chain.
 *
 * @post For a start-to-end or end-to-start traversal of the chain, sites
 * come in the order and orientation that the traversal reads them. For any
 * other traversal, they come start-to-end, in an order the decomposition
 * picks.
 */
std::vector<net_handle_t> child_sites(const SnarlDecomposition& decomposition, const net_handle_t& chain);

/**
 * @brief Get the site that directly encloses a site.
 *
 * @return The enclosing site, start-to-end, or a net handle for which
 * is_root() is true if the site is top-level.
 */
net_handle_t parent_site(const SnarlDecomposition& decomposition, const net_handle_t& site);

/**
 * @brief Count the sites enclosing a site.
 *
 * @return 0 for a top-level site.
 */
size_t depth(const SnarlDecomposition& decomposition, const net_handle_t& site);

/**
 * @brief Return true if a site has no child sites.
 */
bool is_leaf(const SnarlDecomposition& decomposition, const net_handle_t& site);

}

#endif
