#ifndef VG_UNITTEST_SITE_TREE_QUERIES_HPP_INCLUDED
#define VG_UNITTEST_SITE_TREE_QUERIES_HPP_INCLUDED

/** \file
 * Queries about a site's place in a snarl decomposition that only the tests ask.
 */

#include "../../decomposition_sites.hpp"

namespace vg {
namespace unittest {

using namespace std;

/**
 * Get the bounds of a site or chain, read the way the net handle traverses it.
 *
 * For an end-to-start traversal, the start is the end bound read inward.
 *
 * @pre `net` is not the root, and does not start or end at a tip.
 */
inline pair<handle_t, handle_t> site_bounds(const SnarlDecomposition& decomposition,
                                            const HandleGraph& graph, const net_handle_t& net) {
    // get_bound() ignores the traversal, so these come out start-to-end.
    const handle_t start =
        decomposition.get_handle(decomposition.get_bound(net, false, true), &graph);
    const handle_t end =
        decomposition.get_handle(decomposition.get_bound(net, true, false), &graph);
    if (decomposition.starts_at(net) == SnarlDecomposition::END) {
        return {graph.flip(end), graph.flip(start)};
    }
    return {start, end};
}

/// Count the sites enclosing a site. A top-level site has depth 0.
inline size_t depth(const SnarlDecomposition& decomposition, const net_handle_t& site) {
    size_t enclosing = 0;
    for (net_handle_t parent = parent_site(decomposition, site); !decomposition.is_root(parent);
         parent = parent_site(decomposition, parent)) {
        enclosing++;
    }
    return enclosing;
}

/// Return true if a site has no child sites.
inline bool is_leaf(const SnarlDecomposition& decomposition, const net_handle_t& site) {
    return child_chains(decomposition, site).empty();
}

}
}

#endif
