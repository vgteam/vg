#ifndef VG_EDIT_SNARL_DECOMPOSITION_OVERLAY_HPP_INCLUDED
#define VG_EDIT_SNARL_DECOMPOSITION_OVERLAY_HPP_INCLUDED

/** \file
 * A snarl decomposition with edits laid over it.
 */

#include <functional>
#include <mutex>
#include <unordered_map>
#include <vector>

#include <handlegraph/snarl_decomposition.hpp>

#include "grow_only_array.hpp"
#include "../src/handle.hpp"

namespace vg {

using namespace std;
using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

/**
 * A snarl decomposition with edits laid over it.
 *
 * An edit is a site that the decomposition does not have: a stretch of sequence, inside one node or
 * across several nodes, where the reads disagree with the graph. Every site, a snarl or an edit, is
 * a net handle, and has bounds: at each end, a node handle and an offset into it.
 *
 * Edits are leaves of the tree. Each is a child of the lowest node, chain or snarl that holds the
 * nodes of both its bounds, and `for_each_child` lists it after that element's own
 * children. An edit is not a step of any chain, so no net edges lead to it or from it. Edits may
 * overlap each other, and may run across nodes and snarls of the chain they are a child of.
 *
 * An edit is none of the kinds of element `SnarlDecomposition` knows: `is_snarl`, `is_chain`,
 * `is_node`, `is_sentinel` and `is_root` are all false for it, and `is_edit` is true.
 *
 * A merged site is the other kind of site the overlay adds: a run of neighbouring sites of one
 * chain, taken as one site. Like a snarl, it starts after the last base of its start node and ends
 * before the first base of its end node. Its parent is the lowest element that holds both nodes,
 * which is the chain of the sites it is made of. That element lists it after its own children
 * and among its edits, in the order they were added, and still lists the sites the merged site
 * is made of. A merged site has no children and no net edges of its own. `is_merged_site` is
 * true for it, and `is_edit`, like the predicates of `SnarlDecomposition`, is false.
 *
 * Edits and merged sites can only be added, and are never changed, so no net handle is ever
 * invalidated.
 *
 * Several threads may add edits and merged sites and read the overlay at once. Reading an edit through its net
 * handle takes no lock. Listing an element's children takes a lock, and lists the edits added
 * before the listing began. It may or may not list an edit another thread adds during it.
 *
 * The overlay's net handles hold their content shifted up one bit, with the lowest bit telling
 * the two kinds apart: 0 for an element of the inner decomposition, whose net handle is the rest,
 * and 1 for an edit, whose number and orientation are the rest. Overlays laid over one another
 * this way each take the next bit up.
 *
 * @pre The inner decomposition's net handles leave their top bit clear, so that shifting them up
 * loses nothing.
 */
class EditSnarlDecompositionOverlay : public SnarlDecomposition {
public:
    /// Lay an overlay with no edits over `inner`, a decomposition of `graph`.
    EditSnarlDecompositionOverlay(const SnarlDecomposition& inner, const HandleGraph& graph);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Edits
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Add an edit, and return its net handle, which reads it start to end.
     *
     * The edit becomes a child of the lowest node, chain or snarl that holds both `start` and
     * `end`.
     *
     * @param start The node holding the start boundary, oriented to read into the edit.
     * @param start_offset The number of bases of `start` before the start boundary.
     * @param end The node holding the end boundary, oriented to read out of the edit.
     * @param end_offset The number of bases of `end` before the end boundary.
     * @throws std::invalid_argument If no element of the decomposition holds both nodes, or an
     * offset is past the end of its node.
     */
    net_handle_t add_edit(const handle_t& start, size_t start_offset, const handle_t& end,
                          size_t end_offset);

    /// Return true if `net` is an edit.
    bool is_edit(const net_handle_t& net) const;

    /**
     * Add a merged site, and return its net handle, which reads it start to end.
     *
     * @param start The node the site starts after, oriented to read into the site.
     * @param end The node the site ends before, oriented to read out of the site.
     */
    net_handle_t add_merged_site(const handle_t& start, const handle_t& end);

    /// Return true if `net` is a merged site.
    bool is_merged_site(const net_handle_t& net) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Bounds of sites, for snarls, edits and merged sites
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Get the node that holds the start boundary of a site, oriented to read into the site.
     *
     * The site is read the way `site` traverses it. For a snarl, this is the start node.
     *
     * @pre `site` is a start-to-end or end-to-start traversal of a snarl, an edit, or a merged
     * site.
     */
    handle_t get_site_start(const net_handle_t& site) const;

    /// Get the node that holds the end boundary of a site, oriented to read out of the site.
    handle_t get_site_end(const net_handle_t& site) const;

    /**
     * Get the number of bases of `get_site_start(site)` that come before the start boundary.
     *
     * For a snarl or a merged site, the boundary is the inward side of its start node, so this
     * is that node's length.
     */
    size_t get_site_start_offset(const net_handle_t& site) const;

    /**
     * Get the number of bases of `get_site_end(site)` that come before the end boundary.
     *
     * For a snarl or a merged site, the boundary is the inward side of its end node, so this is
     * 0.
     */
    size_t get_site_end_offset(const net_handle_t& site) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // SnarlDecomposition, passed through to the inner decomposition except for edits and merged
    // sites, which the comments below call edits alike
    //
    // Each converts net handles of the inner decomposition to and from the overlay's, going in and
    // coming out, including those it gives to an iteratee.
    ////////////////////////////////////////////////////////////////////////////////////////////////

    net_handle_t get_root() const override;
    bool is_root(const net_handle_t& net) const override;
    bool is_snarl(const net_handle_t& net) const override;
    bool is_chain(const net_handle_t& net) const override;
    bool is_node(const net_handle_t& net) const override;
    bool is_sentinel(const net_handle_t& net) const override;
    net_handle_t get_net(const handle_t& handle, const HandleGraph* node_graph) const override;
    /// May not be called on an edit.
    handle_t get_handle(const net_handle_t& net, const HandleGraph* node_graph) const override;
    /// For an edit, get the element it is a child of.
    net_handle_t get_parent(const net_handle_t& child) const override;
    /// May not be called on an edit, which has no sentinels.
    net_handle_t get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const override;
    /// For an edit, get it read the other way, with its bounds swapped and flipped.
    net_handle_t flip(const net_handle_t& net) const override;
    net_handle_t canonical(const net_handle_t& net) const override;
    endpoint_t starts_at(const net_handle_t& traversal) const override;
    endpoint_t ends_at(const net_handle_t& traversal) const override;
    /// May not be called on edits.
    net_handle_t get_parent_traversal(const net_handle_t& traversal_start,
                                      const net_handle_t& traversal_end) const override;

protected:
    /// List the inner decomposition's children of `traversal`, then its edits and merged sites.
    bool for_each_child_impl(const net_handle_t& traversal,
                             const std::function<bool(const net_handle_t&)>& iteratee) const override;
    /// An edit has two traversals: start to end and end to start.
    bool for_each_traversal_impl(const net_handle_t& item,
                                 const std::function<bool(const net_handle_t&)>& iteratee) const override;
    /// An edit has no net edges.
    bool follow_net_edges_impl(const net_handle_t& here, const HandleGraph* node_graph,
                               bool go_left,
                               const std::function<bool(const net_handle_t&)>& iteratee) const override;

private:
    /// One edit, read start to end.
    struct Edit {
        handle_t start;
        size_t start_offset;
        handle_t end;
        size_t end_offset;
        /// The overlay's canonical net handle for the element the edit is a child of.
        net_handle_t parent;
        /// Whether this is a merged site, and not an edit.
        bool merged;
    };

    /// Return true if `net` is an edit or a merged site: a site this overlay added.
    bool is_added(const net_handle_t& net) const;

    /// Get the edit `net`, read the way `net` reads.
    ///
    /// Read the other way, the end boundary becomes the start. The bases before it along the
    /// flipped node are the ones after it along the node.
    Edit edit_of(const net_handle_t& net) const;

    /// Get the overlay's canonical net handle for the lowest element of the inner decomposition
    /// that holds the nodes of both `a` and `b`.
    net_handle_t lowest_holder(const handle_t& a, const handle_t& b) const;

    /// The decomposition the edits are laid over.
    const SnarlDecomposition& inner;
    /// The graph that `inner` decomposes.
    const HandleGraph& graph;

    /// The edits and the merged sites, by number.
    GrowOnlyArray<Edit> edits;
    /// The numbers of the edits and merged sites that are children of each element, by the
    /// overlay's canonical net handle for the element, in the order they were added.
    unordered_map<net_handle_t, vector<size_t>> children_of;
    /// Guards `children_of`.
    mutable std::mutex children_mutex;
};

}

#endif
