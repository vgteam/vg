#ifndef VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED
#define VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED

/** \file
 * Defines SnarlManagerDecomposition, which presents a SnarlManager through
 * the handlegraph::SnarlDecomposition interface.
 */

#include <handlegraph/snarl_decomposition.hpp>

#include <functional>
#include <mutex>
#include <optional>
#include <unordered_map>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"

namespace vg {

using handlegraph::net_handle_t;

/**
 * A view of a SnarlManager as a SnarlDecomposition. The manager must have
 * had finish() called, must not change while the view exists, and must
 * outlive the view, as must the graph it was built on.
 *
 * The view presents the snarl tree that a SnarlDistanceIndex built from the
 * same snarls presents, with snarls and chains oriented as the manager
 * orients them, except for these differences:
 *
 * - A snarl that holds no nodes and is not an ultrabubble, such as one whose
 *   bounding node has a reversing self loop, is shown as a snarl, whereas a
 *   distance index hides it like a trivial snarl.
 * - A snarl whose start and end Visits are the same, so that it leads out of
 *   a node and around back into it, has that node as both bounds, and makes
 *   a cyclic chain by itself.
 * - A cyclic chain, whose last snarl leads back into its first node, has the
 *   same node as its start and end bound, and lists that node once.
 * - Every chain in the root is a child of the root itself, whereas a
 *   distance index can group root chains into root snarls that have no
 *   bounding nodes.
 * - Going left from a node to the neighbouring node in its chain, across a
 *   hidden trivial snarl, follow_net_edges() reads the node it reaches the
 *   way it travels, whereas a distance index reads it the way the starting
 *   node is read.
 *
 * Trivial snarls, which are ultrabubbles that hold no nodes, are hidden, as
 * in a distance index: a chain lists the two bounding nodes of a trivial
 * snarl as neighbours.
 *
 * A node that bounds no snarl is shown as a trivial chain holding only that
 * node, whose start-to-end traversal reads the node forward.
 *
 * The view cannot present a unary snarl, whose end Visit is its start Visit
 * reversed: leaving its chain through the start bound and through the end
 * bound would read the same handle, so the two bounds could not be told
 * apart. The constructor throws if the manager holds one.
 *
 * The manager does not record internal tips. Apart from the root's
 * tip-to-tip traversal, the view never produces a traversal that starts or
 * ends at an internal tip, and it throws when asked for children that start
 * at internal tips.
 *
 * Methods that take a graph use it to read and make handles, so it must have
 * the manager's node IDs. follow_net_edges() ignores the graph it is
 * given and walks the graph the view was made with.
 */
class SnarlManagerDecomposition : public handlegraph::SnarlDecomposition {
public:
    /// Make a view of the given manager over the given graph.
    SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph);

    ~SnarlManagerDecomposition() override = default;

    /// Get a start-to-end net handle for a snarl the manager owns. For a
    /// hidden trivial snarl, the handle names a snarl that the tree does not
    /// contain.
    net_handle_t get_snarl_net(const Snarl* snarl) const;

    /// Get the snarl the manager owns for a net handle to a snarl or to one
    /// of its bound sentinels.
    const Snarl* get_snarl(const net_handle_t& net) const;

    ////////////////////////////////////////////////////////////////////////
    // SnarlDecomposition interface
    ////////////////////////////////////////////////////////////////////////

    net_handle_t get_root() const override;
    bool is_root(const net_handle_t& net) const override;
    bool is_snarl(const net_handle_t& net) const override;
    bool is_chain(const net_handle_t& net) const override;
    bool is_node(const net_handle_t& net) const override;
    bool is_sentinel(const net_handle_t& net) const override;

    net_handle_t get_net(const handle_t& handle, const HandleGraph* handle_graph) const override;

    /// Get the handle for a node, for a trivial chain's node, or for a snarl
    /// bound sentinel. A sentinel gives its bounding node read into the snarl
    /// if the sentinel faces in, and out of the snarl if it faces out.
    handle_t get_handle(const net_handle_t& net, const HandleGraph* handle_graph) const override;

    net_handle_t get_parent(const net_handle_t& child) const override;
    net_handle_t get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const override;
    net_handle_t flip(const net_handle_t& net) const override;
    net_handle_t canonical(const net_handle_t& net) const override;
    endpoint_t starts_at(const net_handle_t& traversal) const override;
    endpoint_t ends_at(const net_handle_t& traversal) const override;

    /// Get the traversal of a snarl between two of its bound sentinels, or of
    /// a chain between two of its bounding node traversals. Throws if the two
    /// are not both sentinels of one snarl or both bounds of one chain or
    /// trivial chain.
    net_handle_t get_parent_traversal(const net_handle_t& traversal_start,
                                      const net_handle_t& traversal_end) const override;

protected:
    bool for_each_child_impl(const net_handle_t& traversal,
                             const std::function<bool(const net_handle_t&)>& iteratee) const override;
    bool for_each_traversal_impl(const net_handle_t& item,
                                 const std::function<bool(const net_handle_t&)>& iteratee) const override;
    /// Whichever way we go, what we reach is read the way we travel.
    bool follow_net_edges_impl(const net_handle_t& here, const HandleGraph* handle_graph, bool go_left,
                               const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /// Throws, because the manager does not record internal tips.
    bool for_each_tippy_child_impl(const net_handle_t&,
                                   const std::function<bool(const net_handle_t&)>&) const override;

private:

    /**
     * The bounding nodes of the manager's chains, indexed along each chain.
     *
     * A Chain lists (snarl, backward) pairs, where backward means the snarl
     * reads end-to-start going along the chain. A snarl's rank is its
     * position in its chain. Node index i is the node
     * read into the snarl of rank i, going along the chain, and node index
     * chain.size() is the node read out of the last snarl. In a cyclic chain,
     * node index chain.size() is node index 0 again, so the chain has only
     * chain.size() distinct nodes, and node index 0 is both its start bound
     * and its end bound.
     *
     * At a chain's bound nodes, a handle reading along the chain enters it
     * through its START endpoint or leaves through its END endpoint, and one
     * reading against the chain enters through END or leaves through START.
     *
     * A node that bounds any snarls bounds snarls of only one chain.
     */
    class ChainLayout {
    public:
        ChainLayout(const SnarlManager& manager, const HandleGraph& graph);

        /// Where a handle to a bounding node lies in its chain.
        struct Place {
            const Chain* chain;
            /// The node's index in the chain, less than node_count().
            size_t index;
            /// True if the handle reads the same way as the chain.
            bool along;
        };

        /// Get the chain's node at the given index, read the way the chain
        /// reads.
        handle_t node(const Chain& chain, size_t index) const;
        /// Return true if the node read out of the chain's last snarl is its
        /// first node, read the same way.
        bool is_cyclic(const Chain& chain) const;
        /// Count the chain's distinct bounding nodes: chain.size() + 1, or
        /// chain.size() if the chain is cyclic.
        size_t node_count(const Chain& chain) const;
        /// Get the index of the chain's end bound node, which is 0 if the
        /// chain is cyclic.
        size_t last_index(const Chain& chain) const;
        /// Find where a handle's node lies in the chain whose snarls it
        /// bounds, or nothing if it bounds no snarl.
        std::optional<Place> locate(const handle_t& handle) const;
        /// Get the rank of the snarl that a handle at the given place reads
        /// into, or nothing if it reads out of the chain.
        std::optional<size_t> snarl_ahead(const Place& place) const;
        /// Get the endpoint of the chain that a handle at the given place
        /// reads into the chain through, or nothing if the handle does not
        /// read into the chain at an endpoint.
        std::optional<endpoint_t> entry_endpoint(const Place& place) const;
        /// Get the endpoint of the chain that a handle at the given place
        /// reads out of the chain through, or nothing if the handle does not
        /// read out of the chain at an endpoint.
        std::optional<endpoint_t> exit_endpoint(const Place& place) const;

    private:
        const SnarlManager& manager;
        const HandleGraph& graph;
    };

    /**
     * The nodes that bound no snarl, which the view shows as trivial chains,
     * grouped by the snarl directly containing them. Nodes in no snarl are in the
     * root.
     */
    class TrivialChains {
    public:
        TrivialChains(const SnarlManager& manager, const HandleGraph& graph);

        /// Get the trivial chain nodes directly in a snarl, or in the root if
        /// given null, in ID order.
        std::vector<nid_t> nodes_in(const Snarl* parent) const;
        /// Get the snarl directly containing a trivial chain's node, or null
        /// for the root.
        const Snarl* parent_of(nid_t id) const;

    private:
        bool bounds_no_snarl(nid_t id) const;
        /// Fill in parents and root_nodes, on the first call only.
        void index_parents() const;

        const SnarlManager& manager;
        const HandleGraph& graph;

        mutable std::once_flag parents_indexed;
        /// For every trivial chain node, its parent snarl, or null for the
        /// root.
        mutable std::unordered_map<nid_t, const Snarl*> parents;
        /// The trivial chain nodes in the root, in ID order.
        mutable std::vector<nid_t> root_nodes;
    };

    /// Get a snarl's net handle, start-to-end if forward and end-to-start
    /// otherwise.
    net_handle_t snarl_net(const Snarl* snarl, bool forward) const;
    net_handle_t chain_net(const Chain& chain) const;
    const Chain& chain_of_net(const net_handle_t& net) const;
    /// Get the net handle for a parent snarl, or the root if given null.
    net_handle_t parent_snarl_net(const Snarl* parent) const;
    bool is_hidden(const Snarl* snarl) const;

    /// Get the net handle for the chain's node at the given index, read along
    /// the chain or against it.
    net_handle_t chain_node_net(const Chain& chain, size_t index, bool along) const;
    /// Get the net handle for the chain's snarl at the given rank, read along
    /// the chain or against it.
    net_handle_t chain_snarl_net(const Chain& chain, size_t rank, bool along) const;
    /// List the child chains of a snarl, or of the root if given null: the
    /// manager's chains in the manager's order, then the trivial chains in ID
    /// order.
    std::vector<net_handle_t> snarl_children(const Snarl* parent) const;
    /// List a chain's nodes and shown snarls in the order and orientation
    /// that the given traversal of the chain reads them.
    std::vector<net_handle_t> chain_children(const net_handle_t& traversal) const;

    /// Get the next child of a chain, reached by leaving a snarl or node in
    /// it through the end of its traversal, or nothing if that leaves the
    /// chain.
    std::optional<net_handle_t> next_in_chain(const net_handle_t& here) const;
    /// Get what a handle, reached by following an edge from a child of a
    /// snarl, reads into among that snarl's children: the sibling chain it
    /// enters, or the bound sentinel facing out if it reads out of the snarl.
    /// Gets nothing if it reads into the middle of a chain. With a null
    /// snarl, the handle moves among the root's children.
    std::optional<net_handle_t> sibling_entered_by(const handle_t& handle, const Snarl* parent) const;
    /// Get what we reach by following each graph edge out of the end of a
    /// chain's or trivial chain's traversal, or out of a bound sentinel that
    /// faces into its snarl.
    std::vector<net_handle_t> siblings_ahead(const net_handle_t& here) const;

    const SnarlManager& manager;
    const HandleGraph& graph;
    ChainLayout layout;
    TrivialChains trivial_chains;

    /// The manager's snarls, by SnarlManager::snarl_number().
    std::vector<const Snarl*> snarls_by_number;
    /// Whether each snarl, by SnarlManager::snarl_number(), is hidden.
    std::vector<bool> hidden;
};

}

#endif
