#ifndef VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED
#define VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED

/** \file
 * Defines SnarlManagerDecomposition, which presents a SnarlManager through
 * the handlegraph::SnarlDecomposition interface.
 */

#include <handlegraph/snarl_decomposition.hpp>

#include <mutex>
#include <unordered_map>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"

namespace vg {

using namespace std;
using handlegraph::net_handle_t;

/**
 * A view of a finished SnarlManager as a SnarlDecomposition, using only the
 * manager's public interface. The manager and the graph it was built on must
 * outlive the view, and the manager must not change while the view exists.
 *
 * The tree it presents is the one a SnarlDistanceIndex built from the same
 * snarls would present, with these differences:
 *
 * - A trivial snarl (an ultrabubble that holds no nodes and no child snarls)
 *   is hidden: its chain shows its two bounding nodes as neighbours. Every
 *   other snarl the manager holds is shown, including snarls that hold no
 *   nodes but are not ultrabubbles.
 * - A unary snarl, whose start and end are the same node, is shown as a snarl
 *   with that node as both bounds.
 * - A cyclic chain, whose last snarl leads back into its first node, has the
 *   same node as its start and end bound, and lists that node once.
 * - Chains in the root are never grouped under further snarls without bounds.
 * - Snarls and chains are oriented as the manager orients them.
 *
 * Nodes inside a snarl that bound no snarl are shown as trivial chains of one
 * node each, read forward. The manager does not record internal tips, so
 * traversals that start or end at internal tips are never produced, and
 * asking for children that start at internal tips throws.
 */
class SnarlManagerDecomposition : public handlegraph::SnarlDecomposition {
public:
    /// Make a view of the given finished manager over the given graph.
    SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph);

    virtual ~SnarlManagerDecomposition() = default;

    /// Get a start-to-end net handle for a snarl the manager owns. The snarl
    /// must not be a hidden trivial snarl.
    net_handle_t get_snarl_net(const Snarl* snarl) const;

    /// Get the snarl the manager owns for a net handle to a snarl or to one
    /// of its bound sentinels.
    const Snarl* get_snarl(const net_handle_t& net) const;

    ////////////////////////////////////////////////////////////////////////
    // SnarlDecomposition interface
    ////////////////////////////////////////////////////////////////////////

    virtual net_handle_t get_root() const;
    virtual bool is_root(const net_handle_t& net) const;
    virtual bool is_snarl(const net_handle_t& net) const;
    virtual bool is_chain(const net_handle_t& net) const;
    virtual bool is_node(const net_handle_t& net) const;
    virtual bool is_sentinel(const net_handle_t& net) const;

    virtual net_handle_t get_net(const handle_t& handle, const HandleGraph* graph) const;

    /// Get the handle for a node, for a trivial chain's node, or for a snarl
    /// bound sentinel. A sentinel gives its bounding node read the way the
    /// sentinel faces: the start going in or out, or the end going in or
    /// out.
    virtual handle_t get_handle(const net_handle_t& net, const HandleGraph* graph) const;

    virtual net_handle_t get_parent(const net_handle_t& child) const;
    virtual net_handle_t get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const;
    virtual net_handle_t flip(const net_handle_t& net) const;
    virtual net_handle_t canonical(const net_handle_t& net) const;
    virtual endpoint_t starts_at(const net_handle_t& traversal) const;
    virtual endpoint_t ends_at(const net_handle_t& traversal) const;

    /// Get the traversal of a snarl between two of its bound sentinels, or of
    /// a chain between two of its bounding node traversals. Throws if either
    /// is an internal tip.
    virtual net_handle_t get_parent_traversal(const net_handle_t& traversal_start,
                                              const net_handle_t& traversal_end) const;

protected:
    virtual bool for_each_child_impl(const net_handle_t& traversal,
                                     const function<bool(const net_handle_t&)>& iteratee) const;
    virtual bool for_each_traversal_impl(const net_handle_t& item,
                                         const function<bool(const net_handle_t&)>& iteratee) const;
    /// Whichever way we go, what we reach is read the way we travel.
    virtual bool follow_net_edges_impl(const net_handle_t& here, const HandleGraph* graph, bool go_left,
                                       const function<bool(const net_handle_t&)>& iteratee) const;

    /// Throws, because the manager does not record internal tips.
    virtual bool for_each_tippy_child_impl(const net_handle_t& parent,
                                           const function<bool(const net_handle_t&)>& iteratee) const;

private:

    /// What a net handle refers to. The payload is a snarl number for a
    /// snarl or a sentinel, the number of its first snarl for a chain, and a
    /// node ID for a node or a trivial chain.
    enum class NetType : uint64_t {
        ROOT = 0,
        SNARL,
        CHAIN,
        TRIVIAL_CHAIN,
        NODE,
        SENTINEL
    };

    /// Where a node lies in the chain whose snarls it bounds.
    struct ChainPlace {
        const Chain* chain;
        /// Index of the node among the chain's bounding nodes, in chain order.
        size_t index;
        /// True if the located handle reads the same way as the chain.
        bool forward;
    };

    static net_handle_t pack(NetType type, endpoint_t start, endpoint_t end, uint64_t payload);
    static NetType type_of(const net_handle_t& net);
    static uint64_t payload_of(const net_handle_t& net);
    /// Make a net handle to the same thing with a different traversal.
    static net_handle_t with_traversal(const net_handle_t& net, endpoint_t start, endpoint_t end);

    /// Get the handle for a Visit to a node.
    handle_t to_handle(const Visit& visit) const;
    /// Get a node's net handle, read the way the handle reads.
    net_handle_t node_net(const handle_t& handle) const;
    /// Get a snarl's net handle, start-to-end if forward and end-to-start
    /// otherwise.
    net_handle_t snarl_net(const Snarl* snarl, bool forward) const;
    /// Get a chain's net handle with the given traversal.
    net_handle_t chain_net(const Chain& chain, endpoint_t start, endpoint_t end) const;
    /// Get the chain a chain net handle refers to.
    const Chain& chain_of_net(const net_handle_t& net) const;

    /// Return true if the decomposition hides the given snarl.
    bool is_hidden(const Snarl* snarl) const;
    /// Return true if the chain's last snarl leads back into its first node.
    bool is_cyclic(const Chain& chain) const;
    /// Count the distinct bounding nodes of a chain.
    size_t node_count(const Chain& chain) const;
    /// Get a chain's bounding node at the given index, read the way the chain
    /// reads. Index chain.size() is the last node, which for a cyclic chain is
    /// the first node again.
    handle_t chain_node(const Chain& chain, size_t index) const;
    /// List a chain's children in chain order, read the way the chain reads.
    vector<net_handle_t> chain_children(const Chain& chain) const;
    /// Find where a node lies in the chain whose snarls it bounds. Returns
    /// false if the node bounds no snarl.
    bool locate(const handle_t& handle, ChainPlace& place) const;

    /// Get the net handle for a parent snarl, or the root if given null.
    net_handle_t parent_snarl_net(const Snarl* parent) const;
    /// Get the nodes in a snarl that bound no snarl, in ID order.
    vector<nid_t> trivial_chain_nodes(const Snarl* snarl) const;
    /// Build the index of trivial chains to their parents, once.
    void index_trivial_chains() const;
    /// Get the parent snarl of the trivial chain of the given node, or null
    /// if its parent is the root.
    const Snarl* trivial_chain_parent(nid_t id) const;

    /// Call the iteratee on what we reach by leaving here through the end of
    /// its traversal.
    bool follow_right(const net_handle_t& here, const function<bool(const net_handle_t&)>& iteratee) const;

    const SnarlManager& manager;
    const HandleGraph& graph;

    /// The manager's snarls, by snarl number.
    vector<const Snarl*> snarls_by_number;
    /// Whether each snarl, by snarl number, is hidden.
    vector<bool> hidden;

    /// Guards the lazily built trivial chain index.
    mutable once_flag trivial_chains_indexed;
    /// For each node of a trivial chain, its parent snarl, or null if its
    /// parent is the root.
    mutable unordered_map<nid_t, const Snarl*> trivial_chain_parents;
    /// The nodes of trivial chains in the root, in ID order.
    mutable vector<nid_t> root_trivial_chain_nodes;
};

}

#endif
