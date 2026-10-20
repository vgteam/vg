#ifndef VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED
#define VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED

/** \file
 * The handlegraph::SnarlDecomposition interface over the snarls that a
 * SnarlManager holds.
 */

#include <handlegraph/snarl_decomposition.hpp>

#include <functional>
#include <optional>
#include <unordered_map>
#include <vector>

#include "handle.hpp"
#include "snarls.hpp"

namespace vg {

using handlegraph::net_handle_t;

/**
 * A view of the snarls and chains in a SnarlManager, through the
 * SnarlDecomposition interface.
 *
 * A chain alternates chain nodes and snarls, starting and ending with a chain
 * node. Each chain node has a position, counting from 0 at the start bound of
 * the chain. The snarl of rank i lies between the chain nodes at positions i
 * and i + 1.
 *
 * Each snarl is oriented forward along its chain. So a Snarl that is reversed
 * in parent starts at its end() flipped and ends at its start() flipped.
 *
 * In a looping chain, the last snarl leads back into the first chain node.
 * That node is both bounds of the chain, and one child at position 0.
 *
 * Each node that is not a chain node is a free node. It is a trivial chain
 * whose start-to-end traversal visits the node forward.
 *
 * The view keeps a reference to the manager, which must outlive it and must
 * not change. It does not keep the graph. A graph passed to a method must have
 * the same nodes and edges as the graph the view was made with.
 *
 * @note Unlike a SnarlDistanceIndex, the view puts top-level chains directly
 * in the root, and shows snarls that have no children.
 */
class SnarlManagerDecomposition : public handlegraph::SnarlDecomposition {
public:
    /**
     * Make a view of the snarls in a SnarlManager.
     *
     * Walk the graph to find which traversals of each snarl and chain are
     * realizable.
     *
     * @pre finish() has been called on the manager.
     * @throws std::invalid_argument If the manager holds a unary snarl.
     */
    SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph);

    ~SnarlManagerDecomposition() override = default;

    /// Get the start-to-end traversal of a snarl.
    net_handle_t get_snarl_net(const Snarl* snarl) const;

    /**
     * Get the manager's Snarl for a traversal of a snarl.
     *
     * @throws std::invalid_argument If the net handle is not for a snarl.
     */
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

    /**
     * Get the graph handle for a node, a trivial chain, or a sentinel.
     *
     * A sentinel gives its boundary node, oriented the way the sentinel
     * faces.
     *
     * @throws std::invalid_argument If the net handle is for a snarl, the
     * root, or a chain that holds a snarl.
     */
    handle_t get_handle(const net_handle_t& net, const HandleGraph* handle_graph) const override;

    net_handle_t get_parent(const net_handle_t& child) const override;

    /**
     * Get a sentinel of a snarl, or a boundary node of a chain.
     *
     * A looping chain gives its first chain node as both bounds.
     */
    net_handle_t get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const override;

    /// Flip a traversal, or turn a sentinel to face the other way.
    net_handle_t flip(const net_handle_t& net) const override;

    /**
     * Get the canonical traversal of a net graph element.
     *
     * For a snarl, chain or node, pick the first realizable traversal among
     * start-to-end, start-to-tip, end-to-tip, start-to-start, end-to-end and
     * tip-to-tip. For a sentinel, pick its bound facing in.
     *
     * TODO: An element with no realizable traversal gets start-to-end, but
     * the interface requires the canonical traversal to be realizable. Such
     * elements exist: IntegratedSnarlFinder makes a snarl around a cycle that
     * walks from both of its bounds fall into and never leave.
     */
    net_handle_t canonical(const net_handle_t& net) const override;

    endpoint_t starts_at(const net_handle_t& traversal) const override;
    endpoint_t ends_at(const net_handle_t& traversal) const override;

    /**
     * Get the traversal of a snarl or chain from one of its children or bounds
     * to another.
     *
     * @throws std::invalid_argument If the two arguments are in different
     * snarls or chains, or if the traversal between them is not realizable.
     */
    net_handle_t get_parent_traversal(const net_handle_t& traversal_start,
                                      const net_handle_t& traversal_end) const override;

protected:
    /**
     * Call a function on each child of the root, a snarl, a chain or a
     * trivial chain.
     *
     * The root or a snarl gives its chains in the manager's order, then its
     * free nodes in ID order. A chain traversal other than start-to-end and
     * end-to-start gives the children in chain order, each start-to-end.
     */
    bool for_each_child_impl(const net_handle_t& traversal,
                             const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * Call a function on each realizable traversal of a net graph element.
     *
     * Give them ordered by start endpoint, then end endpoint, each in the
     * order START, END, TIP. A sentinel is its own only traversal.
     */
    bool for_each_traversal_impl(const net_handle_t& item,
                                 const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * Call a function on each realizable traversal that comes just after or
     * just before a traversal in its parent.
     *
     * In a looping chain, going right from the last snarl gives the first
     * chain node.
     */
    bool follow_net_edges_impl(const net_handle_t& here, const HandleGraph* handle_graph, bool go_left,
                               const std::function<bool(const net_handle_t&)>& iteratee) const override;

private:

    /// The place of a chain node in its chain.
    struct ChainPosition {
        /// The chain, owned by the manager.
        const Chain* chain;
        /// The position of the node in the chain. The boundary node of a
        /// looping chain is at position 0.
        size_t position;
        /// True if the node is reversed in parent.
        bool is_reverse;
    };

    /**
     * A directory of where each node of a graph sits in the snarl tree.
     *
     * The directory finds chain nodes through the manager's index of the
     * snarls each node side leads into, and lists the free nodes itself.
     */
    class NodeDirectory {
    public:
        /**
         * Find the free nodes of a graph in a manager's snarl tree.
         *
         * The directory keeps a reference to the manager, so it must outlive
         * the directory.
         *
         * @pre finish() has been called on the manager.
         */
        NodeDirectory(const SnarlManager& manager, const HandleGraph& graph);

        /// Get the place of a node in its chain, or nullopt for a free node,
        /// in constant time.
        std::optional<ChainPosition> chain_position(nid_t id) const;

        /**
         * Get the snarl that directly contains a node, or null for the root.
         *
         * @pre The node is in the graph.
         */
        const Snarl* containing_snarl(nid_t id) const;

        /// Get the free nodes in a snarl, or in the root if given null, in ID
        /// order.
        const std::vector<nid_t>& free_nodes(const Snarl* snarl) const;

    private:
        /// The manager whose snarl tree the nodes sit in.
        const SnarlManager& manager;
        /// The snarl that directly contains each free node, or null for the
        /// root.
        std::unordered_map<nid_t, const Snarl*> free_node_snarls;
        /// The free nodes in each snarl that has any, and in the root under
        /// null. Each list is in ID order.
        std::unordered_map<const Snarl*, std::vector<nid_t>> free_nodes_by_snarl;
    };

    /// A set of traversals of one net graph element, each named by the
    /// endpoints where it starts and ends.
    class TraversalSet {
    public:
        /// Add the traversal from one endpoint to another.
        void add(endpoint_t start, endpoint_t end);
        /// Return true if the set holds the traversal from one endpoint to
        /// another.
        bool contains(endpoint_t start, endpoint_t end) const;
        /**
         * Call a function on each traversal in the set, until it returns
         * false.
         *
         * Give the traversals ordered by start endpoint, then by end
         * endpoint, each in the order START, END, TIP.
         */
        template<typename Iteratee>
        bool for_each(const Iteratee& iteratee) const {
            for (endpoint_t start : {START, END, TIP}) {
                for (endpoint_t end : {START, END, TIP}) {
                    if (contains(start, end) && !iteratee(start, end)) {
                        return false;
                    }
                }
            }
            return true;
        }

    private:
        /// One bit for each traversal, at index 3 * start + end, where START,
        /// END and TIP are 0, 1 and 2.
        uint16_t bits = 0;
    };

    /**
     * A survey of the realizable traversals of each snarl and chain in a
     * manager.
     *
     * The constructor does the survey, and nothing changes it afterwards.
     */
    class TraversalSurvey {
    public:
        /// Survey each snarl and chain in a manager.
        TraversalSurvey(const SnarlManager& manager, const HandleGraph& graph, const NodeDirectory& directory);

        /// Get the realizable traversals of a snarl.
        const TraversalSet& of_snarl(const Snarl* snarl) const;
        /// Get the realizable traversals of a chain.
        const TraversalSet& of_chain(const Chain* chain) const;

    private:
        /**
         * A surveyor that finds the realizable traversals of one snarl or
         * chain at a time, by walking the graph.
         *
         * A walk crosses each child as the survey says the traversals of that
         * child allow. So the surveyor must survey each child before its
         * parent.
         */
        class Surveyor {
        public:
            /**
             * Make a surveyor for one survey.
             *
             * The surveyor keeps references to its arguments, which must
             * outlive it.
             */
            Surveyor(const SnarlManager& manager, const HandleGraph& graph, const NodeDirectory& directory,
                     const TraversalSurvey& survey);

            /**
             * Find the realizable traversals of a snarl.
             *
             * @pre Each chain in the snarl has been surveyed.
             */
            TraversalSet survey_snarl(const Snarl* snarl) const;
            /**
             * Find the realizable traversals of a chain.
             *
             * @pre Each snarl in the chain has been surveyed.
             */
            TraversalSet survey_chain(const Chain& chain) const;

        private:
            /// The places in a snarl or chain where a walk can end.
            struct WalkEnds {
                /// True if a walk can leave through the start bound.
                bool start = false;
                /// True if a walk can leave through the end bound.
                bool end = false;
                /// True if a walk can end at a tip.
                bool tip = false;
            };

            /**
             * Find where walks through a snarl from some handles can end.
             *
             * A walk leaves each start through its right side. Each start is
             * on a node in the snarl, or on a bound oriented into the snarl.
             *
             * @pre Each chain in the snarl has been surveyed.
             */
            WalkEnds walk_snarl(const Snarl* snarl, const std::vector<handle_t>& starts) const;

            /**
             * Find where walks along a chain from some chain nodes can end.
             *
             * Each start is the position of a chain node, and true to leave it
             * going in reverse along the chain. Here the end bound has position
             * chain.size() even in a looping chain. A walk that comes to a bound
             * of a looping chain can leave there, and can also go on around the
             * loop, since its boundary node is both bounds.
             *
             * @pre Each snarl in the chain has been surveyed.
             */
            WalkEnds walk_chain(const Chain& chain, const std::vector<std::pair<size_t, bool>>& starts) const;

            /**
             * Find where a walk that enters a child chain or free node can
             * leave it.
             *
             * Add each exit to exits. A walk leaves the child through the
             * right side of an exit. The one exit of a free node is the
             * entering handle.
             *
             * @return True if a walk entering the chain can end at a tip in
             * it.
             *
             * @pre The chain has been surveyed.
             * @throws std::logic_error If the handle enters a chain between its
             * bounds.
             */
            bool enter_child(const handle_t& entering, std::vector<handle_t>& exits) const;

            /**
             * Add the traversals made by walks from the start and end bounds.
             *
             * Also add each of those traversals taken backward that starts at
             * a tip.
             */
            static void add_walks_from_bounds(TraversalSet& traversals, const WalkEnds& from_start,
                                              const WalkEnds& from_end);

            /// The manager whose snarls and chains are surveyed.
            const SnarlManager& manager;
            /// The graph the snarls were found in.
            const HandleGraph& graph;
            /// Tells the surveyor which chain, if any, each node is a chain
            /// node of.
            const NodeDirectory& directory;
            /// Tells the surveyor the traversals of each child already
            /// surveyed.
            const TraversalSurvey& survey;
        };

        /// The realizable traversals of each snarl.
        std::unordered_map<const Snarl*, TraversalSet> snarls;
        /// The realizable traversals of each chain.
        std::unordered_map<const Chain*, TraversalSet> chains;
    };

    /// Get a snarl's start bound oriented into the snarl, or its end bound
    /// oriented out of it.
    Visit snarl_bound(const Snarl* snarl, bool get_end) const;
    /// Get the manager's Snarl for a snarl number.
    const Snarl* snarl_of_number(uint64_t number) const;
    /// Get a snarl's net handle, start-to-end if forward and end-to-start
    /// otherwise.
    net_handle_t snarl_net(const Snarl* snarl, bool forward) const;
    /// Get a chain's start-to-end net handle.
    net_handle_t chain_net(const Chain& chain) const;
    /// Get the chain that a chain's net handle is for.
    const Chain& chain_of_net(const net_handle_t& net) const;
    /// Get the net handle for the chain node at a position, oriented in reverse
    /// along the chain if is_reverse, and forward along it otherwise.
    net_handle_t chain_node_net(const Chain& chain, size_t position, bool is_reverse) const;

    /// Call a function on each child of a snarl, or of the root if given null.
    bool for_each_snarl_child(const Snarl* parent, const std::function<bool(const net_handle_t&)>& iteratee) const;

    /// Call a function on each chain node and snarl of a chain.
    bool for_each_chain_child(const net_handle_t& traversal,
                              const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Get the next child in a chain after a chain node or snarl.
     *
     * Return nullopt when here ends at a tip, leaves its chain at a bound, or
     * is a free node.
     *
     * @pre here is a start-to-end or end-to-start traversal of a node, or any
     * traversal of a snarl.
     */
    std::optional<net_handle_t> next_in_chain(const net_handle_t& here) const;

    /**
     * Get the child of a snarl or the root that a handle enters.
     *
     * If the handle leaves the parent, return the sentinel of the parent for
     * that bound, facing out.
     *
     * @throws std::logic_error If the handle enters a chain between its
     * bounds.
     */
    net_handle_t sibling_entered_by(const HandleGraph& graph, const handle_t& handle, const Snarl* parent) const;

    /**
     * Call a function on each sibling that an edge from the end of a traversal
     * leads into.
     *
     * A sibling is a chain or trivial chain in the same parent, or a sentinel
     * of the parent facing out. A looping chain, a sentinel facing out, and a
     * traversal that ends at a tip have no siblings ahead.
     */
    bool for_each_sibling_ahead(const HandleGraph& graph, const net_handle_t& here,
                                const std::function<bool(const net_handle_t&)>& iteratee) const;

    /// Call a function on each realizable traversal of a net graph element
    /// that starts at the given endpoint.
    bool for_each_traversal_from(const net_handle_t& element, endpoint_t start,
                                 const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Do follow_net_edges() going right.
     *
     * Find the net graph elements that edges from the end of here lead into,
     * and give each realizable traversal of them that starts at the endpoint
     * the edge leads into.
     */
    bool follow_right(const HandleGraph& graph, const net_handle_t& here,
                      const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Find the snarl or chain, and its endpoint, where a traversal of it can
     * start or end at a net graph element or sentinel.
     *
     * @param net One of:
     * - A sentinel facing into its snarl, if is_start, or out of it
     *   otherwise. This gives the snarl and the bound of the sentinel.
     * - A node traversal that enters its chain or trivial chain at a bound,
     *   if is_start, or leaves it at a bound otherwise. This gives the chain
     *   and that bound.
     * - A snarl, chain or trivial chain, which gives its parent and TIP.
     * @return The snarl or chain, start-to-end, and the endpoint.
     *
     * @throws std::invalid_argument If net is the root, a child of the root,
     * a sentinel facing the wrong way, or a node traversal that does not
     * enter or leave its chain at a bound.
     */
    std::pair<net_handle_t, endpoint_t> parent_endpoint(const net_handle_t& net, bool is_start) const;

    /// Get the realizable traversals of a net graph element. A sentinel's
    /// only traversal is itself.
    TraversalSet traversals_of(const net_handle_t& item) const;

    /// The manager whose snarls the view presents.
    const SnarlManager& manager;
    /// Tells the view each node's place in its chain, and the snarl that
    /// contains it.
    NodeDirectory directory;
    /// Tells the view which traversals of each snarl and chain are
    /// realizable.
    TraversalSurvey survey;
};

}

#endif
