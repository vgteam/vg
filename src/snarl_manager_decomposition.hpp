#ifndef VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED
#define VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED

/** \file
 * The handlegraph::SnarlDecomposition interface over the snarls that a
 * SnarlManager holds, for code written against that interface, such as code
 * that also runs on a SnarlDistanceIndex.
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
 * @brief A view of a SnarlManager's snarls through the SnarlDecomposition
 * interface.
 *
 * @details The view shows each snarl and chain the manager holds, in these
 * terms:
 *
 * - A chain alternates *chain nodes* and snarls, starting and ending with a
 *   chain node. Each chain node has a *position*, counted from 0 at the
 *   chain's start bound. The view orients each snarl forward along its chain:
 *   a snarl's start bound is the chain node before it, and its end bound is
 *   the chain node after it. Where the manager's Snarl is backward in its
 *   chain, its start() is the view's end bound.
 * - A *cyclic chain*'s last snarl leads back into its first chain node, in
 *   the same orientation. That node is both of the chain's bounds, and is one
 *   child of the chain.
 * - Each node of the graph is either a chain node of one chain, or a *free
 *   node*. A free node is directly in a snarl, or in the root, as a trivial
 *   chain of its own. The trivial chain's start-to-end traversal is the node
 *   oriented forward.
 * - A *dead end* is a side of a node with no edges. A *tip* of a snarl or
 *   chain is a dead end on a node inside it, other than its bound nodes. A
 *   traversal that starts or ends at a tip has TIP as that endpoint.
 *
 * Methods that take a graph pointer use that graph to convert between
 * handles and node IDs, so it must have the same node IDs as the graph the
 * view was made with.
 *
 * @note A SnarlDistanceIndex built from the same snarl finder presents the
 * same snarl tree, except that:
 * - The view's root holds each root chain directly. A distance index groups
 *   adjacent root chains under root snarls, for which is_root() is true.
 * - The view shows snarls that have no children. A distance index shows the
 *   two bounds of such a snarl as neighbouring chain nodes.
 */
class SnarlManagerDecomposition : public handlegraph::SnarlDecomposition {
public:
    /**
     * @brief Make a view of the snarls in a SnarlManager.
     *
     * @details Finds which traversals of each snarl and chain are
     * realizable, in time linear in the size of the graph.
     *
     * @param manager The manager whose snarls to present. The view keeps a
     * reference to it, so it must outlive the view and must not change while
     * the view exists.
     * @param graph The graph the snarls were found in. The view keeps a
     * reference to it, so it must outlive the view.
     *
     * @pre finish() has been called on the manager.
     * @throws std::invalid_argument If the manager holds a unary snarl, whose
     * end is its start reversed.
     */
    SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph);

    ~SnarlManagerDecomposition() override = default;

    /**
     * @brief Get the start-to-end traversal of a snarl.
     *
     * @param snarl A snarl owned by the manager the view presents.
     * @return The traversal from the snarl's start bound in the view to its
     * end bound. When the Snarl is backward in its chain, this traversal goes
     * from the Snarl's end() to its start().
     */
    net_handle_t get_snarl_net(const Snarl* snarl) const;

    /**
     * @brief Get the manager's Snarl for a traversal of a snarl.
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
     * @brief Get the graph handle for a node, a trivial chain, or a sentinel.
     *
     * @details A trivial chain gives its node, oriented the way the traversal
     * goes. A sentinel gives its bound node, oriented into the snarl if the
     * sentinel faces in, and out of the snarl if the sentinel faces out.
     *
     * @throws std::invalid_argument If the net handle is for a snarl, the
     * root, or a chain that holds a snarl.
     */
    handle_t get_handle(const net_handle_t& net, const HandleGraph* handle_graph) const override;

    net_handle_t get_parent(const net_handle_t& child) const override;

    /// Give a cyclic chain's first chain node as both its start bound and
    /// its end bound. Facing into the chain at the start, or out of it at the
    /// end, that node is oriented forward along the chain.
    net_handle_t get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const override;

    /// Flip a sentinel to the same bound facing the other way.
    net_handle_t flip(const net_handle_t& net) const override;

    /// Pick as the canonical traversal of a snarl, chain or node the first
    /// realizable one of start-to-end, start-to-tip, end-to-tip,
    /// start-to-start, end-to-end and tip-to-tip. Pick as the canonical
    /// traversal of a sentinel its bound facing in.
    /// TODO: An element with no realizable traversal gets start-to-end, but
    /// the interface requires the canonical traversal to be realizable.
    net_handle_t canonical(const net_handle_t& net) const override;

    endpoint_t starts_at(const net_handle_t& traversal) const override;
    endpoint_t ends_at(const net_handle_t& traversal) const override;

    /// Throw std::invalid_argument if the two arguments are in different
    /// snarls or chains, or if their snarl or chain has no realizable
    /// traversal from the first to the second.
    net_handle_t get_parent_traversal(const net_handle_t& traversal_start,
                                      const net_handle_t& traversal_end) const override;

protected:
    /// Give the children of a snarl, or of the root, in this order: its
    /// chains, in the manager's order, then the trivial chains of its free
    /// nodes, by node ID. Give the children of a chain, for a traversal other
    /// than start-to-end and end-to-start, in order from the chain's start.
    bool for_each_child_impl(const net_handle_t& traversal,
                             const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /// Give the realizable traversals in this order, skipping the others:
    /// start-to-start, start-to-end, start-to-tip, end-to-start, end-to-end,
    /// end-to-tip, tip-to-start, tip-to-end, tip-to-tip. Give a sentinel as
    /// its own only traversal.
    bool for_each_traversal_impl(const net_handle_t& item,
                                 const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /// Going right, give each realizable traversal that a walk leaving
    /// through this traversal's end can visit next, starting where the walk
    /// enters it. Going left, give each traversal from which going right
    /// gives this one.
    /// In a cyclic chain, going right from the last snarl gives the first
    /// chain node, and going right from that node gives the first snarl.
    /// Walk the graph the view was made with, and ignore handle_graph.
    bool follow_net_edges_impl(const net_handle_t& here, const HandleGraph* handle_graph, bool go_left,
                               const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /// For a snarl or the root, give each traversal of a child chain that
    /// starts at a tip of that chain, and each traversal of a child chain or
    /// trivial chain that starts at a bound whose outer side is a dead end.
    /// For a chain, give each traversal of a child snarl that starts at a tip
    /// of that snarl. Give them in an arbitrary order.
    bool for_each_tippy_child_impl(const net_handle_t& parent,
                                   const std::function<bool(const net_handle_t&)>& iteratee) const override;

private:

    /// The place of a chain node in its chain.
    struct ChainPosition {
        /// The chain, owned by the manager.
        const Chain* chain;
        /// The node's position, from 0 at the chain's start bound.
        /// Positions are less than chain_node_count(*chain).
        size_t position;
        /// True if the node's reverse orientation is forward along the chain.
        bool is_reverse;
    };

    /// A directory that gives each node's chain position, if it is a chain
    /// node, and the snarl that directly contains it.
    class NodeDirectory {
    public:
        /**
         * @brief Find the place of each node of a graph.
         *
         * @param manager The manager whose chains and snarls hold the nodes.
         * The directory keeps pointers into it, so it must outlive the
         * directory.
         * @param graph The graph whose nodes to place.
         *
         * @pre finish() has been called on the manager.
         */
        NodeDirectory(const SnarlManager& manager, const HandleGraph& graph);

        /**
         * @brief Get a node's place in its chain.
         *
         * @return The place, or nullopt for a free node.
         *
         * @pre The node is in the graph.
         */
        std::optional<ChainPosition> chain_position(nid_t id) const;

        /**
         * @brief Get the snarl that directly contains a node.
         *
         * @return The snarl that contains the node's chain or trivial chain,
         * or null for the root.
         *
         * @pre The node is in the graph.
         */
        const Snarl* containing_snarl(nid_t id) const;

        /// Get the free nodes in a snarl, or in the root if given null, by
        /// ID.
        const std::vector<nid_t>& free_nodes(const Snarl* snarl) const;

    private:
        /// One node's place in the decomposition.
        struct Entry {
            /// The snarl that directly contains the node, or null for the
            /// root.
            const Snarl* containing_snarl;
            /// The node's place in its chain, or nullopt for a free node.
            std::optional<ChainPosition> chain_position;
        };

        /// The entry for each node of the graph, by ID.
        std::unordered_map<nid_t, Entry> entries;
        /// The free nodes in each snarl that has any, and in the root under
        /// null. Each list is in ID order.
        std::unordered_map<const Snarl*, std::vector<nid_t>> free_nodes_by_snarl;
    };

    /// A set of traversal kinds, each named by the endpoints where it starts
    /// and ends.
    class TraversalSet {
    public:
        /// Add the traversal kind from one endpoint to another.
        void add(endpoint_t start, endpoint_t end);
        /// Return true if the set holds the traversal kind.
        bool contains(endpoint_t start, endpoint_t end) const;
        /**
         * @brief Call a function on each traversal kind in the set, until it
         * returns false.
         *
         * @param iteratee A function that takes the start and end endpoints,
         * and returns false to stop.
         * @return False if the function stopped the loop, and true otherwise.
         *
         * @post Kinds come ordered by start endpoint, then by end endpoint,
         * each in the order START, END, TIP.
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
        /// One bit for each traversal kind, at index 3 * start + end, where
        /// START, END and TIP are 0, 1 and 2.
        uint16_t bits = 0;
    };

    /**
     * @brief A survey of the realizable traversals of each snarl and chain.
     *
     * @details The constructor finds them by walking the graph. A walk
     * through a snarl crosses each child chain, turns around in it, or stops
     * at a tip in it, as the child chain's realizable traversals allow. A
     * walk along a chain goes through each of its snarls as that snarl's
     * realizable traversals allow. So the survey needs the traversals of
     * each child before those of its parent.
     */
    class TraversalSurvey {
    public:
        /**
         * @brief Survey each snarl and chain in a manager.
         *
         * @param manager The manager whose snarls and chains to survey.
         * @param graph The graph the snarls were found in.
         * @param directory The places of the graph's nodes.
         */
        TraversalSurvey(const SnarlManager& manager, const HandleGraph& graph, const NodeDirectory& directory);

        /// Get the realizable traversals of a snarl, oriented forward along its
        /// chain.
        const TraversalSet& of_snarl(const Snarl* snarl) const;
        /// Get the realizable traversals of a chain.
        const TraversalSet& of_chain(const Chain* chain) const;

    private:
        /**
         * @brief A walker that finds where walks through a snarl or along a
         * chain can end, from the traversals already surveyed for the
         * snarl's or chain's children.
         *
         * @details The survey's constructor makes one, and asks it about
         * each snarl and chain, innermost first.
         */
        class Walker {
        public:
            /**
             * @brief Make a walker for one survey.
             *
             * @details The walker keeps references to its arguments, which
             * must outlive it.
             */
            Walker(const SnarlManager& manager, const HandleGraph& graph, const NodeDirectory& directory,
                   const TraversalSurvey& survey);

            /**
             * @brief Find the realizable traversals of a snarl.
             *
             * @pre Each chain in the snarl has been surveyed.
             */
            TraversalSet survey_snarl(const Snarl* snarl) const;
            /**
             * @brief Find the realizable traversals of a chain.
             *
             * @pre Each snarl in the chain has been surveyed.
             */
            TraversalSet survey_chain(const Chain& chain) const;

        private:
            /// The endpoints of a snarl or chain where a walk can end.
            struct WalkEnds {
                /// True if the walk can leave through the start bound.
                bool start = false;
                /// True if the walk can leave through the end bound.
                bool end = false;
                /// True if the walk can end at a tip.
                bool tip = false;
            };

            /**
             * @brief Walk through a snarl from some handles, and find where the
             * walk can end.
             *
             * @param starts Handles in the snarl. The walk starts on each one,
             * and leaves it through its right side. A start may be a bound of
             * the snarl, oriented into the snarl.
             * @return Whether the walk can leave the snarl through its start
             * bound, whether it can leave through its end bound, and whether it
             * can end at a tip in the snarl.
             *
             * @pre Each chain in the snarl has been surveyed.
             */
            WalkEnds walk_snarl(const Snarl* snarl, const std::vector<handle_t>& starts) const;

            /**
             * @brief Walk along a chain from some chain nodes, and find where the
             * walk can end.
             *
             * @param starts Places to start: each a chain node's position, and
             * true to leave it going in reverse along the chain, or false to
             * leave it going forward along the chain.
             * @return Whether the walk can leave the chain through its start
             * bound, whether it can leave through its end bound, and whether it
             * can end at a tip in one of the chain's snarls.
             *
             * @pre Each snarl in the chain has been surveyed.
             */
            WalkEnds walk_chain(const Chain& chain, const std::vector<std::pair<size_t, bool>>& starts) const;

            /**
             * @brief Find where a walk that enters a child chain or free node
             * can leave it.
             *
             * @param entering A chain's bound node oriented into the chain, or a
             * handle on a free node.
             * @param exits A list to add the child's exits to. An exit is a
             * handle that a walk leaves the child by, through the handle's right
             * side. A chain's exits are its bound nodes, oriented out of the
             * chain, at each bound where a walk entering the chain can leave it.
             * A free node's one exit is the entering handle.
             * @return True if a walk entering the chain can end at a tip in it.
             *
             * @pre The chain has been surveyed.
             * @throws std::logic_error If the handle enters a chain between
             * its bounds, which a chain does not allow, so the decomposition
             * is broken.
             */
            bool enter_child(const handle_t& entering, std::vector<handle_t>& exits) const;

            /// The manager whose snarls and chains are surveyed.
            const SnarlManager& manager;
            /// The graph the snarls were found in.
            const HandleGraph& graph;
            /// Tells the walker which chain, if any, each node is a chain node
            /// of.
            const NodeDirectory& directory;
            /// The survey being made, which holds the traversals of each child
            /// already surveyed.
            const TraversalSurvey& survey;
        };

        /// The realizable traversals of each snarl, oriented forward along its
        /// chain.
        std::unordered_map<const Snarl*, TraversalSet> snarls;
        /// The realizable traversals of each chain.
        std::unordered_map<const Chain*, TraversalSet> chains;
    };

    /// Get a handle, in the view's graph, on a snarl's start bound oriented
    /// into the snarl, or on its end bound oriented out of the snarl.
    handle_t snarl_bound(const Snarl* snarl, bool get_end) const;
    /// Get a snarl's net handle, start-to-end if forward and end-to-start
    /// otherwise.
    net_handle_t snarl_net(const Snarl* snarl, bool forward) const;
    /// Get a chain's start-to-end net handle.
    net_handle_t chain_net(const Chain& chain) const;
    /// Get the chain that a chain net handle is for.
    const Chain& chain_of_net(const net_handle_t& net) const;
    /// Get the net handle for the chain node at a position, oriented in
    /// reverse along the chain if is_reverse, and forward along it otherwise.
    net_handle_t chain_node_net(const Chain& chain, size_t position, bool is_reverse) const;

    /**
     * @brief Call a function on each child of a snarl, or of the root if
     * given null.
     *
     * @details Gives each child start-to-end: first the snarl's chains, in
     * the manager's order, then the trivial chains of its free nodes, by node
     * ID.
     *
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_snarl_child(const Snarl* parent, const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * @brief Call a function on each chain node and snarl of a chain.
     *
     * @details For a start-to-end or end-to-start traversal, children come in
     * the order the traversal visits them, oriented the way it visits them. For
     * other traversals, they come in order from the chain's start, each
     * start-to-end.
     *
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_chain_child(const net_handle_t& traversal,
                              const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * @brief Get the next child in a chain after a chain node or snarl.
     *
     * @param here A start-to-end or end-to-start traversal of a node, or any
     * traversal of a snarl.
     * @return
     * - The chain node or snarl that a walk leaving through the end of here
     *   enters, oriented the way the walk goes.
     * - nullopt, when here ends at a tip, leaves its chain at a bound, or is
     *   a free node.
     */
    std::optional<net_handle_t> next_in_chain(const net_handle_t& here) const;

    /**
     * @brief Get the child of a snarl or the root that a handle enters.
     *
     * @param handle A handle at the far end of an edge from one of the
     * parent's children or from one of its bounds.
     * @param parent The snarl, or null for the root.
     * @return
     * - The traversal of the child chain or trivial chain that the handle
     *   enters at one of its bounds.
     * - The parent's bound sentinel facing out, if the handle leaves the
     *   parent.
     *
     * @throws std::logic_error If the handle enters a chain between its
     * bounds, which a chain does not allow, so the decomposition is broken.
     */
    net_handle_t sibling_entered_by(const handle_t& handle, const Snarl* parent) const;

    /**
     * @brief Call a function on each sibling that a walk leaving a traversal
     * through its end enters.
     *
     * @details A sibling is a chain or trivial chain in the same parent,
     * given as its start-to-end or end-to-start traversal, or one of the
     * parent's sentinels, facing out. A sentinel facing out, a cyclic chain,
     * and a traversal that ends at a tip have no siblings ahead.
     *
     * @param here A traversal of a chain or trivial chain, or a sentinel.
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_sibling_ahead(const net_handle_t& here,
                                const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * @brief Call a function on each realizable traversal of an element that
     * starts where a walk enters it.
     *
     * @param entered A traversal of the element that starts where the walk
     * enters it.
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_traversal_entered(const net_handle_t& entered,
                                    const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * @brief Call a function on each realizable traversal that a walk leaving
     * a traversal through its end can visit next, starting where the walk
     * enters it.
     *
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool follow_right(const net_handle_t& here, const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * @brief Find the snarl or chain whose traversal an element can start or
     * end, and the endpoint of that snarl or chain.
     *
     * @param net One of:
     * - A sentinel facing into its snarl, if is_start, or out of it
     *   otherwise. This gives the snarl and the sentinel's bound.
     * - A node traversal that enters its chain or trivial chain at a bound,
     *   if is_start, or leaves it at a bound otherwise. This gives the chain
     *   and that bound.
     * - A snarl, chain or trivial chain, which gives its parent and TIP.
     * @param is_start True to find where a traversal starts, and false to
     * find where it ends.
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
    /// The graph the snarls were found in.
    const HandleGraph& graph;
    /// The manager's snarls, by SnarlManager::snarl_number().
    std::vector<const Snarl*> snarls;
    /// Tells the view each node's chain position, and the snarl that
    /// contains it.
    NodeDirectory directory;
    /// Tells the view which traversals of each snarl and chain are
    /// realizable.
    TraversalSurvey survey;
};

}

#endif
