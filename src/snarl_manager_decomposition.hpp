#ifndef VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED
#define VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED

/** \file
 * The handlegraph::SnarlDecomposition interface over the snarls that a
 * SnarlManager holds.
 *
 * Code written for the interface can run on a SnarlManager as well as on a
 * SnarlDistanceIndex.
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
 * The view describes chains in these terms:
 *
 * - A chain alternates *chain nodes* and snarls, starting and ending with a
 *   chain node. Each chain node has a *position* in its chain, counting from 0
 *   at the chain's start bound. The snarl of rank i in the chain lies between
 *   the chain nodes at positions i and i + 1.
 * - Each snarl is oriented forward along its chain. Its start bound is the
 *   chain node before it, and its end bound is the chain node after it. So
 *   for a Snarl that is reversed in parent, the view's start bound is the
 *   Snarl's end() flipped, and the view's end bound is its start() flipped.
 * - In a looping chain, the last snarl leads back into the first chain node,
 *   in the same orientation. That node is both bounds of the chain, and is
 *   one child of the chain, at position 0.
 * - Each node of the graph is either a chain node of one chain, or a *free
 *   node*. A free node is a trivial chain, directly in a snarl or in the
 *   root. The start-to-end traversal of the trivial chain visits the node in
 *   its forward orientation.
 *
 * Methods that take a graph pointer use that graph to convert between
 * handles and node IDs. So it must have the same node IDs as the graph the
 * view was made with.
 *
 * @note A SnarlDistanceIndex built from the same snarl finder presents the
 * same snarl tree, with two differences:
 * - The root of the view holds each top-level chain directly. A distance
 *   index groups top-level chains under root snarls, for which is_root() is
 *   true.
 * - The view shows snarls that have no children. A distance index shows the
 *   two bounds of such a snarl as neighbouring chain nodes.
 */
class SnarlManagerDecomposition : public handlegraph::SnarlDecomposition {
public:
    /**
     * Make a view of the snarls in a SnarlManager.
     *
     * Walk the graph to find which traversals of each snarl and chain are
     * realizable.
     *
     * @param manager The manager whose snarls to present. The view keeps a
     * reference to it, so it must outlive the view and must not change while
     * the view exists.
     * @param graph The graph the snarls were found in. The view keeps a
     * reference to it, so it must outlive the view.
     *
     * @pre finish() has been called on the manager.
     * @throws std::invalid_argument If the manager holds a unary snarl.
     */
    SnarlManagerDecomposition(const SnarlManager& manager, const HandleGraph& graph);

    ~SnarlManagerDecomposition() override = default;

    /**
     * Get the start-to-end traversal of a snarl.
     *
     * For a Snarl that is reversed in parent, the traversal goes from the
     * Snarl's end() to its start().
     *
     * @param snarl A snarl owned by the manager the view presents.
     */
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
     * A trivial chain gives its node, oriented the way the traversal goes. A
     * sentinel gives its boundary node, oriented into the snarl if the
     * sentinel faces in, and out of the snarl if the sentinel faces out.
     *
     * @throws std::invalid_argument If the net handle is for a snarl, the
     * root, or a chain that holds a snarl.
     */
    handle_t get_handle(const net_handle_t& net, const HandleGraph* handle_graph) const override;

    net_handle_t get_parent(const net_handle_t& child) const override;

    /**
     * Get a sentinel of a snarl, or a boundary node of a chain.
     *
     * A looping chain gives its first chain node as both its start bound and
     * its end bound. That node is oriented forward along the chain when it
     * faces into the chain at the start, or out of it at the end.
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
     * the interface requires the canonical traversal to be realizable.
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
     * The root or a snarl gives its chains in the manager's order, then the
     * trivial chains of its free nodes in node ID order. A chain traversal
     * other than start-to-end and end-to-start gives the children in order
     * from the chain's start, each start-to-end.
     */
    bool for_each_child_impl(const net_handle_t& traversal,
                             const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * Call a function on each realizable traversal of a net graph element.
     *
     * Give them in this order: start-to-start, start-to-end, start-to-tip,
     * end-to-start, end-to-end, end-to-tip, tip-to-start, tip-to-end,
     * tip-to-tip. A sentinel is its own only traversal.
     */
    bool for_each_traversal_impl(const net_handle_t& item,
                                 const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * Call a function on each realizable traversal that a walk can visit
     * just after or just before a traversal.
     *
     * Going right, give each realizable traversal that starts where a walk
     * leaving here through its end enters the next net graph element. Going
     * left, give each traversal from which going right gives here. In a
     * looping chain, going right from the last snarl gives the first chain
     * node, and going right from that node gives the first snarl.
     *
     * Walk the graph the view was made with, and ignore handle_graph.
     */
    bool follow_net_edges_impl(const net_handle_t& here, const HandleGraph* handle_graph, bool go_left,
                               const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * Call a function on each realizable traversal of a child that starts at
     * a tip of the parent.
     *
     * For the root or a snarl, give each traversal of a child chain that
     * starts at a tip of that chain. Also give each traversal of a child
     * chain or trivial chain that starts at a bound whose outer side is a
     * dead end. For a chain, give each traversal of a child snarl that starts
     * at a tip of that snarl. Give them in no particular order.
     */
    bool for_each_tippy_child_impl(const net_handle_t& parent,
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
     * The directory gives each node's place in its chain, if it is a chain
     * node, and the snarl that directly contains the node. It also lists the
     * free nodes of each snarl.
     */
    class NodeDirectory {
    public:
        /**
         * Find where each node of a graph sits in a manager's snarl tree.
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
         * Get the place of a node in its chain.
         *
         * @return The place, or nullopt for a free node.
         *
         * @pre The node is in the graph.
         */
        std::optional<ChainPosition> chain_position(nid_t id) const;

        /**
         * Get the snarl that directly contains a node.
         *
         * @return The snarl that contains the chain or trivial chain of the
         * node, or null for the root.
         *
         * @pre The node is in the graph.
         */
        const Snarl* containing_snarl(nid_t id) const;

        /// Get the free nodes in a snarl, or in the root if given null, in ID
        /// order.
        const std::vector<nid_t>& free_nodes(const Snarl* snarl) const;

    private:
        /// The place of one node in the snarl tree.
        struct Entry {
            /// The snarl that directly contains the node, or null for the
            /// root.
            const Snarl* containing_snarl;
            /// The place of the node in its chain, or nullopt for a free node.
            std::optional<ChainPosition> chain_position;
        };

        /// The entry for each node of the graph, by ID.
        std::unordered_map<nid_t, Entry> entries;
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
         *
         * @param iteratee A function that takes the start and end endpoints,
         * and returns false to stop.
         * @return False if the function stopped the loop, and true otherwise.
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
     * Each snarl is oriented forward along its chain. The constructor does
     * the survey, and nothing changes it afterwards.
     */
    class TraversalSurvey {
    public:
        /**
         * Survey each snarl and chain in a manager.
         *
         * @param manager The manager whose snarls and chains to survey.
         * @param graph The graph the snarls were found in.
         * @param directory The directory of the nodes of the graph.
         */
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
         * In a snarl, a walk follows edges from node to node. It crosses each
         * child chain, turns around in it, or stops at a tip in it, as the
         * survey says the traversals of that chain allow. Along a chain, a
         * walk goes through each snarl as the survey says the traversals of
         * that snarl allow. So the surveyor must survey each child before its
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
             * A walk leaves a start through its right side, and goes on until
             * it leaves the snarl through a bound or comes to a dead end.
             *
             * @param starts Handles on nodes in the snarl, or on bounds of the
             * snarl oriented into it.
             *
             * @pre Each chain in the snarl has been surveyed.
             */
            WalkEnds walk_snarl(const Snarl* snarl, const std::vector<handle_t>& starts) const;

            /**
             * Find where walks along a chain from some chain nodes can end.
             *
             * @param starts The places to start. Each is the position of a chain
             * node, and true to leave that node going in reverse along the
             * chain, or false to leave it going forward along the chain.
             * Here the end bound has position chain.size() even in a looping
             * chain, so that a walk going forward out of the last snarl
             * leaves the chain instead of going around again.
             *
             * @pre Each snarl in the chain has been surveyed.
             */
            WalkEnds walk_chain(const Chain& chain, const std::vector<std::pair<size_t, bool>>& starts) const;

            /**
             * Find where a walk that enters a child chain or free node can
             * leave it.
             *
             * @param entering A boundary node of a chain, oriented into the
             * chain, or a handle on a free node.
             * @param exits A list to add the exits of the child to. A walk
             * leaves the child through the right side of an exit. The exits of
             * a chain are its boundary nodes, oriented out of the chain, at
             * each bound where a walk entering the chain can leave it. The one
             * exit of a free node is the entering handle.
             * @return True if a walk entering the chain can end at a tip in
             * it.
             *
             * @pre The chain has been surveyed.
             * @throws std::logic_error If the handle enters a chain between its
             * bounds. A correct snarl decomposition has no edge that does.
             */
            bool enter_child(const handle_t& entering, std::vector<handle_t>& exits) const;

            /**
             * Add the traversals made by walks from the start and end bounds.
             *
             * Also add each traversal from a tip that is one of those
             * traversals taken backward.
             *
             * @param from_start Where walks entering at the start bound can
             * end.
             * @param from_end Where walks entering at the end bound can end.
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
    /// oriented out of it, as a handle in the view's graph.
    handle_t snarl_bound(const Snarl* snarl, bool get_end) const;
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

    /**
     * Call a function on each child of a snarl, or of the root if given null.
     *
     * Give each child start-to-end. Give the chains of the snarl first, in
     * the manager's order, then the trivial chains of its free nodes, in node
     * ID order.
     *
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_snarl_child(const Snarl* parent, const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Call a function on each chain node and snarl of a chain.
     *
     * For a start-to-end or end-to-start traversal, give the children in the
     * order the traversal visits them, oriented the way it visits them. For
     * other traversals, give them in order from the start of the chain, each
     * start-to-end.
     *
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_chain_child(const net_handle_t& traversal,
                              const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Get the next child in a chain after a chain node or snarl.
     *
     * @param here A start-to-end or end-to-start traversal of a node, or any
     * traversal of a snarl.
     * @return The chain node or snarl that a walk leaving here through its end
     * enters, oriented the way the walk goes. Return nullopt when here ends at
     * a tip, leaves its chain at a bound, or is a free node.
     */
    std::optional<net_handle_t> next_in_chain(const net_handle_t& here) const;

    /**
     * Get the child of a snarl or the root that a handle enters.
     *
     * @param handle A handle at the far end of an edge from a child or bound
     * of the parent.
     * @param parent The snarl, or null for the root.
     * @return The traversal of the child chain or trivial chain that the
     * handle enters at one of its bounds. If the handle leaves the parent,
     * return the parent's sentinel for that bound, facing out.
     *
     * @throws std::logic_error If the handle enters a chain between its
     * bounds. A correct snarl decomposition has no edge that does.
     */
    net_handle_t sibling_entered_by(const handle_t& handle, const Snarl* parent) const;

    /**
     * Call a function on each sibling that a walk leaving a traversal through
     * its end enters.
     *
     * A sibling is a chain or trivial chain in the same parent, given as its
     * start-to-end or end-to-start traversal, or one of the parent's
     * sentinels, facing out. A sentinel facing out, a looping chain, and a
     * traversal that ends at a tip have no siblings ahead.
     *
     * @param here A traversal of a chain or trivial chain, or a sentinel.
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_sibling_ahead(const net_handle_t& here,
                                const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Call a function on each realizable traversal of a net graph element
     * that starts where a walk enters it.
     *
     * @param entered A traversal of the element that starts where the walk
     * enters it.
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_traversal_entered(const net_handle_t& entered,
                                    const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * Call a function on each realizable traversal that a walk leaving a
     * traversal through its end can visit next.
     *
     * Each traversal given starts where the walk enters its net graph
     * element.
     *
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool follow_right(const net_handle_t& here, const std::function<bool(const net_handle_t&)>& iteratee) const;

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
    /// The manager's snarls, by SnarlManager::snarl_number(). Net handles for
    /// snarls, sentinels and chains carry snarl numbers as their payloads.
    std::vector<const Snarl*> snarls_by_number;
    /// Tells the view each node's place in its chain, and the snarl that
    /// contains it.
    NodeDirectory directory;
    /// Tells the view which traversals of each snarl and chain are
    /// realizable.
    TraversalSurvey survey;
};

}

#endif
