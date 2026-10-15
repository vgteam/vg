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
 * @details The view shows every snarl and chain the manager holds. These
 * terms describe how:
 *
 * - A chain alternates nodes and snarls, starting and ending with a node.
 *   Its *chain nodes* are numbered by *position* along the chain, from 0 at
 *   its start bound. The view orients each snarl forward along its chain: a snarl's
 *   start bound is the chain node before it, and its end bound is the chain
 *   node after it. Where the manager's Snarl is backward in its chain, its
 *   start() is the view's end bound.
 * - A *cyclic chain*'s last snarl leads back into its first chain node, in
 *   the same orientation. That node is both of the chain's bounds, and the
 *   chain lists it once.
 * - Each node of the graph is either a chain node of one chain, or a *free
 *   node*. A free node sits directly in a snarl, or in the root, as a
 *   trivial chain of its own, whose start-to-end traversal takes the node in
 *   its forward orientation.
 * - A *dead end* is a side of a node with no edges on it. A traversal of a
 *   snarl or chain starts or ends at TIP where it starts or ends at a dead end
 *   inside the snarl or chain.
 *
 * Methods that take a graph use it to look up and make handles, so it must have
 * the manager's node IDs.
 *
 * @note A SnarlDistanceIndex built from the same snarl finder presents the
 * same tree, except that:
 * - The view's root holds every root chain directly. A distance index groups
 *   adjacent root chains under root snarls, for which is_root() is true.
 * - The view shows snarls that have no children. A distance index shows the
 *   two bounds of such a snarl as neighbouring chain nodes.
 */
class SnarlManagerDecomposition : public handlegraph::SnarlDecomposition {
public:
    /**
     * @brief Make a view of the snarls in a SnarlManager.
     *
     * @details Works out which traversals of each snarl and chain are
     * realizable, in time linear in the size of the graph.
     *
     * @param manager The snarls to present. The view keeps a reference to
     * it, so it must outlive the view and must not change while the view
     * exists.
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
     * end bound, which is from the Snarl's end() to its start() when the
     * Snarl is backward in its chain.
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
     * @details Besides nodes, this accepts a trivial chain, giving its node
     * oriented the way the traversal goes, and a sentinel, giving its bound
     * node oriented into the snarl if the sentinel faces in, and out of it
     * if the sentinel faces out.
     *
     * @throws std::invalid_argument If the net handle is for a snarl, a chain
     * of more than one node, or the root.
     */
    handle_t get_handle(const net_handle_t& net, const HandleGraph* handle_graph) const override;

    net_handle_t get_parent(const net_handle_t& child) const override;

    /**
     * @brief Get a bound of a snarl or chain.
     *
     * @details A cyclic chain's start and end bounds are both its first chain
     * node, oriented forward along the chain when it faces into the chain at the start
     * or out of it at the end.
     */
    net_handle_t get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const override;

    /**
     * @brief Get the traversal in the opposite direction.
     *
     * @details A sentinel flips to the same bound facing the other way.
     */
    net_handle_t flip(const net_handle_t& net) const override;

    /**
     * @brief Get the canonical traversal of a net graph element.
     *
     * @details The canonical traversal of a snarl, chain or node is the first
     * realizable one of start-to-end, start-to-tip, end-to-tip,
     * start-to-start, end-to-end and tip-to-tip. A sentinel's is its bound
     * facing in.
     */
    net_handle_t canonical(const net_handle_t& net) const override;

    endpoint_t starts_at(const net_handle_t& traversal) const override;
    endpoint_t ends_at(const net_handle_t& traversal) const override;

    /**
     * @brief Get the traversal of a snarl or chain between two of its
     * children or bounds.
     *
     * @throws std::invalid_argument If the two arguments are not in the same
     * snarl or chain, or if no realizable traversal starts and ends where
     * they say.
     */
    net_handle_t get_parent_traversal(const net_handle_t& traversal_start,
                                      const net_handle_t& traversal_end) const override;

protected:
    /**
     * @brief List the children of a net graph element.
     *
     * @details A snarl's or the root's children are its chains in the
     * manager's order, then the trivial chains of its free nodes in node ID
     * order.
     */
    bool for_each_child_impl(const net_handle_t& traversal,
                             const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * @brief List the realizable traversals of a net graph element.
     *
     * @details Traversals come in this order: start-to-start, start-to-end,
     * start-to-tip, end-to-start, end-to-end, end-to-tip, tip-to-start,
     * tip-to-end, tip-to-tip.
     */
    bool for_each_traversal_impl(const net_handle_t& item,
                                 const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * @brief List the net graph elements next to a traversal.
     *
     * @details Going right gives each realizable traversal that a walk
     * leaving through this traversal's end can visit next, starting where
     * the walk enters it. Going left gives each traversal from which going
     * right gives this one.
     * Going around a cyclic chain, a walk visits the chain's first chain node
     * and goes on.
     *
     * @note Walks the graph the view was made with, and ignores handle_graph.
     */
    bool follow_net_edges_impl(const net_handle_t& here, const HandleGraph* handle_graph, bool go_left,
                               const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * @brief List the traversals of a snarl's or chain's children that start
     * at one of its tips.
     *
     * @details For a snarl, these are each child chain's traversals that
     * start at a tip of the chain, and each traversal of a child chain or
     * trivial chain that starts at a bound facing a dead end. For a chain,
     * they are each child snarl's traversals that start at a tip of the
     * snarl.
     */
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
        /// True if the chain holds the node in reverse orientation.
        bool is_reverse;
    };

    /**
     * @brief A directory of each node's place in the decomposition: a
     * position in a chain, or the snarl or root that it sits free in.
     */
    class NodeDirectory {
    public:
        /**
         * @brief Find the place of every node.
         *
         * @param manager The snarls whose chains and contents place the
         * nodes. The directory keeps pointers into it, so it must outlive
         * the directory.
         * @param graph The graph whose nodes to place.
         *
         * @pre finish() has been called on the manager.
         */
        NodeDirectory(const SnarlManager& manager, const HandleGraph& graph);

        /**
         * @brief Get a node's place in its chain.
         *
         * @return The place, or nullopt for a free node.
         */
        std::optional<ChainPosition> chain_position(nid_t id) const;

        /**
         * @brief Get the snarl that directly contains a node.
         *
         * @return The snarl that holds the node's chain, or that the node is
         * free in, or null if that is the root.
         */
        const Snarl* containing_snarl(nid_t id) const;

        /// Get the free nodes directly in a snarl, or in the root if given
        /// null, in ID order.
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

    /**
     * @brief A set of traversal kinds, each named by the endpoints where it
     * starts and ends.
     */
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
         * @post Kinds come in the order for_each_traversal_impl() promises.
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
        /// One bit for each traversal kind, at index 3 * start + end.
        uint16_t bits = 0;
    };

    /**
     * @brief A survey of which traversals of each snarl and chain are
     * realizable, made by walking the graph.
     *
     * @details The survey walks each snarl among its children, crossing,
     * turning around in, or stopping at a tip in each child chain as that
     * chain's own traversals allow. It walks each chain over its chain
     * nodes, through each snarl as that snarl's traversals allow. So it
     * surveys inner snarls and chains before the ones that hold them.
     */
    class TraversalSurvey {
    public:
        /**
         * @brief Survey every snarl and chain in a manager.
         *
         * @param manager The snarls and chains to survey. The survey keeps
         * pointers into it, so it must outlive the survey.
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
        /// Work out the realizable traversals of a snarl.
        /// @pre Every chain in the snarl has been surveyed.
        TraversalSet survey_snarl(const Snarl* snarl) const;
        /// Work out the realizable traversals of a chain.
        /// @pre Every snarl in the chain has been surveyed.
        TraversalSet survey_chain(const Chain& chain) const;

        /// The endpoints that a walk can reach.
        struct Reach {
            bool start = false;
            bool end = false;
            bool tip = false;
        };

        /**
         * @brief Walk among a snarl's children.
         *
         * @param snarl The snarl.
         * @param starts Handles on nodes inside the snarl. The walk starts on
         * each one, and leaves it through its right side.
         * @return Whether the walk leaves the snarl through its start
         * bound or its end bound, and whether it reaches a dead end.
         */
        Reach walk_snarl(const Snarl* snarl, const std::vector<handle_t>& starts) const;

        /**
         * @brief Walk over a chain's chain nodes.
         *
         * @param chain The chain.
         * @param starts Places to start: each a chain node's position, and
         * true to leave it going in reverse along the chain, or false to
         * leave it going forward along the chain.
         * @return Whether the walk leaves the chain through its start bound or
         * its end bound, and whether it reaches a tip in one of its snarls.
         */
        Reach walk_chain(const Chain& chain, const std::vector<std::pair<size_t, bool>>& starts) const;

        /**
         * @brief Enter a child chain or free node.
         *
         * @param entering A handle that enters a chain at one of its bounds,
         * or a handle on a free node.
         * @param exits Where to add the child's exits. An exit is a handle
         * that a walk leaves the child by, through the handle's right side.
         * A chain's exits are its bound nodes, oriented out of the chain, at
         * each bound where a walk entering the chain can leave it. A free
         * node's one exit is the entering handle.
         * @return True if a walk entering the chain can reach a tip inside it.
         */
        bool enter_child(const handle_t& entering, std::vector<handle_t>& exits) const;

        /// The snarls and chains being surveyed.
        const SnarlManager& manager;
        /// The graph the snarls were found in.
        const HandleGraph& graph;
        /// Tells the survey which chain, if any, each node is a chain node of.
        const NodeDirectory& directory;
        /// The traversals of each snarl, by snarl.
        std::unordered_map<const Snarl*, TraversalSet> snarls;
        /// The traversals of each chain, by chain.
        std::unordered_map<const Chain*, TraversalSet> chains;
    };

    /// Get a snarl's start bound oriented into it, or its end bound oriented
    /// out of it, in the view's orientation, from the view's graph.
    handle_t snarl_bound(const Snarl* snarl, bool get_end) const;
    /// Get a snarl's net handle, start-to-end if forward and end-to-start
    /// otherwise.
    net_handle_t snarl_net(const Snarl* snarl, bool forward) const;
    /// Get a chain's start-to-end net handle.
    net_handle_t chain_net(const Chain& chain) const;
    /// Get the chain that a chain net handle is for.
    const Chain& chain_of_net(const net_handle_t& net) const;
    /// Get the net handle for a chain node, oriented forward along the chain
    /// or in reverse along it.
    net_handle_t chain_node_net(const Chain& chain, size_t position, bool is_reverse) const;

    /**
     * @brief Call a function on each child of a snarl, or of the root if
     * given null, as for_each_child_impl() promises.
     *
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_snarl_child(const Snarl* parent, const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * @brief Call a function on each chain node and snarl of a chain.
     *
     * @details For a start-to-end or end-to-start traversal, children come in
     * the order the traversal visits them, oriented the way it visits them. For any other
     * traversal, they come in start-to-end order, each start-to-end.
     *
     * @return False if the function stopped the loop, and true otherwise.
     */
    bool for_each_chain_child(const net_handle_t& traversal,
                              const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * @brief Get the next child in a chain after a chain node or snarl.
     *
     * @param here A start-to-end or end-to-start traversal of a chain node,
     * or any traversal of a snarl.
     * @return The neighbouring chain node or snarl that a walk leaving
     * through the end of here enters, oriented in the direction of the walk. A
     * walk leaving a snarl through a tip, or leaving a chain at its bound,
     * enters nothing, and gives nullopt.
     */
    std::optional<net_handle_t> next_in_chain(const net_handle_t& here) const;

    /**
     * @brief Get the child of a snarl or the root that a handle enters.
     *
     * @param handle A handle reached by following an edge out of one of the
     * parent's children or out of one of its bounds.
     * @param parent The snarl, or null for the root.
     * @return
     * - The traversal of the sibling chain or trivial chain that the handle
     *   enters at one of its bounds.
     * - The parent's bound sentinel facing out, if the handle leaves the
     *   parent.
     */
    std::optional<net_handle_t> sibling_entered_by(const handle_t& handle, const Snarl* parent) const;

    /**
     * @brief Call a function on each sibling that a walk leaving a traversal
     * through its end enters.
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

    /// Go right from a traversal, as follow_net_edges_impl() does.
    bool follow_right(const net_handle_t& here, const std::function<bool(const net_handle_t&)>& iteratee) const;

    /**
     * @brief Find the parent of a net graph element, and the endpoint of the
     * parent where a traversal starting or ending at the element starts or
     * ends.
     *
     * @param net A sentinel, a bound of a chain or trivial chain, or a child
     * of a chain or snarl.
     * @param is_start True to find where a traversal starts, so that a bound
     * must be oriented into its chain, and false to find where it ends, so
     * that a bound must be oriented out of it.
     * @return The parent, start-to-end, and the endpoint: a sentinel's bound,
     * a chain's or trivial chain's bound, or TIP for a child.
     *
     * @throws std::invalid_argument If the element cannot start or end a
     * traversal of its parent.
     */
    std::pair<net_handle_t, endpoint_t> traversal_bound(const net_handle_t& net, bool is_start) const;

    /// Get the realizable traversals of a net graph element.
    TraversalSet traversals_of(const net_handle_t& item) const;

    /// The snarls being presented.
    const SnarlManager& manager;
    /// The graph the snarls were found in, which follow_net_edges() walks.
    const HandleGraph& graph;
    /// The manager's snarls, by SnarlManager::snarl_number().
    std::vector<const Snarl*> snarls;
    /// Tells the view where each node sits: in which chain, or free in which
    /// snarl.
    NodeDirectory directory;
    /// Tells the view which traversals of each snarl and chain are
    /// realizable.
    TraversalSurvey survey;
};

}

#endif
