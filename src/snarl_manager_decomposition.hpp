#ifndef VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED
#define VG_SNARL_MANAGER_DECOMPOSITION_HPP_INCLUDED

/** \file
 * Defines SnarlManagerDecomposition, which presents a SnarlManager through
 * the handlegraph::SnarlDecomposition interface.
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
 * @brief An adapter that answers SnarlDecomposition questions about the snarls
 * in a SnarlManager.
 *
 * @details The view presents the same snarl tree that a SnarlDistanceIndex
 * built from the same snarl finder presents: the same snarls, the same
 * chains, and the same children in each, except as the notes below say.
 * These terms describe the tree:
 *
 * - A chain alternates nodes and snarls, starting and ending with a node.
 *   Its *chain nodes* are numbered by *position*: position i is the node
 *   read into the chain's snarl i, going along the chain, and the last
 *   position is the node read out of its last snarl. The first and last
 *   chain nodes are the chain's bounds.
 * - A *cyclic chain* returns to where it started: its last snarl leads back
 *   into its first chain node, read the same way. Its first chain node is
 *   both of its bounds and has only one position, 0.
 * - Every node of the graph is either a chain node of exactly one chain, or
 *   a *free node*, which sits alone in its parent snarl, or in the root, and
 *   is shown as a trivial chain holding only that node. A trivial chain's
 *   start-to-end traversal reads its node forward.
 * - A *hidden snarl* has no children: no chains and no free nodes. The view
 *   leaves it out of the tree, so its chain shows its two bounds as
 *   neighbouring chain nodes.
 * - A *dead end* is a side of a node with no edges on it. A *tip node* is a
 *   free node with a dead end. A traversal of a snarl ends at TIP when it
 *   stops at a tip node's dead end. A traversal of a chain ends at TIP when
 *   it crosses the chain's snarls into one where it can reach a tip.
 *
 * Walks inside a snarl, which decide where its traversals can reach tips,
 * cross or turn around in each child chain as the manager's flags on the
 * chain's snarls allow. The traversals a chain lists follow a distance
 * index's rules instead, which the TODOs below describe.
 *
 * @note The view and a SnarlDistanceIndex differ in these ways:
 * - The view refuses unary snarls, whose end is their start reversed.
 * - The view's root holds every root chain directly. A distance index groups
 *   adjacent root chains under extra root snarls, for which is_root() is
 *   true and canonical() gives the root.
 * - Snarls are read as the manager reads them, which can be backward along
 *   their chain. A distance index reads every snarl along its chain.
 * - Snarls and chains that can reach a tip have tip traversals.
 *
 * TODO: These behaviours copy a distance index, and should change once the
 * distance index's [snarl_distance_interface] tests pass:
 * - A snarl with no children is hidden.
 * - A cyclic chain's end bound is its start bound's node read end-to-end
 *   (or start-to-start).
 * - A chain can turn around at its start if any of its shown snarls can turn
 *   around at the side facing the chain's start.
 * - Going left from a chain node to the neighbouring chain node, the node
 *   reached is read the way the traversal reads.
 *
 * Methods that take a graph use it to read and make handles, so it must have
 * the manager's node IDs. follow_net_edges() ignores the graph it is given
 * and walks the graph the view was made with.
 */
class SnarlManagerDecomposition : public handlegraph::SnarlDecomposition {
public:
    /**
     * @brief Make a view of the snarls in a SnarlManager.
     *
     * @details Builds the directory of where each node sits, and works out
     * which traversals of each snarl and chain are realizable, which takes
     * time linear in the size of the graph.
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
     * @brief Get the net handle for a snarl.
     *
     * @param snarl A snarl owned by the manager the view presents.
     * @return A start-to-end traversal of the snarl, even if the snarl cannot
     * be crossed from start to end.
     *
     * @throws std::invalid_argument If the snarl is hidden, so the tree does
     * not contain it.
     */
    net_handle_t get_snarl_net(const Snarl* snarl) const;

    /**
     * @brief Convert a net handle for a snarl to the manager's Snarl.
     *
     * @param net Any traversal of a snarl in the tree.
     * @return The snarl, owned by the manager the view presents.
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
     * @brief Get the graph handle for a node, a trivial chain, or a snarl
     * bound sentinel.
     *
     * @return For a node or trivial chain, the node read the way the
     * traversal leaves it. For a sentinel, its bounding node read into the
     * snarl if the sentinel faces in, and out of the snarl if it faces out.
     *
     * @throws std::invalid_argument If the net handle is for a snarl, a
     * chain with snarls in it, or the root.
     */
    handle_t get_handle(const net_handle_t& net, const HandleGraph* handle_graph) const override;

    /**
     * @brief Get the parent of a net graph element.
     *
     * @return For a chain or a free node's trivial chain at the root, the
     * root itself.
     *
     * @throws std::invalid_argument If given the root.
     */
    net_handle_t get_parent(const net_handle_t& child) const override;

    /**
     * @brief Get the bound of a snarl or chain.
     *
     * @return For a cyclic chain, the node at position 0. Facing out at the
     * end, it is read end-to-end if the chain reads it forward, and
     * start-to-start otherwise; facing in at the end, the other of the two.
     *
     * @throws std::invalid_argument If given the root, a node or a sentinel.
     */
    net_handle_t get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const override;

    /**
     * @brief Get the traversal in the opposite direction.
     *
     * @return For a sentinel, the same bound facing the other way. For
     * anything else, the traversal with its start and end swapped, which
     * leaves start-to-start, end-to-end and tip-to-tip traversals unchanged.
     */
    net_handle_t flip(const net_handle_t& net) const override;

    /**
     * @brief Get the canonical traversal of a net graph element.
     *
     * @return For a snarl, chain, node or trivial chain, the first realizable
     * traversal in this order: start-to-end, start-to-tip, end-to-tip,
     * start-to-start, end-to-end, tip-to-tip. If none is realizable,
     * start-to-end. For a sentinel, its bound facing in.
     */
    net_handle_t canonical(const net_handle_t& net) const override;

    endpoint_t starts_at(const net_handle_t& traversal) const override;
    endpoint_t ends_at(const net_handle_t& traversal) const override;

    /**
     * @brief Get the traversal of a snarl between two of its sentinels, or of
     * a chain or trivial chain between two of its bounds.
     *
     * @details Each argument says where the traversal starts or ends:
     * - a sentinel of a snarl, at its bound;
     * - a traversal of a free node, where its trivial chain's traversal
     *   starts or ends;
     * - a chain's bound read into the chain to start, or out of it to end,
     *   at that bound;
     * - a snarl in a chain, or a chain or trivial chain in a snarl, at a tip
     *   of its parent.
     *
     * @throws std::invalid_argument If an argument is none of these, if the
     * two do not have the same parent, or if the traversal they describe is
     * not realizable.
     */
    net_handle_t get_parent_traversal(const net_handle_t& traversal_start,
                                      const net_handle_t& traversal_end) const override;

protected:
    /**
     * @brief List the children of a net graph element.
     *
     * @post A snarl's or the root's children are its chains in the
     * manager's order, then the trivial chains of its free nodes in node ID
     * order.
     * @post A cyclic chain read start-to-end or end-to-start lists its
     * position 0 node first, and only once.
     * @post A trivial chain's only child is its node. Nodes and sentinels
     * have none.
     */
    bool for_each_child_impl(const net_handle_t& traversal,
                             const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * @brief List the realizable traversals of a net graph element.
     *
     * @post Traversals come in this order, leaving out those that are not
     * realizable: start-to-start, start-to-end, start-to-tip, end-to-start,
     * end-to-end, end-to-tip, tip-to-start, tip-to-end, tip-to-tip.
     * @post Nodes and trivial chains have only start-to-end and end-to-start,
     * the root has only tip-to-tip, and a sentinel has only itself.
     */
    bool for_each_traversal_impl(const net_handle_t& item,
                                 const std::function<bool(const net_handle_t&)>& iteratee) const override;

    /**
     * @brief List the net graph elements reached by leaving a traversal
     * through its end, or through its start if going left.
     *
     * @details Inside a chain, the walk reaches the neighbouring chain node or
     * shown snarl. Among a snarl's children, it follows the graph's edges to
     * the sibling chains, trivial chains and outward-facing bound sentinels
     * they lead into. At the root, it reaches root chains and trivial chains.
     *
     * @post A chain node or snarl reached inside a chain is read in the
     * direction of travel, except as the class's TODOs describe.
     * @post Leaving a traversal through a tip reaches nothing, as does
     * leaving a cyclic chain, an outward-facing sentinel, or the end of a
     * chain from inside it.
     * @note Walks the graph the view was made with, ignoring handle_graph.
     */
    bool follow_net_edges_impl(const net_handle_t& here, const HandleGraph* handle_graph, bool go_left,
                               const std::function<bool(const net_handle_t&)>& iteratee) const override;

private:

    /// @brief The place of a chain node in its chain.
    struct ChainPosition {
        /// @brief The chain, owned by the manager.
        const Chain* chain;
        /// @brief The node's position: 0 for the start bound, counting up
        /// along the chain, and less than chain_node_count(*chain).
        size_t position;
        /// @brief True if the chain reads the node in reverse.
        bool is_reverse;
    };

    /**
     * @brief A directory that lists each node's place in the decomposition:
     * a position in a chain, or a snarl or the root that it sits free in.
     *
     * @details Filled in when made, from every chain and snarl in the manager
     * and every node in the graph.
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
         * @brief Get the snarl directly containing a node.
         *
         * @return The parent snarl of the node's chain, or the snarl the node
         * is free in, or null if that is the root.
         */
        const Snarl* containing_snarl(nid_t id) const;

        /**
         * @brief Get the free nodes directly in a snarl, or in the root if
         * given null, in ID order.
         */
        const std::vector<nid_t>& free_nodes(const Snarl* snarl) const;

    private:
        /// @brief One node's place in the decomposition.
        struct Entry {
            /// @brief The snarl directly containing the node, or null for the
            /// root.
            const Snarl* containing_snarl;
            /// @brief The node's place in its chain, or nullopt for a free node.
            std::optional<ChainPosition> chain_position;
        };

        /// @brief The entry for every node of the graph, by ID.
        std::unordered_map<nid_t, Entry> entries;
        /// @brief The free nodes in each snarl that has any, and in the root
        /// under null, each list in ID order.
        std::unordered_map<const Snarl*, std::vector<nid_t>> free_nodes_by_snarl;
    };

    /**
     * @brief A set of traversal kinds, each a pair of the endpoints where a
     * traversal starts and ends.
     */
    class TraversalSet {
    public:
        /// @brief Add the traversal kind that starts and ends at the given
        /// endpoints.
        void add(endpoint_t start, endpoint_t end);
        /// @brief Add the traversal kind between two endpoints, in both
        /// directions.
        void add_both_ways(endpoint_t one, endpoint_t other);
        /// @brief Return true if the set holds the traversal kind.
        bool contains(endpoint_t start, endpoint_t end) const;
        /**
         * @brief Call a function on each traversal kind in the set, in the
         * order for_each_traversal_impl() promises, until it returns false.
         *
         * @param iteratee A function taking the start and end endpoints and
         * returning false to stop.
         * @return False if the function stopped the loop, and true otherwise.
         */
        template<typename Iteratee>
        bool for_each(const Iteratee& iteratee) const {
            // Counting up the bit index gives the order a distance index
            // lists them in.
            for (endpoint_t start : {START, END, TIP}) {
                for (endpoint_t end : {START, END, TIP}) {
                    if (contains(start, end) && !iteratee(start, end)) {
                        return false;
                    }
                }
            }
            return true;
        }

        /**
         * @brief Work out a chain's traversals from those of its shown snarls.
         *
         * @param shown_snarls Each shown snarl's traversals, and whether the
         * snarl is backward in the chain, in chain order.
         */
        static TraversalSet for_chain(const std::vector<std::pair<TraversalSet, bool>>& shown_snarls);

    private:
        /**
         * @brief Return true if a walk entering a chain at one bound can
         * cross the chain's snarls into one where it reaches a tip.
         *
         * @param shown_snarls As for for_chain().
         * @param from_start True to enter at the chain's start bound, and
         * false to enter at its end bound.
         */
        static bool reaches_tip(const std::vector<std::pair<TraversalSet, bool>>& shown_snarls, bool from_start);

        /**
         * @brief Return true if a walk can go from a tip in one of a chain's
         * snarls, along the chain, to a tip in the same or a later snarl.
         *
         * @param shown_snarls As for for_chain().
         */
        static bool links_tips(const std::vector<std::pair<TraversalSet, bool>>& shown_snarls);

        /// @brief One bit for each traversal kind, at index 3 * start + end.
        uint16_t bits = 0;
    };

    /**
     * @brief A surveyor that works out which traversals of a snarl are
     * realizable, by walking the snarl's net graph.
     *
     * @details The walk follows the graph's edges among the snarl's
     * children. It crosses or turns around in each child chain as the
     * manager's flags on that chain's snarls allow.
     */
    class SnarlSurveyor {
    public:
        /**
         * @brief Make a surveyor for the snarls of a graph.
         *
         * @param graph The graph to walk, which must outlive the surveyor.
         * @param directory The places of the graph's nodes, which must outlive
         * the surveyor.
         */
        SnarlSurveyor(const HandleGraph& graph, const NodeDirectory& directory);

        /// @brief Work out the realizable traversals of a snarl.
        TraversalSet survey(const Snarl* snarl) const;

    private:
        /**
         * @brief Find out whether a walk inside a snarl can stop at a tip
         * node's dead end.
         *
         * @param snarl The snarl to walk in. The walk stops at its bounds.
         * @param starts Handles read in the snarl, each of which the walk
         * leaves through its right side.
         */
        bool reaches_tip(const Snarl* snarl, const std::vector<handle_t>& starts) const;

        /**
         * @brief List the handles a walk can leave a child chain or free node
         * through, after reading into it.
         *
         * @param entering A handle reading into a child chain at one of its
         * bounds, or reading a free node.
         * @return For a free node, the handle itself. For a chain, the
         * chain's bounds read out of it, wherever the walk can cross or turn
         * around to reach them.
         */
        std::vector<handle_t> exits_of_child(const handle_t& entering) const;

        /// @brief The graph to walk.
        const HandleGraph& graph;
        /// @brief The directory of node places, which tells the surveyor
        /// which child of the snarl a handle reads into.
        const NodeDirectory& directory;
    };

    /// @brief A record of one of the manager's snarls: whether the view shows
    /// it, and how it can be traversed.
    struct SnarlRecord {
        /// @brief The snarl, owned by the manager.
        const Snarl* snarl = nullptr;
        /// @brief True if the snarl has no children, so the view leaves it out.
        bool hidden = false;
        /// @brief The snarl's realizable traversals, empty if it is hidden.
        TraversalSet traversals;
    };

    /// @brief Get the record for a snarl or sentinel net handle, or the
    /// record for the first snarl of a chain net handle's chain.
    const SnarlRecord& record_of(const net_handle_t& net) const;
    /// @brief Get a snarl's net handle, start-to-end if forward and
    /// end-to-start otherwise.
    net_handle_t snarl_net(const Snarl* snarl, bool forward) const;
    /// @brief Get a chain's start-to-end net handle.
    net_handle_t chain_net(const Chain& chain) const;
    /// @brief Get the chain a chain net handle is for.
    const Chain& chain_of_net(const net_handle_t& net) const;
    /**
     * @brief Get the net handle for the chain node at a position, read along
     * the chain or in reverse.
     *
     * @param position A position, or for a cyclic chain, chain.size(), which
     * a walk reaches by going around the chain to position 0 again. That
     * gives the chain's end bound facing out.
     */
    net_handle_t chain_node_net(const Chain& chain, size_t position, bool is_reverse) const;
    /// @brief Get the net handle for a chain's snarl, read along the chain or
    /// in reverse.
    net_handle_t chain_snarl_net(const Chain& chain, size_t rank, bool is_reverse) const;

    /**
     * @brief List the children of a snarl, or of the root if given null.
     *
     * @return The manager's child chains in the manager's order, then the
     * trivial chains of the free nodes in ID order, all start-to-end.
     */
    std::vector<net_handle_t> snarl_children(const Snarl* parent) const;

    /**
     * @brief List a chain's nodes and shown snarls.
     *
     * @return For a start-to-end or end-to-start traversal, the children in
     * the order and orientation the traversal reads them. For any other
     * traversal, the children in start-to-end order, each start-to-end.
     */
    std::vector<net_handle_t> chain_children(const net_handle_t& traversal) const;

    /**
     * @brief Get the chain's next child after a chain node or shown snarl.
     *
     * @param here A start-to-end or end-to-start traversal of a chain node
     * or a shown snarl.
     * @return The neighbouring chain node or shown snarl that a walk leaving
     * through the end of here reaches, or nullopt at the end of the chain.
     */
    std::optional<net_handle_t> next_in_chain(const net_handle_t& here) const;

    /**
     * @brief Get the child of a snarl, or the root, that a handle reads into.
     *
     * @param handle A handle reached by following an edge out of one of the
     * parent's children or out of one of its bounds.
     * @param parent The snarl, or null for the root.
     * @return The traversal of the sibling chain or trivial chain that the
     * handle starts, or the parent's bound sentinel facing out if the handle
     * reads out of the parent, or nullopt if the handle reads into the middle
     * of a chain.
     */
    std::optional<net_handle_t> sibling_entered_by(const handle_t& handle, const Snarl* parent) const;

    /**
     * @brief List the siblings reached by following graph edges out of the
     * end of a traversal.
     *
     * @param here A traversal of a chain or trivial chain, or a sentinel
     * facing into its snarl.
     * @return One sibling_entered_by() result for each edge, leaving out
     * edges that read into the middle of a chain. Empty for a traversal
     * ending at a tip, a cyclic chain, or an outward-facing sentinel.
     */
    std::vector<net_handle_t> siblings_ahead(const net_handle_t& here) const;

    /**
     * @brief Find where a traversal of a net graph element's parent starts or
     * ends, if it starts or ends at the element.
     *
     * @param net A sentinel, a bound of a chain or trivial chain, or a child
     * of a chain or snarl.
     * @param is_start True to find where the traversal starts, so that a
     * bound must read into its chain, and false to find where it ends, so
     * that a bound must read out of it.
     * @return The parent, start-to-end, and the endpoint: a sentinel's bound,
     * a chain's or trivial chain's bound, or TIP for a child.
     *
     * @throws std::invalid_argument If the element cannot start or end a
     * traversal of its parent.
     */
    std::pair<net_handle_t, endpoint_t> traversal_bound(const net_handle_t& net, bool is_start) const;

    /// @brief Get the realizable traversals of a net graph element.
    TraversalSet traversals_of(const net_handle_t& item) const;

    /// @brief The snarls being presented.
    const SnarlManager& manager;
    /// @brief The graph the snarls were found in, which follow_net_edges()
    /// walks.
    const HandleGraph& graph;
    /// @brief Answers where each node sits: in which chain, or free in which
    /// snarl.
    NodeDirectory directory;
    /// @brief The manager's snarls, indexed by SnarlManager::snarl_number().
    std::vector<SnarlRecord> snarls;
    /// @brief The realizable traversals of every chain.
    std::unordered_map<const Chain*, TraversalSet> chains;
};

}

#endif
