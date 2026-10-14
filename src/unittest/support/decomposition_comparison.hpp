#ifndef VG_UNITTEST_DECOMPOSITION_COMPARISON_HPP_INCLUDED
#define VG_UNITTEST_DECOMPOSITION_COMPARISON_HPP_INCLUDED

/**
 * \file decomposition_comparison.hpp
 * Provides TraversalLabeler, which names net graph traversals the same way in
 * any SnarlDecomposition of a graph, and require_same_tree(), which checks
 * that two decompositions present the same snarl tree.
 */

#include <string>

#include <handlegraph/snarl_decomposition.hpp>

#include "handle.hpp"

namespace vg::unittest {

using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

/**
 * @brief A label maker that names each traversal of a net graph element by
 * the graph handles at its ends, so that a traversal gets the same label in
 * every SnarlDecomposition of one graph.
 *
 * @details Labels match up traversals between implementations that encode,
 * number or orient elements differently. Tests compare labels made from one
 * decomposition with labels made from another. The labels are:
 *
 * - For the root, or any net handle for which is_root() is true: `root`.
 * - For a node traversal: `node`, the graph handle that get_handle() gives,
 *   and the endpoints where the traversal starts and ends, as two letters
 *   from S for START, E for END and T for TIP. For example, `node 3+ SE`.
 * - For a snarl or chain traversal: `snarl` or `chain`, then the side where
 *   the traversal enters, `>`, and the side where it leaves. A side at a
 *   bound is the graph handle reading in through that bound when entering,
 *   or out through it when leaving; a side at a tip is `tip`. For example,
 *   `snarl 2+>5+` or `chain 4->tip`.
 * - For a sentinel: `sentinel`, its graph handle, `of`, and the bounds label
 *   of its snarl. For example, `sentinel 5- of 2+/5+`.
 *
 * A bounds label names a snarl or chain by its start bound reading in and
 * its end bound reading out, as `start/end`, read whichever way round sorts
 * first. For example, a snarl from 2+ to 5+ is `2+/5+` read either way.
 */
class TraversalLabeler {
public:
    /**
     * @brief Make a labeler for the traversals of a decomposition.
     *
     * @param decomposition The decomposition, which must outlive the labeler.
     * @param graph The decomposed graph, which must outlive the labeler.
     */
    TraversalLabeler(const SnarlDecomposition& decomposition, const HandleGraph& graph);

    /// @brief Label a traversal of any net graph element.
    std::string label(const net_handle_t& net) const;

    /**
     * @brief Label a net graph element, whichever way it is traversed.
     *
     * @return For a snarl or chain, `snarl` or `chain` and its bounds label;
     * for a node, `node` and its ID; and for a sentinel or the root, its
     * label().
     */
    std::string element_label(const net_handle_t& net) const;

    /**
     * @brief Label a snarl or chain by its bounds.
     *
     * @return The same label for every traversal of the element.
     */
    std::string bounds_label(const net_handle_t& net) const;

private:
    /// @brief Label the side of a snarl or chain traversal where it enters or
    /// leaves.
    std::string side(const net_handle_t& net, SnarlDecomposition::endpoint_t endpoint, bool entering) const;

    /// @brief The decomposition whose traversals are labelled.
    const SnarlDecomposition& decomposition;
    /// @brief The decomposed graph, which turns net handles into graph
    /// handles.
    const HandleGraph& graph;
};

/**
 * @brief Require that two decompositions of one graph present the same snarl
 * tree.
 *
 * @details Matches the two decompositions' traversals by their
 * TraversalLabeler labels, and requires that:
 *
 * - Both have the same labels for every realizable traversal of every snarl,
 *   chain and node under the root that starts and ends at a bound, and for
 *   each bound sentinel of every snarl, facing in and facing out.
 * - Each of those traversals has the same children in both: the same
 *   traversals in the same order for a start-to-end or end-to-start
 *   traversal of a chain, and the same set of elements otherwise.
 * - Each of those traversals, other than a sentinel, has the same parent in
 *   both. A child that reads one way along its chain must have the parent
 *   chain read the same way in both; any other parent matches either way
 *   round.
 * - Following net edges from each of those traversals, going left and going
 *   right, reaches the same set of traversals in both. Going left from a
 *   chain node, the chain node it reaches is matched whichever way it is
 *   read.
 *
 * TODO: A distance index lists a chain's children start-to-end for every
 * traversal of the chain, where the interface wants an end-to-start
 * traversal's children listed end-to-start. So the children of an
 * end-to-start traversal are worked out from the start-to-end traversal's
 * children, in both decompositions.
 *
 * TODO: Going left from a chain node to the neighbouring chain node, a
 * distance index reads the node it reaches the way the traversal reads, or
 * the way it travels, depending on how it packed the nodes into records.
 *
 * TODO: A distance index copies a chain's tip traversals from its snarls
 * before it works out the snarls' tip traversals, so its chains never have
 * any, and it can also leave out a snarl's traversal to a tip node that the
 * snarl's start reaches. So traversals that start or end at a tip are left
 * out.
 *
 * @param expected One decomposition.
 * @param actual The decomposition under test, which failure messages show
 * first.
 * @param graph The graph both decompose.
 */
void require_same_tree(const SnarlDecomposition& expected, const SnarlDecomposition& actual,
                       const HandleGraph& graph);

}

#endif
