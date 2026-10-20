#ifndef VG_UNITTEST_DECOMPOSITION_COMPARISON_HPP_INCLUDED
#define VG_UNITTEST_DECOMPOSITION_COMPARISON_HPP_INCLUDED

/**
 * \file decomposition_comparison.hpp
 * Names for the net graph elements of a snarl tree that every
 * SnarlDecomposition implementation shares, so tests can compare
 * implementations.
 */

#include <map>
#include <string>
#include <vector>

#include <handlegraph/snarl_decomposition.hpp>

#include "handle.hpp"

namespace vg::unittest {

using handlegraph::net_handle_t;
using handlegraph::SnarlDecomposition;

/**
 * Name a net graph element the same way in every implementation.
 *
 * Every traversal of one element gets the same name. A node is
 * named `node` and its ID, as in `node 3`.
 *
 * A snarl or chain is named `snarl` or `chain` and its bounds. The bounds are
 * written as the handle of the start bound facing in, a slash, and the handle
 * of the end bound facing out, as in `2+/5+` for a snarl from node 2 forward
 * to node 5 forward. Written from the other end, the same bounds are `5-/2-`.
 * The name uses whichever of the two sorts first.
 *
 * @param element Any traversal of a net graph element.
 */
std::string element_name(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& element);

/**
 * The shape of a snarl tree, as a map from the name of the root and of each
 * snarl and chain to the names of its children. The map also holds, under the
 * name of each snarl and chain followed by ` traversals`, its realizable
 * traversals, each as two letters S, E or T for where it starts and ends,
 * and under the name followed by ` tippy children`, the traversals that
 * for_each_tippy_child() gives that start at a tip, each as the name of its
 * element and its letters.
 *
 * The root is named `root`, and everything else as element_name() names it.
 * The children of a snarl or of the root are sorted by name. The children of
 * a chain are in the order that a traversal of the chain visits them. Of the
 * two traversals, start-to-end and end-to-start, the map keeps the one whose
 * list sorts first. The list for a looping chain starts at its boundary node
 * in both orders.
 */
using TreeShape = std::map<std::string, std::vector<std::string>>;

/**
 * Describe the shape of the snarl tree that a decomposition presents.
 *
 * Snarls with no children are left out, as keys and as children of their
 * chains, because a SnarlDistanceIndex does not present them. The chain of
 * such a snarl lists the two boundary nodes of the snarl next to each other.
 *
 * A child for which is_root() is true is replaced by its own children. A
 * SnarlDistanceIndex groups connected top-level chains under such a root
 * snarl.
 */
TreeShape tree_shape(const SnarlDecomposition& decomposition, const HandleGraph& graph);

}

#endif
