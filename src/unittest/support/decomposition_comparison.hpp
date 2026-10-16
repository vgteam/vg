#ifndef VG_UNITTEST_DECOMPOSITION_COMPARISON_HPP_INCLUDED
#define VG_UNITTEST_DECOMPOSITION_COMPARISON_HPP_INCLUDED

/**
 * \file decomposition_comparison.hpp
 * Comparing the snarl trees that two SnarlDecomposition implementations
 * present for one graph.
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
 * @brief Get a name for a snarl, chain or node, which is the same in each
 * implementation and for each traversal of the element.
 *
 * @details A node is named `node` and its ID, as in `node 3`.
 *
 * A snarl or chain is named `snarl` or `chain` and its bound pair. The bound
 * pair is written as the handle of the start bound facing in, a slash, and
 * the handle of the end bound facing out, as in `2+/5+` for a snarl from
 * node 2 forward to node 5 forward. The same bounds written in the other
 * direction are `5-/2-`. The bound pair is whichever of these two sorts
 * first.
 *
 * @param element Any traversal of a snarl, chain or node.
 */
std::string element_name(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& element);

/**
 * @brief The shape of a snarl tree: the name of the root and of each snarl
 * and chain, mapped to the names of its children.
 *
 * @details The root is named `root`, and the other elements as
 * element_name() names them. A snarl's or the root's children are
 * sorted by name. A chain's children are in the order that a traversal of
 * the chain from one bound to the other visits them; of the two such orders,
 * the one that sorts first is kept. A cyclic chain's children start at its
 * bound node in both orders.
 */
using TreeShape = std::map<std::string, std::vector<std::string>>;

/**
 * @brief Describe the shape of the snarl tree that a decomposition presents.
 *
 * @details A snarl with no children is left out, both as an element and as
 * a child of its chain, so that the shape of an implementation that presents
 * such snarls matches the shape of a SnarlDistanceIndex. The snarl's chain
 * then lists the snarl's two bounds next to each other.
 *
 * A child for which is_root() is true, such as a root snarl that a
 * SnarlDistanceIndex groups root chains under, is replaced by its own
 * children.
 */
TreeShape tree_shape(const SnarlDecomposition& decomposition, const HandleGraph& graph);

}

#endif
