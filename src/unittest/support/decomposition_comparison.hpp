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
 * @brief Name a snarl, chain or node the same way in every implementation.
 *
 * @details The name comes from node IDs and bound handles, so every
 * traversal of the element gives the same name too. A node is named `node` and its ID, as in `node 3`.
 *
 * A snarl or chain is named `snarl` or `chain` and its bound pair. To make
 * the bound pair, take a traversal of the element that goes from one bound
 * to the other. Write the handle that enters the element at the first bound
 * and the handle that leaves at the second, as `2+/5+` for a snarl entered
 * at node 2 forward and left at node 5 forward. The traversal the other way
 * gives `5-/2-`. The bound pair is whichever of the two sorts first, so
 * every traversal of the element gives the same name.
 *
 * @param element Any traversal of a snarl, chain or node.
 */
std::string element_name(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                         const net_handle_t& element);

/**
 * @brief The shape of a snarl tree: each element's name, mapped to the names
 * of its children.
 *
 * @details The root is named `root`. A snarl's or the root's children are
 * sorted by name. A chain's children are in the order that a traversal of
 * the chain from one bound to the other passes them; of the two such orders,
 * the one that sorts first is kept. A cyclic chain's children start at its
 * bound node in both orders.
 */
using TreeShape = std::map<std::string, std::vector<std::string>>;

/**
 * @brief Describe the shape of the snarl tree that a decomposition presents.
 *
 * @details A snarl with no children is left out, both as an element and as
 * a child of its chain, because a SnarlDistanceIndex does not present one.
 * Its chain then lists its two bounds as neighbours. A net handle other than
 * get_root() for which is_root() is true, like the root snarls that a
 * SnarlDistanceIndex groups root chains under, counts as part of the root.
 */
TreeShape tree_shape(const SnarlDecomposition& decomposition, const HandleGraph& graph);

}

#endif
