#ifndef VG_FOLLOW_NET_ELEMENTS_HPP_INCLUDED
#define VG_FOLLOW_NET_ELEMENTS_HPP_INCLUDED

/** \file
 * Walking a snarl decomposition one net graph element at a time, for code
 * that steps along a chain.
 */

#include <functional>

#include <handlegraph/iteratee.hpp>
#include <handlegraph/snarl_decomposition.hpp>

#include "handle.hpp"

namespace vg {

/**
 * @brief Call a function on each net graph element next to a traversal, with
 * one traversal of each element.
 *
 * @details follow_net_edges() gives every realizable traversal of each
 * element that a walk can visit next: crossing it, turning around in it, or
 * stopping at a tip in it. Code that steps along a chain wants one traversal
 * of each element, so this gives the traversal that crosses the element, from
 * one end to the other, when follow_net_edges() gives one. Otherwise it gives
 * one of the traversals follow_net_edges() gives for the element.
 *
 * @param go_left True to step left from the traversal, and false to step
 * right, as for follow_net_edges().
 * @param iteratee A function to call on each element's traversal, which
 * may return false to stop.
 * @return False if the function stopped the loop, and true otherwise.
 */
template<typename Iteratee>
bool follow_net_elements(const handlegraph::SnarlDecomposition& decomposition, const handlegraph::net_handle_t& here,
                         const HandleGraph* graph, bool go_left, const Iteratee& iteratee);

/// Call a bool-returning function on each element, as follow_net_elements()
/// does.
bool follow_net_elements_impl(const handlegraph::SnarlDecomposition& decomposition,
                              const handlegraph::net_handle_t& here, const HandleGraph* graph, bool go_left,
                              const std::function<bool(const handlegraph::net_handle_t&)>& iteratee);

template<typename Iteratee>
bool follow_net_elements(const handlegraph::SnarlDecomposition& decomposition, const handlegraph::net_handle_t& here,
                         const HandleGraph* graph, bool go_left, const Iteratee& iteratee) {
    return follow_net_elements_impl(decomposition, here, graph, go_left,
                                    handlegraph::BoolReturningWrapper<Iteratee>::wrap(iteratee));
}

}

#endif
