#ifndef VG_DECOMPOSITION_NEIGHBORS_HPP_INCLUDED
#define VG_DECOMPOSITION_NEIGHBORS_HPP_INCLUDED

/** \file
 * Stepping from a traversal in a snarl decomposition to each neighboring net
 * graph element.
 */

#include <functional>

#include <handlegraph/iteratee.hpp>
#include <handlegraph/snarl_decomposition.hpp>

#include "handle.hpp"

namespace vg {

/**
 * Call a function on one traversal of each neighboring net graph element of a
 * traversal.
 *
 * Use this instead of follow_net_edges() when you want just one traversal of
 * each next element, as when stepping along a chain. follow_net_edges() gives
 * every realizable traversal of each next element: one that crosses it, one
 * that turns around in it, one that ends at a tip in it. This gives the one
 * that crosses the element if there is one, and otherwise the first one
 * follow_net_edges() gives. The elements come in an arbitrary order.
 *
 * @param decomposition The snarl decomposition of graph.
 * @param here The traversal to step from.
 * @param graph The graph whose edges join here to its neighboring elements.
 * @param go_left True to step out of the left side of the traversal, false
 * for the right side, as for follow_net_edges().
 * @param iteratee A function to call on the traversal of each neighboring
 * element, which may return false to stop.
 * @return False if the function stopped the loop, and true otherwise.
 */
template<typename Iteratee>
bool for_each_neighbor(const handlegraph::SnarlDecomposition& decomposition, const handlegraph::net_handle_t& here,
                       const HandleGraph* graph, bool go_left, const Iteratee& iteratee);

/// Implement for_each_neighbor() for an iteratee that returns bool.
bool for_each_neighbor_impl(const handlegraph::SnarlDecomposition& decomposition,
                            const handlegraph::net_handle_t& here, const HandleGraph* graph, bool go_left,
                            const std::function<bool(const handlegraph::net_handle_t&)>& iteratee);

template<typename Iteratee>
bool for_each_neighbor(const handlegraph::SnarlDecomposition& decomposition, const handlegraph::net_handle_t& here,
                       const HandleGraph* graph, bool go_left, const Iteratee& iteratee) {
    return for_each_neighbor_impl(decomposition, here, graph, go_left,
                                  handlegraph::BoolReturningWrapper<Iteratee>::wrap(iteratee));
}

}

#endif
