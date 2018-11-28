#ifndef VG_ALGORITHMS_TOPOLOGICAL_SORT_HPP_INCLUDED
#define VG_ALGORITHMS_TOPOLOGICAL_SORT_HPP_INCLUDED

/**
 * \file topological_sort.hpp
 *
 * Defines a topological sort algorithm for handle graphs.
 */

#include <unordered_map>

#include "../position.hpp"
#include "../cached_position.hpp"
#include "../vg.pb.h"
#include "../hash_map.hpp"
#include "../handle.hpp"
#include "apply_bulk_modifications.hpp"
#include "is_single_stranded.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// Find all of the nodes with no edges on their left sides.
vector<handle_t> head_nodes(const HandleGraph* g);

/// Find all of the nodes with no edges on their right sides.
vector<handle_t> tail_nodes(const HandleGraph* g);

/**
 * Order and orient the nodes in the graph using a topological sort. The sort is
 * guaranteed to be machine-independent given the initial graph's node and edge
 * ordering. The algorithm is well-defined on non-DAG graphs, but the order is
 * necessarily not a topological order.
 * 
 * We use a bidirected adaptation of Kahn's topological sort (1962), which can handle components with no heads or tails.
 * 
 * L ← Empty list that will contain the sorted and oriented elements
 * S ← Set of nodes which have been oriented, but which have not had their downstream edges examined
 * N ← Set of all nodes that have not yet been put into S
 * 
 * while N is nonempty do
 *     remove a node from N, orient it arbitrarily, and add it to S
 *         (In practice, we use "seeds": the heads all in a batch at the start, and any
 *          nodes we have seen that had too many incoming edges)
 *     while S is non-empty do
 *         remove an oriented node n from S
 *         add n to tail of L
 *         for each node m with an edge e from n to m do
 *             remove edge e from the graph
 *             if m has no other edges to that side then
 *                 orient m such that the side the edge comes to is first
 *                 remove m from N
 *                 insert m into S
 *             otherwise
 *                 put an oriented m on the list of arbitrary places to start when S is empty
 *                     (This helps start at natural entry points to cycles)
 *     return L (a topologically sorted order and orientation)
 */
vector<handle_t> topological_order(const HandleGraph* g);

/**
 * Order the nodes in a graph using a topological sort. The sort is NOT guaranteed
 * to be machine-independent, but it is faster than topological_order(). This algorithm 
 * is invalid in a graph that has any cycles. For safety, consider this property with
 * algorithms::is_directed_acyclic().
 */
vector<handle_t> lazy_topological_order(const HandleGraph* g);
    
/**
 * Order the nodes in a graph using a topological sort. Similar to lazy_topological_order
 * but somewhat faster. The algorithm is invalid in a graph that has any cycles or
 * any reversing edges. For safety, consider these properties with algorithms::is_acyclic()
 * and algorithms::is_single_stranded().
 */
vector<handle_t> lazier_topological_order(const HandleGraph* g);
    
/**
 * Topologically sort the given handle graph, and then apply that sort to re-
 * order the nodes of the graph. The sort is guaranteed to be stable. This sort is well-defined
 * on graphs that are not DAGs, but instead of finding a topological sort ti does a heuristic
 * sort to minimize a feedback arc set.
 */
void topological_sort(MutableHandleGraph* g);
    
/**
 * Topologically sort the given handle graph, and then apply that sort to re-
 * order the nodes of the graph. The sort is NOT guaranteed to be stable or
 * machine-independent, but it is faster than sort(). This algorithm is invalid 
 * in a graph that has any cycles. For safety, consider checking this property with
 * algorithms::is_acyclic().
 */
void lazy_topological_sort(MutableHandleGraph* g);

/**
 * Topologically sort the given handle graph, and then apply that sort to re-
 * order the nodes of the graph. The sort is NOT guaranteed to be stable or
 * machine-independent, but it is faster than sort() and somewhat faster than
 * lazy_sort(). This algorithm is invalid in a graph that has any cycles or reversing
 * edges. For safety, consider checking these properties with algorithms::is_single_stranded()
 * and algorithms::is_acyclic().
 */
void lazier_topological_sort(MutableHandleGraph* g);
    
/**
 * Topologically sort the given handle graph, and then apply that sort to orient
 * all the nodes in the global forward direction. May invalidate any paths
 * stored by the graph. The re-orientation is guaranteed to be stable.
 * Invalidates all handles into the graph (since any node might be flipped).
 */
unordered_set<id_t> orient_nodes_forward(MutableHandleGraph* g);
                                                      
}
}

#endif
