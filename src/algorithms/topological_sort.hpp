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

namespace vg {
namespace algorithms {

using namespace std;

/// Find all of the nodes with no edges on their left sides.
vector<handle_t> head_nodes(const HandleGraph* g);

/// Find all of the nodes with no edges on their right sides.
vector<handle_t> head_nodes(const HandleGraph* g);

/**
 * Order and orient the nodes in the graph using a topological sort.
 * 
 * We use a bidirected adaptation of Kahn's topological sort (1962), which can handle components with no heads or tails.
 * 
 * L ← Empty list that will contain the sorted and oriented elements
 * S ← Set of nodes which have been oriented, but which have not had their downstream edges examined
 * N ← Set of all nodes that have not yet been put into S
 * 
 * while N is nonempty do
 *     remove a node from N, orient it arbitrarily, and add it to S
 *         (In practice, we use "seeds": the heads and any nodes we have seen that had too many incoming edges)
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
vector<handle_t> topological_sort(const HandleGraph* g);
                                                      
}
}

#endif
