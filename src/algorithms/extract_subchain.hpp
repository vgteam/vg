#ifndef VG_ALGORITHMS_EXTRACT_SUBCHAIN_HPP_INCLUDED
#define VG_ALGORITHMS_EXTRACT_SUBCHAIN_HPP_INCLUDED

/**
 * \file extract_subchain.hpp
 *
 * Algorithm that extracts all node ids in a subchain of the snarl decomposition.
 */

#include "../handle.hpp"
#include "../hash_map.hpp"

namespace vg {

//------------------------------------------------------------------------------

/**
 * Extracts all node ids in a subchain of the snarl decomposition.
 *
 * The subchain is defined by two handles, start and end.
 * If the corresponding nodes are on the same chain in the snarl decomposition
 * and end is after start, the function returns all node ids in the subchain
 * between them. This means the identifiers of start and end, as well as all
 * node identifiers in the weakly connected component formed by removing start
 * and end from the graph.
 *
 * If the nodes are not on the same chain or end is before start, the behavior
 * is undefined.
 */
hash_set<nid_t> extract_subchain(const HandleGraph& graph, handle_t start, handle_t end);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_ALGORITHMS_EXTRACT_SUBCHAIN_HPP_INCLUDED
