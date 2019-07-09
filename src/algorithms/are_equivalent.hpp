#ifndef VG_ALGORITHMS_ARE_EQUIVALENT_HPP_INCLUDED
#define VG_ALGORITHMS_ARE_EQUIVALENT_HPP_INCLUDED

/**
 * \file are_equivalent.hpp
 *
 * Defines an algorithm for determining whether two graphs are identical
 */

#include <iostream>

#include "../handle.hpp"

namespace vg {
namespace algorithms {
using namespace std;

    /// Checks whether two graphs are identical in their IDs, sequences, and edges.
    /// Optionally reports why graphs are found non-equivalent to stderr.
    bool are_equivalent(const HandleGraph* graph_1,
                        const HandleGraph* graph_2, bool verbose = false);

    /// Checks whether two graphs are identical in their IDs, sequences, edges, and paths.
    /// Optionally reports why graphs are found non-equivalent to stderr.
    bool are_equivalent_with_paths(const PathHandleGraph* graph_1,
                                   const PathHandleGraph* graph_2, bool verbose = false);
}
}

#endif
