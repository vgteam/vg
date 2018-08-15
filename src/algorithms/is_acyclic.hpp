#ifndef VG_ALGORITHMS_IS_ACYCLIC_HPP_INCLUDED
#define VG_ALGORITHMS_IS_ACYCLIC_HPP_INCLUDED

/**
 * \file is_acyclic.hpp
 *
 * Defines algorithms for deciding if a graph is acyclic or directed acyclic
 */

#include "../handle.hpp"
#include "is_single_stranded.hpp"

#include <unordered_set>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

// Returns true if the graph contains no cycles (i.e. true if no node can reach itself
// along a bidirected walk).
bool is_acyclic(const HandleGraph* graph);
    
/// Returns true if the graph contains no directed cycles. It may contain reversing
/// cycles (i.e. true if no node can reach itself in the same orientation along a
/// bidirected walk, but it might be able to reach itself in the opposite orientation).
bool is_directed_acyclic(const HandleGraph* graph);

}
}

#endif
