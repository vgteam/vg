#ifndef VG_ALGORITHMS_IS_SINGLE_STRANDED_HPP_INCLUDED
#define VG_ALGORITHMS_IS_SINGLE_STRANDED_HPP_INCLUDED

/**
 * \file is_single_stranded.hpp
 *
 * Defines an algorithm for deciding if a graph contains reversing edges.
 */

#include "../handle.hpp"

#include <unordered_set>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

/// Returns true if the graph contains no reversing edges (i.e. edges that connected
/// the locally forward orientation of a node to the locally reverse orientation of
/// of another node).
bool is_single_stranded(const HandleGraph* graph);


}
}

#endif
