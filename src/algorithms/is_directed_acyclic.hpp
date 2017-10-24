#ifndef VG_ALGORITHMS_IS_DIRECTED_ACYCLIC_HPP_INCLUDED
#define VG_ALGORITHMS_IS_DIRECTED_ACYCLIC_HPP_INCLUDED

/**
 * \file is_directed_acyclic.hpp
 *
 * Defines an algorithm for deciding if a graph is directed-acyclic, and
 * contains no way to repeat a handle.
 */

#include "../handle.hpp"

#include <unordered_set>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

/// Returns true if the graph contains no way to repeat a handle.
bool is_directed_acyclic(const HandleGraph* graph);


}
}

#endif
