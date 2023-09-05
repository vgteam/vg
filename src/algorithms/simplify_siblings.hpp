#ifndef VG_ALGORITHMS_SIMPLIFY_SIBLINGS_HPP_INCLUDED
#define VG_ALGORITHMS_SIMPLIFY_SIBLINGS_HPP_INCLUDED


#include "../handle.hpp"
#include "merge.hpp"
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

namespace vg {
namespace algorithms {

using namespace std;

/**
 * Simplify siblings in the given graph.
 *
 * When one base has two successors with the same base value, and those
 * successors have the same set of predecessors, the successors will be merged.
 *
 * Performs only a subset of the possible merges. Can only merge in from one
 * side of a given node in a single invocation. Returns true if it made
 * progress and there may be more merging to do.
 *
 * Preserves paths.
 *
 * Optional can_merge callback will only let nodes get merged together if 
 * this pairwise check returns true. 
 */
bool simplify_siblings(handlegraph::MutablePathDeletableHandleGraph* graph,
                       function<bool(const handle_t&, const handle_t&)> can_merge = nullptr);
    
}
}

#endif
