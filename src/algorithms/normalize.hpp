#ifndef VG_ALGORITHMS_NORMALIZE_HPP_INCLUDED
#define VG_ALGORITHMS_NORMALIZE_HPP_INCLUDED


#include "../handle.hpp"
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

/**
 * Normalize a graph, performing up to the given number of iterations.
 * Simplifies siblings and unchops runs of nodes, in a loop.
 *
 * if "can_merge" specified, it must return true in order for a pair of nodes to get merged
 */
void normalize(handlegraph::MutablePathDeletableHandleGraph* graph, int max_iter = 1,
               bool debug = false,
               function<bool(const handle_t&, const handle_t&)> can_merge = nullptr);

    
}
}

#endif
