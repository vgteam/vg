#ifndef VG_ALGORITHMS_MERGE_HPP_INCLUDED
#define VG_ALGORITHMS_MERGE_HPP_INCLUDED


#include "../handle.hpp"
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

/**
 * Merge the given ranges of bases on the given handles together, rewriting paths.
 * Sequences must match. Handles to a single node may occur no more than once.
 */
void merge(handlegraph::MutablePathDeletableHandleGraph* graph, const vector<pair<handle_t, size_t>>& start, size_t length);
    
}
}

#endif
