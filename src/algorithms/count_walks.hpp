#ifndef VG_ALGORITHMS_COUNT_WALKS_HPP_INCLUDED
#define VG_ALGORITHMS_COUNT_WALKS_HPP_INCLUDED

/**
 * \file count_walks.hpp
 *
 * Defines algorithm for counting the number of distinct walks through a DAG.
 */

#include "../handle.hpp"
#include "topological_sort.hpp"

#include <unordered_map>
#include <vector>

namespace vg {
namespace algorithms {

using namespace std;

    /// Returns the count map through each snarl in a graph. Assumes that
    /// the graph is a single-stranded DAG. Consider checking these properties with
    /// algorithms::is_single_stranded and algorithms::is_directed_acyclic for safety.
    /// Returns numeric_limits<size_t>::max() if the actual number of walks is larger
    /// than this.
    tuple<vector<handle_t>, unordered_map<handle_t, size_t>, bool> count_walks_through_nodes(const HandleGraph* graph);

    /// Returns the number of source-to-sink walks through the graph
    size_t count_walks(const HandleGraph* graph);

}
}

#endif
