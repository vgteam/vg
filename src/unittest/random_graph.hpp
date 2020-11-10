#include "handle.hpp"
#include <vector>

#ifndef VG_UNITTEST_RANDOM_GRAPH_HPP_INCLUDED
#define VG_UNITTEST_RANDOM_GRAPH_HPP_INCLUDED

namespace vg{
namespace unittest{

/// Create a random graph by adding variation to a sequence of length seq_size
/// variant_len is the mean length of a larger variation and variant_count
/// is the number of variations added to the graph
void random_graph(int64_t seq_size, int64_t variant_len, int64_t variant_count,
                  MutablePathMutableHandleGraph* graph);

/// Create a random graph with multiple connected components. 
/// Each connected component created with a sequence length in seq_sizes
void random_graph(vector<int64_t> seq_sizes, int64_t variant_len, int64_t total_variant_count,
                  MutablePathMutableHandleGraph* graph);
                  
              
/// Create a random undirected graph over numbered nodes.
/// Returns a list of adjacency lists for each node.
std::vector<std::vector<size_t>> random_adjacency_list(size_t node_count, size_t edge_count);

}
}

#endif
