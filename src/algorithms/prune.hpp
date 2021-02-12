#ifndef VG_ALGORITHMS_PRUNE_HPP_INCLUDED
#define VG_ALGORITHMS_PRUNE_HPP_INCLUDED

#include "../handle.hpp"

namespace vg {
namespace algorithms {

/// Take all nodes that would introduce paths of > edge_max edge crossings, remove them, and link their neighbors to
/// head_node or tail_node depending on which direction the path extension was stopped.
/// For pruning graph prior to indexing with gcsa2. Returns the number of edges removed.
size_t prune_complex(DeletableHandleGraph& graph,
                   int path_length, int edge_max);

/// Wrap the graph with heads and tails (for GCSA2 indexing) and then prune as with
/// prune_complex. Returns the number of edges removed.
size_t prune_complex_with_head_tail(DeletableHandleGraph& graph,
                                  int path_length, int edge_max);

/// Remove any weakly connected components that have total sequence
/// length under the minimum size. Returns the number of nodes removed.
size_t prune_short_subgraphs(DeletableHandleGraph& graph, int min_size);

/// Remove nodes with >= max_degree total edges on each side. Note that
/// end-to-start self loops count twice. Returns the number of nodes removed.
size_t remove_high_degree_nodes(DeletableHandleGraph& graph, int max_degree);

}
}

#endif
