#pragma once

#include "../handle.hpp"
//#include "../subgraph.hpp"
#include "../types.hpp"
#include "shortest_cycle.hpp"
#include <unordered_map>

namespace vg {
namespace algorithms {

using namespace std;

/// expand the subgraph iteratively for this many steps
void expand_subgraph_by_steps(const HandleGraph& source, MutableHandleGraph& subgraph, const uint64_t& steps, bool forward_only = false);

/// expand the subgraph iteratively until its node count is at least node_count
void expand_subgraph_to_node_count(const HandleGraph& source, MutableHandleGraph& subgraph, const uint64_t& node_count, bool forward_only = false);

/// expand the subgraph iteratively to include at least length new sequence
void expand_subgraph_by_length(const HandleGraph& source, MutableHandleGraph& subgraph, const uint64_t& length, bool forward_only = false);

/// expand the subgraph iterativel until its total sequence length is greater than length
void expand_subgraph_to_length(const HandleGraph& source, MutableHandleGraph& subgraph, const uint64_t& length, bool forward_only = false);

/// expand the context around a single handle position
void extract_context(const HandleGraph& source, MutableHandleGraph& subgraph, const handle_t& handle, const uint64_t& offset, const uint64_t& length, bool go_fwd = true, bool go_rev = true);

/// extract the node id range
void extract_id_range(const HandleGraph& source, const nid_t& id1, const nid_t& id2, MutableHandleGraph& subgraph);

/// extract the path range
/// nodes aren't cut, so the returned graph may start before start and/or end after end
/// if end < 0, then it will walk to the end of the path
void extract_path_range(const PathPositionHandleGraph& source, path_handle_t path_handle, int64_t start, int64_t end, MutableHandleGraph& subgraph);

/// add subpaths to the subgraph, providing a concatenation of subpaths that are discontiguous over the subgraph
/// based on their order in the path position index provided by the source graph
/// will clear any path found in both graphs before writing the new steps into it
/// if subpath_naming is true, a suffix will be added to each path in the subgraph denoting its offset
/// in the source graph (unless the subpath was not cut up at all)
void add_subpaths_to_subgraph(const PathPositionHandleGraph& source, MutablePathHandleGraph& subgraph,
                              bool subpath_naming = false);

/// We can accumulate a subgraph without accumulating all the edges between its nodes
/// this helper ensures that we get the full set
void add_connecting_edges_to_subgraph(const HandleGraph& source, MutableHandleGraph& subgraph);

}
}
