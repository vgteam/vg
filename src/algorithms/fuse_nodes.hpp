#ifndef VG_ALGORITHMS_FUSE_NODES_HPP_INCLUDED
#define VG_ALGORITHMS_FUSE_NODES_HPP_INCLUDED

#include "../handle.hpp"
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

namespace vg {
namespace algorithms {

/// Check whether two nodes can be fused.
bool can_fuse_nodes(const HandleGraph& graph, handle_t left, handle_t right);

/// Fuse two nodes into one node.
handle_t fuse_nodes(handlegraph::MutablePathDeletableHandleGraph* graph,
                    handle_t left, handle_t right);

}
}

#endif
