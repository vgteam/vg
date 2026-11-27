#ifndef VG_ALGORITHMS_FUSE_NODES_HPP_INCLUDED
#define VG_ALGORITHMS_FUSE_NODES_HPP_INCLUDED

#include "../handle.hpp"
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

namespace vg {
namespace algorithms {

/// Check whether fuse_nodes() can join left and right without dropping an edge.
bool can_fuse_nodes(const HandleGraph& graph, handle_t left, handle_t right);

/// Replace left and right, which must pass can_fuse_nodes(), with one node spelling both; paths on either now visit it.
handle_t fuse_nodes(handlegraph::MutablePathDeletableHandleGraph* graph,
                    handle_t left, handle_t right);

}
}

#endif
