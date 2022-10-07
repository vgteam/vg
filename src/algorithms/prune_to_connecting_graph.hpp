#ifndef VG_ALGORITHMS_PRUNE_TO_CONNECTING_GRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_PRUNE_TO_CONNECTING_GRAPH_HPP_INCLUDED

#include "../handle.hpp"

namespace vg {
namespace algorithms {

/// Remove all parts of the graph that are not on some path between the two handles
void prune_to_connecting_graph(DeletableHandleGraph& graph,
                               const handle_t& from, const handle_t& to);
}
}

#endif
