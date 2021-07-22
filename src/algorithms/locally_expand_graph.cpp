/**
 * \file locally_expand_graph.cpp
 *
 * Implementation for the locally_expand_graph algorithm.
 */
 
#include "locally_expand_graph.hpp"
#include <structures/rank_pairing_heap.hpp>

//#define debug_locally_expand_graph

namespace vg {
namespace algorithms {

using namespace structures;

void locally_expand_graph(const HandleGraph& parent, MutableHandleGraph& subgraph,
                          handle_t from, int64_t max_dist) {
    
    RankPairingHeap<handle_t, int64_t, greater<int64_t>> queue;
    
    handle_t start = parent.get_handle(subgraph.get_id(from), subgraph.get_is_reverse(from));
    
    queue.push_or_reprioritize(start, -parent.get_length(start));
    while (!queue.empty()) {
        handle_t handle;
        int64_t dist;
        tie(handle, dist) = queue.top();
        queue.pop();
        
#ifdef debug_locally_expand_graph
        cerr << "at " << parent.get_id(handle) << " " << parent.get_is_reverse(handle) << ", dist " << dist << endl;
#endif
        
        int64_t dist_thru = dist + parent.get_length(handle);
        if (dist_thru < max_dist) {
            parent.follow_edges(handle, false, [&](const handle_t& next) {
#ifdef debug_locally_expand_graph
                cerr << "\ttake edge to " << parent.get_id(next) << " " << parent.get_is_reverse(next) << ", dist " << dist_thru << endl;
#endif
                queue.push_or_reprioritize(next, dist_thru);
                if (!subgraph.has_node(parent.get_id(next))) {
                    subgraph.create_handle(parent.get_sequence(parent.forward(next)), parent.get_id(next));
                }
                subgraph.create_edge(subgraph.get_handle(parent.get_id(handle), parent.get_is_reverse(handle)),
                                     subgraph.get_handle(parent.get_id(next), parent.get_is_reverse(next)));
            });
        }
    }
}

}
}
