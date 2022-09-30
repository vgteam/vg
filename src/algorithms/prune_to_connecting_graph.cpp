#include "algorithms/prune_to_connecting_graph.hpp"

#include <unordered_set>
#include <queue>

//#define debug

namespace vg {
namespace algorithms {

void prune_to_connecting_graph(DeletableHandleGraph& graph,
                               const handle_t& from, const handle_t& to) {
    
    // we use these to remember which hanndles were reached
    unordered_set<handle_t> forward, backward;
    
    // do BFS from both positions
    for (bool fwd : {true, false}) {
        
        auto& reached = fwd ? forward : backward;
        handle_t start = fwd ? from : to;
        
        queue<handle_t> bfs_queue;
        reached.emplace(start);
        bfs_queue.emplace(start);
                
        while (!bfs_queue.empty()) {
            
            handle_t here = bfs_queue.front();
            bfs_queue.pop();
            
#ifdef debug
            cerr << "BFS in direction fwd? " << fwd << " at " << graph.get_id(here) << " " << graph.get_is_reverse(here) << endl;
#endif
            
            graph.follow_edges(here, !fwd, [&](const handle_t& next) {
#ifdef debug
                cerr << "follow edge to " << graph.get_id(next) << " " << graph.get_is_reverse(next) << endl;
#endif
                if (!reached.count(next)) {
                    reached.emplace(next);
                    bfs_queue.emplace(next);
#ifdef debug
                    cerr << "add to queue" << endl;
#endif
                }
            });
        }
    }
    
    // check that each node is on some path between the two handles
    vector<handle_t> to_delete;
    graph.for_each_handle([&](const handle_t& handle) {
        auto flipped = graph.flip(handle);
        if (!((forward.count(handle) && backward.count(handle))
              || (forward.count(flipped) && backward.count(flipped)))) {
#ifdef debug
            cerr << "mark " << graph.get_id(handle) << " for deletion" << endl;
#endif
            to_delete.push_back(handle);
        }
    });
    
    // delete the handles that failed the test
    for (auto handle : to_delete) {
        graph.destroy_handle(handle);
    }
}

}
}
