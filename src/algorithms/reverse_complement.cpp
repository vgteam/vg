#include "reverse_complement.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    unordered_map<id_t, pair<id_t, bool>> reverse_complement_graph(const HandleGraph* source,
                                                                   MutableHandleGraph* into) {
        
        if (into->get_node_count()) {
            cerr << "error:[algorithms] attempted to create reversed graph in a non-empty graph" << endl;
            exit(1);
        }
        
        // the return value, translation from 'into' -> 'source'
        unordered_map<id_t, pair<id_t, bool>> node_translation;
        // for translating b/w the graphs the other direction in the course in the algorithm
        unordered_map<handle_t, handle_t> forward_translation;
        
        node_translation.reserve(source->get_node_count());
        forward_translation.reserve(source->get_node_count());
        
        // make the nodes in reverse orientation
        source->for_each_handle([&](const handle_t& handle) {
            handle_t rev_handle = into->create_handle(source->get_sequence(source->flip(handle)));
            node_translation[into->get_id(rev_handle)] = make_pair(source->get_id(handle), true);
            forward_translation[handle] = rev_handle;
        });
        
        // make the edges
        source->for_each_edge([&](const edge_t& edge) {
            // get the two sides in the correct orientation
            handle_t rev_left, rev_right;
            if (source->get_is_reverse(edge.first)) {
                rev_left = into->flip(forward_translation[source->flip(edge.first)]);
            }
            else {
                rev_left = forward_translation[edge.first];
            }
            
            if (source->get_is_reverse(edge.second)) {
                rev_right = into->flip(forward_translation[source->flip(edge.second)]);
            }
            else {
                rev_right = forward_translation[edge.second];
            }
            
            // actually make the edge
            into->create_edge(rev_right, rev_left);
            
            // always keep going
            return true;
        });
        
        return move(node_translation);
    }
}
}
