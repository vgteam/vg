#include "apply_bulk_modifications.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    unordered_set<id_t> apply_orientations(MutableHandleGraph* graph, const vector<handle_t>& orientations) {
        
        // Track what we flip
        unordered_set<id_t> flipped;
        for (const auto& handle : orientations) {
            if (graph->get_is_reverse(handle)) {
                // This needs to be flipped
                flipped.insert(graph->get_id(handle));
                // Flip it
                graph->apply_orientation(handle);
            }
        }
        return flipped;
    }
    
    void apply_ordering(MutableHandleGraph* graph, const vector<handle_t>& ordering) {
        
        if (graph->node_size() != ordering.size()) {
            cerr << "error:[algorithms] attempting to sort a graph with an incomplete ordering" << endl;
            exit(1);
        }
        
        // TODO: we don't check that all nodes are present only once, which might be nice to do
        
        size_t index = 0;
        graph->for_each_handle([&](const handle_t& at_index) {
            // For each handle in the graph, along with its index
            
            // Swap the handle we observe at this index with the handle that we know belongs at this index.
            // The loop invariant is that all the handles before index are the correct sorted handles in the right order.
            // Note that this ignores orientation
            graph->swap_handles(at_index, ordering.at(index));
            
            // Now we've written the sorted handles through one more space.
            index++;
        });
    }
}
}
