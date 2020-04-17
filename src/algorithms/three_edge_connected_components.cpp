#include "three_edge_connected_components.hpp"

// Grab the Tsin's Algorithm header out of pinchesAndCacti, which installs it into sonLib's include directory
extern "C" {
#include "sonLibList.h"
#include "sonLibTuples.h"
#include "3_Absorb3edge2x.h"
}

//#define debug

namespace vg {
namespace algorithms {

using namespace std;

void three_edge_connected_components_dense(const function<void(const function<void(size_t)>&)>& for_each_node,
    const function<void(size_t, const function<void(size_t)>&)>& for_each_connected_node,
    const function<void(const function<void(const function<void(size_t)>&)>&)>& component_callback) {
    
    // TODO: copy over or reimplement Tsin's Algorithm, or one of its successors in a way that doesn't have to copy everything.
    
    // For now we just call into the properly licensed version that's part of pinchesAndCacti.
    
    // Make the stList of all the vertices, where each vertex is an stList of single element stIntTuple items that point to the ranks of connected nodes.
    // When an item is removed, use the list destructor on it.
    stList* vertices = stList_construct3(0, (void(*)(void *)) stList_destruct);
    
    for_each_node([&](size_t rank) {
        while (rank >= stList_length(vertices)) {
            // Make sure we have an adjacency list allocated for the node
            // When an item in the node's adjacency list is destroyed, run the int tuple destructor.
            stList_append(vertices, stList_construct3(0, (void(*)(void *)) stIntTuple_destruct));
        }
        
        for_each_connected_node(rank, [&](size_t other_rank) {
            // For each edge on the node, represent it as a 1-tuple in the node's list.
            stList_append((stList*) stList_get(vertices, rank), stIntTuple_construct1((int64_t) rank));
            // We don't have to do the back-edge now; we will do it when we visit the other node.
        });
    });
    
    // Now we have the graph in the format Tsin's Algorithm wants, so run it.
    // The components come out as a list of lists, one for each component, with
    // the entries in each component's list being 1-element stIntTuples with
    // ranks in them.
    stList* components = computeThreeEdgeConnectedComponents(vertices);
    
    for(size_t i = 0; i < stList_length(components); i++) {
        // For each component
        stList* component = (stList*) stList_get(components, i);
        // Announce the component
        component_callback([&](const function<void(size_t)>& visit_member) {
            // And when we get the function to feed the members to
            for (size_t j = 0; j < stList_length(component); j++) {
                // Call it with each member
                visit_member(stIntTuple_get((stIntTuple*) stList_get(component, j), 0));
            }
        });
    }

    // Clean up the component result
    stList_destruct(components);

    // Clean up the vertex data
    stList_destruct(vertices);
}

}
}
