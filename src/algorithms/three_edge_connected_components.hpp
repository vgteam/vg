#ifndef VG_ALGORITHMS_THREE_EDGE_CONNECTED_COMPONENTS_HPP_INCLUDED
#define VG_ALGORITHMS_THREE_EDGE_CONNECTED_COMPONENTS_HPP_INCLUDED

#include <functional>
#include <vector>
#include <unordered_map>

namespace vg {
namespace algorithms {

using namespace std;

// Interface

/**
 * Get the three-edge-connected components of an arbitrary graph (not
 * necessarily a handle graph). Only recognizes one kind of edge and one kind
 * of node. Nodes are arbitrary value types (which may need to be hashable).
 *
 * Takes a function that loops an iteratee over all nodes, and a function that,
 * given a node, loops an iteratee over all nodes connected to it.
 *
 * For each component identified, calls the given callback with a function that
 * iterates over all nodes in the component.
 *
 * If you have a graph where you can easily rank the nodes, don't use this. Use
 * three_edge_connected_components_dense() instead. The first thing this
 * function does is asign nodes a dense, 0-based rank space.  
 */
template<typename TECCNode>
void three_edge_connected_components(const function<void(const function<void(TECCNode)>&)>& for_each_node,
    const function<void(TECCNode, const function<void(TECCNode)>&)>& for_each_connected_node,
    const function<void(const function<void(const function<void(TECCNode)>&)>&)>& component_callback);


/**
 * Get the three-edge-connected components of an arbitrary graph (not
 * necessarily a handle graph). Only recognizes one kind of edge and one kind
 * of node. Nodes are dense positive integers starting with 0.
 *
 * Takes a function that loops an iteratee over all nodes, and a function that,
 * given a node, loops an iteratee over all nodes connected to it.
 *
 * For each component identified, calls the given callback with a function that
 * iterates over all nodes in the component.
 */
void three_edge_connected_components_dense(const function<void(const function<void(size_t)>&)>& for_each_node,
    const function<void(size_t, const function<void(size_t)>&)>& for_each_connected_node,
    const function<void(const function<void(const function<void(size_t)>&)>&)>& component_callback);

// Implementation

template<typename TECCNode>
void three_edge_connected_components(const function<void(const function<void(TECCNode)>&)>& for_each_node,
    const function<void(TECCNode, const function<void(TECCNode)>&)>& for_each_connected_node,
    const function<void(const function<void(const function<void(TECCNode)>&)>&)>& component_callback) {
    // Convert to small positive integers
    
    vector<TECCNode> rank_to_node;
    unordered_map<TECCNode, size_t> node_to_rank;
    
    for_each_node([&](TECCNode node) {
        // Populate the rank/node translation.
        // TODO: can we condense this?
        node_to_rank[node] = rank_to_node.size();
        rank_to_node.push_back(node);
    });
    
    three_edge_connected_components_dense([&](const function<void(size_t)>& iteratee) {
        for (size_t i = 0; i < rank_to_node.size(); i++) {
            // All the nodes are just their ranks now
            iteratee(i);
        }
    }, [&](size_t rank, const function<void(size_t)> visit_connected) {
        // Translate the rank we are asked about into a node
        for_each_connected_node(rank_to_node[rank], [&](TECCNode connected) {
            // And translate the node back into a rank
            visit_connected(node_to_rank[connected]);
        });
    }, [&](const function<void(const function<void(size_t)>&)>& for_each_component_member) {
        // When we get a component
        // Call our component callback with a function that takes the iteratee
        component_callback([&](const function<void(TECCNode)>& iteratee) {
            for_each_component_member([&](size_t member) {
                // And for each member of the component we got, translate it and send it off.
                iteratee(rank_to_node[member]);
            });
        });
    });
}


    
}
}

#endif
