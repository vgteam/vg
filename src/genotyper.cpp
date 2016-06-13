#include "genotyper.hpp"

namespace vg {

using namespace std;

vector<Path> Genotyper::get_paths_through_site(VG& graph, NodeTraversal start, NodeTraversal end) {
    // We're going to emit traversals supported by any paths in the graph.
    
    // Put all our subpaths in here to deduplicate them
    set<list<NodeTraversal>> results;
    
    if(graph.paths.has_node_mapping(start.node) && graph.paths.has_node_mapping(end.node)) {
        // If we have some paths that visit both ends (in some orientation)
        
        // Get all the mappings to the end node, by path name
        auto& endMappingsByName = graph.paths.get_node_mapping(end.node);
        
        for(auto& nameAndMappings : graph.paths.get_node_mapping(start.node)) {
            // Go through the paths that visit the start node
            
            if(!endMappingsByName.count(nameAndMappings.first)) {
                // No path by this name has any mappings to the end node. Skip
                // it early.
                continue;
            }
            
            for(auto* mapping : nameAndMappings.second) {
                // Start at each mapping in the appropriate orientation
                
                // How many times have we gone to the next mapping looking for a
                // mapping to the end node in the right orientation?
                size_t traversalCount = 0;
                
                // Do we want to go left (true) or right (false) from this
                // mapping? If start is a forward traversal and we found a
                // forward mapping, we go right. If either is backward we go
                // left, and if both are backward we go right again.
                bool traversalDirection = mapping->position().is_reverse() != start.backward;
                
                // What orientation would we want to find the end node in? If
                // we're traveling backward, we expect to find it in the
                // opposite direction to the one we were given.
                bool expectedEndOrientation = end.backward != traversalDirection;
                
                // We're going to fill in this list with traversals.
                list<NodeTraversal> pathTraversed;
                
                while(mapping != nullptr && traversalCount < max_path_search_steps) {
                    // Traverse along until we hit the end traversal or take too
                    // many steps
                    
                    // Say we visit this node along the path, in this orientation
                    pathTraversed.push_back(NodeTraversal(graph.get_node(mapping->position().node_id()), mapping->position().is_reverse() != traversalDirection));
                    
                    if(mapping->position().node_id() == end.node->id() && mapping->position().is_reverse() == expectedEndOrientation) {
                        // We have stumbled upon the end node in the orientation we wanted it in.
                        // Keep this subpath if it's new
                        results.emplace(std::move(pathTraversed));
                        // Then try the next embedded path
                        break;
                    }
                    
                    // Otherwise just move to the right (or left)
                    if(traversalDirection) {
                        // We're going backwards
                        mapping = graph.paths.traverse_left(mapping);
                    } else {
                        // We're going forwards
                        mapping = graph.paths.traverse_right(mapping);
                    }
                    // Tick the counter so we don't go really far on long paths.
                    traversalCount++;
                    
                }
                
                
            }
        }
        
    }
    
    // Now convert the unique results into Paths
    vector<Path> toReturn;
    
    for(auto& result : results) {
        // Convert each list of node traversals to a proper path
        toReturn.push_back(path_from_node_traversals(result));
    }
    
    return toReturn;
}

}


