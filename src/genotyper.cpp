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
        auto& endmappings_by_name = graph.paths.get_node_mapping(end.node);
        
        for(auto& name_and_mappings : graph.paths.get_node_mapping(start.node)) {
            // Go through the paths that visit the start node
            
            if(!endmappings_by_name.count(name_and_mappings.first)) {
                // No path by this name has any mappings to the end node. Skip
                // it early.
                continue;
            }
            
            for(auto* mapping : name_and_mappings.second) {
                // Start at each mapping in the appropriate orientation
                
                // How many times have we gone to the next mapping looking for a
                // mapping to the end node in the right orientation?
                size_t traversal_count = 0;
                
                // Do we want to go left (true) or right (false) from this
                // mapping? If start is a forward traversal and we found a
                // forward mapping, we go right. If either is backward we go
                // left, and if both are backward we go right again.
                bool traversal_direction = mapping->position().is_reverse() != start.backward;
                
                // What orientation would we want to find the end node in? If
                // we're traveling backward, we expect to find it in the
                // opposite direction to the one we were given.
                bool expected_end_orientation = end.backward != traversal_direction;
                
                // We're going to fill in this list with traversals.
                list<NodeTraversal> path_traversed;
                
                while(mapping != nullptr && traversal_count < max_path_search_steps) {
                    // Traverse along until we hit the end traversal or take too
                    // many steps
                    
                    // Say we visit this node along the path, in this orientation
                    path_traversed.push_back(NodeTraversal(graph.get_node(mapping->position().node_id()), mapping->position().is_reverse() != traversal_direction));
                    
                    if(mapping->position().node_id() == end.node->id() && mapping->position().is_reverse() == expected_end_orientation) {
                        // We have stumbled upon the end node in the orientation we wanted it in.
                        // Keep this subpath if it's new
                        results.emplace(std::move(path_traversed));
                        // Then try the next embedded path
                        break;
                    }
                    
                    // Otherwise just move to the right (or left)
                    if(traversal_direction) {
                        // We're going backwards
                        mapping = graph.paths.traverse_left(mapping);
                    } else {
                        // We're going forwards
                        mapping = graph.paths.traverse_right(mapping);
                    }
                    // Tick the counter so we don't go really far on long paths.
                    traversal_count++;
                    
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

map<Alignment*, vector<int>> Genotyper::get_affinities(VG& graph, const map<string, Alignment*>& reads_by_name,
    vector<id_t>& superbubble_contents, vector<Path>& superbubble_paths) {
    
    // We're going to build this up gradually, appending to all the vectors.
    map<Alignment*, vector<int>> toReturn;
    
    // What reads are relevant to this superbubble?
    set<string> relevant_read_names;
    
    for(auto id : superbubble_contents) {
        // For every node in the superbubble, what paths visit it?
        auto& mappings_by_name = graph.paths.get_node_mapping(id);
        for(auto& name_and_mappings : mappings_by_name) {
            // For each path visiting the node
            if(reads_by_name.count(name_and_mappings.first)) {
                // This path is a read, so add the name to the set if it's not
                // there already
                relevant_read_names.insert(name_and_mappings.first);
            }    
        }
    }
    
    // What IDs are visited by these reads?
    set<id_t> relevant_ids;
    
    for(auto& name : relevant_read_names) {
        // Get the mappings for each read
        auto& mappings = graph.paths.get_path(name);
        for(auto& mapping : mappings) {
            // Add in all the nodes that are visited
            relevant_ids.insert(mapping.position().node_id());
        }
    }
    
    for(auto id : superbubble_contents) {
        // Throw out all the IDs that are also used in the superbubble itself
        relevant_ids.erase(id);
    }
    
    // Make a vg graph with all the nodes used by the reads relevant to the
    // superbubble, but outside the superbubble itself.
    VG surrounding;
    for(auto id : relevant_ids) {
        // Add each node and its edges to the new graph. Ignore dangling edges.
        surrounding.add_node(*graph.get_node(id));
        surrounding.add_edges(graph.edges_of(graph.get_node(id)));
    }
    
    for(auto& path : superbubble_paths) {
        // Now for each superbubble path, make a copy of that graph with it in
        VG allele_graph(surrounding);
        
        for(size_t i = 0; i < path.mapping_size(); i++) {
            // Add in every node on the path
            id_t id = path.mapping(i).position().node_id();
            surrounding.add_node(*graph.get_node(id));
            // And its edges
            surrounding.add_edges(graph.edges_of(graph.get_node(id)));
        }
        
        // Get rid of dangling edges
        allele_graph.remove_orphan_edges();
        
        
        for(auto& name : relevant_read_names) {
            // For every read that touched the superbubble, grab its original
            // Alignment pointer.
            Alignment* read = reads_by_name.at(name);
            
            // Re-align a copy to this graph (using quality-adjusted alignment)
            Alignment aligned = allele_graph.align_qual_adjusted(*read, aligner);
            
            // Grab the score and save it for this read and superbubble path
            toReturn[read].push_back(aligned.score());
            
        }
    }
    
    // After scoring all the reads against all the versions of the superbubble,
    // return the affinities
    return toReturn;
    
    

}

}


