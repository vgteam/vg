#include "genotyper.hpp"
#include "bubbles.hpp"

namespace vg {

using namespace std;

map<pair<NodeTraversal, NodeTraversal>, vector<id_t>> Genotyper::find_sites(VG& graph) {

    // Set up our output map
    map<pair<NodeTraversal, NodeTraversal>, vector<id_t>> to_return;

    // Unfold the graph
    // Copy the graph and unfold the copy. We need to hold this translation
    // table from new node ID to old node and relative orientation.
    map<vg::id_t, pair<vg::id_t, bool>> unfold_translation;
    auto transformed = graph.unfold(unfold_max_length, unfold_translation);
    
    // Fix up any doubly reversed edges
    transformed.flip_doubly_reversed_edges();

    // Now dagify the graph. We need to hold this translation table from new
    // node ID to old node and relative orientation.
    map<vg::id_t, pair<vg::id_t, bool>> dag_translation;
    transformed = transformed.dagify(dagify_steps, dag_translation);
    
    // Compose the complete translation
    map<vg::id_t, pair<vg::id_t, bool>> overall_translation = transformed.overlay_node_translations(dag_translation, unfold_translation);
    dag_translation.clear();
    unfold_translation.clear();
    
    // Find the superbubbles in the DAG
    map<pair<id_t, id_t>, vector<id_t>> superbubbles = vg::superbubbles(transformed);
    
    for(auto& superbubble : superbubbles) {
        
        // Translate the superbubble coordinates into NodeTraversals
        auto& start_translation = overall_translation[superbubble.first.first];
        NodeTraversal start(graph.get_node(start_translation.first), start_translation.second);
        
        auto& end_translation = overall_translation[superbubble.first.second];
        NodeTraversal end(graph.get_node(end_translation.first), end_translation.second);
        
        // Find the vector we want all the nodes in
        auto& superbubble_nodes = to_return[make_pair(start, end)];
        
        for(auto id : superbubble.second) {
            // Translate each ID and put it in the vector
            superbubble_nodes.push_back(overall_translation[id].first);
        }
    }

    // Give back the collections of involved nodes by start and end
    return to_return;    
    
}

vector<Path> Genotyper::get_paths_through_site(VG& graph, NodeTraversal start, NodeTraversal end) {
    // We're going to emit traversals supported by any paths in the graph.
    
    // Put all our subpaths in here to deduplicate them
    set<list<NodeTraversal>> results;

#ifdef debug
    cerr << "Looking for paths between " << start << " and " << end << endl;
#endif
    
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
                
#ifdef debug
                cerr << "Trying mapping of read " << name_and_mappings.first << endl;
#endif
                
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
                    
#ifdef debug
                    cerr << "\tTraversing " << pb2json(*mapping) << endl;
#endif
                    
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

map<Alignment*, vector<double>> Genotyper::get_affinities(VG& graph, const map<string, Alignment*>& reads_by_name,
    vector<id_t>& superbubble_contents, vector<Path>& superbubble_paths) {
    
    // Grab our thread ID, which determines which aligner we get.
    int tid = omp_get_thread_num();
    
    // We're going to build this up gradually, appending to all the vectors.
    map<Alignment*, vector<double>> toReturn;
    
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
        
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "Align to " << pb2json(allele_graph.graph) << endl;
#endif
        
        for(auto& name : relevant_read_names) {
            // For every read that touched the superbubble, grab its original
            // Alignment pointer.
            Alignment* read = reads_by_name.at(name);
            
            Alignment aligned;
            if(read->sequence().size() == read->quality().size()) {
                // Re-align a copy to this graph (using quality-adjusted alignment).
                // TODO: actually use quality-adjusted alignment
                aligned = allele_graph.align(*read);
            } else {
                // If we don't have the right number of quality scores, use un-adjusted alignment instead.
                aligned = allele_graph.align(*read);
            }
            
#ifdef debug
            #pragma omp critical (cerr)
            cerr << "\t" << pb2json(aligned) << endl;
#endif
            
            // Grab the identity and save it for this read and superbubble path
            toReturn[read].push_back(aligned.identity());
            
        }
    }
    
    // After scoring all the reads against all the versions of the superbubble,
    // return the affinities
    return toReturn;
}

Genotype Genotyper::get_genotype(const vector<Path>& superbubble_paths, const map<Alignment*, vector<double>>& affinities) {
    // Compute total affinity for each path
    vector<double> total_affinity(superbubble_paths.size(), 0);
    for(auto& alignment_and_affinities : affinities) {
        for(size_t i = 0; i < alignment_and_affinities.second.size(); i++) {
            // For each affinity of an alignment to a path
            double affinity = alignment_and_affinities.second.at(i);
            // Add it in to the total
            total_affinity.at(i) += affinity;
        }
    }
    
    Genotype genotype;
    
    if(superbubble_paths.empty()) {
        // No paths! Nothing to say!
        throw runtime_error("No paths through superbubble! Can't genotype!");
    }
    
    // If there's only one path, call hom ref
    if(superbubble_paths.size() < 2) {
        *genotype.add_allele() = superbubble_paths.front();
        Support* support = genotype.add_support();
        support->set_quality(total_affinity.front());
        return genotype;
    }
    
    // Find the two best paths
    int best_allele = -1;
    int second_best_allele = -1;
    
    for(int i = 0; i < superbubble_paths.size(); i++) {
        if(best_allele == -1 || total_affinity[best_allele] < total_affinity[i]) {
            // We have a new best allele
            second_best_allele = best_allele;
            best_allele = i;
        } else if(second_best_allele == -1 || total_affinity[second_best_allele] < total_affinity[i]) {
            // We have a new second best allele
            second_best_allele = i;
        }
    }
    
    // Add the best allele with the most support
    *genotype.add_allele() = superbubble_paths[best_allele];
    Support* best_support = genotype.add_support();
    best_support->set_quality(total_affinity[best_allele]);
    
    if(total_affinity[best_allele] > total_affinity[second_best_allele] * max_het_bias) {
        // If we're too biased to one side, call hom that
        return genotype;
    }
    
    // Else call het by adding the second best allele as well
    *genotype.add_allele() = superbubble_paths[second_best_allele];
    Support* second_best_support = genotype.add_support();
    second_best_support->set_quality(total_affinity[second_best_allele]);
    
    // TODO: use more support fields by investigating the Alignment associated
    // with the affinity, as it relates to the allele's path
    
    return genotype;
}

}


