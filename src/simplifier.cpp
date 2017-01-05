#include "simplifier.hpp"

namespace vg {

using namespace std;

Simplifier::Simplifier(VG& graph) : Progressive(),
    graph(graph), site_finder(graph, ""), traversal_finder(graph) {
    // Nothing to do
}

pair<size_t, size_t> Simplifier::simplify_once(size_t iteration) {

    // Set up the deleted node and edge counts
    pair<size_t, size_t> to_return {0, 0};
    auto& deleted_nodes = to_return.first;
    auto& deleted_edges = to_return.second;

    if(!graph.is_valid(true, true, true, true)) {
        // Make sure the graph is valid and not missing nodes or edges
        cerr << "error:[vg::Simplifier] Invalid graph on iteration " << iteration << endl;
        exit(1);
    }

    // Make a list of leaf sites
    list<NestedSite> leaves;
    
    if (show_progress) {
        cerr << "Iteration " << iteration << ": Scanning " << graph.node_count() << " nodes and "
            << graph.edge_count() << " edges for sites..." << endl;
    }
    
    site_finder.for_each_site_parallel([&](NestedSite root) {
        // For every tree of sites
        
        // We keep a queue of sites to process.
        list<NestedSite*> to_check{&root};
        
#ifdef debug
        #pragma omp critical (cerr)
        cerr << "Found root site with " << root.children.size() << " children" << endl;
#endif

        while (!to_check.empty()) {
            // Until we have seen all the sites, check one
            NestedSite& check = *(to_check.front());
            to_check.pop_front();
            
            // Drop any children that have no contents, and subsume them
            // into their parent.
            vector<NestedSite> nonempty_children;
            
            for (auto& child : check.children) {
                if (child.nodes.size() == 2) {
                    // The child hs only a start and an end node, with no
                    // content nodes. It's not legit.
                    
                    for (Node* node : child.nodes) {
                        // Copy up all the nodes
                        check.nodes.insert(node);
                    }
                    
                    for (Edge* edge : child.edges) {
                        // Copy up all the edges
                        check.edges.insert(edge);
                    }
                    
                    // Don't preserve the child
                    
                } else {
                    // The child is legit. Keep it.
                    nonempty_children.emplace_back(std::move(child));
                }
            }

            // Swap the old child list with the new nonempty child list
            swap(check.children, nonempty_children);
            
            if (check.children.empty()) {
                // This is a legitimate leaf. Move it into the leaves queue.
                #pragma omp critical (leaves)
                leaves.emplace_back(std::move(check));
            } else {
                for (auto& child : check.children) {
                    // Check all the children, depth-first.
                    to_check.push_front(&child);
                }
            }
        }
        
    });
    
    if (show_progress) {
        cerr << "Found " << leaves.size() << " leaves" << endl;
    }
    
    // Index all the graph paths
    map<string, unique_ptr<PathIndex>> path_indexes;
    graph.paths.for_each_name([&](const string& name) {
        // For every path name, go index it and put it in this collection
        path_indexes.insert(make_pair(name, move(unique_ptr<PathIndex>(new PathIndex(graph, name)))));
    });
    
    // Now we have a list of all the leaf sites.
    create_progress("simplifying leaves", leaves.size());
    
    for (auto& leaf : leaves) {
        // Look at all the leaves
        
        // For each leaf, calculate its total size.
        size_t total_size = 0;
        for (auto* node : leaf.nodes) {
            // For each node
            if (node == leaf.start.node || node == leaf.end.node) {
                // That isn't a start or end
                continue;
            }
            // Include it in the size figure
            total_size += node->sequence().size();
        }
        
        if (total_size == 0) {
            // This site is just the start and end nodes, so it doesn't make
            // sense to try and remove it.
            continue;
        }
        
        if (total_size >= min_size) {
            // This site is too big to remove
            continue;
        }
        
#ifdef debug
        cerr << "Found " << total_size << " bp leaf" << endl;
        for (auto* node : leaf.nodes) {
            cerr << "\t" << node->id() << ": " << node->sequence() << endl;
        }
#endif
        
        // Otherwise we want to simplify this site away
        
        // Identify the replacement traversal for the bubble
        auto traversals = traversal_finder.find_traversals(leaf);
        
        if (traversals.empty()) {
            // We couldn't find any paths through the site.
            continue;
        }
        
        // Get the collection of visits in the traversal we want to keep.
        // Copy them into a vector for random access.
        vector<SiteTraversal::Visit> visits{traversals.front().visits.begin(), traversals.front().visits.end()};
        
        // Determine the length of the new traversal
        size_t new_site_length = 0;
        for (size_t i = 1; i + 1 < visits.size(); i++) {
            // For every non-anchoring node
            auto& visit = visits[i];
            // Total up the lengths of all the nodes that are newly visited.
            assert(visit.node != nullptr);
            new_site_length += visit.node->sequence().size();
        }

#ifdef debug
        cerr << "Chosen traversal is " << new_site_length << " bp" << endl;
#endif
        
        // Now we have to rewrite paths that visit nodes/edges not on this
        // traversal, or in a different order, or whatever. To be safe we'll
        // just rewrite all paths.
        
        // Find all the paths that traverse this region.
        
        // We start at the start node. Copy out all the mapping pointers on that
        // node, so we can go through them while tampering with them.
        map<string, set<Mapping*> > mappings_by_path = graph.paths.get_node_mapping(leaf.start.node);
        
        // We'll keep a set of the end mappings we managed to find
        set<Mapping*> found_end_mappings;
        
        for (auto& kv : mappings_by_path) {
            // For each path that hits the start node
            
            // Unpack the name
            auto& path_name = kv.first;
            
            // If a path can't be represented after a bubble is popped
            // (because the path reversed and came out the same side as it
            // went in), we just clobber the path entirely. TODO: handle
            // out-the-same-side traversals as valid genotypes somehow..
            bool kill_path = false;
            
            for (Mapping* start_mapping : kv.second) {
                // For each visit to the start node
                
                // Determine what orientation we're going to scan in
                bool backward = start_mapping->position().is_reverse();
                
                // We're going to fill this list with the mappings we need to
                // remove and replace in this path for this traversal. Initially
                // runs from start of site to end of site, but later gets
                // flipped into path-local orientation.
                list<Mapping*> existing_mappings;
                
                // Tracing along forward/backward from each as appropriate, see
                // if the end of the site is found in the expected orientation
                // (or if the path ends first).
                bool found_end = false;
                Mapping* here = start_mapping;
                
                // We want to remember the end mapping when we find it
                Mapping* end_mapping = nullptr;
                
#ifdef debug
                cerr << "Scanning " << path_name << " from " << pb2json(*here)
                    << " for " << leaf.end << " orientation " << backward << endl;
#endif
                
                while (here) {
                    // Until we hit the start/end of the path or the mapping we want
                    
#ifdef debug
                    cerr << "\tat " << pb2json(*here) << endl;
#endif
                    
                    if (here->position().node_id() == leaf.end.node->id() &&
                        here->position().is_reverse() == (leaf.end.backward != backward)) {
                        // We have encountered the end of the site in the
                        // orientation we expect, given the orientation we saw
                        // for the start.
                        
                        found_end = true;
                        end_mapping = here;
                        
                        // Know we got to this mapping at the end from the
                        // start, so we don't need to clobber everything
                        // before it.
                        found_end_mappings.insert(here);
                        
                        // Stop scanning!
                        break;
                    }
                    
                    if (here->position().node_id() == leaf.start.node->id() &&
                        here->position().is_reverse() != (leaf.start.backward != backward)) {
                        // We have encountered the start node with an incorrect orientation.
                        cerr << "warning:[vg simplify] Path " << path_name
                            << " doubles back through start of site "
                            << leaf.start << " - " << leaf.end << "; dropping!" << endl;
                            
                        kill_path = true;
                        break;
                    }
                    
                    if (!leaf.nodes.count(graph.get_node(here->position().node_id()))) {
                        // We really should stay inside the site!
                        cerr << "error:[vg simplify] Path " << path_name
                            << " somehow escapes site " << leaf.start << " - " << leaf.end << endl;
                            
                        exit(1);
                    }
                    
                    if (here != start_mapping) {
                        // Remember the mappings that aren't to the start or
                        // end of the site, so we can remove them later.
                        existing_mappings.push_back(here);
                    }
                    
                    // Scan left along ther path if we found the site start backwards, and right if we found it forwards.
                    Mapping* next = backward ? graph.paths.traverse_left(here) : graph.paths.traverse_right(here);
                    
                    // Make into NodeTraversals
                    NodeTraversal here_traversal(graph.get_node(here->position().node_id()), here->position().is_reverse());
                    NodeTraversal next_traversal(graph.get_node(next->position().node_id()), next->position().is_reverse());
                    
                    if (backward) {
                        // We're scanning the other way
                        swap(here_traversal, next_traversal);
                    }
                    
                    // Make sure we have an edge so we can traverse this node and then the node we're going to.
                    if(graph.get_edge(here_traversal, next_traversal) == nullptr) {
                        cerr << "error:[vg::Simplifier] No edge " << here_traversal << " to " << next_traversal << endl;
                        exit(1);
                    }
                    
                    here = next;
                }
                
                if (kill_path) {
                    // This path can't exist after we pop this bubble.
                    break;
                }
                
                
                if (!found_end) {
                    // This path only partly traverses the site, and is
                    // anchored at the start. Remove the part inside the site.
                    
                    // TODO: let it stay if it matches the one true traversal.
                                 
                    for(auto* mapping : existing_mappings) {
                        // Trim the path out of the site
                        graph.paths.remove_mapping(mapping);
                    }
                    
                    // TODO: update feature positions if we trim off the start of a path

                    // Maybe the next time the path visits the site it will go
                    // all the way through.
                    continue;
                }
                
                // If we found the end, remove all the mappings encountered, in
                // order so that the last one removed is the last one along the
                // path.
                if (backward) {
                    // Make sure the last mapping in the list is the last
                    // mapping to occur along the path.
                    existing_mappings.reverse();
                }
                
                // Where does the variable region of the site start for this
                // traversal of the path? If there are no existing mappings,
                // it's the start mapping's position if we traverse the site
                // backwards and the end mapping's position if we traverse
                // the site forwards. If there are existing mappings, it's
                // the first existing mapping's position in the path. TODO:
                // This is super ugly. Can we view the site in path
                // coordinates or something?
                PathIndex& path_index = *path_indexes.at(path_name).get();
                Mapping* mapping_after_first = existing_mappings.empty() ?
                    (backward ? start_mapping : end_mapping) : existing_mappings.front();
                assert(path_index.mapping_positions.count(mapping_after_first));
                size_t variable_start = path_index.mapping_positions.at(mapping_after_first); 
                
                
                
                // Determine the total length of the old traversal of the site
                size_t old_site_length = 0;
                for (auto* mapping : existing_mappings) {
                    // Add in the lengths of all the mappings that will get
                    // removed.
                    old_site_length += mapping_from_length(*mapping);
                }
#ifdef debug
                cerr << "Replacing " << old_site_length << " bp at " << variable_start
                    << " with " << new_site_length << " bp" << endl;
#endif

                // Actually update any BED features
                features.on_path_edit(path_name, variable_start, old_site_length, new_site_length);
                
                // Where will we insert the new site traversal into the path?
                list<Mapping>::iterator insert_position;
                
                if (!existing_mappings.empty()) {
                    // If there are existing internal mappings, we'll insert right where they were
                    
                    for (auto* mapping : existing_mappings) {
                        // Remove each mapping from left to right along the
                        // path, saving the position after the mapping we just
                        // removed. At the end we'll have the position of the
                        // mapping to the end of the site.
                        
#ifdef debug
                        cerr << path_name << ": Drop mapping " << pb2json(*mapping) << endl;
#endif
                        
                        insert_position = graph.paths.remove_mapping(mapping);
                    }
                } else {
                    // Otherwise we'll insert right before the mapping to
                    // the start or end of the site (whichever occurs last
                    // along the path)
                    insert_position = graph.paths.find_mapping(backward ? start_mapping : here);
                }
                
                // Make sure we're going to insert starting from the correct end of the site.
                if (backward) {
                    assert(insert_position->position().node_id() == leaf.start.node->id());
                } else {
                    assert(insert_position->position().node_id() == leaf.end.node->id());
                }
                
                // Make sure we have a place for internal visits
                assert(visits.size() >= 2);
                
                // Loop through the internal visits in the canonical
                // traversal backwards along the path we are splicing. If
                // it's a forward path this is just right to left, but if
                // it's a reverse path it has to be left to right.
                for (size_t i = 0; i < visits.size() - 2; i++) {
                    // Find the visit we need next, as a function of which
                    // way we need to insert this run of visits. Normally we
                    // go through the visits right to left, but when we have
                    // a backward path we go left to right.
                    auto& visit = backward ? visits[i + 1] : visits[visits.size() - 2 - i];
                    
                    // Make a Mapping to represent it
                    Mapping new_mapping;
                    new_mapping.mutable_position()->set_node_id(visit.node->id());
                    // We hit this node backward if it's backward along the
                    // traversal, xor if we are traversing the traversal
                    // backward
                    new_mapping.mutable_position()->set_is_reverse(visit.backward != backward);
                    
                    // Add an edit
                    Edit* edit = new_mapping.add_edit();
                    edit->set_from_length(visit.node->sequence().size());
                    edit->set_to_length(visit.node->sequence().size());
                    
#ifdef debug
                    cerr << path_name << ": Add mapping " << pb2json(new_mapping) << endl;
#endif
                    
                    // Insert the mapping in the path, moving right to left
                    insert_position = graph.paths.insert_mapping(insert_position, path_name, new_mapping);
                    
                }
                
                // Now we've corrected this site on this path. Update its index.
                // TODO: right now this means retracing the entire path.
                path_indexes[path_name].get()->update_mapping_positions(graph, path_name);
            }
            
            if (kill_path) {
                // Destroy the path completely, because it needs to reverse
                // inside a site that we have popped.
                graph.paths.remove_path(path_name);
            }
            
        }
        
        
        // It's possible a path can enter the site through the end node and
        // never hit the start. So we're going to trim those back before we delete nodes and edges.
        map<string, set<Mapping*> > end_mappings_by_path = graph.paths.get_node_mapping(leaf.end.node);
        
        for (auto& kv : end_mappings_by_path) {
            // Unpack the name
            auto& path_name = kv.first;
            
            // We might have to kill the path, if it reverses inside a
            // bubble we're popping
            bool kill_path = false;
            
            for (Mapping* end_mapping : kv.second) {
                if (found_end_mappings.count(end_mapping)) {
                    // Skip the traversals of the site that we handled.
                    continue;
                }
                
                // Now we're left with paths that leave the site but don't
                // enter. We're going to clobber everything before the path
                // leaves the site.
                
                // Determine what orientation we're going to scan in
                bool backward = end_mapping->position().is_reverse();
                
                // Start at the end
                Mapping* here = end_mapping;
                
                // Keep a list of mappings we need to remove
                list<Mapping*> to_remove;
                
                while (here) {
                    
                    if (here->position().node_id() == leaf.end.node->id() &&
                        here->position().is_reverse() != (leaf.end.backward != backward)) {
                        // We have encountered the end node with an incorrect orientation.
                        cerr << "warning:[vg simplify] Path " << path_name
                            << " doubles back through end of site "
                            << leaf.start << " - " << leaf.end << "; dropping!" << endl;
                            
                        kill_path = true;
                        break;
                    }
                    
                    // Say we should remove the mapping.
                    to_remove.push_back(here);
                    
                    // Scan right along the path if we found the site end backwards, and left if we found it forwards.
                    here = backward ? graph.paths.traverse_right(here) : graph.paths.traverse_left(here);
                    
                    // Eventually we should hit the end of the path, or the
                    // end of the site, since we don't hit the start.
                    
                }
                
                if (kill_path) {
                    // Just go kill the whole path
                    break;
                }
                
                for (auto* mapping: to_remove) {
                    // Get rid of all the mappings once we're done tracing them out.
                    graph.paths.remove_mapping(mapping);
                }
                
            }
            
            if (kill_path) {
                // Destroy the path completely, because it needs to reverse
                // inside a site that we have popped.
                graph.paths.remove_path(path_name);
            }
        }
        
        // Now delete all edges that aren't connecting adjacent nodes on the
        // blessed traversal (before we delete their nodes).
        set<Edge*> blessed_edges;
        for (auto i = visits.begin(); i != --visits.end(); ++i) {
            // For each node and the next node (which won't be the end)
            auto next = i;
            next++;
            
            assert(next != visits.end());
            
            // Find the edge between them
            NodeTraversal here(i->node, i->backward);
            NodeTraversal next_traversal(next->node, next->backward);
            Edge* edge = graph.get_edge(here, next_traversal);
            assert(edge != nullptr);
            
            // Remember we need it
            blessed_edges.insert(edge);
        }
        
        for (auto* edge : leaf.edges) {
            if (!blessed_edges.count(edge)) {
                // Get rid of all the edges not needed for the one true traversal
#ifdef debug
                cerr << leaf.start << " - " << leaf.end << ": Delete edge: " << pb2json(*edge) << endl;
#endif
                graph.destroy_edge(edge);
                deleted_edges++;
            }
        }
       
           
        // Now delete all the nodes that aren't on the blessed traversal.
        
        // What nodes are on it?
        set<Node*> blessed_nodes;
        for (auto& visit : visits) {
            blessed_nodes.insert(visit.node);
        }
        
        for (auto* node : leaf.nodes) {
            // For every node in the site
            if (!blessed_nodes.count(node)) {
                // If we don't need it for the chosen path, destroy it
#ifdef debug
                cerr << leaf.start << " - " << leaf.end << ": Delete node: " << pb2json(*node) << endl;
#endif
                // There may be paths still touching this node, if they
                // managed to get into the site without touching the start
                // node. We'll delete those paths.
                set<string> paths_to_kill;
                for (auto& kv : graph.paths.get_node_mapping(node)) {
                
                    if (mappings_by_path.count(kv.first)) {
                        // We've already actually updated this path; the
                        // node_mapping data is just out of date.
                        continue;
                    }
                
                    paths_to_kill.insert(kv.first);
                }
                for (auto& path : paths_to_kill) {
                    graph.paths.remove_path(path);
                    cerr << "warning:[vg simplify] Path " << path << " removed" << endl;
                }

                graph.destroy_node(node);
                
                deleted_nodes++;
            }
        }
        
        // OK we finished a leaf
        increment_progress();
    }
    
    destroy_progress();
    
    // Reset the ranks in the graph, since we rewrote paths
    graph.paths.clear_mapping_ranks();
    
    // Return the statistics.
    return to_return;

}

void Simplifier::simplify() {
    for (size_t i = 0; i < max_iterations; i++) {
        // Try up to the max number of iterations
        auto deleted_elements = simplify_once(i);
        
        if (show_progress) {
            cerr << "Iteration " << i << ": deleted " << deleted_elements.first
                << " nodes and " << deleted_elements.second << " edges" << endl;
        }
        
        if (deleted_elements.first == 0 && deleted_elements.second == 0) {
            // If nothing gets deleted, stop because trying again won't change
            // things
            break;
        }
    }
}

}
