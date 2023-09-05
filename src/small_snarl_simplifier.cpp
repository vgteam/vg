#include "small_snarl_simplifier.hpp"

#include "cactus_snarl_finder.hpp"

namespace vg {

using namespace std;

SmallSnarlSimplifier::SmallSnarlSimplifier(VG& graph) : Progressive(), graph(graph), traversal_finder(graph) {
    
    // create a SnarlManager using Cactus
    CactusSnarlFinder site_finder(graph);
    site_manager = site_finder.find_snarls();
}

pair<size_t, size_t> SmallSnarlSimplifier::simplify_once(size_t iteration) {

    // Set up the deleted node and edge counts
    pair<size_t, size_t> to_return {0, 0};
    auto& deleted_nodes = to_return.first;
    auto& deleted_edges = to_return.second;

    if(!graph.is_valid(true, true, true, true)) {
        // Make sure the graph is valid and not missing nodes or edges
        cerr << "error:[vg::SmallSnarlSimplifier] Invalid graph on iteration " << iteration << endl;
        exit(1);
    }

    // Make a list of leaf sites
    list<const Snarl*> leaves;
    
    if (show_progress) {
        cerr << "Iteration " << iteration << ": Scanning " << graph.node_count() << " nodes and "
            << graph.edge_count() << " edges for sites..." << endl;
    }
    
    for (const Snarl* top_level_site : site_manager.top_level_snarls()) {
        list<const Snarl*> queue {top_level_site};
        
        while (queue.size()) {
            const Snarl* site = queue.front();
            queue.pop_front();
            
            if (site_manager.is_leaf(site)) {
                // It's a leaf. Filter it out if it is trivial
                
                if (site->type() == ULTRABUBBLE) {
                    auto contents = site_manager.shallow_contents(site, graph, false);
                    if (contents.first.empty()) {
                        // Nothing but the boundary nodes in this snarl
                        continue;
                    }
                }
            
                // Not trivial. Keep it.
                leaves.push_back(site);
            }
            else {
                for (const Snarl* child_site : site_manager.children_of(site)) {
                    queue.push_back(child_site);
                }
            }
        }
    }
    
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
    
    // We can't use the SnarlManager after we modify the graph, so we load the
    // contents of all the leaves we're going to modify first.
    map<const Snarl*, pair<unordered_set<id_t>, unordered_set<edge_t>>> leaf_contents;
    
    // How big is each leaf in bp
    map<const Snarl*, size_t> leaf_sizes;
    
    // We also need to pre-calculate the traversals for the snarls that are the
    // right size, since the traversal finder uses the snarl manager amd might
    // not work if we modify the graph.
    map<const Snarl*, vector<SnarlTraversal>> leaf_traversals;
    
    for (const Snarl* leaf : leaves) {
        // Look at all the leaves
        
        // Get the contents of the bubble, excluding the boundary nodes
        leaf_contents[leaf] = site_manager.deep_contents(leaf, graph, false);
        
        // For each leaf, calculate its total size.
        unordered_set<id_t>& nodes = leaf_contents[leaf].first;
        size_t& total_size = leaf_sizes[leaf];
        for (id_t node_id : nodes) {
            // For each node include it in the size figure
            total_size += graph.get_length(graph.get_handle(node_id));
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
        
        // Identify the replacement traversal for the bubble if it's the right size.
        // We can't necessarily do this after we've modified the graph.
        vector<SnarlTraversal>& traversals = leaf_traversals[leaf];
        traversals = traversal_finder.find_traversals(*leaf);
    }
    
    for (const Snarl* leaf : leaves) {
        // Look at all the leaves
        
        // Get the contents of the bubble, excluding the boundary nodes
        unordered_set<id_t>& nodes = leaf_contents[leaf].first;
        unordered_set<edge_t>& edges = leaf_contents[leaf].second;
        
        // For each leaf, grab its total size.
        size_t& total_size = leaf_sizes[leaf];
        
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
        for (id_t node_id : nodes) {
            cerr << "\t" << node << " = " << node_id << ": " << graph.get_sequence(graph.get_handle(node_id)) << endl;
        }
#endif
        
        // Otherwise we want to simplify this site away
        
        // Grab the replacement traversal for the bubble
        vector<SnarlTraversal>& traversals = leaf_traversals[leaf];
        
        if (traversals.empty()) {
            // We couldn't find any paths through the site.
            continue;
        }
        
        // Get the traversal out of the vector
        SnarlTraversal& traversal = traversals.front();
        
#ifdef debug
        cerr << "Chosen traversal has: " << traversal.visit_size() << " visits" << endl;
#endif

        // Determine the length of the new traversal
        size_t new_site_length = 0;
        for (size_t i = 1; i < traversal.visit_size() - 1; i++) {
            // For every non-anchoring node
            const Visit& visit = traversal.visit(i);

#ifdef debug
            cerr << "Chosen traversal has: " << visit << endl;
#endif

            // Total up the lengths of all the nodes that are newly visited.
            assert(visit.node_id());
            new_site_length += graph.get_node(visit.node_id())->sequence().size();
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
        map<string, set<mapping_t*> > mappings_by_path = graph.paths.get_node_mapping_by_path_name(graph.get_node(leaf->start().node_id()));
        
        // It's possible a path can enter the site through the end node and
        // never hit the start. So we're going to trim those back before we delete nodes and edges.
        map<string, set<mapping_t*> > end_mappings_by_path = graph.paths.get_node_mapping_by_path_name(graph.get_node(leaf->end().node_id()));
        
        if (!drop_hairpin_paths) {
            // We shouldn't drop paths if they hairpin and can't be represented
            // in a simplified bubble. So we instead have to not simplify
            // bubbles that would have that problem.
            bool found_hairpin = false;
            
            for (auto& kv : mappings_by_path) {
                // For each path that hits the start node
                
                if (found_hairpin) {
                    // We only care if there are 1 or more hairpins, not how many
                    break;    
                }
                
                // Unpack the name
                auto& path_name = kv.first;
                
                for (mapping_t* start_mapping : kv.second) {
                    // For each visit to the start node
                
                    if (found_hairpin) {
                        // We only care if there are 1 or more hairpins, not how many
                        break;    
                    }
                
                    // Determine what orientation we're going to scan in
                    bool backward = start_mapping->is_reverse();
                    
                    // Start at the start node
                    mapping_t* here = start_mapping;
                    
                    while (here) {
                        // Until we hit the start/end of the path or the mapping we want
                        if (here->node_id() == leaf->end().node_id() &&
                            here->is_reverse() == (leaf->end().backward() != backward)) {
                            // We made it out.
                            // Stop scanning!
                            break;
                        }
                        
                        if (here->node_id() == leaf->start().node_id() &&
                            here->is_reverse() != (leaf->start().backward() != backward)) {
                            // We have encountered the start node with an incorrect orientation.
                            cerr << "warning:[vg::SmallSnarlSimplifier] Path " << path_name
                                << " doubles back through start of site "
                                << to_node_traversal(leaf->start(), graph) << " - "
                                << to_node_traversal(leaf->end(), graph) << "; skipping site!" << endl;
                                
                            found_hairpin = true;
                            break;
                        }
                        
                        // Scan left along ther path if we found the site start backwards, and right if we found it forwards.
                        here = backward ? graph.paths.traverse_left(here) : graph.paths.traverse_right(here);
                    }
                }
            }
            
            for (auto& kv : end_mappings_by_path) {
                // For each path that hits the end node
                
                if (found_hairpin) {
                    // We only care if there are 1 or more hairpins, not how many
                    break;
                }
                
                // Unpack the name
                auto& path_name = kv.first;
                
                for (mapping_t* end_mapping : kv.second) {
                    
                    if (found_hairpin) {
                        // We only care if there are 1 or more hairpins, not how many
                        break;
                    }
                    
                    // Determine what orientation we're going to scan in
                    bool backward = end_mapping->is_reverse();
                    
                    // Start at the end
                    mapping_t* here = end_mapping;
                    
                    while (here) {
                        
                        if (here->node_id() == leaf->start().node_id() &&
                            here->is_reverse() == (leaf->start().backward() != backward)) {
                            // We made it out.
                            // Stop scanning!
                            break;
                        }
                        
                        if (here->node_id() == leaf->end().node_id() &&
                            here->is_reverse() != (leaf->end().backward() != backward)) {
                            // We have encountered the end node with an incorrect orientation.
                            cerr << "warning:[vg::SmallSnarlSimplifier] Path " << path_name
                                << " doubles back through end of site "
                                << to_node_traversal(leaf->start(), graph) << " - "
                                << to_node_traversal(leaf->end(), graph) << "; dropping site!" << endl;
                            
                            found_hairpin = true;
                            break;
                        }
                        
                        // Scan right along the path if we found the site end backwards, and left if we found it forwards.
                        here = backward ? graph.paths.traverse_right(here) : graph.paths.traverse_left(here);
                        
                    }
                    
                }
                    
            }
            
            if (found_hairpin) {
                // We found a hairpin, so we want to skip the site.
                cerr << "warning:[vg::SmallSnarlSimplifier] Site " << to_node_traversal(leaf->start(), graph) << " - " << to_node_traversal(leaf->end(), graph) << " skipped due to hairpin path." << endl;
                continue;
            }
            
        }
        
        // We'll keep a set of the end mappings we managed to find, starting from the start
        set<mapping_t*> found_end_mappings;
        
        for (auto& kv : mappings_by_path) {
            // For each path that hits the start node
            
            // Unpack the name
            auto& path_name = kv.first;
            
            // If a path can't be represented after a bubble is popped
            // (because the path reversed and came out the same side as it
            // went in), we just clobber the path entirely. TODO: handle
            // out-the-same-side traversals as valid genotypes somehow..
            bool kill_path = false;
            
            for (mapping_t* start_mapping : kv.second) {
                // For each visit to the start node
                
                // Determine what orientation we're going to scan in
                bool backward = start_mapping->is_reverse();
                
                // We're going to fill this list with the mappings we need to
                // remove and replace in this path for this traversal. Initially
                // runs from start of site to end of site, but later gets
                // flipped into path-local orientation.
                list<mapping_t*> existing_mappings;
                
                // Tracing along forward/backward from each as appropriate, see
                // if the end of the site is found in the expected orientation
                // (or if the path ends first).
                bool found_end = false;
                mapping_t* here = start_mapping;
                
                // We want to remember the end mapping when we find it
                mapping_t* end_mapping = nullptr;
                
#ifdef debug
                cerr << "Scanning " << path_name << " from " << *here
                    << " for " << to_node_traversal(leaf->end(), graph) << " orientation " << backward << endl;
#endif
                
                while (here) {
                    // Until we hit the start/end of the path or the mapping we want
                    
#ifdef debug
                    cerr << "\tat " << *here << endl;
#endif
                    
                    if (here->node_id() == leaf->end().node_id() &&
                        here->is_reverse() == (leaf->end().backward() != backward)) {
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
                    
                    if (here->node_id() == leaf->start().node_id() &&
                        here->is_reverse() != (leaf->start().backward() != backward)) {
                        // We have encountered the start node with an incorrect orientation.
                        cerr << "warning:[vg::SmallSnarlSimplifier] Path " << path_name
                            << " doubles back through start of site "
                            << to_node_traversal(leaf->start(), graph) << " - "
                            << to_node_traversal(leaf->end(), graph) << "; dropping!" << endl;
                            
                        assert(drop_hairpin_paths);
                        kill_path = true;
                        break;
                    }
                    
                    if (here->node_id() != leaf->start().node_id() &&
                        here->node_id() != leaf->end().node_id() && 
                        !nodes.count(here->node_id())) {
                        // We aren't the start, the end, or any internal contained node.
                        // That's an error!
                        // We really should stay inside the site!
                        cerr << "error:[vg::SmallSnarlSimplifier] Path " << path_name
                            << " somehow escapes site " << to_node_traversal(leaf->start(), graph)
                            << " - " << to_node_traversal(leaf->end(), graph) << " and reaches non-contained node "
                            << here->node_id() << " at address " << graph.get_node(here->node_id()) << endl;
                            
                        exit(1);
                    }
                    
                    if (here != start_mapping) {
                        // Remember the mappings that aren't to the start or
                        // end of the site, so we can remove them later.
                        existing_mappings.push_back(here);
                    }
                    
                    // Scan left along the path if we found the site start backwards, and right if we found it forwards.
                    mapping_t* next = backward ? graph.paths.traverse_left(here) : graph.paths.traverse_right(here);

#ifdef debug
                    cerr << "Here is " << *here << " at " << here << " and next is " << *next << " at " << next << endl;
#endif
                    
                    if (next == nullptr) {
                        // We hit the end of the path without finding the end of the site.
                        // We've found all the existing mappings, so we can stop.
                        break;
                    }
                    
                    // Make into NodeTraversals
                    NodeTraversal here_traversal(graph.get_node(here->node_id()), here->is_reverse());
                    NodeTraversal next_traversal(graph.get_node(next->node_id()), next->is_reverse());
                    
                    if (backward) {
                        // We're scanning the other way
                        std::swap(here_traversal, next_traversal);
                    }
                    
                    // Make sure we have an edge so we can traverse this node and then the node we're going to.
                    if(graph.get_edge(here_traversal, next_traversal) == nullptr) {
                        cerr << "error:[vg::SmallSnarlSimplifier] No edge " << here_traversal << " to " << next_traversal << endl;
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
                mapping_t* mapping_after_first = existing_mappings.empty() ?
                    (backward ? start_mapping : end_mapping) : existing_mappings.front();
                assert(path_index.mapping_positions.count(mapping_after_first));
                size_t variable_start = path_index.mapping_positions.at(mapping_after_first); 
                
                
                
                // Determine the total length of the old traversal of the site
                size_t old_site_length = 0;
                for (auto* mapping : existing_mappings) {
                    // Add in the lengths of all the mappings that will get
                    // removed.
                    old_site_length += mapping->length;
                }
#ifdef debug
                cerr << "Replacing " << old_site_length << " bp at " << variable_start
                    << " with " << new_site_length << " bp" << endl;
#endif

                // Actually update any BED features
                if (features != nullptr) {
                    features->on_path_edit(path_name, variable_start, old_site_length, new_site_length);
                }
                
                // Where will we insert the new site traversal into the path?
                list<mapping_t>::iterator insert_position;
                
                if (!existing_mappings.empty()) {
                    // If there are existing internal mappings, we'll insert right where they were
                    
                    for (auto* mapping : existing_mappings) {
                        // Remove each mapping from left to right along the
                        // path, saving the position after the mapping we just
                        // removed. At the end we'll have the position of the
                        // mapping to the end of the site.
                        
#ifdef debug
                        cerr << path_name << " forward: Drop mapping " << *mapping << endl;
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
                    assert(insert_position->node_id() == leaf->start().node_id());
                } else {
                    assert(insert_position->node_id() == leaf->end().node_id());
                }
                
#ifdef debug
                cerr << "Chosen traversal has " << traversal.visit_size() << " visits" << endl;
#endif

                // Loop through the internal visits in the canonical
                // traversal backwards along the path we are splicing. If
                // it's a forward path this is just right to left, but if
                // it's a reverse path it has to be left to right.
                for (size_t i = 1; i + 1 < traversal.visit_size(); i++) {
                    // Don't visit the first or last node on the traversal
                    // because they are the snarl start and end which we didn't
                    // remove.

                    // Find the visit we need next, as a function of which
                    // way we need to insert this run of visits. Normally we
                    // go through the visits right to left, but when we have
                    // a backward path we go left to right.
                    const Visit& visit = backward ? traversal.visit(i)
                                                  : traversal.visit(traversal.visit_size() - i - 1);
                    
                    // Make a mapping_t to represent it
                    mapping_t new_mapping;
                    new_mapping.set_node_id(visit.node_id());
                    // We hit this node backward if it's backward along the
                    // traversal, xor if we are traversing the traversal
                    // backward
                    new_mapping.set_is_reverse(visit.backward() != backward);
                    
                    // Add the length
                    size_t node_seq_length = graph.get_node(visit.node_id())->sequence().size();
                    new_mapping.length = node_seq_length;
                    
#ifdef debug
                    cerr << path_name << " backward: Add mapping " << new_mapping << endl;
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
        
        for (auto& kv : end_mappings_by_path) {
            // Now we handle the end mappings not reachable from the start. For each path that touches the end...
        
            // Unpack the name
            auto& path_name = kv.first;
            
            // We might have to kill the path, if it reverses inside a
            // bubble we're popping
            bool kill_path = false;
            
            for (mapping_t* end_mapping : kv.second) {
                if (found_end_mappings.count(end_mapping)) {
                    // Skip the traversals of the site that we handled.
                    continue;
                }
                
                // Now we're left with paths that leave the site but don't
                // enter. We're going to clobber everything before the path
                // leaves the site.
                
                // Determine what orientation we're going to scan in
                bool backward = end_mapping->is_reverse();
                
                // Start at the end
                mapping_t* here = end_mapping;
                
                // Keep a list of mappings we need to remove
                list<mapping_t*> to_remove;
                
                while (here) {
                    
                    if (here->node_id() == leaf->end().node_id() &&
                        here->is_reverse() != (leaf->end().backward() != backward)) {
                        // We have encountered the end node with an incorrect orientation.
                        cerr << "warning:[vg::SmallSnarlSimplifier] Path " << path_name
                            << " doubles back through end of site "
                            << to_node_traversal(leaf->start(), graph) << " - "
                            << to_node_traversal(leaf->end(), graph) << "; dropping!" << endl;
                            
                        assert(drop_hairpin_paths);
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
        for (int i = 0; i < traversal.visit_size() - 1; ++i) {
            // For each node and the next node (which won't be the end)
            
            const Visit visit = traversal.visit(i);
            const Visit next = traversal.visit(i+1);

#ifdef debug
            cerr << "Follow edge from " << visit << " to " << next << endl;
#endif

            // Find the edge between them
            NodeTraversal here(graph.get_node(visit.node_id()), visit.backward());
            NodeTraversal next_traversal(graph.get_node(next.node_id()), next.backward());
            Edge* edge = graph.get_edge(here, next_traversal);
            assert(edge != nullptr);
            
            // Remember we need it
            blessed_edges.insert(edge);
        }

        // The traversal also touches the boundary nodes, so don't do anything special for them.
        
        for (const edge_t& edge_handle : edges) {
            Node* from_node = graph.get_node(graph.get_id(edge_handle.first));
            Node* to_node = graph.get_node(graph.get_id(edge_handle.second));
            Edge* edge = graph.get_edge(NodeTraversal(from_node, graph.get_is_reverse(edge_handle.first)),
                                        NodeTraversal(to_node, graph.get_is_reverse(edge_handle.second)));
            if (!blessed_edges.count(edge)) {
                // Get rid of all the edges not needed for the one true traversal
#ifdef debug
                cerr << to_node_traversal(leaf->start(), graph) << " - "
                     << to_node_traversal(leaf->end(), graph) << ": Delete edge: "
                     << pb2json(*edge) << endl;
#endif
                graph.destroy_edge(edge);
                deleted_edges++;
            }
        }
       
           
        // Now delete all the nodes that aren't on the blessed traversal.
        
        // What nodes are on it?
        set<Node*> blessed_nodes;
        for (int i = 0; i < traversal.visit_size(); i++) {
            const Visit& visit = traversal.visit(i);
            blessed_nodes.insert(graph.get_node(visit.node_id()));
        }
        
        for (id_t node_id : nodes) {
            Node* node = graph.get_node(node_id);
            // For every node in the site
            if (!blessed_nodes.count(node)) {
                // If we don't need it for the chosen path, destroy it
#ifdef debug
                cerr << to_node_traversal(leaf->start(), graph) << " - "
                     << to_node_traversal(leaf->end(), graph) << ": Delete node: "
                     << pb2json(*node) << endl;
#endif
                // There may be paths still touching this node, if they
                // managed to get into the site without touching the start
                // node. We'll delete those paths.
                set<string> paths_to_kill;
                for (auto& kv : graph.paths.get_node_mapping_by_path_name(node)) {
                
                    if (mappings_by_path.count(kv.first)) {
                        // We've already actually updated this path; the
                        // node_mapping data is just out of date.
                        continue;
                    }
                
                    paths_to_kill.insert(kv.first);
                }
                for (auto& path : paths_to_kill) {
                    graph.paths.remove_path(path);
                    cerr << "warning:[vg::SmallSnarlSimplifier] Path " << path << " removed" << endl;
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

void SmallSnarlSimplifier::simplify() {
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
