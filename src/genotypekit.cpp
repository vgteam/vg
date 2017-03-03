#include "genotypekit.hpp"

namespace vg {

using namespace std;

    
CactusUltrabubbleFinder::CactusUltrabubbleFinder(VG& graph,
                                                 const string& hint_path_name,
                                                 bool filter_trivial_bubbles) :
    graph(graph), hint_path_name(hint_path_name), filter_trivial_bubbles(filter_trivial_bubbles) {
    // Make sure the graph is sorted.
    // cactus needs the nodes to be sorted in order to find a source and sink.
    graph.sort();
}

SnarlManager CactusUltrabubbleFinder::find_snarls() {
    
    // Get the bubble tree in Cactus format
    BubbleTree* bubble_tree = ultrabubble_tree(graph);
    
    // Convert to Snarls
    
    vector<Snarl> converted_snarls;
    
    bubble_tree->for_each_preorder([&](BubbleTree::Node* node) {
        
        Bubble& bubble = node->v;
        if (node != bubble_tree->root) {
            // If we aren't the root node of the tree, we need to be a Snarl
            
            if (filter_trivial_bubbles) {
                
                // Check whether the bubble consists of a single edge
                
                set<NodeSide> start_connections = graph.sides_of(bubble.start);
                set<NodeSide> end_connections = graph.sides_of(bubble.end);
                
                if (start_connections.size() == 1
                    && start_connections.count(bubble.end)
                    && end_connections.size() == 1
                    && end_connections.count(bubble.start)) {
                    // This is a single edge bubble, skip it
                    return;
                }
            }
            
            // We're going to fill in this Snarl.
            Snarl snarl;
            
            // Set up the start and end

            // Make sure to preserve original endpoint
            // ordering, because swapping them without flipping their
            // orientation flags will make an inside-out site.
            snarl.mutable_start()->set_node_id(bubble.start.node);
            snarl.mutable_start()->set_backward(!bubble.start.is_end);
            snarl.mutable_end()->set_node_id(bubble.end.node);
            snarl.mutable_end()->set_backward(bubble.end.is_end);
            
            // Mark snarl as an ultrabubble if it's acyclic
            snarl.set_type(bubble.dag ? ULTRABUBBLE : UNCLASSIFIED);
            
            // If not a top level site, add parent info
            if (node->parent != bubble_tree->root) {
                Bubble& bubble_parent = node->parent->v;
                Snarl* snarl_parent = snarl.mutable_parent();
                snarl_parent->mutable_start()->set_node_id(bubble_parent.start.node);
                snarl_parent->mutable_start()->set_backward(!bubble_parent.start.is_end);
                snarl_parent->mutable_end()->set_node_id(bubble_parent.end.node);
                snarl_parent->mutable_end()->set_backward(bubble_parent.end.is_end);
            }
            
            converted_snarls.push_back(snarl);
            
        }
    });
    
    delete bubble_tree;
    
    // Now form the SnarlManager and return
    return SnarlManager(converted_snarls.begin(), converted_snarls.end());
}
    
   
ExhaustiveTraversalFinder::ExhaustiveTraversalFinder(VG& graph, SnarlManager& snarl_manager) :
                                                     graph(graph), snarl_manager(snarl_manager) {
    // nothing more to do
}
    
ExhaustiveTraversalFinder::~ExhaustiveTraversalFinder() {
    // no heap objects
}

void ExhaustiveTraversalFinder::stack_up_valid_walks(NodeTraversal walk_head, vector<NodeTraversal>& stack) {
    
    id_t head_id = walk_head.node->id();
    
    if (walk_head.backward) {
        // we are leaving from the start of the node
        
        // get all edges involving this node so we can filter them down to valid walks
        for (Edge* edge : graph.edges_of(walk_head.node)) {
            if (edge->from() == head_id && edge->from_start()) {
                // the edge is part of a valid walk
                Node* next_node = graph.get_node(edge->to());
                bool next_backward = edge->to_end();
                // add the next traversal in the walk to the stack
                stack.push_back(NodeTraversal(next_node, next_backward));
            }
            else if (edge->to() == head_id && !edge->to_end()) {
                // the edge is part of a valid walk in the opposite orientation
                Node* next_node = graph.get_node(edge->from());
                bool next_backward = !edge->from_start();
                // add the next traversal in the walk to the stack
                stack.push_back(NodeTraversal(next_node, next_backward));
            }
        }
    }
    else {
        // we are leaving from the end of the node
        
        // get all edges involving this node so we can filter them down to valid walks
        for (Edge* edge : graph.edges_of(walk_head.node)) {
            if (edge->from() == head_id && !edge->from_start()) {
                // the edge is part of a valid walk
                Node* next_node = graph.get_node(edge->to());
                bool next_backward = edge->to_end();
                // add the next traversal in the walk to the stack
                stack.push_back(NodeTraversal(next_node, next_backward));
            }
            else if (edge->to() == head_id && edge->to_end()) {
                // the edge is part of a valid walk in the opposite orientation
                Node* next_node = graph.get_node(edge->from());
                bool next_backward = !edge->from_start();
                // add the next traversal in the walk to the stack
                stack.push_back(NodeTraversal(next_node, next_backward));
            }
        }
    }
}

    
vector<SnarlTraversal> ExhaustiveTraversalFinder::find_traversals(const Snarl& site) {

    vector<SnarlTraversal> to_return;
    
    // construct maps that lets us "skip over" child sites
    map<NodeTraversal, const Snarl*> child_site_starts;
    map<NodeTraversal, const Snarl*> child_site_ends;
    for (const Snarl* subsite : snarl_manager.children_of(&site)) {
        child_site_starts[to_node_traversal(subsite->start(), graph)] = subsite;
        // reverse the direction of the end because we want to find it when we're entering
        // the site from that direction
        child_site_ends[to_rev_node_traversal(subsite->end(), graph)] = subsite;
    }
    
    // keeps track of the walk of the DFS traversal
    list<Visit> path;
    
    // these mark the start of the edges out of the node that is on the head of the path
    // they can be used to see how many nodes we need to peel off the path when we're
    // backtracking
    NodeTraversal stack_sentinel(nullptr);
    NodeTraversal site_end = to_node_traversal(site.end(), graph);
    
    // initialize stack for DFS traversal of site
    vector<NodeTraversal> stack;
    stack.push_back(to_node_traversal(site.start(), graph));
    
    while (stack.size()) {
        
        NodeTraversal node_traversal = stack.back();
        stack.pop_back();
        
        // we have traversed all of edges out of the head of the path, so we can pop it off
        if (node_traversal == stack_sentinel) {
            path.pop_back();
            continue;
        }
        
        // have we finished a traversal through the site?
        if (node_traversal == site_end) {
            
            // yield path as a snarl traversal
            SnarlTraversal traversal;
            to_return.push_back(traversal);
            
            // increment past the Snarl's start node, which we don't want in the traversal
            auto iter = path.begin();
            iter++;
            // record the traversal in the return value
            for (; iter != path.end(); iter++) {
                *to_return.back().add_visits() = *iter;
            }
            
            // label which snarl this came from
            *to_return.back().mutable_snarl()->mutable_start() = site.start();
            *to_return.back().mutable_snarl()->mutable_end() = site.end();
            
            // don't proceed to add more onto the DFS stack
            continue;
        }
        
        // mark the beginning of this node/site's edges forward in the stack
        stack.push_back(stack_sentinel);
        
        // initialize empty visit for this iteration
        Visit visit;
        
        if (child_site_starts.count(node_traversal)) {
            // make a visit out of the site
            const Snarl* child_site = child_site_starts[node_traversal];
            transfer_boundary_info(*child_site, *visit.mutable_snarl());
            visit.set_backward(false);
            
            // skip the site and add the other side to the stack
            stack.push_back(to_node_traversal(child_site->end(), graph));
        }
        else if (child_site_ends.count(node_traversal)) {
            // make a visit out of the site
            const Snarl* child_site = child_site_ends[node_traversal];
            transfer_boundary_info(*child_site, *visit.mutable_snarl());
            visit.set_backward(true);
            
            // note: we're traveling through the site backwards, so we reverse the
            // traversal on the start end
            
            // skip the site and add the other side to the stack
            stack.push_back(to_rev_node_traversal(child_site->start(), graph));
        }
        else {
            // make a visit out of the node traversal
            visit.set_node_id(node_traversal.node->id());
            visit.set_backward(node_traversal.backward);
            
            // add all of the node traversals we can reach through valid walks to stack
            stack_up_valid_walks(node_traversal, stack);
        }
        
        // add visit to path
        path.push_back(visit);
    }
    
    return to_return;
}
    
ReadRestrictedTraversalFinder::ReadRestrictedTraversalFinder(VG& graph, SnarlManager& snarl_manager,
                                                             const map<string, Alignment*>& reads_by_name,
                                                             int min_recurrence, int max_path_search_steps) :
                                                             graph(graph), reads_by_name(reads_by_name),
                                                             min_recurrence(min_recurrence),
                                                             max_path_search_steps(max_path_search_steps),
                                                             snarl_manager(snarl_manager) {
    // nothing else to do
}

ReadRestrictedTraversalFinder::~ReadRestrictedTraversalFinder() {
    // no heap variables
}
    
// replaces get_paths_through_site from genotyper
vector<SnarlTraversal> ReadRestrictedTraversalFinder::find_traversals(const Snarl& site) {
    // We're going to emit traversals supported by any paths in the graph.
    
    // Put all our subpaths in here to deduplicate them by sequence they spell
    // out. And to count occurrences. Note that the occurrence count will be
    // boosted to min_recurrence if a non-read path in the graph supports a
    // certain traversal string, so we don't end up dropping unsupported ref
    // alleles.
    map<string, pair<list<Visit>, int>> results;
    
    // construct maps that lets us "skip over" child sites
    map<NodeTraversal, const Snarl*> child_site_starts;
    map<NodeTraversal, const Snarl*> child_site_ends;
    for (const Snarl* subsite : snarl_manager.children_of(&site)) {
        child_site_starts[to_node_traversal(subsite->start(), graph)] = subsite;
        // reverse the direction of the end because we want to find it when we're entering
        // the site from that direction
        child_site_ends[to_rev_node_traversal(subsite->end(), graph)] = subsite;
    }
    
#ifdef debug
#pragma omp critical (cerr)
    cerr << "Looking for paths between " << site.start << " and " << site.end << endl;
#endif
    
    Node* site_start_node = graph.get_node(site.start().node_id());
    Node* site_end_node = graph.get_node(site.end().node_id());
    
    if(graph.paths.has_node_mapping(site_start_node) && graph.paths.has_node_mapping(site_end_node)) {
        // If we have some paths that visit both ends (in some orientation)
        
        // Get all the mappings to the end node, by path name
        auto& endmappings_by_name = graph.paths.get_node_mapping(site_end_node);
        
        for(auto& name_and_mappings : graph.paths.get_node_mapping(site_start_node)) {
            // Go through the paths that visit the start node
            
            // Grab their names
            auto& name = name_and_mappings.first;
            
            if(!endmappings_by_name.count(name_and_mappings.first)) {
                // No path by this name has any mappings to the end node. Skip
                // it early.
                continue;
            }
            
            for(auto* mapping : name_and_mappings.second) {
                // Start at each mapping in the appropriate orientation
                
#ifdef debug
#pragma omp critical (cerr)
                cerr << "Trying mapping of read/path " << name_and_mappings.first << endl;
#endif
                
                // How many times have we gone to the next mapping looking for a
                // mapping to the end node in the right orientation?
                size_t traversal_count = 0;
                
                // Do we want to go left (true) or right (false) from this
                // mapping? If start is a forward traversal and we found a
                // forward mapping, we go right. If either is backward we go
                // left, and if both are backward we go right again.
                bool traversal_direction = mapping->position().is_reverse() != site.start().backward();
                
                // What orientation would we want to find the end node in? If
                // we're traveling backward, we expect to find it in the
                // opposite direction to the one we were given.
                bool expected_end_orientation = site.end().backward() != traversal_direction;
                
                // We're going to fill in this list with traversals.
                list<Visit> path_traversed;
                
                // And we're going to fill this with the sequence
                stringstream allele_stream;
                
                while(mapping != nullptr && traversal_count < max_path_search_steps) {
                    // Traverse along until we hit the end traversal or take too
                    // many steps
                    
#ifdef debug
#pragma omp critical (cerr)
                    cerr << "\tTraversing " << pb2json(*mapping) << endl;
#endif
                    
                    // Say we visit this node along the path, in this orientation
                    NodeTraversal node_traversal = NodeTraversal(graph.get_node(mapping->position().node_id()),
                                                                 mapping->position().is_reverse() != traversal_direction);
                    
                    // Stick the sequence of the node (appropriately oriented) in the stream for the allele sequence
                    string seq = node_traversal.node->sequence();
                    allele_stream << (node_traversal.backward ? reverse_complement(seq) : seq);
                    
                    if(node_traversal.node == site_end_node && node_traversal.backward == expected_end_orientation) {
                        // We have stumbled upon the end node in the orientation we wanted it in.
                        
                        if(results.count(allele_stream.str())) {
                            // It is already there! Increment the observation count.
#ifdef debug
#pragma omp critical (cerr)
                            cerr << "\tFinished; got known sequence " << allele_stream.str() << endl;
#endif
                            
                            if(reads_by_name.count(name)) {
                                // We are a read. Just increment count
                                results[allele_stream.str()].second++;
                            } else {
                                // We are a named path (like "ref")
                                if(results[allele_stream.str()].second < min_recurrence) {
                                    // Ensure that this allele doesn't get
                                    // eliminated, since ref or some other named
                                    // path supports it.
                                    results[allele_stream.str()].second = min_recurrence;
                                } else {
                                    results[allele_stream.str()].second++;
                                }
                            }
                        } else {
                            // Add it in. Give it a count of 1 if we are a read,
                            // and a count of min_recurrence (so it doesn't get
                            // filtered later) if we are a named non-read path
                            // (like "ref").
                            results[allele_stream.str()] = make_pair(path_traversed,
                                                                     reads_by_name.count(name) ? 1 : min_recurrence);
#ifdef debug
#pragma omp critical (cerr)
                            cerr << "\tFinished; got novel sequence " << allele_stream.str() << endl;
#endif
                        }
                        
//                        if(reads_by_name.count(name)) {
//                            // We want to log stats on reads that read all the
//                            // way through sites. But since we may be called
//                            // multiple times we need to send the unique read
//                            // name too.
//                            report_site_traversal(site, name);
//                        }
                        
                        // Then try the next embedded path
                        break;
                    }
                    
                    // We are not yet at the end of the of the site on this path
                    
                    // initialize visit
                    Visit visit;
                    
                    // is this traversal at the start of a nested subsite?
                    Node* site_opposite_side = nullptr;
                    if (child_site_starts.count(node_traversal)) {
                        const Snarl* child_site = child_site_starts[node_traversal];
                        site_opposite_side = graph.get_node(child_site->end().node_id());
                        
                        transfer_boundary_info(*child_site, *visit.mutable_snarl());
                        
                        // add the site into the sequence since we are going to skip it
                        allele_stream << "(" << child_site->start().node_id() << ":" << child_site->end().node_id() << ")";
                        
                    }
                    else if (child_site_ends.count(node_traversal)) {
                        const Snarl* child_site = child_site_starts[node_traversal];
                        site_opposite_side = graph.get_node(child_site->start().node_id());
                        
                        transfer_boundary_info(*child_site, *visit.mutable_snarl());
                        visit.set_backward(true);
                        
                        // add the reverse site into the sequence since we are going to skip it
                        allele_stream << "(" << child_site->end().node_id() << ":" << child_site->start().node_id() << ")";
                    }
                    else {
                        visit = to_visit(node_traversal);
                        allele_stream << node_traversal.node->sequence();
                    }
                    
                    path_traversed.push_back(visit);
                    
                    // Was this node traversal entering a nested subsite?
                    if (site_opposite_side) {
                        // Skip the site
                        if (traversal_direction) {
                            // Go backwards until you hit the other side of the site
                            while (mapping->position().node_id() != site_opposite_side->id()) {
                                mapping = graph.paths.traverse_left(mapping);
                                // Break out of the loop if the path ends before crossing child site
                                if (mapping == nullptr) {
                                    break;
                                }
                                // Tick the counter so we don't go really far on long paths.
                                traversal_count++;
                            }
                        }
                        else {
                            // Go forwards until you hit the other side of the site
                            while (mapping->position().node_id() != site_opposite_side->id()) {
                                mapping = graph.paths.traverse_right(mapping);
                                // Break out of the loop if the path ends before crossing child site
                                if (mapping == nullptr) {
                                    break;
                                }
                                // Tick the counter so we don't go really far on long paths.
                                traversal_count++;
                            }
                        }
                    }
                    else {
                        
                        // Otherwise just move to the right (or left) one position
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
    }
    
    // Now collect the unique results
    vector<SnarlTraversal> to_return;
    
    for(auto& result : results) {
        // Break out each result
        const string& seq = result.first;
        auto& visits = result.second.first;
        auto& count = result.second.second;
        
        if(count < min_recurrence) {
            // We don't have enough initial hits for this sequence to justify
            // trying to re-align the rest of the reads. Skip it. Note that the
            // reference path (and other named paths) will stuff in at least
            // min_recurrence to make sure we don't throw out their alleles.
            continue;
        }
        
        // Send out each list of visits
        to_return.emplace_back();
        for (Visit& visit : visits) {
            *to_return.back().add_visits() = visit;
        }
        
        // label which snarl this came from
        *to_return.back().mutable_snarl()->mutable_start() = site.start();
        *to_return.back().mutable_snarl()->mutable_end() = site.end();
        
    }
    
    return to_return;
}
    
double FixedGenotypePriorCalculator::calculate_log_prior(const Genotype& genotype) {
    // Are all the alleles the same?
    bool all_same = true;
    // What is the common allele number (-1 for unset)
    int allele_value = -1;
    for(size_t i = 0; i < genotype.allele_size(); i++) {
        // For each allele in the genotype
        if(allele_value == -1) {
            // On the first one, grab it. Everyone else will have to match.
            allele_value = genotype.allele(i);
        }
        
        if(allele_value != genotype.allele(i)) {
            // There are two distinct allele values in here
            all_same = false;
            break;
        }
    }
    
    // Return the appropriate prior depending on whether the alleles are all the
    // same (homozygous) or not (heterozygous).
    return all_same ? homozygous_prior_ln : heterozygous_prior_ln;
}

TrivialTraversalFinder::TrivialTraversalFinder(VG& graph) : graph(graph) {
    // Nothing to do!
}

vector<SnarlTraversal> TrivialTraversalFinder::find_traversals(const Snarl& site) {
    assert(site.type() == ULTRABUBBLE);
    
    // We'll fill this in and send it back
    vector<SnarlTraversal> to_return;
    
    // We don't want to be duplicating partial paths, so we store for each
    // NodeTraversal we can reach the previous NodeTraversal we can reach it
    // from.
    map<NodeTraversal, NodeTraversal> previous;
    
    list<NodeTraversal> stack{to_node_traversal(site.start(), graph)};
    
    while (!stack.empty()) { 
        // While there's still stuff on the stack
        
        // Grab the first thing
        NodeTraversal here = stack.front();
        stack.pop_front();
        
        if (here.node->id() == site.end().node_id()) {
            // Trace back a path
            list<NodeTraversal> path;
            
            // Move back one node from the end so it isn't included
            here = previous.at(here);
            
            while (true) {
                // Until we get to the start of the site
                
                if (here.node->id() == site.start().node_id()) {
                    // Stop when we've reached the start of the site
                    break;
                }
                
                // Put this traversal on the front of the path
                path.push_front(here);
                
                // Trace back
                here = previous.at(here);
            }
            
            // Initialize a SnarlTraversal in the return value
            to_return.emplace_back();
            
            // Translate the path into the traversal
            for (NodeTraversal node_traversal : path) {
                *(to_return.back().add_visits()) = to_visit(node_traversal);
            }
            
            // label which snarl this came from
            *to_return.back().mutable_snarl()->mutable_start() = site.start();
            *to_return.back().mutable_snarl()->mutable_end() = site.end();
            
            // Stop early after having found one path
            break;
        } else {
            // We haven't reached the end of the site
            
            for (NodeTraversal next : graph.nodes_next(here)) {
                // Look at all the places we can go from this node
                if (previous.count(next)) {
                    // We already know how to get there.
                    continue;
                }
                
                // Remember how we got there
                previous[next] = here;
                // Explore it, depth first
                stack.push_front(next);
            }
        }
    }
    
    // When we get here, either we found a path, or there isn't one.
    return to_return;
}

RepresentativeTraversalFinder::RepresentativeTraversalFinder(AugmentedGraph& augmented,
    SnarlManager& snarl_manager, PathIndex& index, size_t max_depth,
    size_t max_bubble_paths) : augmented(augmented), snarl_manager(snarl_manager),
    index(index), max_depth(max_depth), max_bubble_paths(max_bubble_paths) {
    
    // Nothing to do!

}

vector<SnarlTraversal> RepresentativeTraversalFinder::find_traversals(const Snarl& site) {
    
    // TODO: we can only do ultrabubbles right now. Other snarls may not have
    // traversals through from end to end.
    assert(site.type() == ULTRABUBBLE);
    
    // Get its nodes and edges (including all child sites)
    pair<unordered_set<Node*>, unordered_set<Edge*>> contents = snarl_manager.deep_contents(&site, augmented.graph, true);
    
    // Copy its node set
    unordered_set<Node*> nodes_left(contents.first);

    // Trace the ref path through the site
    vector<NodeTraversal> ref_path_for_site;
    
    // First figure where the site starts and ends in the primary path
    // TODO: support off-path sites
    size_t site_start = index.by_id.at(site.start().node_id()).first;
    size_t site_end = index.by_id.at(site.end().node_id()).first;
    
    
#ifdef debug
    cerr << "Site starts with " << to_node_traversal(site.start(), augmented.graph)
        << " at " << site_start
        << " and ends with " << to_node_traversal(site.end(), augmented.graph)
        << " at " << site_end << endl;
#endif

    // The primary path may go through the site backward. So get the primary min and max coords
    size_t primary_min = min(site_start, site_end);
    size_t primary_max = max(site_start, site_end);

    // Then walk nodes from min coordinate to max coordinate. This holds the
    // start coordinate of the current node.
    int64_t ref_node_start = primary_min;
    while(ref_node_start <= primary_max) {
    
        // Find the reference node starting here or later.
        auto found = index.by_start.lower_bound(ref_node_start);
        if(found == index.by_start.end()) {
            throw runtime_error("No ref node found when tracing through site!");
        }
        if((*found).first > index.by_id.at(site.end().node_id()).first) {
            // The next reference node we can find is out of the space
            // being replaced. We're done.
            if (verbose) {
                cerr << "Stopping for out-of-bounds node" << endl;
            }
            break;
        }
        
        // Get the corresponding NodeTraversal (by converting through Visit)
        NodeTraversal found_traversal = to_node_traversal(found->second.to_visit(), augmented.graph);
        
        // Add the traversal to the ref path through the site
        ref_path_for_site.push_back(found_traversal);
        
        // Make sure this node is actually in the site
        assert(contents.first.count(found_traversal.node));
        
        // Drop it from the set of nodes in the site we haven't visited.
        nodes_left.erase(found_traversal.node);
        
        // Next iteration look where this node ends.
        ref_node_start = found->first + found_traversal.node->sequence().size();
    }
    
    for(auto node : nodes_left) {
        // Make sure none of the nodes in the site that we didn't visit
        // while tracing along the ref path are on the ref path.
        if(index.by_id.count(node->id())) {
            cerr << "Node " << node->id() << " is on ref path at "
                << index.by_id.at(node->id()).first << " but not traced in site "
                << to_node_traversal(site.start(), augmented.graph) << " to " 
                << to_node_traversal(site.end(), augmented.graph) << endl;
            throw runtime_error("Extra ref node found");
        }
    }
    
    // We need to know all the full-length traversals we're going to consider.
    // XREF states will have to be calculated later, over the whole traversal.
    set<vector<NodeTraversal>> site_traversal_set;
    
    // We have this function to extend a partial traversal into a full
    // traversal and add it as a path. The path must already be rooted on
    // the reference in the correct order and orientation.
    auto extend_into_allele = [&](vector<NodeTraversal> path) {
        // Splice the ref path through the site and the bubble's path
        // through the site together.
        vector<NodeTraversal> extended_path;

        for(auto& traversal : path) {
            // Make sure the site actually has the nodes we're visiting.
            assert(contents.first.count(traversal.node));
#ifdef debug
            if(index.by_id.count(traversal.node->id())) {
                cerr << "Path member " << traversal << " lives on ref at "
                << index.by_id.at(traversal.node->id()).first << endl;
            } else {
                cerr << "Path member " << traversal << " does not live on ref" << endl;
            }
#endif
        }
        
        size_t ref_path_index = 0;
        size_t bubble_path_index = 0;
        
        while(ref_path_for_site.at(ref_path_index) != path.at(bubble_path_index)) {
            // Collect NodeTraversals from the ref path until we hit the one
            // at which the bubble path starts.
            extended_path.push_back(ref_path_for_site[ref_path_index++]);
        }
        while(bubble_path_index < path.size()) {
            // Then take the whole bubble path
            extended_path.push_back(path[bubble_path_index++]);
        }
        while(ref_path_index < ref_path_for_site.size() && ref_path_for_site.at(ref_path_index) != path.back()) {
            // Then skip ahead to the matching point in the ref path
            ref_path_index++;
        }
        if(ref_path_index == ref_path_for_site.size()) {
            // We ran out of ref path before finding the place to leave the alt.
            // This must be a backtracking loop or something; start over from the beginning.
            ref_path_index = 0;
            
            while(ref_path_index < ref_path_for_site.size() && ref_path_for_site.at(ref_path_index) != path.back()) {
                // Then skip ahead to the matching point in the ref path
                ref_path_index++;
            }
            
            if(ref_path_index == ref_path_for_site.size()) {
                // Still couldn't find it!
                stringstream err;
                err << "Couldn't find " << path.back() << " in ref path of site "
                    << to_node_traversal(site.start(), augmented.graph)
                    << " to " << to_node_traversal(site.end(), augmented.graph) << endl;
                throw runtime_error(err.str());
            }
        }
        // Skip the matching NodeTraversal
        ref_path_index++;
        while(ref_path_index < ref_path_for_site.size()) {
            // Then take the entier rest of the ref path
            extended_path.push_back(ref_path_for_site[ref_path_index++]);
        }
        
        // Now add it to the set
        site_traversal_set.insert(extended_path);
    
    };

    for(Node* node : contents.first) {
        // Find the bubble for each node
        
        if(total(augmented.node_supports.at(node)) == 0) {
            // Don't bother with unsupported nodes
            continue;
        }
        
        if(index.by_id.count(node->id())) {
            // Don't try to pathfind to the reference for reference nodes.
            continue;
        }
        
        // TODO: allow bubbles that don't backend into the primary path
        pair<Support, vector<NodeTraversal>> sup_path = find_bubble(node, nullptr);

        vector<NodeTraversal>& path = sup_path.second;
        
        if(path.empty()) {
            // We couldn't find a path back to the primary path. Discard
            // this material.
            if (verbose) {
                cerr << "Warning: No path found for node " << node->id() << endl;
            }
            // TODO: record the node's bases as lost.
            
            // TODO: what if it's already in another bubble/the node is deleted?
            continue;
        }
        
        // Extend it out into an allele
        extend_into_allele(path);
        
    }
    
    for(Edge* edge : contents.second) {
        // Go through all the edges
        
        if(!index.by_id.count(edge->from()) || !index.by_id.count(edge->to())) {
            // Edge doesn't touch reference at both ends. Don't use it
            // because for some reason it makes performance worse
            // overall.
            continue;
        }
        
        // Find a path based around this edge
        pair<Support, vector<NodeTraversal>> sup_path = find_bubble(nullptr, edge);
        vector<NodeTraversal>& path = sup_path.second;
        
#ifdef debug
        cerr << "Edge " << edge->from() << " to " << edge->to() << " yields:" << endl;
        for(auto& traversal : path) {
            cerr << "\t" << traversal << endl;
        }
#endif
        
        if(path.empty()) {
            // We couldn't find a path back to the primary path. Discard
            // this material.
            if (verbose) {
                cerr << "Warning: No path found for edge " << edge->from() << "," << edge->to() << endl;
            }
            // TODO: bases lost
            // TODO: what if it's already in another bubble/the node is deleted?
            continue;
        }
        
        // Extend it out into an allele
        extend_into_allele(path);
    }
    
    
    // Now convert to SnarlTraversals
    vector<SnarlTraversal> unique_traversals;
    
    // Have a function to convert a vector of NodeTraversals including the snarl
    // ends into a SnarlTraversal
    auto emit_traversal = [&](vector<NodeTraversal> node_traversals) {
        // Make it into this SnarlTraversal
        SnarlTraversal trav;
        
        // Label traversal with the snarl
        *trav.mutable_snarl()->mutable_start() = site.start();
        *trav.mutable_snarl()->mutable_end() = site.end();
        
        // Add everything but the first and last nodes as Visits.
        // TODO: think about nested sites?
        
        if (site_start > site_end) {
            // The primary path runs backward, so we need to emit all our lists
            // of NodeTraversals backward so they come out in the snarl's
            // orientation and not the primary path's.
            
            for(size_t i = 1; i + 1 < node_traversals.size(); i++) {
                // Make a Visit for each NodeTraversal but the first and last,
                // but backward and in reverse order.
                *trav.add_visits() = to_visit(node_traversals[node_traversals.size() - 1 - i].reverse());
            }
            
        } else {
            // The primary path and the snarl use the same orientation
            
            for(size_t i = 1; i + 1 < node_traversals.size(); i++) {
                // Make a Visit for each NodeTraversal but the first and last
                *trav.add_visits() = to_visit(node_traversals[i]);
            }
            
        }
        
        // Save the SnarlTraversal
        unique_traversals.push_back(trav);
    };
    
    
    // Do the ref path first
    emit_traversal(ref_path_for_site);
    for(auto& node_traversals : site_traversal_set) {
        // Look at each vector of NodeTraversals
        if (node_traversals != ref_path_for_site) {
            // And do everything other than the ref path
            emit_traversal(node_traversals);            
        }
        
    }
    
    return unique_traversals;
}

pair<Support, vector<NodeTraversal>> RepresentativeTraversalFinder::find_bubble(Node* node, Edge* edge) {

    // What are we going to find our left and right path halves based on?
    NodeTraversal left_traversal;
    NodeTraversal right_traversal;

    if(edge != nullptr) {
        // Be edge-based
        
        // Find the nodes at the ends of the edges. Look at them traversed in the
        // edge's local orientation.
        left_traversal = NodeTraversal(augmented.graph.get_node(edge->from()), edge->from_start());
        right_traversal = NodeTraversal(augmented.graph.get_node(edge->to()), edge->to_end());
        
    } else {
        // Be node-based
        assert(node != nullptr);
        left_traversal = right_traversal = NodeTraversal(node);
    }
    
#ifdef debug
    cerr << "Starting from: " << left_traversal << ", " << right_traversal << endl;
#endif

    // Find paths on both sides, with nodes on the primary path at the outsides
    // and this edge in the middle. Returns path lengths and paths in pairs in a
    // set.
    auto leftPaths = bfs_left(left_traversal);
    auto rightPaths = bfs_right(right_traversal);
    
    // Find a combination of two paths which gets us to the reference in a
    // consistent orientation (meaning that when you look at the ending nodes'
    // Mappings in the reference path, the ones with minimal ranks have the same
    // orientations) and which doesn't use the same nodes on both sides.
    // Track support of up to max_bubble_paths combinations, and return the
    // highest
    pair<Support, vector<NodeTraversal> > bestBubblePath;
    int bubbleCount = 0;
    
    // We need to look in different combinations of lists.
    auto testCombinations = [&](const list<list<NodeTraversal>>& leftList,
        const list<list<NodeTraversal>>& rightList) {

        for(auto leftPath : leftList) {
            // Figure out the relative orientation for the leftmost node.
#ifdef debug        
            cerr << "Left path: " << endl;
            for(auto traversal : leftPath ) {
                cerr << "\t" << traversal << endl;
            }
#endif    
            // Split out its node pointer and orientation
            auto leftNode = leftPath.front().node;
            auto leftOrientation = leftPath.front().backward;
            
            // Get where it falls in the reference as a position, orientation pair.
            auto leftRefPos = index.by_id.at(leftNode->id());
            
            // We have a backward orientation relative to the reference path if we
            // were traversing the anchoring node backwards, xor if it is backwards
            // in the reference path.
            bool leftRelativeOrientation = leftOrientation != leftRefPos.second;
            
            // Make a set of all the nodes in the left path
            set<int64_t> leftPathNodes;
            for(auto visit : leftPath) {
                leftPathNodes.insert(visit.node->id());
            }

            // Get the minimum support in the left path
            Support minLeftSupport = min_support_in_path(leftPath);
            
            for(auto rightPath : rightList) {
                // Figure out the relative orientation for the rightmost node.
#ifdef debug            
                cerr << "Right path: " << endl;
                for(auto traversal : rightPath ) {
                    cerr << "\t" << traversal << endl;
                }
#endif            
                // Split out its node pointer and orientation
                // Remember it's at the end of this path.
                auto rightNode = rightPath.back().node;
                auto rightOrientation = rightPath.back().backward;
                
                // Get where it falls in the reference as a position, orientation pair.
                auto rightRefPos = index.by_id.at(rightNode->id());
                
                // We have a backward orientation relative to the reference path if we
                // were traversing the anchoring node backwards, xor if it is backwards
                // in the reference path.
                bool rightRelativeOrientation = rightOrientation != rightRefPos.second;

                // Get the minimum support in the right path
                Support minRightSupport = min_support_in_path(rightPath);
                
                if(leftRelativeOrientation == rightRelativeOrientation &&
                    ((!leftRelativeOrientation && leftRefPos.first < rightRefPos.first) ||
                    (leftRelativeOrientation && leftRefPos.first > rightRefPos.first))) {
                    // We found a pair of paths that get us to and from the
                    // reference without turning around, and that don't go back to
                    // the reference before they leave.

                    // Get the minimum support of combined left and right paths
                    Support minFullSupport = support_min(minLeftSupport, minRightSupport);
                    
                    // Start with the left path
                    vector<NodeTraversal> fullPath{leftPath.begin(), leftPath.end()};
                    
                    // We need to detect overlap with the left path
                    bool overlap = false;
                    
                    // If we're starting from an edge, we keep the first node on
                    // the right path. If we're starting from a node, we need to
                    // discard it because it's just another copy of our node
                    // we're starting with.
                    for(auto it = (edge != nullptr ? rightPath.begin() : ++(rightPath.begin())); it != rightPath.end(); ++it) {
                        // For all but the first node on the right path, add that in
                        fullPath.push_back(*it);
                        
                        if(leftPathNodes.count((*it).node->id())) {
                            // We already visited this node on the left side. Try
                            // the next right path instead.
                            overlap = true;
                        }
                    }
                    
                    if(overlap) {
                        // Can't combine this right with this left, as they share
                        // nodes and we can't handle the copy number implications.
                        // Try the next right.
                        // TODO: handle the copy number implications.
                        continue;
                    }
                    
                    if(leftRelativeOrientation) {
                        // Turns out our anchored path is backwards.
                        
                        // Reorder everything the other way
                        reverse(fullPath.begin(), fullPath.end());
                        
                        for(auto& traversal : fullPath) {
                            // Flip each traversal
                            traversal = traversal.reverse();
                        }
                    }

                    
#ifdef debug        
                    cerr << "Merged path:" << endl;
                    for(auto traversal : fullPath) {
                        cerr << "\t" << traversal << endl;
                    }                    
#endif
                    // update our best path by seeing if we've found one with higher min support
                    if (total(minFullSupport) > total(bestBubblePath.first)) {
                        bestBubblePath.first = minFullSupport;
                        bestBubblePath.second = fullPath;
                    }

                    // keep things from getting out of hand
                    if (++bubbleCount >= max_bubble_paths) {
                        return bestBubblePath;
                    }
                }
            }
        }
        
        // Return the best path along with its min support
        // (could be empty)
        return bestBubblePath;
        
    };
    
    // Convert sets to lists, which requires a copy again...
    list<list<NodeTraversal>> leftConverted;
    for(auto lengthAndPath : leftPaths) {
        leftConverted.emplace_back(move(lengthAndPath.second));
    }
    list<list<NodeTraversal>> rightConverted;
    for(auto lengthAndPath : rightPaths) {
        rightConverted.emplace_back(move(lengthAndPath.second));
    }
    
    // Look for a valid combination, or return an empty path if one iesn't
    // found.
    return testCombinations(leftConverted, rightConverted);
    
}

Support RepresentativeTraversalFinder::min_support_in_path(const list<NodeTraversal>& path) {
    
    if (path.empty()) {
        return Support();
    }
    auto cur = path.begin();
    auto next = path.begin();
    ++next;
    Support minSupport = augmented.node_supports.count(cur->node) ? augmented.node_supports.at(cur->node) : Support();
    for (; next != path.end(); ++cur, ++next) {
        // check the node support
        Support support = augmented.node_supports.count(next->node) ? augmented.node_supports.at(next->node) : Support();
        minSupport = support_min(minSupport, support);
        
        // check the edge support
        Edge* edge = augmented.graph.get_edge(*cur, *next);
        assert(edge != NULL);
        Support edgeSupport = augmented.edge_supports.count(edge) ? augmented.edge_supports.at(edge) : Support();
        minSupport = support_min(minSupport, edgeSupport);
    }

    return minSupport;
}

set<pair<size_t, list<NodeTraversal>>> RepresentativeTraversalFinder::bfs_left(NodeTraversal node, bool stopIfVisited) {

    // Holds partial paths we want to return, with their lengths in bp.
    set<pair<size_t, list<NodeTraversal>>> toReturn;
    
    // Do a BFS
    
    // This holds the paths to get to NodeTraversals to visit (all of which will
    // end with the node we're starting with).
    list<list<NodeTraversal>> toExtend;
    
    // This keeps a set of all the oriented nodes we already got to and don't
    // need to queue again.
    set<NodeTraversal> alreadyQueued;
    
    // Start at this node at depth 0
    toExtend.emplace_back(list<NodeTraversal> {node});
    // Mark this traversal as already queued
    alreadyQueued.insert(node);
    
    // How many ticks have we spent searching?
    size_t searchTicks = 0;
    // Track how many options we have because size may be O(n).
    size_t stillToExtend = toExtend.size();
    
    while(!toExtend.empty()) {
        // Keep going until we've visited every node up to our max search depth.
        
#ifdef debug
        searchTicks++;
        if(searchTicks % 100 == 0) {
            // Report on how much searching we are doing.
            cerr << "Search tick " << searchTicks << ", " << stillToExtend << " options." << endl;
        }
#endif
        
        // Dequeue a path to extend.
        // Make sure to move out of the list to avoid a useless copy.
        list<NodeTraversal> path(move(toExtend.front()));
        toExtend.pop_front();
        stillToExtend--;
        
        // We can't just throw out longer paths, because shorter paths may need
        // to visit a node twice (in opposite orientations) and thus might get
        // rejected later. Or they might overlap with paths on the other side.
        
        // Look up and see if the front node on the path is on our reference
        // path
        if(index.by_id.count(path.front().node->id())) {
            // This node is on the reference path. TODO: we don't care if it
            // lands in a place that is itself deleted.
            
            // Say we got to the right place
            toReturn.emplace(bp_length(path), move(path));
            
            // Don't bother looking for extensions, we already got there.
        } else if(path.size() <= max_depth) {
            // We haven't hit the reference path yet, but we also haven't hit
            // the max depth. Extend with all the possible extensions.
            
            // Look left
            vector<NodeTraversal> prevNodes;
            augmented.graph.nodes_prev(path.front(), prevNodes);
            
            for(auto prevNode : prevNodes) {
                // For each node we can get to
                Edge* edge = augmented.graph.get_edge(prevNode, path.front());
                assert(edge != NULL);
                
                if((!augmented.node_supports.empty() && (!augmented.node_supports.count(prevNode.node) ||
                    total(augmented.node_supports.at(prevNode.node)) == 0)) ||
                   (!augmented.edge_supports.empty() && (!augmented.edge_supports.count(edge) ||
                    total(augmented.edge_supports.at(edge)) == 0))) {
                    
                    // We have no support at all for visiting this node (but we
                    // do have some node read support data)
                    continue;
                }
                
                if(stopIfVisited && alreadyQueued.count(prevNode)) {
                    // We already have a way to get here.
                    continue;
                }
            
                // Make a new path extended left with the node
                list<NodeTraversal> extended(path);
                extended.push_front(prevNode);
                toExtend.emplace_back(move(extended));
                stillToExtend++;
                
                // Remember we found a way to this node, so we don't try and
                // visit it other ways.
                alreadyQueued.insert(prevNode);
            }
        }
        
    }
    
    return toReturn;
}

set<pair<size_t, list<NodeTraversal>>> RepresentativeTraversalFinder::bfs_right(NodeTraversal node, bool stopIfVisited) {

    // Look left from the backward version of the node.
    auto toConvert = bfs_left(node.reverse(), stopIfVisited);
    
    // Since we can't modify set records in place, we need to do a copy
    set<pair<size_t, list<NodeTraversal>>> toReturn;
    
    for(auto lengthAndPath : toConvert) {
        // Flip every path to run the other way
        lengthAndPath.second.reverse();
        for(auto& traversal : lengthAndPath.second) {
            // And invert the orientation of every node in the path in place.
            traversal = traversal.reverse();
        }
        // Stick it in the new set
        toReturn.emplace(move(lengthAndPath));
    }
    
    return toReturn;
}

size_t RepresentativeTraversalFinder::bp_length(const list<NodeTraversal>& path) {
    size_t length = 0;
    for(auto& traversal : path) {
        // Sum up length of each node's sequence
        length += traversal.node->sequence().size();
    }
    return length;
}

int total(const Support& support) {
    return support.forward() + support.reverse();
}

Support support_min(const Support& a, const Support& b) {
    Support to_return;
    to_return.set_forward(min(a.forward(), b.forward()));
    to_return.set_reverse(min(a.reverse(), b.reverse()));
    return to_return;
}


}
