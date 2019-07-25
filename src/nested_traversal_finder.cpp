/// \file nested_traversal_finder.cpp: Contains implementation for a
/// TraversalFinder that works on non-leaf Snarls and produces covering paths.

#include "nested_traversal_finder.hpp"

namespace vg {

using namespace std;


NestedTraversalFinder::NestedTraversalFinder(SupportAugmentedGraph& augmented,
    SnarlManager& snarl_manager) : augmented(augmented), snarl_manager(snarl_manager) {
    
    // Nothing to do!

}

vector<SnarlTraversal> NestedTraversalFinder::find_traversals(const Snarl& site) {
    
    // We will populate this
    set<SnarlTraversal> to_return;
    
    // We need a function to check supports and convert things to
    // SnarlTraversals
    auto emit_path = [&](const pair<Support, vector<Visit>>& bubble) {
        if (total(bubble.first) == 0) {
            // There's no support for this
            return;
        }
        
        // This is a legit path to consider
        
        // Make a traversal of the path
        SnarlTraversal trav;
        
        // We need to have the anchoring boundary nodes.
        assert(bubble.second.size() >= 2);
        
        // Do we need to reverse the bubble so it hits the start on the left and the end on the right?
        bool flip_path;
        if (bubble.second.front() == site.start()) {
            // We're already forward
            flip_path = false;
        } else {
            // We must be backward
            flip_path = true;
            // Make sure we actually start with the right node for being backward.
            assert(site.end().node_id() == bubble.second.front().node_id());
        }
        
        for (size_t i = 0; i < bubble.second.size(); i++) {
            // For every visit
            
            // Stick the forward or flipped visit in the traversal
            *trav.add_visit() = flip_path ? reverse(bubble.second[bubble.second.size() - 1 - i]) : bubble.second[i];;
        }
            
        // Now emit the actual traversal
        to_return.insert(trav);
    };
    
    // Get our contained nodes and edges
    unordered_set<id_t> nodes;
    unordered_set<edge_t> edges;
    
    // Grab them, including child boundaries but not our boundaries (which we're
    // guaranteed to visit)
    tie(nodes, edges) = snarl_manager.shallow_contents(&site, augmented.graph, false);
    
    for(auto it = nodes.begin(); it != nodes.end(); ) {
        // For each node
        if (snarl_manager.into_which_snarl(*it, false) ||
            snarl_manager.into_which_snarl(*it, true)) {
        
            // If the node is a child boundary, don't use it. Use visits to the
            // child instead.
            it = nodes.erase(it);
            
        } else {
            // Not a child boundary. Try the next one.
            ++it;
        }
    }  
    
    for (id_t node_id : nodes) {
        // Find bubbles for nodes
        Node* node = augmented.graph.get_node(node_id);
        emit_path(find_bubble(node, nullptr, nullptr, site));
    }
    
    for (const edge_t& edge_handle : edges) {
        Edge* edge = augmented.graph.get_edge(
            NodeTraversal(augmented.graph.get_node(augmented.graph.get_id(edge_handle.first)),
                          augmented.graph.get_is_reverse(edge_handle.first)),
            NodeTraversal(augmented.graph.get_node(augmented.graph.get_id(edge_handle.second)),
                          augmented.graph.get_is_reverse(edge_handle.second)));
        // Find bubbles for edges
        emit_path(find_bubble(nullptr, edge, nullptr, site));
    }
    
    for (const Snarl* child : snarl_manager.children_of(&site)) {
        // Find bubbles for children
        emit_path(find_bubble(nullptr, nullptr, child, site));
    }
    
    // Convert to a vector and return
    return vector<SnarlTraversal> {to_return.begin(), to_return.end()};
    
}

pair<Support, vector<Visit>> NestedTraversalFinder::find_bubble(Node* node, Edge* edge, const Snarl* child, const Snarl& site) {

    // What are we going to find our left and right path halves based on?
    Visit left_visit;
    Visit right_visit;
    
    if(edge != nullptr) {
        // Be edge-based
        
        // Find the nodes at the ends of the edges. Look at them traversed in the
        // edge's local orientation.
        left_visit.set_node_id(edge->from());
        left_visit.set_backward(edge->from_start());
        
        right_visit.set_node_id(edge->to());
        right_visit.set_backward(edge->to_end());
    } else if (node != nullptr) {
        // Be node-based. Both roots are the same visit
        left_visit.set_node_id(node->id());
        right_visit = left_visit;
    } else {
        // We must be child-based
        assert(child != nullptr);
        
        // Both roots are the same child snarl visit
        transfer_boundary_info(*child, *left_visit.mutable_snarl());
        right_visit = left_visit;
    }
    
    // Find paths on both sides, anchored into the outside of our snarl. Returns
    // path lengths (in visits) and paths in pairs in a set.
    auto left_paths = search_left(left_visit, site);
    auto right_paths = search_right(right_visit, site);
    
    // Now splice the paths together to get max support.
    
    // TODO: implement that. For now just get any path with nonzero support.
    
    for (auto& left : left_paths) {
        for (auto& right : right_paths) {
            // Splice together the first path from each set
            
            // Start witht he whole left path
            vector<Visit> spliced_path {left.second.begin(), left.second.end()};
            
            auto it = right.second.begin();
            if (it != right.second.end() && !spliced_path.empty() && *it == spliced_path.back()) {
                // If the right path starts with the same visit the left path ended with, skip it.
                ++it;
            }
            
            // Copy over the rest of the right side path
            copy(it, right.second.end(), back_inserter(spliced_path));
            
            // Return this spliced-together path with its support
            return make_pair(min_support_in_path(spliced_path), spliced_path);
            
        }
    }
    
    // If we get here there's no pair of paths to combine. Return the zero value.
    return pair<Support, vector<Visit>>();
    
}

Support NestedTraversalFinder::min_support_in_path(const vector<Visit>& path) {
    
    if (path.empty()) {
        return Support();
    }
    // We look at the current visit and the next visit so we can think about
    // edges.
    auto cur = path.begin();
    auto next = path.begin();
    ++next;
    
    // We have a function to get the support for a visit
    function<Support(const Visit&)> support_for_visit = [&](const Visit& v) { 
        if (v.node_id()) {
            // This is a node visit
            id_t node = v.node_id();
            
            // Return the support for it, or 0 if it's not in the map.
            return augmented.node_supports.count(node) ? augmented.node_supports.at(node) : Support();
        } else {
            // It's a snarl visit. We assume it goes in one side and out the
            // other.
            
            // We don't inspect the whole snarl, but we know you have to visit
            // the start and end nodes, so we look at them.
            return support_min(support_for_visit(v.snarl().start()), support_for_visit(v.snarl().end()));
        }
    }; 
    
    // Start out with the support for the current visit.
    Support min_support = support_for_visit(*cur);
    for (; next != path.end(); ++cur, ++next) {
        // For each subsequent node next
    
        // check the node support
        min_support = support_min(min_support, support_for_visit(*next));
        
        // check the edge support
        NodeSide from_side = to_right_side(*cur);
        NodeSide to_side = to_left_side(*next);
        edge_t edge = augmented.graph.edge_handle(augmented.graph.get_handle(from_side.node, !from_side.is_end),
                                                  augmented.graph.get_handle(to_side.node, to_side.is_end));

        assert(augmented.graph.has_edge(edge.first, edge.second));
        Support edge_support = augmented.edge_supports.count(edge) ? augmented.edge_supports.at(edge) : Support();
        min_support = support_min(min_support, edge_support);
    }

    return min_support;
}

set<pair<size_t, list<Visit>>> NestedTraversalFinder::search_left(const Visit& root, const Snarl& site) {
    
    // Find paths left from the given visit to the start or end of the given
    // ultrabubble.
    
    // Right now just finds the shortest path in nodes.

    // Holds partial paths we want to return, with their lengths in visits.
    set<pair<size_t, list<Visit>>> to_return;
    
    // Do a BFS
    
    // This holds the parent of each Visit in the queue, back to the initial
    // rooting visit.
    map<Visit, Visit> parents;
    
    // This holds all the visits we need to look at.
    list<Visit> queue {root};
    
    const Snarl* managed_site = snarl_manager.manage(site);
    
    while(!queue.empty()) {
        // Keep going until we've emptied the queue
        
        // Dequeue a visit to continue from.
        Visit to_extend_from = queue.front();
        queue.pop_front();
        
        if (to_extend_from == site.start() || to_extend_from == site.end()) {
            // We're done! Do a trace back.
            
            // Fill in this path
            list<Visit> path;           
            
            // With this cursor
            Visit v = to_extend_from;
            
            while (v != root) {
                // Until we get to the root, put this node on the path and get
                // its parent.
                path.push_back(v);
                v = parents.at(v);
            }
            // We hit the root so add it to the path too.
            path.push_back(v);
            
            // Add the path and its length to the list to return
            to_return.insert(make_pair(path.size(), path));
            
            // And go ahead and return it right now (since it's maximally short)
            break;
        } else {
            // We haven't reached a boundary, so just keep searching left from
            // this visit.
            vector<Visit> neighbors = snarl_manager.visits_left(to_extend_from, augmented.graph, managed_site);
            
            for (auto& extension : neighbors) {
                // For each thing off the left
                
                if (parents.count(extension)) {
                    // This next node is already reachable by a path this short or shorter.
                    continue;
                }
                
                // Check the edge to it to make sure it has coverag
                NodeSide from_side = to_right_side(extension);
                NodeSide to_side = to_left_side(to_extend_from);
                edge_t edge = augmented.graph.edge_handle(augmented.graph.get_handle(from_side.node, !from_side.is_end),
                                                  augmented.graph.get_handle(to_side.node, to_side.is_end));
                
                if (!augmented.edge_supports.count(edge) || total(augmented.edge_supports.at(edge)) == 0) {
                    // This edge is not supported, so don't explore this extension.
                    continue;
                }

                // Look up the node we're entering (either the snarl boundary or
                // just the node we're going to visit), so we can check to make
                // sure it has coverage.
                id_t node = to_right_side(extension).node;
                
                if (!augmented.node_supports.count(node) || total(augmented.node_supports.at(node)) == 0) {
                    // This node is not supported, so don't explore this extension.
                    continue;
                }
                
                // If it checks out, queue the neighbor with this visit as its parent
                parents[extension] = to_extend_from;
                queue.push_back(extension);
            }
            
        }
        
    }
    
    // If we get here, either we broke out of search or no path with any support was found.
    // Return the shortest path or our still-empty set.
    return to_return;
}

set<pair<size_t, list<Visit>>> NestedTraversalFinder::search_right(const Visit& root, const Snarl& site) {

    // Make a backwards version of the root
    Visit root_rev = root;
    root_rev.set_backward(!root_rev.backward());

    // Look left from the backward version of the root
    auto to_convert = search_left(root_rev, site);
    
    // Since we can't modify set records in place, we need to do a copy
    set<pair<size_t, list<Visit>>> to_return;
    
    for(auto length_and_path : to_convert) {
        // Flip every path to run the other way
        length_and_path.second.reverse();
        for(auto& visit : length_and_path.second) {
            // And invert the orientation of every visit in the path in place.
            visit.set_backward(!visit.backward());
        }
        // Stick it in the new set
        to_return.emplace(move(length_and_path));
    }
    
    return to_return;
}

}
