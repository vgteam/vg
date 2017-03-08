/// \file nested_traversal_finder.cpp: Contains implementation for a TraversalFinder that works on non-leaf Snarls.

#include "genotypekit.hpp"

namespace vg {

using namespace std;

/**
 * This TraversalFinder emits at least one traversal representing every node,
 * edge, or child Snarl. Only works on ultrabubbles, and so does not handle
 * cycles.
 */
class NestedTraversalFinder : public TraversalFinder {

protected:
    /// The annotated, augmented graph we're finding traversals in
    AugmentedGraph& augmented;
    /// The SnarlManager managiung the snarls we use
    SnarlManager& snarl_manager;
    
    /**
     * Given an edge or node or child snarl in the augmented graph, look out
     * from the edge or node or child in both directions to find a shortest
     * bubble connecting the start and end of the given site.
     *
     * Exactly one of edge and node and child must be non-null.
     *
     * Return the found traversal as a vector of Visits, including anchoring
     * Visits to the site's start and end nodes. Also return the minimum support
     * found on any edge or node in the bubble that is not contained within a
     * child.
     *
     * If there is no path with any support, returns a zero Support and an empty
     * vector.
     */
    pair<Support, vector<Visit>> find_bubble(Node* node, Edge* edge, Snarl* child,
        const Snarl& site, const map<NodeTraversal, const Snarl*>& child_boundary_index);
        
    /**
     * Get the minimum support of all nodes and edges used in the given path
     * that are not inside child snarls.
     */
    Support min_support_in_path(const vector<Visit>& path);
        
    /**
     * Do a breadth-first search left from the given node traversal, and return
     * lengths (in visits) and paths starting at the given node and ending on
     * the given indexed path. Refuses to visit nodes with no support.
     *
     * Lengths are included so that shorter paths sort first.
     */
    set<pair<size_t, list<Visit>>> search_left(const Visit& root, const Snarl& site,
        const map<NodeTraversal, const Snarl*>& child_boundary_index);
        
    /**
     * Do a breadth-first search right from the given node traversal, and return
     * lengths (in visits) and paths starting at the given node and ending on
     * the given indexed path.
     *
     * Lengths are included so that shorter paths sort first.
     */
    set<pair<size_t, list<Visit>>> search_right(const Visit& root, const Snarl& site,
        const map<NodeTraversal, const Snarl*>& child_boundary_index);
        
    /**
     * Get the length of a path through nodes and child sites, in base pairs.
     * Ignores any bases inside child sites.
     */
    size_t bp_length(const list<Visit>& path);
    
public:

    NestedTraversalFinder(AugmentedGraph& augmented, SnarlManager& snarl_manager);
    
    /// Should we emit verbose debugging info?
    bool verbose = false;
    
    virtual ~NestedTraversalFinder() = default;
    
    /**
     * Find traversals to cover the nodes, edges, and children of the snarl.
     * Always emits the primary path traversal first, if applicable.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    
};


NestedTraversalFinder::NestedTraversalFinder(AugmentedGraph& augmented,
    SnarlManager& snarl_manager) : augmented(augmented), snarl_manager(snarl_manager) {
    
    // Nothing to do!

}

vector<SnarlTraversal> NestedTraversalFinder::find_traversals(const Snarl& site) {
    
    // For each node, edge, and child snarl, try and find a traversal that visits it
}

pair<Support, vector<Visit>> NestedTraversalFinder::find_bubble(Node* node, Edge* edge, Snarl* child, const Snarl& site, const map<NodeTraversal, const Snarl*>& child_boundary_index) {

    // What are we going to find our left and right path halves based on?
    Visit left_visit;
    Visit right_visit;

    if(edge != nullptr) {
        // Be edge-based
        
        // TODO: if the edge has no support, don't go looking for a path
        
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
    auto left_paths = search_left(left_visit, site, child_boundary_index);
    auto right_paths = search_right(right_visit, site, child_boundary_index);
    
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
            Node* node = augmented.graph.get_node(v.node_id());
            
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
        Edge* edge = augmented.graph.get_edge(to_left_side(*cur), to_right_side(*next));
        assert(edge != NULL);
        Support edge_support = augmented.edge_supports.count(edge) ? augmented.edge_supports.at(edge) : Support();
        min_support = support_min(min_support, edge_support);
    }

    return min_support;
}

set<pair<size_t, list<Visit>>> NestedTraversalFinder::search_left(const Visit& root, const Snarl& site,
    const map<NodeTraversal, const Snarl*>& child_boundary_index) {
    
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
            return to_return;
        } else {
            // We haven't reached a boundary, so just keep searching left from
            // this visit.
            vector<Visit> neighbors = visits_left(to_extend_from, augmented.graph, child_boundary_index);
            
            for (auto& extension : neighbors) {
                // For each thing off the left
                
                if (parents.count(extension)) {
                    // This next node is already reachable by a path this short or shorter.
                    continue;
                }
                
                // Check the edge to it to make sure it has coverage
                Edge* edge = augmented.graph.get_edge(to_right_side(extension), to_left_side(to_extend_from));
                
                if (!augmented.edge_supports.count(edge) || total(augmented.edge_supports.at(edge)) == 0) {
                    // This edge is not supported, so don't explore this extension.
                    continue;
                }

                // Look up the node we're entering (either the snarl boundary or
                // just the node we're going to visit), so we can check to make
                // sure it has coverage.
                Node* node = augmented.graph.get_node(to_right_side(extension).node);
                
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
    
    // If we get here, no path with any support was found. Return our still-empty set.
    return to_return;
}

set<pair<size_t, list<Visit>>> NestedTraversalFinder::search_right(const Visit& root, const Snarl& site,
    const map<NodeTraversal, const Snarl*>& child_boundary_index) {

    // Make a backwards version of the root
    Visit root_rev = root;
    root_rev.set_backward(!root_rev.backward());

    // Look left from the backward version of the root
    auto to_convert = search_left(root_rev, site, child_boundary_index);
    
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
