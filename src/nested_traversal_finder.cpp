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
     */
    pair<Support, vector<Visit>> find_bubble(Node* node, Edge* edge, Snarl* child, const Snarl& site);
        
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
    set<pair<size_t, list<Visit>>> search_left(Visit root, const Snarl& site,
        const map<NodeTraversal, const Snarl*>& child_boundary_index);
        
    /**
     * Do a breadth-first search right from the given node traversal, and return
     * lengths (in visits) and paths starting at the given node and ending on
     * the given indexed path.
     *
     * Lengths are included so that shorter paths sort first.
     */
    set<pair<size_t, list<Visit>>> search_right(Visit root, const Snarl& site,
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
    
    // TODO: implement
}

pair<Support, vector<Visit>> NestedTraversalFinder::find_bubble(Node* node, Edge* edge, Snarl* child, const Snarl& site) {

    // TODO: implement
    
}

Support NestedTraversalFinder::min_support_in_path(const vector<Visit>& path) {
    
    // TODO: implement
}

set<pair<size_t, list<Visit>>> NestedTraversalFinder::search_left(Visit root, const Snarl& site,
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
            
            
            // TODO: implement
            // For each thing off the left
            
            // Check the edge to it to make sure it has coverage
            
            // If the next thing is a node, check it to make sure it has coverage
            
            // If the next thing is a snarl, maybe check its attached boundary node to make sure it has coverage?
            
            // If it checks out, queue each neighbor with this visit as its parent
            
        }
        
    }
    
    // We should never get here becuase we should eventually reach the end of
    // the site we're searching.
    throw runtime_error("Got lost in site in search_left!");
}

set<pair<size_t, list<Visit>>> NestedTraversalFinder::search_right(Visit root, const Snarl& site,
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
