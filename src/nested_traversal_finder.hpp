#ifndef VG_NESTED_TRAVERSAL_FINDER_H
#define VG_NESTED_TRAVERSAL_FINDER_H

/// \file nested_traversal_finder.hpp: Defines a TraversalFinder that produces
/// covering paths, with an awareness of nested Snarls.

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
     * If there is no path with any support, returns a zero Support and a
     * possibly empty Path.
     */
    pair<Support, vector<Visit>> find_bubble(Node* node, Edge* edge, const Snarl* child,
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

}

#endif
