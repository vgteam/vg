#ifndef VG_TRAVERSAL_FINDER_HPP_INCLUDED
#define VG_TRAVERSAL_FINDER_HPP_INCLUDED
// traversal_finder.hpp: TraversalFinder interfaces

// TraversalFinder is an interface for finding paths through Snarls.  The various
// subclasses implement different algorithms for doing this.  A traversal of a snarl
// can be thought of as an allele in a site.

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <regex>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <list>
#include "vg.pb.h"
#include "vg.hpp"
#include "translator.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "types.hpp"
#include "snarls.hpp"
#include "path_index.hpp"
#include "genotypekit.hpp"

namespace vg {

using namespace std;

class AugmentedGraph;
class SupportAugmentedGraph;

/**
 * Represents a strategy for finding traversals of (nested) sites. Polymorphic
 * base class/interface.
 */
class TraversalFinder {
public:
    virtual ~TraversalFinder() = default;
    
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site) = 0;
};

class ExhaustiveTraversalFinder : public TraversalFinder {
    
    VG& graph;
    SnarlManager& snarl_manager;
    bool include_reversing_traversals;
    
public:
    ExhaustiveTraversalFinder(VG& graph, SnarlManager& snarl_manager,
                              bool include_reversing_traversals = false);
    
    virtual ~ExhaustiveTraversalFinder();
    
    /**
     * Exhaustively enumerate all traversals through the site. Only valid for
     * acyclic Snarls.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    
private:
    void stack_up_valid_walks(NodeTraversal walk_head, vector<NodeTraversal>& stack);
    void add_traversals(vector<SnarlTraversal>& traversals, NodeTraversal traversal_start,
                        set<NodeTraversal>& stop_at, set<NodeTraversal>& yield_at);
    
};
    
class ReadRestrictedTraversalFinder : TraversalFinder {
    
    VG& graph;
    SnarlManager& snarl_manager;
    const map<string, Alignment*>& reads_by_name;
    
    // How many times must a path recur before we try aligning to it? Also, how
    // many times must a node in the graph be visited before we use it in indel
    // realignment for nearby indels? Note that the primary path counts as a
    // recurrence. TODO: novel inserts can't recur, and novel deletions can't be
    // filtered in this way.
    int min_recurrence;
    
    // How many nodes max should we walk when checking if a path runs through a superbubble/site
    int max_path_search_steps;
    
public:
    ReadRestrictedTraversalFinder(VG& graph, SnarlManager& snarl_manager, const map<string,
                                  Alignment*>& reads_by_name, int min_recurrence = 2,
                                  int max_path_search_steps = 100);
    
    virtual ~ReadRestrictedTraversalFinder();
    
    /**
     * For the given site, emit all traversals with unique sequences that run from
     * start to end, out of the paths in the graph. Uses the map of reads by
     * name to determine if a path is a read or a real named path. Paths through
     * the site supported only by reads are subject to a min recurrence count,
     * while those supported by actual embedded named paths are not.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    
};

class PathBasedTraversalFinder : public TraversalFinder{
    vg::VG& graph;
    SnarlManager& snarlmanager;
    public:
    PathBasedTraversalFinder(vg::VG& graph, SnarlManager& sm);
    virtual ~PathBasedTraversalFinder() = default;
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);

};

/**
 * This traversal finder finds one or more traversals through leaf sites with no
 * children. It uses a depth-first search. It doesn't work on non-leaf sites,
 * and is not guaranteed to find all traversals. Only works on ultrabubbles.
 */
class TrivialTraversalFinder : public TraversalFinder {

    // Holds the vg graph we are looking for traversals in.
    VG& graph;

public:
    TrivialTraversalFinder(VG& graph);

    virtual ~TrivialTraversalFinder() = default;
    
    /**
     * Find at least one traversal of the site by depth first search, if any
     * exist. Only works on sites with no children.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
};

/**
 * This TraversalFinder is derived from the old vg call code, and emits at least
 * one traversal representing every node, and one traversal representing every
 * edge.
 */
class RepresentativeTraversalFinder : public TraversalFinder {

protected:
    /// The annotated, augmented graph we're finding traversals in
    SupportAugmentedGraph& augmented;
    /// The SnarlManager managiung the snarls we use
    SnarlManager& snarl_manager;
    
    /// We keep around a function that can be used to get an index for the
    /// appropriate path to use to scaffold a given site, or null if no
    /// appropriate index exists.
    function<PathIndex*(const Snarl&)> get_index;
    
    /// What DFS depth should we search to?
    size_t max_depth;
    /// How many DFS searches should we let there be on the stack at a time?
    size_t max_width;
    /// How many search intermediates can we allow?
    size_t max_bubble_paths;
    
    /**
     * Find a Path that runs from the start of the given snarl to the end, which
     * we can use to backend our traversals into when a snarl is off the primary
     * path.
     */
    Path find_backbone(const Snarl& site);
    
    /**
     * Given an edge or node in the augmented graph, look out from the edge or
     * node or snarl in both directions to find a shortest bubble relative to
     * the path, with a consistent orientation. The bubble may not visit the
     * same node twice.
     *
     * Exactly one of edge and node and snarl must be not null.
     *
     * Takes a max depth for the searches producing the paths on each side.
     * 
     * Return the ordered and oriented nodes in the bubble, with the outer nodes
     * being oriented forward along the path for which an index is provided, and
     * with the first node coming before the last node in the reference.  Also
     * return the minimum support found on any edge or node in the bubble
     * (including the reference node endpoints and their edges which aren't
     * stored in the path).
     */
    pair<Support, vector<Visit>> find_bubble(Node* node, Edge* edge, const Snarl* snarl, PathIndex& index,
                                             const Snarl& site);
        
    /**
     * Get the minimum support of all nodes and edges in path
     */
    Support min_support_in_path(const list<Visit>& path);
        
    /**
     * Do a breadth-first search left from the given node traversal, and return
     * lengths and paths starting at the given node and ending on the given
     * indexed path. Refuses to visit nodes with no support, if support data is
     * available in the augmented graph.
     */
    set<pair<size_t, list<Visit>>> bfs_left(Visit visit, PathIndex& index, bool stopIfVisited = false,
                                            const Snarl* in_snarl = nullptr);
        
    /**
     * Do a breadth-first search right from the given node traversal, and return
     * lengths and paths starting at the given node and ending on the given
     * indexed path. Refuses to visit nodes with no support, if support data is
     * available in the augmented graph.
     */
    set<pair<size_t, list<Visit>>> bfs_right(Visit visit, PathIndex& index, bool stopIfVisited = false,
                                             const Snarl* in_snarl = nullptr);
        
    /**
     * Get the length of a path through nodes, in base pairs.
     */
    size_t bp_length(const list<Visit>& path);
    
public:

    /**
     * Make a new RepresentativeTraversalFinder to find traversals. Uses the
     * given augmented graph as the graph with coverage annotation, and reasons
     * about child snarls with the given SnarlManager. Explores up to max_depth
     * in the BFS search when trying to find its way across snarls, and
     * considers up to max_width search states at a time. When combining search
     * results on either side of a graph element to be represented, thinks about
     * max_bubble_paths combinations.
     *
     * Uses the given get_index function to try and find a PathIndex for a
     * reference path traversing a child snarl.
     */
    RepresentativeTraversalFinder(SupportAugmentedGraph& augmented, SnarlManager& snarl_manager,
        size_t max_depth, size_t max_width, size_t max_bubble_paths,
        function<PathIndex*(const Snarl&)> get_index = [](const Snarl& s) { return nullptr; });
    
    /// Should we emit verbose debugging info?
    bool verbose = false;
    
    virtual ~RepresentativeTraversalFinder() = default;
    
    /**
     * Find traversals to cover the nodes and edges of the snarl. Always emits
     * the primary path traversal first, if applicable.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    
};

}

#endif
