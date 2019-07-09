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

#include <structures/immutable_list.hpp>

#include <vg/vg.pb.h>
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
    
protected:
    void stack_up_valid_walks(NodeTraversal walk_head, vector<NodeTraversal>& stack);
    virtual bool visit_next_node(const Node*, const Edge*) { return true; }
    void add_traversals(vector<SnarlTraversal>& traversals, NodeTraversal traversal_start,
                        set<NodeTraversal>& stop_at, set<NodeTraversal>& yield_at);
};

/** Does exhaustive traversal, but restricting to nodes and edges that meet 
    support thresholds (counts of reads that touch them, taken from augmented graph).
*/
class SupportRestrictedTraversalFinder : public ExhaustiveTraversalFinder {
public:
    AugmentedGraph& aug;
    int min_node_support;
    int min_edge_support;

    SupportRestrictedTraversalFinder(AugmentedGraph& augmented_graph,
                                     SnarlManager& snarl_manager,
                                     int min_node_support = 1,                                     
                                     int min_edge_support = 1,
                                     bool include_reversing_traversals = false);
    virtual ~SupportRestrictedTraversalFinder();
protected:
    virtual bool visit_next_node(const Node*, const Edge*);
};

class ReadRestrictedTraversalFinder : public TraversalFinder {

    AugmentedGraph& aug;
    SnarlManager& snarl_manager;
    
    // How many times must a path recur before we try aligning to it? Also, how
    // many times must a node in the graph be visited before we use it in indel
    // realignment for nearby indels? Note that the primary path counts as a
    // recurrence. TODO: novel inserts can't recur, and novel deletions can't be
    // filtered in this way.
    int min_recurrence;
    
    // How many nodes max should we walk when checking if a path runs through a superbubble/site
    int max_path_search_steps;
    
public:
    ReadRestrictedTraversalFinder(AugmentedGraph& augmented_graph, SnarlManager& snarl_manager,
                                  int min_recurrence = 2,
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

/**
 * Like ReadRestrictedTraversal finder, but works on paths in the graph.  As with
 * the former, it's been cut out of Genotyper and moved here.  
 * 
 * I'm not sure what PathBasedTraversalFinder (see below) does, but it does not work
 * as a drop-in replacement for this class, so keep the two implementations at least 
 * for now.    
 */
class PathRestrictedTraversalFinder : public TraversalFinder {

    VG& graph;
    SnarlManager& snarl_manager;

    // store the reads that are embedded in the augmented graph, by their unique names
    map<string, const Alignment*>& reads_by_name;

    // How many times must a path recur before we try aligning to it? Also, how
    // many times must a node in the graph be visited before we use it in indel
    // realignment for nearby indels? Note that the primary path counts as a
    // recurrence. TODO: novel inserts can't recur, and novel deletions can't be
    // filtered in this way.
    int min_recurrence;
    
    // How many nodes max should we walk when checking if a path runs through a superbubble/site
    int max_path_search_steps;

    // Allow multiple traversals with the same sequence
    bool allow_duplicates;
    
public:
    PathRestrictedTraversalFinder(VG& graph, SnarlManager& snarl_manager,
                                  map<string, const Alignment*>& reads_by_name,
                                  int min_recurrence = 2,
                                  int max_path_search_steps = 100,
                                  bool allow_duplicates = false);
    
    virtual ~PathRestrictedTraversalFinder();

    /**
     * For the given site, emit all subpaths with unique sequences that run from
     * start to end, out of the paths in the graph. Uses the map of reads by
     * name to determine if a path is a read or a real named path. Paths through
     * the snarl supported only by reads are subject to a min recurrence count,
     * while those supported by actual embedded named paths are not.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);

   /**
    * Like above, but return the path name corresponding to each traversal
    */
    virtual pair<vector<SnarlTraversal>, vector<string>> find_named_traversals(const Snarl& site);
    
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
 * This traversal finder finds one or more traversals through leaf sites with
 * no children. It uses a depth-first search. It doesn't work on non-leaf
 * sites, and is not guaranteed to find all traversals. Only works on acyclic
 * sites that are start-end-reachable.
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
    AugmentedGraph& augmented;
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
    /// Minimum support for a node to consider travnersal through it
    size_t min_node_support;
    /// Minimum support for a edge to consider travnersal through it
    size_t min_edge_support;    
    
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
     * Get the minimum support of all nodes and edges in path, in the path's forward orientation.
     */
    Support min_support_in_path(const list<Visit>& path);
        
    /**
     * Do a breadth-first search left from the given node traversal, and return
     * lengths, target-path-relative orientations, and paths starting at the
     * given node and ending on the given indexed path. Refuses to visit nodes
     * with no support, if support data is available in the augmented graph.
     *
     * If in_snarl is not null, restricts the found paths to stay within the
     * given snarl.
     *
     * If both_orientations_distance is not zero, keeps searching up to that
     * many steps after finding the target path to see if it can find a node on
     * the target path in the opposite orientation. This is useful for
     * inversions.
     */
    set<tuple<size_t, bool, structures::ImmutableList<Visit>>> bfs_left(Visit visit, PathIndex& index,
                                                                        bool stop_if_visited = false,
                                                                        const Snarl* in_snarl = nullptr,
                                                                        size_t both_orientations_distance = 0);
        
    /**
     * Do a breadth-first search right from the given node traversal, and return
     * lengths, target-path-relative orientations, and paths starting at the
     * given node and ending on the given indexed path. Refuses to visit nodes
     * with no support, if support data is available in the augmented graph.
     *
     * API is similar to bfs_left().
     */
    set<tuple<size_t, bool, structures::ImmutableList<Visit>>> bfs_right(Visit visit, PathIndex& index,
                                                                         bool stop_if_visited = false,
                                                                         const Snarl* in_snarl = nullptr,
                                                                         size_t both_orientations_distance = 0);
        
    /**
     * Get the length of a path through nodes, in base pairs.
     */
    size_t bp_length(const  structures::ImmutableList<Visit>& path);

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
    RepresentativeTraversalFinder(AugmentedGraph& augmented, SnarlManager& snarl_manager,
        size_t max_depth, size_t max_width, size_t max_bubble_paths,
        size_t min_node_support = 1, size_t min_edge_support = 1,
        function<PathIndex*(const Snarl&)> get_index = [](const Snarl& s) { return nullptr; });
    
    /// Should we emit verbose debugging info?
    bool verbose = false;

    /// Should trivial child snarls have their traversals glommed into ours?
    bool eat_trivial_children = false;
    
    /// What timeout/step limit should we use for finding other orientations of
    /// the reference path after we find one?
    size_t other_orientation_timeout = 10;
    
    virtual ~RepresentativeTraversalFinder() = default;
    
    /**
     * Find traversals to cover the nodes and edges of the snarl. Always emits
     * the primary path traversal first, if applicable.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    
};


/**
 * This TraversalFinder returns a traversals and their corresponding genotypes
 * from an input vcf. It relies on alt-paths in the graph (via construct -a)
 * to map between the vcf and the graph.  
 */
class VCFTraversalFinder : public TraversalFinder {

protected:
    VG& graph;
    
    /// The SnarlManager managiung the snarls we use
    SnarlManager& snarl_manager;

    /// We keep this around from the RepresentativeTraversalFinder for now:
    /// We use a path index for accessing the reference path
    /// We keep around a function that can be used to get an index for the
    /// appropriate path to use to scaffold a given site, or null if no
    /// appropriate index exists.
    function<PathIndex*(const Snarl&)> get_index;

    /// Store variants indexed by an arbitrary node in one of their associated
    /// alt paths. We can then use this to find all variants in a top-level snarl
    unordered_map<id_t, list<vcflib::Variant*>> node_to_variant;

    /// Use this method to prune the search space by selecting alt-alleles
    /// to skip by considering their paths (in SnarlTraversal) format
    function<bool(const SnarlTraversal& alt_path)> skip_alt;

    /// If a snarl has more than this many traversals, return nothing and print
    /// a warning.  Dense and large deletions will make this happen from time
    /// to time.  In practice, skip_alt (above) can be used to prune down
    /// the search space by selecting alleles to ignore.
    size_t max_traversal_cutoff;
    
    /// Include snarl endpoints in traversals
    bool include_endpoints = true;

    /// How far to scan when looking for deletions
    size_t max_deletion_scan_nodes = 100;

public:

    /**
     * Make a new VCFTraversalFinder.  Builds the indexes needed to find all the 
     * variants in a site.
     *
     * The skip_alt() method is defined, it is run on the alt-path of each variant
     * allele in the snarl.  If it returns true, that alt-path will never be included
     * in any traversals returned in find_traversals().  
     * This is used to, for example, use read support to prune the number of traversals
     * that are enumerated.
     */
    VCFTraversalFinder(VG& graph, SnarlManager& snarl_manager, vcflib::VariantCallFile& vcf,
                       function<PathIndex*(const Snarl&)> get_index,
                       FastaReference* fasta_ref = nullptr,
                       FastaReference* ins_ref = nullptr,
                       function<bool(const SnarlTraversal&)> skip_alt = nullptr,
                       size_t max_traversal_cutoff = 500000);
        
    virtual ~VCFTraversalFinder();
    
    /**
     * Find traversals for the site.  Each traversa is returned in a pair with
     * its haplotype.  The haplotype refers to the list of variants (also returned)
     */
    pair<vector<pair<SnarlTraversal, vector<int>>>, vector<vcflib::Variant*>> find_allele_traversals(Snarl site);

    /**
     * Return a list of traversals for the site.  The same traversals as above, only the
     * haplotype information not included
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);

    /** 
     * Get all the variants that are contained in a site */
    vector<vcflib::Variant*> get_variants_in_site(const Snarl& site);

    
protected:

    /** Load up all the variants into our node index
     */
    void create_variant_index(vcflib::VariantCallFile& vcf, FastaReference* ref_fasta = nullptr,
                              FastaReference* ins_fasta = nullptr);
    void delete_variant_index();

    /** Get a traversal for every possible haplotype (but reference)
     * in the most naive way possibe.  This will blow up terribly for sites that contain more than a few
     * variants.  There's an obvious dynamic programming speedup, but the main issue is that
     * the output size is exponential in the number of variants.
     */
    void brute_force_alt_traversals(const Snarl& site,
                                    const vector<vcflib::Variant*>& site_variants,
                                    PathIndex* path_index,
                                    PathIndex::iterator start_it,
                                    PathIndex::iterator end_it,
                                    vector<pair<SnarlTraversal, vector<int> > >& output_traversals);

    /** Get a traversal for a given haplotype.  It gets all the nodes and edges from the alt 
     * paths, and greedily walks over them whenever possible (traversing the reference otherwise).
     * if there is no traversal that can satisfy the haplotype, then the returned bool is set to false
     */
    pair<SnarlTraversal, bool> get_alt_traversal(const Snarl& site,
                                                 const vector<vcflib::Variant*>& site_variants,
                                                 PathIndex* path_index,
                                                 PathIndex::iterator start_it,
                                                 PathIndex::iterator end_it,
                                                 const vector<int>& haplotype);

    /** Get a set of all alt path nodes and deletion edges for a halptype.
     */
    pair<unordered_set<handle_t>, unordered_set<pair<handle_t, handle_t> >>
         get_haplotype_alt_contents(const vector<vcflib::Variant*>& site_variants,
                                    const vector<int>& haplotype,
                                    PathIndex* path_index);

    /** Get one alt-path out of the graph in the form of a snarl traversal.  bool value
     *  returned is true if allele is a deletion, and returned traversal will be a pair of nodes
     *  representing the deletion edge.
     */
    pair<SnarlTraversal, bool> get_alt_path(vcflib::Variant* site_variant, int allele, PathIndex* path_index);

    /**
     * An alt path for a deletion is the deleted reference path.  But sometimes vg construct doesn't
     * write a deletion edge that exactly jumps over the alt path.  In these cases, we need to 
     * search the graph for one. This does a brute-force check of all deletion edges in the vicinity
     * for one that's the same size as the one we're looking for.  It picks the nearest one and
     * returns true if the exact size is found.  
     * Todo: check the sequence as well
     * Also todo: It'd be really nice if construct -fa would make the deletion-edge easily inferrable 
     * from the alt path.  It really shouldn't be necessary to hunt around. 
     * Returns: <size delta between returned deletion and vcf allele,
     *           position delta between returned deletion and alt path>
     */
    pair<int, int> scan_for_deletion(vcflib::Variant* var, int allele, PathIndex* path_index,
                                     PathIndex::iterator& first_path_it, PathIndex::iterator& last_path_it);

    /**
     * Prune our search space using the skip_alt method.  Will return a list of pruned VCF alleles/
     *
     * ex, if the input has A --> T
     *                      G --> C,A
     * there input alleles are <0,1>, <0,1,2>.  If there's no support for the G->C on the second one,
     * the output would be <0,1>, <0,2>.
     *
     */
    vector<vector<int>> get_pruned_alt_alleles(const Snarl& site,
                                               const vector<vcflib::Variant*>& site_variants,
                                               PathIndex* path_index);

    /**
     * Count the possible traversal paths.  Return false if we ever get beyond our cutoff
     */
    bool check_max_trav_cutoff(const vector<vector<int> >& alleles);
                                
};

}

#endif
