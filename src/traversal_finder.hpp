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
#include <limits>

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
#include "gbwt_helper.hpp"

namespace vg {

using namespace std;

class AugmentedGraph;

// some Protobuf replacements
using Traversal = vector<handle_t>;
using PathInterval = pair<step_handle_t, step_handle_t>;
string traversal_to_string(const PathHandleGraph* graph, const Traversal& traversal, int64_t max_steps = 10);
// replaces pb2json(snarl)
string graph_interval_to_string(const HandleGraph* graph, const handle_t& start_handle, const handle_t& end_handle);

/**
 * Represents a strategy for finding traversals of (nested) sites. Polymorphic
 * base class/interface.
 */
class TraversalFinder {
public:
    virtual ~TraversalFinder() = default;
    
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site) = 0;

    // new, protobuf-free interface. hope is to eventually deprecate the old one.  for now it is only supported in a few places
    virtual vector<Traversal> find_traversals(const handle_t& snarl_start, const handle_t& snarl_end) {
        assert(false); return{};
    }
};

class ExhaustiveTraversalFinder : public TraversalFinder {
    
    const HandleGraph& graph;
    SnarlManager& snarl_manager;
    bool include_reversing_traversals;
    
public:
    ExhaustiveTraversalFinder(const HandleGraph& graph, SnarlManager& snarl_manager,
                              bool include_reversing_traversals = false);
    
    virtual ~ExhaustiveTraversalFinder();
    
    /**
     * Exhaustively enumerate all traversals through the site. Only valid for
     * acyclic Snarls.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    
protected:
    void stack_up_valid_walks(handle_t walk_head, vector<Visit>& stack);
    virtual bool visit_next_node(handle_t handle) { return true; }
    void add_traversals(vector<SnarlTraversal>& traversals, handle_t traversal_start,
                        unordered_set<handle_t>& stop_at, unordered_set<handle_t>& yield_at);
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
 *
 * DEPRECATED: Use PathTraversalFinder instead
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
    virtual pair<vector<SnarlTraversal>, vector<string> > find_named_traversals(const Snarl& site);
    
};

class PathBasedTraversalFinder : public TraversalFinder{
    const PathHandleGraph& graph;
    SnarlManager& snarlmanager;
    public:
    PathBasedTraversalFinder(const PathHandleGraph& graph, SnarlManager& sm);
    virtual ~PathBasedTraversalFinder() = default;
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);

};

/** This is a Handle Graph replacement for PathRestrictedTraversalFinder
 * that uses the PathHandleGraph interface instead of the VG-based 
 * path index.  It returns all traversals through a snarl that are contained
 * within paths in the graph.  It can also return a mapping from the traversals
 * to their paths*/
class PathTraversalFinder : public TraversalFinder {
    
protected:
    // our graph with indexed path positions
    const PathHandleGraph& graph;
    
    // restrict to these paths
    unordered_set<path_handle_t> paths;
    
public:
    // if path_names not empty, only those paths will be considered
    PathTraversalFinder(const PathHandleGraph& graph,
                        const vector<string>& path_names = {});

    /**
     * Return all traversals through the site that are sub-paths of embedded paths in the graph
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    virtual vector<Traversal> find_traversals(const handle_t& snarl_start, const handle_t& snarl_end);

    /**
    * Like above, but return the path steps for the for the traversal endpoints
    */
    virtual pair<vector<SnarlTraversal>, vector<PathInterval>> find_path_traversals(const Snarl& site);
    virtual pair<vector<Traversal>, vector<PathInterval>> find_path_traversals(const handle_t& snarl_start, const handle_t& snarl_end);

};    
    

/**
 * This traversal finder finds one or more traversals through leaf sites with
 * no children. It uses a depth-first search. It doesn't work on non-leaf
 * sites, and is not guaranteed to find all traversals. Only works on acyclic
 * sites that are start-end-reachable.
 */
class TrivialTraversalFinder : public TraversalFinder {

    // Holds the vg graph we are looking for traversals in.
    const HandleGraph& graph;

public:
    TrivialTraversalFinder(const HandleGraph& graph);

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
    const PathHandleGraph& graph;

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
    pair<Support, vector<Visit>> find_bubble(id_t node, const edge_t* edge, const Snarl* snarl, PathIndex& index,
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

    /** 
     * Get support
     */
    bool has_supports = false;
    function<Support(id_t)> get_node_support;
    function<Support(edge_t)> get_edge_support;

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
    RepresentativeTraversalFinder(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                                  size_t max_depth, size_t max_width, size_t max_bubble_paths,
                                  size_t min_node_support = 1, size_t min_edge_support = 1,
                                  function<PathIndex*(const Snarl&)> get_index = [](const Snarl& s) { return nullptr; },
                                  function<Support(id_t)> get_node_support = nullptr,
                                  function<Support(edge_t)> get_edge_support = nullptr);
    
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
    const PathHandleGraph& graph;

    /// Use this to check if our snarl runs through a reference path
    /// (may be overkill, but can be used for sanity checking)
    PathTraversalFinder path_finder;
    
    /// The SnarlManager managiung the snarls we use
    SnarlManager& snarl_manager;

    /// Store variants indexed by an arbitrary node in one of their associated
    /// alt paths. We can then use this to find all variants in a top-level snarl
    unordered_map<id_t, list<vcflib::Variant*>> node_to_variant;

    /// Use this method to prune the search space by selecting alt-alleles
    /// to skip by considering their paths (in SnarlTraversal) format
    /// It will try again and again until enough traversals are pruned,
    /// with iteration keeping track of how many tries (so it should become stricter
    /// as iteration increases)
    function<bool(const SnarlTraversal& alt_path, int iteration)> skip_alt;

    /// If a snarl has more than this many traversals, return nothing and print
    /// a warning.  Dense and large deletions will make this happen from time
    /// to time.  In practice, skip_alt (above) can be used to prune down
    /// the search space by selecting alleles to ignore.
    size_t max_traversal_cutoff;

    /// Maximum number of pruning iterations
    size_t max_prune_iterations = 2;
    
    /// Include snarl endpoints in traversals
    bool include_endpoints = true;

    /// How far to scan when looking for deletions
    size_t max_deletion_scan_nodes = 50;

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
    VCFTraversalFinder(const PathHandleGraph& graph, SnarlManager& snarl_manager, vcflib::VariantCallFile& vcf,
                       const vector<string>& ref_path_names = {},
                       FastaReference* fasta_ref = nullptr,
                       FastaReference* ins_ref = nullptr,
                       function<bool(const SnarlTraversal&, int)> skip_alt = nullptr,
                       size_t max_traversal_cutoff = 50000);
        
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
                                    path_handle_t ref_path,
                                    step_handle_t start_step,
                                    step_handle_t end_step,
                                    vector<pair<SnarlTraversal, vector<int> > >& output_traversals);

    /** Get a traversal for a given haplotype.  It gets all the nodes and edges from the alt 
     * paths, and greedily walks over them whenever possible (traversing the reference otherwise).
     * if there is no traversal that can satisfy the haplotype, then the returned bool is set to false
     */
    pair<SnarlTraversal, bool> get_alt_traversal(const Snarl& site,
                                                 const vector<vcflib::Variant*>& site_variants,
                                                 path_handle_t ref_path,
                                                 step_handle_t start_step,
                                                 step_handle_t end_step,
                                                 const vector<int>& haplotype);

    /** Get a set of all alt path nodes and deletion edges for a halptype.
     */
    pair<unordered_set<handle_t>, unordered_set<pair<handle_t, handle_t> >>
         get_haplotype_alt_contents(const vector<vcflib::Variant*>& site_variants,
                                    const vector<int>& haplotype,
                                    path_handle_t ref_path);

    /** Get one alt-path out of the graph in the form of a snarl traversal.  if the path is a deletion,
     *  the edges corresponding to the deletion are also returned.  note that it is indeed possible
     *  for one alt path (and therefore one vcf alleles) to correspond to several deletion edges in the 
     *  graph due to normalization during construction.
     */
    pair<SnarlTraversal, vector<edge_t>> get_alt_path(vcflib::Variant* site_variant, int allele, path_handle_t ref_path);

    /**
     * An alt path for a deletion is the deleted reference path.  But sometimes vg construct doesn't
     * write a deletion edge that exactly jumps over the alt path.  In these cases, we need to 
     * search the graph for one. This does a brute-force check of all deletion edges in the vicinity
     * for one that's the same size as the one we're looking for.  
     * It tries to find a set of nearyby deletions that match the desired length.
     * Todo: check the sequence as well
     * Also todo: It'd be really nice if construct -fa would make the deletion-edge easily inferrable 
     * from the alt path.  It really shouldn't be necessary to hunt around. 
     * Returns: <deletion traversal, list of deletion edges>
     */
    pair<SnarlTraversal, vector<edge_t>> scan_for_deletion(vcflib::Variant* var, int allele, path_handle_t ref_path,
                                                           step_handle_t first_path_step, step_handle_t last_path_step);

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
                                               path_handle_t ref_path);

    /**
     * Count the possible traversal paths.  Return false if we ever get beyond our cutoff
     */
    bool check_max_trav_cutoff(const vector<vector<int> >& alleles);

    /**
     * Lookup a node in the reference path (mimics old PathIndex)
     */
    pair<step_handle_t, bool> step_in_path(handle_t handle, path_handle_t path_handle) const;
                                
};

/** Finds traversals with the most flow.  Node and edge weights are specified 
 * using the callbacks and can be used, ex, to yield read supports.  
 * If one traversal is requested, then the path with the highest flow (whose
 * node or edge with the minimum weight is maximum) is returned.  If K
 * traversals are specified, then the K highest flow traversals are returned.
 * This is designed to be a replacement for RepresentativeTraversalFinder.
 * It should do a better job of enumerating off-reference traversals, and will
 * of course guarantee to return all the optimal traversals (in the context of max flow).
 * Unlike RepresentativeTraversalFinder, it does not currently support nested
 * snarls, so all traversals returned are explicit.  
 * It is possible that it will blow up on massive snarls, espeically for large Ks. 
 */ 
class FlowTraversalFinder : public TraversalFinder {
    
protected:
    const HandleGraph& graph;
    
    SnarlManager& snarl_manager;

    /// The K-best traversals are returned
    size_t K;

    /// Callbacks to get supports
    function<double(handle_t)> node_weight_callback;
    function<double(edge_t)> edge_weight_callback;

    /// Call off the search as soon as a traversal of this length (bp) is encountered
    size_t max_traversal_length;
    
public:
    
    // if path_names not empty, only those paths will be considered
    // if a traversal is found that exceeds max_traversal_length, then the search is called off
    FlowTraversalFinder(const HandleGraph& graph, SnarlManager& snarl_manager,
                        size_t K,
                        function<double(handle_t)> node_weight_callback,
                        function<double(edge_t)> edge_weight_callback,
                        size_t max_traversal_length = numeric_limits<size_t>::max());

    /**
     * Return the K widest (most flow) traversals through the site
     * The reference traversal will be returned first (regardless of its flow).
     * After, the traversals are listed in decreasing order
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);

    /**
     * Return the K widest traversals, along with their flows
     */
    virtual pair<vector<SnarlTraversal>, vector<double>> find_weighted_traversals(const Snarl& site,
                                                                                  bool greedy_avg = false,
                                                                                  const HandleGraph* overlay = nullptr);

    /// Set K
    void setK(size_t k);

};    

/** Rerturn all traversals of a snarl that correspond to haplotypes stored in a GBWT
 */
class GBWTTraversalFinder : public TraversalFinder {

protected:
    
    const HandleGraph& graph;
    const gbwt::GBWT& gbwt;
    
public:
    
    GBWTTraversalFinder(const HandleGraph& graph, const gbwt::GBWT& gbwt);
    
    virtual ~GBWTTraversalFinder();

    /* Return a traversal for every gbwt thread through the snarl 
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    virtual vector<Traversal> find_traversals(const handle_t& snarl_start, const handle_t& snarl_end);

    /** Return the traversals, paired with their path identifiers in the gbwt.  The traversals are 
     *  unique, but there can be more than one path along each one (hence the vector)
     */
    virtual pair<vector<SnarlTraversal>, vector<vector<gbwt::size_type>>>
    find_gbwt_traversals(const Snarl& site, bool return_paths = true);
    virtual pair<vector<Traversal>, vector<vector<gbwt::size_type>>>
    find_gbwt_traversals(const handle_t& snarl_start, const handle_t& snarl_end, bool return_paths = true);

    /** Return traversals paired with path identifiers from the GBWT.  The traversals are *not* unique
     * (which is consistent with PathTraversalFinder)
     * To get the sample name from the path identifier id, use gbwtgraph::get_path_sample_name();
     */
    virtual pair<vector<SnarlTraversal>, vector<gbwt::size_type>> find_path_traversals(const Snarl& site);
    virtual pair<vector<Traversal>, vector<gbwt::size_type>> find_path_traversals(const handle_t& snarl_start, const handle_t& snarl_end);

    const gbwt::GBWT& get_gbwt() { return gbwt; }
    
protected:

    /**
     * Breadth first search from the start to the end, only branching if there's a haplotype 
     * in the GBWT, and returning all unique haplotypes found. 
     */
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > get_spanning_haplotypes(handle_t start, handle_t end);    
};

}

#endif
