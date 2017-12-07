#ifndef VG_GENOTYPEKIT_HPP_INCLUDED
#define VG_GENOTYPEKIT_HPP_INCLUDED
// genotypekit.hpp: defines pluggable modules for building the genotyper

// The basic idea here is we're going to create a few of these classes, fill in
// their public parameter fields, and then wire them up and set them going to
// emit genotypes in a streaming fashion.

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
#include "distributions.hpp"
#include "snarls.hpp"
#include "path_index.hpp"

extern "C" {
#include "sonLib.h"
}

namespace vg {

using namespace std;
    
/**
 * Represents a strategy for finding (nested) sites in a vg graph that can be described
 * by snarls. Polymorphic base class/interface.
 */
class SnarlFinder {
public:
    virtual ~SnarlFinder() = default;
    
    /**
     * Run a function on all root-level NestedSites in parallel. Site trees are
     * passed by value so they have a clear place to live during parallel
     * operations.
     */
    virtual SnarlManager find_snarls() = 0;
};

/**
 * Represents a strategy for finding traversals of (nested) sites. Polymorphic
 * base class/interface.
 */
class TraversalFinder {
public:
    virtual ~TraversalFinder() = default;
    
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site) = 0;
};


/**
 * Represents a strategy for computing consistency between Alignments and
 * SnarlTraversals. Determines whether a read is consistent with a SnarlTraversal
 * or not (but has access to all the SnarlTraversals). Polymorphic base
 * class/interface.
 */
class ConsistencyCalculator {
public:
    virtual ~ConsistencyCalculator() = default;
    
    /**
     * Return true or false for each tarversal of the site, depending on if the
     * read is consistent with it or not.
     */
    virtual vector<bool> calculate_consistency(const Snarl& site,
        const vector<SnarlTraversal>& traversals, const Alignment& read) const = 0;
};


/**
 * Represents a strategy for calculating Supports for SnarlTraversals.
 * Polymorphic base class/interface.
 */ 
class TraversalSupportCalculator {
public:
    virtual ~TraversalSupportCalculator() = default;
    
    /**
     * Return Supports for all the SnarlTraversals, given the reads and their
     * consistency flags.
     */
    virtual vector<Support> calculate_supports(const Snarl& site,
        const vector<SnarlTraversal>& traversals, const vector<Alignment*>& reads,
        const vector<vector<bool>>& consistencies) const = 0;
};



// TODO: This needs to be redesigned vis a vis the Genotype object. Genotypes
// need an accompanying Locus object in order to have the Path of the allele
// and also they are not site tree aware.
/**
 * Represents a strategy for calculating genotype likelihood for a (nested)
 * Site. Polymorphic base class/interface.
 */
class GenotypeLikelihoodCalculator {
public:
    virtual ~GenotypeLikelihoodCalculator() = default;
    
    /**
     * Return the log likelihood of the given genotype.
     */
    virtual double calculate_log_likelihood(const Snarl& site,
        const vector<SnarlTraversal>& traversals, const Genotype& genotype,
        const vector<vector<bool>>& consistencies, const vector<Support>& supports,
        const vector<Alignment*>& reads) = 0;
};

/**
 * Represents a strategy for assigning genotype priors. Polymorphic base
 * class/interface.
 */
class GenotypePriorCalculator {
public:
    virtual ~GenotypePriorCalculator() = default;
    
    /**
     * Return the log prior of the given genotype.
     *
     * TODO: ploidy priors on nested sites???
     */
    virtual double calculate_log_prior(const Genotype& genotype) = 0;
};

/**
 * Represents a strategy for converting Locus objects to VCF records.
 * Polymorphic base class/interface.
 */
class VcfRecordConverter {
public:
    virtual ~VcfRecordConverter() = default;
    
    virtual vcflib::Variant convert(const Locus& locus) = 0;
};


/**
 * Represents a filter that passes or rejects VCF records according to some
 * criteria. Polymorphic base class/interface.
 */
class VcfRecordFilter {
public:
    virtual ~VcfRecordFilter() = default;
    
    /**
     * Returns true if we should keep the given VCF record, and false otherwise.
     */
    virtual bool accept_record(const vcflib::Variant& variant) = 0;
};



/////////////////////////////////
// And now the implementations //
/////////////////////////////////

/// General interface for an augmented graph.  This is a graph that was constructed
/// by adding some read information to an original ("base") graph.  We preserve
/// mappings back to the base graph via translations. Augmented graphs can be
/// made using edit (such as in vg mod -i) or pileup.
/// Todo : further abstract to handle graph interface
struct AugmentedGraph {
    // This holds all the new nodes and edges
    VG graph;

    // This holds the base graph (only required for mapping edges)
    VG* base_graph = NULL;

    // Translations back to the base graph
    Translator translator;
    
    // Map an edge back to the base graph (return NULL if does not exist)
    // In the special case that an edge in the augmented graph represents
    // two abutting positions within a node in the base graph, null, true
    // is returned (todo: less tortured interface?)
    pair<const Edge*, bool> base_edge(const Edge* augmented_edge);

    // Is this node novel?  ie does it not map back to the base graph?
    bool is_novel_node(const Node* augmented_node) {
        Position pos;
        pos.set_node_id(augmented_node->id());
        return !translator.has_translation(pos);
    }
    
    // Is this edge novel?  ie does it not map back to the base graph?
    bool is_novel_edge(const Edge* augmented_edge) {
        auto be_ret = base_edge(augmented_edge);
        return be_ret.first == NULL && be_ret.second == false;
    }
    
    /**
     * Clear the contents.
     */
    virtual void clear();

    /** 
     * Construct an augmented graph using edit() on a set of alignments. Adds
     * the paths of the alignments to the graph.
     * 
     * If unique_names is set, makes sure all the alignments' paths have unique
     * names, which is a requirement for paths in a graph.
     *
     * If leave_edits is set, the alignment's paths are not modified. The
     * alignments' paths will be the paths they were originally aligned to,
     * although the graph will be modified to describe the edits that the
     * alignments found.
     */
    void augment_from_alignment_edits(vector<Alignment>& alignments, bool unique_names = true,
                                      bool leave_edits = false);

    /**
     * Load the translations from a file
     */
    void load_translations(istream& in_file);

    /**
     * Write the translations to a file
     */
    void write_translations(ostream& out_file);
};

/// Augmented Graph that holds some Support annotation data specific to vg call
struct SupportAugmentedGraph : public AugmentedGraph {
        
    // This holds support info for nodes. Note that we discard the "os" other
    // support field from StrandSupport.
    // Supports for nodes are minimum distinct reads that use the node.
    map<Node*, Support> node_supports;
    // And for edges
    map<Edge*, Support> edge_supports;
    
    /**
     * Return true if we have support information, and false otherwise.
     */
    bool has_supports();
    
    /**
     * Get the Support for a given Node, or 0 if it has no recorded support.
     */
    Support get_support(Node* node);
    
    /**
     * Get the Support for a given Edge, or 0 if it has no recorded support.
     */
    Support get_support(Edge* edge);    
    
    /**
     * Clear the contents.
     */
    virtual void clear();

    /**
    * Read the supports from protobuf.
    */
    void load_supports(istream& in_file);

    /**
     * Write the supports to protobuf
     */
    void write_supports(ostream& out_file);
    
};
    


class SimpleConsistencyCalculator : public ConsistencyCalculator{
    public:
    ~SimpleConsistencyCalculator();
    vector<bool> calculate_consistency(const Snarl& site,
        const vector<SnarlTraversal>& traversals, const Alignment& read) const;
};

/**
 * Class for finding all snarls using the base-level Cactus snarl decomposition
 * interface.
 */
class CactusSnarlFinder : public SnarlFinder {
    
    /// Holds the vg graph we are looking for sites in.
    VG& graph;
    
    /// Holds the names of reference path hints
    unordered_set<string> hint_paths;
    
    /// Create a snarl in the given SnarlManager with the given start and end,
    /// containing the given child snarls in the list of chains of children and
    /// the given list of unary children. Recursively creates snarls in the
    /// SnarlManager for the children. Returns a pointer to the finished snarl
    /// in the SnarlManager. Start and end may be empty visits, in which case no
    /// snarl is created, all the child chains are added as root chains, and
    /// null is returned. If parent_start and parent_end are empty Visits, no
    /// parent() is added to the produced snarl.
    const Snarl* recursively_emit_snarls(const Visit& start, const Visit& end,
        const Visit& parent_start, const Visit& parent_end,
        stList* chains_list, stList* unary_snarls_list, SnarlManager& destination);
    
public:
    /**
     * Make a new CactusSnarlFinder to find snarls in the given graph.
     * We can't filter trivial bubbles because that would break our chains.
     *
     * Optionally takes a hint path name.
     */
    CactusSnarlFinder(VG& graph);
    
    /**
     * Make a new CactusSnarlFinder with a single hinted path to base the
     * decomposition on.
     */
    CactusSnarlFinder(VG& graph, const string& hint_path);
    
    /**
     * Find all the snarls with Cactus, and put them into a SnarlManager.
     */
    virtual SnarlManager find_snarls();
    
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

/**
 * This genotype prior calculator has a fixed prior for homozygous genotypes and
 * a fixed prior for hets.
 */
class FixedGenotypePriorCalculator : public GenotypePriorCalculator {
public:
    // These parameters are configurable, but have defaults.
    double homozygous_prior_ln = prob_to_logprob(0.999);
    double heterozygous_prior_ln = prob_to_logprob(0.001);


    virtual ~FixedGenotypePriorCalculator() = default;
    virtual double calculate_log_prior(const Genotype& genotype);
};


class SimpleTraversalSupportCalculator : public TraversalSupportCalculator{
    // A set of traversals through the site
    // A set of alignments to the site
    // And a set of consistencies, one vector for each alignment,
    //    one boolean per traversal.
    public:
    ~SimpleTraversalSupportCalculator();
    vector<Support> calculate_supports(const Snarl& site,
        const vector<SnarlTraversal>& traversals, const vector<Alignment*>& reads,
        const vector<vector<bool>>& consistencies) const;
};

/**
 * TBD
 *
 */
// class StandardVcfRecordConverter {
// private:
//    const ReferenceIndex& index;
//    vcflib::VariantCallFile& vcf;
//    const string& sample_name;
   
// public:
//    StandardVcfRecordConverter();
//    virtual ~StandardVcfRecordConverter() = default;
   
//    virtual vcflib::Variant convert(const Locus& locus) = 0;
// };
    
// We also supply utility functions for working with genotyping Protobuf objects

/**
 * Create a Support for the given forward and reverse coverage and quality.
 */
Support make_support(double forward, double reverse, double quality);

/**
 * Get the total read support in a Support.
 */
double total(const Support& support);

/**
 * Get the minimum support of a pair of Supports, by taking the min in each
 * orientation.
 */
Support support_min(const Support& a, const Support& b);

/**
 * Get the maximum support of a pair of Supports, by taking the max in each
 * orientation.
 */
Support support_max(const Support& a, const Support& b);

/**
 * Add two Support values together, accounting for strand.
 */
Support operator+(const Support& one, const Support& other);

/**
 * Add in a Support to another.
 */
Support& operator+=(Support& one, const Support& other);

/**
 * Scale a Support by a factor.
 */
template<typename Scalar>
Support operator*(const Support& support, const Scalar& scale) {
    Support prod;
    prod.set_forward(support.forward() * scale);
    prod.set_reverse(support.reverse() * scale);
    prod.set_left(support.left() * scale);
    prod.set_right(support.right() * scale);
    
    // log-scaled quality can just be multiplied
    prod.set_quality(support.quality() * scale);
    
    return prod;
}

/**
 * Scale a Support by a factor, in place.
 */
template<typename Scalar>
Support& operator*=(Support& support, const Scalar& scale) {
    support.set_forward(support.forward() * scale);
    support.set_reverse(support.reverse() * scale);
    support.set_left(support.left() * scale);
    support.set_right(support.right() * scale);
    
    // log-scaled quality can just be multiplied
    support.set_quality(support.quality() * scale);
    
    return support;
}

/**
 * Scale a Support by a factor, the other way
 */
template<typename Scalar>
Support operator*(const Scalar& scale, const Support& support) {
    Support prod;
    prod.set_forward(support.forward() * scale);
    prod.set_reverse(support.reverse() * scale);
    prod.set_left(support.left() * scale);
    prod.set_right(support.right() * scale);
    
    // log-scaled quality can just be multiplied
    prod.set_quality(support.quality() * scale);
    
    return prod;
}

/**
 * Divide a Support by a factor.
 */
template<typename Scalar>
Support operator/(const Support& support, const Scalar& scale) {
    Support scaled;
    
    scaled.set_forward(support.forward() / scale);
    scaled.set_reverse(support.reverse() / scale);
    scaled.set_left(support.left() / scale);
    scaled.set_right(support.right() / scale);
    
    // log-scaled quality can just be divided. Maybe.
    scaled.set_quality(support.quality() / scale);
    
    return scaled;
}

/**
 * Divide a Support by a factor, in place.
 */
template<typename Scalar>
Support& operator/=(Support& support, const Scalar& scale) {
    support.set_forward(support.forward() / scale);
    support.set_reverse(support.reverse() / scale);
    support.set_left(support.left() / scale);
    support.set_right(support.right() / scale);
    
    // log-scaled quality can just be divided
    support.set_quality(support.quality() / scale);
    
    return support;
}

/**
 * Support less-than, based on total coverage.
 */
bool operator< (const Support& a, const Support& b);

/**
 * Support greater-than, based on total coverage.
 */
bool operator> (const Support& a, const Support& b);

/**
 * Allow printing a Support.
 */
ostream& operator<<(ostream& stream, const Support& support);

/**
 * Get a VCF-style 1/2, 1|2|3, etc. string from a Genotype.
 */
string to_vcf_genotype(const Genotype& gt);
    
}

#endif

