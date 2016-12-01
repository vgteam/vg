#ifndef VG_GENOTYPEKIT_H
#define VG_GENOTYPEKIT_H
// genotypekit.hpp: defines pluggable modules for building the genotyper

// The basic idea here is we're going to create a few of these classes, fill in
// their public parameter fields, and then wire them up and set them going to
// emit genotypes in a streaming fashion.

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
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

namespace vg {

using namespace std;

// First we need the types we operate on

// Forward declaration
struct NestedSite;

/**
 * Represents a traversal of a (possibly nested) site, going from start to end
 * and visiting nodes, edges, and contained nested sites. Basic component of a
 * genotype.
 */
struct SiteTraversal {
    // We just act like a list of these, which are oriented traversals of
    // either nodes or NestedSites.
    struct Visit {
        Node* node = nullptr;
        NestedSite* child = nullptr;
        // Indicates:
        //   if node != nullptr  : reverse complement of node
        //   if child != nullptr : traversal of child site entering backwards through
        //                          end and leaving backwards through start
        bool backward = false;
        
        /**
         * Make a Visit form a node and an orientation
         */
        inline Visit(Node* node, bool backward = false) : node(node), backward(backward) {
            // Nothing to do!
        }
        
        /**
         * Make a Visit from a child site and an orientation.
         */
        inline Visit(NestedSite* child, bool backward = false) : child(child), backward(backward) {
            // Nothing to do!
        }
        
        /**
         * Make a Visit from a NodeTraversal.
         */
        inline Visit(const NodeTraversal& traversal) : node(traversal.node), backward(traversal.backward) {
            // Nothing to do!
        }
        
        
    };
    
    // We keep a list of these, which we could potentially splice up to flatten
    // down to just nodes.
    //
    // NB: If a Visit contains a child site, the list SHOULD NOT contain a Visit for
    // the node where we entered the NestedSite, but it SHOULD contain a Visit for the
    // node where we leave it (this simplifies the case where two NestedSites are
    // adjacent to each other).
    list<Visit> visits;
};

/**
 * Represents a genotypeable site, with input and output NodeTraversals, that
 * can contain other nested sites within it.
 *
 * Must be understood in relation to some vg graph.
 */
struct NestedSite {
    // Each NestedSite contains all its child sites.
    vector<NestedSite> children;
    // And an index from node traversal to the child that you read into
    // following that traversal (as an index in the above vector). These are
    // traversals of the start and end nodes of children, and the nodes they
    // refer to are conatined within this site's children, rather than this site
    // itself.
    map<NodeTraversal, size_t> child_border_index;
    
    // It also contains pointers to the nodes in it but not in its children.
    set<Node*> nodes;
    // And to the edges in it but not in its children. Edges attaching to nodes
    // not in the node set are attaching to start or end nodes of children,
    // which can be mapped to and from children using child_border_index
    // above.
    set<Edge*> edges;
    
    // And its oriented start and end anchoring node traversals.
    NodeTraversal start; // points into site
    NodeTraversal end; // points out of site
    
    inline bool operator==(const NestedSite& other) const {
        return start == other.start && end == other.end;
    }
    inline bool operator<(const NestedSite& other) const {
        return start < other.start ? true : (end < other.end ? true : false);
    }
    inline bool operator>(const NestedSite& other) const {
        return !(*this < other);
    }
};

// For genotypes we use the protobuf Genotype object, in the context of the
// traversals obtained for a given NestedSite.


// Then the actual pluggable things. We keep the graph and the overall set of
// reads as contextual things we pass in the constructors of the actual
// instances.

/**
 * Represents a strategy for finding (nested) Sites in a vg graph. Polymorphic
 * base class/interface.
 */
class SiteFinder {
public:
    virtual ~SiteFinder() = default;
    
    /**
     * Run a function on all root-level NestedSites in parallel. Site trees are
     * passed by value so they have a clear place to live during parallel
     * operations.
     */
    virtual void for_each_site_parallel(const function<void(NestedSite)>& lambda) = 0;
};

/**
 * Represents a strategy for finding traversals of (nested) sites. Polymorphic
 * base class/interface.
 */
class TraversalFinder {
public:
    virtual ~TraversalFinder() = default;
    
    virtual vector<SiteTraversal> find_traversals(const NestedSite& site) = 0;
};

/**
 * Represents a strategy for computing consistency between Alignments and
 * SiteTraversals. Determines whether a read is consistent with a SiteTraversal
 * or not (but has access to all the SiteTraversals). Polymorphic base
 * class/interface.
 */
class ConsistencyCalculator {
public:
    virtual ~ConsistencyCalculator() = default;
    
    /**
     * Return true or false for each tarversal of the site, depending on if the
     * read is consistent with it or not.
     */
    virtual vector<bool> calculate_consistency(const NestedSite& site, 
        const vector<SiteTraversal>& traversals, const Alignment& read) const = 0;
};

/**
 * Represents a strategy for calculating Supports for SiteTraversals.
 * Polymorphic base class/interface.
 */ 
class TraversalSupportCalculator {
public:
    virtual ~TraversalSupportCalculator() = default;
    
    /**
     * Return Supports for all the SiteTraversals, given the reads and their
     * consistency flags.
     */
    virtual vector<Support> calculate_supports(const NestedSite& site, 
        const vector<SiteTraversal>& traversals, const vector<Alignment*>& reads,
        const vector<vector<bool>>& consistencies) const = 0;
};

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
    virtual double calculate_log_likelihood(const NestedSite& site, 
        const vector<SiteTraversal>& traversals, const Genotype& genotype,
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


// And now the implementations

/**
 * This site finder finds sites with Cactus.
 */
class CactusSiteFinder : public SiteFinder {

    // Holds the vg graph we are looking for sites in.
    VG& graph;

    // Use this path name as a rooting hint, if present.
    string hint_path_name;

public:
    /**
     * Make a new CactusSiteFinder to find sites in the given graph.
     */
    CactusSiteFinder(VG& graph, const string& hint_path_name);
    
    virtual ~CactusSiteFinder() = default;
    
    /**
     * Find all the sites in parallel with Cactus, make the site tree, and call
     * the given function on all the top-level sites.
     */
    virtual void for_each_site_parallel(const function<void(NestedSite)>& lambda);

};
    
class ExhaustiveTraversalFinder : TraversalFinder {
    
    VG& graph;
    
public:
    ExhaustiveTraversalFinder(VG& graph);
    
    virtual ~ExhaustiveTraversalFinder() = default;
    
    /**
     * Exhaustively enumerate all traversals through the site
     *
     * If a traversal includes a NestedSite, the node traversal that enters
     * the site will not be included (instead it will a site traversal), but 
     * the node traversal that leaves the site will be included, unless that
     * node traversal is also enters another site
     */
    virtual vector<SiteTraversal> find_traversals(const NestedSite& site);
    
};
    
class ReadRestrictedTraversalFinder : TraversalFinder {
    
    VG& graph;
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
    ReadRestrictedTraversalFinder(VG& graph, const map<string, Alignment*>& reads_by_name
                                  int min_recurrence = 2, int max_path_search_steps = 100);
    
    virtual ~ReadRestrictedTraversalFinder() = default;
    
    /**
     * For the given site, emit all traversals with unique sequences that run from
     * start to end, out of the paths in the graph. Uses the map of reads by
     * name to determine if a path is a read or a real named path. Paths through
     * the site supported only by reads are subject to a min recurrence count,
     * while those supported by actual embedded named paths are not.
     */
    virtual vector<SiteTraversal> find_traversals(const NestedSite& site);
    
};

/**
 * This traversal finder finds one or more traversals through leaf sites with no
 * children. It uses a depth-first search. It doesn't work on non-leaf sites,
 * and is not guaranteed to find all traversals.
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
    virtual vector<SiteTraversal> find_traversals(const NestedSite& site);
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


}

#endif

