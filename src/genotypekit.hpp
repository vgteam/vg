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
    // We just act like a vector of these, which are oriented traversals of
    // either nodes or NestedSites.
    struct Visit {
        Node* node = nullptr;
        NestedSite* child = nullptr;
        bool backward = false;
    };
    
    // We keep a list of these, which we could potentially splice up to flatten
    // down to just nodes.
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
    // And an index from start or end node pointer to the child that owns that
    // node (as an index in the above vector).
    map<Node*, size_t> child_border_node_index;
    
    // It also contains pointers to the nodes in it but not in its children.
    set<Node*> nodes;
    // And to the edges in it but not in its children. Edges attaching to nodes
    // not in the node set are attaching to start or end nodes of children,
    // which can be mapped to and from children using child_border_node_index
    // above.
    set<Edge*> edges;
    
    // And its oriented start and end anchoring node traversals.
    NodeTraversal start;
    NodeTraversal end;
    
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
     * Run a function on all root-level NestedSites in parallel.
     */
    virtual void for_each_site_parallel(const function<void(const NestedSite&)>& lambda) = 0;
};

/**
 * Represents a strategy for finding traversals of (nested) sites. Polymorphic
 * base class/interface.
 */
class TraversalFinder {
public:
    virtual ~TraversalFinder() = default;
    
    virtual vector<SiteTraversal> find_traversals(const NestedSite& site);
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

