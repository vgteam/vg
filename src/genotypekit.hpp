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
#include "bubbles.hpp"
#include "distributions.hpp"
#include "snarls.hpp"

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


// And now the implementations
    
class CactusUltrabubbleFinder : public SnarlFinder {
    
    // Holds the vg graph we are looking for sites in.
    VG& graph;
    
    // Use this path name as a rooting hint, if present.
    string hint_path_name;
    
public:
    /**
     * Make a new CactusSiteFinder to find sites in the given graph.
     */
    CactusUltrabubbleFinder(VG& graph, const string& hint_path_name);
    
    /**
     * Find all the sites in parallel with Cactus, make the site tree, and call
     * the given function on all the top-level sites.
     */
    virtual SnarlManager find_snarls();
    
};
    
class ExhaustiveTraversalFinder : TraversalFinder {
    
    VG& graph;
    SnarlManager& snarl_manager;
    
public:
    ExhaustiveTraversalFinder(VG& graph, SnarlManager& snarl_manager);
    
    virtual ~ExhaustiveTraversalFinder();
    
    /**
     * Exhaustively enumerate all traversals through the site. Only valid for
     * acyclic Snarls.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    
private:
    void stack_up_valid_walks(NodeTraversal walk_head, vector<NodeTraversal>& stack);
    
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

/**
 * TBD
 *
 */
//class StandardVcfRecordConverter {
//private:
//    const ReferenceIndex& index;
//    vcflib::VariantCallFile& vcf;
//    const string& sample_name;
//    
//public:
//    StandardVcfRecordConverter();
//    virtual ~StandardVcfRecordConverter() = default;
//    
//    virtual vcflib::Variant convert(const Locus& locus) = 0;
//};
    
    
}

#endif

