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
#include <vg/vg.pb.h>
#include "vg.hpp"
#include "translator.hpp"
#include "hash_map.hpp"
#include "types.hpp"
#include "statistics.hpp"
#include "snarls.hpp"
#include "path_index.hpp"
#include "packer.hpp"

namespace vg {

using namespace std;

/**
 * Given a path (which may run either direction through a snarl, or not touch
 * the ends at all), collect a list of NodeTraversals in order for the part
 * of the path that is inside the snarl, in the same orientation as the path.
 */
SnarlTraversal get_traversal_of_snarl(VG& graph, const Snarl* snarl, const SnarlManager& manager, const Path& path);
    
/**
 * Make a SnarlTraversal into the string it represents, including notes for nested child snarls.
 */
string traversal_to_string(VG& graph, const SnarlTraversal& path);

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
    /// This holds all the new nodes and edges
    VG graph;

    /// This holds the base graph (only required for mapping edges)
    VG* base_graph = nullptr;

    /// Translations back to the base graph
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

    // Do we have a base graph in order to run the above methods?
    bool has_base_graph() const {
        return base_graph != nullptr;
    }
    
    /// Get the alignments, if any, embedded in the graph that touch the given
    /// node ID.
    vector<const Alignment*> get_alignments(id_t node_id) const;
    // Or a given edge
    vector<const Alignment*> get_alignments(pair<NodeSide, NodeSide>) const;
    
    /// Get all the embedded alignments.
    vector<const Alignment*> get_alignments() const;

    /**
     * Get the Support for a given Node, or 0 if it has no recorded support.
     * (only forward strand)
     */
    virtual Support get_support(id_t node);
    
    /**
     * Get the Support for a given Edge, or 0 if it has no recorded support.
     * (only forward strand)
     */
    virtual Support get_support(edge_t edge);    

    virtual bool has_supports() const;
    
    /**
     * Clear the contents.
     */
    virtual void clear();

    /** 
     * Construct an augmented graph using edit() on a set of alignments.
     * Modifies the passed-in vector of alignments arbitrarily (in particular,
     * it may be empty after the call returns).
     *
     * Stores a modified version of the reads in the AugmentedGraph. Read
     * alignments that touch a node in the augmented graph can be retrieved with
     * get_alignments(). Reads will have softclips removed.
     *
     * Must only be called ONCE, because all modifications to the graph have
     * to be processed together to update the embedded alignment indexes.
     * 
     * If unique_names is set, makes sure all the alignments have unique names.
     *
     * If leave_edits is set, the alignment's paths are not modified after
     * trimming softclips. The alignments' paths will be the paths they were
     * originally aligned to, although the graph will be modified to describe
     * the edits that the alignments found.
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
    
protected:

    /// Holds all of the alignments that have been embedded in the graph. They
    /// aren't in the graph's paths object, but they are fully embedded without
    /// any edits relative to the graph's sequence.
    vector<Alignment> embedded_alignments;
    
    // Maps from node ID to the alignments that touch that node.
    unordered_map<id_t, vector<Alignment*>> alignments_by_node;

    // Maps from the edge to the alignments that touch that node.
    pair_hash_map<pair<NodeSide, NodeSide>, vector<Alignment*>> alignments_by_edge;
};

/// Augmented Graph that holds some Support annotation data specific to vg call
struct SupportAugmentedGraph : public AugmentedGraph {
        
    // This holds support info for nodes. Note that we discard the "os" other
    // support field from StrandSupport.
    // Supports for nodes are minimum distinct reads that use the node.
    map<id_t, Support> node_supports;
    // And for edges
    unordered_map<edge_t, Support> edge_supports;
    
    /**
     * Return true if we have support information, and false otherwise.
     */
    virtual bool has_supports() const;
    
    /**
     * Get the Support for a given Node, or 0 if it has no recorded support.
     */
    virtual Support get_support(id_t node);
    
    /**
     * Get the Support for a given Edge, or 0 if it has no recorded support.
     */
    virtual Support get_support(edge_t edge);    
    
    /**
     * Clear the contents.
     */
    virtual void clear();

    /**
    * Read the supports from protobuf.
    */
    void load_supports(istream& in_file);

    /**
     * Read the suppors from output of vg pack
     * Everything put in forward support, average used for nodes
     * Graph must implement VectorizableHandleGraph
     */
    void load_pack_as_supports(const string& pack_file_name, const HandleGraph* vectorizable_graph);

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
 * Flip the orientations of a Support.
 */
Support flip(const Support& to_flip);

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

