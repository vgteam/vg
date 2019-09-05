#ifndef VG_SNARL_CALLER_HPP_INCLUDED
#define VG_SNARL_CALLER_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
#include "snarls.hpp"
#include "genotypekit.hpp"
#include "packer.hpp"

namespace vg {

using namespace std;


/**
 * SnarlCaller: Given a list of traversals through a site, come up with a genotype
 * come up with a genotype
 */ 
class SnarlCaller {
public:
    virtual ~SnarlCaller();

    /// Get the genotype of a site
    virtual vector<int> genotype(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 int ref_trav_idx,
                                 int ploidy) = 0;

    /// Update INFO and FORMAT fields of the called variant
    virtual void update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const string& sample_name,
                                 vcflib::Variant& variant) = 0;

    /// Define any header fields needed by the above
    virtual void update_vcf_header(string& header) const = 0;
};

/**
 * Find the genotype of some traversals in a site using read support
 */ 
class SupportBasedSnarlCaller : public SnarlCaller {
public:
   SupportBasedSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager);
    virtual ~SupportBasedSnarlCaller();

    /// Support of an edge
    virtual Support get_edge_support(const edge_t& edge) const = 0;
    virtual Support get_edge_support(id_t from, bool from_reverse, id_t to, bool to_reverse) const = 0;

    /// Effective length of an edge
    virtual int64_t get_edge_length(const edge_t& edge, const unordered_map<id_t, size_t>& ref_offsets) const;

    /// Minimum support of a node
    virtual Support get_min_node_support(id_t node) const = 0;

    /// Average support of a node
    virtual Support get_avg_node_support(id_t node) const = 0;

    /// Use node or edge support as proxy for child support (as was done in original calling code)
    virtual tuple<Support, Support, int> get_child_support(const Snarl& snarl) const;

    /// Get the genotype of a site
    virtual vector<int> genotype(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 int ref_trav_idx,
                                 int ploidy);

    /// Update INFO and FORMAT fields of the called variant
    virtual void update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const string& sample_name,
                                 vcflib::Variant& variant);

    /// Define any header fields needed by the above
    virtual void update_vcf_header(string& header) const;

    /// Get the support of a traversal
    /// Child snarls are handled as in the old call code: their maximum support is used
    virtual Support get_traversal_support(const SnarlTraversal& traversal) const;

    /// Get the support of a set of traversals.  Any support overlapping traversals in shared_travs
    /// will have their support split.  If exclusive_only is true, then any split support gets
    /// rounded down to 0.  if the ref_trav_idx is given, it will be used for computing (deletion) edge lengths
    virtual vector<Support> get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                      const vector<int>& shared_travs,
                                                      bool exclusive_only, int ref_trav_idx = -1) const;

    /// Get the total length of all nodes in the traversal
    virtual vector<int> get_traversal_sizes(const vector<SnarlTraversal>& traversals) const;

protected:

    /// Get the best support out of a list of supports, ignoring skips
    static int get_best_support(const vector<Support>& supports, const vector<int>& skips);

    /// Relic from old code
    static double support_val(const Support& support) { return total(support); };

    /// Get the bias used to for comparing two traversals
    /// (It differrs heuristically depending whether they are alt/ref/het/hom/snp/indel
    ///  see tuning parameters below)
    double get_bias(const vector<int>& traversal_sizes, int best_trav,
                    int second_best_trav, int ref_trav_idx) const;

    /// get a map of the beginning of a node (in forward orientation) on a traversal
    /// used for up-weighting large deletion edges in complex snarls with average support
    unordered_map<id_t, size_t> get_ref_offsets(const SnarlTraversal& ref_trav) const;

    /// Tuning

    /// What's the minimum integer number of reads that must support a call? We
    /// don't necessarily want to call a SNP as het because we have a single
    // supporting read, even if there are only 10 reads on the site.
    int min_total_support_for_call = 1;
    /// What fraction of the reads supporting an alt are we willing to discount?
    /// At 2, if twice the reads support one allele as the other, we'll call
    /// homozygous instead of heterozygous. At infinity, every call will be
    /// heterozygous if even one read supports each allele.
    double max_het_bias = 10;
    /// Like above, but applied to ref / alt ratio (instead of alt / ref)
    double max_ref_het_bias = 4.5;
    /// Like the max het bias, but applies to novel indels.
    double max_indel_het_bias = 6;
    /// As above but used for alt1/alt2 calls
    double max_ma_bias = 3;
    /// what's the minimum ref or alt allele depth to give a PASS in the filter
    /// column? Also used as a min actual support for a second-best allele call
    size_t min_mad_for_filter = 1;
    /// what's the min log likelihood for allele depth assignments to PASS?
    double min_ad_log_likelihood_for_filter = -9;
    /// Use average instead of minimum support when determining a traversal's support
    /// its node and edge supports.
    size_t average_traversal_support_switch_threshold = 50;
    /// Use average instead of minimum support when determining a node's support
    /// its position supports.
    size_t average_node_support_switch_threshold = 50;

    const PathHandleGraph& graph;

    SnarlManager& snarl_manager;

    // todo: background support

    
};

/**
 * Get the read support from a Packer object
 */ 
class PackedSupportSnarlCaller : public SupportBasedSnarlCaller {
public:
    PackedSupportSnarlCaller(const Packer& packer, SnarlManager& snarl_manager);
    virtual ~PackedSupportSnarlCaller();

    /// Support of an edge
    virtual Support get_edge_support(const edge_t& edge) const;
    virtual Support get_edge_support(id_t from, bool from_reverse, id_t to, bool to_reverse) const;

    /// Minimum support of a node
    virtual Support get_min_node_support(id_t node) const;

    /// Average support of a node
    virtual Support get_avg_node_support(id_t node) const;
    
protected:

    /// Derive supports from this pack index
    const Packer& packer;
};

// debug helpers
inline string to_string(const HandleGraph& graph, handle_t handle) {
    return std::to_string(graph.get_id(handle)) + ":" + std::to_string(graph.get_is_reverse(handle));
}
inline string to_string(const HandleGraph& graph, edge_t edge) {
    return to_string(graph, edge.first) + " -> " + to_string(graph, edge.second);
}

}
#endif
