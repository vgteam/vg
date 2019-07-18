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
    virtual ~SnarlCaller() = 0;

    /// Get the genotype of a site
    virtual vector<int> genotype(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 int ref_trav_idx,
                                 int ploidy) = 0;
};

/**
 * Find the genotype of some traversals in a site using read support
 */ 
class SupportBasedSnarlCaller : public SnarlCaller {
public:
    /// use_avg_node_support: A node's support is the average (instead of minimum)
    /// support across its bases
    /// use_avg_trav_support: A traversals support is the average (instead of minimum)
    /// support across its nodes and edges
    SupportBasedSnarlCaller(const PathHandleGraph& graph, bool use_avg_node_support, bool use_avg_trav_support);
    virtual ~SupportBasedSnarlCaller();

    /// Support of an edge
    virtual Support get_edge_support(const edge_t& edge) const = 0;

    /// Effective length of an edge
    virtual int64_t get_edge_length(const edge_t& edge) const;

    /// Minimum support of a node
    virtual Support get_min_node_support(id_t node) const = 0;

    /// Average support of a node
    virtual Support get_avg_node_support(id_t node) const = 0;

    /// Get the support of a node
    function<Support(id_t)> get_node_support;

    /// Use node or edge support as proxy for child support (as was done in original calling code)
    virtual pair<Support, int> get_child_support(const Snarl& snarl) const;

    /// Get the genotype of a site
    virtual vector<int> genotype(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 int ref_trav_idx,
                                 int ploidy);

    /// Get the support of a traversal
    /// Child snarls are handled as in the old call code: their maximum support is used
    virtual Support get_traversal_support(const SnarlTraversal& traversal) const;

    /// Get the support of a set of traversals.  Any support overlapping traversals in shared_travs
    /// will have their support split.  If exclusive_only is true, then any split support gets
    /// rounded down to 0
    virtual vector<Support> get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                      const vector<int>& shared_travs,
                                                      bool exclusive_only) const;

    /// Get the total length of all nodes in the traversal
    virtual vector<int> get_traversal_sizes(const vector<SnarlTraversal>& traversals) const;

protected:

    /// Get the best support out of a list of supports, ignoring skips
    static int get_best_support(const vector<Support>& supports, const vector<int>& skips);

    /// Relic from old code
    static double support_val(const Support& support) { return total(support); };


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
    double max_indel_het_bias =3;
    /// Like the max het bias, but applies to multiallelic indels.
    double max_indel_ma_bias = 6;

    ///

protected:

    const PathHandleGraph& graph;
    bool use_avg_node_support;
    bool use_avg_trav_support;

    // todo: background support
};

/**
 * Get the read support from a Packer object
 */ 
class PackedSupportSnarlCaller : public SupportBasedSnarlCaller {
public:
    PackedSupportSnarlCaller(const Packer& packer, bool use_avg_node_support, bool use_avg_trav_support);
    virtual ~PackedSupportSnarlCaller();

    /// Support of an edge
    virtual Support get_edge_support(const edge_t& edge) const;

    /// Minimum support of a node
    virtual Support get_min_node_support(id_t node) const;

    /// Average support of a node
    virtual Support get_avg_node_support(id_t node) const;
    
protected:

    /// Derive supports from this pack index
    const Packer& packer;
};

}
#endif
