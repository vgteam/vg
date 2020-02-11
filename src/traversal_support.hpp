#ifndef VG_TRAVERSAL_SUPPORT_HPP_INCLUDED
#define VG_TRAVERSAL_SUPPORT_HPP_INCLUDED

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
 * Get the read support of snarl traversals or sets of snarl traversals
 */ 
class TraversalSupportFinder {
public:
    TraversalSupportFinder(const PathHandleGraph& graph, SnarlManager& snarl_manager);
    virtual ~TraversalSupportFinder();

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

    /// Get the support of a traversal
    /// Child snarls are handled as in the old call code: their maximum support is used
    virtual Support get_traversal_support(const SnarlTraversal& traversal) const;
 
    /// wrapper for using get_traversal_set_support to get the support for
    /// some alleles in a genotype, where everything is split evently among them
    /// anything not in the genotype gets a support using "exclusive_count"
    /// where nodes taken by the genotype are counted as 0
    /// stuff not in the genotype is limited to other_trav_subset (or all if empty)
    virtual vector<Support> get_traversal_genotype_support(const vector<SnarlTraversal>& traversals,
                                                           const vector<int>& genotype,
                                                           const set<int>& other_trav_subset,
                                                           int ref_trav_idx = -1);
    
    /// traversals:      get support for each traversal in this set
    /// shared_travs:    if a node appears N times in shared_travs, then it will count as 1 / (N+1) support
    /// tgt_travs:       if not empty, only compute support for these traversals (remaining slots in output vector left 0)
    /// eclusive_only:   shared_travs are completely ignored
    /// exclusive_count_support: anything in shared_travs has this much support subtracted from it
    /// mutual_shared:   shared_travs count as 1/N support (instead of 1/(N+1)).  usefuly for total support
    /// ref_trav_idx:    index of reference traversal if known
    virtual vector<Support> get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                      const vector<int>& shared_travs,
                                                      const set<int>& tgt_travs,
                                                      bool exclusive_only,
                                                      const vector<Support>& exclusive_count_support,
                                                      bool mutual_shared,
                                                      int ref_trav_idx = -1) const;

    /// Get the total length of all nodes in the traversal
    virtual vector<int> get_traversal_sizes(const vector<SnarlTraversal>& traversals) const;

    /// Get the average traversal support thresholdek
    virtual size_t get_average_traversal_support_switch_threshold() const;

    /// Relic from old code
    static double support_val(const Support& support) { return total(support); };

    /// get a map of the beginning of a node (in forward orientation) on a traversal
    /// used for up-weighting large deletion edges in complex snarls with average support
    unordered_map<id_t, size_t> get_ref_offsets(const SnarlTraversal& ref_trav) const;

    /// set the threshold
    virtual void set_support_switch_threshold(size_t trav_thresh, size_t node_thresh);

protected:

    size_t average_traversal_support_switch_threshold = 50;
    /// Use average instead of minimum support when determining a node's support
    /// its position supports.
    size_t average_node_support_switch_threshold = 50;

    const PathHandleGraph& graph;

    SnarlManager& snarl_manager;

};

/**
 * Get the read support from a Packer object
 */ 
class PackedTraversalSupportFinder : public TraversalSupportFinder {
public:
    PackedTraversalSupportFinder(const Packer& packer, SnarlManager& snarl_manager);
    virtual ~PackedTraversalSupportFinder();

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

/**
 * Add a caching overlay to the PackedTravesalSupportFinder to avoid frequent
 * base queries which can become expensive.  Even caching the edges seems
 * to have an impact
 */
class CachedPackedTraversalSupportFinder : public PackedTraversalSupportFinder {
public:
    CachedPackedTraversalSupportFinder(const Packer& packer, SnarlManager& snarl_manager, size_t cache_size = 100000);
    virtual ~CachedPackedTraversalSupportFinder();

    /// Support of an edge
    virtual Support get_edge_support(id_t from, bool from_reverse, id_t to, bool to_reverse) const;
    
    /// Minimum support of a node
    virtual Support get_min_node_support(id_t node) const;

    /// Average support of a node
    virtual Support get_avg_node_support(id_t node) const;
    
protected:

    /// One node cache per threade
    mutable vector<LRUCache<edge_t, Support>*> edge_support_cache;
    mutable vector<LRUCache<nid_t, Support>*> min_node_support_cache;
    mutable vector<LRUCache<nid_t, Support>*> avg_node_support_cache;
};


}

#endif
