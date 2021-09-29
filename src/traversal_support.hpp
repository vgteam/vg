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
    TraversalSupportFinder(const HandleGraph& graph, SnarlManager& snarl_manager);
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

    /// Average MAPQ of reads that map to a node
    virtual size_t get_avg_node_mapq(id_t node) const = 0;

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
                                                           int ref_trav_idx = -1,
                                                           int* max_trav_size = nullptr);
    
    /// traversals:      get support for each traversal in this set
    /// shared_travs:    if a node appears N times in shared_travs, then it will count as 1 / (N+1) support
    /// shared_support:  optional supports for shared_travs.  used to weight support split by traversal support.
    /// tgt_travs:       if not empty, only compute support for these traversals (remaining slots in output vector left 0)
    /// eclusive_only:   shared_travs are completely ignored
    /// exclusive_count_travs: these traversals get subtracted from supports in the target traversals
    /// exclusive_count_support: used with above, to determine amount of support to subtract
    /// ref_trav_idx:    index of reference traversal if known
    /// max_trav_size:   optional input of max trav size.  useful when longest traversral is outside target set
    virtual vector<Support> get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                      const vector<int>& shared_travs,
                                                      const vector<Support>& shared_support,
                                                      const set<int>& tgt_travs,
                                                      bool exclusive_only,
                                                      const vector<int>& exclusive_count_travs,
                                                      const vector<Support>& exclusive_count_support,
                                                      int ref_trav_idx = -1,
                                                      int* max_trav_size = nullptr) const;
    
    /// Get the total length of all nodes in the traversal
    virtual vector<int> get_traversal_sizes(const vector<SnarlTraversal>& traversals) const;

    /// Get the average MAPQ in each traversal
    /// Only consider nodes
    /// Normalize by base coverage (ie avg coverage / node by node length)
    virtual vector<double> get_traversal_mapqs(const vector<SnarlTraversal>& traversals) const;

    /// Get the average traversal support thresholdek
    virtual size_t get_average_traversal_support_switch_threshold() const;

    /// Relic from old code
    static double support_val(const Support& support) { return total(support); };

    /// get a map of the beginning of a node (in forward orientation) on a traversal
    /// used for up-weighting large deletion edges in complex snarls with average support
    unordered_map<id_t, size_t> get_ref_offsets(const SnarlTraversal& ref_trav) const;

    /// set the threshold
    virtual void set_support_switch_threshold(size_t trav_thresh, size_t node_thresh);

    /// set the breakpoint stricter upper override
    virtual void set_min_bp_edge_override(bool bp_override);

    /// apply the override to a set of traversals
    virtual void apply_min_bp_edge_override(const vector<SnarlTraversal>& traversals,
                                            const set<int>& tgt_travs,
                                            vector<Support>& supports, int ref_trav_idx) const;

protected:

    size_t average_traversal_support_switch_threshold = 50;
    /// Use average instead of minimum support when determining a node's support
    /// its position supports.
    size_t average_node_support_switch_threshold = 50;

    /// If on, always apply minimum edge support for breakpoint (ref->offref) edges
    bool min_bp_edge_override = false;

    const HandleGraph& graph;

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

    /// Average MAPQ of reads that map to a node
    virtual size_t get_avg_node_mapq(id_t node) const;
    
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
    // good if cache_size lines up with FlowCaller::max_snarl_edges in graph_caller.hpp
    CachedPackedTraversalSupportFinder(const Packer& packer, SnarlManager& snarl_manager, size_t cache_size = 500000);
    virtual ~CachedPackedTraversalSupportFinder();

    /// Support of an edge
    virtual Support get_edge_support(id_t from, bool from_reverse, id_t to, bool to_reverse) const;
    
    /// Minimum support of a node
    virtual Support get_min_node_support(id_t node) const;

    /// Average support of a node
    virtual Support get_avg_node_support(id_t node) const;

    /// Average MAPQ of reads that map to a node
    virtual size_t get_avg_node_mapq(id_t node) const;
    
protected:

    /// One node cache per threade
    mutable vector<LRUCache<edge_t, Support>*> edge_support_cache;
    mutable vector<LRUCache<nid_t, Support>*> min_node_support_cache;
    mutable vector<LRUCache<nid_t, Support>*> avg_node_support_cache;
    mutable vector<LRUCache<nid_t, size_t>*> avg_node_mapq_cache;
};

/**
 * Add table to keep track of child snarl support that can be maintained by outside logic
 */
class NestedCachedPackedTraversalSupportFinder : public CachedPackedTraversalSupportFinder {
public:
    NestedCachedPackedTraversalSupportFinder(const Packer& packer, SnarlManager& snarl_manager, size_t cache_size = 500000);
    virtual ~NestedCachedPackedTraversalSupportFinder();
    
    virtual tuple<Support, Support, int> get_child_support(const Snarl& snarl) const;

    /**
     * map used for get_child_support().  It's intialized for every snarl so that the values can be
     * updated from different threads. 
     */
    // todo: why can't we use unordered_map -- there's a hash function in snarls.hpp
    //       perhaps we can switch to pointers but not so sure at moment
    struct snarl_less {
        inline bool operator()(const Snarl& s1, const Snarl& s2) const {
            return s1.start() < s2.start() || (s1.start() == s2.start() && s1.end() < s2.end());
        }
    };
    typedef map<Snarl, tuple<Support, Support, int>, snarl_less> SupportMap;
    SupportMap child_support_map;
};
}

#endif
