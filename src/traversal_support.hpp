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

    /// Get the support of a set of traversals.  Any support overlapping traversals in shared_travs
    /// will have their support split.  If exclusive_only is true, then any split support gets
    /// rounded down to 0 (and ignored when computing mins or averages) .
    /// exclusive_count is like exclusive only except shared traversals will be counted (as 0)
    /// when doing average and min support
    /// if the ref_trav_idx is given, it will be used for computing (deletion) edge lengths
    /// if unique is true, then every node or edge will only be counted once
    /// (useful for total support)
    virtual vector<Support> get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                      const vector<int>& shared_travs,
                                                      bool exclusive_only,
                                                      bool exclusive_count,
                                                      bool unique,
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

}

#endif
