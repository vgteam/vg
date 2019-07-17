#ifndef VG_TRAVERSAL_GENOTYPER_HPP_INCLUDED
#define VG_TRAVERSAL_GENOTYPER_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
#include "snarls.hpp"

namespace vg {

using namespace std;


/**
 * TraversalGenotyper: Given a list of traversals through a site, come up with a genotype
 * come up with a genotype
 */ 
class TraversalGenotyper {
public:
    virtual ~TraversalGenotyper() = 0;

    /// Get the genotype of a site
    virtual vector<int> genotype(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 int ref_trav_idx,
                                 int ploidy) = 0;
};

/**
 * Find the genotype of some traversals in a site using read support
 */ 
class SupportBasedTraversalGenotyper : public TraversalGenotyper {
public:
    /// use_avg_node_support: A node's support is the average (instead of minimum)
    /// support across its bases
    /// use_avg_trav_support: A traversals support is the average (instead of minimum)
    /// support across its nodes and edges
    SupportBasedTraversalGenotyper(const PathHandleGraph& graph, bool use_avg_node_support, bool use_avg_trav_support);
    virtual ~SupportBasedTraversalGenotyper();

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
    virtual Support get_traversal_support(const SnarlTraversal& traversal);

    /// Get the support of a set of traversals (desribed by trav_indexes), assuming even split between overlaps
    virtual vector<Support> get_traversal_set_support(const vector<SnarlTraversal>& traversals,
                                                      const vector<int>& trav_indxes,
                                                      bool exclusive_only);

protected:

    const PathHandleGraph& graph;
    bool use_avg_node_support;
    bool use_avg_trav_support;

    // todo: background support


};

}
#endif
