#ifndef VG_GENOME_STATE_HPP_INCLUDED
#define VG_GENOME_STATE_HPP_INCLUDED

/**
 * \file genome_state.hpp
 *
 * MCMC-friendly genome state representation.
 */

#include "handle.hpp"
#include "snarls.hpp"

#include <vector>
#include <unordered_map>
#include <utility>
#include <tuple>



namespace vg {

using namespace std;

/**
 * Represents the state of a snarl: zero or more haplotypes traversing its
 * NetGraph.
 *
 * Only admits full-length traversals of a snarl from start to end.
 */
class SnarlState {

public:
    
    /**
     * How many haplotypes traverse this snarl?
     */
    size_t size();

    // We can add, remove, and swap haplotypes by rank.
    // All these operations are meant to be reversible.
    // TODO: Write real operation descriptions
    void insert(size_t rank, const vector<handle_t>& haplotype);
    vector<handle_t> replace(size_t rank, const vector<handle_t>& haplotype);
    vector<handle_t> erase(size_t rank);
    void swap(size_t rank1, size_t rank2);

protected:
    /**
     * Store a vector of haplotype traversals of this snarl from start to end.
     * Handles are in the Snarl's NetGraph, and represent visits to either the
     * snarl's own nodes or to child snarls. The handles visiting the snarl's
     * start and end nodes are included.
     */
    vector<vector<handle_t>> haplotypes;
    
    
};
 
/**
 * Define a way to represent a phased set of haplotypes on a graph that is under
 * consideration as a variant calling solution.
 *
 * There should only be one of these for any genotyping run; operations all
 * mutate the single copy and, if you don't like the result, can be rolled back
 * by doing the reverse mutation.
 */
class GenomeState {

public:
    
    /**
     * Sample and add a new traversal of the given snarl. Recursively samples
     * and adds new traversals of all the child snarls that need to be
     * traversed. Returns the rank of the new traversal in the given snarl.
     */
    size_t sample_new(const Snarl* snarl);
    
    /**
     * Erase the haplotype of the given rank form the given snarl. Also erases
     * all the corresponding haplotypes in child snarls, recursively. Returns
     * all the insert operations needed to undo the erase.
     */
     vector<tuple<const Snarl*, size_t, vector<handle_t>>> erase(const Snarl* snarl, size_t rank);
     
protected:
    /// We have a state for every snarl
    unordered_map<const Snarl*, SnarlState> state;
    
    /// We precompute all the net graphs and keep them around.
    unordered_map<const Snarl*, NetGraph> net_graphs;
    
    /// We keep a reference to the SnarlManager that knows what children are
    /// where.
    const SnarlManager& manager;

};

}

#endif
