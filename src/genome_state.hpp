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

protected:
    /**
     * Store a vector of haplotype traversals of this snarl from start to end.
     * Handles are in the Snarl's NetGraph, and represent visits to either the
     * snarl's own nodes or to child snarls. The handles visiting the snarl's
     * start and end nodes are included.
     */
    vector<vector<handle_t>> haplotypes;

public:
    
    /// How many haplotypes traverse this snarl?
    size_t size() const;
    
    /// We define a way to iterate over the contained haplotypes
    using const_iterator = decltype(haplotypes)::const_iterator;
    
    /// Get an iterator to the first haplotype
    const_iterator begin() const;
    
    /// Get an iterator to the past-the-end haplotype
    const_iterator end() const;
    
    /// Get the haplotype at the given rank.
    const vector<handle_t>& at(size_t rank) const;

    // We can add, remove, and swap haplotypes by rank.
    // All these operations are meant to be reversible.
    
    /// Insert the given traversal of this snarl from start to end as a
    /// haplotype with the given overall rank. The first rank is 0. Pushes
    /// everything with a higher rank 1 rank up.
    void insert(size_t rank, const vector<handle_t>& haplotype);
    
    /// Replace the traversal of this haplotype at the given rank with the given
    /// new traversal from start to end. Returns the replaced haplotype.
    vector<handle_t> replace(size_t rank, const vector<handle_t>& haplotype);
    
    /// Erase the traversal of this haplotype with the given rank. Shifts
    /// everything with a higher rank 1 rank down. Returns the erased haplotype.
    vector<handle_t> erase(size_t rank);
    
    /// Swap the traversals of this haplotype with the two given ranks.
    void swap(size_t rank1, size_t rank2);
    
};

class GenomeState;

/**
 * Represents a modification of a GenomeState.
 * We use a command pattern to enable undo-ability.
 * Applying a command always returns a command that will undo what you did.
 */
struct GenomeStateCommand {
    virtual ~GenomeStateCommand() = default;
    
    /// Execute this command on the given state and return the reverse command.
    /// Generally ends up calling a command-type-specific method on the GenomeState that does the actual work.
    virtual GenomeStateCommand* execute(GenomeState& state) const = 0;
};

struct InsertHaplotypeCommand : public GenomeStateCommand {

    /// In each snarl's state, at each rank, insert the given haplotype.
    /// Includes all the nested snarls.
    /// To be executed in order from left to right (becasue inserting at ranks changes ranks).
    vector<tuple<const Snarl*, size_t, vector<handle_t>>> insertions;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;
    
    virtual ~InsertHaplotypeCommand() = default;
};

struct DeleteHaplotypeCommand : public GenomeStateCommand {
    /// In each (top level) snarl's state, delete the haplotype at the given rank.
    vector<pair<const Snarl*, size_t>> deletions;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;

    virtual ~DeleteHaplotypeCommand() = default;
};

struct SwapHaplotypesCommand : public GenomeStateCommand {
    /// Swap the haplotypes at the given ranks
    pair<size_t, size_t> to_swap;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;
    
    virtual ~SwapHaplotypesCommand() = default;
};

struct CreateHaplotypeCommand : public GenomeStateCommand {
    /// This is the root snarl to create a new traversal of.
    /// TODO: maybe should be a chain?
    const Snarl* root;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;
    
    virtual ~CreateHaplotypeCommand() = default;
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
    
    /// Make a new GenomeState on the given SnarlManager, manageing snarls in
    /// the given graph, with the given telomere pairs
    GenomeState(const SnarlManager& manager, const HandleGraph* graph,
        const unordered_set<pair<const Snarl*, const Snarl*>> telomeres);
    
    // We execute commands and return the inverse commands
    
    /// Create a haplotype and return a command to delete it
    DeleteHaplotypeCommand create_haplotype(const CreateHaplotypeCommand& c);
    
    /// Insert a haplotype and return a command to delete it
    DeleteHaplotypeCommand insert_haplotype(const InsertHaplotypeCommand& c);
    
    /// Delete a haplotype and return a command to insert it
    InsertHaplotypeCommand delete_haplotype(const DeleteHaplotypeCommand& c);
    
    /// Swap two haplotypes and return a command to swap them back
    SwapHaplotypesCommand swap_haplotypes(const SwapHaplotypesCommand& c);
    
    /// Execute a command. Return a new heap-allocated command that undoes the
    /// command being executed. Frees the passed command. TODO: does that make
    /// sense?
    GenomeStateCommand* execute(GenomeStateCommand* command);
     
protected:
    /// We keep track of pairs of telomere unary snarls. The haplotypes we work
    /// on connect a left telomere and its corresponding right telomere. We can
    /// traverse the snarl decomposition from one telomere to the other either
    /// following chains or pollowing paths inside of some snarl we are in.
    unordered_set<pair<const Snarl*, const Snarl*>> telomeres;

    /// We have a state for every snarl.
    unordered_map<const Snarl*, SnarlState> state;
    
    /// We precompute all the net graphs and keep them around.
    unordered_map<const Snarl*, NetGraph> net_graphs;
    
    /// We keep a reference to the SnarlManager that knows what children are
    /// where and which snarl we should look at next after leaving a previous
    /// snarl.
    const SnarlManager& manager;

};

}

#endif
