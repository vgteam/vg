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

class GenomeState;

/**
 * Represents a modification of a GenomeState.
 * We use a command pattern to enable undo-ability.
 * Applying a command always returns a command that will undo what you did.
 */
struct GenomeStateCommand {
    virtual ~GenomeStateCommand() = default;
    
    /// Execute this command on the given state and return the reverse command.
    virtual GenomeStateCommand* execute(GenomeState& state) const;
};

struct InsertHaplotypeCommand : public GenomeStateCommand {

    /// In each snarl's state, at each rank, insert the given haplotype.
    /// Includes all the nested snarls.
    /// To be executed in order from left to right (becasue inserting at ranks changes ranks).
    vector<tuple<const Snarl*, size_t, vector<handle_t>>> insertions;
    
    virtual ~InsertHaplotypeCommand() = default;
};

struct DeleteHaplotypeCommand : public GenomeStateCommand {
    /// In each (top level) snarl's state, delete the haplotype at the given rank.
    
    vector<pair<const Snarl*, size_t>> deletions;

    virtual ~DeleteHaplotypeCommand() = default;
};

struct SwapHaplotypesCommand : public GenomeStateCommand {
    pair<size_t, size_t> to_swap;
    
    virtual ~SwapHaplotypesCommand() = default;
};

struct CreateHaplotypeCommand : public GenomeStateCommand {
    /// This is the root snarl to create a new traversal of.
    /// TODO: maybe should be a chain?
    const Snarl* root;
    
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
    
    // We execute commands and return the inverse commands
    
    DeleteHaplotypeCommand create_haplotype(const CreateHaplotypeCommand& c);
    
    DeleteHaplotypeCommand insert_haplotype(const InsertHaplotypeCommand& c);
    
    InsertHaplotypeCommand delete_haplotype(const DeleteHaplotypeCommand& c);
    
    SwapHaplotypesCommand swap_haplotypes(const SwapHaplotypesCommand& c);
    
    /// Execute a command.
    /// Return a new heap-allocated command that undoes the command being executed.
    /// Frees the passed command. TODO: does that make sense?
    GenomeStateCommand* execute(GenomeStateCommand* command);
     
protected:
    /// We have a state for every snarl. TODO: What about how top-level snarls
    /// form chains? Mightn't you need to traverse several? We don't really have
    /// a root-level state.
    unordered_map<const Snarl*, SnarlState> state;
    
    /// We precompute all the net graphs and keep them around.
    unordered_map<const Snarl*, NetGraph> net_graphs;
    
    /// We keep a reference to the SnarlManager that knows what children are
    /// where.
    const SnarlManager& manager;

};

}

#endif
