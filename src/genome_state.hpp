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
 *
 * Multiple haplotypes are distinguished by their rank. The rank at which each
 * haplotype visits a child snarl on each visit is recorded.
 */
class SnarlState {

protected:
    // TODO: The rank updates here are going to be O(n^2) in total traversal
    // count, because we're shifting vectors around and updating a bunch of ints
    // on every edit.
    
    // Really we need a dynamic rank/select sort of thing. But for traversal
    // counts on the order of 2, like we expect, this might end up faster.
    
    // We assign ranks to each traversal of each node. The ranks on the end
    // nodes are also the rank of the traversal overall, since each traversal
    // visits each end node only once. We constrain the ranks on the end nodes
    // to always be the same. Internally, rank is more or less arbitrary, but we
    // use it, when we visit child snarls, to say which traversal of that child
    // snarl is being visited.

    /**
     * Store a vector of haplotype traversals of this snarl from start to end.
     * Each traversal is paired with its rank among traversals of that handle
     * (in either orientation) in this SnarlState. Handles are in the Snarl's
     * NetGraph, and represent visits to either the snarl's own nodes or to
     * child snarls. The handles visiting the snarl's start and end nodes are
     * included.
     */
    vector<vector<pair<handle_t, size_t>>> haplotypes;

    /**
     * For each locally forward handle, store a vector of all of the visits to
     * it in either orientation, in order by rank.
     */
    unordered_map<handle_t, vector<decltype(haplotypes)::iterator>> visits_by_rank;
    
    /// We need to keep track of the net graph, because we may need to traverse
    /// haplotypes forward or reverse and we need to flip things.
    const NetGraph* graph;

public:
    
    /// Create a SnarlState that uses the given net graph.
    SnarlState(const NetGraph* graph);
    
    /// How many haplotypes traverse this snarl?
    size_t size() const;
    
    /// Trace the haplotype with the given rank. Call the iteratee on each
    /// handle in the net graph with the rank of this traversal among traversals
    /// of that handle's node or child snarl.
    /// If backward is true, traverses in reverse.
    void trace(size_t rank, bool backward, const function<void(const handle_t&, size_t)>& iteratee) const;

    // We can add, remove, and swap haplotypes by rank.
    // All these operations are meant to be reversible.
    
    /// Insert the given traversal of this snarl from start to end as a
    /// haplotype with the given overall rank. The first rank is 0. Pushes
    /// everything with a higher rank 1 rank up.
    /// 
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
    // No data here; a telomere pair should be randomly chosen.    
    
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
    /// the given graph, with the given telomere pairs. Each telomere can appear
    /// in only one pair.
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
    
    /// Count the number of haplotypes connecting the given pair of telomeres.
    /// It must be a pair of telomeres actually passed during construction.
    size_t count_haplotypes(const pair<const Snarl*, const Snarl*>& telomere_pair);
    
    /// Trace the haplotype starting at the given start telomere snarl with the
    /// given rank. Calls the callback with each backing HandleGraph handle.
    void trace_haplotype(const pair<const Snarl*, const Snarl*>& telomere_pair,
        size_t rank, const function<void(const handle_t&)> iteratee);
     
protected:
    /// We keep track of pairs of telomere unary snarls. The haplotypes we work
    /// on connect a left telomere and its corresponding right telomere. We can
    /// traverse the snarl decomposition from one telomere to the other either
    /// following chains or pollowing paths inside of some snarl we are in.
    unordered_set<pair<const Snarl*, const Snarl*>> telomeres;

    /// We precompute all the net graphs and keep them around.
    unordered_map<const Snarl*, NetGraph> net_graphs;

    /// We have a state for every snarl, which depends on the corresponding net graph.
    unordered_map<const Snarl*, SnarlState> state;
    
    /// We keep a reference to the SnarlManager that knows what children are
    /// where and which snarl we should look at next after leaving a previous
    /// snarl.
    const SnarlManager& manager;

};

}

#endif
