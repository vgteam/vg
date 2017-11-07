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
 * Every traversing haplotype is assigned a "lane" number at which it traverses
 * the snarl. Lane is the same looking backward or forward through the snarl.
 * Within each child snarl or child node, the traversal is also assigned a lane.
 * The lane assignments at the start and end nodes are the same, and define the
 * lane assignments for the overall snarl. Traversals can be inserted at any
 * lane number in any internal node.
 */
class SnarlState {

protected:
    
    // This stores, for each overall lane, the traversal in that lane, annotated with its internal lane assignments
    vector<vector<pair<handle_t, size_t>>> haplotypes;
    
    // This stores, for each forward handle, a vector of all the lanes in order.
    // Each lane is holding an iterator to the entry in the haplotyope that
    // occupies that lane. When we insert into or delete out of the vectors in
    // this map, we update the lane numbers where all the iterators point. TODO:
    // really we need to hold skip lists or something; we need efficient insert
    // at index. But since we still need to pay O(N) fixing up stuff after the
    // insert, it might not be worth it.
    unordered_map<handle_t, vector<decltype(haplotypes)::value_type::iterator>> net_node_lanes;
    
    /// We need to keep track of the net graph, because we may need to traverse
    /// haplotypes forward or reverse and we need to flip things.
    const NetGraph* graph;

public:
    
    /// Create a SnarlState that uses the given net graph.
    SnarlState(const NetGraph* graph);
    
    /// Dump internal state to cerr.
    void dump() const;
    
    /// How many haplotypes traverse this snarl?
    size_t size() const;
    
    /// Trace the haplotype int eh given overall lane in the given orientation.
    /// Yields the oriented handles visited and the per-forward-handle lane
    /// assignments.
    void trace(size_t overall_lane, bool backward, const function<void(const handle_t&, size_t)>& iteratee) const;

    // We can add, remove, and swap haplotypes by rank.
    // All these operations are meant to be reversible.
    
    /// Insert the given traversal of this snarl from start to end, with the
    /// given lane assignments for each oriented handle. If handles to the same
    /// node or child snarl appear more than once, their lane numbers must be
    /// strictly increasing.
    void insert(const vector<pair<handle_t, size_t>>& haplotype);
    
    /// Insert the given traversal of this snarl from start to end, assigning
    /// each visit to a handle to the next available lane. Returns the haplotype
    /// annotated with lane assignments. If handles to the same node or child
    /// snarl appear more than once, their lane numbers will be strictly
    /// increasing.
    const vector<pair<handle_t, size_t>>& append(const vector<handle_t>& haplotype);
    
    /// Insert the given traversal of this snarl from start to end, assigning it
    /// to the given overall lane. Returns the haplotype annotated with lane
    /// assignments for all the internal handles. If the internal handles
    /// represent child snarls, this can be used to recurse down and insert
    /// traversals of them at the right lanes. If handles to the same node or
    /// child snarl appear more than once, their assigned lane numbers will be
    /// strictly increasing. Returns the haplotype annotated with lane
    /// assignments.
    const vector<pair<handle_t, size_t>>& insert(size_t overall_lane, const vector<handle_t>& haplotype);
    
    // TODO: can we do an efficient replace? Or should we just drop and add.
    
    /// Erase the traversal of this haplotype in the given overall lane. Shifts
    /// everything in a higher lane 1 rank down. Returns the erased haplotype
    /// and its old lane assignments.
    vector<pair<handle_t, size_t>> erase(size_t overall_lane);
    
    /// Swap the traversals of this haplotype in the two given overall lanes.
    /// Internal lane assignments (as are used by child snarls) are not
    /// affected.
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
