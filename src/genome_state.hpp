#ifndef VG_GENOME_STATE_HPP_INCLUDED
#define VG_GENOME_STATE_HPP_INCLUDED

/**
 * \file genome_state.hpp
 *
 * MCMC-friendly genome state representation.
 */

#include "handle.hpp"
#include "snarls.hpp"
#include "yeet.hpp"

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
    
    /// Insert the given traversal of this snarl from start to end or end to
    /// start (as determined by the backward flag), assigning each visit to a
    /// handle to the next available lane. Returns the haplotype annotated with
    /// lane assignments. If handles to the same node or child snarl appear more
    /// than once, their lane numbers will be strictly increasing.
    const vector<pair<handle_t, size_t>>& append(const vector<handle_t>& haplotype, bool backward = false);
    
    /// Insert the given traversal of this snarl from start to end or end to
    /// start (as determined by the backward flag), assigning it to the given
    /// overall lane. Returns the haplotype annotated with lane assignments for
    /// all the internal handles. If the internal handles represent child
    /// snarls, this can be used to recurse down and insert traversals of them
    /// at the right lanes. If handles to the same node or child snarl appear
    /// more than once, their assigned lane numbers will be strictly increasing.
    /// Returns the haplotype annotated with lane assignments.
    const vector<pair<handle_t, size_t>>& insert(size_t overall_lane, const vector<handle_t>& haplotype, bool backward = false);
    
    // TODO: can we do an efficient replace? Or should we just drop and add.
    
    /// Erase the traversal of this haplotype in the given overall lane. Shifts
    /// everything in a higher lane 1 rank down. Returns the erased haplotype
    /// and its old lane assignments.
    vector<pair<handle_t, size_t>> erase(size_t overall_lane);
    
    /// Swap the traversals of this haplotype in the two given overall lanes.
    /// Internal lane assignments (as are used by child snarls) are not
    /// affected.
    void swap(size_t lane1, size_t lane2);
    
};

class GenomeState;

// We have all these commands. All the handles in the commands are in the net
// graph of the appropriate snarl, as owned by the GenomeState. If you want to
// populate a command, you need to go get the net graph from the genome state to
// make your handles.

/**
 * Represents a modification of a GenomeState.
 * We use a command pattern to enable undo-ability.
 * Applying a command always returns a command that will undo what you did.
 */
struct GenomeStateCommand {
    virtual ~GenomeStateCommand() = default;
    
    /// Execute this command on the given state and return the reverse command.
    /// Generally ends up calling a command-type-specific method on the
    /// GenomeState that does the actual work.
    virtual GenomeStateCommand* execute(GenomeState& state) const = 0;
};

struct InsertHaplotypeCommand : public GenomeStateCommand {

    /// For each snarl, holds several haplotype traversals. The handles in each
    /// traversal are annotated with their local lane assignments. This
    /// annotation lets the command be basically an undelete, because we can
    /// exactly reverse a delete. Insertions are applied from begin to end
    /// within each vector. Lane numbers for a given forward handle must
    /// strictly increase between vectors and within vectors.
    unordered_map<const Snarl*, vector<vector<pair<handle_t, size_t>>>> insertions;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;
    
    virtual ~InsertHaplotypeCommand() = default;
};

struct DeleteHaplotypeCommand : public GenomeStateCommand {
    
    /// For each snarl, delete the haplotype in each overall lane in the vector.
    /// You must specify out the deletions for all the snarls in a haplotype; we
    /// won't automatically go and find the children if you just list to delete
    /// from the parents. Deletions happen from begin to end through each
    /// vector. Lane numbers for a given snarl must be strictly decreasing.
    unordered_map<const Snarl*, vector<size_t>> deletions;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;

    virtual ~DeleteHaplotypeCommand() = default;
};

struct SwapHaplotypesCommand : public GenomeStateCommand {
    /// Work on the given pair of telomeres
    pair<const Snarl*, const Snarl*> telomere_pair;
    /// Swap the haplotypes at the given ranks
    pair<size_t, size_t> to_swap;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;
    
    virtual ~SwapHaplotypesCommand() = default;
};

struct AppendHaplotypeCommand : public GenomeStateCommand {
    /// We just feed in a full traversal from one end of a telomere pair to the
    /// other. Must start and end on the boundary nodes of telomere snarls.
    /// TODO: unary telomeres will work strangely. Internally the GenomeState
    /// has to work out how to divide this into snarls.
    vector<handle_t> haplotype;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;
    
    virtual ~AppendHaplotypeCommand() = default;
};

struct ReplaceSnarlHaplotypeCommand : public GenomeStateCommand {
    /// Which snarl are we working on?
    const Snarl* snarl;
    /// Which lane in the snarl are we changing?
    size_t lane;
    /// What fully specified haplotype should we replace it with? This gets
    /// around any problems to do with increases or decreases in the copy
    /// numbers of child snarls being underspecified.
    vector<handle_t> haplotype;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;
    
    virtual ~ReplaceSnarlHaplotypeCommand() = default;
    
};

/// We can use this to replace a local haplotype within one or more snarls.
/// We also could just express all the deletions and insertions in terms of
/// this.
struct ReplaceLocalHaplotypeCommand : public GenomeStateCommand {
    /// Holds, for each snarl, the overall lanes that have to be deleted. Each
    /// deleted overall lane for a snarl that is still traversed from its parent
    /// needs to be replaced in insertions. Deletions happen first.
    unordered_map<const Snarl*, vector<size_t>> deletions;
    
    /// Holds, for each Snarl, the visits and their lane assignments that have
    /// to be inserted.
    unordered_map<const Snarl*, vector<vector<pair<handle_t, size_t>>>> insertions;
    
    virtual GenomeStateCommand* execute(GenomeState& state) const;
    
    virtual ~ReplaceLocalHaplotypeCommand() = default;
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
    
    /// Dump internal state to cerr.
    void dump() const;
    
    /// Make a new GenomeState on the given SnarlManager, manageing snarls in
    /// the given graph, with the given telomere pairs. Each telomere can appear
    /// in only one pair.
    GenomeState(const SnarlManager& manager, const HandleGraph* graph,
        const unordered_set<pair<const Snarl*, const Snarl*>> telomeres);
        
    // We can give you a net graph for a snarl so that you can make handles to
    // populate commands.
    const NetGraph* get_net_graph(const Snarl* snarl);
    
    // We execute commands and return the inverse commands
    
    // We have high level commands which are easy to sample
    
    /// Create a haplotype and return a command to delete it
    DeleteHaplotypeCommand append_haplotype(const AppendHaplotypeCommand& c);
    
    /// Swap two haplotypes and return a command to swap them back
    SwapHaplotypesCommand swap_haplotypes(const SwapHaplotypesCommand& c);
    
    /// Replace part(s) of some haplotype(s) with other material.
    ReplaceLocalHaplotypeCommand replace_snarl_haplotype(const ReplaceSnarlHaplotypeCommand& c);
    
    // And low level commands which are very specific and good for being undos
    
    /// Insert a haplotype and return a command to delete it
    DeleteHaplotypeCommand insert_haplotype(const InsertHaplotypeCommand& c);
    
    /// Delete a haplotype and return a command to insert it
    InsertHaplotypeCommand delete_haplotype(const DeleteHaplotypeCommand& c);
    
    /// Replace part(s) of some haplotype(s) with other material.
    ReplaceLocalHaplotypeCommand replace_local_haplotype(const ReplaceLocalHaplotypeCommand& c);
    
    /// Execute a command. Return a new heap-allocated command that undoes the
    /// command being executed. Frees the passed command. TODO: does that make
    /// sense?
    GenomeStateCommand* execute(GenomeStateCommand* command);
    
    /// Count the number of haplotypes connecting the given pair of telomeres.
    /// It must be a pair of telomeres actually passed during construction.
    size_t count_haplotypes(const pair<const Snarl*, const Snarl*>& telomere_pair) const;
    
    /// Count the number of haplotypes within a given snarl
    size_t count_haplotypes(const Snarl* snarl) const;
    
    /// Trace the haplotype starting at the given start telomere snarl with the
    /// given overall lane. Calls the callback with each backing HandleGraph
    /// handle.
    void trace_haplotype(const pair<const Snarl*, const Snarl*>& telomere_pair,
        size_t overall_lane, const function<void(const handle_t&)>& iteratee) const;
     
protected:
    /// We keep track of pairs of telomere snarls. The haplotypes we work on
    /// connect a left telomere and its corresponding right telomere. We can
    /// traverse the snarl decomposition from one telomere to the other either
    /// following chains or pollowing paths inside of some snarl we are in.
    /// The start snarl of a pair must be in its local forward orientation.
    unordered_set<pair<const Snarl*, const Snarl*>> telomeres;

    /// We precompute all the net graphs and keep them around.
    unordered_map<const Snarl*, NetGraph> net_graphs;

    /// We have a state for every snarl, which depends on the corresponding net graph.
    unordered_map<const Snarl*, SnarlState> state;
    
    /// We remember the backing graph, because sometimes we need to translate
    /// from backing graph handles to IDs and orientations to look up snarls.
    const HandleGraph* backing_graph;
    
    /// We keep a reference to the SnarlManager that knows what children are
    /// where and which snarl we should look at next after leaving a previous
    /// snarl.
    const SnarlManager& manager;
    
    /// We have a generic stack-based handle-vector-to-per-snarl-haplotypes
    /// insertion walker function. The handles to add have to span one or more
    /// entire snarls, and we specify the lane in the spanned snarls to put them
    /// in.
    void insert_handles(const vector<handle_t>& to_add,
        unordered_map<const Snarl*, vector<size_t>>& lanes_added,
        size_t top_lane = numeric_limits<size_t>::max());

};

}

#endif
