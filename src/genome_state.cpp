#include "genome_state.hpp"

namespace vg {

using namespace std;

SnarlState::SnarlState(const NetGraph* graph) : graph(graph) {
    // Nothing to do!
}

size_t SnarlState::size() const {
    return haplotypes.size();
}

void SnarlState::trace(size_t rank, bool backward, const function<void(const handle_t&, size_t)>& iteratee) const {
    // Get the haplotype we want to loop over
    auto& haplotype = haplotypes.at(rank);
    
    auto process_traversal = [&](const pair<handle_t, size_t>& handle_and_rank) {
        // For every handle in the haplotype, yield it either forward or
        // backward as determined by our traversal direction.
        iteratee(backward ? graph->flip(handle_and_rank.first) : handle_and_rank.first, handle_and_rank.second);
    };
    
    if (backward) {
        // If we're going backward, go in reverse order.
        // See <https://stackoverflow.com/a/23094303>
        for_each(haplotype.rbegin(), haplotype.rend(), process_traversal);
    } else {
        // Otherwise go in forward order
        for_each(haplotype.begin(), haplotype.end(), process_traversal);
    }
}

void SnarlState::insert(size_t rank, const vector<handle_t>& haplotype) {
    // We're going to insert the given haplotype at the given initial rank.
    
    // For each handle in the haplotype
    
    // If it's the first or the last
    
    // Insert it at the specified rank and push everything else up
    
    // Otherwise just append it at the next available rank.
    
    // Insert it at the last rank for all visits to the forward version of that handle.
    
    // For the first visit (and the last if it isn't equal to the first), actually make it be at the specified rank and move all the other things at or after that rank up.
}


GenomeState::GenomeState(const SnarlManager& manager, const HandleGraph* graph,
    const unordered_set<pair<const Snarl*, const Snarl*>> telomeres) : telomeres(telomeres), manager(manager) {

    manager.for_each_snarl_preorder([&](const Snarl* snarl) {
        // For each snarl
        
        // Make a net graph for it. TODO: we're not considering internal
        // connectivity, but what we really should do is consider internal
        // connectivity but only allowing for start to end traversals (but
        // including in unary snarls)
        net_graphs.emplace(snarl, manager.net_graph_of(snarl, graph, false));
        
        // Make an empty state for it using the net graph
        state.emplace(snarl, SnarlState(&net_graphs.at(snarl)));
        
        // TODO: can we just make the net graph live in the state?
        
        
    });
}

DeleteHaplotypeCommand GenomeState::create_haplotype(const CreateHaplotypeCommand& c) {
    // Sample a new haplotype recursively
    
    // TODO: what RNG are we going to use?
}

DeleteHaplotypeCommand GenomeState::insert_haplotype(const InsertHaplotypeCommand& c) {
    DeleteHaplotypeCommand to_return;
    
    for (const tuple<const Snarl*, size_t, vector<handle_t>>& action : c.insertions) {
        // For every insert action we have to do, in order, do it
        // Insert into this snarl at this rank this haplotype
        state.at(get<0>(action)).insert(get<1>(action), get<2>(action));
        
        // Record how to undo it (delete from this snarl at this rank).
        to_return.deletions.emplace_back(get<0>(action), get<1>(action));
    }
    
    return to_return;
}


GenomeStateCommand* GenomeState::execute(GenomeStateCommand* command) {
    // Just make the command tell us what type it is
    return command->execute(*this);
}

size_t GenomeState::count_haplotypes(const pair<const Snarl*, const Snarl*>& telomere_pair) {
    // We assume all the traversals go through the whole chromosome from telomere to telomere.
    return state.at(telomere_pair.first).size();
}

void GenomeState::trace_haplotype(const pair<const Snarl*, const Snarl*>& telomere_pair,
    size_t rank, const function<void(const handle_t&)> iteratee) {
    
    // Start at the first snarl in the telomere pair
    
    // Start at the rank we were given
    
    // Trace out its traversal
    
    // If we visit a real node, yield it.
    
    // If we visit a child snarl, we need to recurse into it with the appropriate rank.
    
}


}
