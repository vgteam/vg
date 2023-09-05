#include "genome_state.hpp"
#include "utility.hpp"

#include <signal.h>

namespace vg {

using namespace std;

SnarlState::SnarlState(const NetGraph* graph) : graph(graph) {
    // Nothing to do!
}

size_t SnarlState::size() const {
    return haplotypes.size();

}
void SnarlState::dump() const {
    // First dump the haplotypes
    for (size_t i = 0; i < haplotypes.size(); i++) {
        cerr << "Haplotype " << i << ":";
        
        for (auto& record : haplotypes.at(i)) {
            cerr << " " << graph->get_id(record.first) << " " << graph->get_is_reverse(record.first)
                << " at lane " << record.second << ",";
        }
        
        cerr << endl;
    }
    
    
    // Then the lanes index
    for (auto& kv : net_node_lanes) {
        cerr << "Net node " << graph->get_id(kv.first) << " " << graph->get_is_reverse(kv.first) << " lanes:" << endl;
        
        for (size_t i = 0; i < kv.second.size(); i++) {
            cerr << "\tLane " << i << ": " << graph->get_id(kv.second.at(i)->first)
                << " " << graph->get_is_reverse(kv.second.at(i)->first)
                << " at lane " << kv.second.at(i)->second << endl;
        }
        
    }
    
}

void SnarlState::trace(size_t overall_lane, bool backward, const function<void(const handle_t&, size_t)>& iteratee) const {
    // Get the haplotype we want to loop over
    auto& haplotype = haplotypes.at(overall_lane);
    
    auto process_traversal = [&](const pair<handle_t, size_t>& handle_and_lane) {
        // For every handle in the haplotype, yield it either forward or
        // backward as determined by our traversal direction.
        iteratee(backward ? graph->flip(handle_and_lane.first) : handle_and_lane.first, handle_and_lane.second);
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

void SnarlState::insert(const vector<pair<handle_t, size_t>>& haplotype) {
    if (haplotype.empty()) {
        throw runtime_error("Tried to add an empty haplotype to a snarl");
    }
    
    if (haplotype.front().first != graph->get_start() || haplotype.back().first != graph->get_end()) {
        // Fail if it's not actually from start to end.
        stringstream message;
        message << "Tried to add a haplotype to a snarl ("
            << graph->get_id(graph->get_start()) << " " << graph->get_is_reverse(graph->get_start())
            << " -> " << graph->get_id(graph->get_end()) << " " << graph->get_is_reverse(graph->get_end())
            << ") that starts at "
            << graph->get_id(haplotype.front().first) << " " <<  graph->get_is_reverse(haplotype.front().first)
            << " and ends at "
            << graph->get_id(haplotype.back().first) << " " <<  graph->get_is_reverse(haplotype.back().first)
            << " and is not a start-to-end traversal of the snarl";
        throw runtime_error(message.str());
    }
    
    if (haplotype.front().second != haplotype.back().second) {
        // Fail if we try to put something at two different overall lanes
        throw runtime_error("Tried to insert a haplotype with different lanes at the snarl start and end nodes.");
    }

    // TODO: all these inserts at indexes are O(N).

    // Insert the whole traversal into haplotypes at the appropriate index for the overall lane
    size_t overall_lane = haplotype.front().second;
    assert(overall_lane == haplotype.back().second);
    auto inserted = haplotypes.emplace(haplotypes.begin() + overall_lane, haplotype);
    
    for (auto it = inserted->begin(); it != inserted->end(); ++it) {
        // For each handle visit
        auto& handle_visit = *it;
        
        // Insert the iterator record at the right place in net_node_lanes
        auto& node_lanes = net_node_lanes[graph->forward(handle_visit.first)];
        auto lane_iterator = node_lanes.emplace(node_lanes.begin() + handle_visit.second, it);
    
        // Look at whatever is after the lane we just inserted
        ++lane_iterator;
        while (lane_iterator != node_lanes.end()) {
            // Update all the subsequent records in that net node's lane list and bump up their internal lane assignments
            
            // First dereference to get the iterator that points to the actual
            // record, then dereference that, and increment the lane number.
            (*(*lane_iterator)).second++;
            
            ++lane_iterator;
        }
    }
}

const vector<pair<handle_t, size_t>>& SnarlState::append(const vector<handle_t>& haplotype, bool backward) {
    assert(!haplotype.empty());
    
    if (backward) {
        if (haplotype.back() != graph->flip(graph->get_start()) || haplotype.front() != graph->flip(graph->get_end())) {
            // Fail if it's not actually from end to start.
            throw runtime_error("Tried to add a haplotype backward to a snarl that is not a end-to-start traversal of that snarl.");
        }
    } else {
        if (haplotype.front() != graph->get_start() || haplotype.back() != graph->get_end()) {
            // Fail if it's not actually from start to end.
            throw runtime_error("Tried to add a haplotype forward to a snarl that is not a start-to-end traversal of that snarl.");
        }
    }
    
    // Make a new haplotype at the end of our haplotypes vector that's big enough.
    haplotypes.emplace_back(haplotype.size());
    auto& inserted = haplotypes.back();
    
    for (size_t i = backward ? (haplotype.size() - 1) : 0;
        backward ? (i != (size_t) -1) : i < haplotype.size();
        backward ? i-- : i++) {
        
        // For each position in the haplotype, going in the correct direction
        
        // Get the handle we want to insert
        auto handle = haplotype.at(i);
        if (backward) {
            // If we're going backward we need to insert everything the other
            // way around
            handle = graph->flip(handle);
        }
        
        // Work out where we are putting it. We should insert left to right when
        // going forward, and right to left when going backward.
        auto inserted_iterator = inserted.begin() + (backward ? haplotype.size() - 1 - i : i);
        
        // Save the handle
        inserted_iterator->first = handle;
        
        // Find the appropriate node lanes collection
        auto& node_lanes = net_node_lanes[graph->forward(handle)];
        // Save the local lane assignment
        inserted_iterator->second = node_lanes.size();
        // And do the insert
        node_lanes.emplace_back(inserted_iterator);
        
#ifdef debug
        cerr << "At haplotype position " << i << "/" << haplotype.size()
            << " inserted " << graph->get_id(handle) << " " << graph->get_is_reverse(handle)
            << " at lane " << inserted_iterator->second << "/" << node_lanes.size() << endl;
#endif
    }
    
    // Return the completed vector with the lane annotations.
    return inserted;
}

const vector<pair<handle_t, size_t>>& SnarlState::insert(size_t overall_lane, const vector<handle_t>& haplotype, bool backward) {
    assert(!haplotype.empty());
    
    if (backward) {
        if (haplotype.back() != graph->flip(graph->get_start()) || haplotype.front() != graph->flip(graph->get_end())) {
            // Fail if it's not actually from end to start.
            throw runtime_error("Tried to add a haplotype backward to a snarl that is not a end-to-start traversal of that snarl.");
        }
    } else {
        if (haplotype.front() != graph->get_start() || haplotype.back() != graph->get_end()) {
            // Fail if it's not actually from start to end.
            throw runtime_error("Tried to add a haplotype forward to a snarl that is not a start-to-end traversal of that snarl.");
        }
    }
    
    // Insert a haplotype record at the specified overall lane that's big enough.
    auto& inserted = *haplotypes.emplace(haplotypes.begin() + overall_lane, haplotype.size());
    
    for (size_t i = backward ? (haplotype.size() - 1) : 0;
        backward ? (i != (size_t) -1) : i < haplotype.size();
        backward ? i-- : i++) {
        
        // For each position in the haplotype, going in the correct direction
        
        // Get the handle we want to insert
        auto handle = haplotype.at(i);
        if (backward) {
            // If we're going backward we need to insert everything the other
            // way around
            handle = graph->flip(handle);
        }
        
        // Work out where we are putting it. We should insert left to right when
        // going forward, and right to left when going backward.
        auto inserted_iterator = inserted.begin() + + (backward ? haplotype.size() - 1 - i : i);
        
        // Save the handle
        inserted_iterator->first = handle;
        
        // Find the appropriate node lanes collection
        auto& node_lanes = net_node_lanes[graph->forward(handle)];
        
        if (inserted_iterator == inserted.begin() || inserted_iterator + 1 == inserted.end()) {
            // Start and end visits get placed at the predetermined overall_lane
            inserted_iterator->second = overall_lane;
            
            // Insert at the correct offset
            auto lane_iterator = node_lanes.emplace(node_lanes.begin() + overall_lane, inserted_iterator);
            
            // Look at whatever is after the lane we just inserted
            ++lane_iterator;
            while (lane_iterator != node_lanes.end()) {
                // Update all the subsequent records in that net node's lane list and bump up their internal lane assignments
                
                // First dereference to get the iterator that points to the actual
                // record, then dereference that, and increment the lane number.
                (*(*lane_iterator)).second++;
                
                ++lane_iterator;
            }
                
        } else {
            // Interior visits just get appended, which is simplest. No need to bump anything up.
            
            // Save the local lane assignment
            inserted_iterator->second = node_lanes.size();
            // And do the insert
            node_lanes.emplace_back(inserted_iterator);
        }
        
        // Insert the next handle in the next slot in the haplotype
        ++inserted_iterator;
    } 
    
    // Return the annotated haplotype.
    return inserted;
}

vector<pair<handle_t, size_t>> SnarlState::erase(size_t overall_lane) {
    // Copy what we're erasing
    auto copy = haplotypes.at(overall_lane);
    
    for (auto it = copy.rbegin(); it != copy.rend(); ++it) {
        // Trace from end to start and remove from the net node lanes collections.
        // We have to do it backward so we can handle duplicate visits properly.
        auto& node_lanes = net_node_lanes[graph->forward(it->first)];
        auto lane_iterator = node_lanes.erase(node_lanes.begin() + it->second);
        
        while (lane_iterator != node_lanes.end()) {
            // Update all the subsequent records in that net node's lane list and bump down their internal lane assignments
            
            // First dereference to get the iterator that points to the actual
            // record, then dereference that, and decrement the lane number.
            (*(*lane_iterator)).second--;
            
            ++lane_iterator;
        }
    }

    // Drop the actual haplotype
    haplotypes.erase(haplotypes.begin() + overall_lane);
    
    // Return the copy
    return copy;
}

void SnarlState::swap(size_t lane1, size_t lane2) {
    
    // Swap the start and end annotation values
    std::swap(haplotypes.at(lane1).front().second, haplotypes.at(lane2).front().second);
    std::swap(haplotypes.at(lane1).back().second, haplotypes.at(lane2).back().second);
    
    // Swap the start net node index entries
    auto& start_node_lanes = net_node_lanes[graph->forward(graph->get_start())];
    std::swap(start_node_lanes.at(lane1), start_node_lanes.at(lane2));
    
    // Swap the end net node index entries
    auto& end_node_lanes = net_node_lanes[graph->forward(graph->get_end())];
    std::swap(end_node_lanes.at(lane1), end_node_lanes.at(lane2));
    
    // Swap the actual haplotype vectors
    std::swap(haplotypes.at(lane1), haplotypes.at(lane2));
}

GenomeStateCommand* InsertHaplotypeCommand::execute(GenomeState& state) const {
    // Allocate and populate the reverse command.
    return new DeleteHaplotypeCommand(state.insert_haplotype(*this));
}

GenomeStateCommand* DeleteHaplotypeCommand::execute(GenomeState& state) const {
    // Allocate and populate the reverse command.
    return new InsertHaplotypeCommand(state.delete_haplotype(*this));
}

GenomeStateCommand* SwapHaplotypesCommand::execute(GenomeState& state) const {
    // Allocate and populate the reverse command.
    return new SwapHaplotypesCommand(state.swap_haplotypes(*this));
}

GenomeStateCommand* AppendHaplotypeCommand::execute(GenomeState& state) const {
    // Allocate and populate the reverse command.
    return new DeleteHaplotypeCommand(state.append_haplotype(*this));
}

GenomeStateCommand* ReplaceSnarlHaplotypeCommand::execute(GenomeState& state) const {
    // Allocate and populate the reverse command.
    return new ReplaceLocalHaplotypeCommand(state.replace_snarl_haplotype(*this));
}

GenomeStateCommand* ReplaceLocalHaplotypeCommand::execute(GenomeState& state) const {
    // Allocate and populate the reverse command.
    return new ReplaceLocalHaplotypeCommand(state.replace_local_haplotype(*this));
}


GenomeState::GenomeState(const SnarlManager& manager, const HandleGraph* graph,
    const unordered_set<pair<const Snarl*, const Snarl*>> telomeres) : telomeres(telomeres),
    backing_graph(graph), manager(manager) {

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

const NetGraph* GenomeState::get_net_graph(const Snarl* snarl) {
    return &net_graphs.at(snarl);
}

DeleteHaplotypeCommand GenomeState::append_haplotype(const AppendHaplotypeCommand& c) {
    // We'll populate this with all the stuff we added
    DeleteHaplotypeCommand to_return;
    
    // Do a stack-based walk of the specified haplotype and thread it in in the
    // final overall lane.
    insert_handles(c.haplotype, to_return.deletions);
    
    // Reverse all the deletion vectors to delete in reverse insertion order
    for (auto& kv : to_return.deletions) {
        reverse(kv.second.begin(), kv.second.end());
    }
    
    return to_return;
}

ReplaceLocalHaplotypeCommand GenomeState::replace_local_haplotype(const ReplaceLocalHaplotypeCommand& c) {
    ReplaceLocalHaplotypeCommand to_return;
    
    // We need to do all the deletions and create balancing insertions
    // TODO: this code is duplicated with delete_haplotype
    for (auto& kv : c.deletions) {
        // We can handle each snarl independently.
        // TODO: do this in parallel?
        auto& snarl = kv.first;
        auto& overall_lanes = kv.second;
        
        // Find where to log the deletions we need to do
        auto& haplotype_insertions = to_return.insertions[snarl];
        
        for (auto& overall_lane : overall_lanes) {
            // For each haplotype we want to remove from this snarl, in order...
            
#ifdef debug
            cerr << "Delete " << overall_lane << " from " << kv.first->start() << " -> " << kv.first->end() << endl;
#endif
            
            // Remove the haplotype and save a copy
            auto removed = state.at(snarl).erase(overall_lane);  
            
            // Save the insertion to do by logging the haplotype with all its
            // tagged lane assignments.
            haplotype_insertions.emplace_back(removed); 
        }
        
        // Flip the insertions around to happen in reverse order. Things need to
        // get to the lanes we deleted them from.
        reverse(haplotype_insertions.begin(), haplotype_insertions.end());
    }
    
    // Then we need to do all the insertions and create balancing deletions
    // TODO: this is duplicated with insert_haplotype
    for (auto& kv : c.insertions) {
        // We can handle each snarl independently.
        // TODO: do this in parallel?
        auto& snarl = kv.first;
        auto& haplotypes = kv.second;
        
        // Find where to log the deletions we need to do
        auto& haplotype_deletions = to_return.deletions[snarl];
        
        for (auto& haplotype : haplotypes) {
            // For each haplotype we want to add to this snarl, in order...
            
            // Insert the haplotype
            state.at(snarl).insert(haplotype);  
            
            // Save the deletion to do by logging the overall lane used.
            haplotype_deletions.emplace_back(haplotype.front().second); 
        }
        
        // Flip the deletions around to happen in reverse order. Things may not
        // stay in the lane we put them in when we add later things.
        reverse(haplotype_deletions.begin(), haplotype_deletions.end());
    }
    
    return to_return;
}

DeleteHaplotypeCommand GenomeState::insert_haplotype(const InsertHaplotypeCommand& c) {
    DeleteHaplotypeCommand to_return;
    
    for (auto& kv : c.insertions) {
        // We can handle each snarl independently.
        // TODO: do this in parallel?
        auto& snarl = kv.first;
        auto& haplotypes = kv.second;
        
        // Find where to log the deletions we need to do
        auto& haplotype_deletions = to_return.deletions[snarl];
        
        for (auto& haplotype : haplotypes) {
            // For each haplotype we want to add to this snarl, in order...
            
            // Insert the haplotype
            state.at(snarl).insert(haplotype);  
            
            // Save the deletion to do by logging the overall lane used.
            haplotype_deletions.emplace_back(haplotype.front().second); 
        }
        
        // Flip the deletions around to happen in reverse order. Things may not
        // stay in the lane we put them in when we add later things.
        reverse(haplotype_deletions.begin(), haplotype_deletions.end());
    }
    
    return to_return;
}

InsertHaplotypeCommand GenomeState::delete_haplotype(const DeleteHaplotypeCommand& c) {
    InsertHaplotypeCommand to_return;
    
    for (auto& kv : c.deletions) {
        // We can handle each snarl independently.
        // TODO: do this in parallel?
        auto& snarl = kv.first;
        auto& overall_lanes = kv.second;
        
        // Find where to log the deletions we need to do
        auto& haplotype_insertions = to_return.insertions[snarl];
        
        for (auto& overall_lane : overall_lanes) {
            // For each haplotype we want to remove from this snarl, in order...
            
#ifdef debug
            cerr << "Delete " << overall_lane << " from " << kv.first->start() << " -> " << kv.first->end() << endl;
#endif
            
            // Remove the haplotype and save a copy
            auto removed = state.at(snarl).erase(overall_lane);  
            
            // Save the insertion to do by logging the haplotype with all its
            // tagged lane assignments.
            haplotype_insertions.emplace_back(removed); 
        }
        
        // Flip the insertions around to happen in reverse order. Things need to
        // get to the lanes we deleted them from.
        reverse(haplotype_insertions.begin(), haplotype_insertions.end());
    }
    
    return to_return;
}

SwapHaplotypesCommand GenomeState::swap_haplotypes(const SwapHaplotypesCommand& c) {
    
    // We have to walk the chromosome and swap in each top-level snarl.
    
    // Make a visit to track where we are. We start at the start of the forward
    // snarl.
    Visit here = c.telomere_pair.first->start();
    
    // Work out what snarl comes next
    const Snarl* next = manager.into_which_snarl(here);
    
    while (next != nullptr) {
        // Until we run out of snarls
    
        // Work out if we go backward or forward through this one
        bool backward = (here.node_id() != next->start().node_id());
        
        // Swap the lanes in this snarl
        state.at(next).swap(c.to_swap.first, c.to_swap.second);
        
        if (next == c.telomere_pair.second) {
            // We just did the last snarl on the chromosome so stop. Don't go
            // around circular things forever.
            break;
        }
        
        // Now look at the visit out of the snarl we just did
        here = backward ? reverse(next->start()) : next->end();
        
        // See if that puts us in another snarl
        next = manager.into_which_snarl(here);
    }

    // This command is its own inverse
    return c;
}

ReplaceLocalHaplotypeCommand GenomeState::replace_snarl_haplotype(const ReplaceSnarlHaplotypeCommand& c) {
    ReplaceLocalHaplotypeCommand to_return;
    
    // First recursively delete everything that was in that lane.
    // We use this map of lanes to delete. We use sets to make sure we delete lanes from high to low.
    unordered_map<const Snarl*, set<size_t>> lanes_to_delete{make_pair(c.snarl, set<size_t>{c.lane})};
    
    while (!lanes_to_delete.empty()) {
        auto& snarl = lanes_to_delete.begin()->first;
        auto& overall_lanes = lanes_to_delete.begin()->second;
        
        // Find where to log the deletions we need to do
        auto& haplotype_insertions = to_return.insertions[snarl];
        
        for (auto it = overall_lanes.rbegin(); it != overall_lanes.rend(); ++it) {
            // For each haplotype we want to remove from this snarl, in decreasign order by lane number...
            auto& overall_lane = *it;
            
#ifdef debug
            cerr << "Delete " << overall_lane << " from " << snarl->start() << " -> " << snarl->end() << endl;
#endif
            
            // Remove the haplotype and save a copy
            auto removed = state.at(snarl).erase(overall_lane);
            
            for (auto& handle_and_lane : removed) {
                if (net_graphs.at(snarl).is_child(handle_and_lane.first)) {
                    // If any of the handles is to a child chain/unary snarl, we
                    // need to remember to remove this lane from all its snarls.
                    
                    // TODO: Should this get chain from net graph handle logic
                    // be factored out into a utility function somehow?
                    
                    // Get the handle in the backing graph that reads into the child
                    // in the orientation we are visiting it
                    handle_t into = net_graphs.at(snarl).get_inward_backing_handle(handle_and_lane.first);
                
                    // Get the child we are actually reading into from the SnarlManager
                    const Snarl* child = manager.into_which_snarl(to_visit(*backing_graph, into));
                
                    // Get the chain for the child
                    const Chain* child_chain = manager.chain_of(child);
                    
                    for (auto it = chain_begin(*child_chain); it != chain_end(*child_chain); ++it) {
                        // For each snarl in the chain, remember to delete this overall lane
                        lanes_to_delete[it->first].insert(handle_and_lane.second);
                    }
                }
            }
            
            // Save the insertion to do by logging the haplotype with all its
            // tagged lane assignments.
            haplotype_insertions.emplace_back(removed);
        }
        
        // Flip the insertions around to happen in reverse order. Things need to
        // get to the lanes we deleted them from.
        reverse(haplotype_insertions.begin(), haplotype_insertions.end());
        
        // We finished this snarl. Handle another snarl.
        lanes_to_delete.erase(snarl);
    }
    
    
    // Then trace along the new haplotype bit stack-wise like we do for an
    // append. Make sure it ends up in the right lane that we just deleted.
    insert_handles(c.haplotype, to_return.deletions, c.lane);
    
    // Reverse all the deletion vectors to delete in reverse insertion order
    for (auto& kv : to_return.deletions) {
        reverse(kv.second.begin(), kv.second.end());
    }
    
    return to_return;
}

GenomeStateCommand* GenomeState::execute(GenomeStateCommand* command) {
    // Just make the command tell us what type it is
    return command->execute(*this);
}

size_t GenomeState::count_haplotypes(const pair<const Snarl*, const Snarl*>& telomere_pair) const {
    // We assume all the traversals go through the whole chromosome from telomere to telomere.
    return state.at(telomere_pair.first).size();
}

size_t GenomeState::count_haplotypes(const Snarl* snarl) const {
    return state.at(snarl).size();
}

void GenomeState::trace_haplotype(const pair<const Snarl*, const Snarl*>& telomere_pair,
    size_t overall_lane, const function<void(const handle_t&)>& iteratee) const {
    
    // We need to traverse this hierarchy while not emitting visits twice. The
    // hard part is that the same handle represents entering a snarl and the
    // start visit of that snarl. The other hard part is that for back-to-back
    // snarls in a chain be need to not visit the shared node twice.
    
    // So what we do is, we never emit the handle for entering a snarl and
    // always make it be the frist handle from inside the snarl instead. Also,
    // if we visit a child snarl immediately after we did a child snarl, we
    // leave the first handle of the child snarl off because we just did it.
    
    // We define a recursive function to traverse a snarl and all its visited children.
    function<void(const Snarl*, size_t, bool, bool)> recursively_traverse_snarl =
        [&](const Snarl* here, size_t lane, bool orientation, bool skip_first) {

        // Here's the NetGraph we are working on
        auto& net_graph = net_graphs.at(here);
    
#ifdef debug
        cerr << "Tracing snarl " << here->start() << " -> " << here->end() << " lane " << lane
            << " in orientation " << orientation << " and skip first flag " << skip_first << endl;
#endif
    
        // Go through its traversal
        state.at(here).trace(lane, orientation, [&](const handle_t& visit, const size_t child_lane) {
            
#ifdef debug
            cerr << "Snarl " << here->start() << " -> " << here->end() << " has visit "
                << net_graph.get_id(visit) << " " << net_graph.get_is_reverse(visit) << endl;
#endif
            
            if (skip_first) {
                // We aren't supposed to do this visit; it's already been done
                // by the previous snarl in our chain.
                skip_first = false;
                return;
            }
                        
            if (net_graph.is_child(visit)) {
                // If the visit enters a real child chain, we have to do that
                // child chain all the snarls in it.
                
                // Get the handle in the backing graph that reads into the chain
                // in the orientation we are visiting it
                handle_t into = net_graph.get_inward_backing_handle(visit);
                
#ifdef debug
                cerr << "\tInward backing handle is " << backing_graph->get_id(into) << " "
                    << backing_graph->get_is_reverse(into) << endl;
#endif
            
                // Get the snarl we really are entering, because get_inward_backing_handle works.
                const Snarl* entered = manager.into_which_snarl(to_visit(*backing_graph, into));
                
                // Decide if we are entering it through its end
                bool entered_snarl_via_end = entered->start().node_id() != backing_graph->get_id(into);
                
                // Get the chain that is our child
                const Chain* child_chain = manager.chain_of(entered);
                
                // Get this snarl's orientation in that chain
                bool entered_orientation_in_chain = manager.chain_orientation_of(entered);
                
                // If we enter a snarl via its start, and it is forward in the chain, we do the whole chain forward.
                // If we enter a snarl via its end, and it is forward in its chain, we do the whole chain reverse.
                // If we enter a snarl via its start, and it is backward in its chain, we do the whole chain reverse.
                // If we enter a snarl via its end, and it is backward in its chain, we do the whole chain forward.
                // F F F
                // T F T
                // F T T
                // T T F
                // So we use xor here.
                bool chain_orientation = entered_snarl_via_end != entered_orientation_in_chain;
                
#ifdef debug
                cerr << "\tEnters into chain containing " << entered->start() << " -> " << entered->end()
                    << " which we enter in orientation " << entered_snarl_via_end
                    << " and is in orientation " << entered_orientation_in_chain << " in its chain" << endl;
                cerr << "\tOverall, chain should be traversed in orientation " << chain_orientation << endl;
#endif
                
                // Get the iterators to loop over the chain in the orientation we have reached it
                auto it = chain_orientation ? chain_rcbegin(*child_chain) : chain_begin(*child_chain);
                auto end = chain_orientation ? chain_rcend(*child_chain) : chain_end(*child_chain);
           
                // Track if we had a previous snarl that would have emitted any shared nodes.
                bool had_previous = false;
                
                // Make sure we have the end of the chain we were expecting.
                assert(it->first == entered);
                
                for (; it != end; ++it) {
                    // Until we run out of snarls in the chain
                
#ifdef debug
                    cerr << "\t\tHandle " << it->first->start() << " -> " << it->first->end()
                        << " orientation " << it->second << " in chain" << endl;
#endif
                
                    // Handle this one
                    recursively_traverse_snarl(it->first, child_lane, it->second, had_previous);
                    
                    // Say we did a snarl previously
                    had_previous = true;
                }                                                               
            } else {
                // This is a visit to a normal handle in this snarl.
                
                // Emit it
                iteratee(visit);
            }
        });
    };

    // Now we need to walk between the telomeres we were given. It's not quite a
    // chain because the telomeres may be unary snarls.
    
    // Make a visit to track where we are. We start at the start of the forward
    // snarl.
    Visit here = telomere_pair.first->start();
    
    // Work out what snarl comes next
    const Snarl* next = manager.into_which_snarl(here);
    
    // Track if we had a previous snarl that would have emitted any shared nodes.
    bool had_previous = false;
    
    while (next != nullptr) {
        // Until we run out of snarls
    
        // Work out if we go backward or forward through this one
        bool backward = (here.node_id() != next->start().node_id());
        
#ifdef debug
        cerr << "Traverse top level snarl " << pb2json(*next) << " lane " << overall_lane
            << " in orientation " << backward << " and skip first flag " << had_previous << endl;
#endif
        
        // Handle it
        recursively_traverse_snarl(next, overall_lane, backward, had_previous);
        
        // Say we did a snarl previously
        had_previous = true;
        
        if (next == telomere_pair.second) {
            // We just did the last snarl on the chromosome so stop. Don't go
            // around circular things forever.
            break;
        }
        
        // Now look at the visit out of the snarl we just did
        here = backward ? reverse(next->start()) : next->end();
        
        // See if that puts us in another snarl
        next = manager.into_which_snarl(here);
    }
}


void GenomeState::dump() const {
    for (auto& kv : state) {
        cerr << "State of " << kv.first->start() << " -> " << kv.first->end() << ":" << endl;
        kv.second.dump();
    }
}

void GenomeState::insert_handles(const vector<handle_t>& to_add,
    unordered_map<const Snarl*, vector<size_t>>& lanes_added, size_t top_lane) {

    // We can't add an empty haplotype.
    assert(!to_add.empty());
    
    // This holds a stack of all the snarls we are in at a given point in the
    // haplotype we are adding, and the handles we are putting for the
    // traversals of them that we are building. Lane assignments are not
    // necessary since they will always be last.
    list<pair<const Snarl*, vector<handle_t>>> stack;

    // We know we're at the start of a telomere snarl, so we can just jump right
    // into the main loop...
    
    for (auto& next_handle : to_add) {
        // For each handle, look at it as a visit in the base graph
        Visit next_visit = to_visit(*backing_graph, next_handle);
        
#ifdef debug
        cerr << "Stack: ";
        for (auto& frame : stack) {
            cerr << frame.first->start() << " -> " << frame.first->end() << ", ";
        }
        cerr << endl;
        
        cerr << "Encountered visit: " << next_visit << endl;
#endif
        
        // Are we going in and out of snarls?
        auto last_snarl = manager.into_which_snarl(reverse(next_visit));
        auto next_snarl = manager.into_which_snarl(next_visit);
        
        if (last_snarl != nullptr) {
            // If we're leaving a child snarl
            
#ifdef debug
            cerr << "Leaving snarl " << last_snarl->start() << " -> " << last_snarl->end() << endl;
#endif
            
            // Make sure it's the one we have been working on
            assert(!stack.empty());
            assert(stack.front().first == last_snarl);
            
            // Make sure the exit handle is in the haplotype
            stack.front().second.push_back(next_handle);
            
            // What state do we have to work on?
            auto& snarl_state = state.at(last_snarl);
            
            // Did we go backward or forward through this child snarl? This only
            // matters to make sure we get the lane assignments right for its
            // repeated children. We keep going forward through the haplotype
            // either way.
            bool backward = next_visit != last_snarl->end();
            
            // Add in its haplotype, and get the resulting lane assignments.
            // Make sure to insert at the right lane if we are the last thing on
            // the stack (i.e. the top level snarl) and have a particular lane.
            auto& embedded = (stack.size() == 1 && top_lane != numeric_limits<size_t>::max()) ?
                snarl_state.insert(top_lane, stack.front().second, backward) :
                snarl_state.append(stack.front().second, backward);
            
            // Remember to delete the overall lane from this snarl
            assert(!embedded.empty());
            lanes_added[last_snarl].push_back(embedded.front().second);
            
            // Pop from the stack
            stack.pop_front();
            
            // What chain are we in?
            auto chain = manager.chain_of(last_snarl);
            
            if ((next_visit == get_end_of(*chain) || next_visit == reverse(get_start_of(*chain))) && !stack.empty()) {
                // If we exited a chain, record a traversal of the whole chain in
                // the parent snarl's haplotype under construction.
                
#ifdef debug
                cerr << "Snarl was last in its chain of " << chain->size() << endl;
#endif
                
                // Get the parent snarl
                const Snarl* parent = stack.front().first;
                // And its net graph
                auto& net_graph = net_graphs.at(parent);
                
                // Get the backing graph handle reading out of the chain
                handle_t backing_outward = backing_graph->get_handle(next_visit.node_id(), next_visit.backward());
                
                // Make it inward
                handle_t backing_inward = backing_graph->flip(backing_outward);
                
                // Get a handle_t representing the whole chain.
                // We get the representative chain handle for the inward backing handle, and flip it back.
                handle_t chain_handle = net_graph.flip(net_graph.get_handle_from_inward_backing_handle(backing_inward));

#ifdef debug
                cerr << "Represent chain traversal with " << to_visit(net_graph, chain_handle) << endl;
#endif

                // Tack it on to the parent
                stack.front().second.push_back(chain_handle);
            } else {
#ifdef debug
                cerr << "Visit " << next_visit << " was in chain of " << chain->size()
                    << " with " << reverse(get_start_of(*chain)) << " out start and "
                    << get_end_of(*chain) << " out end and " << stack.size() << " on stack" << endl;
#endif
            }
            
        } else if (next_snarl == nullptr && !stack.empty()) {
            // Otherwise, we're an ordinary visit in the snarl we're in. So make
            // sure we're in a snarl (i.e. not the first handle in the whole
            // haplotype).
            
#ifdef debug
            cerr << "In snarl " << stack.front().first->start() << " -> " << stack.front().first->end() << endl;
#endif
            
            // Add this handle in the backing graph, which is going to be used
            // to represent a visit to an ordinary node, to the top haplotype on
            // the stack.
            stack.front().second.push_back(next_handle);
        }
        
        if (next_snarl != nullptr) {
            // When we come to a child snarl, descend into a new stack frame.
            // TODO: relies on the backing graph handles being the end handles
            // in the snarl's net graph.
            
            // Note that we may enter a snarl we don't finish if we aren't
            // inserting across all the snarls in a chain.
            
#ifdef debug
            cerr << "Entering snarl " << next_snarl->start() << " -> " << next_snarl->end() << endl;
#endif
            
            stack.emplace_front(next_snarl, vector<handle_t>{next_handle});
        }
    }

}
    

}
