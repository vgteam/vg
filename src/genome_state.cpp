#include "genome_state.hpp"

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
    assert(!haplotype.empty());
    
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

const vector<pair<handle_t, size_t>>& SnarlState::append(const vector<handle_t>& haplotype) {
    assert(!haplotype.empty());
    
    if (haplotype.front() != graph->get_start() || haplotype.back() != graph->get_end()) {
        // Fail if it's not actually from start to end.
        throw runtime_error("Tried to add a haplotype to a snarl that is not a start-to-end traversal of that snarl.");
    }
    
    // Make a new haplotype at the end of our haplotypes vector that's big enough.
    haplotypes.emplace_back(haplotype.size());
    auto& inserted = haplotypes.back();
    
    // Make an iterator to run through it
    auto inserted_iterator = inserted.begin();
    for (auto& handle : haplotype) {
        // For every handle we need to insert
        
        // Save the handle
        inserted_iterator->first = handle;
        
        // Find the appropriate node lanes collection
        auto& node_lanes = net_node_lanes[graph->forward(handle)];
        // Save the local lane assignment
        inserted_iterator->second = node_lanes.size();
        // And do the insert
        node_lanes.emplace_back(inserted_iterator);
        
        // Insert the next handle in the next slot in the haplotype
        ++inserted_iterator;
    }

    // Return the completed vector with the lane annotations.
    return inserted;
}

const vector<pair<handle_t, size_t>>& SnarlState::insert(size_t overall_lane, const vector<handle_t>& haplotype) {
    assert(!haplotype.empty());
    
    if (haplotype.front() != graph->get_start() || haplotype.back() != graph->get_end()) {
        // Fail if it's not actually from start to end.
        throw runtime_error("Tried to add a haplotype to a snarl that is not a start-to-end traversal of that snarl.");
    }
    
    // Insert a haplotype record at the specified overall lane that's big enough.
    auto& inserted = *haplotypes.emplace(haplotypes.begin() + overall_lane, haplotype.size());
    
    // Make an iterator to run through it
    auto inserted_iterator = inserted.begin();
    for (auto& handle : haplotype) {
        // For every handle we need to insert
        
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

GenomeStateCommand* CreateHaplotypeCommand::execute(GenomeState& state) const {
    // Allocate and populate the reverse command.
    return new DeleteHaplotypeCommand(state.create_haplotype(*this));
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

const NetGraph* GenomeState::get_net_graph(const Snarl* snarl) {
    return &net_graphs.at(snarl);
}

DeleteHaplotypeCommand GenomeState::create_haplotype(const CreateHaplotypeCommand& c) {
    // Sample a new haplotype recursively
    
    // TODO: what RNG are we going to use?
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

GenomeStateCommand* GenomeState::execute(GenomeStateCommand* command) {
    // Just make the command tell us what type it is
    return command->execute(*this);
}

size_t GenomeState::count_haplotypes(const pair<const Snarl*, const Snarl*>& telomere_pair) {
    // We assume all the traversals go through the whole chromosome from telomere to telomere.
    return state.at(telomere_pair.first).size();
}

void GenomeState::trace_haplotype(const pair<const Snarl*, const Snarl*>& telomere_pair,
    size_t overall_lane, const function<void(const handle_t&)>& iteratee) {
    
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

        // We need to remember if the last visit we did was a child snarl.
        bool last_was_child = false;
    
        // Here's the NetGraph we are working on
        auto& net_graph = net_graphs.at(here);
    
        // Go through its traversal
        state.at(here).trace(lane, orientation, [&](const handle_t& visit, const size_t child_lane) {
            
            if (skip_first) {
                // We aren't supposed to do this visit; it's already been done.
                skip_first = false;
                return;
            }
            
            // What child, if any, do we enter on this node?
            const Snarl* child = manager.into_which_snarl(net_graph.get_id(visit), net_graph.get_is_reverse(visit));
            
            if (child != nullptr && child != here) {
                // If the visit enters a snarl (other than this one)
                
                // Work out if we go through it backward
                bool child_backward = net_graph.get_id(visit) != child->start().node_id();
                
                // Recurse and do the snarl. Do it in the lane we have recorded
                // for this child visit, in the orientation appropriate for the
                // handle we are entering with, and skipping the first node if
                // the last thing we visited was also a child snarl (which must
                // be back to back with this one).
                recursively_traverse_snarl(child, child_lane, child_backward, last_was_child);
                
                // Remember that the last thing we did was a child snarl.
                last_was_child = true;
            } else {
                // This is a visit to a normal handle in this snarl.
                
                // Emit it
                iteratee(visit);
                
                // Rember we didn't visit a child last at this level
                last_was_child = false;
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
    
    // Track if we had a previous snarl that would have emitted any shaed nodes.
    bool had_previous = false;
    
    while (next != nullptr) {
        // Until we run out of snarls
    
        // Work out if we go backward or forward through this one
        bool backward = (here.node_id() != next->start().node_id());
        
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


}
