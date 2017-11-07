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
        throw runtime_error("Tried to add a haplotype to a snarl that is not a start-to-end traversal of that snarl.");
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
        inserted_iterator->second = node_lanes.size() - 1;
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
            inserted_iterator->second = node_lanes.size() - 1;
            // And do the insert
            node_lanes.emplace_back(inserted_iterator);
        }
        
        // Insert the next handle in the next slot in the haplotype
        ++inserted_iterator;
    } 
    
    // Return the annotated haplotype.
    return inserted;
}

GenomeStateCommand* DeleteHaplotypeCommand::execute(GenomeState& state) const {
    // Allocate and populate the reverse command.
    return new InsertHaplotypeCommand(state.delete_haplotype(*this));
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

InsertHaplotypeCommand GenomeState::delete_haplotype(const DeleteHaplotypeCommand& c) {
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
