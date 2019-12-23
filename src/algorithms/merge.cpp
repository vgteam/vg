/**
 * \file merge.cpp
 *
 * Defines an algorithm to merge parts of handles together.
 */

#include "merge.hpp"

#include <unordered_set>
#include <iostream>

namespace vg {
namespace algorithms {

using namespace std;

void merge(handlegraph::MutablePathDeletableHandleGraph* graph, const vector<pair<handle_t, size_t>>& start, size_t length) {

#ifdef debug
    cerr << "Merge " << length << " bases on:" << endl;
    for(auto& handle_and_offset : start) {
        cerr << "\t" << graph->get_id(handle_and_offset.first)
            << (graph->get_is_reverse(handle_and_offset.first) ? '-' : '+')
            << "@" << handle_and_offset.second << endl;
    }
#endif

    if (start.empty()) {
        // Nothing to do!
        return;
    }

    // Split out the desired ranges of all the handles and get the middle bits
    vector<handle_t> middles;
    middles.reserve(start.size());
    
    for (auto& handle_and_offset : start) {
        auto& to_split = handle_and_offset.first;
        auto& offset = handle_and_offset.second;
        
        vector<size_t> divisions;
        // By default the part we are interested in will be at index 0.
        size_t output_index = 0;
        if (offset != 0) {
            // We need to break after the start of the node
            divisions.push_back(offset);
            // Our middle handle will end up at index 1
            output_index = 1;
        }
        if (offset + length != graph->get_length(to_split)) {
            // We need to break before the end of the node
            divisions.push_back(offset + length);
        }
        
        if (divisions.empty()) {
            // Just use the whole node
            middles.push_back(to_split);
        } else {
            // Actually split
            auto parts = graph->divide_handle(to_split, divisions);
            middles.push_back(parts.at(output_index));
        }
    }
    
    // Pick one to be the one that survives
    handle_t merged = middles.back();
    middles.pop_back();
    
    // Create edges from everything attached to the other ones to the first one.
    // We collect and then create edges to avoid upsetting iteration.
    unordered_set<handle_t> right_neighbors;
    unordered_set<handle_t> left_neighbors;
    
    for (auto& other : middles) {
        // For each node we merge in
        graph->follow_edges(other, false, [&](const handle_t& h) {
            // Look right and collect neighbors
            right_neighbors.insert(h);
        });
        graph->follow_edges(other, true, [&](const handle_t& h) {
            // Look left and collect neighbors
            left_neighbors.insert(h);
        });
    }
    
    for (auto& h : right_neighbors) {
        // Make all the right edges. May exist already.
        graph->create_edge(merged, h);
    }
    for (auto& h : left_neighbors) {
        // Make all the left edges. May exist already.
        graph->create_edge(h, merged);
    }
    
    // Move all the paths over to the first one
    // Need to aggregate to avoid removing steps as we visit them.
    // Also need to record orientation so we can preserve it
    vector<pair<step_handle_t, bool>> to_move;
    for (auto& other : middles) {
        graph->for_each_step_on_handle(other, [&](const step_handle_t s) {
            // Say we need to rewrite this step, and record an orientation:
            // true if the step runs against the orientartion of other, and
            // false if it runs with it.
            to_move.emplace_back(s, graph->get_is_reverse(other) != graph->get_is_reverse(graph->get_handle_of_step(s)));
        });
    }
    for (auto& step_and_orientation : to_move) {
        // For each thing we are moving, unpack it
        auto& step = step_and_orientation.first;
        auto& flip = step_and_orientation.second;
        // Rewrite the path to go through merged forward if we went through the
        // hendle we're merging in forward, and merged reverse otherwise.
        graph->rewrite_segment(step, step, {flip ? graph->flip(merged) : merged});
    }
    
    for (auto& other : middles) {
        // Delete the other versions of the merged segment.
        // First we have to delete each edge exactly once
        unordered_set<edge_t> to_remove;
        graph->follow_edges(other, false, [&](const handle_t& h) {
            to_remove.insert(graph->edge_handle(other, h));
        });
        graph->follow_edges(other, true, [&](const handle_t& h) {
            to_remove.insert(graph->edge_handle(h, other));
        });
        for (auto& e : to_remove) {
            graph->destroy_edge(e);
        }
        // And then the node itself
        graph->destroy_handle(other);
    }
}


}
}

