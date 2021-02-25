// SPDX-FileCopyrightText: 2014 Erik Garrison
//
// SPDX-License-Identifier: MIT

#include "sorted_id_ranges.hpp"

namespace vg {
namespace algorithms {

using namespace std;

vector<pair<id_t, id_t>> sorted_id_ranges(const HandleGraph* graph) {

    // Build the list of of all the node IDs to operate on
    vector<id_t> graph_ids;
    graph->for_each_handle([&](handle_t handle) {
        // Put all the ids in the list
        graph_ids.push_back(graph->get_id(handle));
    });
    
    // Sort the graph IDs
    std::sort(graph_ids.begin(), graph_ids.end());
    
    // Coalesce them into ranges
    vector<pair<id_t, id_t>> ranges;
    for (auto& id : graph_ids) {
        if (ranges.empty() || ranges.back().second + 1 != id) {
            // We can't glom on to the previous range, so start a new one of just us
            ranges.emplace_back(id, id);
        } else {
            // Extend the previous range to us
            ranges.back().second = id;
        }
    }
    
    return ranges;
}


}
}
