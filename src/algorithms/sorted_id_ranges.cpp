#include "sorted_id_ranges.hpp"

namespace vg {
namespace algorithms {

using namespace std;

std::vector<std::pair<nid_t, nid_t>> sorted_id_ranges(const HandleGraph* graph) {

    // Build the list of of all the node IDs to operate on
    std::vector<nid_t> graph_ids;
    graph->for_each_handle([&](handle_t handle) {
        // Put all the ids in the list
        graph_ids.push_back(graph->get_id(handle));
    });
    
    // Sort the graph IDs
    std::sort(graph_ids.begin(), graph_ids.end());
    
    // Coalesce them into ranges
    std::vector<std::pair<nid_t, nid_t>> ranges;
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

bool is_in_sorted_id_ranges(const nid_t& value, const std::vector<std::pair<nid_t, nid_t>>& ranges) {
    // Find the first interval starting after the value
    auto it = std::upper_bound(ranges.begin(), ranges.end(), std::make_pair(value, value), [&](const std::pair<nid_t, nid_t>& a, const std::pair<nid_t, nid_t>& b) {
        return a.first < b.first;
    });
    if (it == ranges.begin()) {
        // All intervals start after the value, so none can cover it.
        return false;
    }
    // Move back to the range covering the value, if there is one.
    --it;
    // Return whether that range contains the value.
    return it->first <= value && it->second >= value;
}


}
}
