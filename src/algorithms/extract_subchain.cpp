#include "extract_subchain.hpp"

#include <stack>

namespace vg {

//------------------------------------------------------------------------------

hash_set<nid_t> extract_subchain(const HandleGraph& graph, handle_t start, handle_t end) {
    hash_set<nid_t> result;
    nid_t start_id = graph.get_id(start);
    bool start_is_reverse = graph.get_is_reverse(start);
    nid_t end_id = graph.get_id(end);
    bool end_is_reverse = graph.get_is_reverse(end);

    std::stack<nid_t> active;
    active.push(start_id);
    active.push(end_id);
    while (!active.empty()) {
        nid_t node_id = active.top(); active.pop();
        if (result.count(node_id)) {
            continue;
        }
        result.insert(node_id);
        bool go_forward = true, go_backward = true;
        if (node_id == start_id) {
            go_forward &= !start_is_reverse;
            go_backward &= start_is_reverse;

        }
        if (node_id == end_id) {
            go_forward &= end_is_reverse;
            go_backward &= !end_is_reverse;
        }
        handle_t handle = graph.get_handle(node_id, false);
        if (go_forward) {
            graph.follow_edges(handle, false, [&](const handle_t& next) {
                active.push(graph.get_id(next));
            });
        }
        if (go_backward) {
            graph.follow_edges(handle, true, [&](const handle_t& prev) {
                active.push(graph.get_id(prev));
            });
        }
    }

    return result;
}

//------------------------------------------------------------------------------

} // namespace vg
