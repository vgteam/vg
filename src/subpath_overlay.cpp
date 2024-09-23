#include <atomic>
#include "subpath_overlay.hpp"

#include <handlegraph/util.hpp>

namespace vg {

using namespace std;
using namespace handlegraph;

SubpathOverlay::SubpathOverlay(const PathHandleGraph* backing,
                               const step_handle_t& begin, const step_handle_t& end) : backing_graph(backing) {
    
    for (auto step = begin; step != end; step = backing->get_next_step(step)) {
        subpath_handles.emplace_back(backing->get_handle_of_step(step));
    }
}

bool SubpathOverlay::has_node(nid_t node_id) const {
    return node_id <= subpath_handles.size();
}
   
handle_t SubpathOverlay::get_handle(const nid_t& node_id, bool is_reverse) const {
    return handlegraph::number_bool_packing::pack(node_id, is_reverse);
}
    
nid_t SubpathOverlay::get_id(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_number(handle);
}
    
bool SubpathOverlay::get_is_reverse(const handle_t& handle) const {
    return handlegraph::number_bool_packing::unpack_bit(handle);
}

handle_t SubpathOverlay::flip(const handle_t& handle) const {
    return handlegraph::number_bool_packing::toggle_bit(handle);
}
    
size_t SubpathOverlay::get_length(const handle_t& handle) const {
    return backing_graph->get_length(subpath_handles[get_id(handle) - 1]);
}

std::string SubpathOverlay::get_sequence(const handle_t& handle) const {
    return backing_graph->get_sequence(get_underlying_handle(handle));
}
    
size_t SubpathOverlay::get_node_count() const {
    return subpath_handles.size();
}

nid_t SubpathOverlay::min_node_id() const {
    return 1;
}
    
nid_t SubpathOverlay::max_node_id() const {
    return subpath_handles.size();
}

bool SubpathOverlay::follow_edges_impl(const handle_t& handle, bool go_left,
                                       const std::function<bool(const handle_t&)>& iteratee) const {
    if (go_left != get_is_reverse(handle)) {
        if (get_id(handle) == 1) {
            return true;
        }
        else {
            return iteratee(get_handle(get_id(handle) - 1, get_is_reverse(handle)));
        }
    }
    else {
        if (get_id(handle) == subpath_handles.size()) {
            return true;
        }
        else {
            return iteratee(get_handle(get_id(handle) + 1, get_is_reverse(handle)));
        }
    }
}
    
bool SubpathOverlay::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {

    if (!parallel) {
        bool keep_going = true;
        for (size_t i = 1; keep_going && i <= subpath_handles.size(); ++i) {
            keep_going = iteratee(get_handle(i, false));
        }
        return keep_going;
    } else {
        std::atomic<bool> keep_going(true);
#pragma omp parallel for
        for (size_t i = 1; i <= subpath_handles.size(); ++i) {
            keep_going = keep_going && iteratee(get_handle(i, false));
        }
        return keep_going;
    }
}

handle_t SubpathOverlay::get_underlying_handle(const handle_t& handle) const {
    handle_t underlying = subpath_handles[get_id(handle) - 1];
    if (get_is_reverse(handle)) {
        underlying = backing_graph->flip(underlying);
    }
    return underlying;
}

}
