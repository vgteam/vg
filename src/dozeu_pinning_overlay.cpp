/**
 * \file dozeu_pinning_overlay.cpp: contains the implementation of DozeuPinningOverlay
 */


#include "dozeu_pinning_overlay.hpp"


namespace vg {

using namespace std;

DozeuPinningOverlay::DozeuPinningOverlay(const HandleGraph* graph, bool preserve_sinks) : graph(graph), preserve_sinks(preserve_sinks) {
    
    unordered_set<handle_t> empty_nodes;
        
    // find the numeric range of handles in the underlying graph (needed for later bookkeeping)
    // and all nodes with no sequence
    uint64_t min_handle = std::numeric_limits<uint64_t>::max();
    graph->for_each_handle([&](const handle_t& handle) {
        
        min_handle = min<uint64_t>(handlegraph::as_integer(handle), min_handle);
        min_handle = min<uint64_t>(handlegraph::as_integer(graph->flip(handle)), min_handle);
        max_handle = max<uint64_t>(handlegraph::as_integer(handle), max_handle);
        max_handle = max<uint64_t>(handlegraph::as_integer(graph->flip(handle)), max_handle);
        
        if (graph->get_length(handle) == 0) {
            empty_nodes.insert(handle);
        }
    });
    // keep track of these values
    num_null_nodes = empty_nodes.size();
    handle_val_range = max_handle - min_handle + 1;
        
    for (const handle_t& empty : empty_nodes) {
        
        // is this empty node a source/sink (depending on what we're preserving)?
        bool should_preserve = graph->follow_edges(empty, !preserve_sinks, [&](const handle_t& next) { return false; });
                
        if (should_preserve) {
            // check the neighbors of the empty source/sink
            graph->follow_edges(empty, preserve_sinks, [&](const handle_t& next) {
                // walk back in the direction of the empty node and check if there are any different neighbors
                bool must_duplicate = !graph->follow_edges(next, !preserve_sinks, [&](const handle_t& prev) { return prev == empty; });
                
                if (must_duplicate) {
                    // this node has a length 0 path from a source/sink, but it will not once we mask out the
                    // length 0 nodes, so we have to duplicate it so that it retains its ability to be used as
                    // a pinning point
                    duplicated_handles.insert(next);
                }
            });
        }
    }
    num_null_nodes = empty_nodes.size();
}

bool DozeuPinningOverlay::performed_duplications() const {
    return !duplicated_handles.empty();
}

bool DozeuPinningOverlay::has_node(id_t node_id) const {
    if (is_a_duplicate_id(node_id)) {
        id_t under_id = get_underlying_id(node_id);
        if (graph->has_node(under_id)) {
            handle_t handle = graph->get_handle(under_id);
            return duplicated_handles.count(handle);
        }
    }
    else {
        if (graph->has_node(node_id)) {
            handle_t handle = graph->get_handle(node_id);
            return graph->get_length(handle) != 0;
        }
    }
    return false;
}

handle_t DozeuPinningOverlay::get_handle(const id_t& node_id, bool is_reverse) const {
    if (is_a_duplicate_id(node_id)) {
        return get_duplicate_handle(graph->get_handle(get_underlying_id(node_id), is_reverse));
    }
    else {
        return graph->get_handle(node_id, is_reverse);
    }
}

id_t DozeuPinningOverlay::get_id(const handle_t& handle) const {
    if (is_a_duplicate_handle(handle)) {
        return graph->get_id(get_underlying_handle(handle)) + (graph->max_node_id() - graph->min_node_id() + 1);
    }
    else {
        return graph->get_id(handle);
    }
}

bool DozeuPinningOverlay::get_is_reverse(const handle_t& handle) const {
    if (is_a_duplicate_handle(handle)) {
        return graph->get_is_reverse(get_underlying_handle(handle));
    }
    else {
        return graph->get_is_reverse(handle);
    }
}

handle_t DozeuPinningOverlay::flip(const handle_t& handle) const {
    if (is_a_duplicate_handle(handle)) {
        return get_duplicate_handle(graph->flip(get_underlying_handle(handle)));
    }
    else {
        return graph->flip(handle);
    }
}

size_t DozeuPinningOverlay::get_length(const handle_t& handle) const {
    if (is_a_duplicate_handle(handle)) {
        return graph->get_length(get_underlying_handle(handle));
    }
    else {
        return graph->get_length(handle);
    }
}

string DozeuPinningOverlay::get_sequence(const handle_t& handle) const {
    if (is_a_duplicate_handle(handle)) {
        return graph->get_sequence(get_underlying_handle(handle));
    }
    else {
        return graph->get_sequence(handle);
    }
}

bool DozeuPinningOverlay::follow_edges_impl(const handle_t& handle, bool go_left,
                                         const function<bool(const handle_t&)>& iteratee) const {
    handle_t to_iterate = handle;
    if (is_a_duplicate_handle(handle)) {
        // this is a duplicate that we made to preserve sources/sinks
        if (preserve_sinks != (go_left != graph->get_is_reverse(handle))) {
            // we are going in the direction where the duplicate has no edges
            return true;
        }
        to_iterate = get_underlying_handle(handle);
    }
    
    return graph->follow_edges(to_iterate, go_left, [&](const handle_t& next) {
        bool keep_going = true;
        if (graph->get_length(next) > 0) {
            // the node is non-empty, so it hasn't been removed
            keep_going = iteratee(next);
            if (keep_going && duplicated_handles.count(graph->forward(next))) {
                // the node has a duplicate
                if (preserve_sinks != (go_left != graph->get_is_reverse(next))) {
                    // we arrived at this node over an edge that the duplicate also shares
                    keep_going = iteratee(get_duplicate_handle(next));
                }
            }
        }
        return keep_going;
    });
}

bool DozeuPinningOverlay::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    // iterate over the original non empty nodes
    bool keep_going = graph->for_each_handle([&](const handle_t& handle) {
        if (graph->get_length(handle) > 0) {
            return iteratee(handle);
        }
        else {
            return true;
        }
    }, parallel);
    // iterate over the duplicates
    for (auto it = duplicated_handles.begin(), end = duplicated_handles.end(); it != end && keep_going; ++it) {
        keep_going = iteratee(get_duplicate_handle(*it));
    }
    return keep_going;
}

size_t DozeuPinningOverlay::get_node_count() const {
    return graph->get_node_count() - num_null_nodes + duplicated_handles.size();
}

id_t DozeuPinningOverlay::min_node_id() const {
    return graph->min_node_id();
}

id_t DozeuPinningOverlay::max_node_id() const {
    id_t max_id = graph->max_node_id();
    for (const handle_t& handle : duplicated_handles) {
        max_id = max(max_id, get_id(get_duplicate_handle(handle)));
    }
    return max_id;
}

bool DozeuPinningOverlay::is_a_duplicate_handle(const handle_t& handle) const {
    return max_handle < (uint64_t) handlegraph::as_integer(handle);
}

bool DozeuPinningOverlay::is_a_duplicate_id(const id_t& node_id) const {
    return node_id > graph->max_node_id();
}


handle_t DozeuPinningOverlay::get_underlying_handle(const handle_t& handle) const {
    if (!is_a_duplicate_handle(handle)) {
        return handle;
    }
    else {
        return handlegraph::as_handle(uint64_t(handlegraph::as_integer(handle)) - handle_val_range);
    }
}

id_t DozeuPinningOverlay::get_underlying_id(const id_t& node_id) const {
    return node_id - (graph->max_node_id() - graph->min_node_id() + 1);
}
 

handle_t DozeuPinningOverlay::get_duplicate_handle(const handle_t& handle) const {
    return handlegraph::as_handle(uint64_t(handlegraph::as_integer(handle)) + handle_val_range);
}

}
