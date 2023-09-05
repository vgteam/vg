/**
 * \file identity_overlay.cpp: contains the implementation of IdentityOverlay
 */


#include "identity_overlay.hpp"


namespace vg {

using namespace std;

    IdentityOverlay::IdentityOverlay(const HandleGraph* graph) : graph(graph) {

    }
    
    bool IdentityOverlay::has_node(id_t node_id) const {
        return graph->has_node(node_id);
    }
    
    handle_t IdentityOverlay::get_handle(const id_t& node_id, bool is_reverse) const {
        return graph->get_handle(node_id, is_reverse);
    }
    
    id_t IdentityOverlay::get_id(const handle_t& handle) const {
        return graph->get_id(handle);
    }
    
    bool IdentityOverlay::get_is_reverse(const handle_t& handle) const {
        return graph->get_is_reverse(handle);
    }
    
    handle_t IdentityOverlay::flip(const handle_t& handle) const {
        return graph->flip(handle);
    }
    
    size_t IdentityOverlay::get_length(const handle_t& handle) const {
        return graph->get_length(handle);
    }
    
    string IdentityOverlay::get_sequence(const handle_t& handle) const {
        return graph->get_sequence(handle);
    }
    
    bool IdentityOverlay::follow_edges_impl(const handle_t& handle, bool go_left,
                                      const function<bool(const handle_t&)>& iteratee) const {
        return graph->follow_edges(handle, go_left, iteratee);
    }
    
    bool IdentityOverlay::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                          bool parallel) const {
        return graph->for_each_handle(iteratee, parallel);
    }
    
    size_t IdentityOverlay::get_node_count() const {
        return graph->get_node_count();
    }
    
    id_t IdentityOverlay::min_node_id() const {
        return graph->min_node_id();
    }
    
    id_t IdentityOverlay::max_node_id() const {
        return graph->max_node_id();
    }

    handle_t IdentityOverlay::get_underlying_handle(const handle_t& handle) const {
        return handle;
    }
}

