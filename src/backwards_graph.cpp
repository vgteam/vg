/**
 * \file backwards_graph.cpp: contains the implementation of BackwardsGraph
 */


#include "backwards_graph.hpp"


namespace vg {

using namespace std;

    BackwardsGraph::BackwardsGraph(const HandleGraph* forward_graph) : forward_graph(forward_graph) {
        // nothing to do
    }
    
    bool BackwardsGraph::has_node(id_t node_id) const {
        return forward_graph->has_node(node_id);
    }
    
    handle_t BackwardsGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return forward_graph->get_handle(node_id, is_reverse);
    }
    
    id_t BackwardsGraph::get_id(const handle_t& handle) const {
        return forward_graph->get_id(handle);
    }
    
    bool BackwardsGraph::get_is_reverse(const handle_t& handle) const {
        return forward_graph->get_is_reverse(handle);
    }
    
    handle_t BackwardsGraph::flip(const handle_t& handle) const {
        return forward_graph->flip(handle);
    }
    
    size_t BackwardsGraph::get_length(const handle_t& handle) const {
        return forward_graph->get_length(handle);
    }
    
    string BackwardsGraph::get_sequence(const handle_t& handle) const {
        // reverse (non-complementing) the sequence
        string sequence = forward_graph->get_sequence(handle);
        reverse(sequence.begin(), sequence.end());
        return sequence;
    }
    
    bool BackwardsGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                           const function<bool(const handle_t&)>& iteratee) const {
        // the left and right side have been switched, so reverse the direction
        return forward_graph->follow_edges(handle, !go_left, iteratee);
    }
    
    bool BackwardsGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
        return forward_graph->for_each_handle(iteratee, parallel);
    }
    
    size_t BackwardsGraph::node_size() const {
        return forward_graph->node_size();
    }
    
    id_t BackwardsGraph::min_node_id() const {
        return forward_graph->min_node_id();
    }
    
    id_t BackwardsGraph::max_node_id() const {
        return forward_graph->max_node_id();
    }

}

