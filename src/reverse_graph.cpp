/**
 * \file reverse_graph.cpp: contains the implementation of ReverseGraph
 */


#include "reverse_graph.hpp"


namespace vg {

using namespace std;

    ReverseGraph::ReverseGraph(const HandleGraph* forward_graph, bool complement) :
        forward_graph(forward_graph), complement(complement) {
        // nothing to do
    }
    
    bool ReverseGraph::has_node(id_t node_id) const {
        return forward_graph->has_node(node_id);
    }
    
    handle_t ReverseGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return forward_graph->get_handle(node_id, is_reverse);
    }
    
    id_t ReverseGraph::get_id(const handle_t& handle) const {
        return forward_graph->get_id(handle);
    }
    
    bool ReverseGraph::get_is_reverse(const handle_t& handle) const {
        return forward_graph->get_is_reverse(handle);
    }
    
    handle_t ReverseGraph::flip(const handle_t& handle) const {
        return forward_graph->flip(handle);
    }
    
    size_t ReverseGraph::get_length(const handle_t& handle) const {
        return forward_graph->get_length(handle);
    }
    
    string ReverseGraph::get_sequence(const handle_t& handle) const {
        // reverse, possibly complement, the sequence
        string sequence = forward_graph->get_sequence(handle);
        if (complement) {
            reverse_complement_in_place(sequence);
        }
        else {
            reverse(sequence.begin(), sequence.end());
        }
        return sequence;
    }
    
    bool ReverseGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                         const function<bool(const handle_t&)>& iteratee) const {
        // the left and right side have been switched, so reverse the direction
        return forward_graph->follow_edges(handle, !go_left, iteratee);
    }
    
    bool ReverseGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                            bool parallel) const {
        // since the handles are the same we can just execute this
        return forward_graph->for_each_handle(iteratee, parallel);
    }
    
    size_t ReverseGraph::get_node_count() const {
        return forward_graph->get_node_count();
    }
    
    id_t ReverseGraph::min_node_id() const {
        return forward_graph->min_node_id();
    }
    
    id_t ReverseGraph::max_node_id() const {
        return forward_graph->max_node_id();
    }

    handle_t ReverseGraph::get_underlying_handle(const handle_t& handle) const {
        return flip(handle);
    }
}

