/**
 * \file reverse_graph.cpp: contains the implementation of DagifiedGraph
 */


#include "reverse_graph.hpp"


namespace vg {

using namespace std;

    DagifiedGraph::DagifiedGraph(const HandleGraph* graph) :
        graph(graph), complement(complement) {
        // nothing to do
    }
    
    bool DagifiedGraph::has_node(id_t node_id) const {
        return graph->has_node(node_id);
    }
    
    handle_t DagifiedGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return graph->get_handle(node_id, is_reverse);
    }
    
    id_t DagifiedGraph::get_id(const handle_t& handle) const {
        return graph->get_id(handle);
    }
    
    bool DagifiedGraph::get_is_reverse(const handle_t& handle) const {
        return graph->get_is_reverse(handle);
    }
    
    handle_t DagifiedGraph::flip(const handle_t& handle) const {
        return graph->flip(handle);
    }
    
    size_t DagifiedGraph::get_length(const handle_t& handle) const {
        return graph->get_length(handle);
    }
    
    string DagifiedGraph::get_sequence(const handle_t& handle) const {
        // reverse, possibly complement, the sequence
        string sequence = graph->get_sequence(handle);
        if (complement) {
            reverse_complement_in_place(sequence);
        }
        else {
            reverse(sequence.begin(), sequence.end());
        }
        return sequence;
    }
    
    bool DagifiedGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                         const function<bool(const handle_t&)>& iteratee) const {
        // the left and right side have been switched, so reverse the direction
        return graph->follow_edges(handle, !go_left, iteratee);
    }
    
    bool DagifiedGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                            bool parallel) const {
        // since the handles are the same we can just execute this
        return graph->for_each_handle(iteratee, parallel);
    }
    
    size_t DagifiedGraph::get_node_count() const {
        return graph->get_node_count();
    }
    
    id_t DagifiedGraph::min_node_id() const {
        return graph->min_node_id();
    }
    
    id_t DagifiedGraph::max_node_id() const {
        return graph->max_node_id();
    }

    handle_t DagifiedGraph::get_underlying_handle(const handle_t& handle) const {
        return handle;
    }
}

