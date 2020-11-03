/**
 * \file null_masking_graph.cpp: contains the implementation of NullMaskingGraph
 */


#include "null_masking_graph.hpp"


namespace vg {

using namespace std;

NullMaskingGraph::NullMaskingGraph(const HandleGraph* graph) : graph(graph) {
    graph->for_each_handle([&](const handle_t& handle) {
        if (graph->get_length(handle) == 0) {
            num_null_nodes++;
        }
    });
}

bool NullMaskingGraph::has_node(id_t node_id) const {
    bool found_node = false;
    if (graph->has_node(node_id)) {
        if (graph->get_length(graph->get_handle(node_id)) > 0) {
            found_node = true;
        }
    }
    return found_node;
}

handle_t NullMaskingGraph::get_handle(const id_t& node_id, bool is_reverse) const {
    // TODO: should we throw an assert that it's non-empty?
    return graph->get_handle(node_id, is_reverse);
}

id_t NullMaskingGraph::get_id(const handle_t& handle) const {
    return graph->get_id(handle);
}

bool NullMaskingGraph::get_is_reverse(const handle_t& handle) const {
    return graph->get_is_reverse(handle);
}

handle_t NullMaskingGraph::flip(const handle_t& handle) const {
    return graph->flip(handle);
}

size_t NullMaskingGraph::get_length(const handle_t& handle) const {
    return graph->get_length(handle);
}

string NullMaskingGraph::get_sequence(const handle_t& handle) const {
    return graph->get_sequence(handle);
}

bool NullMaskingGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                         const function<bool(const handle_t&)>& iteratee) const {
    
    return graph->follow_edges(handle, go_left, [&](const handle_t& next) {
        bool keep_going = true;
        if (graph->get_length(next) > 0) {
            keep_going = iteratee(next);
        }
        return keep_going;
    });
}

bool NullMaskingGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
    return graph->for_each_handle([&](const handle_t& handle) {
        bool keep_going = true;
        if (graph->get_length(handle) > 0) {
            keep_going = iteratee(handle);
        }
        return keep_going;
    }, parallel);
}

size_t NullMaskingGraph::get_node_count() const {
    return graph->get_node_count() - num_null_nodes;
}

id_t NullMaskingGraph::min_node_id() const {
    return graph->min_node_id();
}

id_t NullMaskingGraph::max_node_id() const {
    return graph->max_node_id();
}

handle_t NullMaskingGraph::get_underlying_handle(const handle_t& handle) const {
    return handle;
}

}

