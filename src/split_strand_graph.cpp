/**
 * \file split_strand_graph.cpp: contains the implementation of StrandSplitGraph
 */


#include "split_strand_graph.hpp"


namespace vg {

using namespace std;

    StrandSplitGraph::StrandSplitGraph(const HandleGraph* graph) : graph(graph){
        // nothing to do
    }
    
    bool StrandSplitGraph::has_node(id_t node_id) const {
        return graph->has_node(node_id >> 1);
    }
    
    handle_t StrandSplitGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return handlegraph::number_bool_packing::pack(node_id, is_reverse);
    }
    
    id_t StrandSplitGraph::get_id(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_number(handle);
    }
    
    bool StrandSplitGraph::get_is_reverse(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_bit(handle);
    }
    
    handle_t StrandSplitGraph::flip(const handle_t& handle) const {
        return handlegraph::number_bool_packing::toggle_bit(handle);
    }
    
    size_t StrandSplitGraph::get_length(const handle_t& handle) const {
        return graph->get_length(graph->get_handle(get_id(handle) >> 1));
    }
    
    string StrandSplitGraph::get_sequence(const handle_t& handle) const {
        return graph->get_sequence(get_underlying_handle(handle));
    }
    
    bool StrandSplitGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                             const function<bool(const handle_t&)>& iteratee) const {
        
        return graph->follow_edges(get_underlying_handle(handle), go_left,
                                   [&] (const handle_t& next) {
            return iteratee(get_handle((graph->get_id(next) << 1) + (graph->get_is_reverse(next) != get_is_reverse(handle)),
                                       get_is_reverse(handle)));
        });
    }
    
    bool StrandSplitGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                                bool parallel) const {
        return graph->for_each_handle([&](const handle_t& underlying_handle) {
            id_t node_id = graph->get_id(underlying_handle);
            // forward version of the node
            bool keep_going = iteratee(get_handle(node_id << 1));
            // reverse version of the node
            if (keep_going) {
                keep_going = iteratee(get_handle((node_id << 1) | 1));
            }
            return keep_going;
        }, parallel);
    }
    
    size_t StrandSplitGraph::get_node_count() const {
        return graph->get_node_count() << 1;
    }
    
    id_t StrandSplitGraph::min_node_id() const {
        return graph->min_node_id() << 1;
    }
    
    id_t StrandSplitGraph::max_node_id() const {
        return (graph->max_node_id() << 1) | 1;
    }

    handle_t StrandSplitGraph::get_underlying_handle(const handle_t& handle) const {
        return graph->get_handle(get_id(handle) >> 1,
                                 ((get_id(handle) & 1) == 1) != get_is_reverse(handle));
    }
    
    bool StrandSplitGraph::has_overlay_node_for(const nid_t& underlying_id) const {
        return this->has_node(underlying_id << 1);
    }
    
    handle_t StrandSplitGraph::get_overlay_handle(const handle_t& underlying_handle) const {
        // Get the ID of the node in the underlyign graph
        id_t underlying_id = graph->get_id(underlying_handle);
        assert(graph->has_node(underlying_id));
        // Get the ID of the named orientation of that node, in our graph.
        id_t overlay_id = (underlying_id << 1) | (graph->get_is_reverse(underlying_handle) ? 1 : 0);
        assert(has_node(overlay_id));
        // Get a handle to the forward orientation of our node
        return get_handle(overlay_id, false);
    }
}

