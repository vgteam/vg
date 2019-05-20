/**
 * \file proto_handle_graph.hpp: implementation of the ProtoHandleGraph
 */


#include "proto_handle_graph.hpp"

#include <handlegraph/util.hpp>
#include <atomic>


namespace vg {

using namespace std;
using namespace handlegraph;

    ProtoHandleGraph::ProtoHandleGraph(const Graph* graph) : graph(graph) {
        // nothing to do
    }
    
    handle_t ProtoHandleGraph::get_handle_by_index(const size_t& i) const {
        return handlegraph::number_bool_packing::pack(i, false);
    }
    
    size_t ProtoHandleGraph::edge_size() const {
        return graph->edge_size();
    }
    
    edge_t ProtoHandleGraph::get_edge_by_index(const size_t& i) const {
        const Edge& edge = graph->edge(i);
        return edge_handle(get_handle(edge.from(), edge.from_start()),
                           get_handle(edge.to(), edge.to_end()));
    }
    
    bool ProtoHandleGraph::has_node(id_t node_id) const {
        bool found = false;
        for (size_t i = 0; i < graph->node_size() && !found; i++) {
            found = graph->node(i).id() == node_id;
        }
        return found;
    }
    
    handle_t ProtoHandleGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        for (size_t i = 0; i < graph->node_size(); i++) {
            if (graph->node(i).id() == node_id) {
                return handlegraph::number_bool_packing::pack(i, is_reverse);
            }
        }
        // tried to find a handle for a node that doesn't exist
        cerr << "error::[ProtoHandleGraph] requested handle for a node that does not exist: " << node_id << endl;
        return as_handle(-1);
    }
    
    id_t ProtoHandleGraph::get_id(const handle_t& handle) const {
        return graph->node(handlegraph::number_bool_packing::unpack_number(handle)).id();
    }
    
    bool ProtoHandleGraph::get_is_reverse(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_bit(handle);
    }
    
    handle_t ProtoHandleGraph::flip(const handle_t& handle) const {
        return handlegraph::number_bool_packing::toggle_bit(handle);
    }
    
    size_t ProtoHandleGraph::get_length(const handle_t& handle) const {
        return graph->node(handlegraph::number_bool_packing::unpack_number(handle)).sequence().size();
    }
    
    string ProtoHandleGraph::get_sequence(const handle_t& handle) const {
        return graph->node(handlegraph::number_bool_packing::unpack_number(handle)).sequence();
    }
    
    bool ProtoHandleGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                             const function<bool(const handle_t&)>& iteratee) const {
        bool keep_going = true;
        bool leftward = (go_left != get_is_reverse(handle));
        for (size_t i = 0; i < graph->edge_size() && keep_going; i++) {
            const Edge& edge = graph->edge(i);
            if (edge.from() == get_id(handle) && leftward == edge.from_start()) {
                keep_going = iteratee(get_handle(edge.to(), go_left != edge.to_end()));
            }
            else if (edge.to() == get_id(handle) && leftward != edge.to_end()) {
                keep_going = iteratee(get_handle(edge.from(), go_left == edge.from_start()));
            }
        }
        return keep_going;
    }
    
    bool ProtoHandleGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel) const {
        if (parallel) {
            size_t num_nodes = graph->node_size();
            atomic<bool> keep_going(true);
#pragma omp parallel shared(num_nodes)
            for (size_t i = 0; i < num_nodes; i++) {
                if (keep_going && !iteratee(handlegraph::number_bool_packing::pack(i, false))) {
                    keep_going = false;
                }
            }
            return keep_going;
        }
        else {
            for (size_t i = 0; i < graph->node_size(); i++) {
                if (!iteratee(handlegraph::number_bool_packing::pack(i, false))) {
                    return false;
                }
            }
            return true;
        }
    }
    
    size_t ProtoHandleGraph::get_node_count() const {
        return graph->node_size();
    }
    
    id_t ProtoHandleGraph::min_node_id() const {
        id_t min_id = numeric_limits<id_t>::max();
        for (size_t i = 0; i < graph->node_size(); i++) {
            min_id = min(min_id, graph->node(i).id());
        }
        return min_id;
    }
    
    id_t ProtoHandleGraph::max_node_id() const {
        id_t max_id = numeric_limits<id_t>::min();
        for (size_t i = 0; i < graph->node_size(); i++) {
            max_id = max(max_id, graph->node(i).id());
        }
        return max_id;
    }
    
    bool ProtoHandleGraph::for_each_edge(const function<bool(const edge_t&)>& iteratee,
                                         bool parallel) const {
        
        // we need to be able to map IDs to indexes to make handles efficiently
        unordered_map<id_t, size_t> node_idx;
        node_idx.reserve(graph->node_size());
        for (size_t i = 0; i < graph->node_size(); i++) {
            node_idx[graph->node(i).id()] = i;
        }
        
        if (parallel) {
            atomic<bool> keep_going(true);
#pragma omp parallel for
            for (size_t i = 0; i < graph->edge_size(); i++) {
                const Edge& edge = graph->edge(i);
                if (keep_going && !iteratee(edge_t(handlegraph::number_bool_packing::pack(node_idx[edge.from()], edge.from_start()),
                                            handlegraph::number_bool_packing::pack(node_idx[edge.to()], edge.to_end())))) {
                
                    keep_going = false;
                }
            }
            return keep_going;
        }
        else {
            bool keep_going = true;
            for (size_t i = 0; i < graph->edge_size() && keep_going; i++) {
                const Edge& edge = graph->edge(i);
                keep_going = iteratee(edge_t(handlegraph::number_bool_packing::pack(node_idx[edge.from()], edge.from_start()),
                                             handlegraph::number_bool_packing::pack(node_idx[edge.to()], edge.to_end())));
            }
            return keep_going;
        }
    }
}

