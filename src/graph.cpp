#include "graph.hpp"

namespace vg {

void from_handle_graph(const HandleGraph& from, Graph& to) {
    from.for_each_handle([&](const handle_t& h) {
        Node* node = to.add_node();
        node->set_id(from.get_id(h));
        node->set_sequence(from.get_sequence(h));
    });
    from.for_each_edge([&](const edge_t& e) {
        Edge* edge = to.add_edge();
        edge->set_from(from.get_id(e.first));
        edge->set_from_start(from.get_is_reverse(e.first));
        edge->set_to(from.get_id(e.second));
        edge->set_to_end(from.get_is_reverse(e.second));
    });
}

void from_path_handle_graph(const PathHandleGraph& from, Graph& to) {
    
    from_handle_graph(from, to);
    
    from.for_each_path_handle([&](const path_handle_t& p) {
        Path* path = to.add_path();
        path->set_name(from.get_path_name(p));
        path->set_is_circular(from.get_is_circular(p));
        int64_t rank = 1;
        for (handle_t step : from.scan_path(p)) {
            Mapping* mapping = path->add_mapping();
            Position* position = mapping->mutable_position();
            position->set_node_id(from.get_id(step));
            position->set_is_reverse(from.get_is_reverse(step));
            mapping->set_rank(rank);
            ++rank;
        }
    });
}

}
