#include "graph.hpp"

namespace vg {

void sort_by_id_dedup_and_clean(Graph& graph) {
    remove_duplicates(graph); // graph is sorted here
    remove_orphan_edges(graph);
}

void remove_duplicates(Graph& graph) {
    remove_duplicate_nodes(graph);
    remove_duplicate_edges(graph);
}

void remove_duplicate_edges(Graph& graph) {
    sort_edges_by_id(graph);
    graph.mutable_edge()->erase(std::unique(graph.mutable_edge()->begin(),
                                            graph.mutable_edge()->end(),
                                            [](const Edge& a, const Edge& b) {
                                                return make_tuple(a.from(), a.to(), a.from_start(), a.to_end())
                                                    == make_tuple(b.from(), b.to(), b.from_start(), b.to_end());
                                            }), graph.mutable_edge()->end());

}

void remove_duplicate_nodes(Graph& graph) {
    sort_nodes_by_id(graph);
    graph.mutable_node()->erase(std::unique(graph.mutable_node()->begin(),
                                            graph.mutable_node()->end(),
                                            [](const Node& a, const Node& b) {
                                                return a.id() == b.id();
                                            }), graph.mutable_node()->end());
}

void remove_orphan_edges(Graph& graph) {
    set<id_t> ids;
    for (auto& node : graph.node()) {
        ids.insert(node.id());
    }
    graph.mutable_edge()->erase(std::remove_if(graph.mutable_edge()->begin(),
                                               graph.mutable_edge()->end(),
                                               [&ids](const Edge& e) {
                                                   return !ids.count(e.from()) || !ids.count(e.to());
                                               }), graph.mutable_edge()->end());
}

void sort_by_id(Graph& graph) {
    sort_nodes_by_id(graph);
    sort_edges_by_id(graph);
}

void sort_nodes_by_id(Graph& graph) {
    std::sort(graph.mutable_node()->begin(),
              graph.mutable_node()->end(),
              [](const Node& a, const Node& b) {
                  return a.id() < b.id();
              });
}

void sort_edges_by_id(Graph& graph) {
    std::sort(graph.mutable_edge()->begin(),
              graph.mutable_edge()->end(),
              [](const Edge& a, const Edge& b) {
                  return make_tuple(a.from(), a.to(), a.from_start(), a.to_end())
                      < make_tuple(b.from(), b.to(), b.from_start(), b.to_end());
              });
}

bool is_id_sortable(const Graph& graph) {
    for (auto& edge : graph.edge()) {
        if (edge.from() >= edge.to()) return false;
    }
    return true;
}

bool has_inversion(const Graph& graph) {
    for (auto& edge : graph.edge()) {
        if (edge.from_start() || edge.to_end()) return true;
    }
    return false;
}

void flip_doubly_reversed_edges(Graph& graph) {
    for (auto& edge : *graph.mutable_edge()) {
        if (edge.from_start() && edge.to_end()) {
            edge.set_from_start(false);
            edge.set_to_end(false);
        }
    }
}
    
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
