#pragma once

#include <string>
#include <vector>
#include <unordered_map>

#include <handlegraph/handle_graph.hpp>
#include "theseus/graph.h"

/**
 * @file handle_graph_theseus_adapter.hpp
 *
 * Builds a theseus::Graph from a HandleGraph and maintains a bidirectional
 * mapping between theseus integer vertex IDs and handle_t values.
 *
 * Why this exists instead of theseus::Graph(HandleGraph):
 *   - The upstream constructor has a naming collision bug: both orientations of
 *     a node get the same name because get_id() strips the orientation bit,
 *     causing silent overwrites in name_to_id_.
 *   - After alignment, theseus reports results as integer vertex IDs with no
 *     built-in path back to handle_t.
 *
 * Vertex naming convention: "<node_id>+" for forward, "<node_id>-" for reverse.
 * Use handle_name() to produce the start_node string required by align().
 *
 * Usage:
 *   HandleGraphTheseusAdapter adapter(my_handle_graph);
 *   TheseusAlignerImpl impl(penalties, adapter.take_graph(), false);
 *   auto aln = impl.align(seq, adapter.handle_name(start_handle));
 *   handle_t result_handle = adapter.vertex_to_handle(aln.some_vertex_id);
 */

namespace vg {

class HandleGraphTheseusAdapter {
public:
    explicit HandleGraphTheseusAdapter(const handlegraph::HandleGraph& hg) {
        const size_t n = hg.get_node_count();
        graph_._vertices.reserve(2 * n);
        graph_.name_to_id_.reserve(2 * n);
        id_to_handle_.reserve(2 * n);
        handle_to_id_.reserve(2 * n);

        // Add a vertex for each handle orientation and record the mapping.
        hg.for_each_handle([&](const handlegraph::handle_t& h) {
            add_vertex(hg, h);
            add_vertex(hg, hg.flip(h));
        });

        // for_each_edge visits each edge once in canonical form (first.id <=
        // second.id). Add the forward traversal and its reverse complement.
        hg.for_each_edge([&](const handlegraph::edge_t& e) {
            add_edge(e.first, e.second);
            // Reverse complement: flip both endpoints and swap direction.
            add_edge(hg.flip(e.second), hg.flip(e.first));
        });
    }

    // Move the built graph out for use with TheseusAlignerImpl.
    // The adapter retains the handle<->vertex mappings; only call this once.
    theseus::Graph&& take_graph() { return std::move(graph_); }

    // The vertex name string to pass as start_node to TheseusAlignerImpl::align().
    static std::string handle_name(const handlegraph::HandleGraph& hg,
                                   const handlegraph::handle_t& h) {
        return std::to_string(hg.get_id(h)) +
               (hg.get_is_reverse(h) ? '-' : '+');
    }

    // Convert a theseus vertex index (from alignment results) back to a handle_t.
    handlegraph::handle_t vertex_to_handle(int vertex_id) const {
        return id_to_handle_.at(static_cast<size_t>(vertex_id));
    }

    // Convert a handle_t to its theseus vertex index.
    int handle_to_vertex(const handlegraph::handle_t& h) const {
        return handle_to_id_.at(h);
    }

private:
    theseus::Graph graph_;
    std::vector<handlegraph::handle_t> id_to_handle_;
    std::unordered_map<handlegraph::handle_t, int> handle_to_id_;

    void add_vertex(const handlegraph::HandleGraph& hg,
                    const handlegraph::handle_t& h) {
        const int idx = static_cast<int>(graph_._vertices.size());
        theseus::Graph::vertex v;
        v.name  = handle_name(hg, h);
        v.value = hg.get_sequence(h);  // sequence in handle's orientation
        graph_._vertices.push_back(std::move(v));
        graph_.name_to_id_[graph_._vertices.back().name] = idx;
        id_to_handle_.push_back(h);
        handle_to_id_[h] = idx;
    }

    void add_edge(const handlegraph::handle_t& from,
                  const handlegraph::handle_t& to) {
        const int from_idx = handle_to_id_.at(from);
        const int to_idx   = handle_to_id_.at(to);
        theseus::Graph::edge e;
        e.from_vertex = from_idx;
        e.to_vertex   = to_idx;
        e.overlap     = 0;
        graph_._vertices[from_idx].out_edges.push_back(e);
        graph_._vertices[to_idx].in_edges.push_back(e);
    }
};

} // namespace vg
