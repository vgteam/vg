#pragma once

#include <string>
#include <vector>
#include <unordered_map>

#include <handlegraph/handle_graph.hpp>
#include <theseus/graph.h>
#include <theseus/alignment.h>
#include <vg/vg.pb.h>

/**
 * @file theseus_interop.hpp
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

/**
 * @brief Convert a theseus::Alignment to a vg::Alignment.
 *
 * theseus edit operations are per-base characters:
 *   'M' = match, 'X' = mismatch/SNP, 'I' = insertion, 'D' = deletion.
 * Consecutive operations of the same type are merged into a single vg Edit.
 * Edits are split at node boundaries according to the reference bases consumed
 * in each vertex of the theseus path.
 *
 * @param theseus_aln   Alignment returned by TheseusAlignerImpl::align().
 * @param sequence      The read sequence that was aligned.
 * @param adapter       Adapter used to map theseus vertex IDs → handle_t.
 * @param graph         The HandleGraph the alignment was made against.
 * @return vg::Alignment with sequence, path, and per-node mappings filled in.
 */
inline vg::Alignment theseus_to_vg_alignment(
        const theseus::Alignment&          theseus_aln,
        const std::string&                 sequence,
        const HandleGraphTheseusAdapter&   adapter,
        const handlegraph::HandleGraph&    graph) {

    vg::Alignment vg_aln;
    vg_aln.set_sequence(sequence);

    if (theseus_aln.path.empty() || theseus_aln.edit_op.empty()) {
        return vg_aln;
    }

    vg::Path* vg_path = vg_aln.mutable_path();

    const size_t n_edits = theseus_aln.edit_op.size();
    const size_t n_nodes = theseus_aln.path.size();
    size_t edit_cursor = 0; // index into edit_op
    size_t seq_cursor  = 0; // index into sequence (read bases consumed)

    for (size_t node_idx = 0; node_idx < n_nodes && edit_cursor < n_edits; ++node_idx) {
        const handlegraph::handle_t h = adapter.vertex_to_handle(theseus_aln.path[node_idx]);
        const size_t node_len = graph.get_length(h);

        // Reference range covered by this alignment on this node.
        // start_offset applies only to the first node;
        // end_offset (exclusive) applies only to the last node.
        const size_t ref_start = (node_idx == 0)
                                 ? static_cast<size_t>(theseus_aln.start_offset)
                                 : 0;
        const size_t ref_end   = (node_idx == n_nodes - 1)
                                 ? static_cast<size_t>(theseus_aln.end_offset)
                                 : node_len;
        size_t ref_remaining   = ref_end - ref_start;

        vg::Mapping* mapping = vg_path->add_mapping();
        mapping->mutable_position()->set_node_id(graph.get_id(h));
        mapping->mutable_position()->set_is_reverse(graph.get_is_reverse(h));
        mapping->mutable_position()->set_offset(static_cast<int64_t>(ref_start));
        mapping->set_rank(static_cast<int64_t>(node_idx + 1));

        while (edit_cursor < n_edits) {
            const char op = theseus_aln.edit_op[edit_cursor];

            // Insertions do not consume reference; attach to the current mapping.
            if (op == 'I') {
                const size_t seq_start = seq_cursor;
                size_t run = 0;
                while (edit_cursor < n_edits && theseus_aln.edit_op[edit_cursor] == 'I') {
                    ++run; ++edit_cursor; ++seq_cursor;
                }
                vg::Edit* e = mapping->add_edit();
                e->set_from_length(0);
                e->set_to_length(static_cast<int32_t>(run));
                e->set_sequence(sequence.substr(seq_start, run));
                continue;
            }

            // Reference-consuming op: stop if this node is exhausted.
            if (ref_remaining == 0) break;

            if (op == 'M') {
                size_t run = 0;
                while (edit_cursor < n_edits
                       && theseus_aln.edit_op[edit_cursor] == 'M'
                       && run < ref_remaining) {
                    ++run; ++edit_cursor; ++seq_cursor;
                }
                vg::Edit* e = mapping->add_edit();
                e->set_from_length(static_cast<int32_t>(run));
                e->set_to_length(static_cast<int32_t>(run));
                // no sequence field for matches
                ref_remaining -= run;
            } else if (op == 'X') {
                const size_t seq_start = seq_cursor;
                size_t run = 0;
                while (edit_cursor < n_edits
                       && theseus_aln.edit_op[edit_cursor] == 'X'
                       && run < ref_remaining) {
                    ++run; ++edit_cursor; ++seq_cursor;
                }
                vg::Edit* e = mapping->add_edit();
                e->set_from_length(static_cast<int32_t>(run));
                e->set_to_length(static_cast<int32_t>(run));
                e->set_sequence(sequence.substr(seq_start, run));
                ref_remaining -= run;
            } else if (op == 'D') {
                size_t run = 0;
                while (edit_cursor < n_edits
                       && theseus_aln.edit_op[edit_cursor] == 'D'
                       && run < ref_remaining) {
                    ++run; ++edit_cursor;
                    // deletions consume no read bases
                }
                vg::Edit* e = mapping->add_edit();
                e->set_from_length(static_cast<int32_t>(run));
                e->set_to_length(0);
                ref_remaining -= run;
            }
        }
    }

    return vg_aln;
}

} // namespace vg
