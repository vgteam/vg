#include "piece_directory.hpp"

#include "gfa.hpp"

#include <handlegraph/algorithms/chop.hpp>

#include <algorithm>
#include <ostream>
#include <stdexcept>

namespace vg {

using handlegraph::handle_t;
using handlegraph::nid_t;

//------------------------------------------------------------------------------

void PieceDirectory::add(nid_t node, const PieceLocation& piece) {
    this->pieces.insert_or_assign(node, piece);
}

const PieceLocation* PieceDirectory::find(nid_t node) const {
    auto iter = this->pieces.find(node);
    return (iter == this->pieces.end() ? nullptr : &iter->second);
}

size_t PieceDirectory::size() const {
    return this->pieces.size();
}

//------------------------------------------------------------------------------

PieceDirectory chop_by_plan(handlegraph::MutablePathDeletableHandleGraph& graph, const ChoppingPlan& plan) {
    PieceDirectory directory;
    // chop() reports each new node with the ID of its original node and the
    // forward offset where it starts in that node. It also reports unchopped
    // nodes that got a new ID, and the plan has no pieces for those.
    handlegraph::algorithms::chop(graph, plan.chop_points, [&](nid_t old_id, size_t offset, size_t, handle_t new_handle) {
        const PieceLocation* piece = plan.find_piece(old_id, offset);
        if (piece == nullptr) {
            return;
        }
        if (graph.get_length(new_handle) != piece->length) {
            throw std::runtime_error(
                "chop_by_plan(): piece at offset " + std::to_string(offset) + " of node " + std::to_string(old_id)
                + " has length " + std::to_string(graph.get_length(new_handle))
                + " in the chopped graph and length " + std::to_string(piece->length) + " in the plan"
            );
        }
        directory.add(graph.get_id(new_handle), *piece);
    });
    if (directory.size() != plan.pieces.size()) {
        throw std::runtime_error(
            "chop_by_plan(): the chopped graph has " + std::to_string(directory.size())
            + " of the " + std::to_string(plan.pieces.size()) + " pieces in the plan"
        );
    }
    return directory;
}

//------------------------------------------------------------------------------

std::string stored_string_name(const std::string& prefix, size_t string_id) {
    return prefix + std::to_string(string_id);
}

void write_stored_strings_fasta(const ChoppingPlan& plan, const std::string& name_prefix, std::ostream& out) {
    constexpr size_t LINE_WIDTH = 80;
    for (size_t string_id = 0; string_id < plan.strings.size(); string_id++) {
        out << ">" << stored_string_name(name_prefix, string_id) << "\n";
        const std::string& sequence = plan.strings[string_id];
        for (size_t offset = 0; offset < sequence.size(); offset += LINE_WIDTH) {
            size_t length = std::min(LINE_WIDTH, sequence.size() - offset);
            out.write(sequence.data() + offset, length);
            out << "\n";
        }
    }
}

namespace {

/// Write a GFA segment line for each node of the graph. A piece in the
/// directory gets tags that give the location of its sequence in a stored
/// string.
void write_gfa_piece_nodes(
    const handlegraph::HandleGraph& graph, const PieceDirectory& directory,
    const std::string& name_prefix, std::ostream& out
) {
    graph.for_each_handle([&](const handle_t& handle) {
        nid_t node_id = graph.get_id(handle);
        out << "S\t" << node_id << "\t";
        const PieceLocation* piece = directory.find(node_id);
        if (piece == nullptr) {
            out << graph.get_sequence(handle);
        } else {
            out << "*"
                << "\tsn:Z:" << stored_string_name(name_prefix, piece->string_id)
                << "\tso:i:" << piece->string_offset
                << "\tLN:i:" << piece->length;
        }
        out << "\n";
    });
}

}

void write_chopped_gfa(
    const handlegraph::PathHandleGraph& graph, const PieceDirectory& directory,
    const std::string& name_prefix, std::ostream& out
) {
    write_gfa_header(&graph, out);
    write_gfa_piece_nodes(graph, directory, name_prefix, out);
    write_gfa_edges(&graph, out);
    write_gfa_paths(&graph, out);
}

//------------------------------------------------------------------------------

}
