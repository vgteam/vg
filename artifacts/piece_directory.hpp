#ifndef VG_PIECE_DIRECTORY_HPP_INCLUDED
#define VG_PIECE_DIRECTORY_HPP_INCLUDED

/**
 * \file piece_directory.hpp
 *
 * Chop a graph by a content-aware chopping plan (see dedup.hpp), and write the
 * plan's stored strings and the chopped graph. The GFA writer gives the
 * sequence of each piece as a location in a stored string.
 */

#include "dedup.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <unordered_map>

#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

namespace vg {

//------------------------------------------------------------------------------

/**
 * A directory that says where the stored strings of a chopping plan hold the
 * sequence of each piece in a chopped graph.
 *
 * It lists the pieces by their node IDs in the chopped graph. It lists only
 * pieces, so unchopped nodes are missing from it.
 *
 * The `node` and `node_offset` fields of each location still give the place of
 * the piece in the original node, before chopping.
 */
class PieceDirectory {
public:
    /// Record that the node with the given ID in the chopped graph is the given
    /// piece. Replace any earlier location for the node.
    void add(handlegraph::nid_t node, const PieceLocation& piece);

    /// Return where the sequence of the node with the given ID is stored, or
    /// `nullptr` if the node is not a piece.
    const PieceLocation* find(handlegraph::nid_t node) const;

    /// Return the number of pieces in the directory.
    size_t size() const;

private:
    std::unordered_map<handlegraph::nid_t, PieceLocation> pieces;
};

/**
 * Chop the graph at the chop points of the plan.
 *
 * Chopping may give new IDs to all nodes, and it invalidates handles into the
 * graph.
 *
 * @pre The plan was made for this graph.
 *
 * @return A directory of the nodes in the chopped graph that are pieces.
 *
 * @throws `std::runtime_error` if the chopped graph does not have each piece of
 * the plan as a node of the right length. The graph is left chopped.
 */
PieceDirectory chop_by_plan(handlegraph::MutablePathDeletableHandleGraph& graph, const ChoppingPlan& plan);

//------------------------------------------------------------------------------

/// Return the name of stored string `string_id`. The name is the prefix
/// followed by the number of the string.
std::string stored_string_name(const std::string& prefix, size_t string_id);

/// Write the stored strings of the plan in FASTA format, named by
/// `stored_string_name()`.
void write_stored_strings_fasta(const ChoppingPlan& plan, const std::string& name_prefix, std::ostream& out);

/**
 * Write the chopped graph in GFA format, with its edges and paths.
 *
 * A node that is a piece in the directory has `*` as its sequence. Tag `sn:Z`
 * on its segment line names the stored string that holds the sequence, by
 * `stored_string_name()`. Tag `so:i` gives the offset of the sequence in the
 * string, and tag `LN:i` its length. Other nodes have their sequences on their
 * segment lines.
 */
void write_chopped_gfa(
    const handlegraph::PathHandleGraph& graph, const PieceDirectory& directory,
    const std::string& name_prefix, std::ostream& out
);

//------------------------------------------------------------------------------

}

#endif
