#ifndef VG_DEDUP_HPP_INCLUDED
#define VG_DEDUP_HPP_INCLUDED

/**
 * \file dedup.hpp
 * Content-aware chopping of long nodes, so that a graph's node sequences can be stored once.
 *
 * A graph that keeps the sequence of every node in one array stores repeated sequence once per
 * copy. If we cut the long nodes at places that depend only on the nearby sequence, two copies of
 * a repeat are cut at the same places, and most of their pieces have identical sequences. Each
 * piece can then point into shared stored text instead of owning its sequence.
 *
 * The planner works in four steps:
 *
 * 1. In each long node, it picks the candidate cuts: the offsets where a k-mer is the strict
 *    minimum, by hash, of the k-mers near it. These offsets depend only on the nearby sequence.
 * 2. It counts how many times the sequence of each candidate segment occurs, over all long
 *    nodes. A candidate segment is the sequence between two candidate cuts at most the maximum
 *    node length apart.
 * 3. In each long node, it chooses the cuts by dynamic programming. A piece of length `L` whose
 *    sequence occurs `c` times costs `piece_weight + L / c`. The `c` copies share the stored
 *    sequence, so each pays `1 / c` of its length. Pieces of repeated sequence are cheap, so the
 *    chosen pieces are long. A variant that one copy of a repeat does not share goes into a short
 *    piece of its own.
 * 4. It packs the sequences of the pieces of the chopped nodes into a list of strings. It places
 *    the longest pieces first. A piece whose sequence is already stored, as another piece or
 *    within a run of stored pieces, points there. Otherwise its sequence is stored.
 *
 * Nodes that are not chopped are left out of the plan. Their sequences can be stored one after
 * another in node order, as before, which keeps most of the genome in a compact layout.
 *
 * Only forward strands are compared. A piece whose reverse complement is stored gets stored again.
 */

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <handlegraph/handle_graph.hpp>

namespace vg {

//------------------------------------------------------------------------------

/// The settings for content-aware chopping.
struct ContentChopParameters {
    /// The maximum length of a node in the chopped graph. Longer nodes are chopped into pieces of
    /// at most this length, and other nodes are left whole.
    size_t max_node_length = 1024;

    /// The cost of each piece, in bases of stored sequence. A larger weight gives fewer and
    /// longer pieces, and stores more sequence.
    double piece_weight = 16.0;

    /// The length of the k-mers that decide the candidate cuts. It may be at most 32.
    size_t kmer_length = 16;

    /// A k-mer is a candidate cut when its hash is smaller than the hash of every other k-mer
    /// starting at most this many bases away. Candidate cuts are then more than this many bases
    /// apart, and about twice this many on average.
    size_t anchor_radius = 32;

    /// Throw std::invalid_argument if the settings cannot be used.
    void validate() const;
};

//------------------------------------------------------------------------------

/// The place of one piece of a node in the graph and in the stored strings.
struct PieceLocation {
    /// The ID of the node, before chopping, that the piece is part of.
    handlegraph::nid_t node;
    /// The offset of the piece in the forward sequence of that node.
    size_t node_offset;
    /// The string that holds the sequence of the piece.
    size_t string_id;
    /// The offset of the sequence of the piece in that string.
    size_t string_offset;
    /// The length of the piece.
    size_t length;
};

/**
 * A ChoppingPlan is a plan for chopping some nodes of a graph and for storing the sequence of
 * each of their pieces once. Nodes that are not chopped are not in the plan.
 *
 * Every piece of every chopped node has a location, and the sequence of each piece equals the
 * substring of the string at that location. Pieces are never empty.
 */
struct ChoppingPlan {
    /// The places to cut each chopped node, as offsets in its forward sequence, in the form that
    /// `handlegraph::algorithms::chop()` takes. Nodes that are not chopped are not listed.
    std::vector<std::pair<handlegraph::nid_t, std::vector<size_t>>> chop_points;

    /// The stored strings. Each one is a substring of the sequence of a chopped node.
    std::vector<std::string> strings;

    /// The location of every piece of every chopped node, in order of node ID and offset in the
    /// node.
    std::vector<PieceLocation> pieces;

    /// Return the location of the piece that starts at this offset in this node, before
    /// chopping, or nullptr if no piece of a chopped node starts there. After `chop()`, its
    /// `record_change` callback gives the node and offset for each new node it reports. Nodes
    /// that it does not report keep their ID and start at offset 0. A new node with no piece in
    /// the plan is a node that was not chopped.
    const PieceLocation* find_piece(handlegraph::nid_t node, size_t node_offset) const;

    /// Return the sequence of the piece.
    std::string_view piece_sequence(const PieceLocation& piece) const;

    /// Return the total length of the stored strings.
    size_t stored_length() const;
};

/**
 * Plan content-aware chopping of the graph, and plan where to store the sequences of the
 * pieces. Chop each node longer than the maximum node length into pieces of at most that
 * length, and leave the other nodes whole.
 *
 * The graph should have its unitigs merged into single nodes, because only the inside of a node
 * is chopped.
 *
 * Finds the cuts on several OpenMP threads, and packs the sequences on one thread.
 *
 * @throws std::invalid_argument if the parameters cannot be used.
 */
ChoppingPlan plan_content_chopping(const handlegraph::HandleGraph& graph,
                                   const ContentChopParameters& parameters);

//------------------------------------------------------------------------------

}

#endif
