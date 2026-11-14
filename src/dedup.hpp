#ifndef VG_DEDUP_HPP_INCLUDED
#define VG_DEDUP_HPP_INCLUDED

/**
 * \file dedup.hpp
 *
 * Content-aware node chopping to enable sequence deduplication in GBZ.
 *
 * Full unclipped graphs contain long unaligned segments such as centromeres
 * that are repetitive and hence compressible. GBZ used to store a naive
 * concatenation of the sequences in memory. From version 4, it can store an
 * arbitrary superstring that can be potentially much shorter.
 *
 * Here we use Content-Defined Chunking (CDC) to identify good positions for
 * chopping long nodes. We simultaneously collect a list of superstrings that
 * cover the sequence labels of all chopped pieces.
 */

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <handlegraph/handle_graph.hpp>

namespace vg {

//------------------------------------------------------------------------------

/// Settings for content-aware node chopping.
struct ContentChopParameters {
    /// Maximum allowed node length. Only nodes longer than this are chopped,
    /// and individual pieces will each be at most this length.
    size_t max_node_length = 1024;

    /// Relative cost (in bp) of each chopped piece. We try to minimize
    /// `L + piece_weight * N`, where `L` is the total length of stored strings
    /// and `N` is the number of chopped pieces.
    double piece_weight = 16.0;

    /// K-mer length for identifying candidate positions. At most 32.
    size_t kmer_length = 16;

    /// A position is a candidate for chopping, if the forward strand k-mer
    /// starting at it has a smaller hash than any k-mer starting within this
    /// radius. Serves as a strict lower bound for the distance between
    /// candidate positions. Similar to minimizers with window length
    /// `2 * anchor_radius + 1`.
    size_t anchor_radius = 32;

    /// Throws `std::invalid_argument` if the settings cannot be used.
    void validate() const;
};

/// Location of a piece (chopped node) in the original node and the stored strings.
struct PieceLocation {
    /// Original node identifier before chopping.
    handlegraph::nid_t node;
    /// Start offset of the piece in the original node.
    size_t node_offset;
    /// String that stores the sequence of the piece.
    size_t string_id;
    /// Start offset of the piece in the string.
    size_t string_offset;
    /// Length of the piece.
    size_t length;
};

/// A plan for content-aware chopping of long nodes. Only stores information for
/// nodes that were selected for chopping.
struct ChoppingPlan {
    /// List of chop points (start positions of pieces) in the original nodes.
    /// Serves as input for `handlegraph::algorithms::chop()`.
    std::vector<std::pair<handlegraph::nid_t, std::vector<size_t>>> chop_points;

    /// Stored strings that are collectively superstrings of the chopped pieces.
    std::vector<std::string> strings;

    /// Location of each chopped piece in the original nodes and stored strings,
    /// sorted by the original location.
    std::vector<PieceLocation> pieces;

    /// Returns the location of a piece starting at the given offset in the
    /// original node, or `nullptr` if no piece starts there.
    const PieceLocation* find_piece(handlegraph::nid_t node, size_t node_offset) const;

    /// Returns the sequence of the piece.
    std::string_view piece_sequence(const PieceLocation& piece) const;

    /// Returns the total length of the stored strings.
    size_t stored_length() const;
};

/**
 * Finds chopping points in the given graph according to the given parameters
 * for content-aware node chopping. Uses OpenMP threads. For best results,
 * long repetitive unitigs should be single nodes in the graph.
 *
 * @throws `std::invalid_argument` if the parameters cannot be used.
 *
 * The algorithm is roughly as follows:
 *
 * 1. Identify candidate positions in each long node by minimum k-mer hash
 *    within the given radius.
 * 2. Compute the frequency of each distinct segment over all long nodes.
 *    A segment is a substring between two candidate positions, at most
 *    the maximum segment length apart.
 * 3. Find the cheapest chopping for each long node using dynamic programming.
 * 4. Pack the sequences corresponding to the selected pieces into stored
 *    strings. We currently consider only substring relationships but not
 *    partial overlaps.
 */
ChoppingPlan plan_content_chopping(const handlegraph::HandleGraph& graph, const ContentChopParameters& parameters);

//------------------------------------------------------------------------------

}

#endif
