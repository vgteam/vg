/** \file
 *
 * Unit tests for recombinator.cpp, specifically for clipping gref fragments to a
 * sampled graph.
 */

#include "../recombinator.hpp"

#include "catch.hpp"

#include <vector>

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

// gbwt::vector_type may store nodes in a narrower type than gbwt::node_type.
typedef gbwt::vector_type::value_type node_t;
node_t fwd(nid_t id) { return gbwt::Node::encode(id, false); }
node_t rev(nid_t id) { return gbwt::Node::encode(id, true); }

// Node length is the node id, so offsets are easy to predict.
size_t id_as_length(gbwt::node_type node) { return gbwt::Node::id(node); }

void check_piece(const PathPiece& piece, size_t start, size_t nodes, size_t bp_offset, size_t bp_length) {
    REQUIRE(piece.start == start);
    REQUIRE(piece.nodes == nodes);
    REQUIRE(piece.bp_offset == bp_offset);
    REQUIRE(piece.bp_length == bp_length);
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Paths are clipped to the nodes and edges of a GBWT index", "[gref][recombinator]") {
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    // Paths 1-2-3, 3-4, 5-6 and 8. Node 7 is inside the alphabet range but unused,
    // and nodes 4 and 5 are both present with no edge between them.
    gbwt::GBWTBuilder builder(gbwt::bit_length(rev(8)), 1024);
    std::vector<gbwt::vector_type> paths {
        { fwd(1), fwd(2), fwd(3) }, { fwd(3), fwd(4) }, { fwd(5), fwd(6) }, { fwd(8) }
    };
    for (const gbwt::vector_type& path : paths) {
        builder.insert(path, true);
    }
    builder.finish();
    const gbwt::DynamicGBWT& index = builder.index;
    REQUIRE(index.contains(fwd(7)));

    SECTION("a path that is entirely in the index is one piece") {
        gbwt::vector_type path { fwd(1), fwd(2), fwd(3) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 3, 0, 6);
    }

    SECTION("the reverse orientation is in the index too") {
        gbwt::vector_type path { rev(3), rev(2), rev(1) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 3, 0, 6);
    }

    SECTION("a missing node splits the path and still counts toward the offset") {
        gbwt::vector_type path { fwd(1), fwd(2), fwd(7), fwd(3), fwd(4) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 2);
        check_piece(pieces[0], 0, 2, 0, 3);
        check_piece(pieces[1], 3, 2, 10, 7);
    }

    SECTION("present nodes without an edge between them are split") {
        gbwt::vector_type path { fwd(3), fwd(4), fwd(5), fwd(6) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 2);
        check_piece(pieces[0], 0, 2, 0, 7);
        check_piece(pieces[1], 2, 2, 7, 11);
    }

    SECTION("a node outside the alphabet range is missing") {
        gbwt::vector_type path { fwd(2), fwd(100) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 1, 0, 2);
    }

    SECTION("nothing survives in an empty index") {
        gbwt::DynamicGBWT empty;
        gbwt::vector_type path { fwd(1), fwd(2) };
        REQUIRE(clip_path_to_index(path, empty, id_as_length).empty());
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg
