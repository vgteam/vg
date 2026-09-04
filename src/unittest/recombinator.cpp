/** \file
 *
 * Unit tests for recombinator.cpp, specifically for clipping gref fragments to a
 * sampled graph.
 */

#include "../recombinator.hpp"

#include "catch.hpp"

#include <set>
#include <utility>
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

// Paths 1-2-3, 3-4, 5-6 and 8. Every node from 1 to 8 except 7 is in the index,
// with 7 inside the alphabet range but unused. Nodes 4 and 5 are both present
// but there is no edge between them.
std::vector<gbwt::vector_type> index_paths() {
    return {
        { fwd(1), fwd(2), fwd(3) },
        { fwd(3), fwd(4) },
        { fwd(5), fwd(6) },
        { fwd(8) },
    };
}

void insert_paths(gbwt::GBWTBuilder& builder, const std::vector<gbwt::vector_type>& paths) {
    for (const gbwt::vector_type& path : paths) {
        builder.insert(path, true);
    }
}

// All (node, successor) pairs in the index, excluding the endmarker.
std::set<std::pair<gbwt::node_type, gbwt::node_type>> edges_in(const gbwt::DynamicGBWT& index) {
    std::set<std::pair<gbwt::node_type, gbwt::node_type>> result;
    for (gbwt::node_type node = index.firstNode(); node < index.sigma(); node++) {
        if (index.empty(node)) { continue; }
        for (gbwt::edge_type edge : index.edges(node)) {
            if (edge.first != gbwt::ENDMARKER) {
                result.emplace(node, edge.first);
            }
        }
    }
    return result;
}

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
    gbwt::GBWTBuilder builder(gbwt::bit_length(rev(8)), 1024);
    insert_paths(builder, index_paths());
    builder.finish();
    const gbwt::DynamicGBWT& index = builder.index;

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

    SECTION("a node inside the alphabet range but not in the index is missing") {
        REQUIRE(index.contains(fwd(7)));
        gbwt::vector_type path { fwd(7) };
        REQUIRE(clip_path_to_index(path, index, id_as_length).empty());
    }

    SECTION("a node outside the alphabet range is missing") {
        gbwt::vector_type path { fwd(2), fwd(100) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 1, 0, 2);
    }

    SECTION("an empty path has no pieces") {
        gbwt::vector_type path;
        REQUIRE(clip_path_to_index(path, index, id_as_length).empty());
    }

    SECTION("nothing survives in an empty index") {
        gbwt::DynamicGBWT empty;
        gbwt::vector_type path { fwd(1), fwd(2) };
        REQUIRE(clip_path_to_index(path, empty, id_as_length).empty());
    }
}

TEST_CASE("A GBWT builder can be finished, extended, and finished again", "[gref][recombinator]") {
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder builder(gbwt::bit_length(rev(8)), 1024);
    std::vector<gbwt::vector_type> paths = index_paths();
    insert_paths(builder, paths);
    builder.finish();
    size_t sequences = builder.index.sequences();
    auto edges = edges_in(builder.index);

    SECTION("a path over existing nodes and edges does not change the topology") {
        gbwt::vector_type path { fwd(1), fwd(2), fwd(3), fwd(4) };
        builder.insert(path, true);
        builder.finish();
        REQUIRE(builder.index.sequences() == sequences + 2);
        REQUIRE(builder.index.bidirectional());
        REQUIRE(edges_in(builder.index) == edges);

        // Everything is still there, and the index can still be compressed.
        paths.push_back(path);
        gbwt::GBWT compressed(builder.index);
        REQUIRE(compressed.sequences() == sequences + 2);
        for (size_t i = 0; i < paths.size(); i++) {
            REQUIRE(compressed.extract(gbwt::Path::encode(i, false)) == paths[i]);
        }
    }

    SECTION("a path with a new edge is detected by the test itself") {
        gbwt::vector_type path { fwd(1), fwd(3) };
        builder.insert(path, true);
        builder.finish();
        auto new_edges = edges_in(builder.index);
        REQUIRE(new_edges.size() == edges.size() + 2);
        REQUIRE(new_edges.count(std::make_pair(fwd(1), fwd(3))) == 1);
        REQUIRE(new_edges.count(std::make_pair(rev(3), rev(1))) == 1);
    }
}

TEST_CASE("Gref fragments are recognized from sample and contig names", "[gref][recombinator]") {
    SECTION("PanSN names") {
        REQUIRE(is_gref_fragment("gref_CHM13", "chr1_7_alt"));
        REQUIRE_FALSE(is_gref_fragment("gref_CHM13", "chr1")); // A gref copy of the reference path.
        REQUIRE_FALSE(is_gref_fragment("GRCh38", "chr1_7_alt")); // A contig that happens to look like a fragment.
        REQUIRE_FALSE(is_gref_fragment("CHM13", "chr1"));
    }

    SECTION("generic paths carry their full name as the contig") {
        REQUIRE(is_gref_fragment(gbwtgraph::GENERIC_PATH_SAMPLE_NAME, "gref_x_1_alt"));
        REQUIRE(is_gref_fragment(PathMetadata::NO_SAMPLE_NAME, "gref_x_1_alt"));
        REQUIRE_FALSE(is_gref_fragment(gbwtgraph::GENERIC_PATH_SAMPLE_NAME, "gref_x"));
        REQUIRE_FALSE(is_gref_fragment(gbwtgraph::GENERIC_PATH_SAMPLE_NAME, "x_1_alt"));
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg
