/** \file
 *
 * Unit tests for SubchainWalk in subchain_walk.cpp, the walk of a GBWT sequence through a
 * subchain. Each walk over a graph of seven nodes is checked against the nodes and the
 * sequences written out in the test.
 */

#include "../../src/recombinator.hpp"
#include "../../src/gbwt_helper.hpp"
#include "../../src/utility.hpp"

#include "catch.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace vg::unittest {

namespace {

using node_t = gbwt::vector_type::value_type;
node_t fwd(nid_t id) { return static_cast<node_t>(gbwt::Node::encode(id, false)); }
node_t rev(nid_t id) { return static_cast<node_t>(gbwt::Node::encode(id, true)); }

// A visit of a walk, and whether it is the first visit of a piece.
using Visit = std::pair<gbwt::node_type, bool>;

// The number of bases a walk keeps of its first node and of the end node. It is k - 1 for
// kmers of 4 bases.
constexpr size_t FLANK = 3;

// A limit that keeps every base of the first node and of the end node.
constexpr size_t ALL = std::numeric_limits<size_t>::max();

// One GBWT path: a fragment of the haplotype of a sample, as the IDs of the nodes it
// visits forward.
struct PathSpec {
    size_t sample;
    size_t fragment;
    std::vector<nid_t> nodes;
};

// The samples, each with one haplotype of contig "chr".
enum Sample : size_t { BROKEN, WHOLE, END_FIRST, GAP, JOINED, LOOP };

// A GBZ with seven nodes and ten paths, and the fragment map of its haplotypes. A path ID
// is the position of the path in `paths`. The haplotype of BROKEN is cut between nodes 3
// and 4. The second fragment of END_FIRST starts with node 6 and has the lower path ID.
// GAP has an empty fragment. JOINED visits the nodes of BROKEN in one fragment, and LOOP
// visits nodes 3 and 4 twice.
struct WalkFixture {
    std::vector<PathSpec> paths {
        { BROKEN, 0, { 1, 2, 3 } },                 // path 0
        { BROKEN, 1, { 4, 6, 7 } },                 // path 1
        { WHOLE, 0, { 1, 2, 5, 6, 7 } },            // path 2
        { END_FIRST, 1, { 6, 7 } },                 // path 3
        { END_FIRST, 0, { 1, 2, 5 } },              // path 4
        { GAP, 0, { 1, 2, 3 } },                    // path 5
        { GAP, 1, { } },                            // path 6
        { GAP, 2, { 4, 6, 7 } },                    // path 7
        { JOINED, 0, { 1, 2, 3, 4, 6, 7 } },        // path 8
        { LOOP, 0, { 1, 2, 3, 4, 3, 4, 6, 7 } }     // path 9
    };
    std::map<nid_t, std::string> sequences {
        { 1, "GATTACA" }, { 2, "CCGGTA" }, { 3, "TTTAAAC" }, { 4, "ACGTACG" },
        { 5, "GGGCCCA" }, { 6, "CATCATG" }, { 7, "TGCATGC" }
    };
    gbwtgraph::GBZ gbz;
    gbwt::FragmentMap fragment_map;

    WalkFixture() {
        std::vector<gbwt::vector_type> gbwt_paths;
        for (const PathSpec& path : this->paths) {
            gbwt::vector_type nodes;
            for (nid_t id : path.nodes) {
                nodes.push_back(fwd(id));
            }
            gbwt_paths.push_back(nodes);
        }
        gbwt::GBWT index = get_gbwt(gbwt_paths);
        index.addMetadata();
        index.metadata.setSamples({ "broken", "whole", "end_first", "gap", "joined", "loop" });
        index.metadata.setHaplotypes(LOOP + 1);
        index.metadata.setContigs({ "chr" });
        for (const PathSpec& path : this->paths) {
            index.metadata.addPath(path.sample, 0, 0, path.fragment);
        }

        gbwtgraph::NaiveGraph source;
        for (const auto& [id, sequence] : this->sequences) {
            source.create_node(id, sequence);
        }
        this->gbz = gbwtgraph::GBZ(index, source);
        this->fragment_map = gbwt::FragmentMap(this->gbz.index.metadata, false);
    }

    // Return the nodes the path visits, in the orientation of its GBWT sequence.
    std::vector<gbwt::node_type> oriented_nodes(size_t path, bool is_reverse) const {
        std::vector<gbwt::node_type> result;
        for (nid_t id : this->paths[path].nodes) {
            result.push_back(fwd(id));
        }
        if (is_reverse) {
            std::reverse(result.begin(), result.end());
            for (gbwt::node_type& node : result) {
                node = gbwt::Node::reverse(node);
            }
        }
        return result;
    }

    // Return the GBWT sequence of the path in the orientation, with its offset in the
    // node of visit `visit`. The offset comes from following the sequence through the
    // GBWT.
    Haplotypes::sequence_type at_visit(size_t path, bool is_reverse, size_t visit) const {
        gbwt::size_type sequence_id = gbwt::Path::encode(path, is_reverse);
        gbwt::edge_type pos = this->gbz.index.start(sequence_id);
        for (size_t i = 0; i < visit; i++) {
            pos = this->gbz.index.LF(pos);
        }
        REQUIRE(pos.first == this->oriented_nodes(path, is_reverse)[visit]);
        return { sequence_id, pos.second };
    }

    // Return the walk of the path in the orientation from `start` to `end`. With a start
    // node, the walk starts at visit `visit` of the path, which must be a visit to `start`.
    SubchainWalk walk(size_t path, bool is_reverse, gbwt::node_type start, gbwt::node_type end, size_t visit) const {
        gbwt::CachedGBWT cache(this->gbz.index);
        Haplotypes::sequence_type sequence(gbwt::Path::encode(path, is_reverse), 0);
        if (start != gbwt::ENDMARKER) {
            REQUIRE(this->oriented_nodes(path, is_reverse)[visit] == start);
            sequence = this->at_visit(path, is_reverse, visit);
        }
        return SubchainWalk(cache, this->fragment_map, start, end, sequence);
    }

    // Return the visits `SubchainWalk::for_each_visit()` reports for the walk of the path.
    std::vector<Visit> reported_visits(size_t path, bool is_reverse, gbwt::node_type start, gbwt::node_type end, size_t visit) const {
        gbwt::CachedGBWT cache(this->gbz.index);
        Haplotypes::sequence_type sequence(gbwt::Path::encode(path, is_reverse), 0);
        if (start != gbwt::ENDMARKER) {
            REQUIRE(this->oriented_nodes(path, is_reverse)[visit] == start);
            sequence = this->at_visit(path, is_reverse, visit);
        }
        std::vector<Visit> result;
        bool complete = SubchainWalk::for_each_visit(cache, this->fragment_map, start, end, sequence, [&](gbwt::node_type node, bool starts_piece) {
            result.emplace_back(node, starts_piece);
            return true;
        });
        REQUIRE(complete);
        return result;
    }
};

// Return the pieces of the walk joined into one string.
std::string joined(const std::vector<std::string>& pieces) {
    std::string result;
    for (const std::string& piece : pieces) {
        result += piece;
    }
    return result;
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("A subchain walk follows a haplotype across its fragments", "[haplotypes][subchain_walk]") {
    WalkFixture fixture;
    const gbwtgraph::GBWTGraph& graph = fixture.gbz.graph;

    SECTION("a walk between two nodes continues past a fragment break") {
        SubchainWalk walk = fixture.walk(0, false, fwd(2), fwd(6), 1);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 2, 3, 4, 6 }));
        REQUIRE(walk.pieces() == 2);
        // Node 4 starts the second piece, and keeps all of its bases.
        REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "GTA" "TTTAAAC", "ACGTACG" "CAT" }));
        REQUIRE(walk.sequence_pieces(graph, ALL, ALL) == std::vector<std::string>({ "CCGGTA" "TTTAAAC", "ACGTACG" "CATCATG" }));
        std::vector<Visit> visits { { fwd(2), true }, { fwd(3), false }, { fwd(4), true }, { fwd(6), false } };
        REQUIRE(fixture.reported_visits(0, false, fwd(2), fwd(6), 1) == visits);
    }

    SECTION("a walk between two nodes in one fragment is one piece") {
        SubchainWalk walk = fixture.walk(2, false, fwd(2), fwd(6), 1);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 2, 5, 6 }));
        REQUIRE(walk.pieces() == 1);
        REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "GTA" "GGGCCCA" "CAT" }));
    }

    SECTION("the walk follows the fragment numbers, and the end node can start a piece") {
        // Fragment 0 of END_FIRST is path 4, and fragment 1 is path 3.
        SubchainWalk walk = fixture.walk(4, false, fwd(2), fwd(6), 1);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 2, 5, 6 }));
        REQUIRE(walk.pieces() == 2);
        REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "GTA" "GGGCCCA", "CAT" }));
    }

    SECTION("the walk skips an empty fragment") {
        SubchainWalk walk = fixture.walk(5, false, fwd(2), fwd(6), 1);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 2, 3, 4, 6 }));
        REQUIRE(walk.pieces() == 2);
        REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "GTATTTAAAC", "ACGTACGCAT" }));
        SubchainWalk broken = fixture.walk(0, false, fwd(2), fwd(6), 1);
        REQUIRE(walk == broken);
    }

    SECTION("a reverse walk continues with the previous fragment") {
        // Path 1 reversed visits 7, 6 and 4 in reverse, and path 0 reversed visits 3, 2 and 1.
        SubchainWalk walk = fixture.walk(1, true, rev(6), rev(2), 1);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 2, 3, 4, 6 }));
        // Path 7 reversed gets there past the empty path 6.
        REQUIRE(fixture.walk(7, true, rev(6), rev(2), 1) == walk);
        REQUIRE(walk.pieces() == 2);
        REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "ATG" "CGTACGT", "GTTTAAA" "TAC" }));
        std::vector<Visit> visits { { rev(6), true }, { rev(4), false }, { rev(3), true }, { rev(2), false } };
        REQUIRE(fixture.reported_visits(1, true, rev(6), rev(2), 1) == visits);
    }

    SECTION("a walk without a start node stops at the end of its fragment") {
        SubchainWalk walk = fixture.walk(0, false, gbwt::ENDMARKER, fwd(6), 0);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 1, 2, 3 }));
        REQUIRE(walk.pieces() == 1);
        REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "ACA" "CCGGTA" "TTTAAAC" }));
        REQUIRE(walk.sequence_pieces(graph, ALL, FLANK) == std::vector<std::string>({ "GATTACA" "CCGGTA" "TTTAAAC" }));

        SubchainWalk to_end = fixture.walk(2, false, gbwt::ENDMARKER, fwd(6), 0);
        REQUIRE(to_end.node_ids() == std::vector<nid_t>({ 1, 2, 5, 6 }));
        REQUIRE(to_end.sequence_pieces(graph, ALL, FLANK) == std::vector<std::string>({ "GATTACA" "CCGGTA" "GGGCCCA" "CAT" }));
    }

    SECTION("a walk without an end node stops at the end of its fragment") {
        SubchainWalk walk = fixture.walk(0, false, fwd(2), gbwt::ENDMARKER, 1);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 2, 3 }));
        REQUIRE(walk.pieces() == 1);
        REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "GTA" "TTTAAAC" }));

        // Node 6 is not the end node, so it keeps all of its bases.
        SubchainWalk whole = fixture.walk(2, false, fwd(2), gbwt::ENDMARKER, 1);
        REQUIRE(whole.node_ids() == std::vector<nid_t>({ 2, 5, 6, 7 }));
        REQUIRE(whole.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "GTA" "GGGCCCA" "CATCATG" "TGCATGC" }));
    }

    SECTION("a walk without start and end nodes covers one fragment") {
        SubchainWalk walk = fixture.walk(0, false, gbwt::ENDMARKER, gbwt::ENDMARKER, 0);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 1, 2, 3 }));
        REQUIRE(walk.pieces() == 1);
        REQUIRE(walk.sequence_pieces(graph, ALL, ALL) == std::vector<std::string>({ "GATTACA" "CCGGTA" "TTTAAAC" }));

        SubchainWalk second = fixture.walk(3, false, gbwt::ENDMARKER, gbwt::ENDMARKER, 0);
        REQUIRE(second.node_ids() == std::vector<nid_t>({ 6, 7 }));
        REQUIRE(second.sequence_pieces(graph, ALL, ALL) == std::vector<std::string>({ "CATCATG" "TGCATGC" }));
    }

    SECTION("a walk whose first node is the end node has no sequence") {
        SubchainWalk walk = fixture.walk(3, false, gbwt::ENDMARKER, fwd(6), 0);
        REQUIRE(walk.node_ids() == std::vector<nid_t>({ 6 }));
        REQUIRE(walk.pieces() == 0);
        REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK).empty());
        std::vector<Visit> visits { { fwd(6), true } };
        REQUIRE(fixture.reported_visits(3, false, gbwt::ENDMARKER, fwd(6), 0) == visits);
    }

    SECTION("an empty or missing sequence has an empty walk") {
        gbwt::CachedGBWT cache(fixture.gbz.index);
        Haplotypes::sequence_type empty(gbwt::Path::encode(6, false), 0);
        Haplotypes::sequence_type missing(fixture.gbz.index.sequences(), 0);
        for (Haplotypes::sequence_type sequence : { empty, missing }) {
            for (gbwt::node_type end : { gbwt::ENDMARKER, gbwt::node_type(fwd(6)) }) {
                SubchainWalk walk(cache, fixture.fragment_map, gbwt::ENDMARKER, end, sequence);
                REQUIRE(walk.node_ids().empty());
                REQUIRE(walk.pieces() == 0);
                REQUIRE(walk.sequence_pieces(graph, FLANK, FLANK).empty());
                REQUIRE(walk == SubchainWalk());
                size_t calls = 0;
                bool complete = SubchainWalk::for_each_visit(cache, fixture.fragment_map, gbwt::ENDMARKER, end, sequence, [&](gbwt::node_type, bool) {
                    calls++;
                    return true;
                });
                REQUIRE(complete);
                REQUIRE(calls == 0);
            }
        }
    }

    SECTION("walks with the same visits in different pieces are unequal") {
        SubchainWalk broken = fixture.walk(0, false, fwd(2), fwd(6), 1);
        SubchainWalk joined_walk = fixture.walk(8, false, fwd(2), fwd(6), 1);
        REQUIRE(joined_walk.pieces() == 1);
        REQUIRE(joined_walk.node_ids() == broken.node_ids());
        REQUIRE(joined(joined_walk.sequence_pieces(graph, FLANK, FLANK)) == joined(broken.sequence_pieces(graph, FLANK, FLANK)));
        REQUIRE_FALSE(joined_walk == broken);
        // Walks that differ are ordered one way, so sequences are grouped by equal walks.
        REQUIRE((broken < joined_walk) != (joined_walk < broken));
    }

    SECTION("a walk starts at the given visit to the start node") {
        // Path 9 visits node 3 twice.
        SubchainWalk first = fixture.walk(9, false, fwd(3), fwd(6), 2);
        REQUIRE(first.node_ids() == std::vector<nid_t>({ 3, 4, 6 }));
        REQUIRE(first.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "AAC" "ACGTACG" "TTTAAAC" "ACGTACG" "CAT" }));
        SubchainWalk second = fixture.walk(9, false, fwd(3), fwd(6), 4);
        REQUIRE(second.node_ids() == std::vector<nid_t>({ 3, 4, 6 }));
        REQUIRE(second.sequence_pieces(graph, FLANK, FLANK) == std::vector<std::string>({ "AAC" "ACGTACG" "CAT" }));
        REQUIRE_FALSE(first == second);
    }

    SECTION("the visitor can stop the walk") {
        gbwt::CachedGBWT cache(fixture.gbz.index);
        Haplotypes::sequence_type sequence = fixture.at_visit(0, false, 1);
        // Stop at the third visit, which is in the second fragment, and at the end node.
        for (size_t stop : { size_t(2), size_t(3) }) {
            std::vector<Visit> visits;
            bool complete = SubchainWalk::for_each_visit(cache, fixture.fragment_map, fwd(2), fwd(6), sequence, [&](gbwt::node_type node, bool starts_piece) {
                visits.emplace_back(node, starts_piece);
                return visits.size() <= stop;
            });
            REQUIRE_FALSE(complete);
            std::vector<Visit> expected { { fwd(2), true }, { fwd(3), false }, { fwd(4), true }, { fwd(6), false } };
            expected.resize(stop + 1);
            REQUIRE(visits == expected);
        }
    }
}

} // namespace vg::unittest
