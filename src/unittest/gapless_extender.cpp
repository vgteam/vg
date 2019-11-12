/** \file
 *
 * Unit tests for gapless_extender.cpp, which implements haplotype-aware gapless seed extension.
 */

#include "../gapless_extender.hpp"
#include "../gbwt_helper.hpp"
#include "../json2pb.h"

#include "catch.hpp"

#include <vector>


namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

/*
  A toy graph for GA(T|GGG)TA(C|A)A with some additional edges.
*/
const std::string gapless_extender_graph = R"(
{
    "node": [
        {"id": 1, "sequence": "G"},
        {"id": 2, "sequence": "A"},
        {"id": 3, "sequence": "T"},
        {"id": 4, "sequence": "GGG"},
        {"id": 5, "sequence": "T"},
        {"id": 6, "sequence": "A"},
        {"id": 7, "sequence": "C"},
        {"id": 8, "sequence": "A"},
        {"id": 9, "sequence": "A"}
    ],
    "edge": [
        {"from": 1, "to": 2},
        {"from": 1, "to": 4},
        {"from": 1, "to": 6},
        {"from": 2, "to": 3},
        {"from": 2, "to": 4},
        {"from": 3, "to": 5},
        {"from": 4, "to": 5},
        {"from": 5, "to": 6},
        {"from": 6, "to": 7},
        {"from": 6, "to": 8},
        {"from": 7, "to": 9},
        {"from": 8, "to": 9}
    ]
}
)";

gbwt::vector_type alt_path {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(2, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(8, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
};

gbwt::vector_type short_path {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(4, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(5, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(6, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(7, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(9, false))
};

// Build GBWT with 2x short_path and alt_path.
gbwt::GBWT build_gbwt_index() {
    std::vector<gbwt::vector_type> gbwt_threads {
        short_path, alt_path, short_path
    };

    return get_gbwt(gbwt_threads);
}

void same_position(const Position& pos, const Position& correct) {
    REQUIRE(pos.node_id() == correct.node_id());
    REQUIRE(pos.is_reverse() == correct.is_reverse());
    REQUIRE(pos.offset() == correct.offset());
}

std::vector<std::pair<pos_t, size_t>> normalize_seeds(std::vector<std::pair<pos_t, size_t>>& seeds) {
    GaplessExtender::cluster_type cluster;
    for (auto seed : seeds) {
        cluster.insert(GaplessExtender::to_seed(seed.first, seed.second));
    }
    std::vector<std::pair<pos_t, size_t>> result;
    for (auto seed : cluster) {
        result.emplace_back(GaplessExtender::get_pos(seed), GaplessExtender::get_read_offset(seed));
    }
    std::sort(result.begin(), result.end());
    return result;
}

void correct_score(const GaplessExtension& extension, const Aligner& aligner) {
    int32_t expected_score = (extension.length() - extension.mismatches()) * aligner.match;
    expected_score -= extension.mismatches() * aligner.mismatch;
    expected_score += extension.left_full * aligner.full_length_bonus;
    expected_score += extension.right_full * aligner.full_length_bonus;
    REQUIRE(extension.score == expected_score);
}

void alignment_matches(const Path& path, const std::vector<std::pair<pos_t, std::string>>& alignment) {
    REQUIRE(path.mapping_size() == alignment.size());
    for (size_t i = 0; i < path.mapping_size(); i++) {
        const Mapping& m = path.mapping(i);
        REQUIRE(make_pos_t(m.position()) == alignment[i].first);
        const std::string& edits = alignment[i].second;
        REQUIRE(m.edit_size() == edits.length());
        for (size_t j = 0; j < m.edit_size(); j++) {
            if (edits[j] > '0' && edits[j] <= '9') {
                int n = edits[j] - '0';
                bool match_length_ok = (m.edit(j).from_length() == n &&
                                        m.edit(j).to_length() == n &&
                                        m.edit(j).sequence().empty());
                REQUIRE(match_length_ok);
            } else {
                std::string s = edits.substr(j, 1);
                bool mismatch_ok = (m.edit(j).from_length() == 1 && m.edit(j).to_length() == 1 && m.edit(j).sequence() == s);
                REQUIRE(mismatch_ok);
            }
        }
    }
}

void full_length_match(const std::vector<std::pair<pos_t, size_t>>& seeds, const std::string& read, const std::vector<std::pair<pos_t, std::string>>& correct_alignment, const GaplessExtender& extender, size_t error_bound, bool check_seeds) {
    GaplessExtender::cluster_type cluster;
    for (auto seed : seeds) {
        cluster.insert(GaplessExtender::to_seed(seed.first, seed.second));
    }
    auto result = extender.extend(cluster, read, error_bound);

    // Empty correct alignment indicates that there should not be a full-length alignment.
    if (correct_alignment.empty()) {
        for (auto& extension : result) {
            if (extension.full()) {
                REQUIRE(extension.mismatches() > error_bound);
            }
        }
    } else {
        REQUIRE(result.size() == 1);
        REQUIRE(!result.front().empty());
        REQUIRE(result.front().full());
        REQUIRE(result.front().mismatches() <= error_bound);
        correct_score(result.front(), *(extender.aligner));
        alignment_matches(result.front().to_path(*(extender.graph), read), correct_alignment);

        // This extension should contain all the seeds. Check that contains() works correctly.
        if (check_seeds) {
            for (auto seed : cluster) {
                REQUIRE(result.front().contains(*(extender.graph), seed));
            }
        }
    }
}

void full_length_matches(const std::vector<std::pair<pos_t, size_t>>& seeds, const std::string& read, const std::vector<std::vector<std::pair<pos_t, std::string>>>& correct_alignments, const GaplessExtender& extender, size_t error_bound) {
    GaplessExtender::cluster_type cluster;
    for (auto seed : seeds) {
        cluster.insert(GaplessExtender::to_seed(seed.first, seed.second));
    }
    auto result = extender.extend(cluster, read, error_bound);

    REQUIRE(result.size() == correct_alignments.size());
    for (size_t i = 0; i < result.size(); i++) {
        REQUIRE(!result[i].empty());
        REQUIRE(result[i].full());
        REQUIRE(result[i].mismatches() <= error_bound);
        correct_score(result[i], *(extender.aligner));
        alignment_matches(result[i].to_path(*(extender.graph), read), correct_alignments[i]);
    }
}

void partial_matches(const std::vector<std::pair<pos_t, size_t>>& seeds, const std::string& read, const std::vector<std::vector<std::pair<pos_t, std::string>>>& correct_extensions, const std::vector<size_t>& correct_offsets, const GaplessExtender& extender, size_t error_bound) {
    GaplessExtender::cluster_type cluster;
    for (auto seed : seeds) {
        cluster.insert(GaplessExtender::to_seed(seed.first, seed.second));
    }
    auto result = extender.extend(cluster, read, error_bound, false);

    REQUIRE(result.size() == correct_extensions.size());
    for (size_t i = 0; i < result.size(); i++) {
        REQUIRE(!(result[i].empty()));
        if (result[i].full()) {
            REQUIRE(result[i].mismatches() > error_bound);
        }
        REQUIRE(result[i].read_interval.first == correct_offsets[i]);
        correct_score(result.front(), *(extender.aligner));
        alignment_matches(result[i].to_path(*(extender.graph), read), correct_extensions[i]);
    }
}

void trimmed_extensions(std::vector<GaplessExtension>& extensions, std::vector<GaplessExtension>& correct, const GaplessExtender& extender, size_t error_bound) {
    for (GaplessExtension& extension : extensions) {
        extension.state = extender.graph->bd_find(extension.path);
    }
    extender.trim(extensions, error_bound);

    REQUIRE(extensions.size() == correct.size());
    for (size_t i = 0; i < extensions.size(); i++) {
        correct_score(extensions[i], *(extender.aligner));
        REQUIRE(extensions[i].path == correct[i].path);
        REQUIRE(extensions[i].offset == correct[i].offset);
        correct[i].state = extender.graph->bd_find(correct[i].path);
        REQUIRE(extensions[i].state == correct[i].state);
        REQUIRE(extensions[i].read_interval == correct[i].read_interval);
        REQUIRE(extensions[i].mismatch_positions == correct[i].mismatch_positions);
        REQUIRE(extensions[i].score == correct[i].score);
        REQUIRE(extensions[i].left_full == correct[i].left_full);
        REQUIRE(extensions[i].right_full == correct[i].right_full);
    }
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Gapless extensions report correct positions", "[gapless_extender]") {

    // Build an XG index.
    Graph graph;
    json2pb(graph, gapless_extender_graph.c_str(), gapless_extender_graph.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph));

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph(gbwt_index, xg_index);

    SECTION("starts and ends at node boundaries") {
        GaplessExtension extension {
            {
                gbwt_graph.get_handle(1, false),
                gbwt_graph.get_handle(4, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 4 }, { },
            0, false, false,
            false, false, 0, 0
        };
        Position correct_start = make_position(1, false, 0);
        Position correct_tail = make_position(4, false, 3);
        same_position(extension.starting_position(gbwt_graph), correct_start);
        same_position(extension.tail_position(gbwt_graph), correct_tail);
    }

    SECTION("starts in the middle") {
        GaplessExtension extension {
            {
                gbwt_graph.get_handle(4, false),
                gbwt_graph.get_handle(5, false)
            },
            1, gbwt::BidirectionalState(),
            { 0, 3 }, { },
            0, false, false,
            false, false, 0, 0
        };
        Position correct_start = make_position(4, false, 1);
        Position correct_tail = make_position(5, false, 1);
        same_position(extension.starting_position(gbwt_graph), correct_start);
        same_position(extension.tail_position(gbwt_graph), correct_tail);
    }

    SECTION("ends in the middle") {
        GaplessExtension extension {
            {
                gbwt_graph.get_handle(1, false),
                gbwt_graph.get_handle(4, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 3 }, { },
            0, false, false,
            false, false, 0, 0
        };
        Position correct_start = make_position(1, false, 0);
        Position correct_tail = make_position(4, false, 2);
        same_position(extension.starting_position(gbwt_graph), correct_start);
        same_position(extension.tail_position(gbwt_graph), correct_tail);
    }

    SECTION("starts and ends in the middle") {
        GaplessExtension extension {
            {
                gbwt_graph.get_handle(4, false)
            },
            1, gbwt::BidirectionalState(),
            { 0, 1 }, { },
            0, false, false,
            false, false, 0, 0
        };
        Position correct_start = make_position(4, false, 1);
        Position correct_tail = make_position(4, false, 2);
        same_position(extension.starting_position(gbwt_graph), correct_start);
        same_position(extension.tail_position(gbwt_graph), correct_tail);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Redundant seeds are removed from a cluster", "[gapless_extender]") {

    SECTION("various types of redundant seeds") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 2), static_cast<size_t>(1) },
            { make_pos_t(4, false, 3), static_cast<size_t>(2) },
            { make_pos_t(5, true, 1), static_cast<size_t>(2) },
            { make_pos_t(5, true, 2), static_cast<size_t>(3) },
            { make_pos_t(6, false, 0), static_cast<size_t>(0) },
            { make_pos_t(6, false, 1), static_cast<size_t>(1) }
        };
        std::vector<std::pair<pos_t, size_t>> correct_seeds {
            { make_pos_t(4, false, 1), static_cast<size_t>(0) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(0) }
        };
        std::vector<std::pair<pos_t, size_t>> extracted_seeds = normalize_seeds(seeds);
        REQUIRE(extracted_seeds.size() == correct_seeds.size());
        REQUIRE(extracted_seeds == correct_seeds);
    }

    SECTION("various types of distinct seeds") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 2), static_cast<size_t>(1) },
            { make_pos_t(4, true, 3), static_cast<size_t>(2) },
            { make_pos_t(5, true, 1), static_cast<size_t>(2) },
            { make_pos_t(5, false, 2), static_cast<size_t>(3) },
            { make_pos_t(6, false, 0), static_cast<size_t>(0) },
            { make_pos_t(6, true, 1), static_cast<size_t>(1) }
        };
        std::vector<std::pair<pos_t, size_t>> correct_seeds {
            { make_pos_t(4, false, 1), static_cast<size_t>(0) },
            { make_pos_t(4, true, 1), static_cast<size_t>(0) },
            { make_pos_t(5, false, 0), static_cast<size_t>(1) },
            { make_pos_t(5, true, 0), static_cast<size_t>(1) },
            { make_pos_t(6, false, 0), static_cast<size_t>(0) },
            { make_pos_t(6, true, 0), static_cast<size_t>(0) }
        };
        std::vector<std::pair<pos_t, size_t>> extracted_seeds = normalize_seeds(seeds);
        REQUIRE(extracted_seeds.size() == correct_seeds.size());
        REQUIRE(extracted_seeds == correct_seeds);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Full-length alignments", "[gapless_extender]") {

    // Build an XG index.
    Graph graph;
    json2pb(graph, gapless_extender_graph.c_str(), gapless_extender_graph.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph));

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph(gbwt_index, xg_index);

    // And finally wrap it in a GaplessExtender with an Aligner.
    Aligner aligner;
    GaplessExtender extender(gbwt_graph, aligner);

    SECTION("read starting in the middle of a node matches exactly") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 2), 0 },
            { make_pos_t(6, false, 0), 2 }
        };
        std::string read = "GTACA";
        std::vector<std::pair<pos_t, std::string>> correct_alignment {
            { make_pos_t(4, false, 2), "1" },
            { make_pos_t(5, false, 0), "1" },
            { make_pos_t(6, false, 0), "1" },
            { make_pos_t(7, false, 0), "1" },
            { make_pos_t(9, false, 0), "1" }
        };
        size_t error_bound = 0;
        full_length_match(seeds, read, correct_alignment, extender, error_bound, true);
    }

    SECTION("read matches with errors") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(5, false, 0), 4 },
            { make_pos_t(4, false, 2), 3 }
        };
        std::string read = "GGAGTAC";
        std::vector<std::pair<pos_t, std::string>> correct_alignment {
            { make_pos_t(1, false, 0), "1" },
            { make_pos_t(4, false, 0), "1A1" },
            { make_pos_t(5, false, 0), "1" },
            { make_pos_t(6, false, 0), "1" },
            { make_pos_t(7, false, 0), "1" }
        };
        size_t error_bound = 1;
        full_length_match(seeds, read, correct_alignment, extender, error_bound, true);
    }

    SECTION("false seeds do not matter") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(5, false, 0), 4 },
            { make_pos_t(4, false, 2), 3 },
            { make_pos_t(2, false, 0), 0 }
        };
        std::string read = "GGAGTAC";
        std::vector<std::pair<pos_t, std::string>> correct_alignment {
            { make_pos_t(1, false, 0), "1" },
            { make_pos_t(4, false, 0), "1A1" },
            { make_pos_t(5, false, 0), "1" },
            { make_pos_t(6, false, 0), "1" },
            { make_pos_t(7, false, 0), "1" }
        };
        size_t error_bound = 1;
        full_length_match(seeds, read, correct_alignment, extender, error_bound, false);
    }

    SECTION("read matches reverse complement and ends within a node") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(5, true, 0), 2 },
            { make_pos_t(6, true, 0), 1 }
        };
        std::string read = "GTACT";
        std::vector<std::pair<pos_t, std::string>> correct_alignment {
            { make_pos_t(7, true, 0), "1" },
            { make_pos_t(6, true, 0), "1" },
            { make_pos_t(5, true, 0), "1" },
            { make_pos_t(4, true, 0), "1T" }
        };
        size_t error_bound = 1;
        full_length_match(seeds, read, correct_alignment, extender, error_bound, true);
    }

    SECTION("there is no full-length alignment") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(5, false, 0), 4 },
            { make_pos_t(4, false, 2), 3 }
        };
        std::string read = "AGAGTAC";
        size_t error_bound = 1;
        full_length_match(seeds, read, { }, extender, error_bound, false);
    }

    SECTION("there is a secondary alignment") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 2), 0 },
            { make_pos_t(6, true, 0), 1 }
        };
        std::string read = "GTAC";
        size_t error_bound = 1;
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_alignments {
            {
                { make_pos_t(4, false, 2), "1" },
                { make_pos_t(5, false, 0), "1" },
                { make_pos_t(6, false, 0), "1" },
                { make_pos_t(7, false, 0), "1" }
            },
            {
                { make_pos_t(7, true, 0), "1" },
                { make_pos_t(6, true, 0), "1" },
                { make_pos_t(5, true, 0), "1" },
                { make_pos_t(4, true, 0), "1" }
            }
        };
        full_length_matches(seeds, read, correct_alignments, extender, error_bound);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Partial alignments without trimming", "[gapless_extender]") {

    // Build an XG index.
    Graph graph;
    json2pb(graph, gapless_extender_graph.c_str(), gapless_extender_graph.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph));

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph(gbwt_index, xg_index);

    // And finally wrap it in a GaplessExtender with an Aligner.
    Aligner aligner;
    GaplessExtender extender(gbwt_graph, aligner);

    SECTION("exact matching") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 0), 1 },  // Match [0, 4); read border / mismatch
            { make_pos_t(2, false, 0), 7 },  // Match [6, 9); graph border / mismatch
            { make_pos_t(5, false, 0), 11 }, // Match [10, 13); mismatch / mismatch
            { make_pos_t(7, false, 0), 15 }, // Match [14, 17); mismatch / graph border
            { make_pos_t(6, false, 0), 20 }  // Match [19, 22); mismatch / read border
        };
        std::string read = "AGGGxCGAGxGTAxACAAxTAA";
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_extensions {
            {
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "3" }
            },
            {
                { make_pos_t(1, false, 0), "1" },
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "1" }
            },
            {
                { make_pos_t(4, false, 2), "1" },
                { make_pos_t(5, false, 0), "1" },
                { make_pos_t(6, false, 0), "1" }
            },
            {
                { make_pos_t(6, false, 0), "1" },
                { make_pos_t(7, false, 0), "1" },
                { make_pos_t(9, false, 0), "1" }
            },
            {
                { make_pos_t(5, false, 0), "1" },
                { make_pos_t(6, false, 0), "1" },
                { make_pos_t(8, false, 0), "1" }
            }
        };
        std::vector<size_t> correct_offsets {
            static_cast<size_t>(0),
            static_cast<size_t>(6),
            static_cast<size_t>(10),
            static_cast<size_t>(14),
            static_cast<size_t>(19)
        };
        size_t error_bound = 0;
        partial_matches(seeds, read, correct_extensions, correct_offsets, extender, error_bound);
    }

    SECTION("exact matching with mismatches in the initial node") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 2), 4 }
        };
        std::string read = "xAGxGTAx";
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_extensions {
            {
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "1x1" },
                { make_pos_t(5, false, 0), "1" },
                { make_pos_t(6, false, 0), "1" }
            }
        };
        std::vector<size_t> correct_offsets {
            static_cast<size_t>(1)
        };
        size_t error_bound = 0;
        partial_matches(seeds, read, correct_extensions, correct_offsets, extender, error_bound);
    }

    SECTION("approximate matching") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 2), 4 }
        };
        std::string read = "xAGGGTxAx";
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_extensions {
            {
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "3" },
                { make_pos_t(5, false, 0), "1" },
                { make_pos_t(6, false, 0), "x" },
                { make_pos_t(8, false, 0), "1" }
            }
        };
        std::vector<size_t> correct_offsets {
            static_cast<size_t>(1)
        };
        size_t error_bound = 1;
        partial_matches(seeds, read, correct_extensions, correct_offsets, extender, error_bound);
    }

    SECTION("removing duplicates") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(2, false, 0), 1 },
            { make_pos_t(4, false, 2), 4 }
        };
        std::string read = "xAGGGTxAx";
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_extensions {
            {
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "3" },
                { make_pos_t(5, false, 0), "1" },
                { make_pos_t(6, false, 0), "x" },
                { make_pos_t(8, false, 0), "1" }
            }
        };
        std::vector<size_t> correct_offsets {
            static_cast<size_t>(1)
        };
        size_t error_bound = 1;
        partial_matches(seeds, read, correct_extensions, correct_offsets, extender, error_bound);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Trimming mismatches", "[gapless_extender]") {

    // Build an XG index.
    Graph graph;
    json2pb(graph, gapless_extender_graph.c_str(), gapless_extender_graph.size());
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph));

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph(gbwt_index, xg_index);

    // And finally wrap it in a GaplessExtender with an Aligner.
    Aligner aligner;
    GaplessExtender extender(gbwt_graph, aligner);

    SECTION("basic trimming") {
        std::vector<GaplessExtension> extensions {
            { // Trim right end.
                {
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false)
                },
                static_cast<size_t>(1), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(5, 9), { static_cast<size_t>(7) },
                static_cast<int32_t>(-1), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            },
            { // Trim left end.
                {
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(0), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(2, 9), { static_cast<size_t>(3) },
                static_cast<int32_t>(2), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            },
            { // Trim both ends.
                {
                    gbwt_graph.get_handle(2, false),
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(0), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(1, 9), { static_cast<size_t>(2), static_cast<size_t>(7) },
                static_cast<int32_t>(-2), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            }
        };
        std::vector<GaplessExtension> correct {
            { // Trim both ends.
                {
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false)
                },
                static_cast<size_t>(1), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(3, 7), { },
                static_cast<int32_t>(4), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            },
            { // Trim left end.
                {
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(2), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(4, 9), { },
                static_cast<int32_t>(5), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            },
            { // Trim right end.
                {
                    gbwt_graph.get_handle(4, false),
                },
                static_cast<size_t>(1), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(5, 7), { },
                static_cast<int32_t>(2), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            }
        };
        size_t error_bound = 0;
        trimmed_extensions(extensions, correct, extender, error_bound);
    }

    SECTION("trimming full-length alignments") {
        std::vector<GaplessExtension> extensions {
            { // Number of mismatches is below the bound.
                {
                    gbwt_graph.get_handle(2, false),
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(0), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(0, 8), { static_cast<size_t>(6), static_cast<size_t>(7) },
                static_cast<int32_t>(8), true, true,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            },
            { // Full-length bonus means more than a single mismatch.
                {
                    gbwt_graph.get_handle(2, false),
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(0), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(0, 8), { static_cast<size_t>(0), static_cast<size_t>(6), static_cast<size_t>(7) },
                static_cast<int32_t>(3), true, true,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            }
        };
        std::vector<GaplessExtension> correct {
            { // Full-length bonus means more than a single mismatch.
                {
                    gbwt_graph.get_handle(2, false),
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false)
                },
                static_cast<size_t>(0), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(0, 6), { static_cast<size_t>(0) },
                static_cast<int32_t>(6), true, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            },
            { // Number of mismatches is below the bound.
                {
                    gbwt_graph.get_handle(2, false),
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(0), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(0, 8), { static_cast<size_t>(6), static_cast<size_t>(7) },
                static_cast<int32_t>(8), true, true,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            }
        };
        size_t error_bound = 2;
        trimmed_extensions(extensions, correct, extender, error_bound);
    }

    SECTION("removing duplicates and empty extensions") {
        std::vector<GaplessExtension> extensions {
            { // First duplicate.
                {
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(0), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(4, 11), { static_cast<size_t>(5) },
                static_cast<int32_t>(2), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            },
            { // Second duplicate.
                {
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(1), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(5, 11), { static_cast<size_t>(5) },
                static_cast<int32_t>(1), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            },
            { // Empty extension.
                {
                },
                static_cast<size_t>(0), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(0, 0), { },
                static_cast<int32_t>(0), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            }
        };
        std::vector<GaplessExtension> correct {
            { // Duplicates.
                {
                    gbwt_graph.get_handle(4, false),
                    gbwt_graph.get_handle(5, false),
                    gbwt_graph.get_handle(6, false),
                    gbwt_graph.get_handle(8, false),
                    gbwt_graph.get_handle(9, false)
                },
                static_cast<size_t>(2), gbwt::BidirectionalState(),
                std::pair<size_t, size_t>(6, 11), { },
                static_cast<int32_t>(5), false, false,
                false, false, static_cast<int32_t>(0), static_cast<int32_t>(0)
            }
        };
        size_t error_bound = 0;
        trimmed_extensions(extensions, correct, extender, error_bound);
    }
}

//------------------------------------------------------------------------------

}
}
