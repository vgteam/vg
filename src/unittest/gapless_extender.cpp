/** \file
 *
 * Unit tests for gapless_extender.cpp, which implements haplotype-aware gapless seed extension.
 */

#include "../gapless_extender.hpp"
#include "../gbwt_helper.hpp"
#include "../json2pb.h"
#include "../utility.hpp"

#include "catch.hpp"

#include <map>
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

// Match: 1-9
// Mismatch: ACGT
// Insertion: + 1-9 str
// Deletion: - 1-9
Path get_path(const std::vector<std::pair<pos_t, std::string>>& mappings) {

    Path result;

    for (const std::pair<pos_t, std::string>& node : mappings) {
        Mapping& mapping = *(result.add_mapping());
        pos_t pos = node.first;
        mapping.mutable_position()->set_node_id(id(pos));
        mapping.mutable_position()->set_offset(offset(pos));
        mapping.mutable_position()->set_is_reverse(is_rev(pos));

        std::string edits = node.second;
        for (size_t i = 0; i < edits.length(); i++) {
            Edit& edit = *(mapping.add_edit());
            if (edits[i] > '0' && edits[i] <= '9') {
                int n = edits[i] - '0';
                edit.set_from_length(n);
                edit.set_to_length(n);
            } else if (edits[i] == '-') {
                i++;
                int n = edits[i] - '0';
                edit.set_from_length(n);
            } else if (edits[i] == '+') {
                i++;
                int n = edits[i] - '0';
                i++;
                edit.set_to_length(n);
                edit.set_sequence(edits.substr(i, n));
                i += n - 1;
            } else {
                edit.set_from_length(1);
                edit.set_to_length(1);
                edit.set_sequence(edits.substr(i, 1));
            }
        }
    }

    return result;
}

Alignment get_alignment(const std::vector<std::pair<pos_t, std::string>>& mappings, const std::string& sequence) {
    Alignment result;
    result.set_sequence(sequence);
    *(result.mutable_path()) = get_path(mappings);
    return result;
}

void paths_match(const Path& path, const Path& correct_path) {
    REQUIRE(path.mapping_size() == correct_path.mapping_size());
    for (size_t i = 0; i < path.mapping_size(); i++) {
        const Mapping& mapping = path.mapping(i);
        const Mapping& correct = correct_path.mapping(i);
        REQUIRE(make_pos_t(mapping.position()) == make_pos_t(correct.position()));
        REQUIRE(mapping.edit_size() == correct.edit_size());
        for (size_t j = 0; j < mapping.edit_size(); j++) {
            REQUIRE(mapping.edit(j).from_length() == correct.edit(j).from_length());
            REQUIRE(mapping.edit(j).to_length() == correct.edit(j).to_length());
            REQUIRE(mapping.edit(j).sequence() == correct.edit(j).sequence());
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
        paths_match(result.front().to_path(*(extender.graph), read), get_path(correct_alignment));

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
        paths_match(result[i].to_path(*(extender.graph), read), get_path(correct_alignments[i]));
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
        paths_match(result[i].to_path(*(extender.graph), read), get_path(correct_extensions[i]));
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

void check_haplotypes(const GaplessExtender& extender, const std::vector<nid_t>& nodes, const std::map<std::vector<handle_t>, std::string>& correct_haplotypes) {
    const gbwtgraph::GBWTGraph& graph = *(extender.graph);
    SubHandleGraph subgraph(extender.graph);
    for (nid_t node : nodes) {
        subgraph.add_handle(graph.get_handle(node, false));
    }

    std::vector<std::vector<handle_t>> haplotype_paths;
    bdsg::HashGraph unfolded;
    extender.unfold_haplotypes(subgraph, haplotype_paths, unfolded);

    // We assume that the correct paths are in canonical orientation (smaller than the reverse).
    std::vector<bool> flipped(haplotype_paths.size(), false);
    for (size_t i = 0; i < haplotype_paths.size(); i++) {
        std::vector<handle_t> reverse = haplotype_paths[i];
        std::reverse(reverse.begin(), reverse.end());
        for (handle_t& handle : reverse) {
            handle = extender.graph->flip(handle);
        }
        if (reverse < haplotype_paths[i]) {
            haplotype_paths[i] = reverse;
            flipped[i] = true;
        }
    }

    REQUIRE(haplotype_paths.size() == correct_haplotypes.size());
    REQUIRE(unfolded.get_node_count() == 2 * correct_haplotypes.size());
    for (size_t i = 0; i < haplotype_paths.size(); i++) {
        auto iter = correct_haplotypes.find(haplotype_paths[i]);
        REQUIRE(iter != correct_haplotypes.end());
        std::string forward = unfolded.get_sequence(unfolded.get_handle(2 * i + 1, false));
        std::string reverse = unfolded.get_sequence(unfolded.get_handle(2 * i + 2, false));
        if (flipped[i]) {
            forward = reverse_complement(forward);
            REQUIRE(forward == iter->second);
            REQUIRE(reverse == iter->second);
        } else {
            REQUIRE(forward == iter->second);
            reverse = reverse_complement(reverse);
            REQUIRE(reverse == iter->second);
        }
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

    // We also test that we avoid finding the same best alignment from multiple seeds.
    SECTION("secondary alignment has more mismatches") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(2, false, 0), 1 }, // First seed for the best alignment (1 mismatch).
            { make_pos_t(4, false, 0), 2 }, // Second seed for the best alignment (1 mismatch).
            { make_pos_t(4, false, 0), 1 }  // Seed for the secondary alignment (2 mismatches).
        };
        std::string read = "GAGGA";
        size_t error_bound = 2;
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_alignments {
            {
                { make_pos_t(1, false, 0), "1" },
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "2A" }
            },
            {
                { make_pos_t(1, false, 0), "1" },
                { make_pos_t(4, false, 0), "A2" },
                { make_pos_t(5, false, 0), "A" }
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

TEST_CASE("Haplotype unfolding", "[gapless_extender]") {

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

    SECTION("normal subgraph") {
        std::vector<nid_t> nodes {
            5, 6, 7, 8
        };
        std::map<std::vector<handle_t>, std::string> correct_haplotypes {
            {
                { gbwt_graph.get_handle(5, false), gbwt_graph.get_handle(6, false), gbwt_graph.get_handle(7, false) },
                "TAC"
            },
            {
                { gbwt_graph.get_handle(5, false), gbwt_graph.get_handle(6, false), gbwt_graph.get_handle(8, false) },
                "TAA"
            }
        };
        check_haplotypes(extender, nodes, correct_haplotypes);
    }

    SECTION("right-maximal extensions") {
        std::vector<nid_t> nodes {
            5, 6, 7
        };
        std::map<std::vector<handle_t>, std::string> correct_haplotypes {
            {
                { gbwt_graph.get_handle(5, false), gbwt_graph.get_handle(6, false) },
                "TA"
            },
            {
                { gbwt_graph.get_handle(5, false), gbwt_graph.get_handle(6, false), gbwt_graph.get_handle(7, false) },
                "TAC"
            }
        };
        check_haplotypes(extender, nodes, correct_haplotypes);
    }

    SECTION("left-maximal extensions") {
        std::vector<nid_t> nodes {
            2, 4, 5, 6, 7
        };
        std::map<std::vector<handle_t>, std::string> correct_haplotypes {
            {
                { gbwt_graph.get_handle(2, false), gbwt_graph.get_handle(4, false), gbwt_graph.get_handle(5, false), gbwt_graph.get_handle(6, false) },
                "AGGGTA"
            },
            {
                { gbwt_graph.get_handle(4, false), gbwt_graph.get_handle(5, false), gbwt_graph.get_handle(6, false), gbwt_graph.get_handle(7, false) },
                "GGGTAC"
            }
        };
        check_haplotypes(extender, nodes, correct_haplotypes);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Alignment transformations", "[gapless_extender]") {

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

    // 0: 1, 2, 4, 5, 6, 8, 9 / GAGGGTAAA
    // 1: 1, 4, 5, 6, 7, 9 / GGGGTACA
    std::vector<std::vector<handle_t>> haplotype_paths(2);
    std::vector<std::string> haplotype_sequences(2);
    for (gbwt::node_type node : alt_path) {
        handle_t handle = gbwtgraph::GBWTGraph::node_to_handle(node);
        haplotype_paths[0].push_back(handle);
        haplotype_sequences[0] += gbwt_graph.get_sequence(handle);
    }
    for (gbwt::node_type node : short_path) {
        handle_t handle = gbwtgraph::GBWTGraph::node_to_handle(node);
        haplotype_paths[1].push_back(handle);
        haplotype_sequences[1] += gbwt_graph.get_sequence(handle);
    }

    SECTION("simple cases") {
        std::vector<std::vector<std::pair<pos_t, std::string>>> unfolded_alignments {
            { // Haplotype 1 forward; substitution in the middle + insertion at node boundary.
                { make_pos_t(1, false, 2), "1C2+1G1" }
            },
            { // Haplotype 1 reverse: substitution over node boundary.
                { make_pos_t(2, false, 2), "1GG2" }
            },
            { // Haplotype 2 forward: start in the middle + deletion over node boundary.
                { make_pos_t(3, false, 2), "1-22" }
            },
            { // Haplotype 2 reverse: substitution of a node + deletion of a node + end in the middle.
                { make_pos_t(4, false, 1), "1C-12" }
            }
        };
        std::vector<std::string> reads {
            "GCGTGA", "TGGCC", "GAC", "CCGG"
        };
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_alignments {
            {
                { make_pos_t(4, false, 0), "1C1" },
                { make_pos_t(5, false, 0), "1" },
                { make_pos_t(6, false, 0), "+1G1" }
            },
            {
                { make_pos_t(6, true, 0), "1" },
                { make_pos_t(5, true, 0), "G" },
                { make_pos_t(4, true, 0), "G2" }
            },
            {
                { make_pos_t(4, false, 1), "1-1" },
                { make_pos_t(5, false, 0), "-1" },
                { make_pos_t(6, false, 0), "1" },
                { make_pos_t(7, false, 0), "1" }
            },
            {
                { make_pos_t(7, true, 0), "1" },
                { make_pos_t(6, true, 0), "C" },
                { make_pos_t(5, true, 0), "-1" },
                { make_pos_t(4, true, 0), "2" }
            }
        };
        for (size_t i = 0; i < unfolded_alignments.size(); i++) {
            Alignment source = get_alignment(unfolded_alignments[i], reads[i]);
            extender.transform_alignment(source, haplotype_paths);
            Path target = get_path(correct_alignments[i]);
            paths_match(source.path(), target);
        }
    }
}

//------------------------------------------------------------------------------

}
}
