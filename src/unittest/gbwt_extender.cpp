/** \file
 *
 * Unit tests for gapless_extender.cpp, which implements haplotype-aware gapless seed extension.
 */

#include "../gbwt_extender.hpp"
#include "../gbwt_helper.hpp"
#include "vg/io/json2pb.h"
#include "../utility.hpp"
#include "../vg.hpp"

#include <bdsg/hash_graph.hpp>

#include "catch.hpp"
#include "randomness.hpp"

#include <map>
#include <unordered_set>
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

// Build a GBWTGraph using the provided GBWT index.
gbwtgraph::GBWTGraph build_gbwt_graph(const gbwt::GBWT& gbwt_index) {
    Graph graph;
    json2pb(graph, gapless_extender_graph.c_str(), gapless_extender_graph.size());
    VG vg_graph(graph);
    return gbwtgraph::GBWTGraph(gbwt_index, vg_graph);
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

void correct_score(const WFAAlignment& alignment, const Aligner& aligner) {
    int32_t expected_score = 0;
    for (auto& edit : alignment.edits) {
        switch (edit.first) {
        case WFAAlignment::match:
            expected_score += edit.second * aligner.match;
            break;
        case WFAAlignment::mismatch:
            expected_score -= edit.second * aligner.mismatch;
            break;
        case WFAAlignment::insertion:
            // Fall-through
        case WFAAlignment::deletion:
            expected_score -= (aligner.gap_open + (edit.second - 1) * aligner.gap_extension);
            break;
        default:
            REQUIRE(false);
        }
    }
    REQUIRE(alignment.score == expected_score);
}

//#define debug_conversion_tests

/// Convert a vg Path into a WFAAlignment, for testing and comparison.
/// Does not set score or sequence offset.
WFAAlignment path_to_wfa_alignment(const Path& path, const HandleGraph& graph) {
    WFAAlignment result {
        {},
        {},
        0,
        0,
        0,
        0,
        true
    };
    
    if (path.mapping_size() == 0) {
        // THis is actually the empty path.
        return result;
    }
    
    result.node_offset = path.mapping(0).position().offset();
    
    for (auto& mapping : path.mapping()) {
        // Remember what node we are on. We assume one mapping per node.
        result.path.push_back(graph.get_handle(mapping.position().node_id(), mapping.position().is_reverse()));
        for (auto& edit : mapping.edit()) {
            result.length += edit.to_length();
            
            if (edit.from_length() == edit.to_length()) {
                if (edit.sequence().empty()) {
                    // A match
                    result.append(WFAAlignment::match, edit.to_length());
                } else {
                    // A mismatch
                    result.append(WFAAlignment::mismatch, edit.to_length());
                }
            } else if (edit.from_length() == 0) {
                // An insertion
                result.append(WFAAlignment::insertion, edit.to_length());
            } else if (edit.to_length() == 0) {
                // A deletion
                result.append(WFAAlignment::deletion, edit.from_length());
            } else {
                // Something weird
                REQUIRE(false);
            }
        }
    }
    
    return result;
}

/// When generating a test alignmentm what happened to each base in the graph?
enum class GraphBaseFate {
    Match = 0,
    Mismatch,
    Deletion,
    InsertionMatch,
    InsertionMismatch,
    MatchInsertion,
    MismatchInsertion,
    End // Represents that we have run out of options.
};

/// What's a plan of base fates that we can turn into an alignment?
using alignment_plan_t = std::vector<GraphBaseFate>;

/// Make a random alignment plan
alignment_plan_t make_random_plan(size_t length, default_random_engine& engine) {
    uniform_int_distribution<int> fate_distribution(0, (int)GraphBaseFate::End - 1);
    alignment_plan_t plan;
    plan.reserve(length);
    for (size_t i = 0; i < length; i++) {
        plan.push_back((GraphBaseFate)fate_distribution(engine));
    }
    return plan;
}

/// Move all insertions at node boundaries into the later node in the path, if possible.
void send_insertions_right(Path& path) {
    for (size_t i = 0; i + 1 < path.mapping_size(); i++) {
        Mapping* here = path.mutable_mapping(i);
        if (here->edit_size() > 0) {
            Edit* last_edit = here->mutable_edit(here->edit_size() - 1);
            if (last_edit->from_length() == 0) {
                // The last edit is an insertion. So move it over to the next mapping
                Mapping* next = path.mutable_mapping(i + 1);
                // Make a space
                next->add_edit();
                // Shift everything up. Make sure to use reverse iterators to avoid clobbering.
                std::copy(next->mutable_edit()->rbegin() + 1, next->mutable_edit()->rend(), next->mutable_edit()->rbegin());
                // Put it in place
                *next->mutable_edit(0) = *last_edit;
                // Merge if needed
                *next = merge_adjacent_edits(*next);
                // Get rid of it on the Mapping we got it from
                here->mutable_edit()->RemoveLast();
            }
        }
    }
}

std::ostream& operator<<(std::ostream& out, const alignment_plan_t& plan) {
    for (auto& fate : plan) {
        switch (fate) {
        case GraphBaseFate::Match:
            out << "|";
            break;
        case GraphBaseFate::Mismatch:
            out << "*";
            break;
        case GraphBaseFate::Deletion:
            out << "-";
            break;
        case GraphBaseFate::InsertionMatch:
            out << "(";
            break;
        case GraphBaseFate::InsertionMismatch:
            out << "<";
            break;
        case GraphBaseFate::MatchInsertion:
            out << ")";
            break;
        case GraphBaseFate::MismatchInsertion:
            out << ">";
            break;
        default:
            throw std::runtime_error("Bad fate");    
        }
    }
    return out;
}


/// Make an alignment plan into a Path
std::pair<Path, std::string> plan_to_path(const alignment_plan_t& plan, const HandleGraph& graph, const vector<handle_t>& base_path, size_t start_offset) {
    
#ifdef debug_conversion_tests
    std::cerr << "Plan: " << plan << "@" << start_offset << std::endl;
#endif
    
    std::pair<Path, std::string> to_return;
    Path& path = to_return.first;
    std::string& sequence = to_return.second;
    
    // What node are we on
    auto current_node = base_path.begin();
    // How far are we done through on the current node?
    size_t node_cursor = start_offset;
    // How far are we done through on the path?
    size_t plan_cursor = 0;
    // What mapping are we working on?
    Mapping* mapping = path.add_mapping();
    // Set up its initial position
    mapping->mutable_position()->set_node_id(graph.get_id(*current_node));
    mapping->mutable_position()->set_is_reverse(graph.get_is_reverse(*current_node));
    mapping->mutable_position()->set_offset(start_offset);
    // Assemble a simulated alignment sequence
    std::stringstream ss;
    
    
    while (plan_cursor < plan.size()) {
        GraphBaseFate fate = plan[plan_cursor];
        
        // Handle this graph base.
        // Just add individual edits and merge them up later.
        
        if (fate == GraphBaseFate::InsertionMatch || fate == GraphBaseFate::InsertionMismatch) {
            // Put an insertion first
            Edit* before = mapping->add_edit();
            before->set_to_length(1);
            before->set_sequence("N");
            ss << "N";
        }
        // Put the main edit that consumes a graph base
        Edit* here = mapping->add_edit();
        here->set_from_length(1);
        if (fate == GraphBaseFate::Mismatch || fate == GraphBaseFate::InsertionMismatch || fate == GraphBaseFate::MismatchInsertion) {
            // Needs a sequence
            here->set_sequence("N");
            ss << "N";
        }
        if (fate == GraphBaseFate::Match || fate == GraphBaseFate::InsertionMatch || fate == GraphBaseFate::MatchInsertion) {
            // Keep the real base.
            // TODO: Use faster base accessor
            ss << graph.get_sequence(*current_node)[node_cursor];
        }
        if (fate != GraphBaseFate::Deletion) {
            // The graph base still exists in the sequence
            here->set_to_length(1);
        }
        if (fate == GraphBaseFate::MatchInsertion || fate == GraphBaseFate::MismatchInsertion) {
            // Put an insertion last
            Edit* after = mapping->add_edit();
            after->set_to_length(1);
            after->set_sequence("N");
            ss << "N";
        }
        
        // Advance cursors
        ++node_cursor;
        ++plan_cursor;
        
        if (node_cursor == graph.get_length(*current_node)) {
            // Advance node
            ++current_node;
            node_cursor = 0;
            if (current_node != base_path.end()) {
                // Merge all the individual edits
                *mapping = merge_adjacent_edits(*mapping);
                // Make a mapping for the next node
                mapping = path.add_mapping();
                mapping->mutable_position()->set_node_id(graph.get_id(*current_node));
                mapping->mutable_position()->set_is_reverse(graph.get_is_reverse(*current_node));
            }
        }
        
    }
    
    // Merge all the adjacent edits in the final mapping
    *mapping = merge_adjacent_edits(*mapping);
    
    // Move all the insertions to the rightmost Mappings they can be in
    send_insertions_right(path);
    
    // Remember to send the sequence
    sequence = ss.str();
    
    return to_return;
}

/// Call the given function with random elaborations of the given path of handles into an alignment  path.
void for_each_random_alignment(const HandleGraph& graph, const vector<handle_t>& base_path, const std::function<void(const Path&, const std::string&)> iteratee) {

    if (base_path.empty()) {
        return;
    }
    
    // Get an RNG
    default_random_engine generator(test_seed_source());
    
    // How many bases of the path could be on the first handle?
    size_t first_length = graph.get_length(base_path.front());
    // How many bases of the path aren't on the first handle?
    size_t middle_length = 0;
    for (auto it = base_path.begin() + 1; it != base_path.end() && (it + 1) != base_path.end(); ++it) {
        middle_length += graph.get_length(*it);
    }
    // How many bases can we leave off from the end of the last handle?
    size_t last_length = graph.get_length(base_path.back());

    for (size_t start_offset = 0; start_offset + 1 < first_length; start_offset++) {
        // For each position we could start at on the first handle
        
        for (size_t last_omitted = 0; last_omitted + 1 < last_length && (base_path.size() > 1 || last_omitted + start_offset < first_length); last_omitted++) {
            // For each number of bases we could leave out of the alignment on the last step (which may also be the first step)
        
            // Get the number of graph bases we should use.
            // All the ones on the first node that are past the offset, plus
            // those in the middle, plus those on the last node that aren't
            // omitted. Except when the first and last nodes are the same we
            // don't double-credit the bases we keep.
            size_t graph_length = first_length - start_offset + middle_length + (base_path.size() > 1 ? last_length : (size_t)0) - last_omitted;
            
            for (size_t i = 0; i < 100; i++) {
                // Generate some random plans at this length
                auto plan = make_random_plan(graph_length, generator);
            
                // Try the alignment from each plan
                auto path_and_sequence = plan_to_path(plan, graph, base_path, start_offset);
                REQUIRE(path_to_length(path_and_sequence.first) == path_and_sequence.second.size());
                iteratee(path_and_sequence.first, path_and_sequence.second);
            }
        }
    }

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

/// Compose some random alignments against the given base path through the
/// given graph and ensure they can be round-tripped from Path to WFAAlignment
/// and back.
void round_trip_versions_of(const std::vector<handle_t>& base_path, const HandleGraph& graph) {
    // For each basic route we want to look at through the graph 
    for_each_random_alignment(graph, base_path, [&](const Path& truth_path, const std::string& truth_sequence) {
        // Consider some alignments and round-trip to WFAAlignment and back.
#ifdef debug_conversion_tests
        std::cerr << "Consider sequence " << truth_sequence << " and Path " << pb2json(truth_path) << std::endl;
#endif
        WFAAlignment converted = path_to_wfa_alignment(truth_path, graph);
#ifdef debug_conversion_tests
        std::cerr << "Becomes WFAAlignment: ";
        converted.print(graph, std::cerr);
        std::cerr << std::endl;
#endif
        Path converted_back = converted.to_path(graph, truth_sequence);
#ifdef debug_conversion_tests
        std::cerr << "Converts back to Path " << pb2json(converted_back) << std::endl;
#endif
        paths_match(converted_back, truth_path);
     });
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

void full_length_match(const std::vector<std::pair<pos_t, size_t>>& seeds, const std::string& read, const std::vector<std::pair<pos_t, std::string>>& correct_alignment, const GaplessExtender& extender, size_t error_bound, bool check_seeds) {
    GaplessExtender::cluster_type cluster;
    for (auto seed : seeds) {
        cluster.insert(GaplessExtender::to_seed(seed.first, seed.second));
    }
    auto result = extender.extend(cluster, read, nullptr, error_bound);

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

void full_length_matches(const std::vector<std::pair<pos_t, size_t>>& seeds, const std::string& read, const std::vector<std::vector<std::pair<pos_t, std::string>>>& correct_alignments, const GaplessExtender& extender, size_t error_bound, double overlap_threshold) {
    GaplessExtender::cluster_type cluster;
    for (auto seed : seeds) {
        cluster.insert(GaplessExtender::to_seed(seed.first, seed.second));
    }
    auto result = extender.extend(cluster, read, nullptr, error_bound, overlap_threshold);

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
    auto result = extender.extend(cluster, read, nullptr, error_bound);

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

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Gapless extensions report correct positions", "[gapless_extender]") {

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph = build_gbwt_graph(gbwt_index);
 
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

TEST_CASE("Overlap detection for gapless extensions", "[gapless_extender]") {

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph = build_gbwt_graph(gbwt_index);
 
    SECTION("unrelated extensions") {
        GaplessExtension a {
            {
                gbwt_graph.get_handle(1, false),
                gbwt_graph.get_handle(4, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 4 }, { },
            0, false, false,
            false, false, 0, 0
        };
        GaplessExtension b {
            {
                gbwt_graph.get_handle(5, false),
                gbwt_graph.get_handle(6, false),
                gbwt_graph.get_handle(8, false),
                gbwt_graph.get_handle(9, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 4 }, { },
            0, false, false,
            false, false, 0, 0
        };
        size_t expected_overlap = 0;
        REQUIRE(a.overlap(gbwt_graph, b) == expected_overlap);
    }

    SECTION("identical extensions") {
        GaplessExtension a {
            {
                gbwt_graph.get_handle(1, false),
                gbwt_graph.get_handle(4, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 4 }, { },
            0, false, false,
            false, false, 0, 0
        };
        GaplessExtension b = a;
        size_t expected_overlap = a.length();
        REQUIRE(a.overlap(gbwt_graph, b) == expected_overlap);
    }

    SECTION("one difference") {
        GaplessExtension a {
            {
                gbwt_graph.get_handle(5, false),
                gbwt_graph.get_handle(6, false),
                gbwt_graph.get_handle(7, false),
                gbwt_graph.get_handle(9, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 4 }, { },
            0, false, false,
            false, false, 0, 0
        };
        GaplessExtension b {
            {
                gbwt_graph.get_handle(5, false),
                gbwt_graph.get_handle(6, false),
                gbwt_graph.get_handle(8, false),
                gbwt_graph.get_handle(9, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 4 }, { },
            0, false, false,
            false, false, 0, 0
        };
        REQUIRE(a.overlap(gbwt_graph, b) == a.length() - 1);
        size_t expected_overlap = a.length() - 1;
    }

    SECTION("partial overlap") {
        GaplessExtension a {
            {
                gbwt_graph.get_handle(4, false),
                gbwt_graph.get_handle(5, false),
                gbwt_graph.get_handle(6, false),
                gbwt_graph.get_handle(8, false)
            },
            2, gbwt::BidirectionalState(),
            { 0, 4 }, { },
            0, false, false,
            false, false, 0, 0
        };
        GaplessExtension b {
            {
                gbwt_graph.get_handle(5, false),
                gbwt_graph.get_handle(6, false),
                gbwt_graph.get_handle(8, false),
                gbwt_graph.get_handle(9, false)
            },
            0, gbwt::BidirectionalState(),
            { 1, 5 }, { },
            0, false, false,
            false, false, 0, 0
        };
        size_t expected_overlap = a.length() - 1;
        REQUIRE(a.overlap(gbwt_graph, b) == expected_overlap);
    }

    SECTION("shifted by one") {
        GaplessExtension a {
            {
                gbwt_graph.get_handle(4, false),
                gbwt_graph.get_handle(5, false),
                gbwt_graph.get_handle(6, false),
                gbwt_graph.get_handle(8, false),
                gbwt_graph.get_handle(9, false)
            },
            2, gbwt::BidirectionalState(),
            { 0, 5 }, { },
            0, false, false,
            false, false, 0, 0
        };
        GaplessExtension b {
            {
                gbwt_graph.get_handle(5, false),
                gbwt_graph.get_handle(6, false),
                gbwt_graph.get_handle(8, false),
                gbwt_graph.get_handle(9, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 4 }, { },
            0, false, false,
            false, false, 0, 0
        };
        size_t expected_overlap = 0;
        REQUIRE(a.overlap(gbwt_graph, b) == expected_overlap);
    }

    SECTION("paths of different lengths") {
        GaplessExtension a {
            {
                gbwt_graph.get_handle(1, false),
                gbwt_graph.get_handle(2, false),
                gbwt_graph.get_handle(4, false),
                gbwt_graph.get_handle(5, false)
            },
            0, gbwt::BidirectionalState(),
            { 0, 6 }, { },
            0, false, false,
            false, false, 0, 0
        };
        GaplessExtension b {
            {
                gbwt_graph.get_handle(1, false),
                gbwt_graph.get_handle(4, false),
                gbwt_graph.get_handle(5, false)
            },
            0, gbwt::BidirectionalState(),
            { 1, 6 }, { },
            0, false, false,
            false, false, 0, 0
        };
        size_t expected_overlap = b.length() - 1;
        REQUIRE(a.overlap(gbwt_graph, b) == expected_overlap);
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

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph = build_gbwt_graph(gbwt_index);

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

    // We also test that we avoid finding the same best alignment from multiple seeds.
    SECTION("secondary alignment has more mismatches") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(2, false, 0), 1 }, // First seed for the best alignment (1 mismatch).
            { make_pos_t(4, false, 0), 2 }, // Second seed for the best alignment (1 mismatch).
            { make_pos_t(4, false, 0), 1 }  // Seed for the secondary alignment (2 mismatches).
        };
        std::string read = "GAGGA";
        size_t error_bound = 2;
        double overlap_threshold = 0.9;
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
        full_length_matches(seeds, read, correct_alignments, extender, error_bound, overlap_threshold);
    }

    SECTION("no secondary alignment found if the overlap is too high") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(2, false, 0), 1 }, // First seed for the best alignment (1 mismatch).
            { make_pos_t(4, false, 0), 2 }, // Second seed for the best alignment (1 mismatch).
            { make_pos_t(4, false, 0), 1 }  // Seed for the secondary alignment (2 mismatches).
        };
        std::string read = "GAGGA";
        size_t error_bound = 2;
        double overlap_threshold = 0.1;
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_alignments {
            {
                { make_pos_t(1, false, 0), "1" },
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "2A" }
            }
        };
        full_length_matches(seeds, read, correct_alignments, extender, error_bound, overlap_threshold);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Local alignments", "[gapless_extender]") {

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph = build_gbwt_graph(gbwt_index);

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

    SECTION("trim left flank") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 2), 4 }
        };
        std::string read = "xAGxGTAx";
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_extensions {
            {
                { make_pos_t(4, false, 2), "1" },
                { make_pos_t(5, false, 0), "1" },
                { make_pos_t(6, false, 0), "1" }
            }
        };
        std::vector<size_t> correct_offsets {
            static_cast<size_t>(4)
        };
        size_t error_bound = 0;
        partial_matches(seeds, read, correct_extensions, correct_offsets, extender, error_bound);
    }

    SECTION("trim right flank") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(4, false, 2), 4 }
        };
        std::string read = "xAGGGTxAx";
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_extensions {
            {
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "3" },
                { make_pos_t(5, false, 0), "1" }
            }
        };
        std::vector<size_t> correct_offsets {
            static_cast<size_t>(1)
        };
        size_t error_bound = 1;
        partial_matches(seeds, read, correct_extensions, correct_offsets, extender, error_bound);
    }

    SECTION("remove duplicates") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(2, false, 0), 1 },
            { make_pos_t(4, false, 2), 4 }
        };
        std::string read = "xAGGGTxAx";
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_extensions {
            {
                { make_pos_t(2, false, 0), "1" },
                { make_pos_t(4, false, 0), "3" },
                { make_pos_t(5, false, 0), "1" }
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

TEST_CASE("Non-ACGT characters do not match", "[gapless_extender]") {

    // Create a single-node GBWTGraph.
    bdsg::HashGraph graph;
    graph.create_handle("NNNGATTACANNN", 1);
    std::vector<gbwt::vector_type> paths = {
        { static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(1, false)) }
    };
    gbwt::GBWT gbwt_index = get_gbwt(paths);
    gbwtgraph::GBWTGraph gbwt_graph(gbwt_index, graph);


    // Wrap it in a GaplessExtender with an Aligner.
    Aligner aligner;
    GaplessExtender extender(gbwt_graph, aligner);

    SECTION("exact matching") {
        std::vector<std::pair<pos_t, size_t>> seeds {
            { make_pos_t(1, false, 5), 4 }
        };
        std::string read = "NNGATTACANN";
        std::vector<std::vector<std::pair<pos_t, std::string>>> correct_extensions {
            {
                { make_pos_t(1, false, 3), "7" }
            }
        };
        std::vector<size_t> correct_offsets {
            static_cast<size_t>(2)
        };
        size_t error_bound = 0;
        partial_matches(seeds, read, correct_extensions, correct_offsets, extender, error_bound);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Gapless extensions can be converted to WFAAlignments and joined", "[wfa_alignment]") {

    // Build a GBWT with three threads including a duplicate.
    gbwt::GBWT gbwt_index = build_gbwt_index();

    // Build a GBWT-backed graph.
    gbwtgraph::GBWTGraph gbwt_graph = build_gbwt_graph(gbwt_index);
    
    // And we need an Aligner.
    Aligner aligner;
    
    // Make an extension
    GaplessExtension a {
        {
            gbwt_graph.get_handle(1, false),
            gbwt_graph.get_handle(4, false)
        },
        0, gbwt::BidirectionalState(),
        { 0, 4 }, { },
        4 * aligner.match, false, false,
        false, false, 0, 0
    };
    correct_score(a, aligner);
    // And turn it into a WFAAlignment
    WFAAlignment a_aln = WFAAlignment::from_extension(a);
    correct_score(a_aln, aligner);
    
    // Make another abutting extension
    GaplessExtension b {
        {
            gbwt_graph.get_handle(5, false),
            gbwt_graph.get_handle(6, false),
            gbwt_graph.get_handle(8, false),
            gbwt_graph.get_handle(9, false)
        },
        0, gbwt::BidirectionalState(),
        { 4, 8 }, { },
        4 * aligner.match, false, false,
        false, false, 0, 0
    };
    // And turn it into a WFAAlignment
    correct_score(b, aligner);
    WFAAlignment b_aln = WFAAlignment::from_extension(b);
    correct_score(b_aln, aligner);
    
    // Now weld them together
    WFAAlignment joined = WFAAlignment::make_empty();
    joined.join(a_aln);
    REQUIRE(joined.length == a_aln.length);
    correct_score(joined, aligner);
    
    joined.join(b_aln);
    REQUIRE(joined.length == a_aln.length + b_aln.length);
    correct_score(joined, aligner);
    
    // And make sure we get the right final alignment
    REQUIRE(joined.edits.size() == 1);
    REQUIRE(joined.edits.at(0).first == WFAAlignment::match);
    REQUIRE(joined.edits.at(0).second == 8);
    REQUIRE(joined.path.size() == 6);
    REQUIRE(joined.path.at(0) == gbwt_graph.get_handle(1, false));
    REQUIRE(joined.path.at(1) == gbwt_graph.get_handle(4, false));
    REQUIRE(joined.path.at(2) == gbwt_graph.get_handle(5, false));
    REQUIRE(joined.path.at(3) == gbwt_graph.get_handle(6, false));
    REQUIRE(joined.path.at(4) == gbwt_graph.get_handle(8, false));
    REQUIRE(joined.path.at(5) == gbwt_graph.get_handle(9, false));
    REQUIRE(joined.score == 8 * aligner.match);
}

namespace {

gbwt::GBWT wfa_linear_gbwt() {
    std::vector<gbwt::vector_type> paths;
    paths.emplace_back();
    paths.back().push_back(gbwt::Node::encode(1, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(3, false));
    paths.back().push_back(gbwt::Node::encode(4, false));
    return get_gbwt(paths);
}

gbwtgraph::GBWTGraph wfa_linear_graph(const gbwt::GBWT& index) {
    gbwtgraph::SequenceSource source;
    source.add_node(1, "CGC");
    source.add_node(2, "GATTACA");
    source.add_node(3, "GATTA");
    source.add_node(4, "TAT");
    return gbwtgraph::GBWTGraph(index, source);
}

gbwt::GBWT wfa_general_gbwt() {
    std::vector<gbwt::vector_type> paths;

    // Path that makes the first choice.
    paths.emplace_back();
    paths.back().push_back(gbwt::Node::encode(1, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(3, false));
    paths.back().push_back(gbwt::Node::encode(5, false));
    paths.back().push_back(gbwt::Node::encode(6, false));
    paths.back().push_back(gbwt::Node::encode(8, false));
    paths.back().push_back(gbwt::Node::encode(9, false));
    paths.back().push_back(gbwt::Node::encode(11, false));

    // Path that makes the second choice.
    paths.emplace_back();
    paths.back().push_back(gbwt::Node::encode(1, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(4, false));
    paths.back().push_back(gbwt::Node::encode(5, false));
    paths.back().push_back(gbwt::Node::encode(7, false));
    paths.back().push_back(gbwt::Node::encode(8, false));
    paths.back().push_back(gbwt::Node::encode(10, false));
    paths.back().push_back(gbwt::Node::encode(11, false));

    return get_gbwt(paths);
}

gbwtgraph::GBWTGraph wfa_general_graph(const gbwt::GBWT& index) {
    gbwtgraph::SequenceSource source;

    // Start.
    source.add_node(1, "CGC");
    source.add_node(2, "GATTACA");

    // Simple bubble.
    source.add_node(3, "G");
    source.add_node(4, "C");

    // Middle.
    source.add_node(5, "ATTA");

    // Nondeterministic bubble.
    source.add_node(6, "TG");
    source.add_node(7, "TC");

    // Middle.
    source.add_node(8, "GAA");

    // Bubble that depends on the nondeterministic choice.
    source.add_node(9, "CAT");
    source.add_node(10, "GTA");

    // End.
    source.add_node(11, "TAT");

    return gbwtgraph::GBWTGraph(index, source);
}
gbwt::GBWT wfa_cycle_gbwt() {
    std::vector<gbwt::vector_type> paths;

    // Path that skips the cycle
    paths.emplace_back();
    paths.back().push_back(gbwt::Node::encode(1, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(3, false));

    // Path that takes the cycle once.
    paths.emplace_back();
    paths.back().push_back(gbwt::Node::encode(1, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(3, false));

    // Path that takes the cycle twice.
    paths.emplace_back();
    paths.back().push_back(gbwt::Node::encode(1, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(2, false));
    paths.back().push_back(gbwt::Node::encode(3, false));

    return get_gbwt(paths);
}

gbwtgraph::GBWTGraph wfa_cycle_graph(const gbwt::GBWT& index) {
    gbwtgraph::SequenceSource source;
    source.add_node(1, "CGC");
    source.add_node(2, "GA");
    source.add_node(3, "TAT");
    return gbwtgraph::GBWTGraph(index, source);
}

void check_score(const WFAAlignment& alignment, const Aligner& aligner, int32_t matches, int32_t mismatches, int32_t gaps, int32_t gap_length) {
    int32_t extensions = gap_length - gaps;
    int32_t expected_score = matches * aligner.match - mismatches * aligner.mismatch - gaps * aligner.gap_open - extensions * aligner.gap_extension;
    REQUIRE(alignment.score == expected_score);
}

void check_unlocalized_insertion(const WFAAlignment& alignment, const std::string& sequence, const Aligner& aligner) {
    REQUIRE(alignment.unlocalized_insertion());

    REQUIRE(alignment.path.empty());

    REQUIRE(alignment.edits.size() == 1);
    REQUIRE(alignment.edits.front().first == WFAAlignment::insertion);
    REQUIRE(alignment.edits.front().second == sequence.length());

    REQUIRE(alignment.seq_offset == 0);
    REQUIRE(alignment.length == sequence.length());

    int32_t gap_score = -aligner.gap_open - (int32_t(alignment.length) - 1) * aligner.gap_extension;
    REQUIRE(alignment.score == gap_score);
}

void check_alignment(const WFAAlignment& alignment, const std::string& sequence, const gbwtgraph::GBWTGraph& graph, const Aligner& aligner, const pos_t* from, const pos_t* to) {

    REQUIRE(alignment);

    // Sequence range is sane and corresponds to the specified endpoints.
    REQUIRE(alignment.seq_offset + alignment.length <= sequence.length());
    REQUIRE((from == nullptr) | (alignment.seq_offset == 0));
    REQUIRE((to == nullptr) | (alignment.seq_offset + alignment.length == sequence.length()));

    // Total length of edits in the sequence.
    uint32_t edit_total = 0;
    for (auto edit : alignment.edits) {
        if (edit.first != WFAAlignment::deletion) {
            edit_total += edit.second;
        }
    }
    REQUIRE(edit_total == alignment.length);

    // Check that the path is valid.
    for (size_t i = 0; i < alignment.path.size(); i++) {
        REQUIRE(graph.has_node(graph.get_id(alignment.path[i])));
        if (i > 0) {
            REQUIRE(graph.has_edge(alignment.path[i - 1], alignment.path[i]));
        }
    }

    // Check that the alignment is between the right positions, if provided, and the start/end nodes are used in the alignment.
    if (!alignment.path.empty()) {
        REQUIRE(alignment.node_offset < graph.get_length(alignment.path.front()));
        if (from != nullptr) {
            handle_t from_handle = graph.get_handle(id(*from), is_rev(*from));
            if (alignment.path.front() == from_handle && alignment.node_offset > 0) {
                REQUIRE(alignment.node_offset == offset(*from) + 1);
            } else {
                REQUIRE(offset(*from) + 1 == graph.get_length(from_handle));
                REQUIRE(graph.has_edge(from_handle, alignment.path.front()));
                REQUIRE(alignment.node_offset == 0);
            }
        }
        uint32_t final_offset = alignment.final_offset(graph);
        REQUIRE(final_offset > 0);
        if (to != nullptr) {
            handle_t to_handle = graph.get_handle(id(*to), is_rev(*to));
            uint32_t final_node_length = graph.get_length(alignment.path.back());
            if (alignment.path.back() == to_handle && final_offset < final_node_length) {
                REQUIRE(final_offset == offset(*to));
            } else {
                REQUIRE(offset(*to) == 0);
                REQUIRE(graph.has_edge(alignment.path.back(), to_handle));
                REQUIRE(final_offset == final_node_length);
            }
        }
    }

    // Check that edits of the same type are merged.
    for (size_t i = 1; i < alignment.edits.size(); i++) {
        REQUIRE(alignment.edits[i - 1].first != alignment.edits[i].first);
    }

    // Compute the score using the parameters from the aligner.
    int32_t score_from_edits = 0;
    for (auto edit : alignment.edits) {
        switch (edit.first)
        {
        case WFAAlignment::match:
            score_from_edits += int32_t(edit.second) * aligner.match;
            break;
        case WFAAlignment::mismatch:
            score_from_edits -= int32_t(edit.second) * aligner.mismatch;
            break;
        case WFAAlignment::insertion: // Fall through.
        case WFAAlignment::deletion:
            // Note that a gap of length n has n - 1 extensions according to VG.
            score_from_edits -= aligner.gap_open + (int32_t(edit.second) - 1) * aligner.gap_extension;
            break;
        }
    }
    REQUIRE(alignment.score == score_from_edits);

    // Check the alignment itself.
    size_t seq_offset = alignment.seq_offset;
    size_t node_offset = alignment.node_offset;
    size_t path_offset = 0;
    for (auto edit : alignment.edits) {
        if (edit.first == WFAAlignment::insertion) {
            seq_offset += edit.second;
            continue;
        }
        size_t edit_end = node_offset + edit.second;
        while (edit_end > node_offset) {
            REQUIRE(path_offset < alignment.path.size());
            std::string node_sequence = graph.get_sequence(alignment.path[path_offset]);
            size_t len = std::min(edit_end, node_sequence.length()) - node_offset;
            if (edit.first == WFAAlignment::match) {
                REQUIRE(sequence.substr(seq_offset, len) == node_sequence.substr(node_offset, len));
                seq_offset += len;
            } else if (edit.first == WFAAlignment::mismatch) {
                for (size_t i = 0; i < len; i++) {
                    REQUIRE(sequence[seq_offset + i] != node_sequence[node_offset + i]);
                }
                seq_offset += len;
            }
            node_offset += len;
            if (node_offset >= node_sequence.length()) {
                node_offset = 0;
                edit_end -= node_sequence.length();
                path_offset++;
            }
        }
    }
    if (!alignment.path.empty()) {
        REQUIRE((path_offset == alignment.path.size() - 1) | (path_offset == alignment.path.size() && node_offset == 0));
    }
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Exact matches in a linear graph", "[wfa_extender]") {
    // Create the structures for graph 1: CGC, 2: GATTACA, 3: GATTA, 4: TAT
    gbwt::GBWT index = wfa_linear_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_linear_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Single node, start to end") {
        std::string sequence("GATTACA");
        pos_t from(1, false, 2); pos_t to(3, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Single node, middle") {
        std::string sequence("ATTAC");
        pos_t from(2, false, 0); pos_t to(2, false, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Multiple nodes, start to end") {
        std::string sequence("GATTACAGATTA");
        pos_t from(1, false, 2); pos_t to(4, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Multiple nodes, middle") {
        std::string sequence("ATTACAGATT");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Multiple nodes, end to start") {
        std::string sequence("AG");
        pos_t from(2, false, 5); pos_t to(3, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Multiple nodes, reverse, start to end") {
        std::string sequence("TAATCTGTAATC");
        pos_t from(4, true, 2); pos_t to(1, true, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Multiple nodes, reverse, middle") {
        std::string sequence("AATCTGTAAT");
        pos_t from(3, true, 0); pos_t to(2, true, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Multiple nodes, reverse, end to start") {
        std::string sequence("CT");
        pos_t from(3, true, 3); pos_t to(2, true, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Mismatches in a linear graph", "[wfa_extender]") {
    // Create the structures for graph 1: CGC, 2: GATTACA, 3: GATTA, 4: TAT
    gbwt::GBWT index = wfa_linear_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_linear_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("In the middle") {
        // MMMXMM|MMMM
        std::string sequence("ATTCCAGATT");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("In the middle, reverse") {
        // MMMM|MMMXMM
        std::string sequence("AATCTGTTAT");
        pos_t from(3, true, 0); pos_t to(2, true, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("At both ends") {
        // XMMMMM|MMMX
        std::string sequence("TTTACAGATA");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 2, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("At both ends, reverse") {
        // XMMM|MMMMMX
        std::string sequence("TATCTGTAAA");
        pos_t from(3, true, 0); pos_t to(2, true, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 2, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Over node boundary") {
        // MMMMMX|XMMM
        std::string sequence("ATTACTTATT");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 2, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Over node boundary, reverse") {
        // MMMX|XMMMMM
        std::string sequence("AATAAGTAAT");
        pos_t from(3, true, 0); pos_t to(2, true, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 2, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Gaps in a linear graph", "[wfa_extender]") {
    // Create the structures for graph 1: CGC, 2: GATTACA, 3: GATTA, 4: TAT
    gbwt::GBWT index = wfa_linear_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_linear_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Deletion in the middle") {
        // MDDMMM|MMMM
        std::string sequence("AACAGATT");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Deletion in the middle, reverse") {
        // MMMM|MMMDDMM
        std::string sequence("AATCTGTT");
        pos_t from(3, true, 0); pos_t to(2, true, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Insertion in the middle") {
        // MMMMMIIM|MMMM
        std::string sequence("ATTATACAGATT");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Insertion in the middle, reverse") {
        // MMMM|MIIMMMMM
        std::string sequence("AATCTCCGTAAT");
        pos_t from(3, true, 0); pos_t to(2, true, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Deletion over node boundary") {
        // MMMMMD|DDMM
        std::string sequence("ATTACTT");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 1, 3);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Deletion over node boundary, reverse") {
        // MMMD|DMMMMM
        std::string sequence("AATGTAAT");
        pos_t from(3, true, 0); pos_t to(2, true, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Insertion at node boundary") {
        // MMMMMMII|MMMM
        std::string sequence("ATTACATTGATT");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Insertion at node boundary, reverse") {
        // MMMMII|MMMMMM
        std::string sequence("AATCAATGTAAT");
        pos_t from(3, true, 0); pos_t to(2, true, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Deletion at the end") {
        // MMMMMM|MMMD
        std::string sequence("ATTACAGAT");
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 1, 1);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Deletion at the end, reverse") {
        // MMMM|MMMMD
        std::string sequence("AATCTGTA");
        pos_t from(3, true, 0); pos_t to(2, true, 5);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 1, 1);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Insertion at the end") {
        // MMMMMM|MMII
        std::string sequence("ATTACAGATT");
        pos_t from(2, false, 0); pos_t to(3, false, 2);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Insertion at the end, reverse") {
        // MMMM|MMMMII
        std::string sequence("AATCTGTAAT");
        pos_t from(3, true, 0); pos_t to(2, true, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Deletion at both ends") {
        // DMMMMMM|MMMMD
        std::string sequence("ATTACAGATT");
        pos_t from(1, false, 2); pos_t to(4, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 2, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Insertion at both ends") {
        // IIMM|MMMMMM|MMMMM|MMII
        std::string sequence("AAGCGATTACAGATTATACC");
        pos_t from(1, false, 0); pos_t to(4, false, 2);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 4, 0, 2, 4);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Special cases in a linear graph", "[wfa_extender]") {
    // Create the structures for graph 1: CGC, 2: GATTACA, 3: GATTA, 4: TAT
    gbwt::GBWT index = wfa_linear_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_linear_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Empty sequence, deletion") {
        // DD
        std::string sequence;
        pos_t from(2, false, 0); pos_t to(2, false, 3);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Empty sequence, failure") {
        // DDDDDD|DDDD
        std::string sequence;
        pos_t from(2, false, 0); pos_t to(3, false, 4);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }

    SECTION("Identical endpoints") {
        std::string sequence;
        pos_t from(2, false, 3); pos_t to(2, false, 3);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }

    SECTION("Adjacent endpoints, empty sequence") {
        std::string sequence;
        pos_t from(2, false, 6); pos_t to(3, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Adjacent endpoints, insertion") {
        std::string sequence("AA");
        pos_t from(2, false, 6); pos_t to(3, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Cannot align within the score bound") {
        // MMXXMXM
        std::string sequence("GAGGAGA");
        pos_t from(1, false, 2); pos_t to(3, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }

    SECTION("Cannot reach target") {
        std::string sequence("GATTA");
        pos_t from(3, false, 0); pos_t to(2, false, 6);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }

    SECTION("Run out of graph") {
        // MM|----
        std::string sequence("ATCCCC");
        pos_t from(4, false, 0); pos_t to(5, false, 3);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }

    SECTION("Mixed edits") {
        // MMMMXM|MMDDM|MM
        std::string sequence("ATTAGAGAATA");
        pos_t from(2, false, 0); pos_t to(4, false, 2);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 1, 1, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Mismatches beat ins + del") {
        // MMXXMMM
        std::string sequence("GACCACA");
        pos_t from(1, false, 2); pos_t to(3, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 2, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Prefixes in a linear graph", "[wfa_extender]") {
    // Create the structures for graph 1: CGC, 2: GATTACA, 3: GATTA, 4: TAT
    gbwt::GBWT index = wfa_linear_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_linear_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Empty sequence") {
        std::string sequence;
        pos_t to(2, false, 1);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Got an unlocalized insertion") {
        // IIII
        std::string sequence("GGGG");
        pos_t to(3, false, 0);
        WFAAlignment result = extender.prefix(sequence, to);
        check_unlocalized_insertion(result, sequence, aligner);
    }

    SECTION("Cannot align within the score bound") {
        // IIIII
        std::string sequence("GGGGGG");
        pos_t to(3, false, 0);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, 0, 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Exact match, middle of a node") {
        std::string sequence("ATTAC");
        pos_t to(2, false, 6);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Exact match, start of a node") {
        std::string sequence("GATTA");
        pos_t to(2, false, 5);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Exact match, end of a node") {
        std::string sequence("CGATTA");
        pos_t to(2, false, 5);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Exact match, reverse") {
        std::string sequence("TAAT");
        pos_t to(2, true, 6);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("With edits") {
        // MMM|MMMMDDM|MMMXM
        std::string sequence("CGCGATTAGATAA");
        pos_t to(4, false, 0);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 1, 1, 1, 2);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Trim the prefix") {
        // ------|MMMMM
        std::string sequence("TTTTTTGATTA");
        pos_t to(4, false, 0);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 6, 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Run out of graph, trim") {
        // ----|MMM|MMMXM
        std::string sequence("ATATCGCGATAA");
        pos_t to(2, false, 5);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 5, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Run out of graph, no trim") {
        // I|MMM|MMMXMMM
        std::string sequence("TCGCGATAACA");
        pos_t to(3, false, 0);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 2, 1, 1, 1);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Suffixes in a linear graph", "[wfa_extender]") {
    // Create the structures for graph 1: CGC, 2: GATTACA, 3: GATTA, 4: TAT
    gbwt::GBWT index = wfa_linear_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_linear_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Empty sequence") {
        std::string sequence;
        pos_t from(2, false, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Got an unlocalized insertion") {
        // IIII
        std::string sequence("GGGG");
        pos_t from(2, false, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_unlocalized_insertion(result, sequence, aligner);
    }

    SECTION("Cannot align within the score bound") {
        // IIIII
        std::string sequence("GGGGGG");
        pos_t from(2, false, 0);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, 0, 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Exact match, middle of a node") {
        std::string sequence("ATTAC");
        pos_t from(2, false, 0);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Exact match, start of a node") {
        std::string sequence("TTACAG");
        pos_t from(2, false, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Exact match, end of a node") {
        std::string sequence("TTACA");
        pos_t from(2, false, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Exact match, reverse") {
        std::string sequence("TAAT");
        pos_t from(2, true, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("With edits") {
        // MMXMMM|MDDMM|MM
        std::string sequence("ATGACAGTATA");
        pos_t from(2, false, 0);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 1, 1, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Trim the suffix") {
        // MMMMM|------
        std::string sequence("GATTAAAAAAA");
        pos_t from(2, false, 6);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 6, 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Run out of graph, trim") {
        // MXMM|MMM|---
        std::string sequence("AATATATCAC");
        pos_t from(3, false, 0);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 4, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Run out of graph, no trim") {
        // MMMMMM|MMMXM|MMM|II
        std::string sequence("ATTACAGATCATATCC");
        pos_t from(2, false, 0);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 3, 1, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Connect in a general graph", "[wfa_extender]") {
    // 1   2         5       8       11
    // CGC|GATTACA|G|ATTA|TG|GAA|CAT|TAT
    // CGC|GATTACA|C|ATTA|TC|GAA|GTA|TAT
    gbwt::GBWT index = wfa_general_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_general_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Exact match") {
        std::string sequence("ACAGATTATGG");
        pos_t from(2, false, 3); pos_t to(8, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Mismatch") {
        // One of these:         X     X
        std::string sequence("ACACATTATGG");
        pos_t from(2, false, 3); pos_t to(8, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Mix of edits") {
        // MMMM|D|MMMM|XM|MMM
        std::string sequence("TACAATTACCGAA");
        pos_t from(2, false, 2); pos_t to(10, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 1, 1, 1, 1);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Nondeterministic choice") {
        // MMM|MM|MXM|MMM|MM, bottom choice
        std::string sequence("TTATCGTAGTATA");
        pos_t from(5, false, 0); pos_t to(11, false, 2);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("No haplotype connection") {
        // Matches 6 -> 8 -> 10
        std::string sequence("GGAAGT");
        pos_t from(6, false, 0); pos_t to(10, false, 2);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }

    SECTION("No haplotype connection that includes the endpoints") {
        // Matches node 8.
        std::string sequence("GAA");
        pos_t from(6, false, 1); pos_t to(10, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }

    SECTION("Not within score bound") {
        // M|XMMMMMM|X|MXXM
        std::string sequence("CCATTACATAGGA");
        pos_t from(1, false, 1); pos_t to(6, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Prefix in a general graph", "[wfa_extender]") {
    // 1   2         5       8       11
    // CGC|GATTACA|G|ATTA|TG|GAA|CAT|TAT
    // CGC|GATTACA|C|ATTA|TC|GAA|GTA|TAT
    gbwt::GBWT index = wfa_general_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_general_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Exact match") {
        std::string sequence("TACACAT");
        pos_t to(5, false, 2);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Mismatch") {
        // MMXM|M|MMM
        std::string sequence("TAGAGATT");
        pos_t to(5, false, 3);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Mix of edits") {
        // M|MMMMMXM|M|MDDM|MM|MMM
        std::string sequence("CGATTAGACAATCGAA");
        pos_t to(10, false, 0);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 1, 1, 1, 2);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Nondeterministic choice") {
        // MMM|MXM|MM|MMM, bottom options on the reverse strand
        std::string sequence("TACTACGATAA");
        pos_t to(5, true, 3);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Trim the prefix") {
        // ------|M|MMMM
        std::string sequence("GGGGGGGATTA");
        pos_t to(6, false, 0);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 6, 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }

    SECTION("Run out of graph, trim") {
        // ------|MMM|MMMM
        std::string sequence("AAAAAACGCGATT");
        pos_t to(2, false, 4);
        WFAAlignment result = extender.prefix(sequence, to);
        check_score(result, aligner, sequence.length() - 6, 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, nullptr, &to);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Suffix in a general graph", "[wfa_extender]") {
    // 1   2         5       8       11
    // CGC|GATTACA|G|ATTA|TG|GAA|CAT|TAT
    // CGC|GATTACA|C|ATTA|TC|GAA|GTA|TAT
    gbwt::GBWT index = wfa_general_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_general_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Exact match") {
        std::string sequence("ATTATGGAA");
        pos_t from(3, false, 0);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Mismatch") {
        // MMM|M|MMMM|MX|M
        std::string sequence("ACACATTATGG");
        pos_t from(2, false, 3);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Mix of edits") {
        // MMMXM|M|MMMD|DM|MM
        std::string sequence("TTAGACATTCGA");
        pos_t from(2, false, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 1, 1, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Nondeterministic choice") {
        // MM|MM|MMX|MMM|M, top options
        std::string sequence("TATGGATCATT");
        pos_t from(5, false, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Trim the suffix") {
        // MMM|MM|------
        std::string sequence("GAAGTCCCCCC");
        pos_t from(7, false, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 6, 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }

    SECTION("Run out of graph, trim") {
        // M|MMM|MMM|-------
        std::string sequence("AGTATATAAAAAAA");
        pos_t from(8, false, 1);
        WFAAlignment result = extender.suffix(sequence, from);
        check_score(result, aligner, sequence.length() - 7, 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, nullptr);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Connect with a cycle", "[wfa_extender]") {
    // CGC GA TAT
    // CGC GA GA TAT
    // CGC GA GA GA TAT
    gbwt::GBWT index = wfa_cycle_gbwt();
    gbwtgraph::GBWTGraph graph = wfa_cycle_graph(index);
    Aligner aligner;
    WFAExtender extender(graph, aligner);

    SECTION("Skip the cycle") {
        std::string sequence("CGAT");
        pos_t from(1, false, 1); pos_t to(3, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Take it once") {
        std::string sequence("CGAGAT");
        pos_t from(1, false, 1); pos_t to(3, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Take it twice") {
        std::string sequence("CGAGAGAT");
        pos_t from(1, false, 1); pos_t to(3, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Take it thrice, insertion") {
        std::string sequence("CGAGAGAGAT");
        pos_t from(1, false, 1); pos_t to(3, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 2, 0, 1, 2);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Take it five times, fail") {
        std::string sequence("CGAGAGAGAGAGAT");
        pos_t from(1, false, 1); pos_t to(3, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        REQUIRE_FALSE(result);
    }

    SECTION("Reach the end twice") {
        // MM|MM|M
        std::string sequence("GCGAG");
        pos_t from(1, false, 0); pos_t to(2, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length(), 0, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Identical endpoints, at start") {
        // M|MX; correct endpoint on the third visit to the node
        std::string sequence("AGG");
        pos_t from(2, false, 0); pos_t to(2, false, 0);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }

    SECTION("Identical endpoints, at end") {
        // MM|X; correct endpoint on the third visit to the node
        std::string sequence("GAT");
        pos_t from(2, false, 1); pos_t to(2, false, 1);
        WFAAlignment result = extender.connect(sequence, from, to);
        check_score(result, aligner, sequence.length() - 1, 1, 0, 0);
        check_alignment(result, sequence, graph, aligner, &from, &to);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("WFAAlignment can be converted to Path", "[wfa_alignment]") {

    SECTION("Generic graph") {
        gbwt::GBWT index = wfa_general_gbwt();
        gbwtgraph::GBWTGraph graph = wfa_general_graph(index);
        // Looks like: 1, 2, {3|4}, 5, {6|7}, 8, {9|10}, 11
        
        std::vector<std::vector<handle_t>> base_paths = {
            {graph.get_handle(1), graph.get_handle(2), graph.get_handle(3), graph.get_handle(5)},
            {graph.get_handle(11, true), graph.get_handle(9, true), graph.get_handle(8, true)}
        };
        
        for (auto& base_path : base_paths) {
            round_trip_versions_of(base_path, graph);
        }
    }

    SECTION("Cyclic graph") {
        gbwt::GBWT index = wfa_cycle_gbwt();
        gbwtgraph::GBWTGraph graph = wfa_cycle_graph(index);
        // Looks like 1, 2*, 3
        
        std::vector<std::vector<handle_t>> base_paths = {
            {graph.get_handle(1), graph.get_handle(2), graph.get_handle(3)},
            {graph.get_handle(3, true), graph.get_handle(2, true), graph.get_handle(2, true), graph.get_handle(2, true), graph.get_handle(1, true)}
        };
        
        for (auto& base_path : base_paths) {
            round_trip_versions_of(base_path, graph);
        }
    }
    
}

//------------------------------------------------------------------------------


}
}
