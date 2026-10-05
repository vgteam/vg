///
/// \file position.cpp
///  
/// Unit tests for Position and pos_t manipulation
///

#include "catch.hpp"
#include "surjector.hpp"
#include "aligner.hpp"
#include "alignment.hpp"
#include "annotation.hpp"

#include "bdsg/hash_graph.hpp"
#include "bdsg/overlays/path_position_overlays.hpp"
#include <vg/vg.pb.h>

namespace vg {
namespace unittest {

class TestSurjector : public Surjector {
public:
    TestSurjector(const PathPositionHandleGraph* graph) : Surjector(graph) {}
    ~TestSurjector() = default;
    
    using Surjector::extract_overlapping_paths;
    using Surjector::filter_redundant_path_chunks;
    using Surjector::prune_and_trim_anchors;
    using Surjector::anchor_has_nearby_repeat;
    using Surjector::choose_primary;
    using Surjector::choose_primary_strand;
    using Surjector::diploid_pair_span;
    using Surjector::score_diploid_pair;
    using Surjector::get_aligner;
    
};


// Read individual SAM fields without assuming SA is the only tag.
static map<string, pair<char, string>> parse_sam_tags_for_test(const Alignment& aln) {
    map<string, pair<char, string>> tags;
    if (has_annotation(aln, "tags")) {
        istringstream input(get_annotation<string>(aln, "tags"));
        string field;
        while (input >> field) {
            REQUIRE(field.size() >= 5);
            REQUIRE(field[2] == ':');
            REQUIRE(field[4] == ':');
            REQUIRE(tags.emplace(field.substr(0, 2),
                                 make_pair(field[3], field.substr(5))).second);
        }
    }
    return tags;
}

TEST_CASE("Diploid pair compatibility is separate from fragment scoring", "[surject][diploid][paired]") {
    bdsg::HashGraph graph;
    auto node = graph.create_handle(string(1000, 'A'));
    auto path = graph.create_path_handle("ref");
    graph.append_step(path, node);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    auto alignment = [](int64_t start, int64_t length, bool reverse) {
        Alignment aln;
        aln.set_score(30);
        auto* pos = aln.add_refpos();
        pos->set_name("ref"); pos->set_offset(start); pos->set_is_reverse(reverse);
        auto* mapping = aln.mutable_path()->add_mapping();
        mapping->mutable_position()->set_node_id(1);
        auto* edit = mapping->add_edit();
        edit->set_from_length(length); edit->set_to_length(length);
        return aln;
    };
    auto first = alignment(100, 20, false), second = alignment(180, 20, true);
    SECTION("Inward facing mates may occur in either mate order") {
        CHECK(surjector.diploid_pair_span(first, second, 100) == optional<int64_t>(100));
        CHECK(surjector.diploid_pair_span(second, first, 100) == optional<int64_t>(100));
        CHECK_FALSE(surjector.diploid_pair_span(first, second, 99));
        CHECK(surjector.diploid_pair_span(first, second, 0) == optional<int64_t>(100));
    }
    SECTION("Other paths, same orientation, and outward facing pairs are incompatible") {
        second.mutable_refpos(0)->set_name("other");
        CHECK_FALSE(surjector.diploid_pair_span(first, second, 0));
        second.mutable_refpos(0)->set_name("ref");
        second.mutable_refpos(0)->set_is_reverse(false);
        CHECK_FALSE(surjector.diploid_pair_span(first, second, 0));
        first.mutable_refpos(0)->set_is_reverse(true);
        CHECK_FALSE(surjector.diploid_pair_span(first, second, 0));
    }
    SECTION("Overlapping and contained mates are allowed") {
        CHECK(surjector.diploid_pair_span(alignment(100, 60, false), alignment(140, 60, true), 100) == optional<int64_t>(100));
        CHECK(surjector.diploid_pair_span(alignment(110, 20, false), alignment(100, 60, true), 60) == optional<int64_t>(60));
    }
    SECTION("Unmapped or missing positions are not valid pair candidates") {
        second.clear_path();
        CHECK_FALSE(surjector.diploid_pair_span(first, second, 0));
        first.clear_refpos();
        CHECK_FALSE(surjector.diploid_pair_span(first, second, 0));
    }
    SECTION("Reference span includes deletions and excludes inserted or clipped bases") {
        auto* deletion = second.mutable_path()->mutable_mapping(0)->add_edit();
        deletion->set_from_length(10);
        auto* match = second.mutable_path()->mutable_mapping(0)->add_edit();
        match->set_from_length(5); match->set_to_length(5);
        auto* clip = second.mutable_path()->mutable_mapping(0)->add_edit();
        clip->set_to_length(30); clip->set_sequence(string(30, 'A'));
        CHECK(surjector.diploid_pair_span(first, second, 0) == optional<int64_t>(115));
    }
    SECTION("Terminal reference gaps do not inflate the fragment length") {
        auto* mapping = second.mutable_path()->mutable_mapping(0);
        mapping->clear_edit();
        // This is terminal on the right in reference order because the mate is reversed.
        mapping->add_edit()->set_from_length(10);
        auto* match = mapping->add_edit(); match->set_from_length(20); match->set_to_length(20);
        CHECK(surjector.diploid_pair_span(first, second, 100) == optional<int64_t>(100));
        mapping->clear_edit();
        mapping->add_edit()->set_from_length(30);
        CHECK_FALSE(surjector.diploid_pair_span(first, second, 0));
    }
    SECTION("Coordinates cannot overflow") {
        second.mutable_refpos(0)->set_offset(numeric_limits<int64_t>::max() - 5);
        CHECK_FALSE(surjector.diploid_pair_span(first, second, 0));
        CHECK_THROWS_AS(surjector.diploid_pair_span(first, second, -1), invalid_argument);
    }
    SECTION("Fixed fragment model changes score but not compatibility") {
        Surjector::DiploidPairingParameters parameters;
        CHECK(surjector.score_diploid_pair(first, second, 100, parameters) == Approx(60.0));
        parameters.fragment_model = Surjector::DiploidFragmentModel{100, 10};
        CHECK(surjector.score_diploid_pair(first, second, 100, parameters) == Approx(60.0));
        const double log_base = surjector.get_aligner()->scorer->get_log_base();
        CHECK(surjector.score_diploid_pair(first, second, 110, parameters) == Approx(60.0 - 0.5 / log_base));
        CHECK(surjector.score_diploid_pair(first, second, 90, parameters) == Approx(60.0 - 0.5 / log_base));
        parameters.fragment_model->mean = 500;
        CHECK(surjector.diploid_pair_span(first, second, 100).has_value());
        CHECK(surjector.score_diploid_pair(first, second, 100, parameters) < 60.0);
        parameters.fragment_model->stddev = 0;
        CHECK_THROWS_AS(surjector.score_diploid_pair(first, second, 100, parameters), invalid_argument);
        parameters.fragment_model->stddev = numeric_limits<double>::quiet_NaN();
        CHECK_THROWS_AS(surjector.score_diploid_pair(first, second, 100, parameters), invalid_argument);
    }
}

TEST_CASE("Diploid pair selection uses joint evidence and preserves source identity", "[surject][diploid][paired]") {
    bdsg::HashGraph graph;
    string sequence;
    for (size_t i = 0; i < 100; ++i) sequence += "ACGT";
    auto a = graph.create_handle(sequence), b = graph.create_handle(sequence);
    auto path_a = graph.create_path_handle("A"), path_b = graph.create_path_handle("B");
    graph.append_step(path_a, a); graph.append_step(path_b, b);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.prune_suspicious_anchors = false;
    auto placement = [&](handle_t node, int64_t start, bool reverse, bool secondary) {
        Alignment aln;
        aln.set_name(reverse ? "read/2" : "read/1");
        aln.set_sequence(sequence.substr(0, 20));
        aln.set_mapping_quality(reverse ? 17 : 90);
        aln.set_is_secondary(secondary);
        (reverse ? aln.mutable_fragment_prev() : aln.mutable_fragment_next())->set_name(reverse ? "read/1" : "read/2");
        auto* mapping = aln.mutable_path()->add_mapping();
        mapping->set_rank(1);
        auto* position = mapping->mutable_position();
        position->set_node_id(graph.get_id(node));
        position->set_is_reverse(reverse);
        position->set_offset(reverse ? sequence.size() - start - 20 : start);
        auto* edit = mapping->add_edit(); edit->set_from_length(20); edit->set_to_length(20);
        return aln;
    };
    vector<pair<Alignment, Alignment>> input{
        {placement(a, 0, false, false), placement(a, 80, true, false)},
        {placement(b, 0, false, true), placement(b, 180, true, true)}};
    const auto original = input;
    Surjector::DiploidPairingParameters parameters;
    parameters.fragment_model = Surjector::DiploidFragmentModel{200, 10};
    unordered_set<path_handle_t> paths{path_a, path_b};
    SECTION("Fragment evidence promotes an alternative pair and preserves each mate's source MAPQ") {
        auto output = surjector.surject_diploid_paired(input, paths, parameters);
        REQUIRE(output.size() == 2);
        CHECK(output.front().source_pair_index == 1);
        CHECK(output.front().compatible);
        for (size_t i = 0; i < output.size(); ++i) {
            const auto& pair = output[i];
            CHECK(pair.first.front().refpos(0).name() == (i == 0 ? "B" : "A"));
            CHECK(pair.second.front().refpos(0).name() == pair.first.front().refpos(0).name());
            CHECK(pair.first.front().is_secondary() == (i != 0));
            CHECK(pair.second.front().is_secondary() == (i != 0));
            CHECK(get_annotation<double>(pair.first.front(), "diploid_source_mapping_quality") == 90);
            CHECK(get_annotation<double>(pair.second.front(), "diploid_source_mapping_quality") == 17);
            CHECK(pair.first.front().mapping_quality() <= 60);
            CHECK(pair.second.front().mapping_quality() <= 17);
            CHECK(pair.first.front().fragment_next().refpos(0).name() == pair.second.front().refpos(0).name());
        }
        CHECK(output.front().score > output.back().score);
        for (size_t i = 0; i < input.size(); ++i) {
            CHECK(input[i].first.SerializeAsString() == original[i].first.SerializeAsString());
            CHECK(input[i].second.SerializeAsString() == original[i].second.SerializeAsString());
        }
    }
    SECTION("Fragment learning excludes ambiguity, unknown quality, and split candidates") {
        input.resize(1);
        input.front().first.set_mapping_quality(60);
        input.front().second.set_mapping_quality(60);
        parameters.fragment_model.reset();
        auto output = surjector.surject_diploid_paired(input, paths, parameters);
        CHECK(surjector.diploid_fragment_length(output) == optional<int64_t>(100));
        output.front().first.front().set_mapping_quality(3);
        CHECK_FALSE(surjector.diploid_fragment_length(output));
        output.front().first.front().set_mapping_quality(60);
        set_annotation(output.front().first.front(), "diploid_source_mapping_quality", 255);
        CHECK_FALSE(surjector.diploid_fragment_length(output));
        set_annotation(output.front().first.front(), "diploid_source_mapping_quality", 60);
        output.front().first.push_back(output.front().first.front());
        CHECK_FALSE(surjector.diploid_fragment_length(output));
    }
    SECTION("Primary may appear last and source_pair_index still refers to caller order") {
        reverse(input.begin(), input.end());
        auto output = surjector.surject_diploid_paired(input, paths, parameters);
        REQUIRE(output.size() == 2);
        CHECK(output.front().source_pair_index == 0);
        CHECK(output.front().first.front().refpos(0).name() == "B");
    }
    SECTION("No fragment model uses alignment scores with deterministic ties") {
        parameters.fragment_model.reset();
        auto output = surjector.surject_diploid_paired(input, paths, parameters);
        REQUIRE(output.size() == 2);
        CHECK(output.front().first.front().refpos(0).name() == "A");
        reverse(input.begin(), input.end());
        output = surjector.surject_diploid_paired(input, paths, parameters);
        CHECK(output.front().first.front().refpos(0).name() == "A");
    }
    SECTION("Duplicate source pairs do not become additional competitors") {
        auto duplicate = input.front();
        duplicate.first.set_is_secondary(true); duplicate.second.set_is_secondary(true);
        input.insert(input.begin(), duplicate);
        CHECK(surjector.surject_diploid_paired(input, paths, parameters).size() == 2);
    }
    SECTION("Unknown and zero source qualities are handled per mate") {
        input.resize(1);
        input[0].first.set_mapping_quality(255);
        input[0].second.set_mapping_quality(0);
        auto output = surjector.surject_diploid_paired(input, paths, parameters);
        REQUIRE(output.size() == 1);
        CHECK(output.front().first.front().mapping_quality() == 60);
        CHECK(output.front().second.front().mapping_quality() == 0);
        CHECK(get_annotation<double>(output.front().first.front(), "diploid_source_mapping_quality") == 255);
    }
    SECTION("A hard distance limit excludes otherwise high scoring pairs before scoring") {
        parameters.maximum_fragment_length = 150;
        auto output = surjector.surject_diploid_paired(input, paths, parameters);
        REQUIRE(output.size() == 1);
        CHECK(output.front().compatible);
        CHECK(output.front().source_pair_index == 0);
    }
    SECTION("No compatible pairs use a documented alignment-only fallback") {
        parameters.maximum_fragment_length = 50;
        auto output = surjector.surject_diploid_paired(input, paths, parameters);
        REQUIRE(output.size() == 1);
        CHECK_FALSE(output.front().compatible);
        CHECK_FALSE(get_annotation<bool>(output.front().first.front(), "proper_pair"));
        CHECK(output.front().source_pair_index == 0);
        CHECK(output.front().first.front().mapping_quality() == 0);
        CHECK_FALSE(has_annotation(output.front().first.front(), "diploid_haplotype_quality"));
    }
    SECTION("Unmapped mates retain a reference placeholder and reciprocal links") {
        input.resize(1);
        input.front().first.clear_path();
        SECTION("Both mates unmapped") { input.front().second.clear_path(); }
        SECTION("Only mate one unmapped") {}
        auto output = surjector.surject_diploid_paired(input, paths, parameters);
        REQUIRE(output.size() == 1);
        CHECK_FALSE(output.front().compatible);
        REQUIRE(output.front().first.front().refpos_size() == 1);
        CHECK(output.front().first.front().refpos(0).name().empty());
        CHECK(output.front().first.front().refpos(0).offset() == -1);
        CHECK(output.front().first.front().fragment_next().name() == "read/2");
        CHECK(output.front().second.front().fragment_prev().name() == "read/1");
    }
    SECTION("Malformed groups and unsupported supplementary input are rejected") {
        SECTION("Bad mate link") { input[0].first.mutable_fragment_next()->set_name("wrong"); }
        SECTION("Inconsistent primary status") { input[0].second.set_is_secondary(true); }
        SECTION("Multiple primaries") { input[1].first.set_is_secondary(false); input[1].second.set_is_secondary(false); }
        SECTION("No primary") { input[0].first.set_is_secondary(true); input[0].second.set_is_secondary(true); }
        SECTION("Different mate sequence") { input[1].second.set_sequence("ACGT"); }
        SECTION("Supplementary annotation") { set_annotation(input[0].first, "supplementary", true); }
        SECTION("Embedded supplementary") { input[0].second.add_supplementary(); }
        SECTION("Invalid MAPQ") { input[0].first.set_mapping_quality(256); }
        SECTION("Zero model spread") { parameters.fragment_model->stddev = 0; }
        CHECK_THROWS_AS(surjector.surject_diploid_paired(input, paths, parameters), invalid_argument);
    }
    SECTION("Empty input has no output") {
        CHECK(surjector.surject_diploid_paired({}, paths, parameters).empty());
    }
}

// Construct one exact graph placement with MAPQ 37 and primary status.
static Alignment diploid_test_placement(int64_t node_id, const string& sequence) {
    Alignment aln;
    aln.set_name("read");
    aln.set_sequence(sequence);
    aln.set_mapping_quality(37);
    auto* mapping = aln.mutable_path()->add_mapping();
    mapping->set_rank(1);
    mapping->mutable_position()->set_node_id(node_id);
    auto* edit = mapping->add_edit();
    edit->set_from_length(sequence.size());
    edit->set_to_length(sequence.size());
    return aln;
}

TEST_CASE("Diploid selection compares graph placements deterministically", "[surject][diploid]") {
    bdsg::HashGraph graph;
    const string sequence = "ACGTTGCACTGATCGATGCA";
    auto a = graph.create_handle(sequence);
    auto b = graph.create_handle(sequence);
    auto path_a = graph.create_path_handle("A");
    auto path_b = graph.create_path_handle("B");
    graph.append_step(path_a, a);
    graph.append_step(path_b, b);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.prune_suspicious_anchors = false;
    auto primary = diploid_test_placement(graph.get_id(b), sequence);
    auto alternative = diploid_test_placement(graph.get_id(a), sequence);
    alternative.set_is_secondary(true);
    vector<Alignment> input{alternative, primary};
    const auto original = input;

    SECTION("Primary input need not be first and tied output is deterministic") {
        auto output = surjector.surject_diploid(input, {path_a, path_b});
        REQUIRE(output.size() == 2);
        CHECK(input[0].SerializeAsString() == original[0].SerializeAsString());
        CHECK(input[1].SerializeAsString() == original[1].SerializeAsString());
        CHECK(output[0].refpos(0).name() == "A");
        CHECK_FALSE(output[0].is_secondary());
        CHECK(output[1].is_secondary());
        for (const auto& aln : output) {
            CHECK(get_annotation<double>(aln, "diploid_source_mapping_quality") == 37);
            CHECK(get_annotation<bool>(aln, "diploid_haplotype_preferred"));
            CHECK(aln.mapping_quality() == 3);
        }
        std::reverse(input.begin(), input.end());
        auto reversed = surjector.surject_diploid(input, {path_b, path_a});
        REQUIRE(reversed.size() == output.size());
        for (size_t i = 0; i < output.size(); ++i) {
            CHECK(reversed[i].refpos(0).name() == output[i].refpos(0).name());
            CHECK(reversed[i].mapping_quality() == output[i].mapping_quality());
        }
    }
    SECTION("A better mapped alternative becomes primary across graph placements") {
        auto short_node = graph.create_handle(sequence.substr(0, 8));
        auto short_path = graph.create_path_handle("short");
        graph.append_step(short_path, short_node);
        primary.clear_path();
        auto* mapping = primary.mutable_path()->add_mapping();
        mapping->set_rank(1);
        mapping->mutable_position()->set_node_id(graph.get_id(short_node));
        auto* edit = mapping->add_edit();
        edit->set_from_length(8);
        edit->set_to_length(8);
        edit = mapping->add_edit();
        edit->set_to_length(sequence.size() - 8);
        edit->set_sequence(sequence.substr(8));
        bdsg::PositionOverlay updated_overlay(&graph);
        TestSurjector updated(&updated_overlay);
        updated.prune_suspicious_anchors = false;
        auto output = updated.surject_diploid({primary, alternative}, {short_path, path_a});
        REQUIRE(output.size() == 2);
        CHECK(output.front().refpos(0).name() == "A");
        CHECK(output.front().score() > output.back().score());
        CHECK_FALSE(output.front().is_secondary());
        CHECK(output.back().is_secondary());
        CHECK(output.front().mapping_quality() <= 37);
    }
}

TEST_CASE("Diploid mapping quality respects known input confidence", "[surject][diploid]") {
    bdsg::HashGraph graph;
    const string sequence = "ACGTTGCACTGATCGATGCA";
    auto b = graph.create_handle(sequence);
    auto path_b = graph.create_path_handle("B");
    graph.append_step(path_b, b);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.prune_suspicious_anchors = false;
    auto primary = diploid_test_placement(graph.get_id(b), sequence);

    SECTION("Original high MAPQ is preserved separately from computed quality") {
        primary.set_mapping_quality(90);
        auto output = surjector.surject_diploid({primary}, {path_b});
        REQUIRE(output.size() == 1);
        CHECK(output.front().mapping_quality() == 60);
        CHECK(get_annotation<double>(output.front(), "diploid_source_mapping_quality") == 90);
    }
    SECTION("Unknown source MAPQ is not treated as a confidence cap") {
        primary.set_mapping_quality(255);
        surjector.max_diploid_mapping_quality = 42;
        auto output = surjector.surject_diploid({primary}, {path_b});
        REQUIRE(output.size() == 1);
        CHECK(output.front().mapping_quality() == 42);
        CHECK(get_annotation<double>(output.front(), "diploid_source_mapping_quality") == 255);
    }
    SECTION("Zero source MAPQ caps global confidence") {
        primary.set_mapping_quality(0);
        auto output = surjector.surject_diploid({primary}, {path_b});
        REQUIRE(output.size() == 1);
        CHECK(output.front().mapping_quality() == 0);
    }
}

TEST_CASE("Diploid haplotype quality compares targets of one graph placement", "[surject][diploid]") {
    bdsg::HashGraph graph;
    const string sequence = "ACGTTGCACTGATCGATGCA";
    auto b = graph.create_handle(sequence);
    auto path_b = graph.create_path_handle("B");
    graph.append_step(path_b, b);
    auto primary = diploid_test_placement(graph.get_id(b), sequence);

    auto other = graph.create_path_handle("C");
    graph.append_step(other, b);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.prune_suspicious_anchors = false;
    surjector.multimap_to_all_paths = false;
    auto output = surjector.surject_diploid({primary}, {path_b, other});
    REQUIRE(output.size() == 2);
    CHECK(get_annotation<bool>(output[0], "diploid_haplotype_preferred"));
    CHECK_FALSE(get_annotation<bool>(output[1], "diploid_haplotype_preferred"));
    CHECK(output[0].refpos(0).name() == "B");
    CHECK(output[1].refpos(0).name() == "C");
    CHECK_FALSE(output[0].is_secondary());
    CHECK(output[1].is_secondary());
    for (const auto& aln : output) {
        CHECK(get_annotation<double>(aln, "diploid_haplotype_quality") == 3);
        CHECK(aln.mapping_quality() == 3);
    }
}

TEST_CASE("Diploid reads without candidates retain unmapped metadata", "[surject][diploid]") {
    bdsg::HashGraph graph;
    const string sequence = "ACGTTGCACTGATCGATGCA";
    auto a = graph.create_handle(sequence);
    auto b = graph.create_handle(sequence);
    auto path_a = graph.create_path_handle("A");
    auto path_b = graph.create_path_handle("B");
    graph.append_step(path_a, a);
    graph.append_step(path_b, b);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.prune_suspicious_anchors = false;
    auto primary = diploid_test_placement(graph.get_id(b), sequence);

    set_annotation(primary, "tags", string("ZZ:Z:keep\tSA:Z:obsolete"));
    unordered_set<path_handle_t> targets;
    SECTION("No selected target paths") {
        targets = {};
    }
    SECTION("Placement does not overlap the selected target") {
        targets = {path_a};
    }
    SECTION("Input placement has an empty path") {
        primary.clear_path();
        targets = {path_a, path_b};
    }
    auto output = surjector.surject_diploid({primary}, targets);
    REQUIRE(output.size() == 1);
    REQUIRE(output.front().refpos_size() == 1);
    CHECK(output.front().refpos(0).name().empty());
    CHECK(output.front().refpos(0).offset() == -1);
    CHECK_FALSE(output.front().refpos(0).is_reverse());
    CHECK(output.front().path().mapping_size() == 0);
    CHECK(output.front().mapping_quality() == 0);
    CHECK(get_annotation<string>(output.front(), "tags") == "ZZ:Z:keep");
    CHECK(get_annotation<double>(output.front(), "diploid_source_mapping_quality") == 37);
}

TEST_CASE("Diploid input requires compatible placements of one unpaired read", "[surject][diploid]") {
    bdsg::HashGraph graph;
    const string sequence = "ACGTTGCACTGATCGATGCA";
    auto a = graph.create_handle(sequence);
    auto b = graph.create_handle(sequence);
    auto path_a = graph.create_path_handle("A");
    auto path_b = graph.create_path_handle("B");
    graph.append_step(path_a, a);
    graph.append_step(path_b, b);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.prune_suspicious_anchors = false;
    auto primary = diploid_test_placement(graph.get_id(b), sequence);
    auto alternative = diploid_test_placement(graph.get_id(a), sequence);
    alternative.set_is_secondary(true);

    SECTION("Malformed groups are rejected") {
        CHECK_THROWS_AS(surjector.surject_diploid({alternative}, {path_a}), invalid_argument);
        alternative.set_is_secondary(false);
        CHECK_THROWS_AS(surjector.surject_diploid({primary, alternative}, {path_a}), invalid_argument);
        alternative.set_is_secondary(true);
        alternative.set_name("different");
        CHECK_THROWS_AS(surjector.surject_diploid({primary, alternative}, {path_a}), invalid_argument);
        alternative = diploid_test_placement(graph.get_id(a), sequence);
        alternative.set_is_secondary(true);
        alternative.set_sequence("ACGT");
        CHECK_THROWS_AS(surjector.surject_diploid({primary, alternative}, {path_a}), invalid_argument);
        alternative = diploid_test_placement(graph.get_id(a), sequence);
        alternative.set_is_secondary(true);
        alternative.set_quality(string(sequence.size(), 30));
        CHECK_THROWS_AS(surjector.surject_diploid({primary, alternative}, {path_a}), invalid_argument);
    }
    SECTION("Unsupported input relationships fail explicitly") {
        auto paired = primary;
        paired.mutable_fragment_next()->set_name("mate");
        CHECK_THROWS_AS(surjector.surject_diploid({paired}, {path_b}), invalid_argument);
        auto supplementary = primary;
        set_annotation(supplementary, "supplementary", true);
        CHECK_THROWS_AS(surjector.surject_diploid({supplementary}, {path_b}), invalid_argument);
        auto embedded = primary;
        embedded.add_supplementary();
        CHECK_THROWS_AS(surjector.surject_diploid({embedded}, {path_b}), invalid_argument);
    }
}

TEST_CASE("Paired diploid alternatives retain candidate-local supplementary mate links",
          "[surject][diploid][paired]") {
    bdsg::HashGraph graph;
    vector<Alignment> input;
    unordered_set<path_handle_t> paths;
    for (string name : {"A", "B"}) {
        auto left = graph.create_handle(string(60, 'A'));
        auto right = graph.create_handle(string(40, 'C'));
        auto gap = graph.create_handle(string(200, 'G'));
        graph.create_edge(left, gap);
        graph.create_edge(gap, right);
        graph.create_edge(left, right);
        auto path = graph.create_path_handle(name);
        graph.append_step(path, left);
        graph.append_step(path, gap);
        graph.append_step(path, right);
        paths.insert(path);
        Alignment source;
        source.set_name("split");
        source.set_sequence(string(60, 'A') + string(40, 'C'));
        source.set_mapping_quality(name == "B" ? 17 : 2);
        source.set_is_secondary(name == "A");
        set_annotation(source, "tags", string("ZZ:Z:keep\tSA:Z:obsolete"));
        for (auto node : {left, right}) {
            auto* mapping = source.mutable_path()->add_mapping();
            mapping->set_rank(source.path().mapping_size());
            mapping->mutable_position()->set_node_id(graph.get_id(node));
            auto* edit = mapping->add_edit();
            edit->set_from_length(graph.get_length(node));
            edit->set_to_length(graph.get_length(node));
        }
        input.push_back(source);
    }
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.prune_suspicious_anchors = false;
    surjector.report_supplementary = true;
    vector<pair<Alignment, Alignment>> pairs;
    for (auto source : input) {
        source.mutable_fragment_next()->set_name("mate");
        set_annotation(source, "mate_info", string("obsolete"));
        Alignment mate;
        mate.set_name("mate");
        mate.set_sequence(string(40, 'G'));
        mate.set_mapping_quality(27);
        mate.set_is_secondary(source.is_secondary());
        mate.mutable_fragment_prev()->set_name(source.name());
        auto* mapping = mate.mutable_path()->add_mapping();
        mapping->set_rank(1);
        mapping->mutable_position()->set_node_id(source.path().mapping(1).position().node_id());
        mapping->mutable_position()->set_is_reverse(true);
        auto* edit = mapping->add_edit(); edit->set_from_length(40); edit->set_to_length(40);
        pairs.emplace_back(std::move(source), std::move(mate));
    }
    Surjector::DiploidPairingParameters parameters;
    auto output = surjector.surject_diploid_paired(pairs, paths, parameters);
    REQUIRE(output.size() == 2);
    for (size_t i = 0; i < output.size(); ++i) {
        const auto& pair = output[i];
        const string name = i == 0 ? "A" : "B";
        CHECK(pair.source_pair_index == i);
        REQUIRE(pair.first.size() == 2);
        REQUIRE(pair.second.size() == 1);
        CHECK_FALSE(has_annotation(pair.first.front(), "mate_info"));
        CHECK(is_supplementary(pair.first.back()));
        CHECK(pair.first.back().is_secondary() == (i != 0));
        CHECK(pair.second.front().is_secondary() == (i != 0));
        REQUIRE(has_annotation(pair.first.back(), "mate_info"));
        auto mate = parse_mate_info(get_annotation<string>(pair.first.back(), "mate_info"));
        CHECK(get<0>(mate) == name);
        CHECK(get<1>(mate) == pair.second.front().refpos(0).offset());
        CHECK(get<2>(mate) == pair.second.front().refpos(0).is_reverse());
        CHECK_FALSE(get<3>(mate));
        for (const auto& piece : pair.first) {
            CHECK(piece.fragment_next().refpos(0).name() == name);
            CHECK(piece.mapping_quality() == pair.first.front().mapping_quality());
            CHECK(get_annotation<double>(piece, "diploid_source_mapping_quality") == 17);
            auto tags = get_annotation<string>(piece, "tags");
            CHECK(tags.find("SA:Z:" + name + ",") != string::npos);
            CHECK(tags.find("obsolete") == string::npos);
            CHECK(count(tags.begin(), tags.end(), ';') == 1);
        }
    }
}

TEST_CASE("Surjection alternatives are classified by read interval",
          "[surject][interval-classification]") {
    bdsg::HashGraph graph;
    auto node = graph.create_handle(string(200, 'A'));
    auto path = graph.create_path_handle("ref");
    auto step = graph.append_step(path, node);
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);

    vector<pair<Alignment, pair<step_handle_t, step_handle_t>>> candidates;
    auto add_candidate = [&](const string& name, size_t begin, size_t end,
                             int32_t score, size_t reference_begin) {
        Alignment aln;
        aln.set_name(name);
        aln.set_sequence(string(100, 'A'));
        aln.set_score(score);
        auto* mapping = aln.mutable_path()->add_mapping();
        mapping->set_rank(1);
        mapping->mutable_position()->set_node_id(graph.get_id(node));
        mapping->mutable_position()->set_offset(reference_begin);
        if (begin != 0) {
            auto* clip = mapping->add_edit();
            clip->set_to_length(begin);
            clip->set_sequence(string(begin, 'A'));
        }
        auto* match = mapping->add_edit();
        match->set_from_length(end - begin);
        match->set_to_length(end - begin);
        if (end != 100) {
            auto* clip = mapping->add_edit();
            clip->set_to_length(100 - end);
            clip->set_sequence(string(100 - end, 'A'));
        }
        candidates.emplace_back(aln, make_pair(step, step));
    };

    // The best candidate is deliberately not first. A and B overlap on
    // the read; C covers the remaining read bases.
    add_candidate("A", 0, 60, 40, 0);
    add_candidate("B", 0, 60, 60, 100);
    add_candidate("C", 60, 100, 30, 160);

    auto by_name = [&]() {
        map<string, Alignment> named;
        for (const auto& candidate : candidates) {
            REQUIRE(named.emplace(candidate.first.name(), candidate.first).second);
        }
        return named;
    };

    SECTION("Supplementary reporting retains the disjoint candidate") {
        surjector.report_supplementary = true;
        surjector.choose_primary(candidates);

        REQUIRE(candidates.size() == 3);
        CHECK(candidates[0].first.name() == "B");
        CHECK(candidates[1].first.name() == "A");
        CHECK(candidates[2].first.name() == "C");
        const auto named = by_name();
        CHECK_FALSE(named.at("B").is_secondary());
        CHECK_FALSE(is_supplementary(named.at("B")));
        CHECK(named.at("B").score() == 60);
        CHECK(named.at("A").is_secondary());
        CHECK_FALSE(is_supplementary(named.at("A")));
        CHECK_FALSE(named.at("C").is_secondary());
        CHECK(is_supplementary(named.at("C")));
    }
    SECTION("Without supplementary reporting the disjoint candidate is omitted") {
        surjector.report_supplementary = false;
        surjector.choose_primary(candidates);

        REQUIRE(candidates.size() == 2);
        CHECK(candidates[0].first.name() == "B");
        CHECK(candidates[1].first.name() == "A");
        const auto named = by_name();
        CHECK_FALSE(named.at("B").is_secondary());
        CHECK_FALSE(is_supplementary(named.at("B")));
        CHECK(named.at("B").score() == 60);
        CHECK(named.at("A").is_secondary());
        CHECK_FALSE(is_supplementary(named.at("A")));
        CHECK(named.count("C") == 0);
    }
}



TEST_CASE("Primary strand selection ignores secondary alternatives",
          "[surject][interval-classification][strand-scoring]") {
    bdsg::HashGraph graph;
    auto node = graph.create_handle(string(1000, 'A'));
    auto path_a = graph.create_path_handle("A");
    auto path_b = graph.create_path_handle("B");
    auto step_a = graph.append_step(path_a, node);
    auto step_b = graph.append_step(path_b, node);
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    surjector.report_supplementary = true;

    using Candidate = pair<Alignment, pair<step_handle_t, step_handle_t>>;
    unordered_map<pair<path_handle_t, bool>, vector<Candidate>> candidates;
    auto add_candidate = [&](path_handle_t path, step_handle_t step,
                             size_t begin, size_t end, int32_t score,
                             size_t reference_begin) {
        Alignment aln;
        aln.set_sequence(string(100, 'A'));
        aln.set_score(score);
        auto* mapping = aln.mutable_path()->add_mapping();
        mapping->set_rank(1);
        mapping->mutable_position()->set_node_id(graph.get_id(node));
        mapping->mutable_position()->set_offset(reference_begin);
        if (begin != 0) {
            auto* clip = mapping->add_edit();
            clip->set_to_length(begin);
            clip->set_sequence(string(begin, 'A'));
        }
        auto* match = mapping->add_edit();
        match->set_from_length(end - begin);
        match->set_to_length(end - begin);
        if (end != 100) {
            auto* clip = mapping->add_edit();
            clip->set_to_length(100 - end);
            clip->set_sequence(string(100 - end, 'A'));
        }
        candidates[make_pair(path, false)].emplace_back(
            aln, make_pair(step, step));
    };

    SECTION("A wins with 110 despite B having three score-105 placements") {
        add_candidate(path_a, step_a, 0, 100, 110, 0);
        add_candidate(path_b, step_b, 0, 100, 105, 0);
        add_candidate(path_b, step_b, 0, 100, 105, 200);
        add_candidate(path_b, step_b, 0, 100, 105, 400);
        surjector.choose_primary(candidates.at(make_pair(path_a, false)));
        surjector.choose_primary(candidates.at(make_pair(path_b, false)));

        const auto& b = candidates.at(make_pair(path_b, false));
        REQUIRE(b.size() == 3);
        CHECK_FALSE(b[0].first.is_secondary());
        CHECK_FALSE(is_supplementary(b[0].first));
        CHECK(b[1].first.is_secondary());
        CHECK_FALSE(is_supplementary(b[1].first));
        CHECK(b[2].first.is_secondary());
        CHECK_FALSE(is_supplementary(b[2].first));
        CHECK(surjector.choose_primary_strand(candidates) == make_pair(path_a, false));
    }
    SECTION("B wins with disjoint scores 65 plus 45 against A's 105") {
        add_candidate(path_a, step_a, 0, 100, 105, 0);
        add_candidate(path_b, step_b, 0, 60, 65, 0);
        add_candidate(path_b, step_b, 60, 100, 45, 200);
        surjector.choose_primary(candidates.at(make_pair(path_a, false)));
        surjector.choose_primary(candidates.at(make_pair(path_b, false)));

        const auto& b = candidates.at(make_pair(path_b, false));
        REQUIRE(b.size() == 2);
        CHECK_FALSE(b[0].first.is_secondary());
        CHECK_FALSE(is_supplementary(b[0].first));
        CHECK(b[0].first.score() == 65);
        CHECK_FALSE(b[1].first.is_secondary());
        CHECK(is_supplementary(b[1].first));
        CHECK(b[1].first.score() == 45);
        CHECK(surjector.choose_primary_strand(candidates) == make_pair(path_b, false));
    }
}


TEST_CASE("Secondary alternatives do not join the supplementary group",
          "[surject][interval-classification]") {
    bdsg::HashGraph graph;
    auto left = graph.create_handle(string(60, 'A'));
    auto right = graph.create_handle(string(40, 'C'));
    auto spacer = graph.create_handle(string(200, 'G'));
    graph.create_edge(left, right);
    graph.create_edge(right, spacer);
    graph.create_edge(spacer, right);
    auto main_path = graph.create_path_handle("main");
    graph.append_step(main_path, left);
    auto repeat_path = graph.create_path_handle("repeat");
    graph.append_step(repeat_path, right);
    graph.append_step(repeat_path, spacer);
    graph.append_step(repeat_path, right);
    bdsg::PositionOverlay pos_graph(&graph);
    Surjector surjector(&pos_graph);
    surjector.report_supplementary = true;

    Alignment read;
    read.set_name("split-read");
    read.set_sequence(string(60, 'A') + string(40, 'C'));
    for (auto node : {left, right}) {
        auto* mapping = read.mutable_path()->add_mapping();
        mapping->set_rank(read.path().mapping_size());
        mapping->mutable_position()->set_node_id(graph.get_id(node));
        auto* edit = mapping->add_edit();
        edit->set_from_length(graph.get_length(node));
        edit->set_to_length(graph.get_length(node));
    }

    // Either repeat copy may win the tie. Verify that SA links use the
    // selected copy's position, not the unused alternative's position.
    auto check_pieces_and_links = [&](const vector<Alignment>& output) {
        REQUIRE(output.size() == 3);
        const Alignment* main = nullptr;
        const Alignment* supplementary = nullptr;
        const Alignment* alternative = nullptr;
        for (const auto& aln : output) {
            REQUIRE(aln.refpos_size() == 1);
            CHECK_FALSE(aln.refpos(0).is_reverse());
            if (aln.refpos(0).name() == "main") {
                REQUIRE(main == nullptr);
                main = &aln;
            } else {
                REQUIRE(aln.refpos(0).name() == "repeat");
                if (is_supplementary(aln)) {
                    REQUIRE(supplementary == nullptr);
                    supplementary = &aln;
                } else {
                    REQUIRE(alternative == nullptr);
                    alternative = &aln;
                }
            }
        }
        REQUIRE(main != nullptr);
        REQUIRE(supplementary != nullptr);
        REQUIRE(alternative != nullptr);
        CHECK(main->refpos(0).offset() == 0);
        CHECK_FALSE(is_supplementary(*main));

        // The 40-base piece occurs at offsets 0 and 240 on the repeat path.
        const auto chosen = supplementary->refpos(0).offset();
        const auto unused = alternative->refpos(0).offset();
        CHECK((chosen == 0 || chosen == 240));
        CHECK((unused == 0 || unused == 240));
        CHECK(chosen != unused);

        const auto main_tags = parse_sam_tags_for_test(*main);
        const auto supplementary_tags = parse_sam_tags_for_test(*supplementary);
        const auto alternative_tags = parse_sam_tags_for_test(*alternative);
        REQUIRE(main_tags.count("SA") == 1);
        CHECK(main_tags.at("SA").first == 'Z');
        CHECK(main_tags.at("SA").second
              == "repeat," + to_string(chosen + 1) + ",+,60S40M,0,0;");
        REQUIRE(supplementary_tags.count("SA") == 1);
        CHECK(supplementary_tags.at("SA").first == 'Z');
        CHECK(supplementary_tags.at("SA").second == "main,1,+,60M40S,0,0;");
        CHECK(alternative_tags.count("SA") == 0);
        return make_tuple(main, supplementary, alternative);
    };

    SECTION("Primary input produces a primary, a supplementary, and a secondary") {
        read.set_is_secondary(false);
        const auto output = surjector.surject(read, {main_path, repeat_path});
        const auto pieces = check_pieces_and_links(output);
        CHECK_FALSE(get<0>(pieces)->is_secondary());
        CHECK_FALSE(get<1>(pieces)->is_secondary());
        CHECK(get<2>(pieces)->is_secondary());
    }
    SECTION("Secondary input passes secondary status to all three outputs") {
        read.set_is_secondary(true);
        const auto output = surjector.surject(read, {main_path, repeat_path});
        const auto pieces = check_pieces_and_links(output);
        CHECK(get<0>(pieces)->is_secondary());
        CHECK(get<1>(pieces)->is_secondary());
        CHECK(get<2>(pieces)->is_secondary());
    }
}



// Add one exact-match mapping with explicit orientation and coordinates.
static void append_anchor_match(path_t& anchor, int64_t node_id,
                                bool reverse, size_t offset, size_t length) {
    auto* mapping = anchor.add_mapping();
    mapping->mutable_position()->set_node_id(node_id);
    mapping->mutable_position()->set_is_reverse(reverse);
    mapping->mutable_position()->set_offset(offset);
    auto* edit = mapping->add_edit();
    edit->set_from_length(length);
    edit->set_to_length(length);
}

TEST_CASE("Target-path repeats respect the slide limit", "[surject][anchor-sliding]") {
    bdsg::HashGraph graph;
    auto node = graph.create_handle("ACGATTACGA");
    auto path = graph.create_path_handle("ref");
    auto step = graph.append_step(path, node);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);

    // ACGA occurs at path offsets 0 and 6. The read has only one occurrence.
    string read = "ACGA";
    path_t mappings;
    SECTION("A duplicate six bases to the right is detected") {
        append_anchor_match(mappings, graph.get_id(node), false, 0, 4);
        Surjector::path_chunk_t chunk{{read.begin(), read.end()}, mappings};
        surjector.max_slide = 6;
        CHECK(surjector.anchor_has_nearby_repeat(read, chunk, {step, step}));
    }
    SECTION("A duplicate six bases to the left is detected") {
        append_anchor_match(mappings, graph.get_id(node), false, 6, 4);
        Surjector::path_chunk_t chunk{{read.begin(), read.end()}, mappings};
        surjector.max_slide = 6;
        CHECK(surjector.anchor_has_nearby_repeat(read, chunk, {step, step}));
    }
    SECTION("A duplicate beyond the radius is ignored") {
        append_anchor_match(mappings, graph.get_id(node), false, 0, 4);
        Surjector::path_chunk_t chunk{{read.begin(), read.end()}, mappings};
        surjector.max_slide = 5;
        CHECK_FALSE(surjector.anchor_has_nearby_repeat(read, chunk, {step, step}));
    }
    SECTION("Zero disables both read and target-path sliding") {
        read = "ACGATTACGA";
        append_anchor_match(mappings, graph.get_id(node), false, 0, 4);
        Surjector::path_chunk_t chunk{{read.begin(), read.begin() + 4}, mappings};
        surjector.max_slide = 0;
        CHECK_FALSE(surjector.anchor_has_nearby_repeat(read, chunk, {step, step}));
    }
}

TEST_CASE("Read repeats use the anchor's read span", "[surject][anchor-sliding]") {
    bdsg::HashGraph graph;
    auto node = graph.create_handle("ACGATTCGCCCC");
    auto path = graph.create_path_handle("ref");
    auto step = graph.append_step(path, node);
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);

    SECTION("A duplicate in the read is detected without a path duplicate") {
        string read = "ACGATTACGA";
        path_t mappings;
        append_anchor_match(mappings, graph.get_id(node), false, 0, 4);
        Surjector::path_chunk_t chunk{{read.begin(), read.begin() + 4}, mappings};
        surjector.max_slide = 6;
        CHECK(surjector.anchor_has_nearby_repeat(read, chunk, {step, step}));
    }
    SECTION("An insertion makes the read search radius larger") {
        // ACGCTAGA aligns as 2M4I2M to ACGA.
        // The second ACGCTAGA begins at read offset 10.
        string read = "ACGCTAGATTACGCTAGA";
        path_t mappings;
        append_anchor_match(mappings, graph.get_id(node), false, 0, 2);
        auto* mapping = mappings.mutable_mapping(0);
        auto* insertion = mapping->add_edit();
        insertion->set_to_length(4);
        insertion->set_sequence("GCTA");
        auto* match = mapping->add_edit();
        match->set_from_length(2);
        match->set_to_length(2);
        Surjector::path_chunk_t chunk{{read.begin(), read.begin() + 8}, mappings};
        surjector.max_slide = 12;
        // Read radius is min(12, 2*8) = 12, not 2*4 = 8.
        CHECK(surjector.anchor_has_nearby_repeat(read, chunk, {step, step}));
    }
}

TEST_CASE("Reverse anchors spanning nodes use the full reference interval",
          "[surject][anchor-sliding]") {
    bdsg::HashGraph graph;
    auto a = graph.create_handle("ACG");
    auto b = graph.create_handle("ATT");
    auto c = graph.create_handle("ACG");
    auto d = graph.create_handle("A");
    graph.create_edge(a, b);
    graph.create_edge(b, c);
    graph.create_edge(c, d);
    auto path = graph.create_path_handle("ref");
    auto step_a = graph.append_step(path, a);
    auto step_b = graph.append_step(path, b);
    graph.append_step(path, c);
    graph.append_step(path, d);

    // Path: ACG | ATT | ACG | A. ACGA occurs at offsets 0 and 6.
    // Read TCGT traverses the first ACGA backward: b's first base, then a.
    string read = "TCGT";
    path_t mappings;
    append_anchor_match(mappings, graph.get_id(b), true, 2, 1);
    append_anchor_match(mappings, graph.get_id(a), true, 0, 3);
    Surjector::path_chunk_t chunk{{read.begin(), read.end()}, mappings};
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.max_slide = 6;
    CHECK(surjector.anchor_has_nearby_repeat(read, chunk, {step_b, step_a}));
}

TEST_CASE("Reverse mappings on reverse path steps use forward path coordinates",
          "[surject][anchor-sliding]") {
    bdsg::HashGraph graph;
    // The reverse-oriented step has path sequence ACGATTACGA.
    auto node = graph.create_handle("TCGTAATCGT");
    auto path = graph.create_path_handle("ref");
    auto step = graph.append_step(path, graph.flip(node));

    string read = "ACGA";
    path_t mappings;
    append_anchor_match(mappings, graph.get_id(node), true, 0, 4);
    Surjector::path_chunk_t chunk{{read.begin(), read.end()}, mappings};
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.max_slide = 6;
    CHECK(surjector.anchor_has_nearby_repeat(read, chunk, {step, step}));
}

TEST_CASE("Repeat pruning keeps anchor and step-range vectors synchronized",
          "[surject][anchor-sliding]") {
    bdsg::HashGraph graph;
    auto repeated = graph.create_handle("ACGATTACGA");
    auto unique = graph.create_handle("TGCC");
    graph.create_edge(repeated, unique);
    auto path = graph.create_path_handle("ref");
    auto repeated_step = graph.append_step(path, repeated);
    auto unique_step = graph.append_step(path, unique);

    // ACGA is ambiguous on the path; TGCC is unique.
    string read = "ACGATGCC";
    path_t first, second;
    append_anchor_match(first, graph.get_id(repeated), false, 0, 4);
    append_anchor_match(second, graph.get_id(unique), false, 0, 4);
    vector<Surjector::path_chunk_t> chunks{
        {{read.begin(), read.begin() + 4}, first},
        {{read.begin() + 4, read.end()}, second}
    };
    vector<pair<step_handle_t, step_handle_t>> ranges{
        {repeated_step, repeated_step}, {unique_step, unique_step}
    };
    bdsg::PositionOverlay overlay(&graph);
    TestSurjector surjector(&overlay);
    surjector.prune_suspicious_anchors = true;
    surjector.prune_tail_region_anchors = false;
    surjector.max_slide = 6;
    surjector.max_tail_anchor_prune = 0;
    surjector.max_low_complexity_anchor_prune = 0;
    surjector.max_low_complexity_anchor_trim = 0;
    surjector.max_anchors = 1000;
    surjector.prune_and_trim_anchors(read, chunks, ranges);

    REQUIRE(chunks.size() == 1);
    REQUIRE(ranges.size() == 1);
    CHECK(chunks.front().first.first == read.begin() + 4);
    CHECK(chunks.front().first.second == read.end());
    CHECK(chunks.front().second.mapping(0).position().node_id() == graph.get_id(unique));
    CHECK(ranges.front() == make_pair(unique_step, unique_step));
}


TEST_CASE("Mapper-declared tails control anchor pruning",
          "[surject][tail-pruning]") {
    bdsg::HashGraph graph;
    auto path = graph.create_path_handle("ref");
    vector<handle_t> nodes{
        graph.create_handle("ACGT"),
        graph.create_handle("TGCA"),
        graph.create_handle("GACT")
    };
    graph.create_edge(nodes[0], nodes[1]);
    graph.create_edge(nodes[1], nodes[2]);

    vector<pair<step_handle_t, step_handle_t>> ranges;
    for (auto node : nodes) {
        auto step = graph.append_step(path, node);
        ranges.emplace_back(step, step);
    }
    const auto original_ranges = ranges;

    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    surjector.prune_suspicious_anchors = false;
    surjector.prune_tail_region_anchors = true;
    surjector.max_anchors = 1000;

    const string sequence = "ACGTTGCAGACT";
    vector<Surjector::path_chunk_t> chunks(3);
    for (size_t i = 0; i < chunks.size(); ++i) {
        chunks[i].first.first = sequence.begin() + 4 * i;
        chunks[i].first.second = sequence.begin() + 4 * (i + 1);
        auto* mapping = chunks[i].second.add_mapping();
        mapping->mutable_position()->set_node_id(graph.get_id(nodes[i]));
        auto* edit = mapping->add_edit();
        edit->set_from_length(4);
        edit->set_to_length(4);
    }

    size_t left_tail = 0;
    size_t right_tail = 0;
    vector<size_t> expected{0, 1, 2};

    SECTION("Disabled flag preserves anchors inside declared tails") {
        surjector.prune_tail_region_anchors = false;
        left_tail = 4;
        right_tail = 4;
    }
    SECTION("Zero tail lengths preserve all anchors") {
    }
    SECTION("Left tail removes its fully contained anchor") {
        left_tail = 4;
        expected = {1, 2};
    }
    SECTION("Right tail removes its fully contained anchor") {
        right_tail = 4;
        expected = {0, 1};
    }
    SECTION("Both tails leave only the middle anchor") {
        left_tail = 4;
        right_tail = 4;
        expected = {1};
    }
    SECTION("Anchors crossing tail boundaries are retained") {
        left_tail = 3;
        right_tail = 3;
    }

    surjector.prune_and_trim_anchors(sequence, chunks, ranges,
                                    left_tail, right_tail);

    REQUIRE(chunks.size() == expected.size());
    REQUIRE(ranges.size() == expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        const auto original = expected[i];
        CHECK(chunks[i].first.first == sequence.begin() + 4 * original);
        CHECK(chunks[i].first.second == sequence.begin() + 4 * (original + 1));
        REQUIRE(chunks[i].second.mapping_size() == 1);
        CHECK(chunks[i].second.mapping(0).position().node_id()
              == graph.get_id(nodes[original]));
        CHECK(ranges[i] == original_ranges[original]);
    }
}


TEST_CASE("Tail pruning follows read intervals regardless of node or target-path orientation",
          "[surject][tail-pruning]") {
    // Three 4-base anchors cover read intervals [0,4), [4,8), [8,12).
    // A 4-base tail contains an anchor; a 3-base tail only overlaps one.
    // Changing node orientations or reversing the read on ref must not swap the tails.
    for (bool reverse_middle_node : {false, true}) {
        for (bool reverse_relative_to_path : {false, true}) {
            bdsg::HashGraph graph;
            vector<handle_t> nodes{graph.create_handle("ACGT"),
                                   graph.create_handle(reverse_middle_node ? "GCTT" : "AAGC"),
                                   graph.create_handle("GACT")};
            if (reverse_middle_node) nodes[1] = graph.flip(nodes[1]);
            auto path = graph.create_path_handle("ref");
            vector<step_handle_t> steps;
            for (size_t i = 0; i < nodes.size(); ++i) {
                steps.push_back(graph.append_step(path, nodes[i]));
                if (i) graph.create_edge(nodes[i - 1], nodes[i]);
            }
            // Reverse the read traversal; the stored target path stays unchanged.
            if (reverse_relative_to_path) {
                reverse(nodes.begin(), nodes.end());
                reverse(steps.begin(), steps.end());
                for (auto& node : nodes) node = graph.flip(node);
            }
            string sequence;
            for (auto node : nodes) sequence += graph.get_sequence(node);
            bdsg::PositionOverlay pos_graph(&graph);
            TestSurjector surjector(&pos_graph);
            surjector.prune_suspicious_anchors = false;
            surjector.prune_tail_region_anchors = true;
            for (bool prune_left : {false, true}) {
                INFO("middle node locally reversed = " << reverse_middle_node
                     << ", read reversed relative to target path = " << reverse_relative_to_path
                     << ", pruning = " << (prune_left ? "left" : "right"));
                vector<Surjector::path_chunk_t> chunks(3);
                vector<pair<step_handle_t, step_handle_t>> ranges;
                for (size_t i = 0; i < nodes.size(); ++i) {
                    chunks[i].first = {sequence.begin() + 4 * i, sequence.begin() + 4 * (i + 1)};
                    auto* mapping = chunks[i].second.add_mapping();
                    mapping->mutable_position()->set_node_id(graph.get_id(nodes[i]));
                    mapping->mutable_position()->set_is_reverse(graph.get_is_reverse(nodes[i]));
                    auto* edit = mapping->add_edit();
                    edit->set_from_length(4);
                    edit->set_to_length(4);
                    ranges.emplace_back(steps[i], steps[i]);
                }
                surjector.prune_and_trim_anchors(sequence, chunks, ranges,
                                                prune_left ? 4 : 3, prune_left ? 3 : 4);
                REQUIRE(chunks.size() == 2);
                REQUIRE(ranges.size() == 2);
                for (size_t i = 0; i < chunks.size(); ++i) {
                    const size_t original = i + (prune_left ? 1 : 0);
                    CHECK(chunks[i].first.first == sequence.begin() + 4 * original);
                    CHECK(chunks[i].second.mapping(0).position().node_id() == graph.get_id(nodes[original]));
                    CHECK(ranges[i].first == steps[original]);
                }
            }
        }
    }
}

TEST_CASE( "Spliced surject algorithm preserves deletions against the path", "[surject]" ) {
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GTCGT");
    handle_t h2 = graph.create_handle("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    handle_t h3 = graph.create_handle("TCCTTGC");
    handle_t h4 = graph.create_handle("A");
    handle_t h5 = graph.create_handle("T");
    handle_t h6 = graph.create_handle("GCCGA");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h3, h5);
    graph.create_edge(h4, h6);
    graph.create_edge(h5, h6);
    
    path_handle_t p = graph.create_path_handle("p");
    graph.append_step(p, h1);
    graph.append_step(p, h2);
    graph.append_step(p, h3);
    graph.append_step(p, h4);
    graph.append_step(p, h6);
    
    bdsg::PositionOverlay pos_graph(&graph);
    Surjector surjector(&pos_graph);
    
    vector<handle_t> read_path{h1, h3, h5, h6};
    
    Alignment read;
    string seq;
    Path* rpath = read.mutable_path();
    for (handle_t h : read_path) {
        Mapping* m = rpath->add_mapping();
        m->set_rank(rpath->mapping_size());
        m->mutable_position()->set_node_id(pos_graph.get_id(h));
        Edit* e = m->add_edit();
        e->set_from_length(pos_graph.get_length(h));
        e->set_to_length(pos_graph.get_length(h));
        
        seq += pos_graph.get_sequence(h);
    }
    read.set_sequence(seq);
    
    read.set_score(Aligner().scorer->score_contiguous_alignment(read));
    
    unordered_set<path_handle_t> paths{p};
    vector<Alignment> surjected_alns = surjector.surject(read, paths, true, true);
    REQUIRE(surjected_alns.size() == 1);
    auto& surjected = surjected_alns.front();
    
    vector<handle_t> surjected_path{h1, h3, h4, h6};
    
    REQUIRE(surjected.path().mapping_size() == read_path.size());
    
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(surjected.path().mapping(i).position().node_id() == graph.get_id(surjected_path[i]));
        REQUIRE(!surjected.path().mapping(i).position().is_reverse());
    }
    
    // should not penalize long deletions (assumed to be splices)
    REQUIRE(surjected.score() == read.score() - Aligner().scorer->mismatch - Aligner().scorer->match);
    REQUIRE(surjected.refpos_size() == 1);
    REQUIRE(surjected.refpos(0).name() == graph.get_path_name(p));
    REQUIRE(surjected.refpos(0).offset() == 0);
    
    
    
    Alignment rev_read;
    rev_read.set_sequence(reverse_complement(seq));
    
    Path* rev_rpath = rev_read.mutable_path();
    for (size_t i = 0; i < read_path.size(); ++i) {
        handle_t h = read_path[read_path.size() - i - 1];
        Mapping* m = rev_rpath->add_mapping();
        m->set_rank(rev_rpath->mapping_size());
        m->mutable_position()->set_node_id(pos_graph.get_id(h));
        m->mutable_position()->set_is_reverse(true);
        Edit* e = m->add_edit();
        e->set_from_length(pos_graph.get_length(h));
        e->set_to_length(pos_graph.get_length(h));
    }
    
    rev_read.set_score(Aligner().scorer->score_contiguous_alignment(rev_read));
    
    vector<Alignment> rev_surjected_alns = surjector.surject(rev_read, paths, true, true);
    REQUIRE(rev_surjected_alns.size() == 1);
    auto& rev_surjected = rev_surjected_alns.front();
    
    REQUIRE(rev_surjected.path().mapping_size() == read_path.size());
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(rev_surjected.path().mapping(i).position().node_id()
                == graph.get_id(surjected_path[surjected_path.size() - i - 1]));
        REQUIRE(rev_surjected.path().mapping(i).position().is_reverse());
    }
    
    // should not penalize long deletions (assumed to be splices)
    REQUIRE(rev_surjected.score() == read.score() - Aligner().scorer->mismatch - Aligner().scorer->match);
    REQUIRE(rev_surjected.refpos_size() == 1);
    REQUIRE(rev_surjected.refpos(0).name() == graph.get_path_name(p));
    REQUIRE(rev_surjected.refpos(0).offset() == 0);
    
}

TEST_CASE( "Spliced surject algorithm works when a read touches the same path in both orientations", "[surject]" ) {
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GGGGGGGGGGGGGGG");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("C");
    handle_t h4 = graph.create_handle("TTTTTTTTT");
    handle_t h5 = graph.create_handle("AAA");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, graph.flip(h4));
    graph.create_edge(h5, h4);
    graph.create_edge(h4, h4);
    graph.create_edge(graph.flip(h5), h1);
    
    path_handle_t p = graph.create_path_handle("p");
    graph.append_step(p, h1);
    graph.append_step(p, h2);
    graph.append_step(p, h4);
    graph.append_step(p, h4);
    graph.append_step(p, graph.flip(h4));
    graph.append_step(p, graph.flip(h5));
    
    bdsg::PositionOverlay pos_graph(&graph);
    Surjector surjector(&pos_graph);
    
    vector<handle_t> read_path{h1, h3, h4, h4};
    
    Alignment read;
    string seq;
    Path* rpath = read.mutable_path();
    for (handle_t h : read_path) {
        Mapping* m = rpath->add_mapping();
        m->set_rank(rpath->mapping_size());
        m->mutable_position()->set_node_id(pos_graph.get_id(h));
        Edit* e = m->add_edit();
        e->set_from_length(pos_graph.get_length(h));
        e->set_to_length(pos_graph.get_length(h));
        
        seq += pos_graph.get_sequence(h);
    }
    read.set_sequence(seq);
    
    read.set_score(Aligner().scorer->score_contiguous_alignment(read));
    
    unordered_set<path_handle_t> paths{p};
    vector<Alignment> surjected_alns = surjector.surject(read, paths, true, true);
    REQUIRE(surjected_alns.size() == 1);
    auto& surjected = surjected_alns.front();
    
    vector<handle_t> surjected_path{h1, h2, h4, h4};
    
    REQUIRE(surjected.path().mapping_size() == read_path.size());
    
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(surjected.path().mapping(i).position().node_id() == graph.get_id(surjected_path[i]));
        REQUIRE(surjected.path().mapping(i).position().is_reverse() == graph.get_is_reverse(surjected_path[i]));
    }
    
    Alignment rev_read;
    rev_read.set_sequence(reverse_complement(seq));
    
    Path* rev_rpath = rev_read.mutable_path();
    for (size_t i = 0; i < read_path.size(); ++i) {
        handle_t h = read_path[read_path.size() - i - 1];
        Mapping* m = rev_rpath->add_mapping();
        m->set_rank(rev_rpath->mapping_size());
        m->mutable_position()->set_node_id(pos_graph.get_id(h));
        m->mutable_position()->set_is_reverse(true);
        Edit* e = m->add_edit();
        e->set_from_length(pos_graph.get_length(h));
        e->set_to_length(pos_graph.get_length(h));
    }
    
    rev_read.set_score(Aligner().scorer->score_contiguous_alignment(rev_read));
    
    vector<Alignment> rev_surjected_alns = surjector.surject(rev_read, paths, true, true);
    REQUIRE(rev_surjected_alns.size() == 1);
    auto& rev_surjected = rev_surjected_alns.front();
    
    REQUIRE(rev_surjected.path().mapping_size() == read_path.size());
    for (size_t i = 0; i < surjected_path.size(); ++i) {
        REQUIRE(rev_surjected.path().mapping(i).position().node_id()
                == graph.get_id(surjected_path[surjected_path.size() - i - 1]));
        REQUIRE(rev_surjected.path().mapping(i).position().is_reverse()
                == !graph.get_is_reverse(surjected_path[surjected_path.size() - i - 1]));
    }
}

TEST_CASE("Path overlapping segments can be identified from multipath alignment",
          "[surject][multipath]"){
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("GTCGT");
    handle_t h2 = graph.create_handle("A");
    handle_t h3 = graph.create_handle("T");
    handle_t h4 = graph.create_handle("TTAGAC");
    handle_t h5 = graph.create_handle("GCA");
    handle_t h6 = graph.create_handle("ATTAGACGCA");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
    graph.create_edge(h4, h5);
    graph.create_edge(h5, h6);
    
    path_handle_t p = graph.create_path_handle("p");
    step_handle_t st0 = graph.append_step(p, h1);
    step_handle_t st1 = graph.append_step(p, h2);
    step_handle_t st2 = graph.append_step(p, h4);
    step_handle_t st3 = graph.append_step(p, h5);
    step_handle_t st4 = graph.append_step(p, h6);
    
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    
    multipath_alignment_t mp_aln;
    mp_aln.set_sequence("CGTATTAGACGC");
    
    auto s0 = mp_aln.add_subpath();
    auto m0 = s0->mutable_path()->add_mapping();
    m0->mutable_position()->set_node_id(graph.get_id(h1));
    m0->mutable_position()->set_is_reverse(false);
    m0->mutable_position()->set_offset(2);
    auto e0 = m0->add_edit();
    e0->set_from_length(3);
    e0->set_to_length(3);
    e0->set_sequence("");
    
    s0->add_next(1);
    s0->add_next(2);
    
    auto s1 = mp_aln.add_subpath();
    auto m1 = s1->mutable_path()->add_mapping();
    m1->mutable_position()->set_node_id(graph.get_id(h2));
    m1->mutable_position()->set_is_reverse(false);
    m1->mutable_position()->set_offset(0);
    auto e1 = m1->add_edit();
    e1->set_from_length(1);
    e1->set_to_length(1);
    e1->set_sequence("");
    
    s1->add_next(3);
    
    auto s2 = mp_aln.add_subpath();
    auto m2 = s2->mutable_path()->add_mapping();
    m2->mutable_position()->set_node_id(graph.get_id(h3));
    m2->mutable_position()->set_is_reverse(false);
    m2->mutable_position()->set_offset(0);
    auto e2 = m2->add_edit();
    e2->set_from_length(1);
    e2->set_to_length(1);
    e2->set_sequence("A");
    
    s2->add_next(3);
    
    auto s3 = mp_aln.add_subpath();
    auto m3 = s3->mutable_path()->add_mapping();
    m3->mutable_position()->set_node_id(graph.get_id(h4));
    m3->mutable_position()->set_is_reverse(false);
    m3->mutable_position()->set_offset(0);
    auto e3 = m3->add_edit();
    e3->set_from_length(6);
    e3->set_to_length(6);
    e3->set_sequence("");
    auto m4 = s3->mutable_path()->add_mapping();
    m4->mutable_position()->set_node_id(graph.get_id(h5));
    m4->mutable_position()->set_is_reverse(false);
    m4->mutable_position()->set_offset(0);
    auto e4 = m4->add_edit();
    e4->set_from_length(2);
    e4->set_to_length(2);
    e4->set_sequence("");
    
    identify_start_subpaths(mp_aln);
    
    unordered_set<path_handle_t> surjection_paths{p};
    
    function<int64_t(int64_t)> node_length = [&](int64_t node_id) {
        return int64_t(graph.get_length(graph.get_handle(node_id)));
    };
    
    SECTION("Forward strand of path") {
        
        unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>> connections;
        
        auto overlaps = surjector.extract_overlapping_paths(&pos_graph, mp_aln,
                                                            surjection_paths,
                                                            connections);
        
        auto fp = make_pair(p, false);
        
        REQUIRE(overlaps.count(fp));
        REQUIRE(overlaps.size() == 1);
        REQUIRE(connections.empty());
        
        auto& p_overlaps = overlaps[fp];
        
        REQUIRE(p_overlaps.first.size() == 1);
        REQUIRE(p_overlaps.first.front().first.first == mp_aln.sequence().begin());
        REQUIRE(p_overlaps.first.front().first.second == mp_aln.sequence().end());
        REQUIRE(p_overlaps.first.front().second.mapping_size() == 4);
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().node_id() == graph.get_id(h1));
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().is_reverse() == false);
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().offset() == 2);
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().node_id() == graph.get_id(h2));
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().is_reverse() == false);
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().offset() == 0);
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().node_id() == graph.get_id(h4));
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().is_reverse() == false);
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().offset() == 0);
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().node_id() == graph.get_id(h5));
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().is_reverse() == false);
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().offset() == 0);
        
        REQUIRE(p_overlaps.second.size() == 1);
        REQUIRE(p_overlaps.second.front().first == st0);
        REQUIRE(p_overlaps.second.front().second == st3);
    }
    
    SECTION("Reverse strand of path"){
        
        unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>> connections;
        
        multipath_alignment_t rev_mp_aln;
        rev_comp_multipath_alignment(mp_aln, node_length, rev_mp_aln);
        
        auto overlaps = surjector.extract_overlapping_paths(&pos_graph, rev_mp_aln,
                                                            surjection_paths,
                                                            connections);
        
        auto rp = make_pair(p, true);
            
        REQUIRE(overlaps.count(rp));
        REQUIRE(overlaps.size() == 1);
        REQUIRE(connections.empty());
        
        auto& p_overlaps = overlaps[rp];
        
        REQUIRE(p_overlaps.first.size() == 1);
        REQUIRE(p_overlaps.first.front().first.first == rev_mp_aln.sequence().begin());
        REQUIRE(p_overlaps.first.front().first.second == rev_mp_aln.sequence().end());
        REQUIRE(p_overlaps.first.front().second.mapping_size() == 4);
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().node_id() == graph.get_id(h5));
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().is_reverse() == true);
        REQUIRE(p_overlaps.first.front().second.mapping(0).position().offset() == 1);
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().node_id() == graph.get_id(h4));
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().is_reverse() == true);
        REQUIRE(p_overlaps.first.front().second.mapping(1).position().offset() == 0);
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().node_id() == graph.get_id(h2));
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().is_reverse() == true);
        REQUIRE(p_overlaps.first.front().second.mapping(2).position().offset() == 0);
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().node_id() == graph.get_id(h1));
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().is_reverse() == true);
        REQUIRE(p_overlaps.first.front().second.mapping(3).position().offset() == 0);
        
        REQUIRE(p_overlaps.second.size() == 1);
        REQUIRE(p_overlaps.second.front().first == st3);
        REQUIRE(p_overlaps.second.front().second == st0);
    }
    
    //ATTAGACGCA
    auto s4 = mp_aln.add_subpath();
    auto m5 = s4->mutable_path()->add_mapping();
    m5->mutable_position()->set_node_id(graph.get_id(h6));
    m5->mutable_position()->set_is_reverse(false);
    m5->mutable_position()->set_offset(0);
    auto e5 = m5->add_edit();
    e5->set_from_length(10);
    e5->set_to_length(10);
    e5->set_sequence("");
    
    auto s5 = mp_aln.add_subpath();
    auto m6 = s5->mutable_path()->add_mapping();
    m6->mutable_position()->set_node_id(graph.get_id(h6));
    m6->mutable_position()->set_is_reverse(false);
    m6->mutable_position()->set_offset(1);
    auto e6 = m6->add_edit();
    e6->set_from_length(9);
    e6->set_to_length(9);
    e6->set_sequence("");
    
    auto c0 = mp_aln.mutable_subpath(0)->add_connection();
    c0->set_next(4);
    c0->set_score(-2);
    auto c1 = mp_aln.mutable_subpath(1)->add_connection();
    c1->set_next(5);
    c1->set_score(-1);
    
    SECTION("Connections break segments and are recorded correctly") {
        
        unordered_map<pair<path_handle_t, bool>, vector<tuple<size_t, size_t, int32_t>>> connections;
        
        auto overlaps = surjector.extract_overlapping_paths(&pos_graph, mp_aln,
                                                            surjection_paths,
                                                            connections);
        
        auto fp = make_pair(p, false);
        
        REQUIRE(overlaps.count(fp));
        REQUIRE(overlaps.size() == 1);
        
        auto& p_overlaps = overlaps[fp];
        
        REQUIRE(p_overlaps.first.size() == 5);
        REQUIRE(p_overlaps.second.size() == 5);
        
        REQUIRE(connections.size() == 1);
        REQUIRE(connections.count(fp));
        
        auto& p_connections = connections[fp];
        
        REQUIRE(p_connections.size() == 2);
        
        for (auto& connection : p_connections) {
            if (get<2>(connection) == -2) {
                REQUIRE(p_overlaps.first[get<0>(connection)].second.mapping(0).position().node_id() == graph.get_id(h1));
                REQUIRE(p_overlaps.first[get<1>(connection)].second.mapping(0).position().node_id() == graph.get_id(h6));
                REQUIRE(p_overlaps.first[get<1>(connection)].second.mapping(0).position().offset() == 0);

            }
            else if (get<2>(connection) == -1) {
                REQUIRE(p_overlaps.first[get<0>(connection)].second.mapping(0).position().node_id() == graph.get_id(h2));
                REQUIRE(p_overlaps.first[get<1>(connection)].second.mapping(0).position().node_id() == graph.get_id(h6));
                REQUIRE(p_overlaps.first[get<1>(connection)].second.mapping(0).position().offset() == 1);
            }
            else {
                REQUIRE(false);
            }
        }
    }
}

TEST_CASE("Multipath alignments can be surjected", "[surject][multipath]") {

    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("A");
    handle_t h2 = graph.create_handle("T");
    handle_t h3 = graph.create_handle("TTAGAC");
    handle_t h4 = graph.create_handle("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    handle_t h5 = graph.create_handle("GCA");
    handle_t h6 = graph.create_handle("G");
    handle_t h7 = graph.create_handle("T");
    handle_t h8 = graph.create_handle("TTAGAC");
    
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h3);
    graph.create_edge(h3, h4);
    graph.create_edge(h3, h5);
    graph.create_edge(h4, h5);
    graph.create_edge(h5, h6);
    graph.create_edge(h5, h7);
    graph.create_edge(h6, h8);
    graph.create_edge(h7, h8);
    
    path_handle_t p = graph.create_path_handle("p");
    step_handle_t st0 = graph.append_step(p, h1);
    step_handle_t st1 = graph.append_step(p, h3);
    step_handle_t st2 = graph.append_step(p, h4);
    step_handle_t st3 = graph.append_step(p, h5);
    step_handle_t st4 = graph.append_step(p, h6);
    step_handle_t st5 = graph.append_step(p, h8);
    
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    
    multipath_alignment_t mp_aln;
    mp_aln.set_sequence("TTTAGACGCTAGA");
    
    auto s0 = mp_aln.add_subpath();
    auto m0 = s0->mutable_path()->add_mapping();
    m0->mutable_position()->set_node_id(graph.get_id(h1));
    m0->mutable_position()->set_is_reverse(false);
    m0->mutable_position()->set_offset(0);
    auto e0 = m0->add_edit();
    e0->set_from_length(1);
    e0->set_to_length(1);
    e0->set_sequence("T");
    
    s0->set_score(1);
    s0->add_next(2);
    
    auto s1 = mp_aln.add_subpath();
    auto m1 = s1->mutable_path()->add_mapping();
    m1->mutable_position()->set_node_id(graph.get_id(h2));
    m1->mutable_position()->set_is_reverse(false);
    m1->mutable_position()->set_offset(0);
    auto e1 = m1->add_edit();
    e1->set_from_length(1);
    e1->set_to_length(1);
    e1->set_sequence("");
    
    s1->set_score(6);
    s1->add_next(2);
    
    auto s2 = mp_aln.add_subpath();
    auto m2 = s2->mutable_path()->add_mapping();
    m2->mutable_position()->set_node_id(graph.get_id(h3));
    m2->mutable_position()->set_is_reverse(false);
    m2->mutable_position()->set_offset(0);
    auto e2 = m2->add_edit();
    e2->set_from_length(6);
    e2->set_to_length(6);
    e2->set_sequence("");
    
    s2->set_score(6);
    s2->add_next(3);
    
    auto s3 = mp_aln.add_subpath();
    auto m3 = s3->mutable_path()->add_mapping();
    m3->mutable_position()->set_node_id(graph.get_id(h5));
    m3->mutable_position()->set_is_reverse(false);
    m3->mutable_position()->set_offset(0);
    auto e3 = m3->add_edit();
    e3->set_from_length(2);
    e3->set_to_length(2);
    e3->set_sequence("");
    
    s3->set_score(2);
    auto c0 = s3->add_connection();
    c0->set_next(4);
    c0->set_score(-1);
    
    auto s4 = mp_aln.add_subpath();
    auto m4 = s4->mutable_path()->add_mapping();
    m4->mutable_position()->set_node_id(graph.get_id(h8));
    m4->mutable_position()->set_is_reverse(false);
    m4->mutable_position()->set_offset(1);
    auto e4 = m4->add_edit();
    e4->set_from_length(4);
    e4->set_to_length(4);
    e4->set_sequence("");
    
    s4->set_score(9);
    
    identify_start_subpaths(mp_aln);
    
    multipath_alignment_t rc_mp_aln;
    rev_comp_multipath_alignment(mp_aln, [&](int64_t node_id) {
        return (int64_t) graph.get_length(graph.get_handle(node_id));
    }, rc_mp_aln);
    
    SECTION("With deletion splices allowed") {
        
        string path_name;
        int64_t path_pos;
        bool path_rev;
        unordered_set<path_handle_t> paths{p};
        vector<tuple<string, int64_t, bool>> positions;
        auto surjected_alns = surjector.surject(mp_aln, paths, positions, true, true);
        REQUIRE(surjected_alns.size() == 1);
        auto& surjected = surjected_alns.front();
        tie(path_name, path_pos, path_rev) = positions.front();
        
        REQUIRE(path_name == graph.get_path_name(p));
        REQUIRE(path_rev == false);
        REQUIRE(path_pos == 0);
        
        // normalize the representation
        merge_non_branching_subpaths(surjected);
        
        REQUIRE(surjected.sequence() == mp_aln.sequence());
        REQUIRE(surjected.subpath_size() == 2);
        REQUIRE(surjected.subpath(0).score() == 1 + 6 + 2);
        REQUIRE(surjected.subpath(0).path().mapping_size() == 3);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h1));
        REQUIRE(surjected.subpath(0).path().mapping(0).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(1).position().node_id() == graph.get_id(h3));
        REQUIRE(surjected.subpath(0).path().mapping(1).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(1).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(2).position().node_id() == graph.get_id(h5));
        REQUIRE(surjected.subpath(0).path().mapping(2).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(2).position().offset() == 0);
        REQUIRE(surjected.subpath(0).next_size() == 0);
        REQUIRE(surjected.subpath(0).connection_size() == 1);
        REQUIRE(surjected.subpath(0).connection(0).next() == 1);
        REQUIRE(surjected.subpath(0).connection(0).score() == -1);
        REQUIRE(surjected.subpath(1).score() == 9);
        REQUIRE(surjected.subpath(1).path().mapping_size() == 1);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().node_id() == graph.get_id(h8));
        REQUIRE(surjected.subpath(1).path().mapping(0).position().is_reverse() == false);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().offset() == 1);
        REQUIRE(surjected.subpath(1).next_size() == 0);
        REQUIRE(surjected.subpath(1).connection_size() == 0);
    }
    
    SECTION("Reverse with deletion splices allowed") {
        
        string path_name;
        int64_t path_pos;
        bool path_rev;
        unordered_set<path_handle_t> paths{p};
        vector<tuple<string, int64_t, bool>> positions;
        auto surjected_alns = surjector.surject(rc_mp_aln, paths, positions, true, true);
        REQUIRE(surjected_alns.size() == 1);
        auto& surjected = surjected_alns.front();
        tie(path_name, path_pos, path_rev) = positions.front();
        
        REQUIRE(path_name == graph.get_path_name(p));
        REQUIRE(path_rev == true);
        REQUIRE(path_pos == 0);
        
        // normalize the representation
        merge_non_branching_subpaths(surjected);
        
        REQUIRE(surjected.sequence() == rc_mp_aln.sequence());
        REQUIRE(surjected.subpath_size() == 2);
        REQUIRE(surjected.subpath(0).score() == 9);
        REQUIRE(surjected.subpath(0).path().mapping_size() == 1);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h8));
        REQUIRE(surjected.subpath(0).path().mapping(0).position().is_reverse() == true);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().offset() == 1);
        REQUIRE(surjected.subpath(0).connection_size() == 1);
        REQUIRE(surjected.subpath(0).connection(0).next() == 1);
        REQUIRE(surjected.subpath(0).connection(0).score() == -1);
        REQUIRE(surjected.subpath(0).score() == 1 + 6 + 2);
        REQUIRE(surjected.subpath(1).path().mapping_size() == 3);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().node_id() == graph.get_id(h5));
        REQUIRE(surjected.subpath(1).path().mapping(0).position().is_reverse() == true);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().offset() == 1);
        REQUIRE(surjected.subpath(1).path().mapping(1).position().node_id() == graph.get_id(h3));
        REQUIRE(surjected.subpath(1).path().mapping(1).position().is_reverse() == true);
        REQUIRE(surjected.subpath(1).path().mapping(1).position().offset() == 0);
        REQUIRE(surjected.subpath(1).path().mapping(2).position().node_id() == graph.get_id(h1));
        REQUIRE(surjected.subpath(1).path().mapping(2).position().is_reverse() == true);
        REQUIRE(surjected.subpath(1).path().mapping(2).position().offset() == 0);
        REQUIRE(surjected.subpath(1).next_size() == 0);
        REQUIRE(surjected.subpath(1).next_size() == 0);
        REQUIRE(surjected.subpath(1).connection_size() == 0);
    }
    
    SECTION("Without deletion splices allowed") {
        
        string path_name;
        int64_t path_pos;
        bool path_rev;
        unordered_set<path_handle_t> paths{p};
        vector<tuple<string, int64_t, bool>> positions;
        auto surjected_alns = surjector.surject(mp_aln, paths, positions, true, false);
        REQUIRE(surjected_alns.size() == 1);
        auto& surjected = surjected_alns.front();
        tie(path_name, path_pos, path_rev) = positions.front();
        
        REQUIRE(path_name == graph.get_path_name(p));
        REQUIRE(path_rev == false);
        REQUIRE(path_pos == 0);
        
        // normalize the representation
        merge_non_branching_subpaths(surjected);
        
        REQUIRE(surjected.sequence() == mp_aln.sequence());
        REQUIRE(surjected.subpath_size() == 2);
        REQUIRE(surjected.subpath(0).score() == 1 + 6 + 2 - 6 - 32);
        REQUIRE(surjected.subpath(0).path().mapping_size() == 4);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().node_id() == graph.get_id(h1));
        REQUIRE(surjected.subpath(0).path().mapping(0).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(0).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(1).position().node_id() == graph.get_id(h3));
        REQUIRE(surjected.subpath(0).path().mapping(1).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(1).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(2).position().node_id() == graph.get_id(h4));
        REQUIRE(surjected.subpath(0).path().mapping(2).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(2).position().offset() == 0);
        REQUIRE(surjected.subpath(0).path().mapping(3).position().node_id() == graph.get_id(h5));
        REQUIRE(surjected.subpath(0).path().mapping(3).position().is_reverse() == false);
        REQUIRE(surjected.subpath(0).path().mapping(3).position().offset() == 0);
        REQUIRE(surjected.subpath(0).next_size() == 0);
        REQUIRE(surjected.subpath(0).connection_size() == 1);
        REQUIRE(surjected.subpath(0).connection(0).next() == 1);
        REQUIRE(surjected.subpath(0).connection(0).score() == -1);
        REQUIRE(surjected.subpath(1).score() == 9);
        REQUIRE(surjected.subpath(1).path().mapping_size() == 1);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().node_id() == graph.get_id(h8));
        REQUIRE(surjected.subpath(1).path().mapping(0).position().is_reverse() == false);
        REQUIRE(surjected.subpath(1).path().mapping(0).position().offset() == 1);
        REQUIRE(surjected.subpath(1).next_size() == 0);
        REQUIRE(surjected.subpath(1).connection_size() == 0);
    }
}

TEST_CASE("Duplicate path chunks can be detected", "[surject][multipath]") {
    
    Alignment aln;
    aln.set_sequence("ACGT");
    
    Surjector::path_chunk_t chunk1;
    Surjector::path_chunk_t chunk2;
    Surjector::path_chunk_t chunk3;
    Surjector::path_chunk_t chunk4;
    
    chunk1.first.first = aln.sequence().begin();
    chunk1.first.second = aln.sequence().begin() + 4;
    auto m00 = chunk1.second.add_mapping();
    m00->mutable_position()->set_node_id(1);
    m00->mutable_position()->set_is_reverse(false);
    m00->mutable_position()->set_offset(2);
    auto e00 = m00->add_edit();
    e00->set_from_length(2);
    e00->set_to_length(2);
    auto m01 = chunk1.second.add_mapping();
    m01->mutable_position()->set_node_id(2);
    m01->mutable_position()->set_is_reverse(false);
    m01->mutable_position()->set_offset(0);
    auto e01 = m01->add_edit();
    e01->set_from_length(2);
    e01->set_to_length(2);
    
    chunk2.first.first = aln.sequence().begin();
    chunk2.first.second = aln.sequence().begin() + 2;
    auto m10 = chunk2.second.add_mapping();
    m10->mutable_position()->set_node_id(1);
    m10->mutable_position()->set_is_reverse(false);
    m10->mutable_position()->set_offset(2);
    auto e10 = m10->add_edit();
    e10->set_from_length(2);
    e10->set_to_length(2);
    
    chunk3.first.first = aln.sequence().begin() + 2;
    chunk3.first.second = aln.sequence().begin() + 4;
    auto m20 = chunk3.second.add_mapping();
    m20->mutable_position()->set_node_id(2);
    m20->mutable_position()->set_is_reverse(false);
    m20->mutable_position()->set_offset(0);
    auto e20 = m20->add_edit();
    e20->set_from_length(2);
    e20->set_to_length(2);
    
    chunk4.first.first = aln.sequence().begin();
    chunk4.first.second = aln.sequence().begin() + 4;
    auto m30 = chunk4.second.add_mapping();
    m30->mutable_position()->set_node_id(1);
    m30->mutable_position()->set_is_reverse(false);
    m30->mutable_position()->set_offset(2);
    auto e30 = m30->add_edit();
    e30->set_from_length(2);
    e30->set_to_length(2);
    auto m31 = chunk4.second.add_mapping();
    m31->mutable_position()->set_node_id(3);
    m31->mutable_position()->set_is_reverse(false);
    m31->mutable_position()->set_offset(0);
    auto e31 = m31->add_edit();
    e31->set_from_length(2);
    e31->set_to_length(2);
    
    bdsg::HashGraph graph;
    handle_t h1 = graph.create_handle("AAAC");
    handle_t h2 = graph.create_handle("GTGT");
    handle_t h3 = graph.create_handle("GTAC");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h1);
    
    path_handle_t p = graph.create_path_handle("path");
    
    step_handle_t s1 = graph.append_step(p, h1);
    step_handle_t s2 = graph.append_step(p, h2);
    step_handle_t s3 = graph.append_step(p, h1);
    step_handle_t s4 = graph.append_step(p, h3);
    
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    
    vector<Surjector::path_chunk_t> path_chunks{chunk1, chunk2, chunk3, chunk4};
    vector<pair<step_handle_t, step_handle_t>> ref_chunks;
    ref_chunks.emplace_back(s1, s2);
    ref_chunks.emplace_back(s1, s1);
    ref_chunks.emplace_back(s2, s2);
    ref_chunks.emplace_back(s3, s4);
    
    vector<tuple<size_t, size_t, int32_t>> connections;
    
    surjector.filter_redundant_path_chunks(false, path_chunks, ref_chunks, connections);
    
    REQUIRE(ref_chunks.size() == path_chunks.size());
    REQUIRE(path_chunks.size() == 2);
    
}

TEST_CASE("Supplementary alignments can be generated", "[surject]") {
    
    bdsg::HashGraph graph;

    path_handle_t p = graph.create_path_handle("p");

    handle_t h1 = graph.create_handle("GTCGT");
    graph.append_step(p, h1);

    handle_t prev = h1;
    for (size_t i = 0; i < 20; ++i) {
        handle_t h = graph.create_handle(string(64, 'A'));
        graph.create_edge(prev, h);
        graph.append_step(p, h);
        prev = h;
    }
    handle_t h2 = graph.create_handle("TCCTTGC");
    graph.create_edge(prev, h2);
    graph.append_step(p, h2);
    
    handle_t h3 = graph.create_handle("TGTC");

    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h3, h2);

    bdsg::PositionOverlay pos_graph(&graph);
    Surjector surjector(&pos_graph);
    surjector.report_supplementary = true;

    unordered_set<path_handle_t> paths{p};
    
    SECTION("In a single-path alignment") {
        // with and without an unaligned middle portion
        vector<vector<handle_t>> read_paths{{h1, h2}, {h1, h3, h2}};
        for (const auto& read_path : read_paths) {
            Alignment read;
            string seq;
            Path* rpath = read.mutable_path();
            for (handle_t h : read_path) {
                Mapping* m = rpath->add_mapping();
                m->set_rank(rpath->mapping_size());
                m->mutable_position()->set_node_id(pos_graph.get_id(h));
                Edit* e = m->add_edit();
                e->set_from_length(pos_graph.get_length(h));
                e->set_to_length(pos_graph.get_length(h));
                
                seq += pos_graph.get_sequence(h);
            }
            read.set_sequence(seq);
            
            read.set_score(Aligner().scorer->score_contiguous_alignment(read));
            
            vector<Alignment> surjected_alns = surjector.surject(read, paths, true, false);
        
            REQUIRE(surjected_alns.size() == 2);
            bool found1 = false, found2 = false;
            for (auto& aln : surjected_alns) {
                const auto& path = aln.path();
                REQUIRE(path.mapping_size() == 1);
                const auto& mapping = path.mapping(0);
                REQUIRE(mapping.edit_size() == 2); // match and soft-clip
                if (mapping.position().node_id() == graph.get_id(h1)) {
                    found1 = true;
                    const auto& clip = mapping.edit(1);
                    REQUIRE(clip.from_length() == 0);
                    REQUIRE(clip.sequence() == aln.sequence().substr(mapping.edit(0).to_length(), string::npos));
                }
                else if (mapping.position().node_id() == graph.get_id(h2)) {
                    found2 = true;
                    const auto& clip = mapping.edit(0);
                    REQUIRE(clip.from_length() == 0);
                    REQUIRE(clip.sequence() == aln.sequence().substr(0, aln.sequence().size() - mapping.edit(1).to_length()));
                }
            }
        
            REQUIRE(found1);
            REQUIRE(found2);
        
            bool suppl1 = is_supplementary(surjected_alns.front());
            bool suppl2 = is_supplementary(surjected_alns.back());
            REQUIRE(suppl1 != suppl2);
        }
    }

    SECTION("In a multipath alignment") {

        vector<handle_t> read_path{h1, h2};
        multipath_alignment_t mp_aln;
        string seq;
        for (handle_t h : read_path) {
            if (mp_aln.subpath_size() != 0) {
                mp_aln.mutable_subpath(mp_aln.subpath_size() - 1)->add_next(mp_aln.subpath_size());
            }
            auto subpath = mp_aln.add_subpath();
            subpath->set_score(pos_graph.get_length(h));
            auto mapping = subpath->mutable_path()->add_mapping();
            auto pos = mapping->mutable_position();
            pos->set_node_id(pos_graph.get_id(h));
            pos->set_is_reverse(pos_graph.get_is_reverse(h));
            pos->set_offset(0);
            auto edit = mapping->add_edit();
            edit->set_from_length(pos_graph.get_length(h));
            edit->set_to_length(pos_graph.get_length(h));
            seq += pos_graph.get_sequence(h);
        }
        mp_aln.add_start(0);
        mp_aln.set_sequence(seq);
        
        {
            // unspliced
            vector<tuple<string, int64_t, bool>> positions;
            vector<multipath_alignment_t> surjected = surjector.surject(mp_aln, paths, positions, true, false);
    
            REQUIRE(surjected.size() == 2);
            bool found1 = false, found2 = false;
            for (auto& surj : surjected) {
                REQUIRE(surj.subpath_size() == 1);
                const auto& path = surj.subpath(0).path();
                REQUIRE(path.mapping_size() == 1);
                const auto& mapping = path.mapping(0);
                REQUIRE(mapping.edit_size() == 2);
                const auto& pos = mapping.position();
                if (pos.node_id() == pos_graph.get_id(h1)) {
                    found1 = true;
                }
                else if (pos.node_id() == pos_graph.get_id(h2)) {
                    found2 = true;
                }
            }

            REQUIRE(found1);
            REQUIRE(found2);
        
            bool suppl1 = is_supplementary(surjected.front());
            bool suppl2 = is_supplementary(surjected.back());
            REQUIRE(suppl1 != suppl2);
        }

        {
            // spliced are not returned as supplementaries
            vector<tuple<string, int64_t, bool>> positions;
            vector<multipath_alignment_t> surjected = surjector.surject(mp_aln, paths, positions, true, true);
    
            REQUIRE(surjected.size() == 1);
            for (auto& surj : surjected) {
                REQUIRE(surj.subpath_size() == 2);
                for (const auto& subpath : surj.subpath()) {
                    REQUIRE(subpath.path().mapping_size() == 1);
                }
                REQUIRE(surj.subpath(0).path().mapping(0).position().node_id() == pos_graph.get_id(h1));
                REQUIRE(surj.subpath(1).path().mapping(0).position().node_id() == pos_graph.get_id(h2));
            }
        }
    }
}

}
}
