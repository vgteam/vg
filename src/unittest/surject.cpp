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
    using Surjector::choose_primary;
    using Surjector::choose_primary_strand;
    using Surjector::add_SA_tag;
    
};


TEST_CASE("Diploid surjection validates and jointly selects read placements",
          "[surject][diploid]") {
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
    auto placement = [&](handle_t node, bool secondary) {
        Alignment aln;
        aln.set_name("read");
        aln.set_sequence(sequence);
        aln.set_mapping_quality(37);
        aln.set_is_secondary(secondary);
        auto* mapping = aln.mutable_path()->add_mapping();
        mapping->set_rank(1);
        mapping->mutable_position()->set_node_id(graph.get_id(node));
        auto* edit = mapping->add_edit();
        edit->set_from_length(sequence.size());
        edit->set_to_length(sequence.size());
        return aln;
    };
    auto primary = placement(b, false);
    auto alternative = placement(a, true);
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
            CHECK(aln.mapping_quality() < 10);
        }
        std::reverse(input.begin(), input.end());
        auto reversed = surjector.surject_diploid(input, {path_b, path_a});
        REQUIRE(reversed.size() == output.size());
        for (size_t i = 0; i < output.size(); ++i) {
            CHECK(reversed[i].refpos(0).name() == output[i].refpos(0).name());
            CHECK(reversed[i].mapping_quality() == output[i].mapping_quality());
        }
    }
    SECTION("A mapped secondary wins over an unmapped input primary") {
        primary.clear_path();
        auto output = surjector.surject_diploid({alternative, primary}, {path_a});
        REQUIRE(output.size() == 1);
        CHECK_FALSE(output.front().is_secondary());
        CHECK(output.front().mapping_quality() == 37);
        CHECK(output.front().refpos(0).name() == "A");
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
    SECTION("Identical graph paths do not create extra competitors") {
        alternative = primary;
        alternative.set_is_secondary(true);
        auto output = surjector.surject_diploid({alternative, primary}, {path_b});
        REQUIRE(output.size() == 1);
        CHECK(output.front().mapping_quality() == 37);
    }
    SECTION("One graph placement compares all target haplotype paths") {
        auto other = graph.create_path_handle("C");
        graph.append_step(other, b);
        bdsg::PositionOverlay updated_overlay(&graph);
        TestSurjector updated(&updated_overlay);
        updated.prune_suspicious_anchors = false;
        auto output = updated.surject_diploid({primary}, {path_b, other});
        REQUIRE(output.size() == 2);
        CHECK(get_annotation<bool>(output[0], "diploid_haplotype_preferred"));
        CHECK_FALSE(get_annotation<bool>(output[1], "diploid_haplotype_preferred"));
        CHECK(get_annotation<double>(output[0], "diploid_haplotype_quality") < 10);
    }
    SECTION("All-unmapped output preserves source quality and unrelated metadata") {
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
    SECTION("Empty input returns no records") {
        CHECK(surjector.surject_diploid({}, {path_a}).empty());
    }
    SECTION("Malformed groups are rejected") {
        CHECK_THROWS_AS(surjector.surject_diploid({alternative}, {path_a}), invalid_argument);
        alternative.set_is_secondary(false);
        CHECK_THROWS_AS(surjector.surject_diploid({primary, alternative}, {path_a}), invalid_argument);
        alternative.set_is_secondary(true);
        alternative.set_name("different");
        CHECK_THROWS_AS(surjector.surject_diploid({primary, alternative}, {path_a}), invalid_argument);
        alternative = placement(a, true);
        alternative.set_sequence("ACGT");
        CHECK_THROWS_AS(surjector.surject_diploid({primary, alternative}, {path_a}), invalid_argument);
        alternative = placement(a, true);
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

TEST_CASE("Diploid supplementary pieces retain their candidate and final SA qualities",
          "[surject][diploid]") {
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
    auto output = surjector.surject_diploid(input, paths);
    REQUIRE(output.size() == 4);
    for (size_t i = 0; i < output.size(); ++i) {
        const auto& aln = output[i];
        const string path = i < 2 ? "A" : "B";
        CHECK(aln.refpos(0).name() == path);
        CHECK(aln.is_secondary() == (i >= 2));
        CHECK(is_supplementary(aln) == (i % 2 == 1));
        CHECK(get_annotation<double>(aln, "diploid_source_mapping_quality") == 17);
        REQUIRE(has_annotation(aln, "tags"));
        const string tags = get_annotation<string>(aln, "tags");
        CHECK(tags.find("ZZ:Z:keep") != string::npos);
        CHECK(tags.find("obsolete") == string::npos);
        CHECK(tags.find("SA:Z:" + path + ",") != string::npos);
        CHECK(std::count(tags.begin(), tags.end(), ';') == 1);
        CHECK(tags.find("," + to_string(aln.mapping_quality()) + ",0;") != string::npos);
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

    SECTION("Supplementary reporting retains the disjoint candidate") {
        surjector.report_supplementary = true;
    }
    SECTION("Without supplementary reporting the disjoint candidate is omitted") {
        surjector.report_supplementary = false;
    }

    surjector.choose_primary(candidates);
    REQUIRE(candidates.size() == (surjector.report_supplementary ? 3 : 2));
    size_t primary_count = 0;
    for (const auto& candidate : candidates) {
        const auto& aln = candidate.first;
        const bool supplementary = has_annotation(aln, "supplementary")
            && get_annotation<bool>(aln, "supplementary");
        CHECK(aln.is_secondary() == (aln.name() == "A"));
        CHECK(supplementary == (aln.name() == "C"));
        if (!aln.is_secondary() && !supplementary) {
            CHECK(aln.name() == "B");
            CHECK(aln.score() == 60);
            ++primary_count;
        }
    }
    CHECK(primary_count == 1);

    if (surjector.report_supplementary) {
        vector<Alignment> output;
        vector<tuple<string, int64_t, bool>> positions;
        for (const auto& candidate : candidates) {
            output.push_back(candidate.first);
            positions.emplace_back("ref", candidate.first.path().mapping(0).position().offset(), false);
        }
        surjector.add_SA_tag(output, positions, pos_graph, false);
        for (const auto& aln : output) {
            if (aln.name() == "A") {
                CHECK_FALSE(has_annotation(aln, "tags"));
            } else {
                CHECK(has_annotation(aln, "tags"));
                if (has_annotation(aln, "tags")) {
                    const string expected = aln.name() == "B"
                        ? "SA:Z:ref,161,+,60S40M,0,0;"
                        : "SA:Z:ref,101,+,60M40S,0,0;";
                    CHECK(get_annotation<string>(aln, "tags") == expected);
                }
            }
        }
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

    path_handle_t expected_primary;
    bool expect_supplementary = false;
    SECTION("Secondary alternatives cannot increase the strand score") {
        add_candidate(path_a, step_a, 0, 100, 110, 0);
        add_candidate(path_b, step_b, 0, 100, 105, 0);
        add_candidate(path_b, step_b, 0, 100, 105, 200);
        add_candidate(path_b, step_b, 0, 100, 105, 400);
        expected_primary = path_a;
    }
    SECTION("Secondary alternatives cannot decrease the strand score") {
        add_candidate(path_a, step_a, 0, 100, 90, 0);
        add_candidate(path_b, step_b, 0, 100, 95, 0);
        add_candidate(path_b, step_b, 0, 100, 95, 200);
        add_candidate(path_b, step_b, 0, 100, 95, 400);
        expected_primary = path_b;
    }
    SECTION("Disjoint supplementary pieces still contribute") {
        add_candidate(path_a, step_a, 0, 100, 105, 0);
        add_candidate(path_b, step_b, 0, 60, 65, 0);
        add_candidate(path_b, step_b, 60, 100, 45, 200);
        expected_primary = path_b;
        expect_supplementary = true;
    }

    for (auto& strand : candidates) {
        surjector.choose_primary(strand.second);
    }
    const auto& pieces = candidates.at(make_pair(path_b, false));
    REQUIRE(pieces.size() == (expect_supplementary ? 2 : 3));
    CHECK_FALSE(pieces.front().first.is_secondary());
    if (expect_supplementary) {
        CHECK_FALSE(pieces[1].first.is_secondary());
        REQUIRE(has_annotation(pieces[1].first, "supplementary"));
        CHECK(get_annotation<bool>(pieces[1].first, "supplementary"));
    } else {
        CHECK(pieces[1].first.is_secondary());
        CHECK(pieces[2].first.is_secondary());
    }
    CHECK(surjector.choose_primary_strand(candidates)
          == make_pair(expected_primary, false));
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

    SECTION("Primary input") {
        read.set_is_secondary(false);
    }
    SECTION("Already-secondary input retains its supplementary piece") {
        read.set_is_secondary(true);
    }

    auto output = surjector.surject(read, {main_path, repeat_path});
    REQUIRE(output.size() == 3);
    size_t supplementary_count = 0, secondary_count = 0, linked_count = 0;
    for (const auto& aln : output) {
        supplementary_count += is_supplementary(aln);
        secondary_count += aln.is_secondary();
        if (has_annotation(aln, "tags")) {
            const auto tags = get_annotation<string>(aln, "tags");
            linked_count += tags.find("SA:Z:") != string::npos;
            // Each member links only to its partner, not the repeat alternative.
            CHECK(std::count(tags.begin(), tags.end(), ';') == 1);
        }
        if (is_supplementary(aln)) {
            REQUIRE(aln.refpos_size() == 1);
            CHECK(aln.refpos(0).name() == "repeat");
            CHECK(aln.is_secondary() == read.is_secondary());
        }
    }
    CHECK(supplementary_count == 1);
    CHECK(secondary_count == (read.is_secondary() ? 3 : 1));
    CHECK(linked_count == 2);
}



TEST_CASE("Anchor sliding checks read and target-path repeats",
          "[surject][anchor-sliding]") {
    string region = "ACGATTACGACCCC";
    size_t anchor_start = 0;
    size_t ref_span = 4;
    int64_t max_slide = 6;
    bool reverse_steps = false, reverse_on_path = false;
    bool split_nodes = false, sentinel_first = false;
    bool insertion = false;
    string read_between;
    bool prune = true;

    SECTION("Target-path duplicate at the positive slide limit") {
    }
    SECTION("Target-path duplicate at the negative slide limit") {
        anchor_start = 6;
    }
    SECTION("Target-path duplicate outside the slide limit") {
        max_slide = 5;
        prune = false;
    }
    SECTION("Zero slide limit disables both searches") {
        max_slide = 0;
        read_between = "TTACGA";
        prune = false;
    }
    SECTION("Unique anchor at the beginning of the path") {
        region = "ACGATTCGCCCC";
        prune = false;
    }
    SECTION("Unique anchor at the end of the path") {
        region = "CCCCTTCGACGA";
        anchor_start = region.size() - ref_span;
        sentinel_first = true;
        prune = false;
    }
    SECTION("Read duplicate without a target-path duplicate") {
        region = "ACGATTCGCCCC";
        read_between = "TTACGA";
    }
    SECTION("Reverse mapping on a reverse step is forward on the path") {
        reverse_steps = true;
    }
    SECTION("Reverse mapping on a forward step is reverse on the path") {
        reverse_on_path = true;
    }
    SECTION("Forward mapping on a reverse step is reverse on the path") {
        reverse_steps = true;
        reverse_on_path = true;
    }
    SECTION("Anchor and duplicate cross node boundaries") {
        split_nodes = true;
    }
    SECTION("Reverse-path anchor spans multiple nodes") {
        split_nodes = true;
        reverse_on_path = true;
    }
    SECTION("Reverse-path anchor spans reverse-oriented nodes") {
        split_nodes = true;
        reverse_steps = true;
        reverse_on_path = true;
    }
    SECTION("Read insertion preserves the read-based search radius") {
        region = "ACGATTCGCCCC";
        insertion = true;
        read_between = "TTACGCTAGA";
        max_slide = 12;
    }

    // A distinct longer anchor prevents the keep-one-anchor fallback from
    // hiding removal of the anchor being tested.
    const string sentinel = "GATCCGTAGTCACTGACCTAGGTC";
    const string path_sequence = sentinel_first
        ? sentinel + region : region + sentinel;
    const size_t sentinel_start = sentinel_first ? 0 : region.size();
    if (sentinel_first) {
        anchor_start += sentinel.size();
    }

    bdsg::HashGraph graph;
    auto path = graph.create_path_handle("ref");
    vector<handle_t> handles;
    vector<step_handle_t> steps;
    vector<size_t> starts;
    for (size_t start = 0; start < path_sequence.size();) {
        size_t length = split_nodes
            ? min<size_t>(3, path_sequence.size() - start)
            : path_sequence.size();
        string sequence = path_sequence.substr(start, length);
        if (reverse_steps) {
            reverse_complement_in_place(sequence);
        }
        auto handle = graph.create_handle(sequence);
        if (reverse_steps) {
            handle = graph.flip(handle);
        }
        if (!handles.empty()) {
            graph.create_edge(handles.back(), handle);
        }
        handles.push_back(handle);
        starts.push_back(start);
        steps.push_back(graph.append_step(path, handle));
        start += length;
    }
    bdsg::PositionOverlay pos_graph(&graph);
    TestSurjector surjector(&pos_graph);
    surjector.prune_suspicious_anchors = true;
    surjector.prune_tail_region_anchors = false;
    surjector.max_tail_anchor_prune = 0;
    surjector.max_low_complexity_anchor_prune = 0;
    surjector.max_low_complexity_anchor_trim = 0;
    surjector.max_anchors = 1000;
    surjector.max_slide = max_slide;

    string anchor_sequence = path_sequence.substr(anchor_start, ref_span);
    string sentinel_sequence = sentinel;
    if (reverse_on_path) {
        reverse_complement_in_place(anchor_sequence);
        reverse_complement_in_place(sentinel_sequence);
    }
    if (insertion) {
        anchor_sequence.insert(2, "GCTA");
    }
    const string sequence = anchor_sequence + read_between + sentinel_sequence;
    vector<Surjector::path_chunk_t> chunks;
    vector<pair<step_handle_t, step_handle_t>> ranges;
    auto add_anchor = [&](size_t path_start, size_t span, size_t read_start,
                          size_t read_span, bool insert_bases) {
        path_t mappings;
        vector<size_t> touched;
        for (size_t j = 0; j < handles.size(); ++j) {
            size_t begin = max(path_start, starts[j]);
            size_t end = min(path_start + span,
                             starts[j] + graph.get_length(handles[j]));
            if (begin < end) {
                touched.push_back(j);
            }
        }
        if (reverse_on_path) {
            std::reverse(touched.begin(), touched.end());
        }
        for (auto j : touched) {
            size_t begin = max(path_start, starts[j]);
            size_t end = min(path_start + span,
                             starts[j] + graph.get_length(handles[j]));
            auto* mapping = mappings.add_mapping();
            auto* position = mapping->mutable_position();
            position->set_node_id(graph.get_id(handles[j]));
            position->set_is_reverse(reverse_steps != reverse_on_path);
            position->set_offset(reverse_on_path
                ? starts[j] + graph.get_length(handles[j]) - end
                : begin - starts[j]);
            auto* edit = mapping->add_edit();
            if (insert_bases) {
                edit->set_from_length(2);
                edit->set_to_length(2);
                edit = mapping->add_edit();
                edit->set_to_length(4);
                edit->set_sequence("GCTA");
                edit = mapping->add_edit();
                edit->set_from_length(2);
                edit->set_to_length(2);
            } else {
                edit->set_from_length(end - begin);
                edit->set_to_length(end - begin);
            }
        }
        chunks.emplace_back(make_pair(sequence.begin() + read_start,
                                      sequence.begin() + read_start + read_span),
                            mappings);
        ranges.emplace_back(steps[touched.front()], steps[touched.back()]);
    };
    add_anchor(anchor_start, ref_span, 0, anchor_sequence.size(), insertion);
    add_anchor(sentinel_start, sentinel.size(),
               anchor_sequence.size() + read_between.size(), sentinel.size(), false);
    const auto original_ranges = ranges;

    surjector.prune_and_trim_anchors(sequence, chunks, ranges, 0, 0);

    REQUIRE(chunks.size() == (prune ? 1 : 2));
    REQUIRE(ranges.size() == chunks.size());
    CHECK(chunks.back().first.first
          == sequence.begin() + anchor_sequence.size() + read_between.size());
    CHECK(ranges.back() == original_ranges.back());
    if (!prune) {
        CHECK(chunks.front().first.first == sequence.begin());
        CHECK(ranges.front() == original_ranges.front());
    }
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


TEST_CASE("Surjection uses tail annotations in either read orientation",
          "[surject][tail-pruning]") {
    bdsg::HashGraph graph;

    const string tail_sequence = "ACGTCAGTGCAT";
    const string bridge_sequence = "GATCTAGC";
    const string core_sequence = "TGCAGATCGTACCTGATGCACTAGGTCAGTAC";

    auto misplaced_tail = graph.create_handle(tail_sequence);
    auto spacer = graph.create_handle(string(80, 'C'));
    auto correct_tail = graph.create_handle(tail_sequence);
    auto reference_bridge = graph.create_handle(bridge_sequence);
    auto core = graph.create_handle(core_sequence);
    auto alternate_bridge = graph.create_handle(bridge_sequence);

    vector<handle_t> reference_nodes{
        misplaced_tail, spacer, correct_tail, reference_bridge, core
    };
    auto reference = graph.create_path_handle("ref");
    for (size_t i = 0; i < reference_nodes.size(); ++i) {
        graph.append_step(reference, reference_nodes[i]);
        if (i != 0) {
            graph.create_edge(reference_nodes[i - 1], reference_nodes[i]);
        }
    }
    graph.create_edge(misplaced_tail, alternate_bridge);
    graph.create_edge(alternate_bridge, core);

    bdsg::PositionOverlay pos_graph(&graph);
    Surjector surjector(&pos_graph);
    surjector.prune_suspicious_anchors = false;
    unordered_set<path_handle_t> paths{reference};

    Alignment forward_read;
    forward_read.set_name("tail-pruning");
    string sequence;
    for (auto node : vector<handle_t>{
             misplaced_tail, alternate_bridge, core}) {
        auto* mapping = forward_read.mutable_path()->add_mapping();
        mapping->set_rank(forward_read.path().mapping_size());
        mapping->mutable_position()->set_node_id(graph.get_id(node));
        auto* edit = mapping->add_edit();
        edit->set_from_length(graph.get_length(node));
        edit->set_to_length(graph.get_length(node));
        sequence += graph.get_sequence(node);
    }
    forward_read.set_sequence(sequence);
    forward_read.set_score(
        Aligner().scorer->score_contiguous_alignment(forward_read));

    auto node_length = [&](nid_t node_id) -> int64_t {
        return graph.get_length(graph.get_handle(node_id));
    };

    bool reverse = false;
    SECTION("Forward read uses its left-tail annotation") {
        reverse = false;
    }
    SECTION("Reverse-complemented read uses its right-tail annotation") {
        reverse = true;
    }

    Alignment read = reverse
        ? reverse_complement_alignment(forward_read, node_length)
        : forward_read;

    // Keep the whole read and use ordinary, unspliced surjection.
    auto project = [&](const Alignment& input, bool prune) {
        surjector.prune_tail_region_anchors = prune;
        auto result = surjector.surject(input, paths, true, false);
        REQUIRE_FALSE(result.empty());
        for (const auto& alignment : result) {
            REQUIRE(alignment.path().mapping_size() > 0);
            REQUIRE(alignment.refpos_size() == 1);
        }
        return result;
    };

    auto check_same_alignment = [&](const Alignment& actual,
                                    const Alignment& expected) {
        CHECK(actual.sequence() == expected.sequence());
        CHECK(actual.score() == expected.score());
        CHECK(actual.path().SerializeAsString()
              == expected.path().SerializeAsString());
        CHECK(actual.refpos(0).name() == expected.refpos(0).name());
        CHECK(actual.refpos(0).offset() == expected.refpos(0).offset());
        CHECK(actual.refpos(0).is_reverse()
              == expected.refpos(0).is_reverse());
    };

    auto check_same_outputs = [&](const vector<Alignment>& actual,
                                  const vector<Alignment>& expected) {
        REQUIRE(actual.size() == expected.size());
        for (size_t i = 0; i < actual.size(); ++i) {
            check_same_alignment(actual[i], expected[i]);
        }
    };

    const auto baseline = project(read, false);
    const auto without_annotations = project(read, true);
    check_same_outputs(without_annotations, baseline);

    // Annotations describe the stored read sequence. Reversing the read
    // moves this tail from left to right. Leave the other annotation absent.
    Alignment annotated = read;
    set_annotation<double>(
        annotated,
        reverse ? "right_tail_length" : "left_tail_length",
        static_cast<double>(tail_sequence.size()));

    const auto disabled = project(annotated, false);
    check_same_outputs(disabled, baseline);
    const auto pruned_outputs = project(annotated, true);
    REQUIRE(pruned_outputs.size() == 1);
    const auto& pruned = pruned_outputs.front();

    CHECK(pruned.refpos(0).name() == "ref");
    CHECK(pruned.refpos(0).is_reverse() == reverse);
    CHECK(pruned.refpos(0).offset()
          == static_cast<int64_t>(tail_sequence.size()
                                  + graph.get_length(spacer)));
    // Detect a disconnected annotation reader: enabling pruning must
    // actually change the placement in this fixture.
    for (const auto& alignment : baseline) {
        CHECK(pruned.path().SerializeAsString()
              != alignment.path().SerializeAsString());
    }

    Alignment normalized = reverse
        ? reverse_complement_alignment(pruned, node_length)
        : pruned;
    CHECK(normalized.sequence() == forward_read.sequence());

    const vector<handle_t> expected_nodes{
        correct_tail, reference_bridge, core
    };
    REQUIRE(normalized.path().mapping_size() == expected_nodes.size());
    for (size_t i = 0; i < expected_nodes.size(); ++i) {
        const auto& mapping = normalized.path().mapping(i);
        const auto node = expected_nodes[i];
        CHECK(mapping.position().node_id() == graph.get_id(node));
        CHECK_FALSE(mapping.position().is_reverse());
        CHECK(mapping.position().offset() == 0);
        REQUIRE(mapping.edit_size() == 1);
        CHECK(mapping.edit(0).from_length() == graph.get_length(node));
        CHECK(mapping.edit(0).to_length() == graph.get_length(node));
        CHECK(mapping.edit(0).sequence().empty());
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
