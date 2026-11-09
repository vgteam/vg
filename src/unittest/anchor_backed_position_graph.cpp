///
/// \file anchor_backed_position_graph.cpp
///
/// Unit tests for AnchorBackedPositionGraph.
///
/// Most tests check positions and steps against ones worked out by hand on a
/// small path, for different sets of anchors. They do not compare against a
/// ReferencePathOverlay, because its step handles are its own and differ from
/// those of the base graph. One test checks that surjecting through this graph
/// and through a ReferencePathOverlay gives the same alignments.
///

#include "catch.hpp"
#include "../anchor_backed_position_graph.hpp"
#include "../surjector.hpp"

#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/reference_path_overlay.hpp>
#include <vg/vg.pb.h>

#include <google/protobuf/util/message_differencer.h>

#include <unordered_set>
#include <vector>

namespace vg {
namespace unittest {

namespace {

/// Build a small graph with one path, y, that visits two nodes twice.
///
/// Nodes:
///   n1: ATCG  (len 4)
///   n2: A     (len 1)
///   n3: G     (len 1)
///   n4: TTC   (len 3)
///   n5: A     (len 1)   [unused by path y]
///   n6: C     (len 1)
///   n7: TAG   (len 3)
///
/// Path y: n1, n3, n4, n1, n2, n4, n6, n7. Length 20 bases, 8 steps.
/// y visits n1 at steps 0 and 3, and n4 at steps 2 and 5.
struct TestGraph {
    bdsg::HashGraph graph;
    handle_t n1, n2, n3, n4, n5, n6, n7;
    path_handle_t y;
    std::vector<step_handle_t> y_steps;          // in path order
    std::vector<size_t>        y_positions;       // offset of the first base of each step
    size_t                     y_total_length = 0;
};

TestGraph build_test_graph() {
    TestGraph t;
    t.n1 = t.graph.create_handle("ATCG");
    t.n2 = t.graph.create_handle("A");
    t.n3 = t.graph.create_handle("G");
    t.n4 = t.graph.create_handle("TTC");
    t.n5 = t.graph.create_handle("A");
    t.n6 = t.graph.create_handle("C");
    t.n7 = t.graph.create_handle("TAG");

    t.graph.create_edge(t.n1, t.n3);
    t.graph.create_edge(t.n3, t.n4);
    t.graph.create_edge(t.n4, t.n1);
    t.graph.create_edge(t.n1, t.n2);
    t.graph.create_edge(t.n2, t.n4);
    t.graph.create_edge(t.n4, t.n6);
    t.graph.create_edge(t.n6, t.n7);

    t.y = t.graph.create_path_handle("y");
    auto append = [&](handle_t h) { t.graph.append_step(t.y, h); };
    append(t.n1); append(t.n3); append(t.n4); append(t.n1);
    append(t.n2); append(t.n4); append(t.n6); append(t.n7);

    size_t pos = 0;
    step_handle_t step = t.graph.path_begin(t.y);
    for (size_t i = 0; i < 8; ++i) {
        t.y_steps.push_back(step);
        t.y_positions.push_back(pos);
        pos += t.graph.get_length(t.graph.get_handle_of_step(step));
        if (i + 1 < 8) {
            step = t.graph.get_next_step(step);
        }
    }
    t.y_total_length = pos;  // 20
    return t;
}

/// Make an anchor whose run is the single step at pos_begin.
PrecomputedAnchor anchor_for_step(const step_handle_t& step, size_t pos_begin) {
    PrecomputedAnchor a;
    a.step_begin = step;
    a.step_end = step;
    a.path_offset_step_begin = pos_begin;
    a.path_offset_step_end   = pos_begin;
    return a;
}

/// Return the index of the step of y that covers a position.
/// Path y bases: [0..3]=n1, 4=n3, [5..7]=n4, [8..11]=n1, 12=n2, [13..15]=n4,
///               16=n6, [17..19]=n7.
size_t expected_step_index_at(size_t pos) {
    if (pos <  4) return 0;
    if (pos <  5) return 1;
    if (pos <  8) return 2;
    if (pos < 12) return 3;
    if (pos < 13) return 4;
    if (pos < 16) return 5;
    if (pos < 17) return 6;
    return 7;
}

}  // namespace

TEST_CASE("AnchorBackedPositionGraph: get_path_length returns the cached length",
          "[surject][anchors][anchor_backed_graph]") {

    TestGraph t = build_test_graph();
    std::vector<PrecomputedAnchor> anchors;
    anchors.push_back(anchor_for_step(t.y_steps.front(), 0));

    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    REQUIRE(adapter.get_path_length(t.y) == 20);
}

TEST_CASE("AnchorBackedPositionGraph: dense anchors — every query resolves via cache",
          "[surject][anchors][anchor_backed_graph]") {

    TestGraph t = build_test_graph();

    // One anchor per step on y.
    std::vector<PrecomputedAnchor> anchors;
    for (size_t i = 0; i < t.y_steps.size(); ++i) {
        anchors.push_back(anchor_for_step(t.y_steps[i], t.y_positions[i]));
    }

    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    SECTION("get_position_of_step matches hand-computed positions for every step") {
        for (size_t i = 0; i < t.y_steps.size(); ++i) {
            INFO("step index " << i);
            REQUIRE(adapter.get_position_of_step(t.y_steps[i]) == t.y_positions[i]);
        }
    }

    SECTION("get_step_at_position returns the right step for every base offset on y") {
        for (size_t pos = 0; pos < t.y_total_length; ++pos) {
            INFO("pos = " << pos);
            step_handle_t s = adapter.get_step_at_position(t.y, pos);
            REQUIRE(s == t.y_steps[expected_step_index_at(pos)]);
        }
    }

    SECTION("get_step_at_position returns the underlying handle of the right node") {
        for (size_t pos = 0; pos < t.y_total_length; ++pos) {
            INFO("pos = " << pos);
            step_handle_t s = adapter.get_step_at_position(t.y, pos);
            handle_t got = adapter.get_handle_of_step(s);
            handle_t want = t.graph.get_handle_of_step(t.y_steps[expected_step_index_at(pos)]);
            REQUIRE(got == want);
        }
    }

    SECTION("position past the end returns path_end") {
        REQUIRE(adapter.get_step_at_position(t.y, t.y_total_length) ==
                t.graph.path_end(t.y));
        REQUIRE(adapter.get_step_at_position(t.y, t.y_total_length + 100) ==
                t.graph.path_end(t.y));
    }
}

TEST_CASE("AnchorBackedPositionGraph: sparse anchors exercise the walk fallback",
          "[surject][anchors][anchor_backed_graph]") {

    TestGraph t = build_test_graph();

    // Only the first and last steps are anchors, so every other position
    // comes from walking the path.
    std::vector<PrecomputedAnchor> anchors;
    anchors.push_back(anchor_for_step(t.y_steps.front(), t.y_positions.front()));
    anchors.push_back(anchor_for_step(t.y_steps.back(),  t.y_positions.back()));

    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    SECTION("get_position_of_step gives correct values for unanchored steps") {
        for (size_t i = 0; i < t.y_steps.size(); ++i) {
            INFO("step index " << i << ", anchored=" << (i == 0 || i + 1 == t.y_steps.size()));
            REQUIRE(adapter.get_position_of_step(t.y_steps[i]) == t.y_positions[i]);
        }
    }

    SECTION("memoization: repeating a query returns the same value") {
        for (size_t i = 1; i + 1 < t.y_steps.size(); ++i) {
            size_t first  = adapter.get_position_of_step(t.y_steps[i]);
            size_t second = adapter.get_position_of_step(t.y_steps[i]);
            REQUIRE(first == second);
            REQUIRE(first == t.y_positions[i]);
        }
    }

    SECTION("get_step_at_position still works with sparse anchors") {
        for (size_t pos = 0; pos < t.y_total_length; ++pos) {
            INFO("pos = " << pos);
            step_handle_t s = adapter.get_step_at_position(t.y, pos);
            REQUIRE(s == t.y_steps[expected_step_index_at(pos)]);
        }
    }
}

TEST_CASE("AnchorBackedPositionGraph: tandem-repeat steps are disambiguated",
          "[surject][anchors][anchor_backed_graph]") {

    // y visits n1 at steps 0 (base 0) and 3 (base 8), and n4 at steps 2
    // (base 5) and 5 (base 13). Only the first visits are anchors. The
    // second visits must get their own positions.

    TestGraph t = build_test_graph();
    std::vector<PrecomputedAnchor> anchors;
    anchors.push_back(anchor_for_step(t.y_steps[0], t.y_positions[0]));   // n1 at base 0
    anchors.push_back(anchor_for_step(t.y_steps[2], t.y_positions[2]));   // n4 at base 5

    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    SECTION("second n1 visit (step 3) resolves to base 8") {
        REQUIRE(adapter.get_position_of_step(t.y_steps[3]) == 8);
    }

    SECTION("second n4 visit (step 5) resolves to base 13") {
        REQUIRE(adapter.get_position_of_step(t.y_steps[5]) == 13);
    }

    SECTION("get_step_at_position(8) returns step 3, not step 0") {
        REQUIRE(adapter.get_step_at_position(t.y, 8) == t.y_steps[3]);
    }

    SECTION("get_step_at_position(13) returns step 5, not step 2") {
        REQUIRE(adapter.get_step_at_position(t.y, 13) == t.y_steps[5]);
    }
}

TEST_CASE("AnchorBackedPositionGraph: PathHandleGraph delegation is correct",
          "[surject][anchors][anchor_backed_graph]") {

    TestGraph t = build_test_graph();
    std::vector<PrecomputedAnchor> anchors;
    for (size_t i = 0; i < t.y_steps.size(); ++i) {
        anchors.push_back(anchor_for_step(t.y_steps[i], t.y_positions[i]));
    }
    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    SECTION("HandleGraph delegation: get_id, get_length, get_sequence, get_node_count") {
        REQUIRE(adapter.get_id(t.n1) == t.graph.get_id(t.n1));
        REQUIRE(adapter.get_length(t.n1) == 4);
        REQUIRE(adapter.get_sequence(t.n1) == "ATCG");
        REQUIRE(adapter.get_node_count() == t.graph.get_node_count());
    }

    SECTION("PathHandleGraph delegation: path_begin, get_next_step, get_handle_of_step") {
        step_handle_t s = adapter.path_begin(t.y);
        REQUIRE(s == t.graph.path_begin(t.y));
        REQUIRE(adapter.get_handle_of_step(s) == t.n1);

        step_handle_t s2 = adapter.get_next_step(s);
        REQUIRE(s2 == t.graph.get_next_step(t.graph.path_begin(t.y)));
        REQUIRE(adapter.get_handle_of_step(s2) == t.n3);
    }

    SECTION("Path metadata: name, step count, circular flag") {
        REQUIRE(adapter.get_path_name(t.y) == "y");
        REQUIRE(adapter.get_step_count(t.y) == 8);
        REQUIRE_FALSE(adapter.get_is_circular(t.y));
    }
}

TEST_CASE("AnchorBackedPositionGraph: out-of-scope queries throw",
          "[surject][anchors][anchor_backed_graph]") {

    TestGraph t = build_test_graph();

    // A second path, which the graph knows no positions on.
    path_handle_t x = t.graph.create_path_handle("x");
    t.graph.append_step(x, t.n1);
    t.graph.append_step(x, t.n2);

    std::vector<PrecomputedAnchor> anchors;
    anchors.push_back(anchor_for_step(t.y_steps.front(), 0));
    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    SECTION("get_path_length on a non-target path throws") {
        REQUIRE_THROWS(adapter.get_path_length(x));
    }
    SECTION("get_step_at_position on a non-target path throws") {
        REQUIRE_THROWS(adapter.get_step_at_position(x, 0));
    }
}

// ---------------------------------------------------------------------------
// Surject the same alignment through a ReferencePathOverlay and through an
// AnchorBackedPositionGraph, and compare the results.
// ---------------------------------------------------------------------------

namespace {

/// Append a mapping that matches len bases of h, from its start, to a path.
void append_match_mapping(::vg::Path* path, const bdsg::HashGraph& g,
                          handle_t h, size_t len) {
    ::vg::Mapping* m = path->add_mapping();
    ::vg::Position* pos = m->mutable_position();
    pos->set_node_id(g.get_id(h));
    pos->set_is_reverse(g.get_is_reverse(h));
    pos->set_offset(0);
    ::vg::Edit* e = m->add_edit();
    e->set_from_length(len);
    e->set_to_length(len);
}

/// Build an alignment of ATCGGTTC to n1, n3 and n4, the first three steps of y.
::vg::Alignment make_test_alignment(const TestGraph& t) {
    ::vg::Alignment aln;
    aln.set_name("read1");
    aln.set_sequence("ATCGGTTC");
    aln.set_quality(std::string(8, 'I'));
    ::vg::Path* p = aln.mutable_path();
    append_match_mapping(p, t.graph, t.n1, 4);
    append_match_mapping(p, t.graph, t.n3, 1);
    append_match_mapping(p, t.graph, t.n4, 3);
    return aln;
}

}  // namespace

TEST_CASE("Surject(adapter) == Surject(overlay) on the same alignment",
          "[surject][anchors][anchor_backed_graph][differential]") {

    TestGraph t = build_test_graph();

    // The expected results come from a position index over y.
    bdsg::ReferencePathOverlay overlay(&t.graph,
                                      /*all_paths=*/false,
                                      std::unordered_set<std::string>{ "y" });

    // Every step of y is an anchor.
    std::vector<PrecomputedAnchor> anchors;
    for (size_t i = 0; i < t.y_steps.size(); ++i) {
        anchors.push_back(anchor_for_step(t.y_steps[i], t.y_positions[i]));
    }
    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    ::vg::Alignment source = make_test_alignment(t);

    Surjector surjector_overlay(&overlay);
    Surjector surjector_adapter(&adapter);

    std::vector<std::tuple<std::string, int64_t, bool>> positions_overlay;
    std::vector<std::tuple<std::string, int64_t, bool>> positions_adapter;

    auto results_overlay = surjector_overlay.surject(
        source, std::unordered_set<path_handle_t>{ t.y },
        positions_overlay, /*allow_negative_scores=*/false,
        /*preserve_deletions=*/false);

    auto results_adapter = surjector_adapter.surject(
        source, std::unordered_set<path_handle_t>{ t.y },
        positions_adapter, /*allow_negative_scores=*/false,
        /*preserve_deletions=*/false);

    SECTION("Both paths return the same number of surjected alignments") {
        REQUIRE(results_overlay.size() == results_adapter.size());
        REQUIRE(positions_overlay.size() == positions_adapter.size());
    }

    SECTION("Both paths return the same positions_out tuples") {
        REQUIRE(positions_overlay == positions_adapter);
    }

    SECTION("Surjected alignments have identical scores") {
        REQUIRE(results_overlay.size() == results_adapter.size());
        for (size_t i = 0; i < results_overlay.size(); ++i) {
            INFO("alignment " << i);
            REQUIRE(results_overlay[i].score() == results_adapter[i].score());
        }
    }

    SECTION("Surjected alignments have identical sequences") {
        REQUIRE(results_overlay.size() == results_adapter.size());
        for (size_t i = 0; i < results_overlay.size(); ++i) {
            INFO("alignment " << i);
            REQUIRE(results_overlay[i].sequence() == results_adapter[i].sequence());
        }
    }

    SECTION("Surjected Path proto messages are byte-equal") {
        REQUIRE(results_overlay.size() == results_adapter.size());
        for (size_t i = 0; i < results_overlay.size(); ++i) {
            INFO("alignment " << i << " — overlay path: " << results_overlay[i].path().DebugString()
                 << "\nadapter path: " << results_adapter[i].path().DebugString());
            REQUIRE(google::protobuf::util::MessageDifferencer::Equals(
                results_overlay[i].path(), results_adapter[i].path()));
        }
    }
}

// ---------------------------------------------------------------------------
// Tests of how much work a query does
// ---------------------------------------------------------------------------

TEST_CASE("AnchorBackedPositionGraph: get_position_of_step walks O(gap) to the "
          "nearest known step, not O(path) from the start [regression]",
          "[surject][anchors][anchor_backed_graph][regression]") {

    // With an anchor at step 5, a query for step 6 walks one step back to it.
    TestGraph t = build_test_graph();

    std::vector<PrecomputedAnchor> anchors;
    anchors.push_back(anchor_for_step(t.y_steps[0], t.y_positions[0]));   // base 0
    anchors.push_back(anchor_for_step(t.y_steps[5], t.y_positions[5]));   // base 13, adjacent to step 6
    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    SECTION("interior query walks backward only to the nearest known step") {
        g_anchor_graph_stats = {};
        REQUIRE(adapter.get_position_of_step(t.y_steps[6]) == t.y_positions[6]);
        REQUIRE(g_anchor_graph_stats.pos_of_step_calls == 1);
        REQUIRE(g_anchor_graph_stats.pos_of_step_walk == 1);   // 6 -> 5
    }

    SECTION("a repeated query is served from the memoized backfill (0 walk)") {
        adapter.get_position_of_step(t.y_steps[6]);            // memoizes step 6
        g_anchor_graph_stats = {};
        REQUIRE(adapter.get_position_of_step(t.y_steps[6]) == t.y_positions[6]);
        REQUIRE(g_anchor_graph_stats.pos_of_step_walk == 0);
    }
}

TEST_CASE("AnchorBackedPositionGraph: for_each_step_on_handle returns only the "
          "target path's steps, not every haplotype's [regression]",
          "[surject][anchors][anchor_backed_graph][regression]") {

    // y visits n1 twice and a second path, x, visits it once. Only y's two
    // visits are reported.
    TestGraph t = build_test_graph();

    path_handle_t x = t.graph.create_path_handle("x");
    t.graph.append_step(x, t.n2);
    t.graph.append_step(x, t.n1);   // x also visits n1
    t.graph.append_step(x, t.n4);   // ...and n4

    // One anchor whose run is all of y.
    PrecomputedAnchor whole;
    whole.step_begin = t.y_steps.front();
    whole.step_end   = t.y_steps.back();
    whole.path_offset_step_begin = t.y_positions.front();
    whole.path_offset_step_end   = t.y_positions.back();
    AnchorBackedPositionGraph adapter(&t.graph, {whole}, t.y, t.y_total_length);

    SECTION("base graph sees all paths on n1; adapter sees only the target's") {
        size_t base_count = 0;
        t.graph.for_each_step_on_handle(t.n1, [&](const step_handle_t&) {
            base_count++; return true; });
        REQUIRE(base_count == 3);   // y@step0, y@step3, x

        std::vector<step_handle_t> got;
        adapter.for_each_step_on_handle(t.n1, [&](const step_handle_t& s) {
            got.push_back(s); return true; });
        REQUIRE(got.size() == 2);          // only y's two visits
        REQUIRE(got.size() < base_count);
    }

    SECTION("returned steps are all on the target path; a repeat node yields both visits") {
        std::vector<step_handle_t> got;
        adapter.for_each_step_on_handle(t.n1, [&](const step_handle_t& s) {
            got.push_back(s); return true; });
        for (const step_handle_t& s : got) {
            REQUIRE(adapter.get_path_handle_of_step(s) == t.y);
        }
        bool has0 = false, has3 = false;
        for (const step_handle_t& s : got) {
            if (s == t.y_steps[0]) has0 = true;
            if (s == t.y_steps[3]) has3 = true;
        }
        REQUIRE(has0);
        REQUIRE(has3);
    }

    SECTION("a node not in the target region yields nothing") {
        std::vector<step_handle_t> got;
        adapter.for_each_step_on_handle(t.n5, [&](const step_handle_t& s) {
            got.push_back(s); return true; });   // n5 is on no path here
        REQUIRE(got.empty());
    }

    SECTION("a step in two anchors' runs is reported once") {
        // The two runs share step 2, on n4. y visits n4 at steps 2 and 5.
        PrecomputedAnchor left, right;
        left.step_begin = t.y_steps[0];
        left.step_end   = t.y_steps[2];
        left.path_offset_step_begin = t.y_positions[0];
        left.path_offset_step_end   = t.y_positions[2];
        right.step_begin = t.y_steps[2];
        right.step_end   = t.y_steps.back();
        right.path_offset_step_begin = t.y_positions[2];
        right.path_offset_step_end   = t.y_positions.back();
        AnchorBackedPositionGraph split(&t.graph, {left, right}, t.y, t.y_total_length);

        std::vector<step_handle_t> got;
        split.for_each_step_on_handle(t.n4, [&](const step_handle_t& s) {
            got.push_back(s); return true; });
        REQUIRE(got.size() == 2);
    }
}

} // namespace unittest
} // namespace vg
