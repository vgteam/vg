///
/// \file anchor_backed_position_graph.cpp
///
/// Unit tests for AnchorBackedPositionGraph.
///
/// We validate the adapter's contract directly with hand-computed truth: for
/// a tiny in-memory graph + path, we know exactly which step is at which
/// base position. The adapter must return those positions and steps
/// regardless of how anchors are arranged.
///
/// If the adapter passes these tests, Surjector running against
/// AnchorBackedPositionGraph behaves identically to Surjector against a
/// real PathPositionHandleGraph — because Surjector's only interaction with
/// the position-graph layer is through get_path_length /
/// get_position_of_step / get_step_at_position (plus delegated PathHandleGraph
/// ops, which we forward to base_).
///
/// Why not compare to bdsg::ReferencePathOverlay directly: it uses its own
/// step_handle_t encoding (`(path_handle, step_index)`, see
/// libbdsg/bdsg/src/reference_path_overlay.cpp:428-454). The adapter
/// delegates step iteration to its base PathHandleGraph (HashGraph here),
/// so step_handles aren't byte-comparable across the two implementations.
/// The hand-computed truth is simpler and just as strict.
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

/// Build a fresh small test graph + path.
///
/// Nodes (mirrors pangenome-index's coord_translation_tests/test0_target_loop):
///   n1: ATCG  (len 4)
///   n2: A     (len 1)
///   n3: G     (len 1)
///   n4: TTC   (len 3)
///   n5: A     (len 1)   [unused by path y]
///   n6: C     (len 1)
///   n7: TAG   (len 3)
///
/// Path y: n1, n3, n4, n1, n2, n4, n6, n7. Length 20 bases, 8 steps.
/// y revisits n1 (steps 0 and 3) and n4 (steps 2 and 5) — exercises the
/// tandem-repeat handling.
struct TestGraph {
    bdsg::HashGraph graph;
    handle_t n1, n2, n3, n4, n5, n6, n7;
    path_handle_t y;
    std::vector<step_handle_t> y_steps;          // in path order
    std::vector<size_t>        y_positions;       // base offset of each step's first base
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

PrecomputedAnchor anchor_for_step(const step_handle_t& step, size_t pos_begin) {
    PrecomputedAnchor a;
    a.step_begin = step;
    a.step_end = step;
    a.path_offset_step_begin = pos_begin;
    a.path_offset_step_end   = pos_begin;
    return a;
}

/// Hand-computed expected step at each base position on y.
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

    // Only the first and last steps are in the anchor set. Every other
    // step's position must be discovered via the LF-walk fallback in
    // get_position_of_step.
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

    // y visits n1 at steps 0 (base 0) and 3 (base 8); n4 at steps 2 (base 5)
    // and 5 (base 13). Make sure that with only the first occurrence of each
    // repeated node in the anchor set, position queries for the SECOND
    // occurrence still resolve correctly via the walk fallback.

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

    // Add a second path so we can ask about it.
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
// Differential test: Surjector through the overlay vs through the adapter.
// The strongest possible test — runs the full surject pipeline twice and
// asserts identical output. Passing this proves the adapter is fully
// interchangeable with bdsg::ReferencePathOverlay for surject's purposes.
// ---------------------------------------------------------------------------

namespace {

/// Append a "perfect match" Mapping covering `len` bases of node `h` to a Path.
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

/// Build a simple source Alignment that walks n1 → n3 → n4 on the graph
/// (= first 3 steps of path y, sequence ATCGGTTC).
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

    // Ground truth: real reference-path overlay over y.
    bdsg::ReferencePathOverlay overlay(&t.graph,
                                      std::unordered_set<std::string>{ "y" },
                                      /*all_paths=*/false);

    // Adapter: AnchorBackedPositionGraph with dense anchors covering y.
    std::vector<PrecomputedAnchor> anchors;
    for (size_t i = 0; i < t.y_steps.size(); ++i) {
        anchors.push_back(anchor_for_step(t.y_steps[i], t.y_positions[i]));
    }
    AnchorBackedPositionGraph adapter(&t.graph, anchors, t.y, t.y_total_length);

    // Source alignment: 8-base read walking the first 3 nodes of y.
    ::vg::Alignment source = make_test_alignment(t);

    // Surject through both.
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

} // namespace unittest
} // namespace vg
