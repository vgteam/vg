/// \file unittest/fill_pins.cpp
///
/// Unit tests for fill anchors: where the fill points go in runs of nodes with no branch, which
/// reads cross them, and which slot each fill anchor gets from the called walks.

#include <cmath>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "../fill_pins.hpp"
#include "../edit_snarl_decomposition_overlay.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "../../src/unittest/catch.hpp"
#include "../../src/utility.hpp"
#include "support/decomposition_pair.hpp"
#include "support/nested_graph.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// Return `length` bases from a fixed pseudo-random sequence, different for each `seed`.
static string bases(size_t length, uint32_t seed) {
    string out;
    uint32_t state = seed * 2654435761U + 1;
    for (size_t i = 0; i < length; ++i) {
        state = state * 1664525U + 1013904223U;
        out.push_back("ACGT"[(state >> 16) & 3]);
    }
    return out;
}

/// Add a node of `length` bases with ID `id`.
static handle_t add_node(bdsg::HashGraph& graph, nid_t id, size_t length) {
    return graph.create_handle(bases(length, (uint32_t)id), id);
}

/// Return the fill points of `graph` at `spacing`, and set `tally`.
static vector<FillPoint> points_of(const HandleGraph& graph, size_t spacing,
                                   FillPointTally& tally) {
    return place_fill_points(graph, spacing, tally);
}

TEST_CASE("Fill points divide a run of nodes evenly", "[fill]") {
    bdsg::HashGraph graph;
    FillPointTally tally;
    SECTION("three nodes of 600 bp at a spacing of 1000 get one point, in the middle node") {
        add_node(graph, 1, 600);
        add_node(graph, 2, 600);
        add_node(graph, 3, 600);
        graph.create_edge(graph.get_handle(1), graph.get_handle(2));
        graph.create_edge(graph.get_handle(2), graph.get_handle(3));
        REQUIRE(points_of(graph, 1000, tally) == vector<FillPoint>{{2, 300}});
        REQUIRE(tally.unitigs == 1);
        REQUIRE(tally.placed == 1);
    }
    SECTION("a node of 2500 bp gets two, a third of the way along and two thirds") {
        add_node(graph, 1, 2500);
        REQUIRE(points_of(graph, 1000, tally) == vector<FillPoint>{{1, 833}, {1, 1667}});
    }
}

TEST_CASE("A run of nodes ends where the graph branches", "[fill]") {
    // A bubble 1 -> (2 or 3) -> 4, then 4 -> 5. Every node is 600 bp, so only the run 4, 5 is
    // longer than the spacing. Its point falls on the join between them and moves into node 5.
    bdsg::HashGraph graph;
    for (nid_t id = 1; id <= 5; ++id) {
        add_node(graph, id, 600);
    }
    for (auto [from, to] : vector<pair<nid_t, nid_t>>{{1, 2}, {1, 3}, {2, 4}, {3, 4}, {4, 5}}) {
        graph.create_edge(graph.get_handle(from), graph.get_handle(to));
    }
    FillPointTally tally;
    REQUIRE(points_of(graph, 1000, tally) == vector<FillPoint>{{5, 1}});
    REQUIRE(tally.unitigs == 1);
}

TEST_CASE("A run no longer than the spacing gets no point", "[fill]") {
    bdsg::HashGraph graph;
    add_node(graph, 1, 1000);
    FillPointTally tally;
    REQUIRE(points_of(graph, 1000, tally).empty());
    REQUIRE(tally.unitigs == 0);
    REQUIRE(points_of(graph, 999, tally) == vector<FillPoint>{{1, 500}});
    REQUIRE(points_of(graph, 0, tally).empty());
}

TEST_CASE("A point on a join between nodes moves forward into a node", "[fill]") {
    bdsg::HashGraph graph;
    FillPointTally tally;
    SECTION("past a node of one base, which can hold no point") {
        add_node(graph, 10, 500);
        add_node(graph, 11, 1);
        add_node(graph, 12, 499);
        graph.create_edge(graph.get_handle(10), graph.get_handle(11));
        graph.create_edge(graph.get_handle(11), graph.get_handle(12));
        REQUIRE(points_of(graph, 500, tally) == vector<FillPoint>{{12, 1}});
        REQUIRE(tally.dropped_at_end == 0);
    }
    SECTION("and is dropped where the move reaches the end of the run") {
        add_node(graph, 20, 2);
        add_node(graph, 21, 1);
        graph.create_edge(graph.get_handle(20), graph.get_handle(21));
        REQUIRE(points_of(graph, 2, tally).empty());
        REQUIRE(tally.dropped_at_end == 1);
    }
}

TEST_CASE("A run gets the same points whichever end it is read from", "[fill]") {
    // A run of two 500 bp nodes has its one point on the join, which moves forward into the
    // second node read, so reading the run from the other end would give another point.
    bdsg::HashGraph graph;
    const handle_t a = add_node(graph, 41, 500);
    const handle_t b = add_node(graph, 42, 500);
    graph.create_edge(a, b);
    const vector<handle_t> forward{a, b};
    const vector<handle_t> backward{graph.flip(b), graph.flip(a)};

    FillPointTally tally;
    vector<FillPoint> raw_forward;
    vector<FillPoint> raw_backward;
    place_unitig_points(graph, forward, 500, raw_forward, tally);
    place_unitig_points(graph, backward, 500, raw_backward, tally);
    REQUIRE(raw_forward == vector<FillPoint>{{42, 1}});
    REQUIRE(raw_backward == vector<FillPoint>{{41, 499}});

    vector<FillPoint> from_forward;
    vector<FillPoint> from_backward;
    place_unitig_points(graph, canonical_unitig(graph, forward), 500, from_forward, tally);
    place_unitig_points(graph, canonical_unitig(graph, backward), 500, from_backward, tally);
    REQUIRE(from_forward == raw_forward);
    REQUIRE(from_backward == raw_forward);
    REQUIRE(points_of(graph, 500, tally) == raw_forward);
}

TEST_CASE("A point that moves off a join past the next point takes its place", "[fill]") {
    // Nodes of 2, 1, 1 and 3 bp make a run of 7 bp, which a spacing of 2 divides at 2, 4 and 5
    // bases. The first point moves off the joins at 2, 3 and 4 bases to 5, passing the second and
    // reaching the third, which are dropped.
    bdsg::HashGraph graph;
    const vector<size_t> lengths{2, 1, 1, 3};
    for (size_t i = 0; i < lengths.size(); ++i) {
        add_node(graph, (nid_t)(i + 1), lengths[i]);
        if (i > 0) {
            graph.create_edge(graph.get_handle((nid_t)i), graph.get_handle((nid_t)(i + 1)));
        }
    }
    FillPointTally tally;
    REQUIRE(points_of(graph, 2, tally) == vector<FillPoint>{{4, 1}});
    REQUIRE(tally.merged == 2);
}

TEST_CASE("A run that closes on itself is read from its node with the smallest ID", "[fill]") {
    // Nodes 7 and 8 joined into a cycle. Its point falls on the join after node 7 and moves one
    // base into node 8.
    bdsg::HashGraph graph;
    const handle_t a = add_node(graph, 7, 600);
    const handle_t b = add_node(graph, 8, 600);
    graph.create_edge(a, b);
    graph.create_edge(b, a);
    FillPointTally tally;
    REQUIRE(points_of(graph, 1000, tally) == vector<FillPoint>{{8, 1}});
    REQUIRE(tally.unitigs == 1);
}

TEST_CASE("A point in a node the run reads in reverse is pinned along the node's forward strand",
          "[fill]") {
    // The run is 51 forward then 52 reverse, joined right side to right side. Its point falls on
    // the join and moves one base into node 52, which the run reads backward, so one base of node
    // 52's forward strand lies after the point and 999 before it.
    bdsg::HashGraph graph;
    const handle_t a = add_node(graph, 51, 1000);
    const handle_t b = add_node(graph, 52, 1000);
    graph.create_edge(a, graph.flip(b));
    FillPointTally tally;
    REQUIRE(points_of(graph, 1500, tally) == vector<FillPoint>{{52, 999}});
}

TEST_CASE("A called genotype marks the fill points of its walks by strand", "[fill]") {
    using Marks = vector<pair<int, uint8_t>>;
    REQUIRE(fill_marks_of_call({3, 5}, -1, false) == Marks{{3, FILL_STRAND_0}, {5, FILL_STRAND_1}});
    REQUIRE(fill_marks_of_call({2, 2}, -1, true) == Marks{{2, FILL_EVERY_STRAND}});
    // One allele at a site that is not nested is in a region of ploidy 1.
    REQUIRE(fill_marks_of_call({4}, -1, false) == Marks{{4, FILL_EVERY_STRAND}});
    REQUIRE(fill_marks_of_call({4}, 1, true) == Marks{{4, FILL_STRAND_1}});
    // One allele at a nested site whose strand is not known marks nothing.
    REQUIRE(fill_marks_of_call({4}, -1, true).empty());
    REQUIRE(fill_marks_of_call({1, 2, 1}, -1, false)
            == Marks{{1, FILL_EVERY_STRAND}, {2, FILL_EVERY_STRAND}});
}

TEST_CASE("A fill anchor's kind and slot come from the marks on its point", "[fill]") {
    auto is = [](FillSlot got, AnchorKind kind, int slot) {
        return got.kind == kind && got.slot == slot;
    };
    REQUIRE(is(classify_fill_point(true, 0), AnchorKind::fill_shared, 0));
    REQUIRE(is(classify_fill_point(false, FILL_STRAND_0 | FILL_STRAND_1), AnchorKind::fill_shared, 0));
    REQUIRE(is(classify_fill_point(false, FILL_EVERY_STRAND), AnchorKind::fill_shared, 0));
    REQUIRE(is(classify_fill_point(false, FILL_STRAND_0), AnchorKind::fill_strand, 0));
    REQUIRE(is(classify_fill_point(false, FILL_STRAND_1), AnchorKind::fill_strand, 1));
    REQUIRE(is(classify_fill_point(false, 0), AnchorKind::fill_off_call, 0));
}

TEST_CASE("A read's fill anchor score is its clamped mismapping probability", "[fill]") {
    const FillScoreParams params{true, 0.02, 0.95};
    REQUIRE(std::abs(fill_read_score(60, params) - (float)(-10.0 * log10(0.02))) < 1e-4);
    REQUIRE(std::abs(fill_read_score(0, params) - (float)(-10.0 * log10(0.95))) < 1e-4);
    REQUIRE(std::abs(fill_read_score(10, params) - 10.0f) < 1e-4);
    const FillScoreParams flat{false, 0.02, 0.95};
    REQUIRE(fill_read_score(0, flat) == fill_read_score(60, params));
}

/// A bubble with long nodes: 1 -> (2 or 3) -> 4, each node 1200 bp, so that each node holds one
/// fill point at a spacing of 1000. Nodes 1 and 4 lie in no site, and 2 and 3 in the site 1..4.
struct FillBubble {
    bdsg::HashGraph graph;
    FillBubble() {
        for (nid_t id = 1; id <= 4; ++id) {
            add_node(graph, id, 1200);
        }
        for (auto [from, to] : vector<pair<nid_t, nid_t>>{{1, 2}, {1, 3}, {2, 4}, {3, 4}}) {
            graph.create_edge(graph.get_handle(from), graph.get_handle(to));
        }
    }
    /// The walk visiting `middle` between nodes 1 and 4.
    Traversal walk(nid_t middle) const {
        return Traversal{graph.get_handle(1), graph.get_handle(middle), graph.get_handle(4)};
    }
};

/// An alignment of `name` matching `graph` along `walk` exactly.
static Alignment along(const HandleGraph& graph, const string& name, const Traversal& walk) {
    Alignment aln;
    aln.set_name(name);
    string sequence;
    for (const handle_t& h : walk) {
        const string node = graph.get_sequence(h);
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(graph.get_id(h));
        m->mutable_position()->set_is_reverse(graph.get_is_reverse(h));
        Edit* e = m->add_edit();
        e->set_from_length((int32_t)node.size());
        e->set_to_length((int32_t)node.size());
        sequence += node;
    }
    aln.set_sequence(sequence);
    aln.set_mapping_quality(60);
    return aln;
}

/// Return the index of the point on node `id` among the points of `fill`.
static size_t point_on(const FillAnchors& fill, nid_t id) {
    const auto& points = fill.points();
    const auto found = std::ranges::find_if(points, [&](const FillPoint& p) { return p.node == id; });
    REQUIRE(found != points.end());
    return (size_t)(found - points.begin());
}

TEST_CASE("A fill anchor's slot comes from the called walks of the site that holds its point",
          "[fill]") {
    FillBubble bubble;
    const HandleGraph& graph = bubble.graph;
    IntegratedSnarlFinder finder(graph);
    DecompositionPair pair(graph, finder);
    EditSnarlDecompositionOverlay sites(pair.adapter, graph);
    FillAnchors fill;
    fill.configure(1000, FillScoreParams{true, 0.02, 0.95}, 0.0);
    fill.place_points(graph, FILL_READ_WINDOW);
    fill.place(sites);
    REQUIRE(fill.points() == vector<FillPoint>{{1, 600}, {2, 600}, {3, 600}, {4, 600}});
    const net_handle_t site = find_snarl(sites, graph, 1, 4);
    auto mark = [&](const vector<int>& genotype, int nested_strand) {
        for (const auto& [allele, bits] :
             fill_marks_of_call(genotype, nested_strand, nested_strand >= 0)) {
            fill.mark_walk(site, bubble.walk(allele == 0 ? 2 : 3), bits);
        }
    };
    auto slot_on = [&](nid_t id) { return fill.slot_of(point_on(fill, id)); };
    auto is = [](FillSlot got, AnchorKind kind, int slot) {
        return got.kind == kind && got.slot == slot;
    };

    SECTION("a heterozygote puts each arm on its own strand, and the nodes outside every site on both") {
        mark({0, 1}, -1);
        REQUIRE(is(slot_on(1), AnchorKind::fill_shared, 0));
        REQUIRE(is(slot_on(2), AnchorKind::fill_strand, 0));
        REQUIRE(is(slot_on(3), AnchorKind::fill_strand, 1));
        REQUIRE(is(slot_on(4), AnchorKind::fill_shared, 0));
    }
    SECTION("the phase order decides which arm is strand 0") {
        mark({1, 0}, -1);
        REQUIRE(is(slot_on(2), AnchorKind::fill_strand, 1));
        REQUIRE(is(slot_on(3), AnchorKind::fill_strand, 0));
    }
    SECTION("a homozygote shares its arm, and leaves the other arm off the call") {
        mark({0, 0}, -1);
        REQUIRE(is(slot_on(2), AnchorKind::fill_shared, 0));
        REQUIRE(is(slot_on(3), AnchorKind::fill_off_call, 0));
    }
    SECTION("a nested haploid site names the strand that carries it") {
        mark({1}, 1);
        REQUIRE(is(slot_on(2), AnchorKind::fill_off_call, 0));
        REQUIRE(is(slot_on(3), AnchorKind::fill_strand, 1));
    }
    SECTION("a site in a region of ploidy 1 shares its arm") {
        mark({0}, -1);
        REQUIRE(is(slot_on(2), AnchorKind::fill_shared, 0));
    }
    SECTION("a walk marks no point that another site holds") {
        fill.mark_walk(sites.get_root(), bubble.walk(2), FILL_STRAND_0);
        REQUIRE(is(slot_on(2), AnchorKind::fill_off_call, 0));
    }
}

TEST_CASE("The innermost genotyped site holding a fill point decides its slot", "[fill]") {
    // A parent site 2..8 with arms 3 and 4..7, where 4..7 is a child site with arms 5 and 6, and
    // flanking nodes 1 and 9. Every node is 1200 bp, so each holds a point at a spacing of 1000.
    bdsg::HashGraph graph;
    for (nid_t id = 1; id <= 9; ++id) {
        add_node(graph, id, 1200);
    }
    for (auto [from, to] : vector<pair<nid_t, nid_t>>{{1, 2}, {2, 3}, {3, 8}, {2, 4}, {4, 5},
                                                      {4, 6}, {5, 7}, {6, 7}, {7, 8}, {8, 9}}) {
        graph.create_edge(graph.get_handle(from), graph.get_handle(to));
    }
    auto walk = [&](const vector<nid_t>& ids) {
        Traversal t;
        for (nid_t id : ids) {
            t.push_back(graph.get_handle(id));
        }
        return t;
    };
    IntegratedSnarlFinder finder(graph);
    DecompositionPair decomposition(graph, finder);
    EditSnarlDecompositionOverlay sites(decomposition.adapter, graph);
    FillAnchors fill;
    fill.configure(1000, FillScoreParams{true, 0.02, 0.95}, 0.0);
    fill.place_points(graph, FILL_READ_WINDOW);
    fill.place(sites);
    const net_handle_t parent = find_snarl(sites, graph, 2, 8);
    const net_handle_t child = find_snarl(sites, graph, 4, 7);
    // The parent is heterozygous, with arm 3 on strand 0 and the child's arm 5 on strand 1.
    auto mark_parent = [&]() {
        fill.mark_walk(parent, walk({2, 3, 8}), FILL_STRAND_0);
        fill.mark_walk(parent, walk({2, 4, 5, 7, 8}), FILL_STRAND_1);
    };
    auto mark_child = [&](const vector<int>& genotype, int nested_strand) {
        for (const auto& [allele, bits] : fill_marks_of_call(genotype, nested_strand, true)) {
            fill.mark_walk(child, walk({4, allele == 0 ? 5 : 6, 7}), bits);
        }
    };
    auto slot_on = [&](nid_t id) { return fill.slot_of(point_on(fill, id)); };
    auto is = [](FillSlot got, AnchorKind kind, int slot) {
        return got.kind == kind && got.slot == slot;
    };

    SECTION("the parent decides the points of a child that was not genotyped") {
        mark_parent();
        REQUIRE(is(slot_on(3), AnchorKind::fill_strand, 0));
        REQUIRE(is(slot_on(4), AnchorKind::fill_strand, 1));
        REQUIRE(is(slot_on(5), AnchorKind::fill_strand, 1));
        REQUIRE(is(slot_on(6), AnchorKind::fill_off_call, 0));
        REQUIRE(is(slot_on(9), AnchorKind::fill_shared, 0));
    }
    SECTION("a genotyped child decides its own points, whatever its parent's walks") {
        // The child, marked first, is homozygous for arm 6, which the parent's walks do not take.
        mark_child({1, 1}, -1);
        mark_parent();
        REQUIRE(is(slot_on(5), AnchorKind::fill_off_call, 0));
        REQUIRE(is(slot_on(6), AnchorKind::fill_shared, 0));
        REQUIRE(is(slot_on(4), AnchorKind::fill_strand, 1));
    }
    SECTION("a nested chain of one copy names the strand that carries it") {
        mark_parent();
        mark_child({0}, 1);
        REQUIRE(is(slot_on(5), AnchorKind::fill_strand, 1));
        REQUIRE(is(slot_on(6), AnchorKind::fill_off_call, 0));
    }
    SECTION("a nested chain of one copy whose strand is not known leaves its points to the parent") {
        mark_child({0}, -1);
        mark_parent();
        REQUIRE(is(slot_on(5), AnchorKind::fill_strand, 1));
    }
}

TEST_CASE("The reads crossing each fill point become a fill anchor's rows", "[fill]") {
    FillBubble bubble;
    const HandleGraph& graph = bubble.graph;
    IntegratedSnarlFinder finder(graph);
    DecompositionPair pair(graph, finder);
    EditSnarlDecompositionOverlay sites(pair.adapter, graph);
    FillAnchors fill;
    fill.configure(1000, FillScoreParams{true, 0.02, 0.95}, 0.0);
    fill.place_points(graph, FILL_READ_WINDOW);
    fill.place(sites);

    // One read along arm 2 forward, and one along arm 3 read backward.
    InMemorySiteReadSource reads;
    reads.add(along(graph, "forward", bubble.walk(2)));
    const Traversal arm_3 = bubble.walk(3);
    Traversal backward;
    for (auto it = arm_3.rbegin(); it != arm_3.rend(); ++it) {
        backward.push_back(graph.flip(*it));
    }
    reads.add(along(graph, "backward", backward));
    fill.resolve_reads(reads);
    REQUIRE(fill.placements().size() == 6);

    const net_handle_t site = find_snarl(sites, graph, 1, 4);
    fill.mark_walk(site, bubble.walk(2), FILL_STRAND_0);
    fill.mark_walk(site, bubble.walk(3), FILL_STRAND_1);

    const string path = temp_file::create();
    AnchorParams params;
    params.min_reads = 1;
    params.fill_spacing = 1000;
    AnchorCollector anchors;
    anchors.configure(path, params, "g", "r", 0.02);
    fill.add_anchors(anchors);
    anchors.write("s");

    // Each A row as its node, snarl, slot and kind, followed by its R rows as (read, strand,
    // offset).
    vector<string> rows;
    map<string, string> read_of;
    ifstream in(path);
    string line;
    while (getline(in, line)) {
        vector<string> fields;
        istringstream columns(line);
        string field;
        while (getline(columns, field, '\t')) {
            fields.push_back(field);
        }
        if (fields[0] == "#read") {
            read_of[fields[1]] = fields[2];
        } else if (fields[0] == "A") {
            rows.push_back(fields[1] + " " + fields[2] + " " + fields[3] + " " + fields[8] + " "
                           + fields[9]);
        } else if (fields[0] == "R") {
            rows.push_back(read_of.at(fields[1]) + " " + fields[2] + " " + fields[3]);
        }
    }
    temp_file::remove(path);
    // The forward read crosses each point after the point's node base 599, which is its own base
    // 599 on node 1, 1799 on node 2 and 2999 on node 4. The backward read visits node 4, then 3,
    // then 1, each backward, so each point lies before that node's base 599, at the read's base
    // 600 of the node.
    REQUIRE(rows == vector<string>{
        "1 fill1.600 0 600 fill-shared", "backward 1 3000", "forward 0 599",
        "2 fill2.600 0 600 fill-strand", "forward 0 1799",
        "3 fill3.600 1 600 fill-strand", "backward 1 1800",
        "4 fill4.600 0 600 fill-shared", "backward 1 600", "forward 0 2999"});
}

TEST_CASE("Reads placed from the windows a source fetches match those resolve_reads fetches",
          "[fill]") {
    FillBubble bubble;
    const HandleGraph& graph = bubble.graph;
    IntegratedSnarlFinder finder(graph);
    DecompositionPair decomposition(graph, finder);
    EditSnarlDecompositionOverlay sites(decomposition.adapter, graph);
    InMemorySiteReadSource reads;
    reads.add(along(graph, "forward", bubble.walk(2)));
    reads.add(along(graph, "other", bubble.walk(3)));
    const SiteReadSource::WindowReads window_reads =
        [&reads](const vector<pair<nid_t, nid_t>>& ranges,
                 const function<void(const SiteRead&)>& iteratee) {
            reads.for_each_read(ranges, iteratee);
        };
    auto placed = [](const FillAnchors& fill) {
        vector<tuple<uint32_t, string, uint32_t, uint32_t, float, uint8_t>> out;
        for (const FillPlacement& p : fill.placements()) {
            out.emplace_back(p.point, read_names().name(p.read), p.alignment, p.offset, p.score,
                             p.direction);
        }
        return out;
    };

    // Windows of 2 node IDs: node 1 in window 0, nodes 2 and 3 in window 1, node 4 in window 2.
    FillAnchors fetched;
    fetched.configure(1000, FillScoreParams{true, 0.02, 0.95}, 0.0);
    fetched.place_points(graph, 2);
    fetched.place(sites);
    fetched.resolve_reads(reads);

    FillAnchors observed;
    observed.configure(1000, FillScoreParams{true, 0.02, 0.95}, 0.0);
    observed.place_points(graph, 2);
    // A range that is not exactly one window is ignored, and a window's reads are placed only
    // once.
    observed.resolve_window_reads(0, 3, window_reads);
    observed.resolve_window_reads(1, 2, window_reads);
    observed.resolve_window_reads(2, 3, window_reads);
    observed.resolve_window_reads(2, 3, window_reads);
    observed.place(sites);
    observed.resolve_reads(reads);

    REQUIRE(placed(fetched).size() == 6);
    REQUIRE(placed(observed) == placed(fetched));
    ostringstream report;
    observed.report(report);
    REQUIRE(report.str().find("reads placed from 1 windows of 2 node IDs as the direct pass "
                              "fetched them, and from 2 windows fetched after it")
            != string::npos);
}

/// An alignment of `name` to node `id` of `graph` from its base `offset`, along the strand
/// `reverse` names, matching `length` bases, with `clip_before` and `clip_after` bases soft clipped
/// before and after.
static Alignment on_node(const HandleGraph& graph, const string& name, nid_t id, bool reverse,
                         size_t offset, size_t length, size_t clip_before, size_t clip_after) {
    Alignment aln;
    aln.set_name(name);
    Mapping* m = aln.mutable_path()->add_mapping();
    m->mutable_position()->set_node_id(id);
    m->mutable_position()->set_is_reverse(reverse);
    m->mutable_position()->set_offset((int64_t)offset);
    auto clip = [&](size_t bases) {
        Edit* e = m->add_edit();
        e->set_from_length(0);
        e->set_to_length((int32_t)bases);
        e->set_sequence(string(bases, 'A'));
    };
    if (clip_before > 0) {
        clip(clip_before);
    }
    Edit* e = m->add_edit();
    e->set_from_length((int32_t)length);
    e->set_to_length((int32_t)length);
    if (clip_after > 0) {
        clip(clip_after);
    }
    aln.set_sequence(string(clip_before, 'A')
                     + graph.get_sequence(graph.get_handle(id, reverse)).substr(offset, length)
                     + string(clip_after, 'A'));
    aln.set_mapping_quality(60);
    return aln;
}

TEST_CASE("A read is placed at a fill point only where it crosses it", "[fill]") {
    // One node of 1200 bp, with its point after base 599.
    bdsg::HashGraph graph;
    add_node(graph, 1, 1200);
    IntegratedSnarlFinder finder(graph);
    DecompositionPair decomposition(graph, finder);
    EditSnarlDecompositionOverlay sites(decomposition.adapter, graph);
    FillAnchors fill;
    fill.configure(1000, FillScoreParams{true, 0.02, 0.95}, 0.0);
    fill.place_points(graph, FILL_READ_WINDOW);
    fill.place(sites);
    REQUIRE(fill.points() == vector<FillPoint>{{1, 600}});

    InMemorySiteReadSource reads;
    // Two reads that do not reach the node base before the point.
    reads.add(on_node(graph, "short", 1, false, 0, 300, 0, 0));
    reads.add(on_node(graph, "late", 1, false, 700, 300, 0, 0));
    // Three reads whose aligned bases end at the point: one ends there, one goes on in a soft
    // clip, and one, along the node's reverse strand, starts there after a soft clip.
    reads.add(on_node(graph, "ends", 1, false, 400, 200, 0, 0));
    reads.add(on_node(graph, "clipped", 1, false, 400, 200, 0, 50));
    reads.add(on_node(graph, "clipped-reverse", 1, true, 600, 200, 50, 0));
    // Two reads that cross it, one each way.
    reads.add(on_node(graph, "crosses", 1, false, 400, 400, 0, 0));
    reads.add(on_node(graph, "crosses-reverse", 1, true, 400, 400, 20, 0));
    fill.resolve_reads(reads);

    vector<pair<uint32_t, uint8_t>> placed;
    for (const FillPlacement& p : fill.placements()) {
        placed.emplace_back(p.offset, p.direction);
    }
    std::ranges::sort(placed);
    // The forward read has node base 599 at its base 199. The reverse read reads node base 599
    // as its 201st aligned base, after the 20 clipped bases.
    REQUIRE(placed == vector<pair<uint32_t, uint8_t>>{{199, 0}, {220, 1}});

    ostringstream report;
    fill.report(report);
    REQUIRE(report.str().find("2 reads on the point's node end before the point or start after it")
            != string::npos);
    REQUIRE(report.str().find("3 stop at the point or go on past it only in a soft clip")
            != string::npos);
}

}
}
