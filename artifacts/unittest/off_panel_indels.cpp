/// \file off_panel_indels.cpp
///
/// Unit tests for off-panel indel detection: a read's indels along its walk, stitched across nodes
/// and left-aligned, the classes a read falls into for one indel, and the key an indel is counted
/// under, which does not depend on where the aligner put it or on which window counts it.

#include <string>
#include <tuple>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "../../src/unittest/catch.hpp"
#include "../../src/alignment.hpp"
#include "../off_panel.hpp"
#include "../site_values.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// CGTCAAA AAGTCCG TTGCA: a run of five A from the last three bases of node 1 into node 2, with
/// the C before it at node 1 offset 3 and the G after it at node 2 offset 2.
static void make_run_graph(bdsg::HashGraph& graph) {
    handle_t h1 = graph.create_handle("CGTCAAA", 1);
    handle_t h2 = graph.create_handle("AAGTCCG", 2);
    handle_t h3 = graph.create_handle("TTGCA", 3);
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
}

/// An edit: from length, to length, and the read's bases where they differ.
using EditSpec = tuple<int, int, string>;

/// A read along nodes 1, 2 and 3 forward, each from offset 0, with the given edits on each.
static Alignment run_read(const string& name, const vector<vector<EditSpec>>& per_node) {
    Alignment aln;
    aln.set_name(name);
    aln.set_mapping_quality(60);
    const string walk = "CGTCAAAAAGTCCGTTGCA";
    string sequence;
    size_t walk_pos = 0;
    for (size_t n = 0; n < per_node.size(); ++n) {
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(n + 1);
        m->mutable_position()->set_offset(0);
        for (const EditSpec& spec : per_node[n]) {
            Edit* e = m->add_edit();
            e->set_from_length(get<0>(spec));
            e->set_to_length(get<1>(spec));
            e->set_sequence(get<2>(spec));
            if (get<0>(spec) == get<1>(spec) && get<2>(spec).empty()) {
                sequence += walk.substr(walk_pos, get<0>(spec));
            } else {
                sequence += get<2>(spec);
            }
            walk_pos += get<0>(spec);
        }
    }
    aln.set_sequence(sequence);
    return aln;
}

static WalkAlignment along_walk(const Alignment& aln) {
    vector<const Mapping*> mappings;
    vector<int64_t> starts{0, 7, 14};
    for (const Mapping& m : aln.path().mapping()) {
        mappings.push_back(&m);
    }
    return walk_alignment(mappings, starts, 0, aln.sequence());
}

static const string WALK = "CGTCAAAAAGTCCGTTGCA";
static char walk_base(int64_t x) {
    return WALK[x];
}

TEST_CASE("A read's deletion split across two nodes is one indel", "[off_panel]") {
    // The last A of node 1 and the first of node 2.
    const Alignment aln = run_read("r", {{{6, 6, ""}, {1, 0, ""}}, {{1, 0, ""}, {6, 6, ""}},
                                         {{5, 5, ""}}});
    const WalkAlignment walked = along_walk(aln);
    REQUIRE(walked.indels.size() == 1);
    REQUIRE(walked.indels[0].at == 6);
    REQUIRE(walked.indels[0].deleted == 2);
    REQUIRE(walked.begin == 0);
    REQUIRE(walked.end == 19);
}

TEST_CASE("A clip is not an insertion", "[off_panel]") {
    const Alignment aln = run_read("r", {{{0, 2, "TT"}, {7, 7, ""}}, {{7, 7, ""}},
                                         {{5, 5, ""}, {0, 3, "GGG"}}});
    REQUIRE(along_walk(aln).indels.empty());
}

TEST_CASE("Indels are left-aligned to the start of their repeat tract", "[off_panel]") {
    SECTION("a deletion in the run moves to its first A, across the node boundary") {
        const Alignment aln = run_read("r", {{{7, 7, ""}}, {{1, 1, ""}, {1, 0, ""}, {5, 5, ""}},
                                             {{5, 5, ""}}});
        WalkAlignment walked = along_walk(aln);
        REQUIRE(walked.indels.size() == 1);
        REQUIRE(walked.indels[0].at == 8);
        left_align(walked, walk_base);
        REQUIRE(walked.indels[0].at == 4);
        REQUIRE(tract_end(walked.indels[0], 19, walk_base) == 9);
    }
    SECTION("an insertion of an A at the end of the run moves to its start") {
        WalkIndel ins;
        ins.at = 9;
        ins.inserted = "A";
        left_align(ins, 1, walk_base);
        REQUIRE(ins.at == 4);
        REQUIRE(ins.inserted == "A");
        REQUIRE(tract_end(ins, 19, walk_base) == 9);
    }
    SECTION("an indel never moves past the read's first base or an earlier indel") {
        WalkIndel del;
        del.at = 8;
        del.deleted = 1;
        left_align(del, 6, walk_base);
        REQUIRE(del.at == 6);
    }
}

TEST_CASE("A read is REF, ALT or other for an indel only if it covers its tract", "[off_panel]") {
    WalkIndel indel;
    indel.at = 4;
    indel.deleted = 1;
    const int64_t after = 9;
    auto observe = [&](const Alignment& aln) {
        WalkAlignment walked = along_walk(aln);
        left_align(walked, walk_base);
        return observe_indel(walked, indel, after);
    };
    SECTION("no indel across the tract is REF") {
        REQUIRE(observe(run_read("r", {{{7, 7, ""}}, {{7, 7, ""}}, {{5, 5, ""}}}))
                == IndelObservation::REF);
    }
    SECTION("the same deletion, wherever the aligner put it, is ALT") {
        REQUIRE(observe(run_read("r", {{{4, 4, ""}, {1, 0, ""}, {2, 2, ""}}, {{7, 7, ""}},
                                       {{5, 5, ""}}}))
                == IndelObservation::ALT);
        REQUIRE(observe(run_read("r", {{{7, 7, ""}}, {{1, 1, ""}, {1, 0, ""}, {5, 5, ""}},
                                       {{5, 5, ""}}}))
                == IndelObservation::ALT);
    }
    SECTION("a deletion of two A is other") {
        REQUIRE(observe(run_read("r", {{{5, 5, ""}, {2, 0, ""}}, {{7, 7, ""}}, {{5, 5, ""}}}))
                == IndelObservation::OTHER);
    }
    SECTION("a read ending inside the run says nothing") {
        Alignment aln = run_read("r", {{{7, 7, ""}}, {{1, 1, ""}}});
        REQUIRE(observe(aln) == IndelObservation::NONE);
    }
}

TEST_CASE("An indel is counted under its left-aligned anchor, whichever window counts it",
          "[off_panel]") {
    bdsg::HashGraph graph;
    make_run_graph(graph);
    vector<Alignment> reads;
    // Three fragments delete one A, each placed differently, one on the other strand.
    reads.push_back(run_read("a", {{{7, 7, ""}}, {{1, 1, ""}, {1, 0, ""}, {5, 5, ""}},
                                   {{5, 5, ""}}}));
    reads.push_back(run_read("b", {{{4, 4, ""}, {1, 0, ""}, {2, 2, ""}}, {{7, 7, ""}},
                                   {{5, 5, ""}}}));
    reads.push_back(reverse_complement_alignment(
        run_read("c", {{{6, 6, ""}, {1, 0, ""}}, {{7, 7, ""}}, {{5, 5, ""}}}),
        [&](nid_t id) { return (int64_t)graph.get_length(graph.get_handle(id)); }));
    // One fragment reads the REF and one deletes two A.
    reads.push_back(run_read("d", {{{7, 7, ""}}, {{7, 7, ""}}, {{5, 5, ""}}}));
    reads.push_back(run_read("e", {{{5, 5, ""}, {2, 0, ""}}, {{7, 7, ""}}, {{5, 5, ""}}}));
    auto each_read = [&](const function<void(const Alignment&)>& visit) {
        for (const Alignment& aln : reads) {
            visit(aln);
        }
    };
    EditCountParams params;
    params.min_fraction = 0.0;

    auto count_on = [&](const vector<nid_t>& owned) {
        return count_edits(graph, each_read,
                           [&](nid_t id) {
                               return std::find(owned.begin(), owned.end(), id) != owned.end();
                           },
                           params);
    };

    SECTION("the deletion is one candidate, anchored on the C before the run") {
        const vector<EditCandidate> found = count_on({1, 2, 3});
        REQUIRE(found.size() == 1);
        const EditCandidate& c = found.front();
        REQUIRE(c.node == 1);
        REQUIRE(c.offset == 3);
        REQUIRE(c.kind == 'D');
        REQUIRE(c.sequence == "A");
        REQUIRE(c.ref == 'C');
        REQUIRE(c.walk == vector<nid_t>{1, 2});
        REQUIRE(c.span == 6);
        REQUIRE(c.homopolymer == 5);
        REQUIRE(c.alt_fragments == 3);
        REQUIRE(c.ref_fragments == 1);
        REQUIRE(c.other_fragments == 1);
    }
    SECTION("the window holding the anchor counts every placement; the next counts none") {
        // Two of the three reads put the deletion on node 2, which the second window owns.
        const vector<EditCandidate> first = count_on({1});
        const vector<EditCandidate> second = count_on({2, 3});
        REQUIRE(first.size() == 1);
        REQUIRE(first.front().alt_fragments == 3);
        REQUIRE(second.empty());
    }
    SECTION("its edit site deletes the first A of the run, inside node 1") {
        const EditCandidate c = count_on({1, 2, 3}).front();
        const EditSite site = edit_site(graph, {graph.get_handle(1), graph.get_handle(2)},
                                        c.offset, c.offset + 1, c.sequence.size(), {""});
        REQUIRE(graph.get_id(site.alleles.walks.front().front()) == 1);
        REQUIRE(graph.get_id(site.alleles.walks.front().back()) == 2);
        REQUIRE(site.alleles.walks[0] == site.alleles.walks[1]);
        REQUIRE(site.alleles.step_sequences(graph, 1) == StepSequences{"CGTCAA", "AAGTCCG"});
        REQUIRE(site.alleles.sequence(graph, 1)
                == "CGTCAAAAGTCCG");
    }
}

}
}
