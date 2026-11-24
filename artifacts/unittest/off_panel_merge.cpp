/// \file off_panel_merge.cpp
///
/// Unit tests for merging touching off-panel candidates on one node into one candidate, whose ALTs
/// are what the reads spell over the bases the candidates change.

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "../../src/unittest/catch.hpp"
#include "../off_panel.hpp"

namespace vg {
namespace unittest {

using namespace std;

static const string NODE1 = "GATCCTGACGTTAGCATGCA";
static const string NODE2 = "CTGTCAGCTA";

/// `count` reads along nodes 1 and 2, named from `prefix`, each with the walk bases in `changes`
/// read as the bases given: one for a substitution, none for a deletion, and more for an insertion
/// after the first.
static vector<Alignment> reads_of(const string& prefix, size_t count,
                                  const map<size_t, string>& changes) {
    const string walk = NODE1 + NODE2;
    Alignment aln;
    aln.set_mapping_quality(60);
    string sequence;
    for (size_t x = 0; x < walk.size(); ++x) {
        if (x == 0 || x == NODE1.size()) {
            aln.mutable_path()->add_mapping()->mutable_position()->set_node_id(x == 0 ? 1 : 2);
        }
        Mapping* m = aln.mutable_path()->mutable_mapping(aln.path().mapping_size() - 1);
        auto found = changes.find(x);
        const string bases = found == changes.end() ? walk.substr(x, 1) : found->second;
        Edit* e = m->add_edit();
        e->set_from_length(1);
        e->set_to_length(bases.empty() ? 0 : 1);
        if (!bases.empty() && bases[0] != walk[x]) {
            e->set_sequence(bases.substr(0, 1));
        }
        if (bases.size() > 1) {
            e = m->add_edit();
            e->set_to_length(bases.size() - 1);
            e->set_sequence(bases.substr(1));
        }
        sequence += bases;
    }
    aln.set_sequence(sequence);
    vector<Alignment> reads(count, aln);
    for (size_t i = 0; i < count; ++i) {
        reads[i].set_name(prefix + to_string(i));
    }
    return reads;
}

TEST_CASE("Touching off-panel changes on one node are one candidate, with the ALTs the reads spell",
          "[off_panel]") {
    bdsg::HashGraph graph;
    graph.create_edge(graph.create_handle(NODE1, 1), graph.create_handle(NODE2, 2));
    EditCallBar bar;
    auto count = [&](const vector<vector<Alignment>>& haplotypes, bool merge = true) {
        vector<Alignment> reads;
        for (const auto& haplotype : haplotypes) {
            reads.insert(reads.end(), haplotype.begin(), haplotype.end());
        }
        EditCountParams params;
        params.bar = merge ? &bar : nullptr;
        return count_edits(graph,
                           [&](const function<void(const Alignment&)>& visit) {
                               for (const Alignment& aln : reads) {
                                   visit(aln);
                               }
                           },
                           [](nid_t) { return true; }, params);
    };
    // Node 1 is GATCCTGA C G TTAGCATGCA: C at 8 and G at 9.
    const vector<Alignment> ref = reads_of("r", 4, {});

    SECTION("two touching substitutions, one haplotype with both and one with the second alone") {
        const auto both = reads_of("x", 8, {{8, "A"}, {9, "T"}});
        const auto second = reads_of("y", 8, {{9, "T"}});
        const vector<EditCandidate> found = count({both, second, ref});
        REQUIRE(found.size() == 1);
        const EditCandidate& c = found.front();
        REQUIRE(c.kind == 'M');
        REQUIRE(c.node == 1);
        REQUIRE(c.offset == 8);
        REQUIRE(c.span == 2);
        REQUIRE(c.sequence == "CG");
        REQUIRE(c.alts == vector<string>{"AT", "CT"});
        REQUIRE(c.alt_fragments == 16);
        REQUIRE(c.ref_fragments == 4);
        REQUIRE(c.walk == vector<nid_t>{1});

        // Whatever order the reads come in.
        const vector<EditCandidate> reordered = count({ref, second, both});
        REQUIRE(reordered.size() == 1);
        REQUIRE(reordered.front().alts == c.alts);
        REQUIRE(reordered.front().alt_fragments == c.alt_fragments);
        // Without the bar, each is a candidate of its own.
        REQUIRE(count({both, second, ref}, false).size() == 2);
    }
    SECTION("a fragment whose reads spell two strings is counted under neither") {
        auto discordant = reads_of("d", 1, {{8, "A"}, {9, "T"}});
        discordant.push_back(reads_of("d", 1, {}).front());
        const vector<EditCandidate> found =
            count({reads_of("x", 8, {{8, "A"}, {9, "T"}}), reads_of("y", 8, {{9, "T"}}), ref,
                   discordant});
        REQUIRE(found.size() == 1);
        REQUIRE(found.front().alt_fragments == 16);
        REQUIRE(found.front().ref_fragments == 4);
        REQUIRE(found.front().discordant_fragments == 1);
    }
    SECTION("a substitution and an insertion after it are one complex change") {
        const vector<EditCandidate> found = count({reads_of("x", 8, {{8, "AA"}}), ref});
        REQUIRE(found.size() == 1);
        REQUIRE(found.front().kind == 'M');
        REQUIRE(found.front().sequence == "C");
        REQUIRE(found.front().alts == vector<string>{"AA"});
    }
    SECTION("two substitutions of one base are one change with two ALTs") {
        const vector<EditCandidate> found =
            count({reads_of("x", 8, {{8, "A"}}), reads_of("y", 8, {{8, "T"}}), ref});
        REQUIRE(found.size() == 1);
        REQUIRE(found.front().alts == vector<string>{"A", "T"});
    }
    SECTION("a change the call bar does not admit is not merged") {
        const vector<EditCandidate> found =
            count({reads_of("x", 8, {{8, "A"}}), reads_of("y", 4, {{9, "T"}}), ref,
                   reads_of("s", 4, {})});
        REQUIRE(found.size() == 2);
        REQUIRE(found[0].kind == 'S');
        REQUIRE(found[1].kind == 'S');
    }
    SECTION("changes that touch across a node end are not merged") {
        const vector<EditCandidate> found = count({reads_of("x", 8, {{19, "T"}, {20, "A"}}), ref});
        REQUIRE(found.size() == 2);
        REQUIRE(found[0].node == 1);
        REQUIRE(found[1].node == 2);
    }
}

}
}
