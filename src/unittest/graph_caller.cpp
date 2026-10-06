/// \file graph_caller.cpp
/// Tests for the VCF output buffer's record ordering.
///
/// Records sharing a position are rare in real output, so the order of ties is checked here.

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "catch.hpp"
#include "../graph_caller.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("The buffered record order is total, so two runs cannot disagree", "[graph_caller]") {
    using BufferedRecordKey = VCFOutputCaller::BufferedRecordKey;
    const auto buffered_record_key_less = &VCFOutputCaller::buffered_record_key_less;

    // Two records at the same position, distinguished only by the snarl they came from. This is the
    // case (contig, POS) cannot order, and it is not rare: a nested site sits at or near its
    // parent's position, and flattening can put two snarls on one anchor base.
    const BufferedRecordKey a{"chr20", 1000, ">1>5"};
    const BufferedRecordKey b{"chr20", 1000, ">6>9"};

    // Antisymmetry: without the ID in the key, the comparator would return false in both directions
    // for this pair, so `std::sort` could leave them in input order.
    REQUIRE(buffered_record_key_less(a, b));
    REQUIRE_FALSE(buffered_record_key_less(b, a));

    // Irreflexivity, on every field, since a strict weak ordering needs it and a `<=` typo here
    // would make std::sort's behaviour undefined rather than merely wrong.
    REQUIRE_FALSE(buffered_record_key_less(a, a));

    SECTION("the earlier field always dominates the later one") {
        // Position beats ID: a later ID at an earlier position still sorts first, so the tie-break
        // cannot reorder the file.
        REQUIRE(buffered_record_key_less(BufferedRecordKey{"chr20", 999, ">9>9"},
                                         BufferedRecordKey{"chr20", 1000, ">1>1"}));
        // Contig beats both.
        REQUIRE(buffered_record_key_less(BufferedRecordKey{"chr1", 5000, ">9>9"},
                                         BufferedRecordKey{"chr20", 1, ">1>1"}));
    }

    SECTION("every input permutation sorts to one output") {
    // Sorting each permutation of a set containing a tie must give the same sequence every time.
        vector<BufferedRecordKey> keys{
            {"chr20", 1000, ">6>9"},
            {"chr20", 1000, ">1>5"},
            {"chr20", 900, ">2>3"},
            {"chr21", 10, ">4>7"},
        };
        sort(keys.begin(), keys.end(), buffered_record_key_less);
        const vector<string> want{">2>3", ">1>5", ">6>9", ">4>7"};

        // Start from the first permutation in id order, so next_permutation walks all of them
        // rather than reporting exhaustion on its first call.
        vector<BufferedRecordKey> permuted = keys;
        sort(permuted.begin(), permuted.end(),
             [](const BufferedRecordKey& x, const BufferedRecordKey& y) { return x.id < y.id; });
        size_t checked = 0;
        do {
            vector<BufferedRecordKey> copy = permuted;
            sort(copy.begin(), copy.end(), buffered_record_key_less);
            for (size_t i = 0; i < want.size(); ++i) {
                REQUIRE(copy[i].id == want[i]);
            }
            ++checked;
        } while (next_permutation(permuted.begin(), permuted.end(),
                                  [](const BufferedRecordKey& x, const BufferedRecordKey& y) {
                                      return x.id < y.id;
                                  }));
        // All 24 permutations of four distinct records, so the claim is exhaustive rather than
        // sampled.
        REQUIRE(checked == 24);
    }
}

/// A traversal over the given node ids, as plain node visits.
static SnarlTraversal make_trav(const vector<nid_t>& nodes) {
    SnarlTraversal t;
    for (nid_t n : nodes) {
        t.add_visit()->set_node_id(n);
    }
    return t;
}

/// A child chain with the given boundary nodes.
static Snarl make_child(nid_t start, nid_t end) {
    Snarl s;
    s.mutable_start()->set_node_id(start);
    s.mutable_end()->set_node_id(end);
    return s;
}

TEST_CASE("offset_of_child reports where a traversal enters a chain", "[graph_caller]") {
    const SnarlTraversal t = make_trav({1, 2, 3, 4, 5});
    SECTION("the entry index, not the exit") {
        REQUIRE(FlowCaller::offset_of_child(t, make_child(2, 4)) == 1);
    }
    SECTION("entering from either boundary is the same crossing") {
        REQUIRE(FlowCaller::offset_of_child(t, make_child(4, 2)) == 1);
    }
    SECTION("a chain the traversal does not cross has no offset") {
        REQUIRE(FlowCaller::offset_of_child(t, make_child(7, 9)) == -1);
    }
    SECTION("touching one boundary only is not a crossing") {
        REQUIRE(FlowCaller::offset_of_child(t, make_child(3, 99)) == -1);
    }
}


TEST_CASE("ChildOffsets gives base_offset_of_child's answer by lookup", "[graph_caller]") {
    // Node n is n bases long.
    bdsg::HashGraph graph;
    for (nid_t id = 1; id <= 9; ++id) {
        graph.create_handle(string((size_t)id, 'A'), id);
    }
    // A traversal that revisits nodes, with a child-snarl visit in the middle, which counts for
    // neither the entry nor the bases.
    SnarlTraversal t;
    for (nid_t n : {1, 2, 3, 2, 5}) {
        t.add_visit()->set_node_id(n);
    }
    Visit* snarl_visit = t.add_visit();
    snarl_visit->mutable_snarl()->mutable_start()->set_node_id(5);
    snarl_visit->mutable_snarl()->mutable_end()->set_node_id(6);
    for (nid_t n : {6, 3, 7, 2, 9, 9}) {
        t.add_visit()->set_node_id(n);
    }

    // base_offset_of_child's definition: offset_of_child's entry, then the bases of the node
    // visits before it.
    auto expected = [&](const Snarl& child) -> int64_t {
        const int entry = FlowCaller::offset_of_child(t, child);
        if (entry < 0) {
            return -1;
        }
        int64_t bases = 0;
        for (int i = 0; i < entry; ++i) {
            if (!t.visit(i).has_snarl()) {
                bases += (int64_t)graph.get_length(graph.get_handle(t.visit(i).node_id()));
            }
        }
        return bases;
    };

    // Every pair of boundary nodes, including one the traversal never visits (10), both
    // orientations, and a chain that starts and ends on one node.
    const FlowCaller::ChildOffsets offsets(graph, t);
    size_t crossing = 0;
    for (nid_t start = 1; start <= 10; ++start) {
        for (nid_t end = 1; end <= 10; ++end) {
            const Snarl child = make_child(start, end);
            REQUIRE(offsets.base_offset(child) == expected(child));
            crossing += expected(child) >= 0 ? 1 : 0;
        }
    }
    // The comparison is not vacuous: many pairs cross, and some do not.
    REQUIRE(crossing > 20);
    REQUIRE(crossing < 100);
}


/// Exposes the anchor path's strand lookup and the tables it reads.
class StrandLookup : public VCFOutputCaller {
public:
    StrandLookup() : VCFOutputCaller("S") {
        render_lambda_temper = 1.0;
        render_lambda_ceiling = 1.0;
    }
    void add_read(const string& name, size_t phase_set, bool multi) {
        ReadLambda read;
        read.lambda = 2.0;
        read.sites = 1;
        read.phase_set = phase_set;
        read.multi_phase_set = multi;
        render_lambda[(uint64_t)std::hash<string>{}(name)] = read;
    }
    void set_site_phase_set(size_t record_key, size_t phase_set) {
        render_lambda_phase_set[record_key] = phase_set;
    }
    using VCFOutputCaller::read_strand_log_odds;
};

TEST_CASE("A strand from another phase set is NaN in the anchor path, and no strand is 0",
          "[graph_caller]") {
    // --anchors-hom-split drops a NaN read from a split site and places a 0 read by a coin, so the
    // two must stay distinct here; re-genotyping, which does not use this lookup, gives both 0.
    StrandLookup caller;
    caller.add_read("here", 7, false);
    caller.add_read("elsewhere", 9, false);
    caller.add_read("both", 7, true);
    caller.set_site_phase_set(1, 7);

    const double here = caller.read_strand_log_odds(1, "here");
    REQUIRE(std::isfinite(here));
    REQUIRE(here > 0.0);
    REQUIRE(std::isnan(caller.read_strand_log_odds(1, "elsewhere")));
    REQUIRE(std::isnan(caller.read_strand_log_odds(1, "both")));
    REQUIRE(caller.read_strand_log_odds(1, "unseen") == 0.0);

    // A site with no phase set takes any read found in one phase set, but not one found in two.
    REQUIRE(caller.read_strand_log_odds(2, "elsewhere") > 0.0);
    REQUIRE(std::isnan(caller.read_strand_log_odds(2, "both")));
}

TEST_CASE("A parent's phase swap carries its nested strands and their haplotypes with it",
          "[graph_caller]") {
    const size_t W = LinkageModel::WILDCARD;
    // Diploid parent 1, swapped by read phasing. Its ploidy-1 child 2 is on strand 1 with
    // haplotype 7, and 2's own ploidy-1 child 3 is on strand 0 with haplotype 5. Diploid child 4 is
    // not swapped, so its ploidy-1 child 5 keeps its strand.
    auto call = [](size_t key, size_t ploidy, int strand, size_t first, size_t second) {
        LinkageCollector::PhaseCall pc;
        pc.record_key = key;
        pc.ploidy = ploidy;
        pc.nested_strand = (int8_t)strand;
        pc.hap_first = first;
        pc.hap_second = second;
        return pc;
    };
    vector<LinkageCollector::PhaseCall> phased = {
        call(1, 2, -1, 3, 4), call(2, 1, 1, W, 7), call(3, 1, 0, 5, W),
        call(4, 2, -1, 3, 4), call(5, 1, 0, 6, W)};
    std::unordered_map<size_t, size_t> index;
    for (size_t i = 0; i < phased.size(); ++i) {
        index[phased[i].record_key] = i;
    }
    // Out of level order, as the staged sites can be.
    vector<FlowCaller::NestedLink> links = {
        {3, 2, 2}, {5, 4, 2}, {2, 1, 1}, {4, 1, 1}, {1, 0, 0}};
    const unordered_set<size_t> flips = {1};

    SECTION("each strand moves, and its haplotype moves to the slot the strand names") {
        REQUIRE(FlowCaller::cascade_nested_strands(phased, index, links, flips) == 2);
        REQUIRE(phased[1].nested_strand == 0);
        REQUIRE(phased[1].hap_first == 7);
        REQUIRE(phased[1].hap_second == W);
        REQUIRE(phased[2].nested_strand == 1);
        REQUIRE(phased[2].hap_first == W);
        REQUIRE(phased[2].hap_second == 5);
        REQUIRE(phased[4].nested_strand == 0);
        REQUIRE(phased[4].hap_first == 6);
    }
    SECTION("a chain left out of the links stops the swap reaching its children") {
        links.erase(links.begin() + 2);   // chain 2
        REQUIRE(FlowCaller::cascade_nested_strands(phased, index, links, flips) == 0);
        REQUIRE(phased[2].nested_strand == 0);
    }
}

}
}
