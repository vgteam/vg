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
#include "../flow_caller.hpp"
#include "../vcf_output_caller.hpp"
#include "../child_placer.hpp"
#include "../phase_table.hpp"
#include "../read_strand_table.hpp"
#include "../staged_site.hpp"
#include "../vcf_record.hpp"

namespace vg {
namespace unittest {

using namespace vg::multipass;

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

/// A graph with nodes 1 to 10 and 99, node n being n bases long.
static void add_numbered_nodes(bdsg::HashGraph& graph) {
    for (nid_t id : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 99}) {
        graph.create_handle(string((size_t)id, 'A'), id);
    }
}

/// A walk over the given node ids, each read forward.
static Traversal make_walk(const HandleGraph& graph, const vector<nid_t>& nodes) {
    Traversal walk;
    for (nid_t n : nodes) {
        walk.push_back(graph.get_handle(n));
    }
    return walk;
}

/// A child site with the given boundary nodes.
static SiteBounds make_child(const HandleGraph& graph, nid_t start, nid_t end) {
    return SiteBounds{graph.get_handle(start), graph.get_handle(end)};
}

TEST_CASE("offset_of_child reports where a walk enters a chain", "[graph_caller]") {
    bdsg::HashGraph graph;
    add_numbered_nodes(graph);
    const Traversal t = make_walk(graph, {1, 2, 3, 4, 5});
    SECTION("the entry index, not the exit") {
        REQUIRE(ChildPlacer::offset_of_child(graph, t, make_child(graph, 2, 4)) == 1);
    }
    SECTION("entering from either boundary is the same crossing") {
        REQUIRE(ChildPlacer::offset_of_child(graph, t, make_child(graph, 4, 2)) == 1);
    }
    SECTION("a chain the walk does not cross has no offset") {
        REQUIRE(ChildPlacer::offset_of_child(graph, t, make_child(graph, 7, 9)) == -1);
    }
    SECTION("touching one boundary only is not a crossing") {
        REQUIRE(ChildPlacer::offset_of_child(graph, t, make_child(graph, 3, 99)) == -1);
    }
}


TEST_CASE("ChildOffsets gives base_offset_of_child's answer by lookup", "[graph_caller]") {
    bdsg::HashGraph graph;
    add_numbered_nodes(graph);
    // A walk that revisits nodes.
    const Traversal t = make_walk(graph, {1, 2, 3, 2, 5, 6, 3, 7, 2, 9, 9});

    // base_offset_of_child's definition: offset_of_child's entry, then the bases of the handles
    // before it.
    auto expected = [&](const SiteBounds& child) -> int64_t {
        const int entry = ChildPlacer::offset_of_child(graph, t, child);
        if (entry < 0) {
            return -1;
        }
        int64_t bases = 0;
        for (int i = 0; i < entry; ++i) {
            bases += (int64_t)graph.get_length(t[i]);
        }
        return bases;
    };

    // Every pair of boundary nodes, including one the walk never visits (10), both orientations,
    // and a chain that starts and ends on one node.
    const ChildPlacer::ChildOffsets offsets(graph, t);
    size_t crossing = 0;
    for (nid_t start = 1; start <= 10; ++start) {
        for (nid_t end = 1; end <= 10; ++end) {
            const SiteBounds child = make_child(graph, start, end);
            REQUIRE(offsets.base_offset(graph, child) == expected(child));
            crossing += expected(child) >= 0 ? 1 : 0;
        }
    }
    // The comparison is not vacuous: many pairs cross, and some do not.
    REQUIRE(crossing > 20);
    REQUIRE(crossing < 100);
}


TEST_CASE("A site inside a node is crossed by each visit to it, and a cyclic site by a second visit",
          "[graph_caller]") {
    bdsg::HashGraph graph;
    add_numbered_nodes(graph);
    // Both are bounded by node 2 read forward at each end. The cyclic site lies between the
    // node's sides, outside it; the other lies inside it, between the sides of its bases.
    const SiteBounds inside{graph.get_handle(2), graph.get_handle(2), true};
    const SiteBounds cyclic = make_child(graph, 2, 2);
    REQUIRE(inside != cyclic);
    auto crossings = [&](const Traversal& walk, const SiteBounds& child) {
        return ChildPlacer::crossings_of_child(graph, ChildPlacer::index_traversal_nodes(graph, walk),
                                               child);
    };

    SECTION("one visit crosses the site inside the node, and not the cyclic site") {
        const Traversal once = make_walk(graph, {1, 2, 3});
        REQUIRE(crossings(once, inside) == 1);
        REQUIRE(crossings(once, cyclic) == 0);
        REQUIRE(ChildPlacer::offset_of_child(graph, once, inside) == 1);
        REQUIRE(ChildPlacer::offset_of_child(graph, once, cyclic) == -1);
        const ChildPlacer::ChildOffsets offsets(graph, once);
        REQUIRE(offsets.base_offset(graph, inside) == 1);
        REQUIRE(offsets.base_offset(graph, cyclic) == -1);
    }
    SECTION("two visits cross the site inside the node twice, and the cyclic site once") {
        const Traversal twice = make_walk(graph, {1, 2, 3, 2, 5});
        REQUIRE(crossings(twice, inside) == 2);
        REQUIRE(crossings(twice, cyclic) == 1);
        // Both are entered at the first visit.
        REQUIRE(ChildPlacer::offset_of_child(graph, twice, inside) == 1);
        REQUIRE(ChildPlacer::offset_of_child(graph, twice, cyclic) == 1);
    }
    SECTION("a visit read backward crosses the site inside the node too") {
        Traversal backward = make_walk(graph, {1, 3});
        backward.insert(backward.begin() + 1, graph.get_handle(2, true));
        REQUIRE(crossings(backward, inside) == 1);
        REQUIRE(ChildPlacer::offset_of_child(graph, backward, inside) == 1);
    }
    SECTION("a walk that misses the node does not cross it") {
        const Traversal missing = make_walk(graph, {1, 3, 4});
        REQUIRE(crossings(missing, inside) == 0);
        REQUIRE(ChildPlacer::offset_of_child(graph, missing, inside) == -1);
        REQUIRE(ChildPlacer::ChildOffsets(graph, missing).base_offset(graph, inside) == -1);
    }
    SECTION("the crossing mask sets the bit of each walk that visits the node") {
        const vector<ChildPlacer::TraversalNodeIndex> visits{
            ChildPlacer::index_traversal_nodes(graph, make_walk(graph, {1, 2, 3})),
            ChildPlacer::index_traversal_nodes(graph, make_walk(graph, {1, 3})),
            ChildPlacer::index_traversal_nodes(graph, make_walk(graph, {1, 2, 3, 2, 5}))};
        bool known = false;
        REQUIRE(ChildPlacer::child_crossing_mask(graph, visits, inside, &known) == 0b101);
        REQUIRE(known);
        REQUIRE(ChildPlacer::child_crossing_mask(graph, visits, cyclic, &known) == 0b100);
    }
}


TEST_CASE("The staged-site tree places each site under the sites enclosing it, once each",
          "[graph_caller]") {
    bdsg::HashGraph graph;
    add_numbered_nodes(graph);
    auto bounds = [&](nid_t start, bool start_back, nid_t end, bool end_back) {
        return SiteBounds{graph.get_handle(start, start_back), graph.get_handle(end, end_back)};
    };
    auto name = [&](const SiteBounds& site) {
        return site_name(graph, site.start, site.end, nullptr, false);
    };
    StagedSiteTable table;
    table.start(1);
    auto stage = [&](const SiteBounds& site, const string& id, vector<SiteBounds> enclosing) {
        StagedSite staged;
        staged.bounds = site;
        staged.id = id;
        staged.enclosing = std::move(enclosing);
        table.add_top_level(std::move(staged));
    };
    // A top-level site 1..8, and inside it 2..7, staged read backward as the reference might
    // orient it. Inside that, 3..5 is not staged, and inside 3..5 a site inside node 4 has an ID
    // of its own.
    const SiteBounds outer = bounds(1, false, 8, false);
    const SiteBounds middle = bounds(2, false, 7, false);
    const SiteBounds unstaged = bounds(3, false, 5, false);
    SiteBounds inside = bounds(4, false, 4, false);
    inside.inside_node = true;
    stage(outer, name(outer), {});
    stage(bounds(7, true, 2, true), "<7<2", {outer});
    stage(inside, "4:1:G", {unstaged, middle, outer});

    const StagedSiteTree tree(table, graph, name);
    vector<SiteTree::site_t> order;
    tree.for_each_site([&](SiteTree::site_t site) { order.push_back(site); }, true);
    // The middle site is one node, though it was staged one way round and encloses the other.
    REQUIRE(order.size() == 4);
    for (size_t i = 0; i < order.size(); ++i) {
        const SiteTree::site_t parent = tree.parent_of(order[i]);
        // Preorder: each site comes after the site enclosing it.
        REQUIRE((parent == nullptr
                 || std::find(order.begin(), order.begin() + i, parent) != order.begin() + i));
    }
    SiteTree::site_t own = nullptr;
    for (SiteTree::site_t site : order) {
        if (tree.id_of(site) != nullptr) {
            REQUIRE(own == nullptr);
            own = site;
        }
    }
    REQUIRE(own != nullptr);
    REQUIRE(*tree.id_of(own) == "4:1:G");
    // Its ancestors, innermost first, are named by their boundary visits. The middle site keeps
    // the orientation it was first added in, which the nesting tags read either way round.
    vector<string> ancestors;
    for (SiteTree::site_t site = tree.parent_of(own); site != nullptr; site = tree.parent_of(site)) {
        REQUIRE(tree.id_of(site) == nullptr);
        const SiteEnds ends = tree.ends_of(site);
        ancestors.push_back(site_name(ends.start_id, ends.start_backward, ends.end_id,
                                      ends.end_backward, nullptr, false));
    }
    REQUIRE(ancestors == vector<string>({">3>5", "<7<2", ">1>8"}));
}


/// Sets the tables the anchor path's strand lookup reads.
class StrandLookup : public ReadStrandTable {
public:
    StrandLookup() {
        lambda_temper = 1.0;
        lambda_ceiling = 1.0;
    }
    void add_read(const string& name, size_t phase_set, bool multi) {
        ReadLambda read;
        read.lambda = 2.0;
        read.sites = 1;
        read.phase_set = phase_set;
        read.multi_phase_set = multi;
        lambda[(uint64_t)std::hash<string>{}(name)] = read;
    }
    void set_site_phase_set(size_t record_key, size_t phase_set) {
        lambda_phase_set[record_key] = phase_set;
    }
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
    const size_t W = LinkageModel::NO_HAPLOTYPE;
    // Diploid parent 1, swapped by read phasing. Its ploidy-1 child 2 is on strand 1 with
    // haplotype 7, and 2's own ploidy-1 child 3 is on strand 0 with haplotype 5. Diploid child 4 is
    // not swapped, so its ploidy-1 child 5 keeps its strand.
    auto call = [](size_t key, size_t ploidy, int strand, size_t first, size_t second) {
        PhaseTable::PhaseCall pc;
        pc.record_key = key;
        pc.ploidy = ploidy;
        pc.nested_strand = (int8_t)strand;
        pc.hap_first = first;
        pc.hap_second = second;
        return pc;
    };
    PhaseTable table;
    vector<PhaseTable::PhaseCall>& phased = table.calls();
    phased = {
        call(1, 2, -1, 3, 4), call(2, 1, 1, W, 7), call(3, 1, 0, 5, W),
        call(4, 2, -1, 3, 4), call(5, 1, 0, 6, W)};
    // Out of level order, as the staged sites can be.
    vector<PhaseTable::NestedLink> links = {
        {3, 2, 2}, {5, 4, 2}, {2, 1, 1}, {4, 1, 1}, {1, 0, 0}};
    const unordered_set<size_t> flips = {1};

    SECTION("each strand moves, and its haplotype moves to the slot the strand names") {
        REQUIRE(table.swap_strands(flips, links) == 2);
        // The swapped parent's own haplotypes change strands.
        REQUIRE(phased[0].hap_first == 4);
        REQUIRE(phased[0].hap_second == 3);
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
        REQUIRE(table.swap_strands(flips, links) == 0);
        REQUIRE(phased[2].nested_strand == 0);
    }
}

}
}
