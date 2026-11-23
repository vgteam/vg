/// \file allele_splicer.cpp
///
/// Unit tests for the sequences AlleleSplicer resolves for the strands of a site and for its
/// reference allele.

#include "../../src/unittest/catch.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "../allele_splicer.hpp"
#include "support/decomposition_pair.hpp"
#include "support/nested_graph.hpp"
#include "support/staged_nested_sites.hpp"
#include "../site_values.hpp"

namespace vg {
namespace unittest {

using namespace vg::multipass;

TEST_CASE("AlleleSplicer resolves each strand through the sites nested under a site",
          "[allele_splicer]") {
    StagedNestedSites fixture;
    SECTION("A homozygous reference site takes each strand's alleles at the sites under it") {
        // The top site carries its reference walk on both strands. The middle site carries the
        // reference on strand 0 and node 6 on strand 1. The inner site is on strand 0 alone,
        // where it leaves out node 4.
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 1);
        fixture.phase(fixture.inner_key, 1, -1, 1, 0);
        fixture.index();
        const AlleleSplicer& splicer = fixture.splicer;

        REQUIRE(splicer.reference_sequence(fixture.site(fixture.top_key)) == "GCATGCTGAGCAGCTGA");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 0) == "GCATGGCAGCTGA");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 1) == "GCATTGCTGA");
        REQUIRE(splicer.strand_traversal(fixture.site(fixture.top_key), 1) == 0);

        // The nested sites resolve the same way.
        REQUIRE(splicer.reference_sequence(fixture.site(fixture.middle_key)) == "TGCTGAGCAG");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.middle_key), 0) == "TGGCAG");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.middle_key), 1) == "TTG");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.inner_key), 0) == "GGCA");
        // Strand 1 does not cross the inner site.
        REQUIRE(!splicer.strand_sequence(fixture.site(fixture.inner_key), 1));
        REQUIRE(!splicer.strand_traversal(fixture.site(fixture.inner_key), 1));

        // The record of the top site spells no allele differently, since both strands carry the
        // reference.
        AlleleSplicer::SplicedAlleles spliced;
        REQUIRE(!splicer.splice(fixture.site(fixture.top_key), spliced));
    }
    SECTION("A strand that skips a nested site takes nothing from it") {
        // Strand 0 of the top site goes from node 1 to node 8. The middle site is on strand 1.
        fixture.phase(fixture.top_key, 1, 0);
        fixture.phase(fixture.middle_key, 1, -1, 1, 1);
        fixture.index();
        const AlleleSplicer& splicer = fixture.splicer;
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 0) == "GCACTGA");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 1) == "GCATTGCTGA");
        REQUIRE(splicer.reference_sequence(fixture.site(fixture.top_key)) == "GCATGCTGAGCAGCTGA");
        REQUIRE(!splicer.strand_sequence(fixture.site(fixture.middle_key), 0));
    }
    SECTION("A change inside a node reaches the strand that carries it") {
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 0);
        fixture.phase(fixture.inner_key, 0, 0);
        fixture.phase(fixture.edit_key, 1, 0);
        fixture.index();
        const AlleleSplicer& splicer = fixture.splicer;
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 0) == "GCATGCAGAGCAGCTGA");
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.top_key), 1) == "GCATGCTGAGCAGCTGA");
        REQUIRE(splicer.reference_sequence(fixture.site(fixture.top_key)) == "GCATGCTGAGCAGCTGA");
    }
    SECTION("A site with no phase has no strands") {
        fixture.phase(fixture.middle_key, 0, 1);
        fixture.index();
        REQUIRE(!fixture.splicer.strand_sequence(fixture.site(fixture.top_key), 0));
        REQUIRE(fixture.splicer.reference_sequence(fixture.site(fixture.top_key))
                == "GCATGCTGAGCAGCTGA");
    }
}

/// A site, 1 to 9, that holds a chain of two sites that share a bound node: 2 to 5, which takes
/// node 3 or node 4, and 5 to 8, which takes node 6 or node 7. An edge from 1 to 9 leaves the
/// chain out.
///
/// The site has the reference walk through nodes 3 and 6, the walk from 1 straight to 9, and the
/// walk through nodes 4 and 7. Each site of the chain has its reference walk, through node 3 or
/// node 6, and its walk through the other node.
///
/// Node 2 is GGAGGAGACT, and two sites inside it are staged under the top site. One has an
/// allele without the node's first three bases, and the other an allele without its last four.
struct StagedChainSites {
    VG graph;
    vector<Node*> nodes;
    unique_ptr<IntegratedSnarlFinder> finder;
    unique_ptr<DecompositionPair> pair;
    unique_ptr<EditSnarlDecompositionOverlay> overlay;
    PhaseTable phases;
    StagedSiteTable table;
    unique_ptr<AlleleSplicer> splicer;

    const size_t top_key = 1;
    const size_t first_key = 2;
    const size_t second_key = 3;
    const size_t front_edit_key = 4;
    const size_t back_edit_key = 5;

    StagedChainSites() {
        for (const string& sequence :
             {"GCATTACGGATCA", "GGAGGAGACT", "A", "C", "GG", "T", "A", "CC", "TGACCATGCATGA"}) {
            nodes.push_back(graph.create_node(sequence));
        }
        for (const auto& [from, to] : vector<std::pair<int, int>>{
                 {1, 2}, {1, 9}, {2, 3}, {2, 4}, {3, 5}, {4, 5}, {5, 6}, {5, 7}, {6, 8}, {7, 8},
                 {8, 9}}) {
            graph.create_edge(nodes[from - 1], nodes[to - 1]);
        }
        finder = make_unique<IntegratedSnarlFinder>(graph);
        pair = make_unique<DecompositionPair>(graph, *finder);
        overlay = make_unique<EditSnarlDecompositionOverlay>(pair->adapter, graph);
        splicer = make_unique<AlleleSplicer>(graph, *overlay, phases);
        table.start(1);
        table.add_top_level(
            snarl(top_key, 0, {{1, 2, 3, 5, 6, 8, 9}, {1, 9}, {1, 2, 4, 5, 7, 8, 9}}));
        table.add_nested(snarl(first_key, top_key, {{2, 3, 5}, {2, 4, 5}}));
        table.add_nested(snarl(second_key, top_key, {{5, 6, 8}, {5, 7, 8}}));
        table.add_nested(deletion(front_edit_key, 0, 3));
        table.add_nested(deletion(back_edit_key, 6, 4));
    }

    /// Give the site with record key `key` a phase with traversal `first` on strand 0 and
    /// `second` on strand 1, or at `ploidy` 1 `first` alone, on strand `nested_strand` of its
    /// parent.
    void phase(size_t key, int first, int second, size_t ploidy = 2, int nested_strand = -1) {
        PhaseTable::PhaseCall call;
        call.record_key = key;
        call.contig = "ref";
        call.trav_first = first;
        call.trav_second = second;
        call.ploidy = ploidy;
        call.nested_strand = (int8_t)nested_strand;
        phases.calls().push_back(call);
    }

    /// Hand the nested sites to the render queue and index them, as the render does.
    void index() {
        phases.freeze_for_render(true);
        table.gather_nested();
        table.hand_off();
        splicer->index(table, true);
    }

    /// Get the staged site with record key `key`, once `index` has run.
    StagedSite& site(size_t key) {
        StagedSite* found = nullptr;
        table.for_each([&](StagedSite& staged) {
            if (staged.record_key == key) {
                found = &staged;
            }
        });
        if (found == nullptr) {
            throw std::out_of_range("no staged site has the record key");
        }
        return *found;
    }

    /// Return the bases of the nodes `ids`, in order.
    string spell(const vector<nid_t>& ids) const {
        string bases;
        for (const nid_t id : ids) {
            bases += graph.get_sequence(graph.get_handle(id));
        }
        return bases;
    }

private:
    /// Make the staged site inside node 2 whose one ALT is the node without the `removed` bases
    /// from offset `at` on.
    StagedSite deletion(size_t key, size_t at, size_t removed) {
        const handle_t node = graph.get_handle(2);
        StagedSite site;
        site.id = "edit" + std::to_string(key);
        site.record_key = key;
        site.placement.parent_record_key = top_key;
        site.net = overlay->add_edit(node, at, node, at + removed);
        site.alleles.walks = {{node}, {node}};
        site.alleles.ref_allele = 0;
        site.alleles.edits = {{}, {AlleleEdit{0, at, removed, ""}}};
        return site;
    }

    /// Make the staged site of the snarl that the walks `walks` cross, the first being the
    /// reference.
    StagedSite snarl(size_t key, size_t parent_key, const vector<vector<nid_t>>& walks) {
        StagedSite site;
        const nid_t start = walks.front().front();
        const nid_t end = walks.front().back();
        site.id = ">" + std::to_string(start) + ">" + std::to_string(end);
        site.record_key = key;
        site.placement.parent_record_key = parent_key;
        site.net = find_snarl(*overlay, graph, start, end);
        for (const vector<nid_t>& walk : walks) {
            Traversal visits;
            for (const nid_t id : walk) {
                visits.push_back(graph.get_handle(id));
            }
            site.alleles.walks.push_back(std::move(visits));
        }
        site.alleles.ref_allele = 0;
        return site;
    }
};

TEST_CASE("AlleleSplicer splices both of two nested sites that share a bound node",
          "[allele_splicer]") {
    StagedChainSites fixture;
    SECTION("A strand takes the chosen allele of each site of a chain") {
        // Strand 1 of the top site carries an allele given as a sequence, which crosses neither
        // site of the chain. Strand 0 carries the walk through nodes 4 and 7, and the two sites
        // of the chain, which are on that strand alone, carry their reference alleles, through
        // nodes 3 and 6.
        fixture.phase(fixture.top_key, 2, 3);
        fixture.phase(fixture.first_key, 0, -1, 1, 0);
        fixture.phase(fixture.second_key, 0, -1, 1, 0);
        fixture.index();
        StagedSite& top = fixture.site(fixture.top_key);
        REQUIRE(add_sequence_allele(fixture.graph, top.alleles, "TTTTTT") == 3);

        REQUIRE(fixture.splicer->strand_sequence(top, 0) == fixture.spell({1, 2, 3, 5, 6, 8, 9}));
        REQUIRE(fixture.splicer->strand_sequence(top, 1)
                == fixture.spell({1}) + "TTTTTT" + fixture.spell({9}));
        // The record spells the allele of strand 0 the same way, each node once.
        AlleleSplicer::SplicedAlleles spliced;
        REQUIRE(fixture.splicer->splice(top, spliced));
        REQUIRE(spliced.travs[2] == top.alleles.walks[0]);
        REQUIRE(allele_sequence(fixture.graph, spliced.travs[2], &spliced.sequences[2])
                == fixture.spell({1, 2, 3, 5, 6, 8, 9}));
    }
    SECTION("Only the second site differs from its parent's walk") {
        fixture.phase(fixture.top_key, 2, 1);
        fixture.phase(fixture.first_key, 1, -1, 1, 0);
        fixture.phase(fixture.second_key, 0, -1, 1, 0);
        fixture.index();
        REQUIRE(fixture.splicer->strand_sequence(fixture.site(fixture.top_key), 0)
                == fixture.spell({1, 2, 4, 5, 6, 8, 9}));
        REQUIRE(fixture.splicer->strand_sequence(fixture.site(fixture.top_key), 1)
                == fixture.spell({1, 9}));
    }
    SECTION("Each site of the chain carries an allele given as a sequence") {
        // The bases of each sequence allele follow the site's start bound, which for the second
        // site is the node the two share.
        fixture.phase(fixture.top_key, 2, 1);
        fixture.phase(fixture.first_key, 2, -1, 1, 0);
        fixture.phase(fixture.second_key, 2, -1, 1, 0);
        fixture.index();
        CandidateAlleles& first = fixture.site(fixture.first_key).alleles;
        CandidateAlleles& second = fixture.site(fixture.second_key).alleles;
        REQUIRE(add_sequence_allele(fixture.graph, first, "AAAA") == 2);
        REQUIRE(add_sequence_allele(fixture.graph, second, "CCCCC") == 2);
        REQUIRE(fixture.splicer->strand_sequence(fixture.site(fixture.top_key), 0)
                == fixture.spell({1, 2}) + "AAAA" + fixture.spell({5}) + "CCCCC"
                       + fixture.spell({8, 9}));
    }
    SECTION("Two changes inside one node that a repeat could slide together are both made") {
        // Without its first three bases, node 2 is GGAGACT, which is also the node without its
        // fifth to seventh bases, and those reach the bases the other change removes.
        fixture.phase(fixture.top_key, 2, 1);
        fixture.phase(fixture.first_key, 1, -1, 1, 0);
        fixture.phase(fixture.second_key, 1, -1, 1, 0);
        fixture.phase(fixture.front_edit_key, 1, -1, 1, 0);
        fixture.phase(fixture.back_edit_key, 1, -1, 1, 0);
        fixture.index();
        REQUIRE(fixture.splicer->strand_sequence(fixture.site(fixture.top_key), 0)
                == fixture.spell({1}) + "GGA" + fixture.spell({4, 5, 7, 8, 9}));
        REQUIRE(fixture.splicer->strand_sequence(fixture.site(fixture.top_key), 1)
                == fixture.spell({1, 9}));
    }
    SECTION("Both strands cross the chain, each with its own alleles") {
        fixture.phase(fixture.top_key, 2, 2);
        fixture.phase(fixture.first_key, 0, 1);
        fixture.phase(fixture.second_key, 1, 0);
        fixture.index();
        REQUIRE(fixture.splicer->strand_sequence(fixture.site(fixture.top_key), 0)
                == fixture.spell({1, 2, 3, 5, 7, 8, 9}));
        REQUIRE(fixture.splicer->strand_sequence(fixture.site(fixture.top_key), 1)
                == fixture.spell({1, 2, 4, 5, 6, 8, 9}));
    }
}

TEST_CASE("AlleleSplicer makes an edit that overlaps only an edit that was skipped",
          "[allele_splicer]") {
    // Three sites change one node of two bases: one deletes both, and two change one base
    // each. The edits are made from the right, and one is skipped where it takes a base that an
    // edit already made took.
    StepSequences steps{"AA"};
    vector<AlleleSplicer::StepEdit> edits{{0, 0, 2, "", 1}, {0, 0, 1, "G", 2}, {0, 1, 1, "C", 3}};
    REQUIRE(AlleleSplicer::apply_step_edits(edits, steps));
    // The deletion overlaps the change of the second base and is skipped. The change of the
    // first base overlaps the deletion alone, which was not made.
    REQUIRE(steps == StepSequences{"GC"});

    SECTION("of two edits that share a base, the one further right is made") {
        StepSequences shared{"ACGT"};
        REQUIRE(AlleleSplicer::apply_step_edits({{0, 0, 2, "X", 1}, {0, 1, 2, "Y", 2}}, shared));
        REQUIRE(shared == StepSequences{"AYT"});
    }
    SECTION("two insertions at one offset are both made") {
        StepSequences inserted{"ACGT"};
        REQUIRE(AlleleSplicer::apply_step_edits({{0, 2, 0, "TT", 1}, {0, 2, 0, "GG", 2}},
                                                inserted));
        REQUIRE(inserted[0].size() == 8);
    }
}

}
}
