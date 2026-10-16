/// \file allele_sequences.cpp
///
/// Unit tests for alleles that their walks do not spell: an allele given as a walk plus the
/// sequence of each of its steps. It is spelled by those sequences, and no panel haplotype carries
/// it, even where its walk is another allele's.

#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>
#include <gbwt/dynamic_gbwt.h>

#include "catch.hpp"
#include "../panel_lookup.hpp"
#include "../site_values.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A SNP between nodes 1 and 4: node 2 (T) or node 3 (G).
static void make_snp_graph(bdsg::HashGraph& graph) {
    handle_t h1 = graph.create_handle("AAAA", 1);
    handle_t h2 = graph.create_handle("T", 2);
    handle_t h3 = graph.create_handle("G", 3);
    handle_t h4 = graph.create_handle("CCCC", 4);
    graph.create_edge(h1, h2);
    graph.create_edge(h1, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h3, h4);
}

/// A walk over the given node IDs, each read forward.
static Traversal forward_walk(const HandleGraph& graph, const vector<nid_t>& ids) {
    Traversal walk;
    for (nid_t id : ids) {
        walk.push_back(graph.get_handle(id));
    }
    return walk;
}

TEST_CASE("An allele is spelled by its step sequences where it has them, and by its walk otherwise",
          "[allele_sequences]") {
    bdsg::HashGraph graph;
    make_snp_graph(graph);
    const Traversal ref = forward_walk(graph, {1, 2, 4});
    // The SNP's G on the reference's walk, rather than through node 3.
    const AlleleSequences sequences{{}, {"AAAA", "G", "CCCC"}};

    SECTION("an allele with no step sequences is spelled by its walk") {
        REQUIRE(own_sequences(&sequences, 0) == nullptr);
        REQUIRE(allele_sequence(graph, ref, own_sequences(&sequences, 0)) == "AAAATCCCC");
    }
    SECTION("an allele with step sequences is spelled by them, though its walk is the reference's") {
        REQUIRE(own_sequences(&sequences, 1) != nullptr);
        REQUIRE(allele_sequence(graph, ref, own_sequences(&sequences, 1)) == "AAAAGCCCC");
        // As the graph's own SNP allele is.
        REQUIRE(allele_sequence(graph, forward_walk(graph, {1, 3, 4}), nullptr) == "AAAAGCCCC");
    }
    SECTION("a walk read backward is spelled backward") {
        const Traversal backward{graph.get_handle(4, true), graph.get_handle(2, true),
                                 graph.get_handle(1, true)};
        REQUIRE(allele_sequence(graph, backward, nullptr) == "GGGGATTTT");
    }
    SECTION("no sequences, or an index past them, means the walk spells the allele") {
        REQUIRE(own_sequences(nullptr, 1) == nullptr);
        REQUIRE(own_sequences(&sequences, 2) == nullptr);
        const AlleleSequences none;
        REQUIRE(own_sequences(&none, 0) == nullptr);
    }
}

TEST_CASE("No panel haplotype carries an allele with step sequences, whatever its walk",
          "[allele_sequences]") {
    bdsg::HashGraph graph;
    make_snp_graph(graph);
    // Haplotype 0 takes node 2 and haplotype 1 node 3, each stored in both orientations, so that
    // GBWT sequences 0 and 1 are haplotype 0 and sequences 2 and 3 haplotype 1.
    gbwt::vector_type text;
    for (const vector<nid_t>& path : vector<vector<nid_t>>{{1, 2, 4}, {1, 3, 4}}) {
        for (nid_t id : path) {
            text.push_back(gbwt::Node::encode(id, false));
        }
        text.push_back(gbwt::ENDMARKER);
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            text.push_back(gbwt::Node::encode(*it, true));
        }
        text.push_back(gbwt::ENDMARKER);
    }
    gbwt::DynamicGBWT dynamic_gbwt;
    dynamic_gbwt.insert(text, true);
    const gbwt::GBWT index(dynamic_gbwt);
    const vector<size_t> haplotype_of_sequence{0, 0, 1, 1};
    const PanelLookup lookup(&graph, &index, &haplotype_of_sequence, 2);

    // The reference, the graph's SNP allele, and the SNP's G again on the reference's walk.
    const vector<Traversal> travs{forward_walk(graph, {1, 2, 4}), forward_walk(graph, {1, 3, 4}),
                                  forward_walk(graph, {1, 2, 4})};

    SECTION("without step sequences, a haplotype carries whichever allele follows its walk") {
        // Haplotype 0 follows the walk of alleles 0 and 2, and the later one wins.
        REQUIRE(lookup.alleles(travs) == vector<int>({2, 1}));
    }
    SECTION("with step sequences, the allele its walk does not spell is carried by none") {
        const AlleleSequences sequences{{}, {}, {"AAAA", "G", "CCCC"}};
        REQUIRE(lookup.alleles(travs, &sequences) == vector<int>({0, 1}));
    }
}

TEST_CASE("An edit site changes one base inside a node, read either way along the node",
          "[allele_sequences][off_panel]") {
    bdsg::HashGraph graph;
    make_snp_graph(graph);
    const handle_t node = graph.get_handle(1);   // AAAA

    SECTION("Read forward, the ALT replaces the base at the forward offset") {
        EditSite site = snv_edit_site(graph, node, 1, 'G');
        REQUIRE(site.bounds.inside_node);
        REQUIRE(site.bounds.offset == 1);
        REQUIRE(site.travs.size() == 2);
        REQUIRE(own_sequences(&site.sequences, 0) == nullptr);
        REQUIRE(allele_sequence(graph, site.travs[0], own_sequences(&site.sequences, 0)) == "AAAA");
        REQUIRE(allele_sequence(graph, site.travs[1], own_sequences(&site.sequences, 1)) == "AGAA");
        REQUIRE(bases_before_site(graph, site.bounds) == 1);
    }

    SECTION("Read backward, the ALT is complemented and the offset counted from the other end") {
        EditSite site = snv_edit_site(graph, graph.flip(node), 1, 'G');
        REQUIRE(allele_sequence(graph, site.travs[0], own_sequences(&site.sequences, 0)) == "TTTT");
        REQUIRE(allele_sequence(graph, site.travs[1], own_sequences(&site.sequences, 1)) == "TTCT");
        REQUIRE(bases_before_site(graph, site.bounds) == 2);
    }

    SECTION("Its ID has no underscore and names the node, offset and ALT") {
        REQUIRE(edit_site_id(12, 3, 'S', "T") == "op12.3.ST");
    }
}

}
}
