/// \file unittest/allele_likelihood_scoring.cpp
///
/// Unit tests for scoring reads against alleles from their existing graph
/// alignment, on hand-built graphs with hand-built alignments.
///
/// These exist because the model tests in allele_likelihood.cpp deliberately use
/// hand-built matrices, so they cannot see anything wrong with the scoring that
/// *produces* those matrices. Two bugs got through exactly that gap: reverse
/// strand reads failed to anchor and were scored against the wrong allele, and
/// reads traversing a deletion edge were discarded as uninformative. Both are
/// pinned below.
///

#include <vector>

#include <bdsg/hash_graph.hpp>

#include "allele_likelihood.hpp"
#include "alignment_scorer.hpp"
#include "catch.hpp"
#include "site_read_source.hpp"
#include "snarls.hpp"
#include "utility.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A site with a SNP and a deletion of the SNP-bearing node:
///
///        2 (T)
///      /      \
///   1 --- 3 (G) --- 4          and the deletion edge 1 -> 4
///      \__________/
///
/// Node 1 and 4 are the snarl boundaries; 2 and 3 are the SNP alleles.
struct SnpAndDeletionSite {
    bdsg::HashGraph graph;
    Snarl snarl;
    vector<SnarlTraversal> traversals;   // ref (1,2,4), alt (1,3,4), deletion (1,4)
    unique_ptr<SnarlManager> manager;

    SnpAndDeletionSite() {
        handle_t h1 = graph.create_handle("AAAACCCC", 1);
        handle_t h2 = graph.create_handle("T", 2);
        handle_t h3 = graph.create_handle("G", 3);
        handle_t h4 = graph.create_handle("GGGGTTTT", 4);

        graph.create_edge(h1, h2);
        graph.create_edge(h2, h4);
        graph.create_edge(h1, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h1, h4);   // the deletion

        snarl.mutable_start()->set_node_id(1);
        snarl.mutable_end()->set_node_id(4);
        snarl.set_type(ULTRABUBBLE);
        snarl.set_start_end_reachable(true);

        vector<Snarl> snarls{snarl};
        manager.reset(new SnarlManager(snarls.begin(), snarls.end()));

        vector<vector<nid_t>> allele_paths{{1, 2, 4}, {1, 3, 4}, {1, 4}};
        traversals.resize(allele_paths.size());
        for (size_t i = 0; i < allele_paths.size(); ++i) {
            for (nid_t id : allele_paths[i]) {
                Visit* v = traversals[i].add_visit();
                v->set_node_id(id);
                v->set_backward(false);
            }
        }
    }
};

/// Build an all-match alignment over the given (node, is_reverse) steps.
/// The read sequence is the concatenation of the visited node sequences, so every
/// edit is a perfect match.
static Alignment make_matching_alignment(const HandleGraph& graph, const string& name,
                                         const vector<pair<nid_t, bool>>& steps) {
    Alignment aln;
    aln.set_name(name);
    string seq;
    for (auto& step : steps) {
        string node_seq = graph.get_sequence(graph.get_handle(step.first, step.second));
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(step.first);
        m->mutable_position()->set_is_reverse(step.second);
        m->mutable_position()->set_offset(0);
        Edit* e = m->add_edit();
        e->set_from_length(node_seq.size());
        e->set_to_length(node_seq.size());
        seq += node_seq;
    }
    aln.set_sequence(seq);
    aln.set_quality(string(seq.size(), (char)30));
    aln.set_mapping_quality(60);
    return aln;
}

/// Run the calculator over one site with the given reads.
static AlleleReadLikelihoods score_site(SnpAndDeletionSite& site, const vector<Alignment>& reads) {
    InMemorySiteReadSource source;
    for (const Alignment& aln : reads) {
        source.add(aln);
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, *site.manager, source, qual_scorer,
                                                      plain_scorer);
    return calculator.compute(site.snarl, site.traversals);
}

TEST_CASE("A reverse-strand read scores the same as its forward equivalent",
          "[allele_likelihood][scoring]") {
    // The bug this pins: alleles impose a reading direction on the site, and a read
    // aligned to the other strand visits the same nodes with the opposite
    // orientation flag. Anchoring on the raw flag meant reverse-strand reads
    // matched nothing, fell through to the substitution path, and were scored
    // against the wrong allele -- roughly half of all reads, at every site.
    SnpAndDeletionSite site;

    Alignment forward = make_matching_alignment(site.graph, "fwd", {{1, false}, {2, false}, {4, false}});
    // The same underlying fragment sequenced the other way round: the path runs
    // backwards through the site and every step is flipped.
    Alignment reverse = make_matching_alignment(site.graph, "rev", {{4, true}, {2, true}, {1, true}});

    AlleleReadLikelihoods matrix = score_site(site, {forward, reverse});
    REQUIRE(matrix.num_reads() == 2);

    SECTION("both reads prefer the reference allele they actually traverse") {
        for (size_t r = 0; r < matrix.num_reads(); ++r) {
            REQUIRE(matrix.rel(r, 0) == Approx(1.0));
            REQUIRE(matrix.rel(r, 1) < 1.0);
        }
    }

    SECTION("the two reads are scored identically") {
        for (size_t a = 0; a < matrix.num_alleles(); ++a) {
            REQUIRE(matrix.rel(0, a) == Approx(matrix.rel(1, a)));
        }
    }

    SECTION("so a matrix of one strand genotypes the same as a mix of both") {
        AlleleReadLikelihoods fwd_only = score_site(site, {forward, forward});
        REQUIRE(fwd_only.genotype_likelihood({0, 0}) == Approx(matrix.genotype_likelihood({0, 0})));
        REQUIRE(fwd_only.genotype_likelihood({1, 1}) == Approx(matrix.genotype_likelihood({1, 1})));
    }
}

TEST_CASE("A read spanning a deletion is kept and prefers the deletion allele",
          "[allele_likelihood][scoring]") {
    // The bug this pins: a read traversing straight from one boundary node to the
    // other touches no interior node, so a node-based "is it informative" test
    // discarded it. That read is the only direct evidence the deletion allele ever
    // gets, so dropping it silently destroyed deletion genotyping -- the caller
    // saw reference-supporting reads only, called hom-ref, and emitted no record.
    SnpAndDeletionSite site;

    Alignment deletion_read =
        make_matching_alignment(site.graph, "del", {{1, false}, {4, false}});

    AlleleReadLikelihoods matrix = score_site(site, {deletion_read});

    SECTION("it is not discarded") {
        REQUIRE(matrix.num_reads() == 1);
    }

    SECTION("it prefers the deletion allele over both spanning alleles") {
        REQUIRE(matrix.rel(0, 2) == Approx(1.0));
        REQUIRE(matrix.rel(0, 0) < 1.0);
        REQUIRE(matrix.rel(0, 1) < 1.0);
    }

    SECTION("and it makes the homozygous deletion the best genotype") {
        auto scored = matrix.score_genotypes(2);
        size_t best = 0;
        for (size_t i = 1; i < scored.size(); ++i) {
            if (scored[i].second > scored[best].second) {
                best = i;
            }
        }
        REQUIRE(scored[best].first == vector<int>({2, 2}));
    }
}

TEST_CASE("A read entirely inside one boundary node is still dropped",
          "[allele_likelihood][scoring]") {
    // The counterpart to the test above: widening informativeness to include
    // boundary-to-boundary reads must not accidentally admit reads that genuinely
    // cannot discriminate. A read sitting inside a single boundary node uses no
    // edge inside the site and touches no interior node, so every allele explains
    // it identically.
    SnpAndDeletionSite site;

    Alignment inside_boundary = make_matching_alignment(site.graph, "flank", {{1, false}});

    AlleleReadLikelihoods matrix = score_site(site, {inside_boundary});
    REQUIRE(matrix.num_reads() == 0);
}

TEST_CASE("A read over the SNP discriminates between the two SNP alleles",
          "[allele_likelihood][scoring]") {
    SnpAndDeletionSite site;

    Alignment ref_read = make_matching_alignment(site.graph, "ref", {{1, false}, {2, false}, {4, false}});
    Alignment alt_read = make_matching_alignment(site.graph, "alt", {{1, false}, {3, false}, {4, false}});

    AlleleReadLikelihoods matrix = score_site(site, {ref_read, alt_read});
    REQUIRE(matrix.num_reads() == 2);

    // Each read fits the allele it traverses best, and they disagree.
    REQUIRE(matrix.rel(0, 0) == Approx(1.0));
    REQUIRE(matrix.rel(0, 1) < 1.0);
    REQUIRE(matrix.rel(1, 1) == Approx(1.0));
    REQUIRE(matrix.rel(1, 0) < 1.0);

    // One read each way is the textbook heterozygote.
    auto scored = matrix.score_genotypes(2);
    size_t best = 0;
    for (size_t i = 1; i < scored.size(); ++i) {
        if (scored[i].second > scored[best].second) {
            best = i;
        }
    }
    REQUIRE(scored[best].first == vector<int>({0, 1}));
}

TEST_CASE("Every allele is scored over the same span of read bases",
          "[allele_likelihood][scoring]") {
    // The window invariant. A read placeable over more bases on one allele than
    // another must not gain from the length difference alone: bases an allele
    // cannot place are charged, never omitted. Here the deletion allele cannot
    // place the SNP base that the reference and alt alleles can.
    SnpAndDeletionSite site;

    Alignment ref_read = make_matching_alignment(site.graph, "ref", {{1, false}, {2, false}, {4, false}});
    AlleleReadLikelihoods matrix = score_site(site, {ref_read});
    REQUIRE(matrix.num_reads() == 1);

    // The deletion allele must be penalised relative to the allele the read
    // actually traverses, not rewarded for being shorter.
    REQUIRE(matrix.rel(0, 0) == Approx(1.0));
    REQUIRE(matrix.rel(0, 2) < 1.0);
}

}
}
