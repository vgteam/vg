/// \file sequence_allele.cpp
///
/// Unit tests for an allele of a snarl that is given as a sequence between the snarl's bounds:
/// how it is spelled, how reads are scored at a site that has one, and what it means for the
/// sites nested under the snarl.

#include <cstdio>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "../allele_likelihood.hpp"
#include "../../src/alignment_scorer.hpp"
#include "../child_crossing_finder.hpp"
#include "../linkage_model.hpp"
#include "../../src/unittest/catch.hpp"
#include "../site_read_source.hpp"
#include "../site_values.hpp"
#include "support/snarl_fixtures.hpp"
#include "support/staged_nested_sites.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// A snarl from node 1 to node 4 with two branches of 20 bases, nodes 2 and 3, that differ at
/// two bases.
struct TwoBranchSite {
    bdsg::HashGraph graph;
    Snarl snarl;
    vector<SnarlTraversal> traversals;

    TwoBranchSite() {
        const handle_t h1 = graph.create_handle("ACGTTGCAAGGCTTAACCGG", 1);
        const handle_t h2 = graph.create_handle("GATTACAGATTACAGGCTCA", 2);
        const handle_t h3 = graph.create_handle("GATTACAGTTTACAGCCTCA", 3);
        const handle_t h4 = graph.create_handle("TTGACCGTAGGCATCGATCA", 4);
        graph.create_edge(h1, h2);
        graph.create_edge(h2, h4);
        graph.create_edge(h1, h3);
        graph.create_edge(h3, h4);
        snarl.mutable_start()->set_node_id(1);
        snarl.mutable_end()->set_node_id(4);
        for (const vector<nid_t>& path : vector<vector<nid_t>>{{1, 2, 4}, {1, 3, 4}}) {
            SnarlTraversal& traversal = traversals.emplace_back();
            for (nid_t id : path) {
                traversal.add_visit()->set_node_id(id);
            }
        }
    }

    /// Return the alleles of the two branches, the first being the reference.
    CandidateAlleles alleles() const {
        return CandidateAlleles{walks_of(graph, traversals), 0};
    }
};

/// One edit of a read's alignment to a node: `from` bases of the node become the read's `bases`,
/// or, where `bases` is empty and `to` is `from`, are matched.
struct ReadEdit {
    size_t from;
    size_t to;
    string bases;
};

/// Return the read that visits each node forward with the given edits.
Alignment read_of(const HandleGraph& graph, const string& name,
                  const vector<pair<nid_t, vector<ReadEdit>>>& steps) {
    Alignment aln;
    aln.set_name(name);
    string sequence;
    for (const auto& [id, edits] : steps) {
        const string node = graph.get_sequence(graph.get_handle(id));
        Mapping* mapping = aln.mutable_path()->add_mapping();
        mapping->mutable_position()->set_node_id(id);
        size_t at = 0;
        for (const ReadEdit& edit : edits) {
            Edit* out = mapping->add_edit();
            out->set_from_length((int32_t)edit.from);
            out->set_to_length((int32_t)edit.to);
            if (edit.bases.empty() && edit.from == edit.to) {
                sequence += node.substr(at, edit.from);
            } else {
                out->set_sequence(edit.bases);
                sequence += edit.bases;
            }
            at += edit.from;
        }
    }
    aln.set_sequence(sequence);
    aln.set_mapping_quality(60);
    return aln;
}

/// Return the read turned to the other strand.
Alignment turned(const HandleGraph& graph, const Alignment& read) {
    return reverse_complement_alignment(read, [&](nid_t id) {
        return (int64_t)graph.get_length(graph.get_handle(id));
    });
}

/// Score `reads` at the snarl of `site` with `alleles`. Row r of the result is `reads[r]`.
AlleleReadLikelihoods score(const TwoBranchSite& site, const CandidateAlleles& alleles,
                            const vector<Alignment>& reads,
                            AlleleLikelihoodParams params = AlleleLikelihoodParams()) {
    InMemorySiteReadSource source;
    for (size_t i = 0; i < reads.size(); ++i) {
        Alignment named = reads[i];
        char prefix[16];
        snprintf(prefix, sizeof(prefix), "%06zu:", i);
        named.set_name(prefix + reads[i].name());
        source.add(named);
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, source, qual_scorer,
                                                      plain_scorer, params);
    return calculator.compute(materialized_of(site.graph, site.snarl, alleles), {}, 2);
}

/// The 60 bases the sequence allele of the tests puts in the middle of node 2's.
const string NOVEL = "CCATGGTACGTTAGCAGTCAATCGGATCCTTAGGCAATCGTACCGATTGACCTGAAGTCC";

/// Return the reads of the test of a site with no sequence allele: one of each branch, one with
/// three more bases and a changed base, one that stops inside the snarl, and one on the other
/// strand.
vector<Alignment> walk_reads(const TwoBranchSite& site) {
    vector<Alignment> reads;
    reads.push_back(read_of(site.graph, "first", {{1, {{20, 20, ""}}}, {2, {{20, 20, ""}}},
                                                  {4, {{20, 20, ""}}}}));
    reads.push_back(read_of(site.graph, "second", {{1, {{20, 20, ""}}}, {3, {{20, 20, ""}}},
                                                   {4, {{20, 20, ""}}}}));
    reads.push_back(read_of(site.graph, "changed",
                            {{1, {{20, 20, ""}}},
                             {2, {{4, 4, ""}, {1, 1, "C"}, {5, 5, ""}, {0, 3, "TTG"}, {10, 10, ""}}},
                             {4, {{20, 20, ""}}}}));
    reads.push_back(read_of(site.graph, "partial", {{1, {{20, 20, ""}}}, {3, {{12, 12, ""}}}}));
    reads.push_back(turned(site.graph, reads[1]));
    return reads;
}

}

TEST_CASE("A site with no sequence allele is scored as it was", "[sequence_allele]") {
    const TwoBranchSite site;
    const vector<Alignment> reads = walk_reads(site);
    // The best ln-likelihood of each read, and its relative likelihood under each allele, from
    // the calculator before it scored sequence alleles.
    const vector<vector<double>> expected{
        {82.999516124318689, 1, 9.8241433544953682e-07},
        {82.999516124318689, 9.8241433544953682e-07, 1},
        {65.916287630716312, 1, 9.2326444237392154e-28},
        {44.266408599636634, 1.5350779870416275e-11, 1},
        {82.999516124318689, 9.8241433544953682e-07, 1}};
    for (bool optimal : {false, true}) {
        AlleleLikelihoodParams params;
        params.optimal_pairing = optimal;
        params.insertion_gap_nats = 0.9;
        const AlleleReadLikelihoods matrix = score(site, site.alleles(), reads, params);
        REQUIRE(matrix.num_reads() == reads.size());
        for (size_t r = 0; r < reads.size(); ++r) {
            REQUIRE(matrix.best_ln_likelihood(r) == Approx(expected[r][0]).epsilon(1e-13));
            REQUIRE(matrix.rel(r, 0) == Approx(expected[r][1]).epsilon(1e-13));
            REQUIRE(matrix.rel(r, 1) == Approx(expected[r][2]).epsilon(1e-13));
        }
        REQUIRE(matrix.traversal_length(0) == 20);
        REQUIRE(matrix.traversal_length(1) == 20);
    }
}

TEST_CASE("A sequence allele spells its sequence between the bounds of its snarl",
          "[sequence_allele]") {
    const TwoBranchSite site;
    CandidateAlleles alleles = site.alleles();
    REQUIRE(add_sequence_allele(site.graph, alleles, NOVEL) == 2);
    REQUIRE(alleles.walks.size() == 3);
    REQUIRE(alleles.walks[2] == Traversal{site.graph.get_handle(1), site.graph.get_handle(4)});
    REQUIRE(alleles.ref_allele == 0);
    REQUIRE(alleles.spelled_by_walk(0));
    REQUIRE(alleles.spelled_by_walk(1));
    REQUIRE(!alleles.spelled_by_walk(2));
    REQUIRE(alleles.sequence(site.graph, 2)
            == "ACGTTGCAAGGCTTAACCGG" + NOVEL + "TTGACCGTAGGCATCGATCA");
    REQUIRE(alleles.sequence(site.graph, 1)
            == "ACGTTGCAAGGCTTAACCGGGATTACAGTTTACAGCCTCATTGACCGTAGGCATCGATCA");
    // An allele that deletes everything between the bounds is still not spelled by its walk.
    REQUIRE(add_sequence_allele(site.graph, alleles, "") == 3);
    REQUIRE(!alleles.spelled_by_walk(3));
    REQUIRE(alleles.sequence(site.graph, 3) == "ACGTTGCAAGGCTTAACCGGTTGACCGTAGGCATCGATCA");

    // A snarl read the other way takes its bounds the way its walks read them.
    CandidateAlleles turned_alleles;
    turned_alleles.walks.push_back({site.graph.get_handle(4, true), site.graph.get_handle(2, true),
                                    site.graph.get_handle(1, true)});
    turned_alleles.ref_allele = 0;
    REQUIRE(add_sequence_allele(site.graph, turned_alleles, "AC") == 1);
    REQUIRE(turned_alleles.sequence(site.graph, 1)
            == "TGATCGATGCCTACGGTCAAACCCGGTTAAGCCTTGCAACGT");

    CandidateAlleles none;
    REQUIRE(!add_sequence_allele(site.graph, none, NOVEL));
    REQUIRE(none.walks.empty());
}

TEST_CASE("A walk allele spells its sequence with changes to the steps of the walk",
          "[sequence_allele]") {
    bdsg::HashGraph graph;
    const Traversal walk{graph.create_handle("ACGTACGT", 1), graph.create_handle("GGCC", 2),
                         graph.create_handle("TTAA", 3), graph.create_handle("CAGTCAGT", 4)};
    for (size_t k = 0; k + 1 < walk.size(); ++k) {
        graph.create_edge(walk[k], walk[k + 1]);
    }
    CandidateAlleles alleles{{walk}, 0};
    const StringAlignerParams aligner{.wavefront = true};
    for (const string interior : {"GGCATTAA", "TAA", "GGCCTTAAGG", "TTGGCCTTAA", ""}) {
        const optional<size_t> allele = add_walk_allele(graph, alleles, walk, interior, aligner);
        REQUIRE(allele == alleles.walks.size() - 1);
        REQUIRE(alleles.walks[*allele] == walk);
        REQUIRE(!alleles.spelled_by_walk(*allele));
        REQUIRE(alleles.sequence(graph, *allele) == "ACGTACGT" + interior + "CAGTCAGT");
        // The bound nodes keep their bases.
        const StepSequences steps = alleles.step_sequences(graph, *allele);
        REQUIRE(steps.front() == "ACGTACGT");
        REQUIRE(steps.back() == "CAGTCAGT");
    }
    // A run of mismatches changes the step that holds it.
    REQUIRE(alleles.edits[1] == vector<AlleleEdit>{{1, 3, 1, "A"}});
    // A deletion across two nodes changes both, and the node it takes whole stays in the walk.
    REQUIRE(alleles.edits[2] == vector<AlleleEdit>{{1, 0, 4, ""}, {2, 0, 1, ""}});
    REQUIRE(alleles.step_sequences(graph, 2)[1].empty());
    // An insertion goes after the base before it, or before the first base between the bounds.
    REQUIRE(alleles.edits[3] == vector<AlleleEdit>{{2, 4, 0, "GG"}});
    REQUIRE(alleles.edits[4] == vector<AlleleEdit>{{1, 0, 0, "TT"}});

    REQUIRE(!add_walk_allele(graph, alleles, {walk.front()}, "A", aligner));
}

TEST_CASE("Reads that span a snarl with a sequence allele are scored by aligning them to every allele",
          "[sequence_allele]") {
    const TwoBranchSite site;
    const string branch = site.graph.get_sequence(site.graph.get_handle(2));
    CandidateAlleles alleles = site.alleles();
    REQUIRE(add_sequence_allele(site.graph, alleles,
                                branch.substr(0, 10) + NOVEL + branch.substr(10)) == 2);

    // The novel bases as a read has them, with one base changed and one lost.
    string noisy = NOVEL;
    noisy[20] = noisy[20] == 'A' ? 'C' : 'A';
    noisy.erase(41, 1);

    vector<Alignment> reads = walk_reads(site);
    const size_t first = 0;
    const size_t second = 1;
    const size_t partial = 3;
    // A read with the novel bases, aligned as an insertion in node 2.
    const size_t inserted = reads.size();
    reads.push_back(read_of(site.graph, "inserted",
                            {{1, {{20, 20, ""}}},
                             {2, {{10, 10, ""}, {0, noisy.size(), noisy}, {10, 10, ""}}},
                             {4, {{20, 20, ""}}}}));
    // The same read on the other strand.
    const size_t reversed = reads.size();
    reads.push_back(turned(site.graph, reads[inserted]));
    // A read with the novel bases that the aligner put through node 3, with the insertion at
    // the end of node 1's mapping and the first ten bases of the branch as an insertion too.
    const size_t rerouted = reads.size();
    reads.push_back(read_of(site.graph, "rerouted",
                            {{1, {{20, 20, ""}, {0, 10 + NOVEL.size(), branch.substr(0, 10) + NOVEL}}},
                             {3, {{10, 0, ""}, {5, 5, ""}, {1, 1, "G"}, {4, 4, ""}}},
                             {4, {{20, 20, ""}}}}));
    // A read with five more bases than the sequence allele has.
    const size_t longer = reads.size();
    reads.push_back(read_of(site.graph, "longer",
                            {{1, {{20, 20, ""}}},
                             {2, {{10, 10, ""}, {0, NOVEL.size() + 5, NOVEL + "ACGTA"}, {10, 10, ""}}},
                             {4, {{20, 20, ""}}}}));

    for (bool optimal : {false, true}) {
        AlleleLikelihoodParams params;
        params.optimal_pairing = optimal;
        params.insertion_gap_nats = 0.9;
        const AlleleReadLikelihoods matrix = score(site, alleles, reads, params);
        const AlleleReadLikelihoods walks_only = score(site, site.alleles(), reads, params);
        REQUIRE(matrix.num_reads() == reads.size());
        REQUIRE(matrix.num_alleles() == 3);

        SECTION(optimal ? "Optimal pairing: each read fits its own allele best"
                        : "Greedy pairing: each read fits its own allele best") {
            for (size_t r : {inserted, reversed, rerouted, longer}) {
                REQUIRE(matrix.rel(r, 2) == Approx(1.0));
                REQUIRE(matrix.rel(r, 0) < 1e-6);
                REQUIRE(matrix.rel(r, 1) < 1e-6);
            }
            REQUIRE(matrix.rel(first, 0) == Approx(1.0));
            REQUIRE(matrix.rel(first, 2) < 1e-6);
            REQUIRE(matrix.rel(second, 1) == Approx(1.0));
            REQUIRE(matrix.rel(second, 2) < 1e-6);
            // The two strands of a read are scored alike.
            for (size_t a = 0; a < 3; ++a) {
                REQUIRE(matrix.rel(reversed, a) == Approx(matrix.rel(inserted, a)));
            }
            REQUIRE(matrix.best_ln_likelihood(reversed)
                    == Approx(matrix.best_ln_likelihood(inserted)));
        }
        SECTION(optimal ? "Optimal pairing: an aligned read is scored in the units of a paired one"
                        : "Greedy pairing: an aligned read is scored in the units of a paired one") {
            // A read that matches a walk allele scores on it what its visits score, and loses to
            // the other walk allele what two mismatches cost.
            for (size_t r : {first, second}) {
                REQUIRE(matrix.best_ln_likelihood(r)
                        == Approx(walks_only.best_ln_likelihood(r)));
                REQUIRE(matrix.rel(r, 0) == Approx(walks_only.rel(r, 0)));
                REQUIRE(matrix.rel(r, 1) == Approx(walks_only.rel(r, 1)));
            }
            // The read with the novel bases loses to the sequence allele only its two errors:
            // a mismatch, and a gap of one base in the allele that the read lacks.
            const MatrixAlignmentScorer scorer;
            const double log_base = scorer.get_log_base();
            const double perfect = log_base * scorer.match * (40.0 + 80.0);
            const double errors = log_base * (2.0 * scorer.match + scorer.mismatch
                                              + scorer.gap_open);
            REQUIRE(matrix.best_ln_likelihood(inserted) == Approx(perfect - errors));
        }
        SECTION(optimal ? "Optimal pairing: a read that does not span says nothing of the allele"
                        : "Greedy pairing: a read that does not span says nothing of the allele") {
            REQUIRE(matrix.best_ln_likelihood(partial)
                    == Approx(walks_only.best_ln_likelihood(partial)));
            REQUIRE(matrix.rel(partial, 0) == Approx(walks_only.rel(partial, 0)));
            REQUIRE(matrix.rel(partial, 1) == Approx(walks_only.rel(partial, 1)));
            REQUIRE(matrix.rel(partial, 2) == Approx(1.0));
        }
        SECTION(optimal ? "Optimal pairing: the depth length is the allele's own"
                        : "Greedy pairing: the depth length is the allele's own") {
            REQUIRE(matrix.traversal_length(0) == 20);
            REQUIRE(matrix.traversal_length(1) == 20);
            REQUIRE(matrix.traversal_length(2) == 80);
        }
        SECTION(optimal ? "Optimal pairing: past a limit the sequence allele is not offered"
                        : "Greedy pairing: past a limit the sequence allele is not offered") {
            AlleleLikelihoodParams few_bases = params;
            AlleleLikelihoodParams few_alignments = params;
            // Eight reads span the snarl, and each is aligned to three alleles.
            few_alignments.sequence_allele_max_alignments = 23;
            const AlleleReadLikelihoods capped = score(site, alleles, reads, few_alignments);
            REQUIRE(capped.num_reads() == reads.size());
            for (size_t r = 0; r < reads.size(); ++r) {
                REQUIRE(capped.rel(r, 2) == 0.0);
                REQUIRE(capped.rel(r, 0) == Approx(walks_only.rel(r, 0)));
                REQUIRE(capped.rel(r, 1) == Approx(walks_only.rel(r, 1)));
                REQUIRE(capped.best_ln_likelihood(r)
                        == Approx(walks_only.best_ln_likelihood(r)));
            }
            few_alignments.sequence_allele_max_alignments = 24;
            const AlleleReadLikelihoods within = score(site, alleles, reads, few_alignments);
            REQUIRE(within.rel(inserted, 2) == Approx(1.0));
            REQUIRE(within.rel(inserted, 0) < 1e-6);
            // A read with more bases between the bounds than the limit is scored as one that
            // does not span.
            few_bases.sequence_allele_max_bases = 80;
            const AlleleReadLikelihoods mixed = score(site, alleles, reads, few_bases);
            REQUIRE(mixed.rel(first, 2) < 1e-6);
            REQUIRE(mixed.rel(inserted, 2) == Approx(1.0));
            REQUIRE(mixed.rel(inserted, 0) < 1e-6);
            REQUIRE(mixed.rel(longer, 2) == Approx(1.0));
            REQUIRE(mixed.rel(longer, 0) == Approx(walks_only.rel(longer, 0)));
            REQUIRE(mixed.rel(longer, 1) == Approx(walks_only.rel(longer, 1)));
            REQUIRE(mixed.best_ln_likelihood(longer)
                    == Approx(walks_only.best_ln_likelihood(longer)));
        }
        SECTION(optimal ? "Optimal pairing: an allele too long to align a read to leaves the "
                          "others as they are"
                        : "Greedy pairing: an allele too long to align a read to leaves the "
                          "others as they are") {
            // An allele with 9,000 bases between the bounds, where the reads have under 100.
            string long_bases;
            for (size_t i = 0; i < 9000; ++i) {
                long_bases += "ACGT"[(i * i + i / 7) % 4];
            }
            CandidateAlleles with_long = alleles;
            REQUIRE(add_sequence_allele(site.graph, with_long, long_bases) == 3);
            const AlleleReadLikelihoods scored = score(site, with_long, reads, params);
            REQUIRE(scored.num_reads() == reads.size());
            for (size_t r = 0; r < reads.size(); ++r) {
                for (size_t a = 0; a < 3; ++a) {
                    REQUIRE(scored.rel(r, a) == Approx(matrix.rel(r, a)));
                }
                REQUIRE(scored.best_ln_likelihood(r) == Approx(matrix.best_ln_likelihood(r)));
            }
            // A read that spans the snarl does not fit the long allele.
            for (size_t r : {inserted, reversed, rerouted, longer}) {
                REQUIRE(scored.rel(r, 3) < 1e-6);
            }
        }
    }
}

TEST_CASE("A strand with a sequence allele crosses no site nested in its snarl",
          "[sequence_allele]") {
    using namespace vg::multipass;
    StagedNestedSites fixture;
    const string novel = "ACGTACGTAC";
    // The top site's alleles are the reference walk, the walk from 1 straight to 8, and the walk
    // through node 6. The sequence allele is the fourth.
    const size_t sequence_allele = 3;

    SECTION("No child is crossed by the sequence allele's walk") {
        CandidateAlleles alleles;
        alleles.walks = {{fixture.nested.handle(1), fixture.nested.handle(2),
                          fixture.nested.handle(3), fixture.nested.handle(4),
                          fixture.nested.handle(5), fixture.nested.handle(7),
                          fixture.nested.handle(8)}};
        alleles.ref_allele = 0;
        REQUIRE(add_sequence_allele(fixture.nested.graph, alleles, novel) == 1);
        vector<ChildCrossingFinder> walks;
        for (const Traversal& walk : alleles.walks) {
            walks.emplace_back(fixture.overlay, fixture.nested.graph, walk);
        }
        for (const pair<nid_t, nid_t>& bounds : vector<pair<nid_t, nid_t>>{{2, 7}, {3, 5}}) {
            const net_handle_t child = find_snarl(fixture.overlay, fixture.nested.graph,
                                                  bounds.first, bounds.second);
            REQUIRE(walks[1].crossings(child).empty());
            bool known = false;
            const uint64_t mask = crossing_mask(walks, child, &known);
            REQUIRE(known);
            REQUIRE(mask == 0b01);
            // Genotype linking gives the child one copy, on the reference's strand, where the
            // other strand has the sequence allele, and none where both have it.
            const LinkageCollector::Relation het = LinkageCollector::relate_to_parent(mask, 1, 0);
            REQUIRE(het.known);
            REQUIRE(het.copies == 1);
            REQUIRE(het.carrying_trav == 0);
            REQUIRE(LinkageCollector::relate_to_parent(mask, 1, 1).copies == 0);
        }
        // The same holds for the site inside node 4.
        fixture.index();
        REQUIRE(crossing_mask(walks, fixture.site(fixture.edit_key).net) == 0b01);
    }
    SECTION("The splicer spells the allele as it is, and the other strand through its children") {
        // Strand 0 of the top site has the sequence allele. Strand 1 has the walk through node
        // 6, where the middle site, on that strand alone, has its reference.
        fixture.phase(fixture.top_key, (int)sequence_allele, 2);
        fixture.phase(fixture.middle_key, 0, -1, 1, 1);
        fixture.index();
        StagedSite& top = fixture.site(fixture.top_key);
        REQUIRE(add_sequence_allele(fixture.nested.graph, top.alleles, novel) == sequence_allele);
        const AlleleSplicer& splicer = fixture.splicer;

        REQUIRE(splicer.strand_traversal(top, 0) == (int)sequence_allele);
        REQUIRE(splicer.strand_sequence(top, 0) == "GCA" + novel + "CTGA");
        REQUIRE(splicer.strand_sequence(top, 1) == "GCATGCTGAGCAGCTGA");
        REQUIRE(!splicer.strand_sequence(fixture.site(fixture.middle_key), 0));
        REQUIRE(splicer.strand_sequence(fixture.site(fixture.middle_key), 1) == "TGCTGAGCAG");

        // The record's alleles keep the sequence allele's spelling beside the spliced one.
        AlleleSplicer::SplicedAlleles spliced;
        REQUIRE(splicer.splice(top, spliced));
        REQUIRE(spliced.travs[sequence_allele] == top.alleles.walks[sequence_allele]);
        REQUIRE(allele_sequence(fixture.nested.graph, spliced.travs[sequence_allele],
                                own_sequences(&spliced.sequences, sequence_allele))
                == "GCA" + novel + "CTGA");
        REQUIRE(allele_sequence(fixture.nested.graph, spliced.travs[2],
                                own_sequences(&spliced.sequences, 2))
                == "GCATGCTGAGCAGCTGA");
    }
}

}
}
