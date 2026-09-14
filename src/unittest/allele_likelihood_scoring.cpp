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

    /// `tail` is node 4's sequence. It is a parameter only so a test can vary the
    /// flank length while holding the variant fixed; every existing caller gets the
    /// original graph.
    SnpAndDeletionSite(const string& tail = "GGGGTTTT") {
        handle_t h1 = graph.create_handle("AAAACCCC", 1);
        handle_t h2 = graph.create_handle("T", 2);
        handle_t h3 = graph.create_handle("G", 3);
        handle_t h4 = graph.create_handle(tail, 4);

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

/// Run the calculator over one site with the given reads, at the given ploidy.
///
/// `realign` selects the walk: false is the greedy default, true the optimal one that
/// `--realign` and `--preset ont` turn on. The invariants below must hold for BOTH -- they are
/// properties of the scoring model, not of how the correspondence is searched for.
static AlleleReadLikelihoods score_site(SnpAndDeletionSite& site, const vector<Alignment>& reads,
                                        int ploidy = 2, bool realign = false) {
    InMemorySiteReadSource source;
    for (const Alignment& aln : reads) {
        source.add(aln);
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    AlleleLikelihoodParams params;
    params.realign = realign;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, *site.manager, source, qual_scorer,
                                                      plain_scorer, params);
    return calculator.compute(site.snarl, site.traversals, ploidy);
}

TEST_CASE("The depth rate is per haplotype, so it follows the site's ploidy",
          "[allele_likelihood][scoring]") {
    // The same reads over the same site, genotyped haploid and diploid. The window's
    // read density is a property of the data and does not change; the *per-haplotype*
    // rate does, by exactly the ploidy ratio, and lambda with it.
    //
    // This is a regression test. The ploidy was fixed at 2 inside the rate calculation
    // and the calculator was never told the site's own, so every haploid region -- chrY,
    // and chrX under --ploidy-regex -- got a lambda wrong by a factor of two while the
    // observed read count was right. It survived because the tier-2 evaluation is
    // autosomes only, and because the two components either side of the seam were both
    // tested on their own: the depth term is exercised by passing a rate in directly,
    // and haploid genotyping is exercised without the depth term.
    SnpAndDeletionSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 12; ++i) {
        reads.push_back(make_matching_alignment(site.graph, "r" + std::to_string(i),
                                  {{(nid_t)1, false}, {(nid_t)2, false}, {(nid_t)4, false}}));
    }

    AlleleReadLikelihoods diploid = score_site(site, reads, 2);
    AlleleReadLikelihoods haploid = score_site(site, reads, 1);

    if (diploid.uses_depth_term()) {
        vector<int> hom{0, 0};
        double lambda_diploid = diploid.expected_reads(hom);
        double lambda_haploid = haploid.expected_reads(hom);
        REQUIRE(lambda_diploid > 0.0);
        // Halving the ploidy doubles the per-haplotype rate, and lambda with it.
        REQUIRE(lambda_haploid == Approx(2.0 * lambda_diploid));
    }
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

TEST_CASE("A one-base indel costs the same whichever side carries it",
          "[allele_likelihood][scoring]") {
    // The bug this pins: the walk in score_read_against_allele searches *ahead in the
    // allele* for each read node, but when it finds none it assumes the read node
    // SUBSTITUTES for the allele's current node and consumes that node. When the read
    // node is instead a pure insertion -- the allele simply lacks it -- consuming the
    // allele node desynchronises the walk. The allele node the read would have matched
    // is burned against the wrong read node, and the read's real visit to it then falls
    // through to the allele-exhausted branch and is charged a second time.
    //
    // Concretely, the read (1,2,4) against the deletion allele (1,4) was scored as
    // mismatch(T vs G) + gap(7) + gap(8) -- node 4's whole length charged twice, once as
    // a length difference against node 2 and once as an unplaceable read node -- instead
    // of the single one-base gap the event actually is. On real ONT data the same 1 bp
    // event cost 1 score unit in one direction and 67-68 in the other, and 22% of read
    // rows have a best allele carrying a gap.
    //
    // The invariant asserted here is direction symmetry, which needs no knowledge of
    // gap_open or the log base: one inserted base and one deleted base are the same
    // event seen from the two sides, so they must carry the same penalty.
    SnpAndDeletionSite site;

    // Read takes the SNP node; the deletion allele (index 2) lacks it -> 1 bp insertion.
    Alignment spanning = make_matching_alignment(site.graph, "spanning",
                                                 {{1, false}, {2, false}, {4, false}});
    // Read skips it; the reference allele (index 0) carries it -> 1 bp deletion.
    Alignment deleting = make_matching_alignment(site.graph, "deleting", {{1, false}, {4, false}});

    AlleleReadLikelihoods matrix = score_site(site, {spanning, deleting});
    REQUIRE(matrix.num_reads() == 2);

    SECTION("each read matches its own allele exactly") {
        REQUIRE(matrix.rel(0, 0) == Approx(1.0));   // spanning read vs reference
        REQUIRE(matrix.rel(1, 2) == Approx(1.0));   // deleting read vs deletion
    }

    SECTION("the insertion costs a few score units, not tens") {
        // A bare `> 0.0` would NOT gate this: before the fix the value was a denormal
        // around 1e-23, which prints as 0.0 but is not zero. 1e-10 is about 17 score
        // units at the model's 1.3833 nats per unit -- far above the handful a one-base
        // event can justify, and far below the ~38 the bug charged.
        REQUIRE(matrix.rel(0, 2) > 1e-10);
    }

    SECTION("and costs within one match unit of the deletion, not the flank's length") {
        // The two are not exactly equal: an inserted base exists in the read and could
        // have been matched, so it forgoes `match` credit that a deleted base never had.
        // That residual is one score unit and belongs to the score model, not the walk.
        // What the walk must not do is charge the FLANK, which is what the bug did.
        REQUIRE(matrix.rel(0, 2) < matrix.rel(1, 0));
        REQUIRE(matrix.rel(0, 2) > 0.1 * matrix.rel(1, 0));
    }

    SECTION("and does not grow with the length of the flanking node") {
        // The sharpest statement of the bug, and parameter-free. It charged node 4's
        // whole length twice -- once as a length difference against node 2, once as an
        // unplaceable read node -- so the cost of a ONE BASE insertion scaled with the
        // flank. A walk that charges the event itself cannot care how long the flank is.
        SnpAndDeletionSite long_site("GGGGTTTT" + string(32, 'A'));
        Alignment long_spanning = make_matching_alignment(
            long_site.graph, "spanning", {{1, false}, {2, false}, {4, false}});
        Alignment long_deleting = make_matching_alignment(
            long_site.graph, "deleting", {{1, false}, {4, false}});
        AlleleReadLikelihoods long_matrix = score_site(long_site, {long_spanning, long_deleting});

        REQUIRE(long_matrix.rel(0, 2) == Approx(matrix.rel(0, 2)));
        REQUIRE(long_matrix.rel(1, 0) == Approx(matrix.rel(1, 0)));
    }

    SECTION("equal-length substituted nodes are still charged as a substitution") {
        // Regression guard for the fix itself: the discriminator must not divert the
        // genuine substitution case, where the allele's node really is the read node's
        // counterpart and consuming it is right.
        Alignment other_snp = make_matching_alignment(site.graph, "othersnp",
                                                      {{1, false}, {3, false}, {4, false}});
        AlleleReadLikelihoods snp_matrix = score_site(site, {other_snp});
        REQUIRE(snp_matrix.rel(0, 1) == Approx(1.0));   // its own allele
        // One mismatched base against the other SNP allele, and both flanks still matched,
        // so it must score well above the deletion allele, which differs by a whole node.
        REQUIRE(snp_matrix.rel(0, 0) > snp_matrix.rel(0, 2));
    }
}

/// A site with several alternative interior paths between two boundary nodes, and a
/// configurable tail. Each allele is a list of interior node sequences; an empty list is a
/// bypass. Node ids run in construction order and the tail is last.
struct MultiAlleleSite {
    bdsg::HashGraph graph;
    Snarl snarl;
    vector<SnarlTraversal> traversals;
    unique_ptr<SnarlManager> manager;
    vector<vector<pair<nid_t, bool>>> read_paths;

    MultiAlleleSite(const vector<vector<string>>& alleles, const string& tail) {
        handle_t head = graph.create_handle("AAAACCCC", 1);
        nid_t next = 2;
        vector<vector<nid_t>> interiors;
        vector<handle_t> lasts;
        for (const vector<string>& a : alleles) {
            vector<nid_t> ids;
            handle_t prev = head;
            for (const string& seq : a) {
                handle_t h = graph.create_handle(seq, next);
                graph.create_edge(prev, h);
                prev = h;
                ids.push_back(next);
                ++next;
            }
            interiors.push_back(ids);
            lasts.push_back(prev);
        }
        nid_t tail_id = next;
        handle_t tail_h = graph.create_handle(tail, tail_id);
        for (handle_t l : lasts) {
            graph.create_edge(l, tail_h);
        }

        snarl.mutable_start()->set_node_id(1);
        snarl.mutable_end()->set_node_id(tail_id);
        snarl.set_type(ULTRABUBBLE);
        snarl.set_start_end_reachable(true);
        vector<Snarl> snarls{snarl};
        manager.reset(new SnarlManager(snarls.begin(), snarls.end()));

        traversals.resize(alleles.size());
        read_paths.resize(alleles.size());
        for (size_t i = 0; i < alleles.size(); ++i) {
            vector<nid_t> path{1};
            path.insert(path.end(), interiors[i].begin(), interiors[i].end());
            path.push_back(tail_id);
            for (nid_t id : path) {
                Visit* v = traversals[i].add_visit();
                v->set_node_id(id);
                v->set_backward(false);
                read_paths[i].push_back({id, false});
            }
        }
    }
};

TEST_CASE("Allele sequence between two anchors is charged, outside them is not",
          "[allele_likelihood][scoring]") {
    // Whether an allele node must be paid for depends on where it sits relative to the read's
    // matched node visits, and getting that wrong is silent. A draft of this walk let the
    // correspondence stop early and leave allele nodes unconsumed anywhere, so a read spanning
    // a deletion scored the reference allele at rel = 1.0 -- preferring neither -- and the
    // homozygous deletion stopped being callable at all. Nothing else in the suite caught it.
    //
    // Between two matched visits the allele's extra nodes are sequence the read skipped: a
    // deletion, and charged. Before the first match or after the last they are simply beyond
    // the read's window, and charging them would penalise a read for being short.
    SnpAndDeletionSite site;

    SECTION("an internal skip is a deletion and is charged") {
        // Reads node 1 then node 4, skipping the SNP node the reference allele carries. Both
        // flanks ARE anchors here, so node 2 is strictly internal.
        Alignment deleting = make_matching_alignment(site.graph, "del", {{1, false}, {4, false}});
        AlleleReadLikelihoods matrix = score_site(site, {deleting});
        REQUIRE(matrix.rel(0, 2) == Approx(1.0));   // the deletion allele, matched exactly
        REQUIRE(matrix.rel(0, 0) < 1.0);            // reference: one node deleted, charged
        REQUIRE(matrix.rel(0, 1) < 1.0);
    }

    SECTION("and allele sequence past the last anchor is not") {
        // Stops after the SNP node, so the reference allele's node 4 lies beyond the read's
        // last anchor. That is outside the window and must cost nothing -- the reference
        // allele has to stay a perfect fit. (A read touching only a boundary node is dropped
        // as uninformative before it reaches scoring, so the read has to enter the site.)
        Alignment partial = make_matching_alignment(site.graph, "partial",
                                                    {{1, false}, {2, false}});
        AlleleReadLikelihoods matrix = score_site(site, {partial});
        REQUIRE(matrix.num_reads() == 1);
        REQUIRE(matrix.rel(0, 0) == Approx(1.0));   // reference: node 4 unreached, not charged
        REQUIRE(matrix.rel(0, 2) < 1.0);            // deletion allele lacks node 2: charged
    }
}

TEST_CASE("No read's allele preference depends on the flank's length",
          "[allele_likelihood][scoring]") {
    // The systematic version of the anchor-desync regression above. That test pins one
    // configuration; this one asserts the invariant across many, because the walk is
    // greedy and single-pass and the fixed desync was only one way for a greedy walk to
    // pick a bad correspondence.
    //
    // The invariant: rel is normalised by each read's own best allele, and lengthening a
    // node EVERY allele shares adds the same match credit to all of them. So no rel value
    // may move. Any walk that charges shared flank against one allele and not another --
    // which is exactly what the desync did -- breaks it.
    const vector<vector<vector<string>>> configurations = {
        {{"T"}, {"G"}},                          // SNP
        {{"T"}, {}},                             // one-base insertion against a bypass
        {{"T", "C"}, {}},                        // two adjacent inserted nodes
        {{"T", "C"}, {"T"}},                     // one extra node beside a shared one
        {{"TTTT"}, {"T"}},                       // unequal-length substituted nodes
        {{"T", "C", "G"}, {"T", "G"}},           // an extra node in the middle
        {{"T", "C"}, {"C", "T"}},                // same nodes, different order
        {{}, {"A"}, {"AA"}, {"AAA"}},            // a homopolymer ladder, four alleles
    };

    for (size_t c = 0; c < configurations.size(); ++c) {
        MultiAlleleSite shortf(configurations[c], "GGGGTTTT");
        MultiAlleleSite longf(configurations[c], "GGGGTTTT" + string(40, 'A'));

        vector<Alignment> short_reads, long_reads;
        for (size_t i = 0; i < configurations[c].size(); ++i) {
            short_reads.push_back(make_matching_alignment(
                shortf.graph, "r" + std::to_string(i), shortf.read_paths[i]));
            long_reads.push_back(make_matching_alignment(
                longf.graph, "r" + std::to_string(i), longf.read_paths[i]));
        }

        InMemorySiteReadSource short_src, long_src;
        for (const Alignment& a : short_reads) short_src.add(a);
        for (const Alignment& a : long_reads) long_src.add(a);
        QualAdjAlignmentScorer qs;
        MatrixAlignmentScorer ps;
        // Both walks: greedy is the default and the optimal one is what --realign selects.
        for (bool realign : {false, true}) {
        AlleleLikelihoodParams params;
        params.realign = realign;
        GraphAlignedAlleleLikelihoodCalculator short_calc(shortf.graph, *shortf.manager,
                                                          short_src, qs, ps, params);
        GraphAlignedAlleleLikelihoodCalculator long_calc(longf.graph, *longf.manager,
                                                         long_src, qs, ps, params);
        AlleleReadLikelihoods sm = short_calc.compute(shortf.snarl, shortf.traversals, 2);
        AlleleReadLikelihoods lm = long_calc.compute(longf.snarl, longf.traversals, 2);

        INFO("configuration " << c << (realign ? " (--realign)" : " (greedy)"));
        REQUIRE(sm.num_reads() == lm.num_reads());
        REQUIRE(sm.num_alleles() == lm.num_alleles());
        for (size_t r = 0; r < sm.num_reads(); ++r) {
            for (size_t a = 0; a < sm.num_alleles(); ++a) {
                INFO("config " << c << " read " << r << " allele " << a
                                << (realign ? " (--realign)" : " (greedy)"));
                REQUIRE(sm.rel(r, a) == Approx(lm.rel(r, a)));
            }
        }
        }
    }
}

TEST_CASE("Every read is placeable against every allele, whatever the node layout",
          "[allele_likelihood][scoring]") {
    // An unplaceable read contributes -inf, which normalises to a relative likelihood of
    // exactly 0 -- and 0 is indistinguishable from "scored, and hopeless". So a walk that
    // cannot reach an allele at all fails silently: the genotype simply never considers it.
    //
    // Every defect in the walk's state machine has surfaced here first, and in several cases
    // only here. Restricting which allele columns a read step may reach -- an optimisation
    // tried and reverted -- produced exactly this, as did forbidding the transitions that let
    // a read cross a deleted node. Neither moved any other assertion in this file.
    //
    // So: every read must reach every allele. A read may of course prefer one strongly, but a
    // relative likelihood of 0 means the walk could not get there at all.
    const vector<vector<vector<string>>> configurations = {
        {{"T"}, {"G"}},                          // SNP
        {{"T"}, {}},                             // one-base insertion against a bypass
        {{"T", "C"}, {}},                        // two adjacent inserted nodes
        {{"T", "C"}, {"T"}},                     // one extra node beside a shared one
        {{"TTTT"}, {"T"}},                       // unequal-length substituted nodes
        {{"T", "C", "G"}, {"T", "G"}},           // an extra node in the middle
        {{"T", "C"}, {"C", "T"}},                // same nodes, different order
        {{}, {"A"}, {"AA"}, {"AAA"}},            // a homopolymer ladder, four alleles
        {{"ACGTACGTAC"}, {}},                    // a ten-base deletion to walk across
        {{"AC", "GT", "AC", "GT"}, {"AC", "GT"}},// two nodes deleted from a run of four
        {{"A", "C", "G", "T"}, {"T", "G", "C", "A"}},   // four nodes, reversed order
    };

    for (size_t c = 0; c < configurations.size(); ++c) {
        MultiAlleleSite site(configurations[c], "GGGGTTTT");
        vector<Alignment> reads;
        for (size_t i = 0; i < configurations[c].size(); ++i) {
            reads.push_back(make_matching_alignment(site.graph, "r" + std::to_string(i),
                                                    site.read_paths[i]));
        }
        InMemorySiteReadSource src;
        for (const Alignment& a : reads) src.add(a);
        QualAdjAlignmentScorer qs;
        MatrixAlignmentScorer ps;
        for (bool realign : {false, true}) {
        AlleleLikelihoodParams params;
        params.realign = realign;
        GraphAlignedAlleleLikelihoodCalculator calc(site.graph, *site.manager, src, qs, ps, params);
        AlleleReadLikelihoods m = calc.compute(site.snarl, site.traversals, 2);

        INFO("configuration " << c << (realign ? " (--realign)" : " (greedy)"));
        REQUIRE(m.num_reads() == configurations[c].size());
        for (size_t r = 0; r < m.num_reads(); ++r) {
            // A read built to follow allele r's own path matches it exactly, so it is that
            // read's best allele and rel is 1 by construction. If the bounds were to forbid
            // the very pairing the read was built from, this is where it shows.
            INFO("config " << c << " read " << r << " against its own allele");
            REQUIRE(m.rel(r, r) == Approx(1.0));
            for (size_t a = 0; a < m.num_alleles(); ++a) {
                INFO("config " << c << " read " << r << " allele " << a
                                << (realign ? " (--realign)" : " (greedy)"));
                REQUIRE(m.rel(r, a) > 0.0);
            }
        }
        }
    }
}

TEST_CASE("The optimal walk keeps the indel invariants the greedy one has",
          "[allele_likelihood][scoring]") {
    // --realign changes how the read-to-allele correspondence is SEARCHED for, not what a
    // correspondence costs. So the properties pinned for the greedy walk have to survive it,
    // and the anchor-desync regression above -- direction symmetry of a one-base indel -- is
    // the sharpest of them, being parameter-free.
    SnpAndDeletionSite site;
    Alignment spanning = make_matching_alignment(site.graph, "spanning",
                                                 {{1, false}, {2, false}, {4, false}});
    Alignment deleting = make_matching_alignment(site.graph, "deleting", {{1, false}, {4, false}});

    AlleleReadLikelihoods greedy = score_site(site, {spanning, deleting}, 2, false);
    AlleleReadLikelihoods exact = score_site(site, {spanning, deleting}, 2, true);

    for (const auto& named : {make_pair("greedy", &greedy), make_pair("--realign", &exact)}) {
        const AlleleReadLikelihoods& m = *named.second;
        INFO(named.first);
        // Each read matches its own allele exactly.
        REQUIRE(m.rel(0, 0) == Approx(1.0));
        REQUIRE(m.rel(1, 2) == Approx(1.0));
        // A one-base insertion costs a few score units, not tens, and stays within one match
        // unit of the one-base deletion: the same event seen from the two sides.
        REQUIRE(m.rel(0, 2) > 1e-10);
        REQUIRE(m.rel(0, 2) < m.rel(1, 0));
        REQUIRE(m.rel(0, 2) > 0.1 * m.rel(1, 0));
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
