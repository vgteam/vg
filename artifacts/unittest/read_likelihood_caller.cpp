/// \file read_likelihood_caller.cpp
///
/// Unit tests for ReadLikelihoodSnarlCaller, which turns a matrix of read and allele likelihoods
/// into a genotype and its quality fields.
///
/// The matrix itself is tested under `[allele_likelihood]`. These tests check the quality
/// fields. GQ is multiplied by the explained share and by the depth discount. The tests also
/// check the keys of the genotype likelihoods.

#include <cmath>
#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include "../allele_likelihood.hpp"
#include "../../src/alignment_scorer.hpp"
#include "../../src/unittest/catch.hpp"
#include "support/snarl_fixtures.hpp"
#include "../read_likelihood_caller.hpp"
#include "../site_read_source.hpp"
#include "../../src/snarls.hpp"
#include "../../src/traversal_support.hpp"
#include "../../src/utility.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// A snarl from node 1 to node 4 with three traversals: the reference through node 2, a
/// substitution through node 3, and a deletion along the edge from node 1 to node 4.
///
/// Nodes 2 and 3 hold 60 bp each, so the substitution has the reference's length and the
/// deletion is 60 bp shorter. The depth discount applies to a call only when a called allele
/// differs in length from the reference allele by at least the minimum length given to the
/// caller. The lengths leave out the boundary nodes 1 and 4. With a minimum of 50 bp, the
/// discount can apply to a call of the deletion and never to a call of the substitution.
struct CallerSite {
    bdsg::HashGraph graph;
    Snarl snarl;
    /// The traversals by index: 0 is the reference (1,2,4), 1 the substitution (1,3,4), and 2
    /// the deletion (1,4).
    vector<SnarlTraversal> traversals;

    CallerSite() {
        string interior(60, 'T');
        string interior_alt(60, 'G');
        graph.create_handle("AAAACCCCAAAACCCC", 1);
        graph.create_handle(interior, 2);
        graph.create_handle(interior_alt, 3);
        graph.create_handle("GGGGTTTTGGGGTTTT", 4);
        graph.create_edge(graph.get_handle(1), graph.get_handle(2));
        graph.create_edge(graph.get_handle(1), graph.get_handle(3));
        graph.create_edge(graph.get_handle(1), graph.get_handle(4));
        graph.create_edge(graph.get_handle(2), graph.get_handle(4));
        graph.create_edge(graph.get_handle(3), graph.get_handle(4));

        snarl.mutable_start()->set_node_id(1);
        snarl.mutable_end()->set_node_id(4);
        snarl.set_type(ULTRABUBBLE);

        auto make_trav = [&](const vector<nid_t>& nodes) {
            SnarlTraversal t;
            for (nid_t n : nodes) {
                Visit* v = t.add_visit();
                v->set_node_id(n);
                v->set_backward(false);
            }
            return t;
        };
        traversals.push_back(make_trav({1, 2, 4}));
        traversals.push_back(make_trav({1, 3, 4}));
        traversals.push_back(make_trav({1, 4}));
    }
};

/// An all-match alignment along the given nodes.
Alignment matching_read(const HandleGraph& graph, const string& name,
                        const vector<nid_t>& nodes, int mapq = 60) {
    Alignment aln;
    aln.set_name(name);
    string seq;
    for (nid_t n : nodes) {
        string node_seq = graph.get_sequence(graph.get_handle(n, false));
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(n);
        m->mutable_position()->set_is_reverse(false);
        m->mutable_position()->set_offset(0);
        Edit* e = m->add_edit();
        e->set_from_length(node_seq.size());
        e->set_to_length(node_seq.size());
        seq += node_seq;
    }
    aln.set_sequence(seq);
    aln.set_quality(string(seq.size(), (char)30));
    aln.set_mapping_quality(mapq);
    return aln;
}

/// The result of genotyping one site with `call_site`.
struct Called {
    vector<int> genotype;
    /// The call info that `genotype_at` returned with the genotype.
    unique_ptr<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo> info;
    /// The caller's `gq_factor` for `info`, under the settings it was called with.
    double gq_factor = 0.0;
};

/// Genotype `site` from `reads` with a new caller. `configure`, where given, changes the
/// caller's settings first.
Called call_site(CallerSite& site, const vector<Alignment>& reads,
                 const function<void(ReadLikelihoodSnarlCaller&)>& configure = nullptr,
                 int ploidy = 2, int region_ploidy = 0, bool also_score_other = false) {
    InMemorySiteReadSource source;
    for (const Alignment& aln : reads) {
        source.add(aln);
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, source, qual_scorer,
                                                      plain_scorer);
    ReadLikelihoodSnarlCaller caller(site.graph, calculator);
    if (configure) {
        configure(caller);
    }
    auto result = caller.genotype_at(
        materialized_of(site.graph, site.snarl,
            CandidateAlleles{walks_of(site.graph, site.traversals), 0}),
        Ploidies{.ploidy = ploidy,
                 .region_ploidy = region_ploidy,
                 .also_score_other = also_score_other},
        /* top-level: nothing encloses it */ {});
    Called out;
    out.genotype = result.first;
    out.info = std::move(result.second);
    if (out.info != nullptr) {
        out.gq_factor = caller.gq_factor(*out.info);
    }
    return out;
}

}  // namespace

TEST_CASE("GQ rises with the evidence and never exceeds the undiscounted value",
          "[read_likelihood_caller]") {
    CallerSite site;

    vector<Alignment> few;
    for (int i = 0; i < 4; ++i) {
        few.push_back(matching_read(site.graph, "f" + std::to_string(i), {1, 2, 4}));
    }
    vector<Alignment> many;
    for (int i = 0; i < 40; ++i) {
        many.push_back(matching_read(site.graph, "m" + std::to_string(i), {1, 2, 4}));
    }

    Called weak = call_site(site, few);
    Called strong = call_site(site, many);
    REQUIRE(weak.info != nullptr);
    REQUIRE(strong.info != nullptr);

    // Same unanimous answer either way; ten times the reads must not make it less sure.
    REQUIRE(weak.genotype == vector<int>({0, 0}));
    REQUIRE(strong.genotype == vector<int>({0, 0}));
    REQUIRE(strong.info->gq >= weak.info->gq);

    // The discounts are factors in [0, 1], so GQ never exceeds GQI. Rounding can push the
    // share just above 1, which the clamp on it prevents.
    REQUIRE(strong.info->gq <= strong.info->gq_undiscounted + 1e-9);
    REQUIRE(weak.info->gq <= weak.info->gq_undiscounted + 1e-9);
    REQUIRE(strong.info->explained_share <= 1.0);
}

TEST_CASE("GQ is GQI times the explained share, or GQI with the share discount off",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    // The share falls below 1 only when some reads prefer an allele the call does not
    // contain. Ten reference reads, ten substitution reads and four deletion reads give the
    // call 0/1, which leaves the four deletion reads unexplained.
    for (int i = 0; i < 10; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 2, 4}));
        reads.push_back(matching_read(site.graph, "a" + std::to_string(i), {1, 3, 4}));
    }
    for (int i = 0; i < 4; ++i) {
        reads.push_back(matching_read(site.graph, "d" + std::to_string(i), {1, 4}));
    }

    Called discounted = call_site(site, reads);
    Called raw = call_site(site, reads,
                           [](ReadLikelihoodSnarlCaller& c) { c.set_share_discount(false); });
    REQUIRE(discounted.info != nullptr);
    REQUIRE(raw.info != nullptr);

    REQUIRE(raw.info->gq == Approx(raw.info->gq_undiscounted));
    // GQI is the undiscounted value in both cases: turning the discount off must change
    // GQ, not the record of what GQ would have been.
    REQUIRE(discounted.info->gq_undiscounted == Approx(raw.info->gq_undiscounted));
    // Strict, not <=: with unexplained reads present the discount must actually bite,
    // so an accidentally disconnected multiplication fails here rather than passing.
    REQUIRE(discounted.info->explained_share < 1.0);
    REQUIRE(discounted.info->gq < raw.info->gq);
    REQUIRE(discounted.info->gq == Approx(raw.info->gq * discounted.info->explained_share));
}

TEST_CASE("GQ and GQN can be recomputed from gq_factor and achievable_gap",
          "[read_likelihood_caller]") {
    // Where the linkage model writes a genotype other than the direct call, the record's GQ and
    // GQN are rewritten from `gq_factor` and `achievable_gap`. So GQ must be GQI times
    // `gq_factor`, and GQN must follow from GQI, `achievable_gap` and the explained share. Both
    // must hold under every setting that changes GQ.
    CallerSite site;
    vector<Alignment> mixed;
    for (int i = 0; i < 10; ++i) {
        mixed.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 2, 4}));
        mixed.push_back(matching_read(site.graph, "a" + std::to_string(i), {1, 3, 4}));
    }
    for (int i = 0; i < 4; ++i) {
        mixed.push_back(matching_read(site.graph, "d" + std::to_string(i), {1, 4}));
    }
    vector<Alignment> deletion;
    for (int i = 0; i < 12; ++i) {
        deletion.push_back(matching_read(site.graph, "x" + std::to_string(i), {1, 4}));
    }
    auto no_share = [](ReadLikelihoodSnarlCaller& c) { c.set_share_discount(false); };
    auto depth = [](ReadLikelihoodSnarlCaller& c) { c.set_depth_quality(1.0, 50); };
    auto both = [](ReadLikelihoodSnarlCaller& c) {
        c.set_share_discount(false);
        c.set_depth_quality(1.0, 50);
    };

    Called shared = call_site(site, mixed);
    Called raw = call_site(site, mixed, no_share);
    Called deep = call_site(site, deletion, depth);
    Called deep_raw = call_site(site, deletion, both);
    for (const Called* c : {&shared, &raw, &deep, &deep_raw}) {
        REQUIRE(c->info != nullptr);
        REQUIRE(c->info->gq == Approx(c->info->gq_undiscounted * c->gq_factor));
        // GQN is the gap in nats over the achievable gap, held at 1, times the share.
        REQUIRE(c->info->achievable_gap > 0.0);
        const double gap_nats = c->info->gq_undiscounted * log(10.0) / 10.0;
        REQUIRE(c->info->gq_fraction
                == Approx(min(1.0, gap_nats / c->info->achievable_gap)
                          * c->info->explained_share));
    }
    // The share is in the factor only while the share discount is on.
    REQUIRE(shared.info->explained_share < 1.0);
    REQUIRE(shared.gq_factor == Approx(shared.info->explained_share));
    REQUIRE(raw.gq_factor == Approx(1.0));
    // The depth discount is in it either way.
    REQUIRE(deep_raw.info->depth_discount < 1.0);
    REQUIRE(deep_raw.gq_factor == Approx(deep_raw.info->depth_discount));
}

TEST_CASE("The depth discount applies only to calls that change length by at least the minimum",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 20; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 2, 4}));
    }

    // Every read is on the reference, so the call is 0/0 and no called allele changes length.
    // The discount is meant for ranking structural variants, and its length condition keeps it
    // off small variants. With a minimum length of 50 bp, it must leave this call's GQ alone.
    Called min_50 = call_site(site, reads, [](ReadLikelihoodSnarlCaller& c) {
        c.set_depth_quality(1.0, 50);
    });
    REQUIRE(min_50.info != nullptr);
    REQUIRE(min_50.genotype == vector<int>({0, 0}));

    Called discount_off = call_site(site, reads);
    REQUIRE(discount_off.info != nullptr);
    REQUIRE(min_50.info->gq == Approx(discount_off.info->gq));

    // With a minimum length of 0, the discount applies to every call, and it lowers this call's
    // GQ. So the length of the call is what kept the discount off it above.
    Called min_0 = call_site(site, reads, [](ReadLikelihoodSnarlCaller& c) {
        c.set_depth_quality(1.0, 0);
    });
    REQUIRE(min_0.info != nullptr);
    REQUIRE(min_0.info->gq < discount_off.info->gq);
}

TEST_CASE("A depth discount with an exponent of 0 leaves GQ unchanged",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 12; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 4}));
    }

    // Every read carries the deletion, a length change of 60 bp, which is at least the 50 bp
    // minimum. With an exponent of 1 the discount lowers this call's GQ. With an exponent of 0
    // the discount is 1, so GQ is GQI times the explained share alone.
    Called one = call_site(site, reads, [](ReadLikelihoodSnarlCaller& c) {
        c.set_depth_quality(1.0, 50);
    });
    Called zero = call_site(site, reads, [](ReadLikelihoodSnarlCaller& c) {
        c.set_depth_quality(0.0, 50);
    });
    REQUIRE(one.info != nullptr);
    REQUIRE(zero.info != nullptr);
    REQUIRE(one.info->gq < zero.info->gq);
    REQUIRE(zero.info->gq == Approx(zero.info->gq_undiscounted * zero.info->explained_share));
}

TEST_CASE("Genotype likelihoods are keyed by the genotype that was scored",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 15; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 3, 4}));
    }

    Called called = call_site(site, reads);
    REQUIRE(called.info != nullptr);
    REQUIRE(called.genotype == vector<int>({1, 1}));

    // Every genotype scored is present, keyed by its sorted traversal indices, not by
    // VCF allele indices.
    REQUIRE(called.info->genotype_lls.count(vector<int>({1, 1})) == 1);
    REQUIRE(called.info->genotype_lls.count(vector<int>({0, 0})) == 1);
    REQUIRE(called.info->genotype_lls.count(vector<int>({0, 1})) == 1);

    // The called genotype is the argmax over what was scored, by construction.
    double best = called.info->genotype_lls.at(vector<int>({1, 1}));
    for (const auto& entry : called.info->genotype_lls) {
        REQUIRE(entry.second <= best + 1e-9);
    }
}

TEST_CASE("Haploid sites are genotyped with one allele and still get a quality",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 15; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 3, 4}));
    }

    Called called = call_site(site, reads, nullptr, 1);
    REQUIRE(called.info != nullptr);
    REQUIRE(called.genotype.size() == 1);
    REQUIRE(called.genotype[0] == 1);
    REQUIRE(called.info->gq >= 0.0);
    REQUIRE(called.info->gq <= called.info->gq_undiscounted + 1e-9);
}

TEST_CASE("A nested site's depth rate is per haplotype of the region, not of the site",
          "[read_likelihood_caller]") {
    // A nested site that only one of a diploid parent's alleles crosses is genotyped at
    // ploidy 1, but the reads its depth rate is measured from come from both haplotypes.
    // With the region's ploidy, one copy of the allele predicts half the reads that two
    // copies do, so the same reads give twice the DR. With the site's own ploidy, the rate
    // would double and the two DRs would be equal.
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 15; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 3, 4}));
    }

    Called diploid = call_site(site, reads);
    Called nested = call_site(site, reads, nullptr, 1, 2, true);

    REQUIRE(diploid.info != nullptr);
    REQUIRE(nested.info != nullptr);
    REQUIRE(diploid.genotype == vector<int>({1, 1}));
    REQUIRE(nested.genotype == vector<int>({1}));
    REQUIRE(diploid.info->depth_ratio > 0.0);
    REQUIRE(nested.info->depth_ratio == Approx(2.0 * diploid.info->depth_ratio));

    // The site at ploidy 2, which the linkage pass takes if both parent alleles turn out to cross
    // it, predicts what the diploid call does.
    REQUIRE(nested.info->alt_ploidy_info != nullptr);
    REQUIRE(nested.info->alt_ploidy_info->depth_ratio == Approx(diploid.info->depth_ratio));
}

TEST_CASE("DR can be computed for a genotype other than the direct call",
          "[read_likelihood_caller]") {
    // The linkage model can write a genotype other than the direct call, and DR must then
    // describe the written one. The CallInfo keeps the depth context of the matrix, which gives DR
    // for any genotype, and gives back exactly the direct call's DR for the direct call.
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 15; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 3, 4}));
    }
    Called called = call_site(site, reads, nullptr, 2, 0, true);
    REQUIRE(called.info != nullptr);
    REQUIRE(called.genotype == vector<int>({1, 1}));
    REQUIRE(called.info->depth_ratio > 0.0);
    REQUIRE(called.info->depth.ratio({1, 1}) == called.info->depth_ratio);

    // Another genotype, against lambda's formula.
    const auto& info = *called.info;
    const DepthContext& depth = info.depth;
    REQUIRE(!depth.lengths.empty());
    double lambda_01 = depth.rate * ((double)depth.lengths[0] + depth.read_length - 1.0
                                     + (double)depth.lengths[1] + depth.read_length - 1.0);
    REQUIRE(info.depth.ratio({0, 1}) == Approx(depth.observed / lambda_01));

    // The alternate ploidy's CallInfo, which the linkage pass can promote, carries the same context.
    REQUIRE(info.alt_ploidy_info != nullptr);
    REQUIRE(info.alt_ploidy_info->depth.ratio({1, 1}) == info.depth_ratio);
}

TEST_CASE("Recomputed GQ takes the explained share of the new best genotype",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    // As above: the call is 0/1 and the four deletion reads are unexplained.
    for (int i = 0; i < 10; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 2, 4}));
        reads.push_back(matching_read(site.graph, "a" + std::to_string(i), {1, 3, 4}));
    }
    for (int i = 0; i < 4; ++i) {
        reads.push_back(matching_read(site.graph, "d" + std::to_string(i), {1, 4}));
    }
    Called called = call_site(site, reads);
    REQUIRE(called.info != nullptr);
    REQUIRE(called.genotype == vector<int>({0, 1}));
    auto& info = *called.info;

    // Re-genotyping changes the likelihoods in place; here 1/1 is made the best genotype, with
    // 0/1 second.
    const double het_ll = info.genotype_lls.at(vector<int>({0, 1}));
    info.genotype_lls[vector<int>({1, 1})] = het_ll + 2.0;
    const double gap = logprob_to_phred(het_ll) - logprob_to_phred(het_ll + 2.0);
    // Only the substitution reads have allele 1 as their best allele.
    const double share = info.allele_support[1] / (double)info.n_reads;
    REQUIRE(share < 1.0);

    InMemorySiteReadSource source;
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, source, qual_scorer,
                                                      plain_scorer);
    ReadLikelihoodSnarlCaller caller(site.graph, calculator);

    SECTION("with the share discount") {
        caller.recompute_gq(info);
        REQUIRE(info.gq == Approx(gap * share));
        REQUIRE(info.gq < gap);
    }
    SECTION("without the share discount") {
        caller.set_share_discount(false);
        caller.recompute_gq(info);
        REQUIRE(info.gq == Approx(gap));
    }
}

TEST_CASE("A record with a star allele still gets QUAL from the reference column",
          "[read_likelihood_caller]") {
    // A nested site that one strand crosses is genotyped at ploidy 1 and written with a star
    // allele for the other strand. The star allele's empty traversal matches no scored column,
    // and QUAL needs only the reference allele's.
    CallerSite site;
    InMemorySiteReadSource source;
    for (int i = 0; i < 15; ++i) {
        source.add(matching_read(site.graph, "r" + std::to_string(i), {1, 3, 4}));
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, source, qual_scorer,
                                                      plain_scorer);
    ReadLikelihoodSnarlCaller caller(site.graph, calculator);
    auto result = caller.genotype_at(
        materialized_of(site.graph, site.snarl,
            CandidateAlleles{walks_of(site.graph, site.traversals), 0}),
        Ploidies{.ploidy = 1}, /* top-level: nothing encloses it */ {});
    REQUIRE(result.first == vector<int>({1}));
    const unique_ptr<SnarlCaller::CallInfo> call_info(std::move(result.second));

    // Record alleles: reference, the called alternative, and the star allele, which was not
    // scored.
    const vector<int> scored_allele{0, 1, -1};
    vcflib::Variant variant;
    caller.write_vcf_fields(scored_allele, {1, 2}, call_info.get(), "s", variant);
    REQUIRE(variant.quality > 10.0);
    // GL needs every genotype of the record's alleles, so it is left out.
    REQUIRE(variant.samples["s"].count("GL") == 0);
}

TEST_CASE("A moved record's quality fields come from the stored direct call",
          "[read_likelihood_caller]") {
    const auto rewrite = &ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype;
    // Chosen 1/1, so the chosen genotype's GL entry is -3 against a best other of -1: a margin of
    // -20 phred. The printed GQI is capped and the printed GQN is 0, so a divisor recovered from
    // them would be wrong or missing; the stored achievable gap is 100 phred.
    vcflib::Variant line;
    line.sequenceName = "chr1";
    line.position = 100;
    line.id = ">1>4";
    line.ref = "A";
    line.alt = {"G"};
    line.quality = 30;
    line.filter = "PASS";
    line.format = {"GT", "GL", "GQ", "GQI", "GQN"};
    line.samples["s"] = {{"GT", {"1/1"}}, {"GL", {"-5.000000", "-1.000000", "-3.000000"}},
                         {"GQ", {"40"}}, {"GQI", {"256"}}, {"GQN", {"0.000"}}};
    auto field = [](vcflib::Variant& record, const string& key) {
        return record.samples["s"][key].at(0);
    };
    auto filter = [](vcflib::Variant& record) {
        return record.filter;
    };
    LinkageCollector::MovedQuality moved;
    moved.posterior = 0.995;   // -10 log10(0.005) = 23.01
    moved.direct.explained_share = 0.5;
    moved.direct.achievable_gap = 100.0 * log(10.0) / 10.0;

    SECTION("GQ takes the direct call's GQ factor, whatever it is made of") {
        vcflib::Variant shared = line, unshared = line;
        moved.direct.gq_factor = 0.5;   // the share, as by default
        REQUIRE(rewrite(shared, "s", moved, 0.0));
        REQUIRE(field(shared, "GQ") == "11");
        moved.direct.gq_factor = 1.0;   // no share discount and no depth discount
        REQUIRE(rewrite(unshared, "s", moved, 0.0));
        REQUIRE(field(unshared, "GQ") == "23");
    }
    SECTION("GQN divides by the stored achievable gap and multiplies by the share") {
        vcflib::Variant rewritten = line;
        REQUIRE(rewrite(rewritten, "s", moved, 0.05));
        REQUIRE(field(rewritten, "GQN") == "-0.100");
        REQUIRE(filter(rewritten) == "lowconf");
    }
    SECTION("A record with no FORMAT fields for the sample is left alone") {
        vcflib::Variant other = line;
        REQUIRE(!rewrite(other, "t", moved, 0.05));
        REQUIRE(field(other, "GQ") == "40");
    }
    SECTION("GQN is missing where the direct call had no achievable gap") {
        vcflib::Variant rewritten = line;
        moved.direct.achievable_gap = 0.0;
        REQUIRE(rewrite(rewritten, "s", moved, 0.05));
        REQUIRE(field(rewritten, "GQN") == ".");
        REQUIRE(filter(rewritten) == "PASS");
    }
}

}  // namespace unittest
}  // namespace vg
