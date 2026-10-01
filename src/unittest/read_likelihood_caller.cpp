/// \file read_likelihood_caller.cpp
///
/// Unit tests for the layer that turns a reads x alleles matrix into a genotype and
/// its quality fields.
///
/// The matrix itself is tested under `[allele_likelihood]`. These tests check what the
/// quality fields contain: the two multiplicative discounts on GQ, one of them gated on
/// allele size, and the keys of the genotype likelihoods.

#include <vector>

#include <bdsg/hash_graph.hpp>

#include "allele_likelihood.hpp"
#include "alignment_scorer.hpp"
#include "catch.hpp"
#include "read_likelihood_caller.hpp"
#include "site_read_source.hpp"
#include "snarls.hpp"
#include "traversal_support.hpp"
#include "utility.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// A site with a reference path, a SNP alternative, and a deletion, so that both a
/// same-length and a length-changing call can be exercised. Node 4 is long enough that
/// dropping nodes 2/3 is a >= 50 bp change, which is what arms the depth discount.
struct CallerSite {
    bdsg::HashGraph graph;
    Snarl snarl;
    vector<SnarlTraversal> traversals;   // 0: ref (1,2,4)  1: alt (1,3,4)  2: del (1,4)
    unique_ptr<SnarlManager> manager;

    CallerSite() {
        // 60 bp of interior sequence, so ref -> deletion is a 60 bp length change.
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

        vector<Snarl> snarls{snarl};
        manager.reset(new SnarlManager(snarls.begin(), snarls.end()));
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

/// Genotype one site, returning the call and its info together so the fields can be
/// inspected. `configure` runs against the caller before genotyping.
struct Called {
    vector<int> genotype;
    const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo* info;
    unique_ptr<SnarlCaller::CallInfo> owned;
    /// The caller's `gq_factor` for `info`, under the settings it was called with.
    double gq_factor = 0.0;
};

Called call_site(CallerSite& site, const vector<Alignment>& reads,
                 const function<void(ReadLikelihoodSnarlCaller&)>& configure = nullptr,
                 int ploidy = 2) {
    InMemorySiteReadSource source;
    for (const Alignment& aln : reads) {
        source.add(aln);
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, *site.manager, source,
                                                      qual_scorer, plain_scorer);
    NullTraversalSupportFinder support(site.graph, *site.manager);
    ReadLikelihoodSnarlCaller caller(site.graph, *site.manager, support, calculator);
    if (configure) {
        configure(caller);
    }
    auto result = caller.genotype(site.snarl, site.traversals, 0, ploidy, "", {0, 0});
    Called out;
    out.genotype = result.first;
    out.owned = std::move(result.second);
    out.info = dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(
        out.owned.get());
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

TEST_CASE("--no-share-quality makes GQ the raw ratio, and the share still reports",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    // The share falls below 1 only when some reads prefer an allele the call does not
    // contain. Ten reference reads, ten SNP reads and four deletion reads give the call
    // 0/1, which leaves the four deletion reads unexplained.
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

TEST_CASE("gq_factor and achievable_gap are what GQ and GQN were computed with",
          "[read_likelihood_caller]") {
    // A moved record's quality fields are rewritten from these two, so each must reproduce the
    // per-site value it stands for, under every setting that changes GQ.
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

TEST_CASE("The depth discount is gated on the called allele's length change",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 20; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 2, 4}));
    }

    // Reference against SNP is a 0 bp change, so however implausible the depth, the
    // discount must not fire. This is the gate that keeps a ranking signal aimed at
    // structural variants from touching the SNVs, which are the bulk of every call set.
    Called snp_site = call_site(site, reads, [](ReadLikelihoodSnarlCaller& c) {
        c.set_depth_quality(1.0, 50);
    });
    REQUIRE(snp_site.info != nullptr);
    REQUIRE(snp_site.genotype == vector<int>({0, 0}));

    Called undiscounted = call_site(site, reads);
    REQUIRE(undiscounted.info != nullptr);
    // Identical GQ: the only difference between the runs is a discount that is gated
    // off. If the gate were removed or its comparison inverted, these would diverge.
    REQUIRE(snp_site.info->gq == Approx(undiscounted.info->gq));
}

TEST_CASE("A zero depth-quality exponent is inert", "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 12; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 4}));
    }

    Called off = call_site(site, reads);
    Called zero = call_site(site, reads, [](ReadLikelihoodSnarlCaller& c) {
        c.set_depth_quality(0.0, 50);
    });
    REQUIRE(off.info != nullptr);
    REQUIRE(zero.info != nullptr);
    REQUIRE(zero.info->gq == Approx(off.info->gq));
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
    ReadLikelihoodSnarlCaller::set_region_ploidy(2);
    ReadLikelihoodSnarlCaller::set_want_alt_ploidy(true);
    Called nested = call_site(site, reads, nullptr, 1);
    ReadLikelihoodSnarlCaller::set_want_alt_ploidy(false);
    ReadLikelihoodSnarlCaller::set_region_ploidy(0);

    REQUIRE(diploid.info != nullptr);
    REQUIRE(nested.info != nullptr);
    REQUIRE(diploid.genotype == vector<int>({1, 1}));
    REQUIRE(nested.genotype == vector<int>({1}));
    REQUIRE(diploid.info->depth_ratio > 0.0);
    REQUIRE(nested.info->depth_ratio == Approx(2.0 * diploid.info->depth_ratio));

    // The site at ploidy 2, which the barrier takes if both parent alleles turn out to cross
    // it, predicts what the diploid call does.
    REQUIRE(nested.info->alt_ploidy_info != nullptr);
    REQUIRE(nested.info->alt_ploidy_info->depth_ratio == Approx(diploid.info->depth_ratio));
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
    auto& info = dynamic_cast<ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo&>(*called.owned);

    // Re-genotyping changes the likelihoods in place; here 1/1 is made the best genotype, with
    // 0/1 second.
    const double het_ll = info.genotype_lls.at(vector<int>({0, 1}));
    info.genotype_lls[vector<int>({1, 1})] = het_ll + 2.0;
    const double gap = logprob_to_phred(het_ll) - logprob_to_phred(het_ll + 2.0);
    // Only the SNP reads have allele 1 as their best allele.
    const double share = info.allele_support[1] / (double)info.n_reads;
    REQUIRE(share < 1.0);

    InMemorySiteReadSource source;
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, *site.manager, source,
                                                      qual_scorer, plain_scorer);
    NullTraversalSupportFinder support(site.graph, *site.manager);
    ReadLikelihoodSnarlCaller caller(site.graph, *site.manager, support, calculator);

    SECTION("with the share discount") {
        caller.recompute_gq(info);
        REQUIRE(info.gq == Approx(gap * share));
        REQUIRE(info.gq < gap);
    }
    SECTION("with --no-share-quality") {
        caller.set_share_discount(false);
        caller.recompute_gq(info);
        REQUIRE(info.gq == Approx(gap));
    }
}

}  // namespace unittest
}  // namespace vg
