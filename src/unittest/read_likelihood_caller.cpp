/// \file read_likelihood_caller.cpp
///
/// Unit tests for the layer that turns a reads x alleles matrix into a genotype and
/// its quality fields.
///
/// The matrix itself is covered by `[allele_likelihood]`, and the whole pipeline by
/// `t/18_vg_call.t`. Neither covers this: the integration tests assert that GQ, GL and
/// GQI are *present* and well-formed, not what they contain, and the matrix tests stop
/// before the discounts. So the arithmetic that decides how confident a call looks --
/// two multiplicative discounts, one of them gated on allele size -- had no test that
/// could fail if it were wrong.
///
/// That gap matters more here than the line count suggests. The two worst defects in
/// this work were both index confusions in exactly this kind of glue: genotype
/// likelihoods keyed by traversal index and read as VCF allele indices, twice, once
/// corrupting the heap and once silently changing 34% of confident genotypes. Tests
/// that pin the *meaning* of each field are the cheapest guard against a third.

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

    // The discounts are multiplicative factors in [0,1], so this direction is an
    // invariant rather than a property of these reads. It is asserted because the
    // share is a ratio that floating-point accumulation can push a hair above 1, which
    // would *raise* GQ -- the clamp exists for that and this is what notices if it goes.
    REQUIRE(strong.info->gq <= strong.info->gq_undiscounted + 1e-9);
    REQUIRE(weak.info->gq <= weak.info->gq_undiscounted + 1e-9);
    REQUIRE(strong.info->explained_share <= 1.0);
}

TEST_CASE("--no-share-quality makes GQ the raw ratio, and the share still reports",
          "[read_likelihood_caller]") {
    CallerSite site;
    vector<Alignment> reads;
    // The share only falls below 1 when reads prefer an allele the *call* does not
    // contain -- a two-way split calls the heterozygote, which explains everything and
    // leaves the discount a no-op. So: ten reference reads, ten SNP reads, and four that
    // support the deletion. The call is 0/1 and those four are unexplained.
    //
    // Checked rather than assumed: the first version of this test used the two-way split
    // and passed with share exactly 1, asserting nothing about the discount at all.
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

    // Every genotype scored is present, keyed by its *sorted traversal indices*. The
    // two heap-corrupting bugs in this work both came from reading one index space as
    // another, so the key's identity is worth asserting rather than assuming.
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

}  // namespace unittest
}  // namespace vg
