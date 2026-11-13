/// \file read_likelihood_caller.cpp
///
/// Unit tests for ReadLikelihoodSnarlCaller, which turns a matrix of read and allele likelihoods
/// into a genotype and its quality fields.
///
/// The matrix itself is tested under `[allele_likelihood]`. These tests check the quality
/// fields. GQ is multiplied by the explained share and by the depth discount. The tests also
/// check the keys of the genotype likelihoods.

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <set>
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
#include "../../src/statistics.hpp"
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
                 int ploidy = 2, int region_ploidy = 0, bool also_score_other = false,
                 const vector<int>* group_of = nullptr) {
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
    auto [genotype, call_info] = caller.genotype_at(
        materialized_of(site.graph, site.snarl,
            CandidateAlleles{walks_of(site.graph, site.traversals), 0}),
        Ploidies{.ploidy = ploidy,
                 .region_ploidy = region_ploidy,
                 .also_score_other = also_score_other},
        /* top-level: nothing encloses it */ {}, group_of);
    Called out;
    out.genotype = genotype;
    out.info = std::move(call_info);
    if (out.info != nullptr) {
        out.gq_factor = caller.gq_factor(*out.info);
    }
    return out;
}

/// Build the matrix that `call_site` genotypes `site` from at ploidy 2, from `reads`.
AlleleReadLikelihoods aligned_matrix(const CallerSite& site, const vector<Alignment>& reads) {
    InMemorySiteReadSource source;
    for (const Alignment& aln : reads) {
        source.add(aln);
    }
    QualAdjAlignmentScorer qual_scorer;
    MatrixAlignmentScorer plain_scorer;
    GraphAlignedAlleleLikelihoodCalculator calculator(site.graph, source, qual_scorer,
                                                      plain_scorer);
    return calculator.compute(
        materialized_of(site.graph, site.snarl,
            CandidateAlleles{walks_of(site.graph, site.traversals), 0}),
        /* top-level: nothing encloses it */ {}, 2);
}

/// A calculator that gives every site the matrix that `make` builds.
class GivenMatrix : public AlleleLikelihoodCalculator {
public:
    explicit GivenMatrix(function<AlleleReadLikelihoods()> make) : make(std::move(make)) {}

    AlleleReadLikelihoods compute(const MaterializedSite&, const vector<nid_t>&,
                                  int) override {
        return make();
    }

private:
    function<AlleleReadLikelihoods()> make;
};

/**
 * Return a function that builds a matrix of reads that each fit a set of alleles.
 *
 * Each entry of `reads` is a set of alleles and a number of reads. Each of those reads fits every
 * allele of the set equally, and every other allele e^20 times worse. Every read has e_r 0.02.
 * With `lengths`, the mixture is weighted by the length of each allele, for reads of 100 bp.
 */
function<AlleleReadLikelihoods()> reads_fitting(size_t alleles,
                                               const vector<pair<set<int>, int>>& reads,
                                               const vector<size_t>& lengths = {}) {
    return [alleles, reads, lengths]() {
        AlleleReadLikelihoodsBuilder builder(alleles, 0.02);
        builder.set_allele_lengths(lengths);
        int named = 0;
        for (const auto& [fits, count] : reads) {
            vector ln(alleles, -20.0);
            for (int a : fits) {
                ln[a] = 0.0;
            }
            for (int i = 0; i < count; ++i) {
                builder.add_read(ln, 1e-6, "r" + std::to_string(named++),
                                 lengths.empty() ? 0 : 100);
            }
        }
        return builder.build();
    };
}

/// Genotype a site whose matrix `make` builds, with `group_of` as the groups of its alleles.
///
/// Every walk of the site is the SNP walk of `CallerSite`, since the matrix alone decides the
/// genotype. Walk 0 is the reference.
Called call_matrix(const function<AlleleReadLikelihoods()>& make, int ploidy,
                   const vector<int>* group_of) {
    CallerSite site;
    GivenMatrix calculator(make);
    ReadLikelihoodSnarlCaller caller(site.graph, calculator);
    const size_t alleles = make().num_alleles();
    auto [genotype, call_info] = caller.genotype_at(
        materialized_of(site.graph, site.snarl,
            CandidateAlleles{
                vector<Traversal>(alleles, walks_of(site.graph, site.traversals)[1]), 0}),
        Ploidies{.ploidy = ploidy}, /* top-level: nothing encloses it */ {}, group_of);
    Called out;
    out.genotype = genotype;
    out.info = std::move(call_info);
    out.gq_factor = caller.gq_factor(*out.info);
    return out;
}

/// Return the best shared score of `genotypes` less their best own score, with the depth term when
/// `with_depth` is true.
double gain_of(const AlleleReadLikelihoods& matrix, const vector<vector<int>>& genotypes,
               const vector<int>& group_of, bool with_depth) {
    double shared = -numeric_limits<double>::infinity();
    double own = -numeric_limits<double>::infinity();
    for (const auto& [genotype, ll] : matrix.score_genotypes(genotypes, &group_of, with_depth)) {
        shared = max(shared, ll);
    }
    for (const auto& [genotype, ll] : matrix.score_genotypes(genotypes, nullptr, with_depth)) {
        own = max(own, ll);
    }
    return shared - own;
}

/// Return the sorted multiset of the groups `group_of` gives the alleles of `genotype`.
vector<int> groups_of(vector<int> genotype, const vector<int>& group_of) {
    for (int& allele : genotype) {
        allele = group_of[allele];
    }
    std::ranges::sort(genotype);
    return genotype;
}

/// The best shared and own scores among the genotypes of a multiset of groups, and the first
/// genotype in VCF GL order with each.
struct MultisetBest {
    double shared;
    double own;
    vector<int> by_shared;
    vector<int> by_own;
};

/// Return the MultisetBest of each multiset of groups at `ploidy`, with or without the depth term.
map<vector<int>, MultisetBest> best_by_multiset(const AlleleReadLikelihoods& matrix, int ploidy,
                                                const vector<int>& group_of, bool with_depth) {
    const auto genotypes = AlleleReadLikelihoods::enumerate_genotypes(matrix.num_alleles(), ploidy);
    const auto own = matrix.score_genotypes(genotypes, nullptr, with_depth);
    const auto shared = matrix.score_genotypes(genotypes, &group_of, with_depth);
    map<vector<int>, MultisetBest> best;
    for (size_t i = 0; i < own.size(); ++i) {
        const vector<int>& genotype = own[i].first;
        MultisetBest& b = best.try_emplace(groups_of(genotype, group_of), shared[i].second,
                                           own[i].second, genotype, genotype)
                              .first->second;
        if (shared[i].second > b.shared) {
            b.shared = shared[i].second;
            b.by_shared = genotype;
        }
        if (own[i].second > b.own) {
            b.own = own[i].second;
            b.by_own = genotype;
        }
    }
    return best;
}

/**
 * Return the gain of the multiset `plus`, taken over genotypes grown from best genotypes one copy
 * smaller.
 *
 * For each group g of two or more alleles that `plus` holds at least twice, `plus` less one copy of
 * g is a multiset in `reads`. Its best genotype by shared score and its best by own score, each
 * with any one allele of g added, are among the genotypes the gain is taken over.
 *
 * @param reads The MultisetBest of each multiset one copy smaller than `plus`, without the depth
 * term.
 */
double gain_above(const AlleleReadLikelihoods& matrix, const vector<int>& group_of,
                  const map<vector<int>, MultisetBest>& reads, const vector<int>& plus) {
    vector<vector<int>> extended;
    for (int g : plus) {
        if (std::ranges::count(group_of, g) < 2 || std::ranges::count(plus, g) < 2) {
            continue;
        }
        vector<int> less = plus;
        less.erase(std::ranges::find(less, g));
        const MultisetBest& from = reads.at(less);
        for (int allele = 0; allele < (int)group_of.size(); ++allele) {
            if (group_of[allele] != g) {
                continue;
            }
            for (const vector<int>& start : {from.by_shared, from.by_own}) {
                extended.push_back(start);
                extended.back().push_back(allele);
                std::ranges::sort(extended.back());
            }
        }
    }
    return gain_of(matrix, extended, group_of, false);
}

/**
 * Return the score of every genotype of `matrix` at `ploidy` when the groups `group_of` choose.
 *
 * A multiset of groups scores its shared score, the best ln likelihood of its genotypes with each
 * read taking its best relative likelihood within each group. With `charge`, it scores that less
 * the largest drop in its gain, shared less own score, that one more copy of one of its groups of
 * two or more alleles brings. Each genotype scores the score of its multiset less its shortfall
 * from the best genotype of the multiset. The gains, and the best genotypes they start from, come
 * from the scores without the depth term. `gain_above` gives the gain of a multiset at
 * `ploidy` + 1.
 */
map<vector<int>, double> group_scores(const AlleleReadLikelihoods& matrix, int ploidy,
                                      const vector<int>& group_of, bool charge) {
    const auto best = best_by_multiset(matrix, ploidy, group_of, true);
    const auto reads = best_by_multiset(matrix, ploidy, group_of, false);
    map<vector<int>, double> scores;
    for (const auto& [genotype, ll] : matrix.score_genotypes(ploidy)) {
        const vector<int> groups = groups_of(genotype, group_of);
        const MultisetBest& b = best.at(groups);
        const MultisetBest& r = reads.at(groups);
        double drop = 0.0;
        for (int group : groups) {
            if (charge && std::ranges::count(group_of, group) > 1) {
                vector<int> plus = groups;
                plus.push_back(group);
                std::ranges::sort(plus);
                drop = max(drop, (r.shared - r.own) - gain_above(matrix, group_of, reads, plus));
            }
        }
        scores[genotype] = b.shared - drop - (b.own - ll);
    }
    return scores;
}

/// Require each genotype of `called` to score as `expected` gives.
void require_scores(const Called& called, const map<vector<int>, double>& expected) {
    REQUIRE(called.info->genotype_lls.size() == expected.size());
    for (const auto& [genotype, ll] : expected) {
        REQUIRE(called.info->genotype_lls.at(genotype) == Approx(ll));
    }
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

TEST_CASE("Groups of a site's walks choose its genotype from aligned reads",
          "[read_likelihood_caller][group_candidates]") {
    // The reference and the SNP walk form one group, and the deletion a group of its own. In a
    // graph, `candidate_groups` joins two walks that differ only inside a child chain. The
    // haplotypes are the reference and the SNP, and a few reads fit the deletion.
    // Sharing lets one copy of the group carry the reads of both its walks, beside a copy of the
    // deletion. One more copy of the group takes over the reads of the SNP walk, so the multiset
    // of the group and the deletion is charged, and the groups call the reference and the SNP.
    CallerSite site;
    vector<Alignment> reads;
    for (int i = 0; i < 6; ++i) {
        reads.push_back(matching_read(site.graph, "r" + std::to_string(i), {1, 2, 4}));
    }
    for (int i = 0; i < 4; ++i) {
        reads.push_back(matching_read(site.graph, "s" + std::to_string(i), {1, 3, 4}));
    }
    for (int i = 0; i < 3; ++i) {
        reads.push_back(matching_read(site.graph, "d" + std::to_string(i), {1, 4}));
    }
    const vector group_of{0, 0, 1};

    Called walks = call_site(site, reads);
    Called groups = call_site(site, reads, nullptr, 2, 0, false, &group_of);
    REQUIRE(walks.genotype == vector<int>({0, 1}));
    REQUIRE(groups.genotype == vector<int>({0, 1}));
    REQUIRE(groups.info->gq > 0.0);
    // Without the charge, the group and the deletion would score above two copies of the group.
    const AlleleReadLikelihoods matrix = aligned_matrix(site, reads);
    const auto uncharged = group_scores(matrix, 2, group_of, false);
    REQUIRE(uncharged.at({0, 2}) > uncharged.at({0, 1}));
    require_scores(groups, group_scores(matrix, 2, group_of, true));
    // Genotypes with the same multiset of groups keep the differences between their likelihoods
    // from the walks.
    const auto gap = [](const Called& c) {
        return c.info->genotype_lls.at({0, 1}) - c.info->genotype_lls.at({0, 0});
    };
    REQUIRE(gap(groups) == Approx(gap(walks)));
}

TEST_CASE("Groups choose the genotype where one haplotype is a mosaic of a group's walks",
          "[read_likelihood_caller][group_candidates]") {
    // Walks 0 to 8 form group 0, and walks 9 and 10 are groups of their own. One haplotype is a
    // mosaic of the nine walks of group 0: each of its reads fits one of them. The other haplotype
    // is walk 10. Its reads fit walk 10, or walks 9 and 10, and some of them fit walk 9 alone.
    // No one walk explains the mosaic, so the walks alone call 9 and 10. Group 0 explains it.
    const auto make = reads_fitting(
        11, {{{0}, 4}, {{1}, 3}, {{2}, 3}, {{3}, 3}, {{4}, 3}, {{5}, 3}, {{6}, 3}, {{7}, 3},
             {{8}, 3}, {{9}, 6}, {{10}, 10}, {{9, 10}, 10}});
    const vector group_of{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2};

    Called walks = call_matrix(make, 2, nullptr);
    Called groups = call_matrix(make, 2, &group_of);
    REQUIRE(walks.genotype == vector<int>({9, 10}));
    // Within group 0, walk 0 has the most reads.
    REQUIRE(groups.genotype == vector<int>({0, 10}));
    // A second copy of group 0 lets the own score explain the reads of one more of its nine walks,
    // but lowers the mixture weight of every copy from 1/2 to 1/3. The own score loses more than it
    // gains, so the gain from sharing does not drop, and no multiset of group 0 and another group
    // is charged.
    const auto uncharged = group_scores(make(), 2, group_of, false);
    REQUIRE(groups.info->genotype_lls.at({0, 10}) == Approx(uncharged.at({0, 10})));
    REQUIRE(groups.info->genotype_lls.at({3, 9}) == Approx(uncharged.at({3, 9})));
    require_scores(groups, group_scores(make(), 2, group_of, true));
    // Genotypes with the same multiset of groups keep the differences between their likelihoods
    // from the walks.
    const auto gap = [](const Called& c) {
        return c.info->genotype_lls.at({0, 10}) - c.info->genotype_lls.at({1, 10});
    };
    REQUIRE(gap(groups) == Approx(gap(walks)));
}

TEST_CASE("One copy of a group cannot carry the reads of two haplotypes",
          "[read_likelihood_caller][group_candidates]") {
    // Walks 0 and 1 are 1,000 bp long and form one group. Walk 2 is a deletion of them, with a
    // group of its own. The haplotypes are walks 0 and 1, and a few reads fit the deletion. With
    // the deletion, a copy of the group takes nearly all the weight of the mixture. Sharing would
    // let that one copy carry the reads of both walks, and fit the deletion's reads as well.
    const auto make = reads_fitting(3, {{{0}, 20}, {{1}, 20}, {{2}, 5}}, {1000, 1000, 0});
    const vector group_of{0, 0, 1};

    const auto uncharged = group_scores(make(), 2, group_of, false);
    REQUIRE(uncharged.at({0, 2}) > uncharged.at({0, 1}));
    // A second copy of the group takes over the reads of the second walk, so the copy with the
    // deletion is charged, and two copies of the group win.
    Called groups = call_matrix(make, 2, &group_of);
    REQUIRE(groups.genotype == vector<int>({0, 1}));
    REQUIRE(groups.info->genotype_lls.at({0, 2}) < uncharged.at({0, 2}) - 10.0);
    // Two copies of the deletion hold no group of two or more walks, so they are not charged.
    REQUIRE(groups.info->genotype_lls.at({2, 2}) == uncharged.at({2, 2}));
    // GL and GQ come from the charged scores. GQ before its discounts is the phred-scaled gap
    // between the two best of them.
    require_scores(groups, group_scores(make(), 2, group_of, true));
    vector<double> lls;
    for (const auto& [genotype, ll] : groups.info->genotype_lls) {
        lls.push_back(ll);
    }
    std::ranges::sort(lls, std::greater<>());
    REQUIRE(groups.info->gq_undiscounted
            == Approx(logprob_to_phred(lls[1]) - logprob_to_phred(lls[0])));
}

TEST_CASE("A haploid copy of a group is charged against two copies of it",
          "[read_likelihood_caller][group_candidates]") {
    // Walks 0 and 1 form one group, and walk 2 is a group of its own. At ploidy 1, sharing lets
    // one copy of the group carry the reads of both its walks, which outnumber the reads of walk
    // 2. A diploid genotype of the group explains them with two walks, so the haploid copy is
    // charged, and walk 2 wins.
    const auto make = reads_fitting(3, {{{0}, 10}, {{1}, 10}, {{2}, 16}});
    const vector group_of{0, 0, 1};

    const auto uncharged = group_scores(make(), 1, group_of, false);
    REQUIRE(uncharged.at({0}) > uncharged.at({2}));
    Called groups = call_matrix(make, 1, &group_of);
    REQUIRE(groups.genotype == vector<int>({2}));
    require_scores(groups, group_scores(make(), 1, group_of, true));
}

TEST_CASE("Where a site's walks form one group, the walks alone choose its genotype",
          "[read_likelihood_caller][group_candidates]") {
    // Every genotype has the one multiset of groups, so the groups have nothing to choose between.
    const auto make = reads_fitting(3, {{{0}, 10}, {{1}, 10}, {{2}, 3}});
    const vector group_of{0, 0, 0};

    Called walks = call_matrix(make, 2, nullptr);
    Called groups = call_matrix(make, 2, &group_of);
    REQUIRE(groups.genotype == walks.genotype);
    REQUIRE(groups.info->genotype_lls == walks.info->genotype_lls);
}

TEST_CASE("A multiset of groups is charged from the fit of its reads, without the depth term",
          "[read_likelihood_caller][group_candidates]") {
    // Walks 0 and 1 form one group, and are 1,000 and 20,000 bp long. Walk 2 is a deletion, with a
    // group of its own. Many reads fit both walks of the group. With the depth term, the shared
    // score could take whichever walk of the group has the length that best fits the number of
    // reads, and the own score could not.
    const vector<size_t> lengths{1000, 20000, 0};
    const auto plain =
        reads_fitting(3, {{{0}, 40}, {{1}, 40}, {{2}, 30}, {{0, 1}, 100}}, lengths);
    const auto with_depth = [&plain, &lengths]() {
        AlleleReadLikelihoods matrix = plain();
        matrix.set_depth_context(lengths, 0.01, 100, 0.1);
        return matrix;
    };
    const vector group_of{0, 0, 1};
    REQUIRE(with_depth().uses_depth_term());

    // Return the charge on each genotype: its score before the charge less its score.
    const auto charges = [&](const function<AlleleReadLikelihoods()>& make) {
        const Called called = call_matrix(make, 2, &group_of);
        map<vector<int>, double> charge;
        for (const auto& [genotype, ll] : group_scores(make(), 2, group_of, false)) {
            charge[genotype] = ll - called.info->genotype_lls.at(genotype);
        }
        return charge;
    };
    const auto without = charges(plain);
    const auto with = charges(with_depth);
    REQUIRE(without.at({0, 2}) > 1.0);
    for (const auto& [genotype, charge] : without) {
        REQUIRE(with.at(genotype) == Approx(charge).margin(1e-6));
    }
    require_scores(call_matrix(with_depth, 2, &group_of),
                   group_scores(with_depth(), 2, group_of, true));
    // With the depth term, the gain of the group and the deletion would differ by more than a nat.
    const vector<vector<int>> group_and_deletion{{0, 2}, {1, 2}};
    REQUIRE(std::abs(gain_of(with_depth(), group_and_deletion, group_of, true)
                     - gain_of(with_depth(), group_and_deletion, group_of, false))
            > 1.0);
}

TEST_CASE("One more copy of a group is added to the best genotypes of the multiset",
          "[read_likelihood_caller][group_candidates]") {
    // Walks 0, 1 and 2 form group 0, and walk 3 is a group of its own. Reads fit walks 0 and 1,
    // walks 0 and 2, walk 1, walk 2, and walk 3. Walks 0, 1 and 2 are 300, 200 and 100 bp long.
    // Walk 0 is the longest, so a copy of group 0 has the most weight in the mixture as walk 0, and
    // the best genotypes by shared score have no ties.
    const auto make = reads_fitting(
        4, {{{0, 1}, 10}, {{0, 2}, 10}, {{1}, 6}, {{2}, 6}, {{3}, 10}}, {300, 200, 100, 100});
    const vector group_of{0, 0, 0, 1};
    const AlleleReadLikelihoods matrix = make();
    const auto gain = [&](const vector<vector<int>>& genotypes) {
        return gain_of(matrix, genotypes, group_of, false);
    };

    // Walks 0 and 3 are the best genotype of group 0 and walk 3 by both scores.
    for (const vector<int>* groups : {static_cast<const vector<int>*>(nullptr), &group_of}) {
        const auto scored = matrix.score_genotypes({{0, 3}, {1, 3}, {2, 3}}, groups);
        REQUIRE(scored[0].second > scored[1].second);
        REQUIRE(scored[0].second > scored[2].second);
    }
    // With one more copy of group 0, the best genotype by own score holds walks 1, 2 and 3. Adding
    // a walk to walks 0 and 3 cannot give that genotype, so the gain over the genotypes made by
    // adding a walk is larger than the gain over every genotype, and the charge is smaller.
    const double before = gain({{0, 3}, {1, 3}, {2, 3}});
    const double extended = gain({{0, 0, 3}, {0, 1, 3}, {0, 2, 3}});
    const double every = gain({{0, 0, 3}, {0, 1, 3}, {0, 2, 3}, {1, 1, 3}, {1, 2, 3}, {2, 2, 3}});
    REQUIRE(extended > every + 1.0);
    REQUIRE(before - extended > 1.0);

    // The charge on walks 0 and 3 is the gain before the extra copy less the gain over the
    // genotypes made by adding a walk.
    const Called groups = call_matrix(make, 2, &group_of);
    const auto uncharged = group_scores(make(), 2, group_of, false);
    REQUIRE(uncharged.at({0, 3}) - groups.info->genotype_lls.at({0, 3})
            == Approx(before - extended));
    require_scores(groups, group_scores(make(), 2, group_of, true));
}

TEST_CASE("Tied best genotypes of a multiset go to the first in VCF order",
          "[read_likelihood_caller][group_candidates]") {
    // Walks 0, 1 and 2 form group 0, and walk 3 is a group of its own. Reads fit walk 0, walk 1,
    // walk 2, walks 0 and 2, walks 1 and 2, and walk 3. At ploidy 1, every walk of group 0 has the
    // same shared score without the depth term, so walk 0 is the best genotype of group 0 by that
    // score. Walk 2 fits the most reads, so it is the best by own score. The walks have different
    // lengths, so their depth terms differ, and the tie must still go to walk 0.
    const vector<size_t> lengths{2000, 3000, 8500, 3000};
    const auto make = [&lengths]() {
        AlleleReadLikelihoods matrix = reads_fitting(
            4, {{{0}, 10}, {{1}, 10}, {{2}, 6}, {{0, 2}, 6}, {{1, 2}, 6}, {{3}, 30}})();
        matrix.set_depth_context(lengths, 0.01, 100, 0.1);
        return matrix;
    };
    const vector group_of{0, 0, 0, 1};
    const AlleleReadLikelihoods matrix = make();
    REQUIRE(matrix.uses_depth_term());
    const vector<vector<int>> group{{0}, {1}, {2}};
    const auto shared = matrix.score_genotypes(group, &group_of, false);
    const auto own = matrix.score_genotypes(group, nullptr, false);
    REQUIRE(shared[1].second == shared[0].second);
    REQUIRE(shared[2].second == shared[0].second);
    REQUIRE(own[2].second > own[0].second);
    REQUIRE(own[2].second > own[1].second);

    // With two copies of group 0, walks 0 and 1 are the best genotype by own score. Adding a walk
    // to walk 0 gives it, and adding a walk to walk 2 does not, so the charge on group 0 depends on
    // starting from walk 0 as well as from walk 2.
    const double before = gain_of(matrix, group, group_of, false);
    const double extended =
        gain_of(matrix, {{0, 0}, {0, 1}, {0, 2}, {1, 2}, {2, 2}}, group_of, false);
    const double from_own = gain_of(matrix, {{0, 2}, {1, 2}, {2, 2}}, group_of, false);
    REQUIRE(from_own > extended + 1.0);
    REQUIRE(before - extended > 1.0);

    const Called groups = call_matrix(make, 1, &group_of);
    const auto uncharged = group_scores(matrix, 1, group_of, false);
    REQUIRE(uncharged.at({2}) - groups.info->genotype_lls.at({2}) == Approx(before - extended));
    require_scores(groups, group_scores(matrix, 1, group_of, true));
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

TEST_CASE("A record with an unscored allele still gets QUAL from the reference column",
          "[read_likelihood_caller]") {
    // A written allele that was not scored has no column, and QUAL needs only the reference
    // allele's.
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

    // Record alleles: reference, the called alternative, and an allele that was not scored.
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
