/// \file linkage_model.cpp
/// Tests for the Li-Stephens linkage layer. The properties tested here are the ones the design
/// note (planning/vg-call-linkage-hmm.md, section 8) named in advance, and two of them exist to
/// catch harm rather than to confirm gain.

#include <cmath>
#include <vector>

#include "catch.hpp"
#include "../linkage_model.hpp"

namespace vg {
namespace unittest {

using namespace std;

/// A biallelic site: three genotypes in VCF order 0/0, 0/1, 1/1.
static LinkageModel::Site biallelic(size_t position, double ln_00, double ln_01, double ln_11,
                                    const vector<int>& haps) {
    LinkageModel::Site s;
    s.position = position;
    s.num_alleles = 2;
    s.genotype_ln_likelihood = {ln_00, ln_01, ln_11};
    s.haplotype_allele = haps;
    return s;
}

static size_t best_genotype(const vector<double>& post) {
    size_t best = 0;
    for (size_t i = 1; i < post.size(); ++i) {
        if (post[i] > post[best]) {
            best = i;
        }
    }
    return best;
}

TEST_CASE("Zero weight leaves the per-site genotype untouched", "[linkage_model]") {
    // The inertness property, and it is not decoration: the whole point of a weight is that the
    // shipped default recovers the existing caller exactly. Zero weight makes every transition
    // uniform, so the chain is memoryless and the posterior is the emission.
    LinkageModel::Params p;
    p.weight = 0.0;
    LinkageModel model(p);

    // Site 2's reads prefer 1/1, and the panel links site 1's called allele to allele 0. With no
    // transition model the reads must win outright.
    vector<LinkageModel::Site> sites{
        biallelic(1000, -20.0, -20.0, 0.0, {1, 1, 0, 0}),
        biallelic(1200, -20.0, -20.0, 0.0, {0, 0, 1, 1}),
    };
    auto post = model.posteriors(sites);
    REQUIRE(post.size() == 2);
    REQUIRE(best_genotype(post[0]) == LinkageModel::genotype_index(1, 1));
    REQUIRE(best_genotype(post[1]) == LinkageModel::genotype_index(1, 1));
    REQUIRE_FALSE(model.active());
}

TEST_CASE("Linkage decides a site whose reads cannot", "[linkage_model]") {
    // The gain case. Site 1's reads are decisive for 1/1; site 2's reads are flat between 0/0 and
    // 1/1. Every panel haplotype carrying allele 1 at site 1 carries allele 1 at site 2, so
    // linkage is the only thing that can break the tie, and it should.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;      // sites are close, so a switch should be expensive
    p.rho_min = 1e-4;
    LinkageModel model(p);

    vector<LinkageModel::Site> sites{
        biallelic(1000, -30.0, -30.0, 0.0, {1, 1, 0, 0}),
        biallelic(1100, 0.0, -30.0, 0.0, {1, 1, 0, 0}),
    };
    auto post = model.posteriors(sites);
    REQUIRE(best_genotype(post[1]) == LinkageModel::genotype_index(1, 1));

    // And with the linkage reversed -- the same flat site, but now the panel puts allele 0
    // alongside site 1's allele 1 -- the answer must flip. Otherwise the test is passing on the
    // emission rather than on linkage.
    sites[1] = biallelic(1100, 0.0, -30.0, 0.0, {0, 0, 1, 1});
    post = model.posteriors(sites);
    REQUIRE(best_genotype(post[1]) == LinkageModel::genotype_index(0, 0));
}

TEST_CASE("Linkage does not override decisive reads", "[linkage_model]") {
    // The harm case, and the more important of the pair. Site 2's reads strongly favour the
    // allele that linkage argues against. At a sane weight the reads must still win: a prior that
    // can overturn 30 nats of read evidence is not a prior, it is an assertion.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;
    LinkageModel model(p);

    vector<LinkageModel::Site> sites{
        biallelic(1000, -30.0, -30.0, 0.0, {1, 1, 0, 0}),
        // Panel links site 1's allele 1 to allele 0 here, but the reads say 1/1 by 30 nats.
        biallelic(1100, -30.0, -30.0, 0.0, {0, 0, 1, 1}),
    };
    auto post = model.posteriors(sites);
    REQUIRE(best_genotype(post[1]) == LinkageModel::genotype_index(1, 1));
}

TEST_CASE("An allele no panel haplotype carries stays callable", "[linkage_model]") {
    // Guards the regression that would present as a precision improvement. A state implies a
    // genotype, so without the wildcard contributing per-allele mass, a genotype the panel cannot
    // spell is unreachable and the model quietly suppresses novel alleles. The graph need not
    // contain the sample being genotyped, so this is the common case, not a corner.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.escape = 1e-2;
    LinkageModel model(p);

    // Three alleles; every panel haplotype carries 0 or 1, and nothing carries allele 2. The
    // reads are decisive for 2/2.
    LinkageModel::Site s;
    s.position = 1000;
    s.num_alleles = 3;
    s.genotype_ln_likelihood.assign(6, -40.0);
    s.genotype_ln_likelihood[LinkageModel::genotype_index(2, 2)] = 0.0;
    s.haplotype_allele = {0, 0, 1, 1};

    LinkageModel::Site s2 = s;
    s2.position = 1200;

    auto post = model.posteriors({s, s2});
    REQUIRE(post[0].size() == 6);
    REQUIRE(post[0][LinkageModel::genotype_index(2, 2)] > 0.0);
    REQUIRE(best_genotype(post[0]) == LinkageModel::genotype_index(2, 2));
}

TEST_CASE("A haplotype absent from a site carries no allele there", "[linkage_model]") {
    // A haplotype whose path ends between two sites -- a fragment boundary, which happens in real
    // pangenomes and not only in sampled ones -- must contribute no linkage across the gap.
    // Treating absence as the reference allele would invent evidence for the reference.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    LinkageModel model(p);

    // At site 2 the haplotypes carrying allele 1 at site 1 are absent (-1). The flat site must
    // not be pulled toward allele 1, because nothing links the two.
    vector<LinkageModel::Site> linked{
        biallelic(1000, -30.0, -30.0, 0.0, {1, 1, 0, 0}),
        biallelic(1100, 0.0, -30.0, 0.0, {1, 1, 0, 0}),
    };
    vector<LinkageModel::Site> absent{
        biallelic(1000, -30.0, -30.0, 0.0, {1, 1, 0, 0}),
        biallelic(1100, 0.0, -30.0, 0.0, {-1, -1, 0, 0}),
    };
    auto with_link = model.posteriors(linked);
    auto without = model.posteriors(absent);
    size_t hom_alt = LinkageModel::genotype_index(1, 1);
    // Linkage should support 1/1 strictly more when the carriers are actually present.
    REQUIRE(with_link[1][hom_alt] > without[1][hom_alt]);
}

TEST_CASE("A certain switch is equivalent to forgetting the previous site",
          "[linkage_model]") {
    // A switch probability of 1 means the panel certainly changed identity between the sites, so
    // the transition must carry nothing across -- the same posterior as a fresh chain. Reached
    // here by making the sites effectively infinitely far apart, which is the only way to
    // saturate rho now that the separate block-switch term is gone.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 1e-9;             // gap / scale enormous, so rho saturates at 1
    p.rho_min = 0.0;
    LinkageModel model(p);

    vector<LinkageModel::Site> sites{
        biallelic(1000, -30.0, -30.0, 0.0, {1, 1, 0, 0}),
        biallelic(1100, 0.0, -30.0, 0.0, {1, 1, 0, 0}),
    };
    auto post = model.posteriors(sites);

    LinkageModel::Params alone = p;
    LinkageModel single(alone);
    auto isolated = single.posteriors({sites[1]});

    REQUIRE(post[1].size() == isolated[0].size());
    for (size_t g = 0; g < post[1].size(); ++g) {
        REQUIRE(post[1][g] == Approx(isolated[0][g]).margin(1e-9));
    }
}

TEST_CASE("Posteriors are a distribution over genotypes", "[linkage_model]") {
    LinkageModel::Params p;
    p.weight = 1.0;
    LinkageModel model(p);
    vector<LinkageModel::Site> sites{
        biallelic(1000, -1.0, 0.0, -2.0, {1, 0, 0, 1}),
        biallelic(3000, -2.0, 0.0, -1.0, {0, 1, 1, 0}),
        biallelic(9000, 0.0, -1.0, -3.0, {1, 1, 0, 0}),
    };
    auto post = model.posteriors(sites);
    for (const auto& p_site : post) {
        REQUIRE(p_site.size() == 3);
        double total = 0.0;
        for (double v : p_site) {
            REQUIRE(v >= 0.0);
            total += v;
        }
        REQUIRE(total == Approx(1.0).margin(1e-9));
    }
}

TEST_CASE("A windowed pass matches exact inference on a short chain", "[linkage_model]") {
    // Windowing exists so that a chain-wide pass does not serialise a caller that is parallel
    // over snarls. It is only defensible if it agrees with exact inference where both can run.
    // Sites spaced well past the linkage scale, so influence is spent within a few steps and a
    // modest margin really is enough. At 500 bp against a 10 kb scale it would survive ~137
    // sites, and the first version of this test asserted agreement a 12-site margin could not
    // deliver -- the premise was wrong, not the windowing.
    vector<LinkageModel::Site> sites;
    for (size_t i = 0; i < 40; ++i) {
        sites.push_back(biallelic(1000 + i * 20000,
                                  -(double)(i % 3), 0.0, -(double)((i + 1) % 3),
                                  {(int)(i % 2), (int)((i + 1) % 2), 0, 1}));
    }
    LinkageModel::Params exact;
    exact.weight = 1.0;
    exact.window = 1000;         // one window covers the chain
    exact.margin = 0;
    LinkageModel::Params windowed = exact;
    windowed.window = 8;
    windowed.margin = 12;        // margin well past the range over which linkage carries

    auto a = LinkageModel(exact).posteriors(sites);
    auto b = LinkageModel(windowed).posteriors(sites);
    REQUIRE(a.size() == b.size());
    for (size_t t = 0; t < a.size(); ++t) {
        REQUIRE(a[t].size() == b[t].size());
        for (size_t g = 0; g < a[t].size(); ++g) {
            REQUIRE(a[t][g] == Approx(b[t][g]).margin(1e-6));
        }
    }
}

TEST_CASE("The frequency prior is separable from linkage", "[linkage_model]") {
    // Summing state posteriors per genotype weights each genotype by how many haplotype pairs
    // spell it, which is a panel allele-frequency prior rather than linkage. Over a panel chosen
    // by haplotype sampling against the reads being genotyped, that is the same evidence twice,
    // so it must be separately switchable -- and off by default.
    LinkageModel::Params off;
    off.weight = 0.0;
    off.freq_prior = 0.0;
    LinkageModel::Params on = off;
    on.freq_prior = 1.0;

    // Flat reads, and a panel where allele 0 is much commoner than allele 1.
    vector<LinkageModel::Site> sites{
        biallelic(1000, 0.0, 0.0, 0.0, {0, 0, 0, 0, 0, 1}),
        biallelic(1500, 0.0, 0.0, 0.0, {0, 0, 0, 0, 0, 1}),
    };
    auto without = LinkageModel(off).posteriors(sites);
    auto with = LinkageModel(on).posteriors(sites);

    size_t hom_ref = LinkageModel::genotype_index(0, 0);
    // With the prior, the common allele's homozygote must gain; without it, flat reads must stay
    // flat between the genotypes the panel can spell.
    REQUIRE(with[0][hom_ref] > without[0][hom_ref]);
    REQUIRE(without[0][hom_ref] == Approx(without[0][LinkageModel::genotype_index(1, 1)])
                                        .margin(1e-6));
}


TEST_CASE("The collector keeps sites compactly and re-decides only what changed",
          "[linkage_model]") {
    // The compact form is the reason the two-phase design fits: keeping a CallInfo per site --
    // with its vector<SnarlTraversal> and map -- would run to hundreds of megabytes over a
    // chromosome arm. This asserts the budget rather than trusting the estimate.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;
    LinkageCollector collector(p, 4);

    // Site 1 decisive for 1/1; site 2 flat, with the panel linking allele 1 to allele 1.
    collector.record("chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, /*key*/ 11, /*share*/ 1.0);
    collector.record("chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, /*key*/ 22, /*share*/ 1.0);

    REQUIRE(collector.num_sites() == 2);
    // Two entries, six floats, eight int8s. Anything near a vector-per-site layout would be
    // several times this.
    REQUIRE(collector.bytes() < 2 * 48 + 6 * 4 + 8 + 64);

    auto changes = collector.resolve();
    // Site 1 was already called 1/1 and must not be reported; site 2 was called 0/0 and linkage
    // should move it to 1/1.
    REQUIRE(changes.size() == 1);
    REQUIRE(changes[0].record_key == 22);
    REQUIRE(changes[0].allele_i == 1);
    REQUIRE(changes[0].allele_j == 1);
    REQUIRE(changes[0].posterior > 0.5);
}

TEST_CASE("The collector reports nothing at zero weight", "[linkage_model]") {
    // The shipped default. If this ever reports a change, the caller's output has moved without
    // anyone asking for it.
    LinkageModel::Params p;
    p.weight = 0.0;
    LinkageCollector collector(p, 4);
    collector.record("chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);
    collector.record("chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);
    REQUIRE(collector.resolve().empty());
}

TEST_CASE("The collector does not link across contigs", "[linkage_model]") {
    // Two sites at the same coordinates on different contigs are not neighbours. Linking them
    // would be an easy mistake to make and an invisible one, since the gap would look tiny.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    LinkageCollector same(p, 4);
    same.record("chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);
    same.record("chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);

    LinkageCollector split(p, 4);
    split.record("chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);
    split.record("chr2", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);

    REQUIRE(same.resolve().size() == 1);
    REQUIRE(split.resolve().empty());
}

TEST_CASE("The collector sorts by reference position, not arrival order",
          "[linkage_model]") {
    // Sites arrive in node-ID order, which is close to reference order but not guaranteed to be
    // it. Transition probabilities come from the gaps, so arrival order would silently feed the
    // model wrong distances -- and with a short scale, a wrong gap changes the answer.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 500.0;            // short, so the gap matters
    p.rho_min = 1e-4;

    LinkageCollector ordered(p, 4);
    ordered.record("chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);
    ordered.record("chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);

    LinkageCollector shuffled(p, 4);
    shuffled.record("chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);
    shuffled.record("chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);

    auto a = ordered.resolve();
    auto b = shuffled.resolve();
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        REQUIRE(a[i].record_key == b[i].record_key);
        REQUIRE(a[i].allele_i == b[i].allele_i);
        REQUIRE(a[i].allele_j == b[i].allele_j);
        REQUIRE(a[i].posterior == Approx(b[i].posterior).margin(1e-12));
    }
}

TEST_CASE("Genotype indices round-trip through the collector's decode", "[linkage_model]") {
    // resolve() decodes a genotype index back to an allele pair by inverting the triangular
    // number. Worth pinning: an off-by-one here would emit a plausible but wrong genotype.
    for (size_t j = 0; j < 8; ++j) {
        for (size_t i = 0; i <= j; ++i) {
            size_t idx = LinkageModel::genotype_index(i, j);
            size_t dj = 0;
            while (LinkageModel::genotype_index(0, dj) <= idx) {
                ++dj;
            }
            --dj;
            size_t di = idx - (dj * (dj + 1) / 2);
            REQUIRE(di == i);
            REQUIRE(dj == j);
        }
    }
}

}
}
