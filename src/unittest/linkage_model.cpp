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

/// These tests were written against the collector's earlier argument shape: a dense likelihood vector
/// in genotype-index order and a panel in VCF allele numbering. The collector now takes the
/// genotyper's own space -- likelihoods keyed by candidate traversal pairs, the panel as traversal
/// indices -- because symbolic collapsing maps distinct traversals onto one VCF allele, and a parent
/// whose haplotypes differ only inside its child chains is homozygous in the collapsed numbering and
/// heterozygous in the real one.
///
/// The tests are kept in the old shape on purpose: they exercise the model and the arenas, not the
/// numbering, and the identity case -- traversal index equals VCF allele -- is where the two agree,
/// so translating them would add noise without adding coverage. Construction of the compact space
/// itself is covered separately, against `compact_allele_space`.
static void record_dense(LinkageCollector& c, const string& contig, size_t position,
                         size_t num_alleles, const vector<double>& dense_gls,
                         const vector<int>& panel, size_t called_i, size_t called_j,
                         size_t record_key, double share, size_t ploidy = 2,
                         int64_t start_node = 0, int64_t end_node = 0,
                         bool nested = false, size_t parent_record_key = 0, uint64_t parent_crossing = 0, size_t generation = 0,
                         bool emitted = true) {
    map<vector<int>, double> gls;
    if (ploidy == 1) {
        for (size_t a = 0; a < num_alleles && a < dense_gls.size(); ++a) {
            gls[vector<int>{(int)a}] = dense_gls[a];
        }
    } else {
        for (size_t j = 0; j < num_alleles; ++j) {
            for (size_t i = 0; i <= j; ++i) {
                size_t g = LinkageModel::genotype_index(i, j);
                if (g < dense_gls.size()) {
                    gls[vector<int>{(int)i, (int)j}] = dense_gls[g];
                }
            }
        }
    }
    vector<int> ident(num_alleles);
    for (size_t i = 0; i < num_alleles; ++i) {
        ident[i] = (int)i;
    }
    c.record(contig, position, gls, panel, (int)called_i, (int)called_j, ident, record_key,
             share, ploidy, start_node, end_node,
             LinkageCollector::SiteContext{
                 .nested = nested,
                 .parent_record_key = parent_record_key,
                 .parent_crossing = parent_crossing,
                 .generation = generation,
                 .emitted = emitted,
             });
}

static bool respecify_dense(LinkageCollector& c, size_t record_key,
                            const string& contig, size_t position, size_t num_alleles,
                            const vector<double>& dense_gls, const vector<int>& panel,
                            size_t called_i, size_t called_j, size_t ploidy,
                            bool nested, uint64_t parent_crossing,
                            bool emitted = true) {
    map<vector<int>, double> gls;
    if (ploidy == 1) {
        for (size_t a = 0; a < num_alleles && a < dense_gls.size(); ++a) {
            gls[vector<int>{(int)a}] = dense_gls[a];
        }
    } else {
        for (size_t j = 0; j < num_alleles; ++j) {
            for (size_t i = 0; i <= j; ++i) {
                size_t g = LinkageModel::genotype_index(i, j);
                if (g < dense_gls.size()) {
                    gls[vector<int>{(int)i, (int)j}] = dense_gls[g];
                }
            }
        }
    }
    vector<int> ident(num_alleles);
    for (size_t i = 0; i < num_alleles; ++i) {
        ident[i] = (int)i;
    }
    return c.respecify(record_key, contig, position, gls, panel, (int)called_i, (int)called_j,
                       ident, ploidy, nested, parent_crossing, emitted);
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
    record_dense(collector, "chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, /*key*/ 11, /*share*/ 1.0);
    record_dense(collector, "chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, /*key*/ 22, /*share*/ 1.0);

    REQUIRE(collector.num_sites() == 2);
    // Two entries, six floats, eight int8s, all in flat arenas, plus the record-key index at two
    // hash nodes a site. The index is counted here because `bytes()` counts it: it is real memory,
    // and a reported figure that omitted it would drift from the measured one by a third.
    //
    // A bound on per-site overhead rather than the exact figure, which is 128 bytes an entry and 48
    // of index today: an exact assertion fails whenever a field is added, which is the wrong reason
    // to fail. What this catches is a vector-per-site layout -- three vector headers, 72 bytes, plus
    // three heap allocations for every site -- which would put this at 496 and over the bound.
    // Raised from +64 to +128 when the traversal-derived ordering added four fields to Entry --
    // align_rank, chain_index, chain_backward, unpositioned -- taking this from 400 to 428. That is
    // about 12 bytes a site, 2.6 MB over chr20's 220k entries, against a 4.5 GB peak.
    //
    // The bound still does its job, which the comment above states: a vector-per-site layout adds
    // 72 bytes a site and would put this at 572, far over 480. What it must NOT become is an exact
    // assertion that fails on every added field.
    REQUIRE(collector.bytes() < 2 * (128 + 48) + 128);

    const size_t moved = collector.resolve();
    // Site 1 was already called 1/1 and must not be counted; site 2 was called 0/0 and linkage
    // should move it to 1/1.
    //
    // Asserted on the settled genotype rather than on a patch's contents. The record is built from
    // the settled pair now, so that pair is the observable and the patch it used to describe does
    // not exist.
    REQUIRE(moved == 1);
    int a = -1, b = -1;
    size_t settled_ploidy = 0;
    REQUIRE(collector.settled_traversals(22, &a, &b, &settled_ploidy));
    REQUIRE(settled_ploidy == 2);
    REQUIRE(a == 1);
    REQUIRE(b == 1);
    // And the site that was already right stayed where it was.
    REQUIRE(collector.settled_traversals(11, &a, &b, &settled_ploidy));
    REQUIRE(a == 1);
    REQUIRE(b == 1);
}

TEST_CASE("The collector reports nothing at zero weight", "[linkage_model]") {
    // The shipped default. If this ever reports a change, the caller's output has moved without
    // anyone asking for it.
    LinkageModel::Params p;
    p.weight = 0.0;
    LinkageCollector collector(p, 4);
    record_dense(collector, "chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);
    record_dense(collector, "chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);
    REQUIRE(collector.resolve() == 0);
}

TEST_CASE("The collector does not link across contigs", "[linkage_model]") {
    // Two sites at the same coordinates on different contigs are not neighbours. Linking them
    // would be an easy mistake to make and an invisible one, since the gap would look tiny.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    LinkageCollector same(p, 4);
    record_dense(same, "chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);
    record_dense(same, "chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);

    LinkageCollector split(p, 4);
    record_dense(split, "chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);
    record_dense(split, "chr2", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);

    REQUIRE(same.resolve() == 1);
    REQUIRE(split.resolve() == 0);
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
    record_dense(ordered, "chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);
    record_dense(ordered, "chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);

    LinkageCollector shuffled(p, 4);
    record_dense(shuffled, "chr1", 1100, 2, {0.0, -30.0, 0.0}, {1, 1, 0, 0}, 0, 0, 22, /*share*/ 1.0);
    record_dense(shuffled, "chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1, 0, 0}, 1, 1, 11, /*share*/ 1.0);

    // Order-independence asserted on the settled genotypes, which is what the caller reads, rather
    // than on the order of a patch list that no longer exists.
    REQUIRE(ordered.resolve() == shuffled.resolve());
    for (size_t key : {(size_t)11, (size_t)22}) {
        int ai = -1, aj = -1, bi = -1, bj = -1;
        size_t ap = 0, bp = 0;
        REQUIRE(ordered.settled_traversals(key, &ai, &aj, &ap));
        REQUIRE(shuffled.settled_traversals(key, &bi, &bj, &bp));
        REQUIRE(ai == bi);
        REQUIRE(aj == bj);
        REQUIRE(ap == bp);
    }
}

TEST_CASE("respecify moves a site to the ploidy its settled parent implies", "[linkage_model]") {
    // The coherence guarantee in one function. Descent records a child at the ploidy its parent's
    // *pre-linkage* genotype implied; the barrier then learns what the settled genotype implies and
    // moves the entry before its own generation resolves. Without this the child keeps a ploidy its
    // own parent contradicts, which is what the nested_diploid FILTER used to label.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    const size_t PARENT = 9, CHILD = 91;
    LinkageCollector collector(p, 2);
    record_dense(collector, "chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1}, 1, 1, PARENT,
                     1.0, 2, 10, 20);
    // Recorded haploid, as descent would: one allele, one likelihood per allele.
    record_dense(collector, "chr1", 1010, 2, {0.0, -30.0}, {0, 0}, 0, 0, CHILD,
                     1.0, 1, 11, 12, true, PARENT, /*crossing*/ (uint64_t)1 << 1);

    // The settled parent turns out to cross on both haplotypes, so the child is diploid. Its
    // likelihoods are the triangular vector now, and it is no longer a nested haploid site.
    REQUIRE(respecify_dense(collector, CHILD, "chr1", 1010, 2, {-30.0, -30.0, 0.0}, {1, 1}, 1, 1, 2,
                                /*nested*/ false, 0, 0));
    // An unknown key must say so rather than silently doing nothing.
    REQUIRE_FALSE(respecify_dense(collector, 12345, "chr1", 1010, 2, {0.0, -30.0, -30.0}, {0, 0},
                                  0, 0, 1, true, 0, 0));

    vector<LinkageCollector::PhaseCall> phased;
    collector.resolve(&phased);

    // Nothing to flag: the ploidy the record carries is the one its parent implies.
        const LinkageCollector::PhaseCall* child = nullptr;
    for (const auto& pc : phased) {
        if (pc.record_key == CHILD) {
            child = &pc;
        }
    }
    REQUIRE(child != nullptr);
    REQUIRE(child->ploidy == 2);
    // Asserted on the traversal, not the compact allele index. The collector's space is the called
    // pair plus the panel-carried traversals, so a site whose panel names one traversal has a
    // one-element space and compact 0 *is* traversal 1. The genome fact -- which traversal each
    // strand is on -- is what this test is about, and it is what trav_* carries.
    REQUIRE(child->trav_first == 1);
    REQUIRE(child->trav_second == 1);
}

TEST_CASE("retract drops a site the settled parent does not carry", "[linkage_model]") {
    // The other half. Where the settled genotype crosses the chain on neither haplotype the sample
    // has no copy of it, so there is no site: the entry leaves the chains, the phasing and every
    // count. Marked rather than erased, because the arenas are flat and every other entry holds
    // offsets into them -- so the neighbours must survive it untouched.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    const size_t PARENT = 9, CHILD = 91, SIBLING = 92;
    LinkageCollector collector(p, 2);
    record_dense(collector, "chr1", 1000, 2, {-30.0, -30.0, 0.0}, {1, 1}, 1, 1, PARENT,
                     1.0, 2, 10, 20);
    record_dense(collector, "chr1", 1010, 2, {0.0, -30.0}, {0, 0}, 0, 0, CHILD,
                     1.0, 1, 11, 12, true, PARENT, /*crossing*/ (uint64_t)1 << 0);
    // A neighbour recorded *after* the retracted one, so its arena offsets sit beyond the hole.
    record_dense(collector, "chr1", 1020, 2, {-30.0, -30.0, 0.0}, {1, 1}, 1, 1, SIBLING,
                     1.0, 2, 13, 14);

    size_t before = collector.num_sites_at(0);
    REQUIRE(collector.retract(CHILD));
    REQUIRE_FALSE(collector.retract(CHILD));        // already gone
    REQUIRE_FALSE(collector.retract(999));          // never existed

    vector<LinkageCollector::PhaseCall> phased;
    collector.resolve(&phased);

    for (const auto& pc : phased) {
        REQUIRE(pc.record_key != CHILD);
    }
    // The neighbour past the hole is still read correctly, which is the thing marking rather than
    // erasing is for.
    const LinkageCollector::PhaseCall* sibling = nullptr;
    for (const auto& pc : phased) {
        if (pc.record_key == SIBLING) {
            sibling = &pc;
        }
    }
    REQUIRE(sibling != nullptr);
    REQUIRE(sibling->position == 1020);
    REQUIRE(sibling->trav_first == 1);
    REQUIRE(collector.num_sites_at(0) < before);
}

TEST_CASE("A nested site takes its strand from the parent traversal that carries it",
          "[linkage_model]") {
    // The single derivation, and the reason it replaced three. A nested chain is carried by exactly
    // one of its parent's settled traversals; that it has one copy and that it sits on that
    // traversal's strand are the same statement, so nothing can disagree about them.
    //
    // What this pins is the identity match. The collector sorts the called pair and the Viterbi then
    // orients it against the panel, so the *index* a traversal sat at during descent means nothing
    // by the time the child is placed -- which is why the traversal itself is recorded. Both strands
    // are exercised, because getting this backwards puts every child on the wrong haplotype and
    // still produces perfectly well-formed output.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    const size_t PARENT = 9, CHILD = 91;

    // Which parent traversal the child hangs off -- stated as the CROSSING MASK, which is what
    // the placement derives from -- and whether that traversal is in the settled pair.
    struct Case { int carrying_trav; bool placed; };
    const Case cases[] = {{0, true}, {1, true}, {7, false}};
    for (const Case& c : cases) {
        LinkageCollector collector(p, 2);
        // A het parent, decisively 0/1, with one panel haplotype on each allele so the Viterbi has
        // a strand to put each on.
        record_dense(collector, "chr1", 1000, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1, PARENT,
                         /*share*/ 1.0, /*ploidy*/ 2, /*start*/ 10, /*end*/ 20);
        record_dense(collector, "chr1", 1010, 2, {0.0, -30.0}, {0, 0}, 0, 0, CHILD,
                         /*share*/ 1.0, /*ploidy*/ 1, /*start*/ 11, /*end*/ 12,
                         /*nested*/ true, PARENT,
                         /*crossing*/ (uint64_t)1 << c.carrying_trav);

        vector<LinkageCollector::PhaseCall> phased;
        collector.resolve(&phased);

        const LinkageCollector::PhaseCall* child = nullptr;
        const LinkageCollector::PhaseCall* parent = nullptr;
        for (const auto& pc : phased) {
            if (pc.record_key == CHILD) {
                child = &pc;
            } else if (pc.record_key == PARENT) {
                parent = &pc;
            }
        }
        REQUIRE(child != nullptr);
        REQUIRE(parent != nullptr);
        // In the parent's block either way: the record is at that locus whether or not a strand
        // could be named for it, and starting a block of its own would fragment the phasing.
        REQUIRE(child->phase_set == parent->phase_set);

        if (!c.placed) {
            // The parent settled on a pair that does not contain this traversal, so there is no
            // haplotype to name. Claiming one would put a variant in the emitted genome that the
            // parent record does not carry.
            REQUIRE(child->nested_strand == -1);
            REQUIRE(child->hap_first == LinkageModel::WILDCARD);
            REQUIRE(child->hap_second == LinkageModel::WILDCARD);
            continue;
        }
        // Placed on the strand the carrying traversal was phased onto -- found by asking the
        // parent's own phased pair, not by trusting an index recorded earlier.
        REQUIRE(child->nested_strand >= 0);
        const int want = parent->trav_first == c.carrying_trav ? 0 : 1;
        REQUIRE(child->nested_strand == want);
        if (want == 0) {
            REQUIRE(child->hap_first == parent->hap_first);
            REQUIRE(child->hap_second == LinkageModel::WILDCARD);
        } else {
            REQUIRE(child->hap_second == parent->hap_second);
            REQUIRE(child->hap_first == LinkageModel::WILDCARD);
        }
    }
}

TEST_CASE("The phasing comes back in reference order even with nested sites in it",
          "[linkage_model]") {
    // The mosaic output reads the phasing as one ordered sweep per contig and closes a segment only
    // where the haplotype changes, so order is a contract and not a convenience.
    //
    // Nested sites break it by construction: placing one needs its parent already phased, so they
    // are appended after every chain. Unsorted, a nested site early on a contig shares a segment
    // with one far along it -- chr20 emitted five segments spanning tens of megabases, one claiming
    // 284 sites between ref_start 451,374 and ref_end 65,512,343. The site totals still added up,
    // which is what the harness checks, so nothing looked wrong.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    LinkageCollector collector(p, 2);
    // Two ordinary sites far apart, and a nested child of the first that belongs between them.
    record_dense(collector, "chr1", 1000, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1, /*key*/ 1,
                     /*share*/ 1.0, /*ploidy*/ 2, /*start*/ 10, /*end*/ 20);
    record_dense(collector, "chr1", 900000, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1, /*key*/ 2,
                     /*share*/ 1.0, /*ploidy*/ 2, /*start*/ 90, /*end*/ 100);
    record_dense(collector, "chr1", 1010, 2, {0.0, -30.0}, {0, 0}, 0, 0, /*key*/ 3,
                     /*share*/ 1.0, /*ploidy*/ 1, /*start*/ 11, /*end*/ 12,
                     /*nested*/ true, /*parent*/ 1, /*crossing*/ 1);

    vector<LinkageCollector::PhaseCall> phased;
    collector.resolve(&phased);

    REQUIRE(phased.size() == 3);
    for (size_t i = 1; i < phased.size(); ++i) {
        REQUIRE(phased[i - 1].contig <= phased[i].contig);
        if (phased[i - 1].contig == phased[i].contig) {
            REQUIRE(phased[i - 1].position < phased[i].position);
        }
    }
    // Specifically: the nested site sits between its parent and the far site, not after both.
    REQUIRE(phased[1].record_key == 3);
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

//------------------------------------------------------------------------------
// phasing()

/// Every site free, sized to the chain.
static vector<size_t> unconstrained(size_t n) {
    return vector<size_t>(n, LinkageModel::NO_CONSTRAINT);
}

static size_t count_switches(const vector<LinkageModel::Phase>& ph) {
    size_t n = 0;
    for (size_t t = 1; t < ph.size(); ++t) {
        n += (ph[t].first != ph[t - 1].first) + (ph[t].second != ph[t - 1].second);
    }
    return n;
}

TEST_CASE("Phasing recovers a planted mosaic", "[linkage_model]") {
    // Haplotypes 0 and 1 spell what the reads want over the first half, 2 and 3 over the second,
    // so the only path explaining the whole chain switches once per strand in the middle. This is
    // the property the mosaic output rests on and it cannot be checked against posteriors(),
    // which never forms a path at all.
    LinkageModel::Params p;
    p.weight = 2.0;
    LinkageModel model(p);

    vector<LinkageModel::Site> sites;
    for (size_t i = 0; i < 6; ++i) {
        vector<int> haps = (i < 3) ? vector<int>{0, 1, 0, 0} : vector<int>{0, 0, 0, 1};
        sites.push_back(biallelic(1000 + i * 100, -20.0, 0.0, -20.0, haps));
    }
    auto ph = model.phasing(sites, unconstrained(sites.size()));
    REQUIRE(ph.size() == sites.size());

    for (size_t t = 0; t < sites.size(); ++t) {
        const vector<int>& h = sites[t].haplotype_allele;
        int a = ph[t].first == LinkageModel::WILDCARD ? -1 : h[ph[t].first];
        int b = ph[t].second == LinkageModel::WILDCARD ? -1 : h[ph[t].second];
        REQUIRE(((a == 0 && b == 1) || (a == 1 && b == 0) || a < 0 || b < 0));
    }
    // One switch per strand, not one per site: the mosaic is piecewise, which is the whole basis
    // for run-length encoding the output.
    REQUIRE(count_switches(ph) <= 2);
}

TEST_CASE("Constrained phasing spells the required genotype everywhere", "[linkage_model]") {
    // The consistency guarantee the mosaic output exists to provide. If this can fail, the
    // emitted genome and the emitted VCF can disagree, which is the one thing it must not do.
    LinkageModel::Params p;
    p.weight = 2.0;
    LinkageModel model(p);

    // Haplotype 1 carries allele 1 throughout and haplotype 2 alternates, so the pair (1,2)
    // spells the constraint at every site with no switch at all. Every site's *reads* are set to
    // prefer 0/0, so only the constraint can be producing the answer.
    //
    // The constraints have to be jointly explicable by some pair, which is not a weakness of the
    // test but the model working: constraints that flip genotype faster than any panel pair can
    // follow are better explained by the wildcard than by paying a switch per site, and the model
    // will correctly say so. Real calls come from the panel, so they do not look like that.
    size_t n = 6;
    vector<LinkageModel::Site> sites;
    vector<size_t> want;
    for (size_t i = 0; i < n; ++i) {
        vector<int> haps{0, 1, (int)(i % 2), 0};
        sites.push_back(biallelic(1000 + i * 100, 0.0, -20.0, -20.0, haps));
        want.push_back(LinkageModel::genotype_index(1, i % 2));
    }
    auto ph = model.phasing(sites, want);
    REQUIRE(ph.size() == n);
    for (size_t t = 0; t < n; ++t) {
        const vector<int>& h = sites[t].haplotype_allele;
        REQUIRE(ph[t].first != LinkageModel::WILDCARD);
        REQUIRE(ph[t].second != LinkageModel::WILDCARD);
        size_t got = LinkageModel::genotype_index((size_t)h[ph[t].first],
                                                  (size_t)h[ph[t].second]);
        REQUIRE(got == want[t]);
    }
    // And it should find the pair that needs no switching, not merely *a* consistent pair.
    REQUIRE(count_switches(ph) == 0);
}

TEST_CASE("A constraint no panel pair can follow routes through the wildcard",
          "[linkage_model]") {
    // The other side of the same behaviour, pinned because it looks like a bug when first seen.
    // Here the constraint flips 1/1, 0/0, 0/1 at 100 bp spacing. A panel explanation would have to
    // switch both strands twice, at about 10.6 nats a switch, against 4.6 nats per free strand
    // for the wildcard -- so "explained by nothing in the panel" is genuinely the better answer
    // and the model returns it rather than forcing an implausible mosaic.
    LinkageModel::Params p;
    p.weight = 2.0;
    LinkageModel model(p);
    vector<LinkageModel::Site> sites{
        biallelic(1000, 0.0, -20.0, -20.0, {0, 1, 0, 1}),
        biallelic(1100, -20.0, 0.0, -20.0, {0, 1, 0, 1}),
        biallelic(1200, -20.0, -20.0, 0.0, {0, 1, 0, 1}),
    };
    vector<size_t> want{
        LinkageModel::genotype_index(1, 1),
        LinkageModel::genotype_index(0, 0),
        LinkageModel::genotype_index(0, 1),
    };
    auto ph = model.phasing(sites, want);
    REQUIRE(ph.size() == 3);
    bool any_wildcard = false;
    for (const auto& e : ph) {
        any_wildcard |= (e.first == LinkageModel::WILDCARD
                         || e.second == LinkageModel::WILDCARD);
    }
    REQUIRE(any_wildcard);
}

TEST_CASE("Phasing stays feasible where the panel cannot spell the call", "[linkage_model]") {
    // Panel enumeration makes this unreachable, but --enumerate-support does not, and returning
    // nothing there would drop a whole chain over one site. The wildcard is what keeps the
    // constrained problem solvable.
    LinkageModel::Params p;
    p.weight = 2.0;
    LinkageModel model(p);

    vector<LinkageModel::Site> sites{
        biallelic(1000, 0.0, -20.0, -20.0, {0, 0, 0, 0}),
        biallelic(1100, 0.0, -20.0, -20.0, {0, 0, 0, 0}),
    };
    vector<size_t> want{LinkageModel::genotype_index(1, 1),
                        LinkageModel::genotype_index(1, 1)};
    auto ph = model.phasing(sites, want);
    REQUIRE(ph.size() == 2);
    REQUIRE(ph[0].first == LinkageModel::WILDCARD);
    REQUIRE(ph[0].second == LinkageModel::WILDCARD);
}

TEST_CASE("Window seams do not manufacture switches", "[linkage_model]") {
    // The failure guarded against here would look like a result rather than a bug: decoded
    // independently, two windows choose the state at their shared seam twice, so the join shows a
    // switch at every window boundary -- and that count scales with the site count, which is
    // exactly the shape a real biological signal would have.
    //
    // A chain one haplotype pair explains throughout must phase to zero switches whatever the
    // window size, so run it with a window far shorter than the chain to force seams.
    LinkageModel::Params p;
    p.weight = 2.0;
    p.window = 5;
    p.margin = 2;
    LinkageModel model(p);

    vector<LinkageModel::Site> sites;
    for (size_t i = 0; i < 40; ++i) {
        sites.push_back(biallelic(1000 + i * 100, -20.0, 0.0, -20.0, {0, 1, 0, 1}));
    }
    auto windowed = model.phasing(sites, unconstrained(sites.size()));

    LinkageModel::Params q = p;
    q.window = 1000;   // one window, no seams
    LinkageModel exact(q);
    auto whole = exact.phasing(sites, unconstrained(sites.size()));

    REQUIRE(windowed.size() == whole.size());
    REQUIRE(count_switches(windowed) == 0);
    REQUIRE(count_switches(whole) == 0);
}

/// A haploid site: one likelihood per allele, indexed by allele rather than by genotype pair.
static LinkageModel::Site haploid_site(size_t position, double ln_0, double ln_1,
                                       const vector<int>& haps) {
    LinkageModel::Site s;
    s.position = position;
    s.num_alleles = 2;
    s.ploidy = 1;
    s.genotype_ln_likelihood = {ln_0, ln_1};
    s.haplotype_allele = haps;
    return s;
}

TEST_CASE("A haploid chain gets a mosaic", "[linkage_model]") {
    // chrY and non-pseudoautosomal chrX are haploid, and before this existed they were dropped
    // from the linkage pass entirely -- no transition model and no mosaic for about 5% of a
    // genome. The diploid model cannot express them: its state is a pair.
    //
    // Haplotype 1 carries the alleles the reads want over the first half and haplotype 3 over the
    // second, so the one path explaining the chain switches once.
    LinkageModel::Params p;
    p.weight = 2.0;
    LinkageModel model(p);

    vector<LinkageModel::Site> sites;
    for (size_t i = 0; i < 6; ++i) {
        vector<int> haps = (i < 3) ? vector<int>{0, 1, 0, 0} : vector<int>{0, 0, 0, 1};
        sites.push_back(haploid_site(1000 + i * 100, -20.0, 0.0, haps));
    }
    vector<size_t> free_all(sites.size(), LinkageModel::NO_CONSTRAINT);
    auto path = model.haploid_phasing(sites, free_all);
    REQUIRE(path.size() == sites.size());

    // Every site must be explained by a haplotype carrying the allele the reads chose.
    for (size_t t = 0; t < sites.size(); ++t) {
        if (path[t] == LinkageModel::WILDCARD) {
            continue;
        }
        REQUIRE(sites[t].haplotype_allele[path[t]] == 1);
    }
    size_t switches = 0;
    for (size_t t = 1; t < path.size(); ++t) {
        switches += (path[t] != path[t - 1]);
    }
    REQUIRE(switches <= 1);
}

TEST_CASE("Haploid posteriors are a distribution over alleles", "[linkage_model]") {
    LinkageModel::Params p;
    p.weight = 2.0;
    LinkageModel model(p);
    vector<LinkageModel::Site> sites{
        haploid_site(1000, 0.0, -5.0, {0, 1, 0, 1}),
        haploid_site(1100, -5.0, 0.0, {0, 1, 0, 1}),
    };
    auto post = model.haploid_posteriors(sites);
    REQUIRE(post.size() == 2);
    for (const auto& row : post) {
        // One entry per allele, not per genotype pair -- getting that wrong would index a
        // triangular vector with an allele and silently mis-call.
        REQUIRE(row.size() == 2);
        double sum = 0.0;
        for (double v : row) {
            REQUIRE(v >= 0.0);
            sum += v;
        }
        REQUIRE(sum == Approx(1.0).margin(1e-9));
    }
}

TEST_CASE("Constrained haploid phasing spells the called allele", "[linkage_model]") {
    // The same consistency guarantee the diploid path gives: the emitted mosaic must agree with
    // the emitted VCF.
    LinkageModel::Params p;
    p.weight = 2.0;
    LinkageModel model(p);
    vector<LinkageModel::Site> sites;
    vector<size_t> want;
    for (size_t i = 0; i < 6; ++i) {
        // Haplotype 2 carries allele 1 throughout; the reads prefer allele 0 everywhere, so only
        // the constraint can produce the answer.
        sites.push_back(haploid_site(1000 + i * 100, 0.0, -20.0, {0, 0, 1, 0}));
        want.push_back(1);
    }
    auto path = model.haploid_phasing(sites, want);
    REQUIRE(path.size() == 6);
    for (size_t t = 0; t < 6; ++t) {
        REQUIRE(path[t] != LinkageModel::WILDCARD);
        REQUIRE(sites[t].haplotype_allele[path[t]] == 1);
    }
}

TEST_CASE("Haploid window seams do not manufacture switches", "[linkage_model]") {
    LinkageModel::Params p;
    p.weight = 2.0;
    p.window = 5;
    p.margin = 2;
    LinkageModel model(p);
    vector<LinkageModel::Site> sites;
    for (size_t i = 0; i < 40; ++i) {
        sites.push_back(haploid_site(1000 + i * 100, -20.0, 0.0, {0, 1, 0, 1}));
    }
    vector<size_t> free_all(sites.size(), LinkageModel::NO_CONSTRAINT);
    auto path = model.haploid_phasing(sites, free_all);
    size_t switches = 0;
    for (size_t t = 1; t < path.size(); ++t) {
        switches += (path[t] != path[t - 1]);
    }
    REQUIRE(switches == 0);
}

TEST_CASE("Phasing is deterministic", "[linkage_model]") {
    // Two runs over one input must agree exactly. Ties broken by iteration order are
    // deterministic; ties broken by anything else would make the emitted genome irreproducible
    // without moving any accuracy metric.
    LinkageModel::Params p;
    p.weight = 2.0;
    LinkageModel model(p);
    vector<LinkageModel::Site> sites;
    for (size_t i = 0; i < 20; ++i) {
        sites.push_back(biallelic(1000 + i * 137, -3.0, 0.0, -2.0,
                                  {(int)(i % 2), 1, 0, (int)((i + 1) % 2)}));
    }
    auto a = model.phasing(sites, unconstrained(sites.size()));
    auto b = model.phasing(sites, unconstrained(sites.size()));
    REQUIRE(a.size() == b.size());
    for (size_t t = 0; t < a.size(); ++t) {
        REQUIRE(a[t].first == b[t].first);
        REQUIRE(a[t].second == b[t].second);
    }
}


TEST_CASE("A site below depth 1 inherits its parent's strand, not strand 0",
          "[linkage_model]") {
    // A nested parent occupies one haplotype, so everything inside it is on that haplotype. The
    // identity match cannot discover which: a nested haploid parent has trav_first == trav_second
    // and ploidy 1, so the match can only ever answer strand 0 -- and on chr20 that put all 448
    // depth->=2 sites under a strand-1 parent onto the wrong haplotype, while the mosaic read the
    // parent's wildcard hap_first and simultaneously called them unexplained.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    const size_t TOP = 7, MID = 71, DEEP = 711;

    // Try both parent orientations, since which traversal the Viterbi puts on which haplotype is
    // its choice and the grandchild must follow whichever it made.
    for (int which : {0, 1}) {
        LinkageCollector collector(p, 2);
        // A het top-level parent with one panel haplotype on each allele.
        record_dense(collector, "chr1", 1000, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1, TOP,
                     /*share*/ 1.0, /*ploidy*/ 2, /*start*/ 10, /*end*/ 40);
        // A nested haploid child hanging off one of the parent's traversals.
        record_dense(collector, "chr1", 1010, 2, {0.0, -30.0}, {0, 0}, 0, 0, MID,
                     /*share*/ 1.0, /*ploidy*/ 1, /*start*/ 11, /*end*/ 30,
                     /*nested*/ true, TOP,
                     /*crossing*/ (uint64_t)1 << which, /*generation*/ 1);
        // And a grandchild hanging off the child's own single traversal.
        record_dense(collector, "chr1", 1020, 2, {0.0, -30.0}, {0, 0}, 0, 0, DEEP,
                     /*share*/ 1.0, /*ploidy*/ 1, /*start*/ 12, /*end*/ 20,
                     /*nested*/ true, MID,
                     /*crossing*/ (uint64_t)1 << 0, /*generation*/ 2);

        // Each generation in turn, accumulating -- a nested site is only produced by the pass for
        // its own generation, and a deeper site needs its parent already phased.
        vector<LinkageCollector::PhaseCall> phased;
        for (size_t gen = 0; gen <= 2; ++gen) {
            collector.resolve_generation(gen, gen == 2, &phased);
        }

        const LinkageCollector::PhaseCall* mid = nullptr;
        const LinkageCollector::PhaseCall* deep = nullptr;
        for (const auto& pc : phased) {
            if (pc.record_key == MID) {
                mid = &pc;
            } else if (pc.record_key == DEEP) {
                deep = &pc;
            }
        }
        REQUIRE(mid != nullptr);
        REQUIRE(deep != nullptr);
        REQUIRE(mid->nested_strand >= 0);
        // The property: the grandchild is on the same haplotype as the site that contains it. Before
        // this it was always strand 0, so it agreed only when the middle site happened to be there.
        REQUIRE(deep->nested_strand == mid->nested_strand);
        // And it names a haplotype rather than inheriting the parent's wildcard slot.
        const size_t named = deep->nested_strand == 0 ? deep->hap_first : deep->hap_second;
        REQUIRE(named != LinkageModel::WILDCARD);
    }
}

TEST_CASE("A child of a haploid locus gets no strand, so it is not written as a half-call",
          "[linkage_model]") {
    // A strand is a claim that the locus has two haplotypes and this allele sits on one of them.
    // Where the parent has one copy because the whole contig is haploid -- chrX outside the
    // pseudoautosomal regions, or chrY -- that claim is false, and the renderer turns a strand into
    // "a|.", which says the other haplotype carries nothing here. On a haploid contig there is no
    // other haplotype to be empty and the record should be a bare `a`.
    //
    // This is the case the depth->=2 test above cannot reach: its parent is a het DIPLOID site, so
    // the strand it hands down is real. The guard that was missing applied only to `strand = 0`
    // (`strand = 1` was already conditioned on ploidy 2), so a haploid parent matched trav_first and
    // handed down strand 0 unconditionally. Measured on chrX: 8,056 "1|." records in the haploid
    // interior against 1,965 before, and chrX was the one contig whose F1 fell while all 22
    // autosomes rose.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    const size_t TOP = 9, KID = 91;

    LinkageCollector collector(p, 2);
    // A haploid top-level site: one copy because the contig has one, not because a sibling allele
    // deleted anything. It has no parent and therefore no strand of its own.
    record_dense(collector, "chrX", 1000, 2, {0.0, -30.0}, {0, 0}, 0, 0, TOP,
                 /*share*/ 1.0, /*ploidy*/ 1, /*start*/ 10, /*end*/ 40);
    // A nested child hanging off its single traversal.
    record_dense(collector, "chrX", 1010, 2, {0.0, -30.0}, {0, 0}, 0, 0, KID,
                 /*share*/ 1.0, /*ploidy*/ 1, /*start*/ 11, /*end*/ 30,
                 /*nested*/ true, TOP,
                 /*crossing*/ (uint64_t)1 << 0, /*generation*/ 1);

    vector<LinkageCollector::PhaseCall> phased;
    for (size_t gen = 0; gen <= 1; ++gen) {
        collector.resolve_generation(gen, gen == 1, &phased);
    }

    const LinkageCollector::PhaseCall* kid = nullptr;
    for (const auto& pc : phased) {
        if (pc.record_key == KID) {
            kid = &pc;
        }
    }
    REQUIRE(kid != nullptr);
    REQUIRE(kid->ploidy == 1);
    // The whole point: no strand, because the locus has no second haplotype for one to mean
    // anything against. `nested_strand >= 0` is what makes the renderer write "a|." instead of "a".
    REQUIRE(kid->nested_strand < 0);
}

TEST_CASE("The per-strand transition reduces to the single-rho form when both strands agree",
          "[linkage_model]") {
    // The generalisation has to be a generalisation: with one value for both strands it must compute
    // what the scalar version computed. Asserted to 1e-12 rather than bit-identity, and deliberately
    // so -- the old grouping `stay * jump * (row[a] + col[b])` cannot be preserved once the two
    // coefficients differ, so the split into two separately-coefficiented terms is a guaranteed
    // re-association. Demanding bit-identity here would produce a gate that has to be waived.
    const size_t m = 5;
    vector<double> in(m * m);
    // A fixed, uneven state vector. Deterministic rather than random: a failure has to be
    // reproducible, and there is nothing about randomness this test needs.
    for (size_t i = 0; i < m * m; ++i) {
        in[i] = 0.001 + 0.37 * ((double)((i * 7919) % 101) / 101.0);
    }

    // The single-rho form as it was, written out locally so the comparison is against the old
    // arithmetic and not against the new code's own factorisation.
    auto reference = [&](double rho, vector<double>& out) {
        double stay = 1.0 - rho;
        double jump = rho / (double)m;
        vector<double> row(m, 0.0), col(m, 0.0);
        double total = 0.0;
        for (size_t a = 0; a < m; ++a) {
            for (size_t b = 0; b < m; ++b) {
                double v = in[a * m + b];
                row[a] += v;
                col[b] += v;
                total += v;
            }
        }
        out.assign(m * m, 0.0);
        for (size_t a = 0; a < m; ++a) {
            for (size_t b = 0; b < m; ++b) {
                out[a * m + b] = stay * stay * in[a * m + b]
                                 + stay * jump * (row[a] + col[b])
                                 + jump * jump * total;
            }
        }
    };

    for (double rho : {0.0, 1e-4, 0.013, 0.5, 0.97, 1.0}) {
        vector<double> want, got;
        reference(rho, want);
        transition_apply(in, m, rho, rho, got);
        REQUIRE(got.size() == want.size());
        for (size_t i = 0; i < want.size(); ++i) {
            REQUIRE(got[i] == Approx(want[i]).margin(1e-12));
        }
    }
}

TEST_CASE("The per-strand transition can move one strand and leave the other",
          "[linkage_model]") {
    // The point of the pair form, and the case the scalar version could not express at all: one
    // strand certain to stay, the other certain to jump. The closed form is exact -- strand a keeps
    // its index, strand b is uniform over m -- so this is checkable without a reference
    // implementation, and it is what a per-haplotype distance will actually ask for.
    const size_t m = 4;
    vector<double> in(m * m, 0.0);
    // All mass on (2, 1).
    in[2 * m + 1] = 1.0;

    vector<double> out;
    transition_apply(in, m, /*rho_a*/ 0.0, /*rho_b*/ 1.0, out);
    for (size_t a = 0; a < m; ++a) {
        for (size_t b = 0; b < m; ++b) {
            // Row 2 spreads uniformly across b; every other row is empty.
            const double want = (a == 2) ? 1.0 / (double)m : 0.0;
            REQUIRE(out[a * m + b] == Approx(want).margin(1e-12));
        }
    }

    // And the mirror image, so neither axis is privileged by accident.
    transition_apply(in, m, /*rho_a*/ 1.0, /*rho_b*/ 0.0, out);
    for (size_t a = 0; a < m; ++a) {
        for (size_t b = 0; b < m; ++b) {
            const double want = (b == 1) ? 1.0 / (double)m : 0.0;
            REQUIRE(out[a * m + b] == Approx(want).margin(1e-12));
        }
    }
}

TEST_CASE("Two strands with different deletion content get different switch probabilities",
          "[linkage_model]") {
    // Stage 13's capability, kept after stage 15(b) was reverted for failing its own accuracy
    // criterion. What 15(b) tried to do -- feed the two strands genuinely different distances,
    // derived from the indel content each carries -- measured WORSE (JointIndel -0.00061 on chr20
    // against a required +0.0005), so the plumbing is gone. The arithmetic that would carry it is
    // still here and still has to be right, because a later stage may supply better distances than
    // the called-allele approximation could.
    //
    // switch_probability is monotone in the gap and the transition treats its two axes
    // independently: those are the two properties any per-strand distance scheme rests on.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 10000.0;
    p.rho_min = 1e-3;
    LinkageModel model(p);

    const size_t ref_gap = 5000;
    // switch_probability is monotone in the gap, so the strand that travelled less must be stickier.
    REQUIRE(model.switch_probability(ref_gap - 4000) < model.switch_probability(ref_gap));

    // And the transition has to use the two differently rather than symmetrising them. All mass on
    // one state, strand a pinned and strand b loosened: row 1 must keep everything while column 1
    // must not.
    const size_t m = 4;
    vector<double> in(m * m, 0.0);
    in[1 * m + 1] = 1.0;
    vector<double> out;
    transition_apply(in, m, /*rho_a*/ 0.0, /*rho_b*/ 0.8, out);
    double row1 = 0.0, col1 = 0.0;
    for (size_t k = 0; k < m; ++k) {
        row1 += out[1 * m + k];
        col1 += out[k * m + 1];
    }
    REQUIRE(row1 == Approx(1.0).margin(1e-12));
    REQUIRE(col1 < 0.9);
}

TEST_CASE("Every generation is resolved, not only the first", "[linkage_model]") {
    // Chain construction skips entries above the generation being resolved, so resolving generation
    // 0 alone drops every nested site from linkage, from phasing and from the mosaic. The caller
    // that runs the deferred-descent barrier loops the generations itself, which is why this is
    // latent rather than live -- but the invariant is that a recorded site is always settled, and it
    // should not depend on which caller got there first.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    const size_t TOP = 3, DEEP = 31;

    LinkageCollector collector(p, 2);
    record_dense(collector, "chr1", 1000, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1, TOP,
                 /*share*/ 1.0, /*ploidy*/ 2, /*start*/ 10, /*end*/ 20);
    record_dense(collector, "chr1", 1010, 2, {0.0, -30.0}, {0, 0}, 0, 0, DEEP,
                 /*share*/ 1.0, /*ploidy*/ 1, /*start*/ 11, /*end*/ 12,
                 /*nested*/ true, TOP, /*crossing*/ (uint64_t)1 << 0,
                 /*generation*/ 1);

    // resolve() is generation 0 only, by contract. The generation-1 site must not appear.
    vector<LinkageCollector::PhaseCall> gen0;
    collector.resolve(&gen0);
    bool deep_in_gen0 = false;
    for (const auto& pc : gen0) {
        deep_in_gen0 = deep_in_gen0 || pc.record_key == DEEP;
    }
    REQUIRE_FALSE(deep_in_gen0);
    REQUIRE(collector.max_generation() == 1);

    // Resolving the generation it belongs to is what produces it, and that is what any caller
    // reaching the writer must do for every generation the collector holds.
    vector<LinkageCollector::PhaseCall> all = gen0;
    collector.resolve_generation(1, true, &all);
    bool deep_now = false;
    for (const auto& pc : all) {
        deep_now = deep_now || pc.record_key == DEEP;
    }
    REQUIRE(deep_now);
}

TEST_CASE("A nested site drops a haplotype that does not carry the allele it settled on",
          "[linkage_model]") {
    // A nested site takes its haplotype from whichever of its parent's strands carries the chain.
    // That is a claim about the parent's allele, not about this site's: nothing checked that the
    // named panel haplotype carries what the child actually settled on, and the per-strand pass can
    // move that allele afterwards.
    //
    // Naming the wrong haplotype is worse than naming none. A consumer walking it reads a different
    // sequence than the record states, and the two outputs then disagree about the same site with
    // nothing to say which is wrong.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    const size_t PARENT = 4, CHILD = 41;

    LinkageCollector collector(p, 2);
    // A het parent, decisively 0/1, one panel haplotype on each allele so each strand gets one.
    record_dense(collector, "chr1", 1000, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1, PARENT,
                 /*share*/ 1.0, /*ploidy*/ 2, /*start*/ 10, /*end*/ 20);
    // The child hangs off the parent's traversal 0, so it inherits that strand's haplotype -- panel
    // haplotype 1, per the parent's {1, 0} above. At the child, BOTH panel haplotypes carry allele 0
    // while the reads decide allele 1, so whichever haplotype it inherits demonstrably carries
    // something else.
    record_dense(collector, "chr1", 1010, 2, {-30.0, 0.0}, {0, 0}, 1, 1, CHILD,
                 /*share*/ 1.0, /*ploidy*/ 1, /*start*/ 11, /*end*/ 12,
                 /*nested*/ true, PARENT, /*crossing*/ (uint64_t)1 << 0);

    vector<LinkageCollector::PhaseCall> phased;
    collector.resolve(&phased);

    const LinkageCollector::PhaseCall* child = nullptr;
    for (const auto& pc : phased) {
        if (pc.record_key == CHILD) {
            child = &pc;
        }
    }
    REQUIRE(child != nullptr);
    // It is still placed -- the strand is known, the phase set is the parent's -- and it still names
    // an allele. What it no longer does is name a haplotype it contradicts.
    //
    // Which strand index that is, is deliberately not asserted: the Viterbi decides which of the
    // parent's traversals lands on which haplotype, so pinning traversal 0 to strand 0 would be
    // testing an orientation the design specifically does not promise. (It landed on strand 1 here.)
    REQUIRE(child->nested_strand >= 0);
    const size_t named = child->nested_strand == 0 ? child->hap_first : child->hap_second;
    REQUIRE(named == LinkageModel::WILDCARD);
}

TEST_CASE("A revised site stops being unemitted when the revision writes a line",
          "[linkage_model]") {
    // The whole point of recording unemitted sites is that a collapsed parent can still be phased.
    // But "unemitted" is a property of the *current* record, not of the site, and the barrier is
    // exactly where it changes: a chain no called parent allele reached during the sweep is recorded
    // with no line, and becomes a real record once the settled parent turns out to carry it.
    //
    // respecify() did not update the flag, so such a chain stayed unemitted for the rest of the run.
    // Both the genotype patch and the phase patch skip an unemitted entry -- deliberately, there is
    // normally no line to patch -- so the record came out with neither: no linkage correction and no
    // phase set, on 75 of chr20's 117,097 records. Nothing in the output said so; the site simply
    // had a slash where every other record had a bar.
    LinkageModel::Params p;
    p.weight = 1.0;
    p.scale = 100000.0;
    p.rho_min = 1e-4;

    const size_t A = 1, B = 2;

    LinkageCollector collector(p, 2);
    // An ordinary neighbour, so there is a chain for the revised site to be phased within.
    record_dense(collector, "chr1", 1000, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1, A,
                 /*share*/ 1.0, /*ploidy*/ 2, /*start*/ 10, /*end*/ 20);
    // Recorded with no line of its own, as a chain nothing was written for during the sweep.
    record_dense(collector, "chr1", 1010, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1, B,
                 /*share*/ 1.0, /*ploidy*/ 2, /*start*/ 11, /*end*/ 12,
                 /*nested*/ false, /*parent*/ 0, /*crossing*/ 0,
                 /*generation*/ 0, /*emitted*/ false);

    // The barrier revises it and this time a line is written.
    // Re-emitted 3 bp along, because changing the emitted allele set moves POS: this is exactly
    // the case where an entry left at the sweep-time position becomes unreachable by both patches.
    REQUIRE(respecify_dense(collector, B, "chr1", 1013, 2, {-30.0, 0.0, -30.0}, {1, 0}, 0, 1,
                            /*ploidy*/ 2, /*nested*/ false,
                            /*crossing*/ 0, /*emitted*/ true));

    vector<LinkageCollector::PhaseCall> phased;
    collector.resolve(&phased);

    const LinkageCollector::PhaseCall* revised = nullptr;
    for (const auto& pc : phased) {
        if (pc.record_key == B) {
            revised = &pc;
        }
    }
    REQUIRE(revised != nullptr);
    // The assertion that matters: it is patchable. Without it the record is emitted and then
    // skipped by every patch, which is indistinguishable in the output from never having been
    // phased at all.
    REQUIRE(revised->emitted);
    // And it is patchable at the position the replacement line was actually written to. The patch
    // indices are keyed on (contig, POS), so an entry still filed at the sweep-time position is
    // never looked up -- the patches are not declined, they are never offered.
    REQUIRE(revised->position == 1013);
}

}
}
