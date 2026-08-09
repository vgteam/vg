/// \file unittest/allele_likelihood.cpp
///
/// Unit tests for the read-level genotype likelihood model.
///
/// These build the reads x alleles matrix by hand, so the genotyping maths is
/// tested independently of the scoring that normally produces it. Several of
/// these are cases where a subtly wrong implementation still produces
/// plausible-looking VCF, which is exactly why they are pinned explicitly.
///

#include <cmath>
#include <limits>
#include <vector>

#include "allele_likelihood.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

using namespace std;

static const double NEG_INF = -numeric_limits<double>::infinity();

/// Find the highest scoring genotype among those scored.
static vector<int> best_genotype(const vector<pair<vector<int>, double>>& scored) {
    REQUIRE(!scored.empty());
    size_t best = 0;
    for (size_t i = 1; i < scored.size(); ++i) {
        if (scored[i].second > scored[best].second) {
            best = i;
        }
    }
    return scored[best].first;
}

TEST_CASE("Row normalisation puts every row's maximum at exactly 1", "[allele_likelihood]") {
    // The central invariant of the matrix: each read's likelihoods are relative
    // to that read's own best explanation, so the background term in the
    // genotype likelihood is dimensionless and exactly 1.
    AlleleReadLikelihoodsBuilder builder(3);
    builder.add_read({-5.0, -12.0, -30.0}, 1e-6, "read_a");
    builder.add_read({-100.0, -100.5, -101.0}, 1e-6, "read_b");
    builder.add_read({0.0, NEG_INF, -1.0}, 1e-6, "read_c");
    AlleleReadLikelihoods matrix = builder.build();

    REQUIRE(matrix.num_reads() == 3);
    REQUIRE(matrix.num_alleles() == 3);

    for (size_t r = 0; r < matrix.num_reads(); ++r) {
        double row_max = 0.0;
        for (size_t a = 0; a < matrix.num_alleles(); ++a) {
            REQUIRE(matrix.rel(r, a) >= 0.0);
            REQUIRE(matrix.rel(r, a) <= 1.0);
            row_max = max(row_max, matrix.rel(r, a));
        }
        REQUIRE(row_max == Approx(1.0));
    }

    SECTION("the absolute best fit is retained separately") {
        // Not used by the model, but it is the only surviving record of absolute
        // fit, so a read that fits everything badly stays distinguishable.
        REQUIRE(matrix.best_ln_likelihood(0) == Approx(-5.0));
        REQUIRE(matrix.best_ln_likelihood(1) == Approx(-100.0));
        REQUIRE(matrix.best_ln_likelihood(2) == Approx(0.0));
    }

    SECTION("no valid placement becomes exactly zero, not a small number") {
        REQUIRE(matrix.rel(2, 1) == 0.0);
    }
}

TEST_CASE("A read placing on no allele at all is dropped and counted", "[allele_likelihood]") {
    // Reachable in practice: retrieval fetches by node ID range, so a read can
    // overlap the range yet place on no traversal. Normalising it would divide
    // by zero and quietly poison every genotype at the site with NaN.
    AlleleReadLikelihoodsBuilder builder(2);
    builder.add_read({0.0, -3.0}, 1e-6, "good");
    builder.add_read({NEG_INF, NEG_INF}, 1e-6, "hopeless");
    builder.add_read({-1.0, 0.0}, 1e-6, "also_good");
    AlleleReadLikelihoods matrix = builder.build();

    REQUIRE(matrix.num_reads() == 2);
    REQUIRE(matrix.num_unplaceable() == 1);

    // Nothing downstream sees a NaN.
    for (auto& scored : matrix.score_genotypes(2)) {
        REQUIRE(!std::isnan(scored.second));
        REQUIRE(std::isfinite(scored.second));
    }
}

TEST_CASE("Clear heterozygous evidence calls the heterozygote", "[allele_likelihood]") {
    // Half the reads fit allele 0 and not allele 1, half the reverse. Allele
    // balance has to emerge from the 1/|G| mixture weights, with no het bias
    // parameter anywhere.
    AlleleReadLikelihoodsBuilder builder(2);
    for (int i = 0; i < 10; ++i) {
        builder.add_read({0.0, -20.0}, 1e-6);
    }
    for (int i = 0; i < 10; ++i) {
        builder.add_read({-20.0, 0.0}, 1e-6);
    }
    AlleleReadLikelihoods matrix = builder.build();

    auto scored = matrix.score_genotypes(2);
    REQUIRE(best_genotype(scored) == vector<int>({0, 1}));
}

TEST_CASE("Clear homozygous evidence calls the homozygote", "[allele_likelihood]") {
    AlleleReadLikelihoodsBuilder builder(2);
    for (int i = 0; i < 20; ++i) {
        builder.add_read({0.0, -20.0}, 1e-6);
    }
    AlleleReadLikelihoods matrix = builder.build();

    auto scored = matrix.score_genotypes(2);
    REQUIRE(best_genotype(scored) == vector<int>({0, 0}));

    SECTION("a homozygous genotype reduces to the sum of that allele's column") {
        // The two 1/2 weights are identical and sum to 1, so hom needs no
        // special casing. With a negligible mismap probability the total should
        // be the plain sum of ln rel(r, 0), which here is ln(1) = 0 per read.
        double hom_ref = matrix.genotype_likelihood({0, 0});
        REQUIRE(hom_ref == Approx(0.0).margin(1e-4));
    }
}

TEST_CASE("set_max_allele removes the heterozygote's mixture penalty", "[allele_likelihood]") {
    // The heterozygous-deletion failure in miniature. Allele 0 is a long allele,
    // allele 1 a deletion. Twenty "interior" reads fit only the long allele --
    // on a real heterozygous deletion they all come from the intact haplotype --
    // and three "junction" reads fit only the deletion.
    AlleleReadLikelihoodsBuilder builder(2);
    for (int i = 0; i < 20; ++i) {
        builder.add_read({0.0, -30.0}, 1e-6);   // interior: long only
    }
    for (int i = 0; i < 3; ++i) {
        builder.add_read({-30.0, 0.0}, 1e-6);   // junction: deletion only
    }

    SECTION("under the mixture the interior reads outvote the junction reads") {
        AlleleReadLikelihoods matrix = builder.build();
        // Each interior read prefers hom-long by ln 2; each junction read prefers
        // the het by ln(0.5/e_r), but there are far fewer of them.
        REQUIRE(matrix.genotype_likelihood({0, 0}) > matrix.genotype_likelihood({0, 1}));
        REQUIRE(best_genotype(matrix.score_genotypes(2)) == vector<int>({0, 0}));
    }

    SECTION("under a maximum the interior reads cancel and the het is recovered") {
        builder.set_max_allele(true);
        AlleleReadLikelihoods matrix = builder.build();
        REQUIRE(matrix.genotype_likelihood({0, 1}) > matrix.genotype_likelihood({0, 0}));
        REQUIRE(best_genotype(matrix.score_genotypes(2)) == vector<int>({0, 1}));
    }
}

TEST_CASE("set_max_allele cannot score a heterozygote below a homozygote",
          "[allele_likelihood]") {
    // Why this is a diagnostic and not a shipping default. max over the genotype
    // is monotone in the allele set, so adding an allele never lowers any read's
    // term: 0/1 dominates both 0/0 and 1/1 at every site. Evidence that should
    // call a confident homozygote calls a heterozygote instead.
    AlleleReadLikelihoodsBuilder builder(2);
    for (int i = 0; i < 30; ++i) {
        builder.add_read({0.0, -30.0}, 1e-6);   // overwhelming hom-ref evidence
    }
    builder.add_read({-30.0, 0.0}, 1e-6);       // one stray read, e.g. an error
    builder.set_max_allele(true);
    AlleleReadLikelihoods matrix = builder.build();

    REQUIRE(matrix.genotype_likelihood({0, 1}) >= matrix.genotype_likelihood({0, 0}));
    REQUIRE(matrix.genotype_likelihood({0, 1}) >= matrix.genotype_likelihood({1, 1}));
    REQUIRE(best_genotype(matrix.score_genotypes(2)) == vector<int>({0, 1}));
}

TEST_CASE("A flat matrix yields flat genotype likelihoods", "[allele_likelihood]") {
    // No evidence must produce no preference. This is the correct behaviour for
    // a depth-agnostic model, and it is why a site no read spans reports a
    // no-call rather than a confident hom-ref.
    AlleleReadLikelihoodsBuilder builder(3);
    for (int i = 0; i < 5; ++i) {
        builder.add_read({0.0, 0.0, 0.0}, 1e-6);
    }
    AlleleReadLikelihoods matrix = builder.build();

    auto scored = matrix.score_genotypes(2);
    REQUIRE(scored.size() == 6);
    for (auto& entry : scored) {
        REQUIRE(entry.second == Approx(scored.front().second));
    }
}

TEST_CASE("The mismapping term bounds how much one read can penalise a genotype",
          "[allele_likelihood]") {
    // A read that fits allele 0 perfectly and allele 1 not at all is strong
    // evidence against {1,1} -- but bounded, where without the background term
    // it would be unboundedly negative.
    AlleleReadLikelihoodsBuilder builder(2);
    builder.add_read({0.0, NEG_INF}, 1e-3);
    AlleleReadLikelihoods matrix = builder.build();

    double against = matrix.genotype_likelihood({1, 1});
    double e_r = matrix.mismap_prob(0);

    REQUIRE(std::isfinite(against));
    // With rel == 0 for the only allele in the genotype, the bracket collapses
    // to exactly e_r.
    REQUIRE(against == Approx(log(e_r)));

    SECTION("the per-read term always lies within [ln e_r, 0]") {
        for (auto& scored : matrix.score_genotypes(2)) {
            REQUIRE(scored.second <= 0.0 + 1e-12);
            REQUIRE(scored.second >= log(e_r) - 1e-12);
        }
    }
}

TEST_CASE("Lower MAPQ shrinks a read's influence without flipping the ranking",
          "[allele_likelihood]") {
    // e_r is genotype-independent, so a wrong or pessimistic MAPQ can only
    // compress a read's contribution. It must never reverse which genotype is
    // preferred: that is what makes taking MAPQ at face value safe.
    auto build_with = [](double mismap) {
        AlleleReadLikelihoodsBuilder builder(2, 1e-8, 0.4);
        for (int i = 0; i < 6; ++i) {
            builder.add_read({0.0, -15.0}, mismap);
        }
        return builder.build();
    };

    AlleleReadLikelihoods confident = build_with(1e-6);   // ~MAPQ 60
    AlleleReadLikelihoods unsure = build_with(0.4);       // very low MAPQ, clamped

    double confident_gap =
        confident.genotype_likelihood({0, 0}) - confident.genotype_likelihood({1, 1});
    double unsure_gap = unsure.genotype_likelihood({0, 0}) - unsure.genotype_likelihood({1, 1});

    // Same direction...
    REQUIRE(confident_gap > 0.0);
    REQUIRE(unsure_gap > 0.0);
    // ...but the low-confidence reads discriminate much less.
    REQUIRE(unsure_gap < confident_gap);

    REQUIRE(best_genotype(confident.score_genotypes(2)) == vector<int>({0, 0}));
    REQUIRE(best_genotype(unsure.score_genotypes(2)) == vector<int>({0, 0}));
}

TEST_CASE("The mismapping probability is clamped at both ends", "[allele_likelihood]") {
    // The upper clamp is load-bearing rather than hygiene. Many mappers use
    // MAPQ 0 to mean "multi-mapping" rather than P(wrong) = 1, and an unclamped
    // e_r of 1 would collapse the read's term to ln(1) = 0 for every genotype,
    // so the read would silently contribute nothing at all.
    AlleleReadLikelihoodsBuilder builder(2, 1e-8, 0.1);
    builder.add_read({0.0, -10.0}, 1.0);   // as if from MAPQ 0
    builder.add_read({0.0, -10.0}, 0.0);   // as if from an impossibly good MAPQ
    AlleleReadLikelihoods matrix = builder.build();

    REQUIRE(matrix.mismap_prob(0) == Approx(0.1));
    REQUIRE(matrix.mismap_prob(1) == Approx(1e-8));

    // Strictly inside (0,1), which is what keeps every per-read log finite.
    for (size_t r = 0; r < matrix.num_reads(); ++r) {
        REQUIRE(matrix.mismap_prob(r) > 0.0);
        REQUIRE(matrix.mismap_prob(r) < 1.0);
    }
}

TEST_CASE("Reads that cannot tell two alleles apart still discriminate against a third",
          "[allele_likelihood]") {
    // A read whose window spans a region where alleles 0 and 1 agree must score
    // identically on both, staying neutral between them while still ruling out
    // allele 2. This is how a partially overlapping read contributes partial
    // information rather than none or too much.
    AlleleReadLikelihoodsBuilder builder(3);
    for (int i = 0; i < 8; ++i) {
        builder.add_read({0.0, 0.0, -25.0}, 1e-6);
    }
    AlleleReadLikelihoods matrix = builder.build();

    for (size_t r = 0; r < matrix.num_reads(); ++r) {
        REQUIRE(matrix.rel(r, 0) == Approx(matrix.rel(r, 1)));
        REQUIRE(matrix.rel(r, 2) < matrix.rel(r, 0));
    }

    // {0,0}, {0,1} and {1,1} are indistinguishable; anything involving 2 is worse.
    double hom0 = matrix.genotype_likelihood({0, 0});
    REQUIRE(matrix.genotype_likelihood({0, 1}) == Approx(hom0));
    REQUIRE(matrix.genotype_likelihood({1, 1}) == Approx(hom0));
    REQUIRE(matrix.genotype_likelihood({2, 2}) < hom0);
    REQUIRE(matrix.genotype_likelihood({0, 2}) < hom0);
}

TEST_CASE("Reads with very different amounts of evidence combine correctly",
          "[allele_likelihood]") {
    // Reads contributing unequal amounts of information is correct, and must not
    // be "fixed" by normalising per base: the row normalisation is per read, so
    // it is insensitive to how long that read's window was.
    AlleleReadLikelihoodsBuilder builder(2);
    builder.add_read({0.0, -40.0}, 1e-6, "long_read");    // lots of evidence
    builder.add_read({0.0, -1.0}, 1e-6, "short_read");    // barely any
    AlleleReadLikelihoods matrix = builder.build();

    REQUIRE(best_genotype(matrix.score_genotypes(2)) == vector<int>({0, 0}));
    // The weakly informative read has not been discarded, and it does pull a
    // little: {0,1} must beat {1,1}.
    REQUIRE(matrix.genotype_likelihood({0, 1}) > matrix.genotype_likelihood({1, 1}));
}

TEST_CASE("Genotypes are enumerated in VCF GL order", "[allele_likelihood]") {
    // A genotype's position in the enumeration is its GL field index, so getting
    // this wrong silently mislabels every likelihood in the VCF.
    SECTION("diploid, three alleles") {
        auto genotypes = AlleleReadLikelihoods::enumerate_genotypes(3, 2);
        REQUIRE(genotypes.size() == 6);
        REQUIRE(genotypes[0] == vector<int>({0, 0}));
        REQUIRE(genotypes[1] == vector<int>({0, 1}));
        REQUIRE(genotypes[2] == vector<int>({1, 1}));
        REQUIRE(genotypes[3] == vector<int>({0, 2}));
        REQUIRE(genotypes[4] == vector<int>({1, 2}));
        REQUIRE(genotypes[5] == vector<int>({2, 2}));
    }

    SECTION("haploid") {
        auto genotypes = AlleleReadLikelihoods::enumerate_genotypes(3, 1);
        REQUIRE(genotypes.size() == 3);
        REQUIRE(genotypes[0] == vector<int>({0}));
        REQUIRE(genotypes[1] == vector<int>({1}));
        REQUIRE(genotypes[2] == vector<int>({2}));
    }

    SECTION("the count matches the standard formula for any ploidy") {
        // (n + p - 1) choose p
        REQUIRE(AlleleReadLikelihoods::enumerate_genotypes(4, 2).size() == 10);
        REQUIRE(AlleleReadLikelihoods::enumerate_genotypes(2, 3).size() == 4);
        REQUIRE(AlleleReadLikelihoods::enumerate_genotypes(1, 2).size() == 1);
    }

    SECTION("every genotype is sorted non-decreasing") {
        for (auto& genotype : AlleleReadLikelihoods::enumerate_genotypes(4, 3)) {
            for (size_t i = 1; i < genotype.size(); ++i) {
                REQUIRE(genotype[i - 1] <= genotype[i]);
            }
        }
    }
}

TEST_CASE("Haploid genotyping works without special casing", "[allele_likelihood]") {
    AlleleReadLikelihoodsBuilder builder(3);
    for (int i = 0; i < 6; ++i) {
        builder.add_read({-30.0, 0.0, -30.0}, 1e-6);
    }
    AlleleReadLikelihoods matrix = builder.build();

    REQUIRE(best_genotype(matrix.score_genotypes(1)) == vector<int>({1}));
}

TEST_CASE("Multi-allelic sites need no extra machinery", "[allele_likelihood]") {
    // A five-allele site is fifteen genotypes; exhaustive enumeration is cheap
    // once the matrix exists, so there is no candidate pruning to get wrong.
    AlleleReadLikelihoodsBuilder builder(5);
    for (int i = 0; i < 7; ++i) {
        builder.add_read({-20.0, 0.0, -20.0, -20.0, -20.0}, 1e-6);
    }
    for (int i = 0; i < 7; ++i) {
        builder.add_read({-20.0, -20.0, -20.0, 0.0, -20.0}, 1e-6);
    }
    AlleleReadLikelihoods matrix = builder.build();

    auto scored = matrix.score_genotypes(2);
    REQUIRE(scored.size() == 15);
    REQUIRE(best_genotype(scored) == vector<int>({1, 3}));
}

TEST_CASE("Negative allele markers in a genotype are ignored rather than read out of bounds",
          "[allele_likelihood]") {
    // The VCF layer uses negative sentinels for star and missing alleles in
    // nested calling mode, and they reach update_vcf_info inside the genotype
    // vector. Indexing the matrix with them would read out of bounds.
    AlleleReadLikelihoodsBuilder builder(2);
    builder.add_read({0.0, -10.0}, 1e-6);
    AlleleReadLikelihoods matrix = builder.build();

    double with_marker = matrix.genotype_likelihood({0, -1});
    REQUIRE(std::isfinite(with_marker));

    double only_marker = matrix.genotype_likelihood({-2, -1});
    REQUIRE(std::isfinite(only_marker));
}

TEST_CASE("An empty matrix is handled without dividing by zero", "[allele_likelihood]") {
    // A site no read overlaps. Every genotype scores 0 (an empty product), which
    // is a flat likelihood: correct, and distinct from a confident call.
    AlleleReadLikelihoodsBuilder builder(2);
    AlleleReadLikelihoods matrix = builder.build();

    REQUIRE(matrix.num_reads() == 0);
    auto scored = matrix.score_genotypes(2);
    REQUIRE(scored.size() == 3);
    for (auto& entry : scored) {
        REQUIRE(entry.second == Approx(0.0));
    }
}

}
}
