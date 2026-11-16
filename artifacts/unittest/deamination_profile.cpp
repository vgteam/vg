/// \file unittest/deamination_profile.cpp
///
/// Unit tests for reading deamination profiles and for the posterior that SAFARI uses to accept
/// a graph k-mer as the undamaged copy of a read k-mer.

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "../deamination_profile.hpp"
#include "../../src/unittest/catch.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

const string PROFILE_HEADER = "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n";

/// A profile row in which C becomes T and G becomes A at the given rates, and every other change
/// has rate 0.001.
string damage_row(double c_to_t, double g_to_a) {
    return "0.001\t0.001\t0.001\t0.001\t0.001\t" + to_string(c_to_t) + "\t"
         + to_string(g_to_a) + "\t0.001\t0.001\t0.001\t0.001\t0.001\n";
}

}

TEST_CASE("A deamination profile is read from a .prof table", "[safari][deamination]") {
    vector<SubstitutionMatrix> rows = DeaminationProfile::parse_profile(
        PROFILE_HEADER + damage_row(0.3, 0.02) + damage_row(0.1, 0.01), "test");

    REQUIRE(rows.size() == 2);
    // Row 0, true base C (1): C>A, C>G, C>T, and C read unchanged with what is left.
    REQUIRE(rows[0][1][0] == Approx(0.001));
    REQUIRE(rows[0][1][2] == Approx(0.001));
    REQUIRE(rows[0][1][3] == Approx(0.3));
    REQUIRE(rows[0][1][1] == Approx(1.0 - 0.302));
    // Row 1, true base G (2) becomes A (0).
    REQUIRE(rows[1][2][0] == Approx(0.01));
    REQUIRE(rows[1][2][2] == Approx(1.0 - 0.012));
    // Each true base's probabilities sum to 1.
    for (const SubstitutionMatrix& row : rows) {
        for (size_t true_base = 0; true_base < 4; true_base++) {
            double total = 0;
            for (size_t read_base = 0; read_base < 4; read_base++) {
                total += row[true_base][read_base];
            }
            REQUIRE(total == Approx(1.0));
        }
    }
}

TEST_CASE("Invalid deamination profiles are rejected", "[safari][deamination]") {
    SECTION("A row with the wrong number of fields") {
        REQUIRE_THROWS_AS(DeaminationProfile::parse_profile(PROFILE_HEADER + "0.1\t0.2\n", "test"), runtime_error);
    }
    SECTION("A header with no rows") {
        REQUIRE_THROWS_AS(DeaminationProfile::parse_profile(PROFILE_HEADER, "test"), runtime_error);
    }
    SECTION("Rates for one base summing to more than 1") {
        REQUIRE_THROWS_AS(DeaminationProfile::parse_profile(PROFILE_HEADER + damage_row(0.999, 0.0), "test"), runtime_error);
    }
    SECTION("A field that is not a number") {
        string row = "x\t0.001\t0.001\t0.001\t0.001\t0.1\t0.1\t0.001\t0.001\t0.001\t0.001\t0.001\n";
        REQUIRE_THROWS_AS(DeaminationProfile::parse_profile(PROFILE_HEADER + row, "test"), runtime_error);
    }
}

TEST_CASE("Substitution probabilities use the end more likely to change the base", "[safari][deamination]") {
    // C becomes T often at the 5' end, and G becomes A often at the 3' end.
    vector<SubstitutionMatrix> five = DeaminationProfile::parse_profile(
        PROFILE_HEADER + damage_row(0.3, 0.001) + damage_row(0.01, 0.001), "five");
    vector<SubstitutionMatrix> three = DeaminationProfile::parse_profile(
        PROFILE_HEADER + damage_row(0.001, 0.25) + damage_row(0.001, 0.02), "three");
    DeaminationProfile profile(five, three);

    // At the first base of a 50-base fragment the 5' rates apply to C.
    REQUIRE(profile.substitution_probability(50, 0, 'C', 'T') == Approx(0.3));
    // At the last base the 3' rates apply to G.
    REQUIRE(profile.substitution_probability(50, 49, 'G', 'A') == Approx(0.25));
    // At the first base, G changes more at the 3' end's 50th base, which takes the last 3' row.
    REQUIRE(profile.substitution_probability(50, 0, 'G', 'A') == Approx(0.02));
    // In the middle both ends take their last rows.
    REQUIRE(profile.substitution_probability(50, 25, 'C', 'T') == Approx(0.01));
    // Fragments longer than 1000 bases work too.
    REQUIRE(profile.substitution_probability(5000, 4999, 'G', 'A') == Approx(0.25));
    // An N on either side has probability one quarter.
    REQUIRE(profile.substitution_probability(50, 0, 'N', 'T') == 0.25);
    REQUIRE(profile.substitution_probability(50, 0, 'C', 'N') == 0.25);
}

TEST_CASE("The damage posterior separates damage from other differences", "[safari][deamination]") {
    vector<SubstitutionMatrix> rows = DeaminationProfile::parse_profile(PROFILE_HEADER + damage_row(0.1, 0.1), "test");
    DeaminationProfile profile(rows, rows);
    DamageModelParameters parameters;

    // 25 bases, with 6 Cs and 6 Gs.
    string graph_kmer = "ACGTACGTACGTACGTACGTACGTA";

    SECTION("Equal k-mers have posterior 1") {
        REQUIRE(damage_posterior(graph_kmer, graph_kmer, 60, 0, profile, parameters) == 1.0);
    }

    SECTION("A C read as T is likely damage") {
        string read_kmer = graph_kmer;
        read_kmer[1] = 'T';
        double posterior = damage_posterior(graph_kmer, read_kmer, 60, 0, profile, parameters);
        REQUIRE(posterior > 0.5);
        REQUIRE(posterior < 1.0);
    }

    SECTION("An A read as C is not likely damage") {
        string read_kmer = graph_kmer;
        read_kmer[0] = 'C';
        REQUIRE(damage_posterior(graph_kmer, read_kmer, 60, 0, profile, parameters) < 0.5);
    }

    SECTION("A higher prior on spurious matches lowers the posterior") {
        string read_kmer = graph_kmer;
        read_kmer[1] = 'T';
        DamageModelParameters skeptical;
        skeptical.spurious_prior = 0.99;
        REQUIRE(damage_posterior(graph_kmer, read_kmer, 60, 0, profile, skeptical)
                < damage_posterior(graph_kmer, read_kmer, 60, 0, profile, parameters));
    }

    SECTION("A change the profile never makes has posterior 0") {
        SubstitutionMatrix exact{};
        for (size_t base = 0; base < 4; base++) {
            exact[base][base] = 1.0;
        }
        DeaminationProfile no_damage({exact}, {exact});
        string read_kmer = graph_kmer;
        read_kmer[1] = 'T';
        REQUIRE(damage_posterior(graph_kmer, read_kmer, 60, 0, no_damage, parameters) == 0.0);
    }

    SECTION("Short k-mers use the binomial model of spurious matches") {
        string short_graph = "ACGTACGTACGTACGTAC";
        string short_read = short_graph;
        short_read[1] = 'T';
        double posterior = damage_posterior(short_graph, short_read, 60, 0, profile, parameters);
        REQUIRE(posterior > 0.0);
        REQUIRE(posterior < 1.0);
        REQUIRE(std::isfinite(posterior));
    }
}

}
}
