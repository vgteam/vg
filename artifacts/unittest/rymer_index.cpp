/// \file unittest/rymer_index.cpp
///
/// Unit tests for the RYmer index and for the seed rescue that SAFARI builds on it.

#include <algorithm>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>
#include <gbwtgraph/minimizer.h>

#include "../deamination_profile.hpp"
#include "../rymer_index.hpp"
#include "../rymer_seed_rescue.hpp"
#include "../../src/utility.hpp"
#include "../../src/unittest/catch.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

constexpr size_t TEST_K = 25;
constexpr size_t TEST_W = 5;

/// A fixed pseudo-random DNA sequence.
string random_dna(size_t length, unsigned seed) {
    mt19937 generator(seed);
    string result;
    for (size_t i = 0; i < length; i++) {
        result.push_back("ACGT"[generator() % 4]);
    }
    return result;
}

/// A GBZ graph of one node holding `sequence`, with one path along it.
gbwtgraph::GBZ single_path_gbz(const string& sequence) {
    string filename = temp_file::create();
    {
        ofstream gfa(filename);
        gfa << "H\tVN:Z:1.0\n"
            << "S\t1\t" << sequence << "\n"
            << "P\tref\t1+\t*\n";
    }
    auto [index, source] = gbwtgraph::gfa_to_gbwt(filename);
    temp_file::remove(filename);
    return gbwtgraph::GBZ(index, source);
}

/// An RYmer index of the graph's haplotypes.
gbwtgraph::DefaultMinimizerIndex rymer_index_of(const gbwtgraph::GBZ& gbz) {
    gbwtgraph::DefaultMinimizerIndex index(TEST_K, TEST_W, RYMER_PAYLOAD_SIZE, false);
    index_rymers(gbz, index);
    return index;
}

const string PROFILE_HEADER = "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n";
const string PROFILE_ROW = "0.001\t0.001\t0.001\t0.001\t0.001\t0.1\t0.1\t0.001\t0.001\t0.001\t0.001\t0.001\n";

}

TEST_CASE("RYmer keys keep one bit per base", "[safari][rymer]") {
    // A and G are purines (0), C and T pyrimidines (1).
    gbwtgraph::Key64 dna = gbwtgraph::Key64::encode("ACGTTG");
    REQUIRE(rymer_key(dna, 6).get_key() == 0b010110u);
    // Changing C to T or G to A leaves the RYmer key alone.
    REQUIRE(rymer_key(gbwtgraph::Key64::encode("ATATTA"), 6) == rymer_key(dna, 6));
    // Changing A to C does not.
    REQUIRE(rymer_key(gbwtgraph::Key64::encode("CCGTTG"), 6) != rymer_key(dna, 6));
    REQUIRE(rymer_string("ACGTNacgt") == "RYRYNRYRY");
}

TEST_CASE("The RYmer index finds the graph k-mer of each minimizer", "[safari][rymer]") {
    string sequence = random_dna(150, 12345);
    gbwtgraph::GBZ gbz = single_path_gbz(sequence);
    gbwtgraph::DefaultMinimizerIndex index = rymer_index_of(gbz);

    REQUIRE(is_rymer_index(index));
    REQUIRE(index.k() == TEST_K);
    REQUIRE(index.payload_size() == RYMER_PAYLOAD_SIZE);
    REQUIRE(!index.empty());

    // Every minimizer of the indexed sequence is found under its RYmer, with its DNA key.
    auto minimizers = index.minimizers(sequence);
    REQUIRE(!minimizers.empty());
    for (const auto& minimizer : minimizers) {
        vector<gbwtgraph::Key64> found = graph_kmers_for_rymer(index, to_rymer_minimizer(minimizer, TEST_K));
        bool has_own_key = false;
        for (gbwtgraph::Key64 key : found) {
            has_own_key |= (key == minimizer.key);
        }
        REQUIRE(has_own_key);
    }

    SECTION("An ordinary minimizer index is not taken for an RYmer index") {
        gbwtgraph::DefaultMinimizerIndex plain(TEST_K, TEST_W, 0, false);
        REQUIRE(!is_rymer_index(plain));
    }
}

TEST_CASE("RYmer seed rescue recovers graph k-mers that deamination changed", "[safari][rymer]") {
    string sequence = random_dna(150, 6789);
    gbwtgraph::GBZ gbz = single_path_gbz(sequence);
    gbwtgraph::DefaultMinimizerIndex index = rymer_index_of(gbz);
    vector<SubstitutionMatrix> rows = DeaminationProfile::parse_profile(PROFILE_HEADER + PROFILE_ROW, "test");
    DeaminationProfile profile(rows, rows);
    RymerSeedRescue rescue(index, profile, 0.5, DamageModelParameters());

    string read = sequence.substr(20, 80);

    // Each stand-in sequence holds bases of the graph at their read offsets and N elsewhere.
    auto check_stand_ins = [&](const vector<string>& stand_ins) {
        for (const string& stand_in : stand_ins) {
            REQUIRE(stand_in.size() == read.size());
            size_t bases = 0;
            for (size_t i = 0; i < read.size(); i++) {
                if (stand_in[i] != 'N') {
                    REQUIRE(stand_in[i] == read[i]);
                    bases++;
                }
            }
            REQUIRE(bases == TEST_K);
        }
    };

    SECTION("An undamaged read gets stand-ins equal to its own k-mers") {
        vector<string> stand_ins = rescue.stand_in_sequences(read);
        REQUIRE(!stand_ins.empty());
        check_stand_ins(stand_ins);
    }

    SECTION("A C read as T is restored in some stand-in") {
        size_t restored = 0;
        for (size_t p = 0; p < read.size(); p++) {
            if (read[p] != 'C') {
                continue;
            }
            string damaged = read;
            damaged[p] = 'T';
            vector<string> stand_ins = rescue.stand_in_sequences(damaged);
            check_stand_ins(stand_ins);
            auto restores_base = [p](const string& stand_in) { return stand_in[p] == 'C'; };
            if (ranges::any_of(stand_ins, restores_base)) {
                restored++;
            }
        }
        REQUIRE(restored > 0);
    }

    SECTION("A change of purine to pyrimidine is never restored") {
        for (size_t p = 0; p < read.size(); p++) {
            if (read[p] != 'A') {
                continue;
            }
            string changed = read;
            changed[p] = 'C';
            for (const string& stand_in : rescue.stand_in_sequences(changed)) {
                REQUIRE(stand_in[p] != 'A');
            }
        }
    }

    SECTION("Settings out of range are rejected") {
        REQUIRE_THROWS_AS(RymerSeedRescue(index, profile, 1.5, DamageModelParameters()), invalid_argument);
        gbwtgraph::DefaultMinimizerIndex plain(TEST_K, TEST_W, 0, false);
        REQUIRE_THROWS_AS(RymerSeedRescue(plain, profile, 0.5, DamageModelParameters()), invalid_argument);
    }
}

}
}
