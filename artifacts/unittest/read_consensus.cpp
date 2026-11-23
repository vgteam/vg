/// \file read_consensus.cpp
///
/// Unit tests for the consensus of reads that span a window, and a hidden test that compares the
/// consensus with recorded ones over a file of cases.

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>

#include "../../src/unittest/catch.hpp"
#include "../read_consensus.hpp"

namespace vg {
namespace unittest {

/// Return `length` random bases.
static string random_bases(size_t length, std::mt19937& generator) {
    string bases;
    for (size_t i = 0; i < length; ++i) {
        bases += "ACGT"[generator() % 4];
    }
    return bases;
}

/// Return `sequence` with about one base in `every` changed, removed or doubled.
static string with_errors(const string& sequence, size_t every, std::mt19937& generator) {
    string read;
    for (char base : sequence) {
        if (generator() % every != 0) {
            read += base;
            continue;
        }
        switch (generator() % 3) {
        case 0:
            read += base == 'A' ? 'C' : 'A';
            break;
        case 1:
            break;
        default:
            read += base;
            read += "ACGT"[generator() % 4];
            break;
        }
    }
    return read;
}

TEST_CASE("A read consensus takes the base most reads have", "[read_consensus]") {
    std::mt19937 generator(1);
    const string backbone = random_bases(200, generator);
    string changed = backbone;
    changed[100] = changed[100] == 'G' ? 'T' : 'G';
    vector<string> reads(5, changed);
    reads.resize(8, backbone);

    REQUIRE(read_consensus(reads, backbone) == changed);
    // Without a majority the backbone stays.
    reads.resize(5);
    reads.resize(11, backbone);
    REQUIRE(read_consensus(reads, backbone) == backbone);
    // Too few reads give nothing.
    reads.resize(2);
    REQUIRE(!read_consensus(reads, backbone));
}

TEST_CASE("A read consensus spells a long insertion that reads with errors carry",
          "[read_consensus]") {
    std::mt19937 generator(2);
    const string left = random_bases(300, generator);
    const string novel = random_bases(400, generator);
    const string right = random_bases(300, generator);
    const string backbone = left + right;
    const string truth = left + novel + right;
    vector<string> reads;
    for (size_t i = 0; i < 12; ++i) {
        reads.push_back(with_errors(truth, 40, generator));
    }

    ReadConsensusParams params;
    for (bool seed : {false, true}) {
        params.seed_from_reads = seed;
        const optional<string> consensus = read_consensus(reads, backbone, params);
        REQUIRE(consensus);
        REQUIRE(edit_distance(*consensus, truth) <= 2);
    }
}

TEST_CASE("A read consensus changes the length of a tandem repeat", "[read_consensus]") {
    std::mt19937 generator(3);
    const string left = random_bases(150, generator);
    const string right = random_bases(150, generator);
    const string unit = "CAGGTTCA";
    string six;
    for (size_t i = 0; i < 6; ++i) {
        six += unit;
    }
    const string backbone = left + six + right;
    const string longer = left + six + unit + unit + unit + right;
    const string shorter = left + unit + unit + right;

    vector<string> reads;
    for (size_t i = 0; i < 9; ++i) {
        reads.push_back(with_errors(longer, 60, generator));
    }
    optional<string> consensus = read_consensus(reads, backbone);
    REQUIRE(consensus);
    REQUIRE(edit_distance(*consensus, longer) <= 1);

    reads.clear();
    for (size_t i = 0; i < 9; ++i) {
        reads.push_back(with_errors(shorter, 60, generator));
    }
    consensus = read_consensus(reads, backbone);
    REQUIRE(consensus);
    REQUIRE(edit_distance(*consensus, shorter) <= 1);
}

TEST_CASE("A read consensus keeps the backbone's lengths of runs of one base when asked",
          "[read_consensus]") {
    const string backbone = "ACGTACGGATTTTTTTCAGCATGCAGGCTAGCATCGA";
    // One T fewer in the run, and a C in place of a G.
    const string changed = "ACGTACGGATTTTTTCAGCATCCAGGCTAGCATCGA";
    const vector<string> reads(6, changed);

    ReadConsensusParams params;
    params.keep_homopolymer_lengths = false;
    REQUIRE(read_consensus(reads, backbone, params) == changed);
    params.keep_homopolymer_lengths = true;
    REQUIRE(read_consensus(reads, backbone, params) == "ACGTACGGATTTTTTTCAGCATCCAGGCTAGCATCGA");

    REQUIRE(collapse_homopolymers("AAACCGTTTA") == "ACGTA");
    // A deleted base that is a whole run is not a change of length.
    REQUIRE(keep_homopolymer_lengths("ACTG", "ACGTG") == "ACTG");
    REQUIRE(keep_homopolymer_lengths("ACGGGTG", "ACGTG") == "ACGTG");
}

/// Return the string value of `key` in the one-line JSON object `line`, which has no escapes.
static string json_string(const string& line, const string& key) {
    const string marker = "\"" + key + "\": \"";
    const size_t start = line.find(marker);
    REQUIRE(start != string::npos);
    const size_t from = start + marker.size();
    return line.substr(from, line.find('"', from) - from);
}

/// Return the strings of the array `key` in the one-line JSON object `line`.
static vector<string> json_strings(const string& line, const string& key) {
    const string marker = "\"" + key + "\": [";
    const size_t start = line.find(marker);
    REQUIRE(start != string::npos);
    vector<string> values;
    size_t at = start + marker.size();
    const size_t end = line.find(']', at);
    while (true) {
        const size_t open = line.find('"', at);
        if (open == string::npos || open > end) {
            break;
        }
        const size_t close = line.find('"', open + 1);
        values.push_back(line.substr(open + 1, close - open - 1));
        at = close + 1;
    }
    return values;
}

/// Return the middle value of `values`.
static int64_t median_of(vector<int64_t> values) {
    std::ranges::sort(values);
    return values[values.size() / 2];
}

// Reads the file named by VG_CONSENSUS_CASES, one JSON object a line with `backbone`, `reads`
// and `abpoa`, a recorded consensus. Writes, for each case and in sum, the edit distance from
// the engine's consensus and from the backbone to the recorded one, with and without runs of one
// base collapsed, and the summed edit distance from the reads to the engine's consensus and to
// the recorded one. VG_CONSENSUS_SEED=0 turns off seeding from the reads.
TEST_CASE("The read consensus is compared with recorded consensuses", "[.][read_consensus_cases]") {
    const char* path = std::getenv("VG_CONSENSUS_CASES");
    REQUIRE(path != nullptr);
    std::ifstream in(path);
    REQUIRE(in.good());
    ReadConsensusParams params;
    params.keep_homopolymer_lengths = false;
    const char* seed = std::getenv("VG_CONSENSUS_SEED");
    params.seed_from_reads = seed == nullptr || string(seed) != "0";

    vector<int64_t> engine;
    vector<int64_t> backbone_distance;
    vector<int64_t> engine_collapsed;
    vector<int64_t> backbone_collapsed;
    size_t failed = 0;
    size_t nearer_reads = 0;
    size_t as_near_reads = 0;
    double total_seconds = 0;
    double worst_seconds = 0;
    string line;
    std::cout << "case\tbackbone_len\treads\tabpoa_len\tengine_len\tengine\tbackbone\tengine_hp"
              << "\tbackbone_hp\tseconds\treads_to_engine\treads_to_abpoa" << std::endl;
    for (size_t number = 0; std::getline(in, line); ++number) {
        const string backbone = json_string(line, "backbone");
        const string abpoa = json_string(line, "abpoa");
        const vector<string> reads = json_strings(line, "reads");
        const auto start = std::chrono::steady_clock::now();
        const optional<string> consensus = read_consensus(reads, backbone, params);
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        total_seconds += seconds;
        worst_seconds = std::max(worst_seconds, seconds);
        if (!consensus) {
            ++failed;
            std::cout << number << "\t" << backbone.size() << "\t" << reads.size() << "\t"
                      << abpoa.size() << "\tnone" << std::endl;
            continue;
        }
        const string abpoa_collapsed = collapse_homopolymers(abpoa);
        engine.push_back(*edit_distance(*consensus, abpoa));
        backbone_distance.push_back(*edit_distance(backbone, abpoa));
        engine_collapsed.push_back(
            *edit_distance(collapse_homopolymers(*consensus), abpoa_collapsed));
        backbone_collapsed.push_back(
            *edit_distance(collapse_homopolymers(backbone), abpoa_collapsed));
        std::cout << number << "\t" << backbone.size() << "\t" << reads.size() << "\t"
                  << abpoa.size() << "\t" << consensus->size() << "\t" << engine.back() << "\t"
                  << backbone_distance.back() << "\t" << engine_collapsed.back() << "\t"
                  << backbone_collapsed.back() << "\t" << seconds;
        int64_t to_engine = 0;
        int64_t to_abpoa = 0;
        if (*consensus != abpoa) {
            for (const string& read : reads) {
                to_engine += *edit_distance(read, *consensus);
                to_abpoa += *edit_distance(read, abpoa);
            }
        }
        nearer_reads += to_engine < to_abpoa;
        as_near_reads += to_engine == to_abpoa;
        std::cout << "\t" << to_engine << "\t" << to_abpoa << std::endl;
    }
    REQUIRE(!engine.empty());
    auto summarise = [&](const string& name, const vector<int64_t>& ours,
                         const vector<int64_t>& theirs) {
        size_t exact = 0;
        size_t within5 = 0;
        size_t within20 = 0;
        int64_t worst_excess = std::numeric_limits<int64_t>::min();
        for (size_t i = 0; i < ours.size(); ++i) {
            exact += ours[i] == 0;
            within5 += ours[i] <= 5;
            within20 += ours[i] <= 20;
            worst_excess = std::max(worst_excess, ours[i] - theirs[i]);
        }
        std::cout << name << ": cases " << ours.size() << " exact " << exact << " within5 "
                  << within5 << " within20 " << within20 << " median " << median_of(ours)
                  << " worst " << *std::ranges::max_element(ours) << " backbone_median "
                  << median_of(theirs) << " worst_excess_over_backbone " << worst_excess
                  << std::endl;
    };
    summarise("raw", engine, backbone_distance);
    summarise("collapsed", engine_collapsed, backbone_collapsed);
    std::cout << "reads nearer the engine's consensus " << nearer_reads << " as near "
              << as_near_reads << " nearer the recorded " << engine.size() - nearer_reads
              - as_near_reads << std::endl;
    std::cout << "failed " << failed << " seconds " << total_seconds << " worst_seconds "
              << worst_seconds << std::endl;
}

}
}
