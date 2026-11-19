/** \file
 *
 * Unit tests for ParallelKFFReader, which reads kmers and their counts from a
 * KFF file in many threads at once.
 */

#include "../../src/kff.hpp"
#include "../../src/utility.hpp"

#include "catch.hpp"

#include <algorithm>
#include <array>
#include <random>

namespace vg::unittest {

//------------------------------------------------------------------------------

namespace {

using kmer_count_type = std::pair<ParallelKFFReader::kmer_type, size_t>;

// A block of one or more overlapping kmers: their sequence and the count of
// each kmer.
struct TestBlock {
    std::string sequence;
    std::vector<uint8_t> counts;
};

// A raw section of blocks and the variables written before it.
struct TestSection {
    size_t k;
    size_t max;
    std::vector<TestBlock> blocks;
};

// Write the sections to a KFF file with the given encoding, with one count
// byte per kmer.
void write_kff(const std::string& filename, const std::vector<TestSection>& sections, const std::array<uint8_t, 4>& encoding) {
    Kff_file file(filename, "w");
    file.write_encoding(encoding[0], encoding[1], encoding[2], encoding[3]);
    file.set_uniqueness(true);
    file.set_canonicity(false);
    for (const TestSection& section : sections) {
        Section_GV variables(&file);
        variables.write_var("k", section.k);
        variables.write_var("max", section.max);
        variables.write_var("data_size", 1);
        variables.close();
        Section_Raw raw(&file);
        for (const TestBlock& block : section.blocks) {
            std::vector<uint8_t> packed = kff_encode(block.sequence, encoding.data());
            std::vector<uint8_t> counts = block.counts;
            raw.write_compacted_sequence(packed.data(), block.sequence.length(), counts.data());
        }
        raw.close();
    }
    file.close();
}

// Return the kmers and the counts in the sections in minimizer index encoding,
// sorted.
std::vector<kmer_count_type> expected_kmers(const std::vector<TestSection>& sections) {
    std::vector<kmer_count_type> result;
    for (const TestSection& section : sections) {
        for (const TestBlock& block : section.blocks) {
            for (size_t i = 0; i + section.k <= block.sequence.length(); i++) {
                std::string kmer = block.sequence.substr(i, section.k);
                result.emplace_back(gbwtgraph::Key64::encode(kmer).get_key(), block.counts[i]);
            }
        }
    }
    std::ranges::sort(result);
    return result;
}

// Read the file in batches of `n` kmers with the given number of threads, and
// return the kmers and the counts, sorted. Each batch except the last must have
// `n` kmers.
std::vector<kmer_count_type> read_kmers(const std::string& filename, size_t n, size_t threads) {
    ParallelKFFReader reader(filename);
    std::vector<kmer_count_type> result;
    size_t short_batches = 0;
    #pragma omp parallel num_threads(threads)
    {
        while (true) {
            std::vector<kmer_count_type> batch = reader.read(n);
            #pragma omp critical
            {
                result.insert(result.end(), batch.begin(), batch.end());
                if (batch.size() < n) {
                    short_batches++;
                }
            }
            if (batch.empty()) {
                break;
            }
        }
    }
    // Each thread stops at a short batch, and only the last batch with kmers
    // can be short.
    REQUIRE(short_batches <= threads + 1);
    std::ranges::sort(result);
    return result;
}

// Return a random sequence of the given length.
std::string random_sequence(size_t length, std::mt19937& rng) {
    std::string result;
    for (size_t i = 0; i < length; i++) {
        result.push_back("ACGT"[rng() & 3]);
    }
    return result;
}

// Return sections with blocks of up to 3 kmers, then sections of single kmers,
// with k = 5.
std::vector<TestSection> mixed_sections() {
    std::mt19937 rng(0x5EED);
    std::vector<TestSection> result;
    TestSection multi { 5, 3, {} };
    multi.blocks.push_back({ "GATTACA", { 3, 1, 4 } });
    multi.blocks.push_back({ "CATTA", { 7 } });
    multi.blocks.push_back({ "TTACCA", { 2, 9 } });
    result.push_back(multi);
    for (size_t section = 0; section < 2; section++) {
        TestSection single { 5, 1, {} };
        for (size_t i = 0; i < 25; i++) {
            single.blocks.push_back({ random_sequence(5, rng), { static_cast<uint8_t>(rng() & 0xFF) } });
        }
        result.push_back(single);
    }
    return result;
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("ParallelKFFReader reads the kmers and the counts in the file", "[kff]") {
    std::vector<TestSection> sections = mixed_sections();
    std::vector<kmer_count_type> truth = expected_kmers(sections);
    std::array<uint8_t, 4> kff_encoding { 0, 1, 3, 2 };
    std::array<uint8_t, 4> trivial_encoding { 0, 1, 2, 3 };
    for (const std::array<uint8_t, 4>& encoding : { kff_encoding, trivial_encoding }) {
        std::string filename = temp_file::create("parallel-kff-reader");
        write_kff(filename, sections, encoding);
        for (size_t threads : { 1, 4 }) {
            for (size_t n : { 1, 2, 4, 7, 1000 }) {
                std::vector<kmer_count_type> kmers = read_kmers(filename, n, threads);
                REQUIRE(kmers == truth);
            }
        }
        temp_file::remove(filename);
    }
}

TEST_CASE("ParallelKFFReader takes the kmer length from the first section", "[kff]") {
    std::vector<TestSection> sections = mixed_sections();
    std::string filename = temp_file::create("parallel-kff-reader");
    std::array<uint8_t, 4> encoding { 0, 1, 3, 2 };

    SECTION("a section of shorter kmers") {
        sections.push_back({ 4, 1, { { "GATT", { 1 } } } });
        write_kff(filename, sections, encoding);
        ParallelKFFReader reader(filename);
        REQUIRE(reader.k == 5);
        REQUIRE_THROWS_AS(reader.read(1000), std::runtime_error);
    }

    SECTION("a section of longer blocks") {
        TestSection longer { 5, 4, {} };
        longer.blocks.push_back({ "GATTACAT", { 1, 2, 3, 4 } });
        sections.push_back(longer);
        write_kff(filename, sections, encoding);
        ParallelKFFReader reader(filename);
        REQUIRE(reader.max_kmers_per_block == 3);
        REQUIRE_THROWS_AS(reader.read(1000), std::runtime_error);
    }

    temp_file::remove(filename);
}

//------------------------------------------------------------------------------

} // namespace vg::unittest
