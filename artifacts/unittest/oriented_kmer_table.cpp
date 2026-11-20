/** \file
 *
 * Unit tests for OrientedKmerTable and OrientedKmerCounts in subchain_walk.cpp. The table
 * is checked against a map filled from the kmer list in order, and the counts against
 * occurrence lists written out in the test.
 */

#include "../../src/recombinator.hpp"
#include "../../src/utility.hpp"

#include "catch.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace vg::unittest {

namespace {

using kmer_type = Haplotypes::Subchain::kmer_type;

// Return a string of random bases.
std::string random_bases(std::mt19937_64& rng, size_t length) {
    static const char BASES[] = "ACGT";
    std::string result(length, 'A');
    for (char& c : result) {
        c = BASES[rng() & 3];
    }
    return result;
}

kmer_type encode(const std::string& kmer) {
    return gbwtgraph::Key64::encode(kmer).get_key();
}

kmer_type reverse_complement_key(kmer_type kmer, size_t k) {
    return gbwtgraph::Key64(kmer).reverse_complement(k).get_key();
}

// Return a list of kmers of length k. The first `n` are the kmer of all A and random
// kmers. The reverse complements of every tenth of them follow, and then, for even k, five
// palindromes.
std::vector<kmer_type> test_kmers(std::mt19937_64& rng, size_t k, size_t n) {
    std::vector<kmer_type> result { encode(std::string(k, 'A')) };
    while (result.size() < n) {
        result.push_back(encode(random_bases(rng, k)));
    }
    for (size_t i = 0; i < n; i += 10) {
        result.push_back(reverse_complement_key(result[i], k));
    }
    if (k % 2 == 0) {
        for (size_t i = 0; i < 5; i++) {
            std::string half = random_bases(rng, k / 2);
            result.push_back(encode(half + reverse_complement(half)));
        }
    }
    return result;
}

// Return the position of each kmer of the list and of its reverse complement. The map is
// filled in list order, so a later position replaces an earlier one.
std::map<kmer_type, std::uint32_t> expected_positions(const std::vector<kmer_type>& kmers, size_t k) {
    std::map<kmer_type, std::uint32_t> result;
    for (size_t i = 0; i < kmers.size(); i++) {
        result[kmers[i]] = static_cast<std::uint32_t>(i);
        result[reverse_complement_key(kmers[i], k)] = static_cast<std::uint32_t>(i);
    }
    return result;
}

// Check the table against the map for each kmer in the map and for random kmers.
void check_table(const OrientedKmerTable& table, const std::map<kmer_type, std::uint32_t>& expected, std::mt19937_64& rng, size_t k) {
    for (const auto& [kmer, position] : expected) {
        REQUIRE(table.find(kmer) == position);
    }
    size_t absent = 0;
    for (size_t i = 0; i < 1000; i++) {
        kmer_type kmer = encode(random_bases(rng, k));
        auto iter = expected.find(kmer);
        std::uint32_t position = (iter == expected.end() ? OrientedKmerTable::NONE : iter->second);
        absent += (position == OrientedKmerTable::NONE);
        REQUIRE(table.find(kmer) == position);
    }
    REQUIRE(absent > 0);
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("An oriented kmer table finds each kmer in both orientations at its last position", "[haplotypes][oriented_kmer_table]") {
    for (size_t k : { size_t(29), size_t(30), size_t(31) }) {
        DYNAMIC_SECTION("kmers of length " << k) {
            std::mt19937_64 rng(k);
            // Hundreds of kmers make keys share slots in the table.
            for (size_t n : { size_t(1), size_t(2), size_t(3), size_t(300) }) {
                INFO(n << " random kmers");
                std::vector<kmer_type> kmers = test_kmers(rng, k, n);
                std::map<kmer_type, std::uint32_t> expected = expected_positions(kmers, k);
                // The reverse complement of the kmer of all A comes later in the list.
                REQUIRE(expected.at(kmers.front()) == n);
                OrientedKmerTable table(kmers, k);
                check_table(table, expected, rng, k);
            }
        }
    }

    SECTION("the key reverse complement agrees with the string reverse complement") {
        std::mt19937_64 rng(1);
        for (size_t k : { size_t(29), size_t(30), size_t(31) }) {
            std::string kmer = random_bases(rng, k);
            REQUIRE(reverse_complement_key(encode(kmer), k) == encode(reverse_complement(kmer)));
        }
    }

    SECTION("an empty list has no kmers") {
        std::mt19937_64 rng(2);
        OrientedKmerTable table({}, 31);
        REQUIRE(table.find(encode(std::string(31, 'A'))) == OrientedKmerTable::NONE);
        check_table(table, {}, rng, 31);
    }
}

TEST_CASE("Oriented kmer counts count the kmers of one sequence at a time", "[haplotypes][oriented_kmer_table]") {
    typedef std::vector<std::pair<std::uint32_t, std::uint32_t>> occurrence_list;

    SECTION("kmers in both orientations, in pieces, and in consecutive sequences") {
        size_t k = 29;
        std::mt19937_64 rng(100);
        std::string a = random_bases(rng, k);
        std::string b = random_bases(rng, k);
        std::string c = random_bases(rng, k);
        std::vector<kmer_type> kmers { encode(a), encode(b), encode(c) };
        OrientedKmerTable table(kmers, k);
        OrientedKmerCounts counts(kmers.size());

        // Occurrences are listed in the order of the first one, and a reverse
        // complement counts for its kmer.
        counts.count({ b + "T" + reverse_complement(a) + b }, table, k);
        REQUIRE(counts.occurrence_list() == occurrence_list({ { 1, 2 }, { 0, 1 } }));
        REQUIRE(counts[2] == 0);

        // A kmer cut by the end of a piece or holding an N is not counted.
        std::string with_n = a;
        with_n[k / 2] = 'N';
        counts.count({ a.substr(0, k - 1), a.substr(k - 1) + with_n, "N" + a }, table, k);
        REQUIRE(counts.occurrence_list() == occurrence_list({ { 0, 1 } }));

        // Counting a sequence clears the counts of the previous one.
        counts.count({ a + c }, table, k);
        counts.count({ b }, table, k);
        REQUIRE(counts.occurrence_list() == occurrence_list({ { 1, 1 } }));
        REQUIRE(counts[0] == 0);
        counts.count({}, table, k);
        REQUIRE(counts.occurrence_list().empty());
    }

    SECTION("a kmer after a run of N counts once it has k bases") {
        size_t k = 29;
        std::string kmer(k, 'A');
        std::vector<kmer_type> kmers { encode(kmer) };
        OrientedKmerTable table(kmers, k);
        OrientedKmerCounts counts(kmers.size());
        // The first piece has k - 1 bases after the N, and the second has k + 1.
        std::vector<std::string> pieces { "N" + kmer.substr(1), "CCN" + kmer + "A" };
        counts.count(pieces, table, k);
        REQUIRE(counts[0] == 2);
        REQUIRE(counts.occurrence_list() == std::vector<std::pair<std::uint32_t, std::uint32_t>>({ { 0, 2 } }));
        // The reverse complement of the kmer counts for it too.
        counts.count({ std::string(k, 'T') }, table, k);
        REQUIRE(counts[0] == 1);
    }

    SECTION("an empty list counts nothing") {
        std::mt19937_64 rng(3);
        OrientedKmerTable table({}, 31);
        OrientedKmerCounts counts(0);
        counts.count({ random_bases(rng, 200), "", random_bases(rng, 50) }, table, 31);
        REQUIRE(counts.occurrence_list().empty());
    }
}

} // namespace vg::unittest
