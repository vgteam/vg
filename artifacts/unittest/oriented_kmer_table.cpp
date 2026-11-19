/** \file
 *
 * Unit tests for OrientedKmerTable and OrientedKmerCounts in recombinator.cpp. The table is
 * checked against a map filled from the kmer list in order, and the counts against a count
 * of every window of k bases in the pieces of a sequence.
 */

#include "../../src/recombinator.hpp"
#include "../../src/utility.hpp"

#include "catch.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <random>
#include <set>
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

// Return the pieces of a sequence with `copies` copies of listed kmers, in either
// orientation, among random bases. Some copies follow a run of N, some have an N inside,
// and some are cut by the end of a piece.
std::vector<std::string> test_pieces(std::mt19937_64& rng, const std::vector<kmer_type>& kmers, size_t k, size_t copies) {
    std::string text;
    std::vector<size_t> cuts;
    for (size_t i = 0; i < copies; i++) {
        text += random_bases(rng, rng() % 8);
        std::string kmer = gbwtgraph::Key64(kmers[rng() % kmers.size()]).decode(k);
        if (rng() & 1) {
            kmer = reverse_complement(kmer);
        }
        switch (rng() % 6) {
        case 0:
            text += std::string(1 + rng() % 4, 'N');
            break;
        case 1:
            kmer[rng() % k] = 'N';
            break;
        case 2:
            cuts.push_back(text.size() + 1 + rng() % (k - 1));
            break;
        default:
            break;
        }
        text += kmer;
    }
    cuts.push_back(text.size());
    std::sort(cuts.begin(), cuts.end());

    // Cuts at the same place leave an empty piece.
    std::vector<std::string> result;
    size_t from = 0;
    for (size_t cut : cuts) {
        result.push_back(text.substr(from, cut - from));
        from = cut;
    }
    return result;
}

// The counts of the listed kmers in a sequence, found by checking every window of k bases.
struct ExpectedCounts {
    std::vector<std::uint32_t> counts;

    // (position, occurrences) for the kmers that occur, in the order of first occurrence.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> occurrences;

    // Windows that are the reverse complement of a listed kmer and not a listed kmer.
    size_t reverse_only = 0;

    size_t total() const {
        size_t result = 0;
        for (std::uint32_t count : this->counts) {
            result += count;
        }
        return result;
    }
};

// Count the listed kmers in the pieces by checking every window of k bases.
ExpectedCounts count_windows(
    const std::vector<std::string>& pieces, const std::vector<kmer_type>& kmers,
    const std::map<kmer_type, std::uint32_t>& positions, size_t k
) {
    ExpectedCounts result;
    result.counts.resize(kmers.size(), 0);
    std::set<kmer_type> listed(kmers.begin(), kmers.end());
    std::vector<std::uint32_t> order;
    for (const std::string& piece : pieces) {
        for (size_t start = 0; start + k <= piece.size(); start++) {
            std::string window = piece.substr(start, k);
            if (window.find_first_not_of("ACGT") != std::string::npos) {
                continue;
            }
            kmer_type key = encode(window);
            auto iter = positions.find(key);
            if (iter == positions.end()) {
                continue;
            }
            result.reverse_only += (listed.contains(key) ? 0 : 1);
            if (result.counts[iter->second]++ == 0) {
                order.push_back(iter->second);
            }
        }
    }
    for (std::uint32_t position : order) {
        result.occurrences.emplace_back(position, result.counts[position]);
    }
    return result;
}

// Check the counts in the object against the counts of the windows.
void check_counts(const OrientedKmerCounts& counts, const ExpectedCounts& expected) {
    for (size_t position = 0; position < expected.counts.size(); position++) {
        REQUIRE(counts[position] == expected.counts[position]);
    }
    REQUIRE(counts.occurrence_list() == expected.occurrences);
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

TEST_CASE("Oriented kmer counts match a count of every window", "[haplotypes][oriented_kmer_table]") {
    for (size_t k : { size_t(29), size_t(30), size_t(31) }) {
        DYNAMIC_SECTION("kmers of length " << k) {
            std::mt19937_64 rng(100 + k);
            std::vector<kmer_type> kmers = test_kmers(rng, k, 300);
            std::map<kmer_type, std::uint32_t> positions = expected_positions(kmers, k);
            OrientedKmerTable table(kmers, k);

            std::vector<std::string> first = test_pieces(rng, kmers, k, 60);
            std::vector<std::string> second = test_pieces(rng, kmers, k, 60);
            ExpectedCounts first_expected = count_windows(first, kmers, positions, k);
            ExpectedCounts second_expected = count_windows(second, kmers, positions, k);

            // The pieces have kmers in both orientations, runs of N, and kmers cut by a
            // piece boundary, which the joined pieces would count.
            auto check_pieces = [&](const std::vector<std::string>& pieces, const ExpectedCounts& expected) {
                REQUIRE(expected.occurrences.size() > 1);
                REQUIRE(expected.reverse_only > 0);
                REQUIRE(std::any_of(pieces.begin(), pieces.end(), [](const std::string& piece) {
                    return piece.find('N') != std::string::npos;
                }));
                std::string joined;
                for (const std::string& piece : pieces) {
                    joined += piece;
                }
                REQUIRE(count_windows({ joined }, kmers, positions, k).total() > expected.total());
            };
            check_pieces(first, first_expected);
            check_pieces(second, second_expected);
            // Some kmer of the first sequence is missing from the second, so stale counts
            // would show.
            bool first_only = false;
            for (size_t position = 0; position < kmers.size(); position++) {
                first_only |= (first_expected.counts[position] > 0 && second_expected.counts[position] == 0);
            }
            REQUIRE(first_only);

            // One object counts the sequences in turn.
            OrientedKmerCounts counts(kmers.size());
            counts.count(first, table, k);
            check_counts(counts, first_expected);
            counts.count(second, table, k);
            check_counts(counts, second_expected);
            counts.count({}, table, k);
            check_counts(counts, count_windows({}, kmers, positions, k));
            counts.count(first, table, k);
            check_counts(counts, first_expected);
        }
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
