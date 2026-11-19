/** \file
 *
 * Unit tests for DiploidPairScorer, which scores pairs of haplotype sequences in
 * a subchain as the two haplotypes of a diploid sample.
 */

#include "../diploid_pair_scorer.hpp"

#include "catch.hpp"

#include <random>
#include <string>
#include <vector>

namespace vg::unittest {

//------------------------------------------------------------------------------

namespace {

using kmer_types_type = std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>;

// Return a subchain with `clusters * size` sequences over `num_kmers` kmers.
// Each cluster holds a random row and copies of it with a few bits flipped, so
// the compressed presence matrix stores most sequences as differences from a
// parent.
Haplotypes::Subchain random_subchain(size_t num_kmers, size_t clusters, size_t size, std::mt19937& rng) {
    Haplotypes::Subchain result;
    result.type = Haplotypes::Subchain::normal;
    result.start = gbwt::Node::encode(1, false);
    result.end = gbwt::Node::encode(2, false);
    for (size_t kmer_id = 0; kmer_id < num_kmers; kmer_id++) {
        result.kmers.push_back(kmer_id + 1);
    }
    result.kmer_counts = sdsl::int_vector<0>(num_kmers, 0, 64);

    std::vector<KmerPresenceMatrix::sequence_type> sequences;
    sdsl::bit_vector matrix(clusters * size * num_kmers, 0);
    for (size_t cluster = 0; cluster < clusters; cluster++) {
        std::vector<bool> parent(num_kmers);
        for (size_t kmer_id = 0; kmer_id < num_kmers; kmer_id++) {
            parent[kmer_id] = rng() & 1;
        }
        for (size_t member = 0; member < size; member++) {
            size_t row = cluster * size + member;
            sequences.push_back({ row, 0 });
            for (size_t kmer_id = 0; kmer_id < num_kmers; kmer_id++) {
                bool flipped = (member > 0 && rng() % 16 == 0);
                matrix[row * num_kmers + kmer_id] = (parent[kmer_id] != flipped);
            }
        }
    }
    result.kmers_present = KmerPresenceMatrix(sequences, result.sequences, num_kmers, matrix);
    return result;
}

// Return a random class for each kmer, with every class present.
kmer_types_type random_kmer_types(size_t num_kmers, std::mt19937& rng) {
    kmer_types_type result;
    for (size_t kmer_id = 0; kmer_id < num_kmers; kmer_id++) {
        auto type = static_cast<KmerPresenceMatrix::kmer_classification>(rng() % 5);
        result.emplace_back(type, 1.0);
    }
    return result;
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("DiploidPairScorer agrees with the presence matrix", "[haplotypes]") {
    std::mt19937 rng(0xD1F);
    // Several trials of each length, so that the kmers next to the word
    // boundaries take every class.
    std::vector<size_t> lengths;
    for (size_t num_kmers : { 0, 1, 63, 64, 65, 127, 128, 130, 1000 }) {
        lengths.insert(lengths.end(), 8, num_kmers);
    }
    for (size_t num_kmers : lengths) {
        Haplotypes::Subchain subchain = random_subchain(num_kmers, 4, 5, rng);
        if (num_kmers >= 64) {
            // Most sequences are stored as differences from a parent.
            REQUIRE(subchain.kmers_present.get_num_parents() < subchain.num_sequences());
        }
        kmer_types_type kmer_types = random_kmer_types(num_kmers, rng);

        // Candidates in an arbitrary order, with one sequence listed twice.
        std::vector<size_t> candidates;
        for (size_t i = 0; i < subchain.num_sequences(); i += 2) {
            candidates.push_back(subchain.num_sequences() - 1 - i);
        }
        candidates.push_back(candidates.front());

        DiploidPairScorer scorer(subchain, candidates, kmer_types);
        for (size_t first = 0; first < candidates.size(); first++) {
            for (size_t second = 0; second < candidates.size(); second++) {
                std::int64_t expected = subchain.diploid_score(candidates[first], candidates[second], kmer_types);
                REQUIRE(scorer.score(first, second) == expected);
            }
        }
    }
}

//------------------------------------------------------------------------------

} // namespace vg::unittest
