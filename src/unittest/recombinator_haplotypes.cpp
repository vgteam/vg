/** \file
 *
 * Unit tests for `KmerPresenceMatrix` and for the kmer presence iterators in
 * `Haplotypes::Subchain`.
 */

#include "catch.hpp"

#include "../recombinator_haplotypes.hpp"

#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace vg {
namespace unittest {

//------------------------------------------------------------------------------

namespace {

// Converts a matrix of `0` and `1` characters into a concatenated bitvector.
sdsl::bit_vector build_bitvector(const std::vector<std::string>& rows, size_t num_kmers) {
    sdsl::bit_vector result(rows.size() * num_kmers, 0);
    for (size_t i = 0; i < rows.size(); i++) {
        REQUIRE(rows[i].length() == num_kmers);
        for (size_t j = 0; j < num_kmers; j++) {
            if (rows[i][j] == '1') {
                result[i * num_kmers + j] = 1;
            }
        }
    }
    return result;
}

// Sequence identifiers for the compressed constructor. Sequence `i` is `(i, i + 1)`.
std::vector<KmerPresenceMatrix::sequence_type> sequence_ids(size_t num_sequences) {
    std::vector<KmerPresenceMatrix::sequence_type> result;
    for (size_t i = 0; i < num_sequences; i++) {
        result.push_back({ i, i + 1 });
    }
    return result;
}

// Returns the rows in the order given by the reordered sequences, assuming that
// the sequences were identified using `sequence_ids`. Checks that the reordering
// is a permutation of the original sequences.
std::vector<std::string> reorder_rows(
    const std::vector<std::string>& rows,
    const std::vector<KmerPresenceMatrix::compact_sequence_type>& reordered
) {
    REQUIRE(reordered.size() == rows.size());
    std::vector<bool> seen(rows.size(), false);
    std::vector<std::string> result;
    for (auto sequence : reordered) {
        REQUIRE(sequence.first < rows.size());
        REQUIRE(sequence.second == sequence.first + 1);
        REQUIRE(!seen[sequence.first]);
        seen[sequence.first] = true;
        result.push_back(rows[sequence.first]);
    }
    return result;
}

// TODO: Later tests may require the other fields as well.
// Builds a `Subchain` with the given kmer presence matrix, given as a vector of
// strings of '0' and '1' characters. Uses dummy values for the other fields.
Haplotypes::Subchain build_subchain(const std::vector<std::string>& rows, size_t kmer_count) {
    Haplotypes::Subchain result;
    result.type = Haplotypes::Subchain::normal;
    result.start = gbwt::Node::encode(1, false);
    result.end = gbwt::Node::encode(2, false);

    for (size_t kmer_id = 0; kmer_id < kmer_count; kmer_id++) {
        result.kmers.push_back(kmer_id + 1);
    }
    result.kmer_counts = sdsl::int_vector<0>(kmer_count, 0, 64);
    std::vector<KmerPresenceMatrix::sequence_type> sequences;
    for (size_t i = 0; i < rows.size(); i++) {
        sequences.push_back({ i, 0 });
    }

    // Sequence `i` in the subchain is row `result.sequences[i].first`.
    sdsl::bit_vector buffer = build_bitvector(rows, kmer_count);
    result.kmers_present = KmerPresenceMatrix(sequences, result.sequences, kmer_count, buffer);

    return result;
}

// Generates a random matrix of `0` and `1` values. There are `num_kmers`
// columns and `clusters * children` rows. Each cluster consists of `children`
// rows with `mutations` random mutations relative to the parent, which is
// not in the matrix.
std::vector<std::string> clustered_rows(
    size_t clusters, size_t children, size_t num_kmers, size_t mutations, std::uint64_t seed
) {
    std::mt19937_64 rng(seed);
    std::vector<std::string> result;
    for (size_t cluster = 0; cluster < clusters; cluster++) {
        std::string parent(num_kmers, '0');
        for (size_t kmer_id = 0; kmer_id < num_kmers; kmer_id++) {
            parent[kmer_id] = (rng() % 2 == 0 ? '0' : '1');
        }
        for (size_t child = 0; child < children; child++) {
            std::string sequence = parent;
            for (size_t i = 0; i < mutations; i++) {
                size_t kmer_id = rng() % num_kmers;
                sequence[kmer_id] = (sequence[kmer_id] == '0' ? '1' : '0');
            }
            result.push_back(sequence);
        }
    }

    std::shuffle(result.begin(), result.end(), rng);
    return result;
}

// Gives an arbitrary classification for a kmer based on its index.
KmerPresenceMatrix::kmer_classification kmer_class(size_t kmer_id) {
    switch (kmer_id % 4) {
        case 0: return KmerPresenceMatrix::kmer_classification::absent;
        case 1: return KmerPresenceMatrix::kmer_classification::heterozygous;
        case 2: return KmerPresenceMatrix::kmer_classification::present;
        case 3: return KmerPresenceMatrix::kmer_classification::frequent;
    }
    return KmerPresenceMatrix::kmer_classification::absent;
}

// Arbitrary integer weights for the kmers. Integers are exact in `double`, so we can
// compare the scores for equality even though the matrix may sum them in another order.
double kmer_weight(size_t kmer_id) {
    return static_cast<double>(kmer_id % 7) - 3.0;
}

// Checks `KmerPresenceMatrix` queries and iterators agains the original rows.
void check_queries(
    const std::string& name, const KmerPresenceMatrix& matrix,
    const std::vector<std::string>& rows, size_t num_kmers
) {
    REQUIRE(matrix.get_num_sequences() == rows.size());
    REQUIRE(matrix.get_num_kmers() == num_kmers);
    REQUIRE(matrix.size() == rows.size() * num_kmers);
    REQUIRE(matrix.empty() == (rows.size() * num_kmers == 0));
    REQUIRE(matrix.get_num_parents() <= rows.size());

    // Number of present kmers in each sequence, and in total.
    size_t total_present = 0;
    for (size_t i = 0; i < rows.size(); i++) {
        size_t expected = 0;
        for (size_t j = 0; j < num_kmers; j++) {
            expected += (rows[i][j] == '1');
        }
        if (matrix.get_num_present(i) != expected) {
            std::cerr << name << ": wrong number of present kmers in sequence " << i << std::endl;
        }
        REQUIRE(matrix.get_num_present(i) == expected);
        total_present += expected;
    }
    REQUIRE(matrix.get_total_present() == total_present);
    REQUIRE(matrix.get_num_present(rows.size()) == 0);

    // Iterating over the kmers in a single sequence.
    for (size_t i = 0; i < rows.size(); i++) {
        std::vector<bool> found;
        matrix.for_each_kmer(i, [&](size_t kmer_id, bool is_present) {
            REQUIRE(kmer_id == found.size());
            found.push_back(is_present);
        });
        if (found.size() != num_kmers) {
            std::cerr << name << ": wrong number of kmers in sequence " << i << std::endl;
        }
        REQUIRE(found.size() == num_kmers);
        for (size_t j = 0; j < num_kmers; j++) {
            if (found[j] != (rows[i][j] == '1')) {
                std::cerr << name << ": wrong value for kmer " << j << " in sequence " << i << std::endl;
            }
            REQUIRE(found[j] == (rows[i][j] == '1'));
        }
    }

    // Out-of-range sequences do not produce any calls.
    {
        size_t calls = 0;
        matrix.for_each_kmer(rows.size(), [&](size_t, bool) { calls++; });
        REQUIRE(calls == 0);
    }

    // Kmers with arbitrary classifications and scores.
    std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>> kmers;
    for (size_t j = 0; j < num_kmers; j++) {
        kmers.emplace_back(kmer_class(j), kmer_weight(j));
    }

    // Iterating over pairs of sequences and computing their diploid scores.
    for (size_t i = 0; i < rows.size(); i++) {
        for (size_t j = i; j < rows.size(); j++) {
            std::int64_t score = matrix.diploid_score(i, j, kmers);
            std::int64_t expected_score = 0;
            for (size_t kmer_id = 0; kmer_id < num_kmers; kmer_id++) {
                std::int64_t num_present = (rows[i][kmer_id] == '1') + (rows[j][kmer_id] == '1');
                if (kmers[kmer_id].first == KmerPresenceMatrix::absent) {
                    expected_score += 1 - num_present;
                } else if (kmers[kmer_id].first == KmerPresenceMatrix::heterozygous) {
                    if (num_present == 1) {
                        expected_score += 1;
                    }
                } else if (kmers[kmer_id].first == KmerPresenceMatrix::present) {
                    expected_score += num_present - 1;
                }
            }
            REQUIRE(score == expected_score);
        }
    }

    // Scoring the sequences.
    {
        std::vector<size_t> reported(rows.size(), 0);
        std::vector<double> scores(rows.size(), 0.0);
        matrix.score_sequences(
            kmers,
            [&](size_t sequence, double score) {
                REQUIRE(sequence < rows.size());
                reported[sequence]++;
                scores[sequence] = score;
            }
        );
        for (size_t i = 0; i < rows.size(); i++) {
            REQUIRE(reported[i] == 1);
            double expected = 0.0;
            for (size_t j = 0; j < num_kmers; j++) {
                expected += (rows[i][j] == '1' ? 1.0 : -1.0) * kmer_weight(j);
            }
            if (scores[i] != expected) {
                std::cerr << name << ": wrong score for sequence " << i << std::endl;
            }
            REQUIRE(scores[i] == expected);
        }
    }
}

// Builds and returns a `KmerPresenceMatrix` from the rows. Checks queries, iterators,
// and serialization. A compressed matrix reorders the sequences, and the queries are
// checked against the rows in the stored order.
KmerPresenceMatrix check_matrix(
    const std::string& name, const std::vector<std::string>& rows, size_t num_kmers, bool compress
) {
    KmerPresenceMatrix matrix;
    if (compress) {
        std::vector<KmerPresenceMatrix::compact_sequence_type> reordered;
        matrix = KmerPresenceMatrix(sequence_ids(rows.size()), reordered, num_kmers, build_bitvector(rows, num_kmers));
        check_queries(name, matrix, reorder_rows(rows, reordered), num_kmers);
    } else {
        matrix = KmerPresenceMatrix(rows.size(), num_kmers, build_bitvector(rows, num_kmers));
        check_queries(name, matrix, rows, num_kmers);
    }

    std::stringstream stream;
    matrix.simple_sds_serialize(stream);
    REQUIRE(stream.str().length() == 8 * matrix.simple_sds_size());

    KmerPresenceMatrix loaded;
    loaded.simple_sds_load(stream);
    bool success = (loaded == matrix);
    if (!success) {
        std::cerr << name << ": loaded matrix does not match original" << std::endl;
    }
    REQUIRE(success);

    return matrix;
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Empty kmer presence matrix", "[haplotypes]") {
    SECTION("default constructor") {
        KmerPresenceMatrix matrix;
        check_queries("default", matrix, {}, 0);
        REQUIRE(matrix.get_num_parents() == 0);
    }

    SECTION("no sequences, no kmers") {
        std::string base_name = "no sequences, no kmers";
        for (bool compress : { false, true }) {
            std::string name = base_name + (compress ? " (compressed)" : " (uncompressed)");
            KmerPresenceMatrix matrix = check_matrix(name, {}, 0, compress);
            REQUIRE(matrix.get_num_parents() == 0);
        }
    }

    SECTION("no sequences") {
        std::string base_name = "no sequences";
        for (bool compress : { false, true }) {
            std::string name = base_name + (compress ? " (compressed)" : " (uncompressed)");
            // The number of kmers is normalized to 0, as it cannot be recovered.
            KmerPresenceMatrix matrix;
            if (compress) {
                std::vector<KmerPresenceMatrix::compact_sequence_type> reordered(1);
                matrix = KmerPresenceMatrix(sequence_ids(0), reordered, 5, sdsl::bit_vector());
                REQUIRE(reordered.empty());
            } else {
                matrix = KmerPresenceMatrix(0, 5, sdsl::bit_vector());
            }
            REQUIRE(matrix.get_num_kmers() == 0);
            check_queries(name, matrix, {}, 0);
        }
    }

    SECTION("no kmers") {
        std::string base_name = "no kmers";
        for (bool compress : { false, true }) {
            std::string name = base_name + (compress ? " (compressed)" : " (uncompressed)");
            KmerPresenceMatrix matrix = check_matrix(name, { "", "", "" }, 0, compress);
            // Without kmers, all sequences are identical and form a single cluster.
            REQUIRE(matrix.get_num_parents() == (compress ? 1 : 3));
        }
    }
}

TEST_CASE("Uncompressed kmer presence matrix", "[haplotypes]") {
    std::vector<std::string> rows {
        "00000",
        "11111",
        "10110",
        "01001"
    };

    SECTION("every sequence is a parent") {
        KmerPresenceMatrix matrix = check_matrix("uncompressed", rows, 5, false);
        REQUIRE(matrix.get_num_parents() == rows.size());
    }

    SECTION("a single sequence") {
        KmerPresenceMatrix matrix = check_matrix("single sequence", { "10110" }, 5, false);
        REQUIRE(matrix.get_num_parents() == 1);
    }

    SECTION("a single kmer") {
        KmerPresenceMatrix matrix = check_matrix("single kmer", { "0", "1", "1" }, 1, false);
        REQUIRE(matrix.get_num_parents() == 3);
    }
}

TEST_CASE("Compressed kmer presence matrix", "[haplotypes]") {
    SECTION("identical sequences") {
        std::vector<std::string> rows(8, std::string(128, '0'));
        for (size_t i = 0; i < 128; i += 3) {
            for (std::string& row : rows) { row[i] = '1'; }
        }
        KmerPresenceMatrix matrix = check_matrix("identical", rows, 128, true);
        REQUIRE(matrix.get_num_parents() == 1);

        KmerPresenceMatrix uncompressed(rows.size(), 128, build_bitvector(rows, 128));
        REQUIRE(matrix.simple_sds_size() < uncompressed.simple_sds_size());
    }

    SECTION("clustered sequences") {
        std::vector<std::string> rows = clustered_rows(3, 10, 256, 3, 0xdeadbeef);
        KmerPresenceMatrix matrix = check_matrix("clustered", rows, 256, true);
        REQUIRE(matrix.get_num_parents() == 3);

        KmerPresenceMatrix uncompressed(rows.size(), 256, build_bitvector(rows, 256));
        REQUIRE(matrix.simple_sds_size() < uncompressed.simple_sds_size());
    }

    SECTION("distinct sequences") {
        // Random rows differ in around half of the kmers, which is far too many.
        std::vector<std::string> rows = clustered_rows(16, 1, 64, 0, 0x12345678);
        KmerPresenceMatrix matrix = check_matrix("distinct", rows, 64, true);
        REQUIRE(matrix.get_num_parents() == rows.size());
    }

    SECTION("a single sequence") {
        KmerPresenceMatrix matrix = check_matrix("single sequence", { "10110" }, 5, true);
        REQUIRE(matrix.get_num_parents() == 1);
    }

    SECTION("a single kmer") {
        // With one kmer, only identical sequences can be clustered together.
        KmerPresenceMatrix matrix = check_matrix("single kmer", { "0", "1", "1" }, 1, true);
        REQUIRE(matrix.get_num_parents() == 2);
    }

    SECTION("duplicates of distinct sequences") {
        std::vector<std::string> distinct = clustered_rows(4, 1, 64, 0, 0x87654321);
        std::vector<std::string> rows;
        for (size_t i = 0; i < 3; i++) {
            rows.insert(rows.end(), distinct.begin(), distinct.end());
        }
        KmerPresenceMatrix matrix = check_matrix("duplicates", rows, 64, true);
        REQUIRE(matrix.get_num_parents() == distinct.size());
    }
}

TEST_CASE("Kmer presence matrix dimensions", "[haplotypes]") {
    // These exercise the 64-bit block boundaries in the queries and in the compression.
    std::vector<size_t> kmer_counts { 1, 2, 63, 64, 65, 127, 128, 129, 200 };

    for (size_t num_kmers : kmer_counts) {
        std::vector<std::string> rows = clustered_rows(2, 3, num_kmers, 1, num_kmers);
        for (bool compress : { false, true }) {
            std::string name = "kmers " + std::to_string(num_kmers) + (compress ? " (compressed)" : "");
            KmerPresenceMatrix matrix = check_matrix(name, rows, num_kmers, compress);
            REQUIRE(matrix.get_num_parents() >= 1);
            REQUIRE(matrix.get_num_parents() <= rows.size());
        }
    }
}

TEST_CASE("Invalid kmer presence matrix", "[haplotypes]") {
    SECTION("wrong matrix size (uncompressed)") {
        REQUIRE_THROWS_AS(
            KmerPresenceMatrix(3, 5, sdsl::bit_vector(14, 0)),
            std::runtime_error
        );
    }

    SECTION("wrong matrix size (compressed)") {
        std::vector<KmerPresenceMatrix::compact_sequence_type> reordered;
        REQUIRE_THROWS_AS(
            KmerPresenceMatrix(sequence_ids(3), reordered, 5, sdsl::bit_vector(14, 0)),
            std::runtime_error
        );
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Iterating over the kmers in a sequence", "[haplotypes]") {
    // The last two rows become children of the first two, which reorders the sequences.
    std::vector<std::string> matrix {
        "00000",
        "11111",
        "10110",
        "01001",
        "10111",
        "01000"
    };
    Haplotypes::Subchain subchain = build_subchain(matrix, 5);

    SECTION("every kmer is visited once in order") {
        REQUIRE(subchain.sequences.size() == matrix.size());
        for (size_t i = 0; i < matrix.size(); i++) {
            const std::string& row = matrix[subchain.sequences[i].first];
            std::vector<bool> found;
            subchain.for_each_kmer(i, [&](size_t kmer_id, bool is_present) {
                REQUIRE(kmer_id == found.size());
                found.push_back(is_present);
            });
            REQUIRE(found.size() == subchain.kmers.size());
            for (size_t j = 0; j < found.size(); j++) {
                bool should_be_present = (row[j] == '1');
                REQUIRE(found[j] == should_be_present);
            }
        }
    }

    SECTION("no kmers") {
        Haplotypes::Subchain empty = build_subchain({ "", "", "" }, 0);
        for (size_t i = 0; i < 3; i++) {
            size_t calls = 0;
            empty.for_each_kmer(i, [&](size_t, bool) { calls++; });
            REQUIRE(calls == 0);
        }
    }
}

TEST_CASE("Scoring the haplotypes in a subchain", "[haplotypes]") {
    // The last two rows become children of the first two, which reorders the sequences.
    std::vector<std::string> matrix {
        "00000",
        "11111",
        "10110",
        "01001",
        "10111",
        "01000"
    };
    Haplotypes::Subchain subchain = build_subchain(matrix, 5);

    // Arbitrary weights for the kmers. (Presence does not matter.)
    std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>> kmers {
        { KmerPresenceMatrix::kmer_classification::present, 1.0 },
        { KmerPresenceMatrix::kmer_classification::present, -2.0 },
        { KmerPresenceMatrix::kmer_classification::present, 0.5 },
        { KmerPresenceMatrix::kmer_classification::present, 3.0 },
        { KmerPresenceMatrix::kmer_classification::present, -1.5 }
    };

    SECTION("every haplotype is scored once") {
        std::vector<size_t> reported(subchain.sequences.size(), 0);
        std::vector<double> scores(subchain.sequences.size(), 0.0);
        subchain.score_sequences(
            kmers,
            [&](size_t haplotype_id, double score) {
                REQUIRE(haplotype_id < subchain.sequences.size());
                reported[haplotype_id]++;
                scores[haplotype_id] = score;
            }
        );
        for (size_t i = 0; i < scores.size(); i++) {
            REQUIRE(reported[i] == 1);
            const std::string& row = matrix[subchain.sequences[i].first];
            double expected = 0.0;
            for (size_t j = 0; j < kmers.size(); j++) {
                expected += (row[j] == '1' ? 1.0 : -1.0) * kmers[j].second;
            }
            REQUIRE(scores[i] == expected);
        }
    }

    SECTION("no kmers") {
        Haplotypes::Subchain empty = build_subchain({ "", "", "" }, 0);
        std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>> kmers;
        std::vector<double> scores;
        empty.score_sequences(
            kmers,
            [&](size_t haplotype_id, double score) {
                REQUIRE(haplotype_id == scores.size());
                scores.push_back(score);
            }
        );
        REQUIRE(scores.size() == 3);
        for (double score : scores) {
            REQUIRE(score == 0.0);
        }
    }

    SECTION("no haplotypes") {
        Haplotypes::Subchain empty = build_subchain({}, 5);
        size_t calls = 0;
        empty.score_sequences(
            kmers, // We reuse the earlier kmers with actual scores.
            [&](size_t, double) { calls++; }
        );
        REQUIRE(calls == 0);
    }
}

//------------------------------------------------------------------------------

TEST_CASE("Kmer presence rows", "[haplotypes]") {
    // Clustered rows make a compressed matrix with children.
    for (size_t num_kmers : { size_t(0), size_t(1), size_t(63), size_t(64), size_t(130) }) {
        std::vector<std::string> rows(12, "");
        if (num_kmers > 0) {
            rows = clustered_rows(3, 4, num_kmers, 2, 0x1234 + num_kmers);
        }
        Haplotypes::Subchain subchain = build_subchain(rows, num_kmers);
        if (num_kmers >= 63) {
            REQUIRE(subchain.kmers_present.get_num_parents() < subchain.num_sequences());
        }
        for (size_t i = 0; i < subchain.num_sequences(); i++) {
            std::vector<std::uint64_t> row;
            subchain.presence_row(i, row);
            REQUIRE(row.size() == (num_kmers + 63) / 64);
            const std::string& source = rows[subchain.sequences[i].first];
            for (size_t kmer_id = 0; kmer_id < num_kmers; kmer_id++) {
                bool bit = (row[kmer_id / 64] >> (kmer_id % 64)) & 1;
                REQUIRE(bit == (source[kmer_id] == '1'));
            }
        }
    }
}

namespace {

// Adds kmer copies to a subchain. `copies[i]` lists (kmer id, copies) for sequence `i`.
void add_copies(Haplotypes::Subchain& subchain, const std::vector<std::vector<std::pair<size_t, size_t>>>& copies) {
    size_t total = 0;
    for (const auto& list : copies) {
        total += list.size();
    }
    subchain.copy_offsets = sdsl::int_vector<0>(copies.size() + 1, 0, 64);
    subchain.copy_kmers = sdsl::int_vector<0>(total, 0, 64);
    subchain.copy_values = sdsl::int_vector<0>(total, 0, 64);
    size_t offset = 0;
    for (size_t i = 0; i < copies.size(); i++) {
        subchain.copy_offsets[i] = offset;
        for (const auto& [kmer, value] : copies[i]) {
            subchain.copy_kmers[offset] = kmer;
            subchain.copy_values[offset] = value;
            offset++;
        }
    }
    subchain.copy_offsets[copies.size()] = offset;
}

std::vector<std::pair<size_t, size_t>> get_copies(const Haplotypes::Subchain& subchain, size_t sequence) {
    std::vector<std::pair<size_t, size_t>> result;
    subchain.for_each_copy(sequence, [&result](size_t kmer_id, size_t copies) {
        result.emplace_back(kmer_id, copies);
    });
    return result;
}

void check_same_kmers(const Haplotypes::Subchain& subchain, const Haplotypes::Subchain& loaded) {
    REQUIRE(loaded.type == subchain.type);
    REQUIRE(loaded.start == subchain.start);
    REQUIRE(loaded.end == subchain.end);
    REQUIRE(loaded.kmers == subchain.kmers);
    REQUIRE(loaded.sequences == subchain.sequences);
    for (size_t i = 0; i < subchain.num_sequences(); i++) {
        std::vector<std::uint64_t> expected;
        std::vector<std::uint64_t> found;
        subchain.presence_row(i, expected);
        loaded.presence_row(i, found);
        REQUIRE(found == expected);
    }
}

} // anonymous namespace

TEST_CASE("Kmer copies in a subchain", "[haplotypes]") {
    std::vector<std::string> matrix {
        "1101",
        "0111",
        "1100"
    };
    Haplotypes::Subchain subchain = build_subchain(matrix, 4);
    std::vector<std::vector<std::pair<size_t, size_t>>> copies {
        { { 0, 2 }, { 3, 5 } },
        {},
        { { 1, 3 } }
    };

    SECTION("a subchain without copies") {
        REQUIRE(!subchain.has_copies());
        for (size_t i = 0; i < subchain.num_sequences(); i++) {
            REQUIRE(get_copies(subchain, i).empty());
        }
    }

    SECTION("copies are listed by sequence") {
        add_copies(subchain, copies);
        REQUIRE(subchain.has_copies());
        for (size_t i = 0; i < subchain.num_sequences(); i++) {
            REQUIRE(get_copies(subchain, i) == copies[i]);
        }
    }

    SECTION("serialization") {
        add_copies(subchain, copies);
        std::stringstream buffer;
        buffer.exceptions(std::ios::failbit | std::ios::badbit);
        subchain.simple_sds_serialize(buffer);
        REQUIRE(buffer.str().size() == subchain.simple_sds_size() * sizeof(sdsl::simple_sds::element_type));
        Haplotypes::Subchain loaded;
        loaded.simple_sds_load(buffer, Haplotypes::Header::VERSION);
        REQUIRE(buffer.peek() == EOF);
        check_same_kmers(subchain, loaded);
        REQUIRE(loaded.has_copies());
        for (size_t i = 0; i < subchain.num_sequences(); i++) {
            REQUIRE(get_copies(loaded, i) == copies[i]);
        }
    }

    SECTION("an earlier version has no copies") {
        std::uint32_t version = Haplotypes::Header::COPIES_VERSION - 1;
        std::stringstream buffer;
        buffer.exceptions(std::ios::failbit | std::ios::badbit);
        subchain.simple_sds_serialize(buffer, version);
        REQUIRE(buffer.str().size() == subchain.simple_sds_size(version) * sizeof(sdsl::simple_sds::element_type));
        // Loading replaces the copies the subchain had.
        Haplotypes::Subchain loaded = subchain;
        add_copies(loaded, copies);
        loaded.simple_sds_load(buffer, version);
        REQUIRE(buffer.peek() == EOF);
        check_same_kmers(subchain, loaded);
        REQUIRE(!loaded.has_copies());
    }

    SECTION("copies that do not match the sequences are invalid") {
        // Each case breaks one of the conditions: copies for each sequence, as many
        // values as kmers, the last offset, offsets that do not decrease, kmer ids of
        // the subchain, kmer ids in increasing order, and the range of the values.
        std::vector<Haplotypes::Subchain> invalid(8, subchain);
        add_copies(invalid[0], { { { 0, 2 } }, {} });
        add_copies(invalid[1], copies);
        invalid[1].copy_values.resize(invalid[1].copy_kmers.size() - 1);
        add_copies(invalid[2], copies);
        invalid[2].copy_offsets[subchain.num_sequences()] = invalid[2].copy_kmers.size() + 1;
        add_copies(invalid[3], copies);
        invalid[3].copy_offsets[1] = 3;
        add_copies(invalid[4], { { { 4, 2 } }, {}, {} });
        add_copies(invalid[5], { { { 3, 2 }, { 0, 2 } }, {}, {} });
        add_copies(invalid[6], { { { 0, 1 } }, {}, {} });
        add_copies(invalid[7], { { { 0, Haplotypes::Subchain::MAX_COPIES + 1 } }, {}, {} });
        for (const Haplotypes::Subchain& broken : invalid) {
            std::stringstream buffer;
            broken.simple_sds_serialize(buffer);
            Haplotypes::Subchain loaded;
            REQUIRE_THROWS_AS(loaded.simple_sds_load(buffer, Haplotypes::Header::VERSION), sdsl::simple_sds::InvalidData);
        }
    }
}

TEST_CASE("Kmer presence must match the kmers of a subchain", "[haplotypes]") {
    Haplotypes::Subchain subchain = build_subchain({ "1101", "0111", "1100" }, 4);
    subchain.kmers.push_back(subchain.kmers.back() + 1);
    subchain.kmer_counts.resize(subchain.kmers.size());
    std::stringstream buffer;
    subchain.simple_sds_serialize(buffer);
    Haplotypes::Subchain loaded;
    REQUIRE_THROWS_AS(loaded.simple_sds_load(buffer, Haplotypes::Header::VERSION), sdsl::simple_sds::InvalidData);
}

TEST_CASE("Haplotype index file version", "[haplotypes]") {
    Haplotypes haplotypes;
    haplotypes.header.top_level_chains = 1;
    haplotypes.header.construction_jobs = 1;
    haplotypes.chains.resize(1);
    haplotypes.chains[0].offset = 0;
    haplotypes.chains[0].job_id = 0;
    haplotypes.chains[0].placeholder = 0;
    haplotypes.chains[0].subchains.push_back(build_subchain({ "1101", "0111", "1100" }, 4));
    haplotypes.header.total_subchains = 1;
    haplotypes.header.total_kmers = 4;

    // Serializes the haplotypes and returns the version in the file, after loading it.
    auto round_trip = [&haplotypes]() {
        std::stringstream buffer;
        buffer.exceptions(std::ios::failbit | std::ios::badbit);
        haplotypes.simple_sds_serialize(buffer);
        std::uint32_t version = sdsl::simple_sds::load_value<Haplotypes::Header>(buffer).version;
        buffer.seekg(0);
        Haplotypes loaded;
        loaded.simple_sds_load(buffer);
        REQUIRE(buffer.peek() == EOF);
        REQUIRE(loaded.chains.size() == 1);
        REQUIRE(loaded.chains[0].subchains.size() == 1);
        check_same_kmers(haplotypes.chains[0].subchains[0], loaded.chains[0].subchains[0]);
        REQUIRE(loaded.chains[0].subchains[0].has_copies() == haplotypes.chains[0].subchains[0].has_copies());
        return version;
    };

    SECTION("an index without copies uses the earlier version") {
        REQUIRE(haplotypes.file_version() == Haplotypes::Header::COPIES_VERSION - 1);
        REQUIRE(round_trip() == Haplotypes::Header::COPIES_VERSION - 1);
    }

    SECTION("an index with copies uses the current version") {
        add_copies(haplotypes.chains[0].subchains[0], { { { 0, 2 } }, {}, { { 1, 3 } } });
        REQUIRE(haplotypes.file_version() == Haplotypes::Header::VERSION);
        REQUIRE(round_trip() == Haplotypes::Header::VERSION);
    }
}

//------------------------------------------------------------------------------

} // namespace unittest
} // namespace vg
