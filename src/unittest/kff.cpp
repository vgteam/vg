/** \file
 *
 * Unit tests for kff.cpp, which provides tools for working with KFF files.
 */

#include "../kff.hpp"
#include "../utility.hpp"

#include "catch.hpp"

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

std::vector<std::string> kmers {
    "GATTA",
    "CATTAC",
    "GATTACA",
    "ATTACCAT",
};

void encode_decode(const std::string& kmer, const uint8_t* encoding) {
    std::vector<uint8_t> encoded = kff_encode(kmer, encoding);
    std::string decoding = kff_invert(encoding);
    std::string decoded = kff_decode(encoded.data(), kmer.length(), decoding);
    REQUIRE(decoded == kmer);
}

void minimizer_recode(const std::string& kmer, const uint8_t* encoding) {
    gbwtgraph::Key64 encoded = gbwtgraph::Key64::encode(kmer);
    std::vector<uint8_t> recoded = kff_recode(encoded.get_key(), kmer.length(), encoding);
    std::string decoding = kff_invert(encoding);
    std::string decoded = kff_decode(recoded.data(), kmer.length(), decoding);
    REQUIRE(decoded == kmer);
}

void rev_comp(const std::string& kmer, const uint8_t* encoding) {
    std::vector<uint8_t> encoded = kff_encode(kmer, encoding);
    std::vector<uint8_t> complemented = kff_reverse_complement(encoded.data(), kmer.length(), encoding);
    std::string decoding = kff_invert(encoding);
    std::string decoded = kff_decode(complemented.data(), kmer.length(), decoding);
    std::string expected = reverse_complement(kmer);
    REQUIRE(decoded == expected);
}

} // Anonymous namespace.

//------------------------------------------------------------------------------

TEST_CASE("Encode and decode", "[kff]") {
    SECTION("default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        for (auto& kmer : kmers) {
            encode_decode(kmer, encoding);
        }
    }

    SECTION("kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        for (auto& kmer : kmers) {
            encode_decode(kmer, encoding);
        }
    }
}

TEST_CASE("Recode minimizers", "[kff]") {
    SECTION("default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        for (auto& kmer : kmers) {
            minimizer_recode(kmer, encoding);
        }
    }

    SECTION("kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        for (auto& kmer : kmers) {
            minimizer_recode(kmer, encoding);
        }
    }
}

TEST_CASE("Parse integers", "[kff]") {
    SECTION("single byte") {
        uint8_t bytes[1] = { 0x12 };
        uint64_t expected = 0x12;
        uint64_t parsed = kff_parse(bytes, sizeof(bytes));
        REQUIRE(parsed == expected);
    }

    SECTION("multiple bytes") {
        uint8_t bytes[5] = { 0x12, 0x34, 0x56, 0x78, 0x9A };
        uint64_t expected = 0x123456789A;
        uint64_t parsed = kff_parse(bytes, sizeof(bytes));
        REQUIRE(parsed == expected);
    }

    SECTION("maximum length") {
        uint8_t bytes[8] = { 0x12, 0x34, 0x56, 0x78, 0x9A, 0xBC, 0xDE, 0xF0 };
        uint64_t expected = 0x123456789ABCDEF0;
        uint64_t parsed = kff_parse(bytes, sizeof(bytes));
        REQUIRE(parsed == expected);
    }
}

TEST_CASE("Reverse complements", "[kff]") {
    SECTION("default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        for (auto& kmer : kmers) {
            rev_comp(kmer, encoding);
        }
    }

    SECTION("kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        for (auto& kmer : kmers) {
            rev_comp(kmer, encoding);
        }
    }
}

//------------------------------------------------------------------------------

}
}
