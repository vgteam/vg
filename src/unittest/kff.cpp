/** \file
 *
 * Unit tests for kff.cpp, which provides tools for working with KFF files.
 */

#include "../kff.hpp"

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

//------------------------------------------------------------------------------

}
}
