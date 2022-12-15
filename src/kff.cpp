#include "kff.hpp"

#include <algorithm>

namespace vg {

//------------------------------------------------------------------------------

// Encode up to 4 characters in one byte.
uint8_t kff_encode(const std::string& kmer, size_t start, size_t limit, const uint8_t* encoding) {
    uint8_t val = 0;
    for (size_t i = start; i < limit; i++) {
        val <<= 2;
        auto packed = gbwtgraph::CHAR_TO_PACK[static_cast<uint8_t>(kmer[i])];
        if (packed < 4) {
            val |= encoding[packed];
        }
    }
    return val;
}

std::vector<uint8_t> kff_encode(const std::string& kmer, const uint8_t* encoding) {
    std::vector<uint8_t> result;
    result.reserve((kmer.length() + 3) / 4);

    // If k is not a multiple of 4, KFF adds the padding to the high-order bits
    // of the first byte.
    size_t remainder = kmer.length() & 3;
    if (remainder > 0) {
        result.push_back(kff_encode(kmer, 0, remainder, encoding));
    }
    for (size_t i = remainder; i < kmer.length(); i += 4) {
        result.push_back(kff_encode(kmer, i, i + 4, encoding));
    }

    return result;
}

//------------------------------------------------------------------------------

std::string kff_invert(const uint8_t* encoding) {
    std::string result(4, ' ');
    result[encoding[0]] = 'A';
    result[encoding[1]] = 'C';
    result[encoding[2]] = 'G';
    result[encoding[3]] = 'T';
    return result;
}

// Decode up to 4 characters from one byte
void kff_decode(uint8_t byte, size_t chars, const std::string& decoding, std::string& output) {
    size_t offset = 2 * chars;
    for (size_t i = 0; i < chars; i++) {
        offset -= 2;
        output.push_back(decoding[(byte >> offset) & 3]);
    }
}

std::string kff_decode(const uint8_t* kmer, size_t k, const std::string& decoding) {
    std::string result;
    result.reserve(k);

    size_t byte_len = (k + 3) / 4;
    size_t chars = k & 3;
    if (chars == 0) {
        chars = 4;
    }
    for (size_t i = 0; i < byte_len; i++) {
        kff_decode(kmer[i], chars, decoding, result);
        chars = 4;
    }

    return result;
}

std::string kff_decode(const uint8_t* kmer, Kff_reader& reader) {
    return kff_decode(kmer, reader.k, kff_invert(reader.get_encoding()));
}

//------------------------------------------------------------------------------

// Recode up to 4 characters in one byte.
uint8_t kff_recode(gbwtgraph::Key64::value_type kmer, size_t k, size_t chars, const uint8_t* encoding) {
    size_t offset = 2 * k;
    uint8_t val = 0;
    for (size_t i = 0; i < chars; i++) {
        offset -= 2;
        val = (val << 2) | encoding[(kmer >> offset) & 3];
    }
    return val;
}

std::vector<uint8_t> kff_recode(gbwtgraph::Key64::value_type kmer, size_t k, const uint8_t* encoding) {
    std::vector<uint8_t> result;
    result.reserve((k + 3) / 4);

    size_t remainder = k & 3;
    if (remainder > 0) {
        result.push_back(kff_recode(kmer, k, remainder, encoding));
    }
    for (size_t i = remainder; i < k; i += 4) {
        result.push_back(kff_recode(kmer, k - i, 4, encoding));
    }

    return result;
}

std::vector<uint8_t> kff_recode(gbwtgraph::Key64::value_type kmer, Kff_reader& reader) {
    return kff_recode(kmer, reader.k, reader.get_encoding());
}

//------------------------------------------------------------------------------

} // namespace vg