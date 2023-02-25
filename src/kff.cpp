#include "kff.hpp"

namespace vg {

//------------------------------------------------------------------------------

bool kff_is_trivial(const uint8_t* encoding) {
    for (size_t i = 0; i < 4; i++) {
        if (encoding[i] != i) {
            return false;
        }
    }
    return true;
}

std::string kff_invert(const uint8_t* encoding) {
    std::string result(4, ' ');
    result[encoding[0]] = 'A';
    result[encoding[1]] = 'C';
    result[encoding[2]] = 'G';
    result[encoding[3]] = 'T';
    return result;
}

kff_recoding_t kff_recoding(const uint8_t* encoding) {
    kff_recoding_t result;
    for (size_t i = 0; i < 4; i++) {
        result.data[encoding[i]] = i;
    }
    return result;
}

uint64_t kff_parse(const uint8_t* data, size_t bytes) {
    uint64_t value = 0;
    size_t shift = 8 * bytes;
    for (size_t i = 0; i < bytes; i++) {
        shift -= 8;
        value |= static_cast<uint64_t>(data[i]) << shift;
    }
    return value;
}

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
    result.reserve(kff_bytes(kmer.length()));

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

    size_t bytes = kff_bytes(k);
    size_t chars = k & 3;
    if (chars == 0) {
        chars = 4;
    }
    for (size_t i = 0; i < bytes; i++) {
        kff_decode(kmer[i], chars, decoding, result);
        chars = 4;
    }

    return result;
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
    result.reserve(kff_bytes(k));

    size_t remainder = k & 3;
    if (remainder > 0) {
        result.push_back(kff_recode(kmer, k, remainder, encoding));
    }
    for (size_t i = remainder; i < k; i += 4) {
        result.push_back(kff_recode(kmer, k - i, 4, encoding));
    }

    return result;
}

gbwtgraph::Key64::value_type kff_recode(const uint8_t* kmer, size_t k, kff_recoding_t recoding) {
    gbwtgraph::Key64::value_type result = 0;

    size_t bytes = kff_bytes(k);
    size_t chars = k & 3;
    if (chars == 0) {
        chars = 4;
    }
    for (size_t i = 0; i < bytes; i++) {
        size_t offset = 2 * chars;
        for (size_t j = 0; j < chars; j++) {
            offset -= 2;
            result = (result << 2) | recoding.data[(kmer[i] >> offset) & 3];
        }
        chars = 4;
    }

    return result;
}

gbwtgraph::Key64::value_type kff_recode_trivial(const uint8_t* kmer, size_t k, size_t bytes) {
    gbwtgraph::Key64::value_type result = 0;
    for (size_t i = 0; i < bytes; i++) {
        result = (result << 8) | kmer[i];
    }
    return result & sdsl::bits::lo_set[2 * k];
}

std::vector<gbwtgraph::Key64::value_type> kff_recode(const uint8_t* kmers, size_t n, size_t k, kff_recoding_t recoding) {
    std::vector<gbwtgraph::Key64::value_type> result;
    result.reserve(n);

    size_t total_chars = n + k - 1;
    size_t bytes = kff_bytes(total_chars);
    size_t chars = total_chars & 3;
    if (chars == 0) {
        chars = 4;
    }

    gbwtgraph::Key64::value_type curr = 0;
    for (size_t i = 0, processed = 0; i < bytes; i++) {
        size_t offset = 2 * chars;
        for (size_t j = 0; j < chars; j++) {
            offset -= 2;
            curr = (curr << 2) | recoding.data[(kmers[i] >> offset) & 3];
            processed++;
            if (processed >= k) {
                result.push_back(curr & sdsl::bits::lo_set[2 * k]);
            }
        }
        chars = 4;
    }

    return result;
}

//------------------------------------------------------------------------------

uint8_t kff_get(const uint8_t* kmer, size_t i) {
    size_t byte = i / 4;
    size_t offset = 3 - (i & 3);
    return (kmer[byte] >> (2 * offset)) & 3;
}

void kff_set(std::vector<uint8_t>& kmer, size_t i, uint8_t value) {
    size_t byte = i / 4;
    size_t offset = 3 - (i & 3);
    kmer[byte] |= value << (2 * offset);
}

std::vector<uint8_t> kff_reverse_complement(const uint8_t* kmer, size_t k, const uint8_t* encoding) {
    uint8_t complement[4];
    complement[encoding[0]] = encoding[3];
    complement[encoding[1]] = encoding[2];
    complement[encoding[2]] = encoding[1];
    complement[encoding[3]] = encoding[0];

    size_t offset = (4 - (k & 3)) & 3;
    std::vector<uint8_t> result(kff_bytes(k), 0);
    for (size_t i = 0; i < k; i++) {
        kff_set(result, 4 * result.size() - 1 - i, complement[kff_get(kmer, i + offset)]);
    }
    return result;
}

gbwtgraph::Key64::value_type minimizer_reverse_complement(gbwtgraph::Key64::value_type kmer, size_t k) {
    gbwtgraph::Key64::value_type result = 0;
    for (size_t i = 0; i < k; i++) {
        result = (result << 2) | ((kmer & 3) ^ 3);
        kmer >>= 2;
    }
    return result;
}

//------------------------------------------------------------------------------

ParallelKFFReader::ParallelKFFReader(const std::string& filename) :
    reader(filename)
{
    this->k = this->reader.get_var("k");
    if (this->k > gbwtgraph::Key64::KMER_MAX_LENGTH) {
        throw std::runtime_error("ParallelKFFReader: file " + filename + " contains " + std::to_string(this->k) +
            "-mers; cannot use k > " + std::to_string(gbwtgraph::Key64::KMER_MAX_LENGTH));
    }

    this->max_kmers_per_block = this->reader.get_var("max");
    this->data_bytes = this->reader.get_var("data_size");

    std::uint8_t* buf = this->reader.get_encoding();
    for (size_t i = 0; i < 4; i++) {
        this->encoding[i] = buf[i];
    }
    this->recoding = kff_recoding(this->encoding);
}

std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> ParallelKFFReader::read(size_t n) {
    std::vector<std::pair<kmer_type, size_t>> result;
    result.reserve(n);

    std::lock_guard<std::mutex> lock(this->mtx);

    while (!this->buffer.empty() && result.size() < n) {
        result.push_back(this->buffer.front());
        this->buffer.pop_front();
    }
    if (result.size() >= n) {
        return result;
    }

    // Because we read kmers by blocks, we have to preallocate the buffers.
    uint8_t* block = new uint8_t[kff_bytes(this->max_kmers_per_block + this->k - 1)];
    uint8_t* data = new uint8_t[this->max_kmers_per_block * this->data_bytes];
    size_t kmer_bytes = kff_bytes(this->k);
    bool trivial_encoding(kff_is_trivial(this->encoding));
    while (this->reader.has_next() && result.size() < n) {
        size_t block_size = this->reader.next_block(block, data);
        if (block_size > 1) {
            std::vector<kmer_type> kmers = kff_recode(block, block_size, this->k, this->recoding);
            for (size_t i = 0; i < block_size; i++) {
                std::pair<kmer_type, size_t> kmer(kmers[i], kff_parse(data + i * data_bytes, data_bytes));
                if (result.size() < n) {
                    result.push_back(kmer);
                } else {
                    buffer.push_back(kmer);
                }
            }
        } else {
            kmer_type kmer = (trivial_encoding ? kff_recode_trivial(block, this->k, kmer_bytes) : kff_recode(block, this->k, this->recoding));
            result.push_back({ kmer, kff_parse(data, data_bytes) });
        }
    }
    delete[] block; block = nullptr;
    delete[] data; data = nullptr;

    return result;
}

//------------------------------------------------------------------------------

} // namespace vg