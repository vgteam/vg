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
        auto packed = gbwtgraph::KmerEncoding::CHAR_TO_PACK[static_cast<uint8_t>(kmer[i])];
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

//------------------------------------------------------------------------------

namespace {

// Return the value of a variable that the latest variable section of the file
// sets. Throw `std::runtime_error` if it does not set it.
uint64_t kff_variable(const Kff_file& file, const std::string& name) {
    auto found = file.global_vars.find(name);
    if (found == file.global_vars.end()) {
        throw std::runtime_error("ParallelKFFReader: file " + file.filename + " does not set variable " + name);
    }
    return found->second;
}

} // anonymous namespace

ParallelKFFReader::ParallelKFFReader(const std::string& filename) :
    file(filename, "r")
{
    // The variables that the first block section uses come before it.
    this->open_next_section();

    this->k = kff_variable(this->file, "k");
    if (this->k > gbwtgraph::Key64::KMER_MAX_LENGTH) {
        throw std::runtime_error("ParallelKFFReader: file " + filename + " contains " + std::to_string(this->k) +
            "-mers; cannot use k > " + std::to_string(gbwtgraph::Key64::KMER_MAX_LENGTH));
    }

    this->max_kmers_per_block = kff_variable(this->file, "max");
    this->data_bytes = kff_variable(this->file, "data_size");
    this->block_buffer.resize(kff_bytes(this->max_kmers_per_block + this->k - 1));
    this->data_buffer.resize(this->max_kmers_per_block * this->data_bytes);

    for (size_t i = 0; i < 4; i++) {
        this->encoding[i] = this->file.encoding[i];
    }
    this->recoding = kff_recoding(this->encoding);
    this->trivial_encoding = kff_is_trivial(this->encoding);
}

void ParallelKFFReader::open_next_section() {
    this->section.reset();
    while (this->section == nullptr && this->file.tellp() < this->file.end_position) {
        char type = this->file.read_section_type();
        if (type == 'v') {
            // Opening the section reads its variables into `file.global_vars`.
            Section_GV variables(&this->file);
        } else if (type == 'i') {
            Section_Index index(&this->file);
            index.close();
        } else {
            this->section.reset(Block_section_reader::construct_section(&this->file));
            if (this->section == nullptr) {
                throw std::runtime_error("ParallelKFFReader: file " + this->file.filename +
                    " has a section of unknown type " + std::string(1, type));
            }
        }
    }
}

Block_section_reader* ParallelKFFReader::next_section() {
    while (this->section != nullptr && this->section->remaining_blocks == 0) {
        this->open_next_section();
        if (this->section != nullptr && (this->section->k != this->k || this->section->data_size != this->data_bytes ||
            this->section->max > this->max_kmers_per_block)) {
            throw std::runtime_error("ParallelKFFReader: file " + this->file.filename +
                " has sections with different kmer lengths, count sizes, or block lengths");
        }
    }
    return this->section.get();
}

bool ParallelKFFReader::has_single_kmer_records(const Block_section_reader& blocks) {
    // A raw section stores the number of kmers in a block only if a block can
    // hold more than one kmer. Otherwise a block is the packed kmer followed by
    // its count.
    return (blocks.max == 1 && dynamic_cast<const Section_Raw*>(&blocks) != nullptr);
}

std::pair<ParallelKFFReader::kmer_type, size_t> ParallelKFFReader::decode_single(const uint8_t* block, const uint8_t* data) const {
    kmer_type kmer = (this->trivial_encoding ? kff_recode_trivial(block, this->k, kff_bytes(this->k)) : kff_recode(block, this->k, this->recoding));
    return { kmer, kff_parse(data, this->data_bytes) };
}

void ParallelKFFReader::read_block(Block_section_reader& from, std::vector<std::pair<kmer_type, size_t>>& result, size_t space) {
    size_t block_size = from.read_compacted_sequence(this->block_buffer.data(), this->data_buffer.data());
    if (block_size == 1) {
        result.push_back(this->decode_single(this->block_buffer.data(), this->data_buffer.data()));
        return;
    }
    std::vector<kmer_type> kmers = kff_recode(this->block_buffer.data(), block_size, this->k, this->recoding);
    for (size_t i = 0; i < kmers.size(); i++) {
        std::pair<kmer_type, size_t> kmer(kmers[i], kff_parse(this->data_buffer.data() + i * this->data_bytes, this->data_bytes));
        if (i < space) {
            result.push_back(kmer);
        } else {
            this->buffer.push_back(kmer);
        }
    }
}

std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> ParallelKFFReader::read(size_t n) {
    std::vector<std::pair<kmer_type, size_t>> result;
    result.reserve(n);

    // Single-kmer blocks are copied from the file under the lock and decoded
    // after it. Allocate the space for them before taking the lock.
    size_t kmer_bytes = kff_bytes(this->k);
    size_t record_bytes = kmer_bytes + this->data_bytes;
    auto records = std::make_unique_for_overwrite<uint8_t[]>(n * record_bytes);
    size_t num_records = 0;

    {
        std::scoped_lock lock(this->mtx);

        while (!this->buffer.empty() && result.size() < n) {
            result.push_back(this->buffer.front());
            this->buffer.pop_front();
        }

        Block_section_reader* current = this->next_section();
        while (current != nullptr && result.size() + num_records < n) {
            size_t space = n - result.size() - num_records;
            if (has_single_kmer_records(*current)) {
                // The blocks are fixed-size records, so a run of them is one read.
                size_t run = std::min(static_cast<size_t>(current->remaining_blocks), space);
                this->file.read(records.get() + num_records * record_bytes, run * record_bytes);
                current->remaining_blocks -= run;
                num_records += run;
            } else {
                this->read_block(*current, result, space);
            }
            current = this->next_section();
        }
    }

    for (size_t i = 0; i < num_records; i++) {
        const uint8_t* record = records.get() + i * record_bytes;
        result.push_back(this->decode_single(record, record + kmer_bytes));
    }

    return result;
}

//------------------------------------------------------------------------------

} // namespace vg