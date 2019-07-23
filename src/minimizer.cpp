#include "minimizer.hpp"
#include "wang_hash.hpp"

#include <algorithm>

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t MinimizerIndex::KMER_LENGTH;
constexpr size_t MinimizerIndex::WINDOW_LENGTH;
constexpr size_t MinimizerIndex::KMER_MAX_LENGTH;
constexpr size_t MinimizerIndex::INITIAL_CAPACITY;
constexpr double MinimizerIndex::MAX_LOAD_FACTOR;
constexpr MinimizerIndex::key_type MinimizerIndex::NO_KEY;
constexpr MinimizerIndex::code_type MinimizerIndex::NO_VALUE;

constexpr std::uint32_t MinimizerIndex::Header::TAG;
constexpr std::uint32_t MinimizerIndex::Header::VERSION;
constexpr std::uint32_t MinimizerIndex::Header::MIN_VERSION;

constexpr size_t MinimizerIndex::PACK_WIDTH;
constexpr MinimizerIndex::key_type MinimizerIndex::PACK_MASK;
constexpr size_t MinimizerIndex::ID_OFFSET;
constexpr size_t MinimizerIndex::REV_OFFSET;
constexpr MinimizerIndex::code_type MinimizerIndex::REV_MASK;
constexpr MinimizerIndex::code_type MinimizerIndex::OFF_MASK;

//------------------------------------------------------------------------------

// Other class variables.

const std::vector<unsigned char> MinimizerIndex::CHAR_TO_PACK = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

const std::vector<char> MinimizerIndex::PACK_TO_CHAR = { 'A', 'C', 'G', 'T' };

const std::vector<MinimizerIndex::key_type> MinimizerIndex::KMER_MASK = {
    0x0000000000000000ull,
    0x0000000000000003ull,
    0x000000000000000Full,
    0x000000000000003Full,
    0x00000000000000FFull,
    0x00000000000003FFull,
    0x0000000000000FFFull,
    0x0000000000003FFFull,
    0x000000000000FFFFull,
    0x000000000003FFFFull,
    0x00000000000FFFFFull,
    0x00000000003FFFFFull,
    0x0000000000FFFFFFull,
    0x0000000003FFFFFFull,
    0x000000000FFFFFFFull,
    0x000000003FFFFFFFull,
    0x00000000FFFFFFFFull,
    0x00000003FFFFFFFFull,
    0x0000000FFFFFFFFFull,
    0x0000003FFFFFFFFFull,
    0x000000FFFFFFFFFFull,
    0x000003FFFFFFFFFFull,
    0x00000FFFFFFFFFFFull,
    0x00003FFFFFFFFFFFull,
    0x0000FFFFFFFFFFFFull,
    0x0003FFFFFFFFFFFFull,
    0x000FFFFFFFFFFFFFull,
    0x003FFFFFFFFFFFFFull,
    0x00FFFFFFFFFFFFFFull,
    0x03FFFFFFFFFFFFFFull,
    0x0FFFFFFFFFFFFFFFull,
    0x3FFFFFFFFFFFFFFFull
};

//------------------------------------------------------------------------------

MinimizerIndex::Header::Header() :
    tag(TAG), version(VERSION),
    flags(0),
    k(KMER_LENGTH), w(WINDOW_LENGTH),
    keys(0), capacity(INITIAL_CAPACITY), max_keys(INITIAL_CAPACITY * MAX_LOAD_FACTOR),
    values(0), unused1(0),
    unique(0), unused2(0)
{
}

MinimizerIndex::Header::Header(size_t kmer_length, size_t window_length) :
    tag(TAG), version(VERSION),
    flags(0),
    k(kmer_length), w(window_length),
    keys(0), capacity(INITIAL_CAPACITY), max_keys(INITIAL_CAPACITY * MAX_LOAD_FACTOR),
    values(0), unused1(0),
    unique(0), unused2(0)
{
    this->sanitize();
}

void MinimizerIndex::Header::sanitize() {
    if (this->k > KMER_MAX_LENGTH) {
        std::cerr << "warning: [MinimizerIndex] Adjusting k from " << this->k << " to " << KMER_MAX_LENGTH << std::endl;
        this->k = KMER_MAX_LENGTH;
    }
    if (this->k == 0) {
        std::cerr << "warning: [MinimizerIndex] Adjusting k from " << this->k << " to " << 1 << std::endl;
        this->k = 1;
    }

    if (this->w == 0) {
        std::cerr << "warning: [MinimizerIndex] Adjusting w from " << this->w << " to " << 1 << std::endl;
        this->w = 1;
    }
}

bool MinimizerIndex::Header::check() const {
    return (this->tag == TAG && this->version >= MIN_VERSION && this->version <= VERSION && this->flags == 0);
}

bool MinimizerIndex::Header::operator==(const Header& another) const {
    return (this->tag == another.tag && this->version == another.version &&
            this->flags == another.flags &&
            this->k == another.k && this->w == another.w &&
            this->keys == another.keys && this->capacity == another.capacity && this->max_keys == another.max_keys &&
            this->values == another.values && this->unused1 == another.unused1 &&
            this->unique == another.unique && this->unused2 == another.unused2);
}

//------------------------------------------------------------------------------

MinimizerIndex::MinimizerIndex() :
    header(),
    hash_table(this->header.capacity, empty_cell()),
    is_pointer(this->header.capacity, false)
{
}

MinimizerIndex::MinimizerIndex(size_t kmer_length, size_t window_length) :
    header(kmer_length, window_length),
    hash_table(this->header.capacity, empty_cell()),
    is_pointer(this->header.capacity, false)
{
}

MinimizerIndex::MinimizerIndex(const MinimizerIndex& source) {
    this->copy(source);
}

MinimizerIndex::MinimizerIndex(MinimizerIndex&& source) {
    *this = std::move(source);
}

MinimizerIndex::~MinimizerIndex() {
    this->clear();
}

void MinimizerIndex::swap(MinimizerIndex& another) {
    if (&another == this) {
        return;
    }

    std::swap(this->header, another.header);
    this->hash_table.swap(another.hash_table);
    this->is_pointer.swap(another.is_pointer);
}

MinimizerIndex& MinimizerIndex::operator=(const MinimizerIndex& source) {
    if (&source != this) {
        this->copy(source);
    }
    return *this;
}

MinimizerIndex& MinimizerIndex::operator=(MinimizerIndex&& source) {
    if (&source != this) {
        this->header = std::move(source.header);
        this->hash_table = std::move(source.hash_table);
        this->is_pointer = std::move(source.is_pointer);
    }
    return *this;
}

namespace mi {

constexpr static size_t BLOCK_SIZE = 4 * 1024 * 1024;
constexpr static size_t WORD_BITS  = 64;

// Serialize a simple element.
template<typename Element>
size_t serialize(std::ostream& out, const Element& element, bool& ok) {
    out.write(reinterpret_cast<const char*>(&element), sizeof(element));
    if (out.fail()) {
        ok = false;
        return 0;
    }
    return sizeof(element);
}

// Load a simple element and return true if successful.
template<typename Element>
bool load(std::istream& in, Element& element) {
    in.read(reinterpret_cast<char*>(&element), sizeof(element));
    return (in.gcount() == sizeof(element));
}

// Serialize the size of a container.
template<class Container>
size_t serialize_size(std::ostream& out, const Container& c, bool &ok) {
    size_t size = c.size();
    return serialize(out, size, ok);
}

// Resize the container to the serialized size.
template<class Container>
bool load_size(std::istream& in, Container& c) {
    size_t size = 0;
    if (!load(in, size)) {
        return false;
    }
    c.resize(size);
    return true;
}

// Serialize a vector of simple elements in blocks.
template<typename Element>
size_t serialize_vector(std::ostream& out, const std::vector<Element>& v, bool& ok) {
    size_t bytes = 0;

    bytes += serialize_size(out, v, ok);

    // Data in blocks of BLOCK_SIZE elements.
    for (size_t i = 0; i < v.size(); i += BLOCK_SIZE) {
        size_t block_size = std::min(v.size() - i, BLOCK_SIZE);
        size_t byte_size = block_size * sizeof(Element);
        out.write(reinterpret_cast<const char*>(v.data() + i), byte_size);
        if (out.fail()) {
            ok = false;
            return bytes;
        }
        bytes += byte_size;
    }

    return bytes;
}

// Load a serialized vector of simple elements.
template<typename Element>
bool load_vector(std::istream& in, std::vector<Element>& v) {
    if (!load_size(in, v)) {
        return false;
    }

    // Data in blocks of BLOCK_SIZE elements.
    for (size_t i = 0; i < v.size(); i += BLOCK_SIZE) {
        size_t block_size = std::min(v.size() - i, BLOCK_SIZE);
        size_t byte_size = block_size * sizeof(Element);
        in.read(reinterpret_cast<char*>(v.data() + i), byte_size);
        if (in.gcount() != byte_size) {
            return false;
        }
    }

    return true;
}

// Serialize a hash table, replacing pointers with empty values.
// The hash table can be loaded with load_vector().
size_t serialize_hash_table(std::ostream& out, const std::vector<MinimizerIndex::cell_type>& hash_table,
                            const std::vector<bool>& is_pointer, bool& ok) {
    size_t bytes = 0;

    bytes += serialize_size(out, hash_table, ok);

    // Data in blocks of BLOCK_SIZE elements. Replace pointers with NO_VALUE to ensure
    // that the file contents are deterministic.
    for (size_t i = 0; i < hash_table.size(); i += BLOCK_SIZE) {
        size_t block_size = std::min(hash_table.size() - i, BLOCK_SIZE);
        size_t byte_size = block_size * sizeof(MinimizerIndex::cell_type);
        std::vector<MinimizerIndex::cell_type> buffer(hash_table.begin() + i, hash_table.begin() + i + block_size);
        for (size_t j = 0; j < buffer.size(); j++) {
            if (is_pointer[i + j]) {
                buffer[j].second.value = MinimizerIndex::NO_VALUE;
            }
        }
        out.write(reinterpret_cast<const char*>(buffer.data()), byte_size);
        if (out.fail()) {
            ok = false;
            return bytes;
        }
        bytes += byte_size;
    }

    return bytes;
}

// Serialize a boolean vector in blocks.
size_t serialize_bool_vector(std::ostream& out, const std::vector<bool>& v, bool& ok) {
    size_t bytes = 0;

    bytes += serialize_size(out, v, ok);

    // Data in blocks of BLOCK_SIZE words.
    for (size_t i = 0; i < v.size(); i += BLOCK_SIZE * WORD_BITS) {
        size_t block_size = std::min(v.size() - i, BLOCK_SIZE * WORD_BITS);
        size_t word_size = (block_size + WORD_BITS - 1) / WORD_BITS;
        size_t byte_size = word_size * sizeof(std::uint64_t);
        std::vector<std::uint64_t> buffer(word_size, 0);
        for (size_t j = 0; j < block_size; j++) {
            if (v[i + j]) {
                buffer[j / WORD_BITS] |= static_cast<std::uint64_t>(1) << (j % WORD_BITS);
            }
        }
        out.write(reinterpret_cast<const char*>(buffer.data()), byte_size);
        if (out.fail()) {
            ok = false;
            return bytes;
        }
        bytes += byte_size;
    }

    return bytes;
}

// Load a serialized boolean vector.
bool load_bool_vector(std::istream& in, std::vector<bool>& v) {
    if (!load_size(in, v)) {
        return false;
    }

    // Data in blocks of BLOCK_SIZE words.
    for (size_t i = 0; i < v.size(); i += BLOCK_SIZE * WORD_BITS) {
        size_t block_size = std::min(v.size() - i, BLOCK_SIZE * WORD_BITS);
        size_t word_size = (block_size + WORD_BITS - 1) / WORD_BITS;
        size_t byte_size = word_size * sizeof(std::uint64_t);
        std::vector<std::uint64_t> buffer(word_size, 0);
        in.read(reinterpret_cast<char*>(buffer.data()), byte_size);
        if (in.gcount() != byte_size) {
            return false;
        }
        for (size_t j = 0; j < block_size; j++) {
            v[i + j] = static_cast<bool>(buffer[j / WORD_BITS] & (static_cast<std::uint64_t>(1) << (j % WORD_BITS)));
        }
    }

    return true;
}

} // namespace mi

std::pair<size_t, bool> MinimizerIndex::serialize(std::ostream& out) const {
    size_t bytes = 0;
    bool ok = true;

    bytes += mi::serialize(out, this->header, ok);
    bytes += mi::serialize_hash_table(out, this->hash_table, this->is_pointer, ok);
    bytes += mi::serialize_bool_vector(out, this->is_pointer, ok);

    // Serialize the occurrence lists.
    for (size_t i = 0; i < this->capacity(); i++) {
        if (this->is_pointer[i]) {
            bytes += mi::serialize_vector(out, *(this->hash_table[i].second.pointer), ok);
        }
    }

    if (!ok) {
        std::cerr << "error: [MinimizerIndex] serialization failed" << std::endl;
    }

    return std::make_pair(bytes, ok);
}

bool MinimizerIndex::load(std::istream& in) {
    bool ok = true;

    // Load and check the header.
    ok &= mi::load(in, this->header);
    if (!(this->header.check())) {
        std::cerr << "error: [MinimizerIndex] invalid or old index file" << std::endl;
        std::cerr << "error: [MinimizerIndex] index version is " << this->header.version << "; expected " << Header::VERSION << std::endl;
        return false;
    }

    // Load the hash table.
    if (ok) {
        ok &= mi::load_vector(in, this->hash_table);
    }
    if (ok) {
        ok &= mi::load_bool_vector(in, this->is_pointer);
    }

    // Load the occurrence lists.
    if (ok) {
        for (size_t i = 0; i < this->capacity(); i++) {
            if (this->is_pointer[i]) {
                this->hash_table[i].second.pointer = new std::vector<code_type>();
                ok &= mi::load_vector(in, *(this->hash_table[i].second.pointer));
            }
        }
    }

    if (!ok) {
        std::cerr << "error: [MinimizerIndex] index loading failed" << std::endl;
    }

    return ok;
}

bool MinimizerIndex::operator==(const MinimizerIndex& another) const {
    if (this->header != another.header || this->is_pointer != another.is_pointer) {
        return false;
    }

    for (size_t i = 0; i < this->capacity(); i++) {
        cell_type a = this->hash_table[i], b = another.hash_table[i];
        if (a.first != b.first) {
            return false;
        }
        if (this->is_pointer[i]) {
            if (*(a.second.pointer) != *(b.second.pointer)) {
                return false;
            }
        } else {
            if (a.second.value != b.second.value) {
                return false;
            }
        }
    }

    return true;
}

void MinimizerIndex::copy(const MinimizerIndex& source) {
    this->clear();
    this->header = source.header;
    this->hash_table = source.hash_table;
    this->is_pointer = source.is_pointer;
}

void MinimizerIndex::clear(size_t i) {
    if (this->is_pointer[i]) {
        delete this->hash_table[i].second.pointer;
        this->hash_table[i].second.value = NO_VALUE;
        this->is_pointer[i] = false;
    }
}

void MinimizerIndex::clear() {
    for (size_t i = 0; i < this->hash_table.size(); i++) {
        this->clear(i);
    }    
}

//------------------------------------------------------------------------------

namespace mi {

/*
  A circular buffer of size 2^i for all minimizer candidates. The candidates are sorted
  by both key and position. The candidate at position i is removed when we reach position
  i + w. Candidates at the tail are purged when we advance with a smaller key.
*/
struct CircularBuffer {
    typedef MinimizerIndex::key_type key_type;
    typedef MinimizerIndex::offset_type offset_type;
    typedef MinimizerIndex::minimizer_type minimizer_type;

    std::vector<minimizer_type> buffer;
    size_t head, tail;
    size_t w;

    constexpr static size_t BUFFER_SIZE = 16;

    CircularBuffer(size_t capacity) :
        buffer(),
        head(0), tail(0), w(capacity)
    {
        size_t buffer_size = BUFFER_SIZE;
        while (buffer_size < this->w) {
            buffer_size *= 2;
        }
        this->buffer.resize(buffer_size);
    }

    bool empty() const {
        return (this->head >= this->tail);
    }

    minimizer_type& front() {
        return this->buffer[this->head & (this->buffer.size() - 1)];
    }

    minimizer_type& back() {
        return this->buffer[(this->tail - 1) & (this->buffer.size() - 1)];
    }

    // Advance to the next starting position with a valid k-mer.
    void advance(offset_type pos, key_type forward_key, key_type reverse_key) {
        if (!(this->empty()) && this->front().offset + this->w <= pos) {
            this->head++;
        }
        size_t forward_hash = wang_hash_64(forward_key), reverse_hash = wang_hash_64(reverse_key);
        size_t hash = std::min(forward_hash, reverse_hash);
        while (!(this->empty()) && this->back().hash > hash) {
            this->tail--;
        }
        this->tail++;
        if (reverse_hash < forward_hash) {
            this->back() = { reverse_key, reverse_hash, pos, true };
        } else {
            this->back() = { forward_key, forward_hash, pos, false };
        }
    }

    // Advance to the next starting position without a valid k-mer.
    void advance(std::uint32_t pos) {
        if (!(this->empty()) && this->front().offset + this->w <= pos) {
            this->head++;
        }
    }
};

void
update_forward_key(MinimizerIndex::key_type& key, size_t k, unsigned char c, size_t& valid_chars) {
    MinimizerIndex::key_type packed = MinimizerIndex::CHAR_TO_PACK[c];
    if (packed > MinimizerIndex::PACK_MASK) {
        key = 0;
        valid_chars = 0;
    } else {
        key = ((key << MinimizerIndex::PACK_WIDTH) | packed) & MinimizerIndex::KMER_MASK[k];
        valid_chars++;
    }
}

void
update_reverse_key(MinimizerIndex::key_type& key, size_t k, unsigned char c) {
    MinimizerIndex::key_type packed = MinimizerIndex::CHAR_TO_PACK[c];
    if (packed > MinimizerIndex::PACK_MASK) {
        key = 0;
    } else {
        packed ^= MinimizerIndex::PACK_MASK; // The complement of the base.
        key = (packed << ((k - 1) * MinimizerIndex::PACK_WIDTH)) | (key >> MinimizerIndex::PACK_WIDTH);
    }
}

} // namespace mi


MinimizerIndex::minimizer_type
MinimizerIndex::minimizer(std::string::const_iterator begin, std::string::const_iterator end) const {
    minimizer_type result { NO_KEY, static_cast<size_t>(0), static_cast<std::uint32_t>(0), false };
    if (end - begin < this->k()) {
        return result;
    }

    size_t valid_chars = 0, best_hash = 0;
    bool found = false;
    key_type forward_key = 0, reverse_key = 0;
    for (std::string::const_iterator iter = begin; iter != end; ++iter) {
        mi::update_forward_key(forward_key, this->k(), *iter, valid_chars);
        mi::update_reverse_key(reverse_key, this->k(), *iter);
        if (valid_chars >= this->k()) {
            size_t forward_hash = wang_hash_64(forward_key), reverse_hash = wang_hash_64(reverse_key);
            size_t hash = std::min(forward_hash, reverse_hash);
            if (!found || hash < result.hash) {
                offset_type pos = (iter - begin) + 1 - this->k();
                if (reverse_hash < forward_hash) {
                    result = { reverse_key, reverse_hash, pos, true };
                } else {
                    result = { forward_key, forward_hash, pos, false };
                }
                found = true;
            }
        }
    }

    // It was more convenient to use the first offset of the kmer, regardless of the orientation.
    // If the minimizer is a reverse complement, we must return the last offset instead.
    if (result.is_reverse) {
        result.offset += this->k() - 1;
    }

    return result;
}

std::vector<MinimizerIndex::minimizer_type>
MinimizerIndex::minimizers(std::string::const_iterator begin, std::string::const_iterator end) const {

    std::vector<minimizer_type> result;
    size_t window_length = this->k() + this->w() - 1, total_length = end - begin;
    if (total_length < window_length) {
        return result;
    }

    // Find the minimizers.
    mi::CircularBuffer buffer(this->w());
    size_t valid_chars = 0, start_pos = 0;
    key_type forward_key = 0, reverse_key = 0;
    std::string::const_iterator iter = begin;
    while (iter != end) {
        mi::update_forward_key(forward_key, this->k(), *iter, valid_chars);
        mi::update_reverse_key(reverse_key, this->k(), *iter);
        if (valid_chars >= this->k()) {
            buffer.advance(start_pos, forward_key, reverse_key);
        } else {
            buffer.advance(start_pos);
        }
        ++iter;
        if (iter - begin >= this->k()) {
            start_pos++;
        }
        // We have a full window with a minimizer.
        if (iter - begin >= window_length && !buffer.empty()) {
            if (result.empty() || result.back().offset != buffer.front().offset) {
                result.emplace_back(buffer.front());
            }
        }
    }

    // It was more convenient to use the first offset of the kmer, regardless of the orientation.
    // If the minimizer is a reverse complement, we must return the last offset instead.
    for (minimizer_type& minimizer : result) {
        if (minimizer.is_reverse) {
            minimizer.offset += this->k() - 1;
        }
    }
    std::sort(result.begin(), result.end());

    return result;
}

void MinimizerIndex::insert(const minimizer_type& minimizer, const pos_t& pos) {
    if (minimizer.empty() || is_empty(pos)) {
        return;
    }

    size_t offset = this->find_offset(minimizer.key, minimizer.hash);
    code_type code = encode(pos);
    if (this->hash_table[offset].first == NO_KEY) {
        this->insert(minimizer.key, code, offset);
    } else if (this->hash_table[offset].first == minimizer.key) {
        this->append(minimizer.key, code, offset);
    }
}

std::vector<pos_t> MinimizerIndex::find(const minimizer_type& minimizer) const {
    std::vector<pos_t> result;
    if (minimizer.empty()) {
        return result;
    }

    size_t offset = this->find_offset(minimizer.key, minimizer.hash);
    cell_type cell = this->hash_table[offset];
    if (cell.first == minimizer.key) {
        if (this->is_pointer[offset]) {
            result.reserve(cell.second.pointer->size());
            for (code_type pos : *(cell.second.pointer)) {
                result.emplace_back(decode(pos));
            }
        } else {
            result.emplace_back(decode(cell.second.value));
        }
    }

    return result;
}

size_t MinimizerIndex::count(const minimizer_type& minimizer) const {
    if (minimizer.empty()) {
        return 0;
    }

    size_t offset = this->find_offset(minimizer.key, minimizer.hash);
    if (this->hash_table[offset].first == minimizer.key) {
        return (this->is_pointer[offset] ? this->hash_table[offset].second.pointer->size() : 1);
    }

    return 0;
}

size_t MinimizerIndex::find_offset(key_type key, size_t hash) const {
    size_t offset = hash & (this->capacity() - 1);
    for (size_t attempt = 0; attempt < this->capacity(); attempt++) {
        if (this->hash_table[offset].first == NO_KEY || this->hash_table[offset].first == key) {
            return offset;
        }

        // Quadratic probing with triangular numbers.
        offset = (offset + attempt + 1) & (this->capacity() - 1);
    }

    // This should not happen.
    std::cerr << "error: [MinimizerIndex] Cannot find the offset for key " << key << std::endl;
    return 0;
}

void MinimizerIndex::insert(key_type key, code_type pos, size_t offset) {
    this->hash_table[offset].first = key;
    this->hash_table[offset].second.value = pos;
    this->header.keys++;
    this->header.values++;
    this->header.unique++;

    if (this->size() > this->max_keys()) {
        this->rehash();
    }
}

void MinimizerIndex::append(key_type key, code_type pos, size_t offset) {
    if (this->contains(offset, pos)) {
        return;
    }

    if (this->is_pointer[offset]) {
        std::vector<code_type>* occs = this->hash_table[offset].second.pointer;
        occs->push_back(pos);
        size_t offset = occs->size() - 1;
        while(offset > 0 && occs->at(offset - 1) > occs->at(offset)) {
            std::swap(occs->at(offset - 1), occs->at(offset));
            offset--;
        }
    } else {
        std::vector<code_type>* occs = new std::vector<code_type>(2);
        occs->at(0) = this->hash_table[offset].second.value;
        occs->at(1) = pos;
        if (occs->at(0) > occs->at(1)) {
            std::swap(occs->at(0), occs->at(1));
        }
        this->hash_table[offset].second.pointer = occs;
        this->is_pointer[offset] = true;
        this->header.unique--;
    }
    this->header.values++;
}

bool MinimizerIndex::contains(size_t offset, code_type pos) const {
    if (this->is_pointer[offset]) {
        const std::vector<code_type>* occs = this->hash_table[offset].second.pointer;
        return std::binary_search(occs->begin(), occs->end(), pos);
    } else {
        return (this->hash_table[offset].second.value == pos);
    }
}

void MinimizerIndex::rehash() {
    // Reinitialize with a larger hash table.
    std::vector<cell_type> old_hash_table(2 * this->capacity(), empty_cell());
    std::vector<bool> old_is_pointer(2 * this->capacity(), false);
    this->hash_table.swap(old_hash_table);
    this->is_pointer.swap(old_is_pointer);
    this->header.capacity = this->hash_table.size();
    this->header.max_keys = this->capacity() * MAX_LOAD_FACTOR;

    // Move the keys to the new hash table.
    for (size_t i = 0; i < old_hash_table.size(); i++) {
        key_type key = old_hash_table[i].first;
        if (key == NO_KEY) {
            continue;
        }

        size_t offset = this->find_offset(key, wang_hash_64(key));
        this->hash_table[offset] = old_hash_table[i];
        this->is_pointer[offset] = old_is_pointer[i];
    }
}

//------------------------------------------------------------------------------

} // namespace vg
