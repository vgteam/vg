#ifndef VG_MINIMIZER_HPP_INCLUDED
#define VG_MINIMIZER_HPP_INCLUDED

/** \file 
 * A minimizer index prototype with minimal VG dependencies.
 */

#include <cstdint>
#include <iostream>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include "types.hpp"

namespace vg {

//------------------------------------------------------------------------------

/**
 * A class that implements the minimizer index as a hash table mapping kmers to sets of pos_t.
 * The hash table uses quadratic probing with power-of-two size.
 * A minimizer is the lexicographically smallest kmer in a window of w consecutive kmers.
 * There is an option to specify an upper bound for the number of occurrences of each
 * minimizer. If the actual number is higher, the occurrences will not be stored.
 */
class MinimizerIndex {
public:
    typedef std::uint64_t key_type;
    typedef std::uint64_t code_type;

    union value_type {
        code_type               value;
        std::vector<code_type>* pointer;
    };

    typedef std::pair<key_type, value_type> cell_type;
    typedef std::pair<key_type, size_t>     minimizer_type;

    // Public constants.
    // Any key is smaller than NO_KEY and NO_VALUE maps to an empty pos_t.
    constexpr static size_t    KMER_LENGTH      = 21;
    constexpr static size_t    WINDOW_LENGTH    = 11;
    constexpr static size_t    KMER_MAX_LENGTH  = 31;
    constexpr static size_t    INITIAL_CAPACITY = 1024;
    constexpr static double    MAX_LOAD_FACTOR  = 0.77;
    constexpr static size_t    MAX_OCCS         = std::numeric_limits<size_t>::max();
    constexpr static key_type  NO_KEY           = std::numeric_limits<key_type>::max();
    constexpr static code_type NO_VALUE         = 0;

    static cell_type empty_cell() { return cell_type(NO_KEY, { NO_VALUE }); }

    struct Header {
        std::uint32_t tag, version;
        std::uint64_t flags;
        size_t        k, w;
        size_t        keys, capacity, max_keys;
        size_t        values, max_occs;
        size_t        unique, frequent;

        constexpr static std::uint32_t TAG = 0x31513151;
        constexpr static std::uint32_t VERSION = 1;
        constexpr static std::uint32_t MIN_VERSION = 1;

        Header();
        Header(size_t kmer_length, size_t window_length, size_t max_occs_per_key);
        void sanitize();
        bool check() const;

        bool operator==(const Header& another) const;
        bool operator!=(const Header& another) const { return !(this->operator==(another)); }
    };

//------------------------------------------------------------------------------

    /// Constructs an index with the default parameters.
    MinimizerIndex();

    /// Constructs an index with the specified parameter values.
    MinimizerIndex(size_t kmer_length, size_t window_length, size_t max_occs_per_key = MAX_OCCS);

    /// Copy constructor.
    MinimizerIndex(const MinimizerIndex& source);

    /// Move constructor.
    MinimizerIndex(MinimizerIndex&& source);

    /// Destructor.
    ~MinimizerIndex();

    /// Swaps the contents of the indexes.
    void swap(MinimizerIndex& another);

    /// Copy assignment.
    MinimizerIndex& operator=(const MinimizerIndex& source);

    /// Move assignment.
    MinimizerIndex& operator=(MinimizerIndex&& source);

    /// Serialize the index to the ostream. Returns the number of bytes written and
    /// true if the serialization was successful.
    std::pair<size_t, bool> serialize(std::ostream& out) const;

    /// Load the index from the istream and return true if successful.
    bool load(std::istream& in);

    /// Equality comparison for testing.
    bool operator==(const MinimizerIndex& another) const;

    /// Inequality comparison for testing.
    bool operator!=(const MinimizerIndex& another) const { return !(this->operator==(another)); }

//------------------------------------------------------------------------------

    /// Returns the minimizer and its starting offset in the window specified by the
    /// iterators. If no minimizer exists (e.g. because all kmers contain invalid
    /// characters), returns (NO_KEY, 0).
    minimizer_type minimizer(std::string::const_iterator begin, std::string::const_iterator end) const;

    /// Returns all minimizers in the string specified by the iterators. The return
    /// value is a vector of (key, offset) pairs. A minimizer cannot contain invalid
    /// characters.
    std::vector<minimizer_type> minimizers(std::string::const_iterator begin, std::string::const_iterator end) const;

    /// Returns all minimizers in the string. The return value is a vector of
    /// (key, offset) pairs. A minimizer cannot contain invalid characters.
    std::vector<minimizer_type> minimizers(const std::string& str) const {
        return this->minimizers(str.begin(), str.end());
    }

    /// Inserts the minimizer encoded in the key at the given position into the index.
    /// Minimizers with key NO_KEY or a position encoded as NO_VALUE are not inserted.
    /// Use minimizer() or minimizers() to get the key.
    void insert(key_type key, pos_t pos);

    /// Returns the sorted set of occurrences of the kmer encoded in the key.
    /// If the occurrence limit has been exceeded, returns a vector containing an
    /// empty position.
    /// Use minimizer() or minimizers() to get the key.
    std::vector<pos_t> find(key_type key) const;

    /// Returns the occurrence count of the minimizer with the given key.
    /// If the occurrence limit has been exceeded, returns 0.
    /// Use minimizer() or minimizers() to get the key.
    size_t count(key_type key) const;

//------------------------------------------------------------------------------

    /// Length of the kmers in the index.
    size_t k() const { return this->header.k; }

    /// Window length for the minimizers.
    size_t w() const { return this->header.w; }

    /// Number of keys in the index.
    size_t size() const { return this->header.keys; }

    /// Is the index empty.
    bool empty() const { return (this->size() == 0); }

    /// Number of values (minimizer occurrences) in the index.
    size_t values() const { return this->header.values; }

    /// Size of the hash table.
    size_t capacity() const { return this->header.capacity; }

    /// Actual capacity of the hash table. Exceeding it will initiate rehashing.
    size_t max_keys() const { return this->header.max_keys; }

    /// Current load factor of the hash table.
    double load_factor() const { return static_cast<double>(this->size()) / static_cast<double>(this->capacity()); }

    /// Number of minimizers with a single occurrence.
    size_t unique_keys() const { return this->header.unique; }

    /// Number of minimizers with too many occurrences.
    size_t frequent_keys() const { return this->header.frequent; }

//------------------------------------------------------------------------------

private:
    Header                 header;
    std::vector<cell_type> hash_table;
    std::vector<bool>      is_pointer;

//------------------------------------------------------------------------------

public:
    // Constants for the encoding between std::string and key_type.
    constexpr static size_t   PACK_WIDTH = 2;
    constexpr static key_type PACK_MASK  = 0x3;

    // Arrays for the encoding between std::string and key_type.
    const static std::vector<unsigned char> CHAR_TO_PACK;
    const static std::vector<char>          PACK_TO_CHAR;
    const static std::vector<key_type>      KMER_MASK;

    // Constants for the encoding between pos_t and code_type.
    constexpr static size_t    ID_OFFSET  = 11;
    constexpr static size_t    REV_OFFSET = 10;
    constexpr static code_type REV_MASK   = 0x400;
    constexpr static code_type OFF_MASK   = 0x3FF;

    /// Encode pos_t as code_type.
    static code_type encode(pos_t pos) {
        return (static_cast<code_type>(std::get<0>(pos)) << ID_OFFSET) |
               (static_cast<code_type>(std::get<1>(pos)) << REV_OFFSET) |
               (static_cast<code_type>(std::get<2>(pos)) & OFF_MASK);
    }

    /// Copied from position.cpp to avoid excessive dependencies.
    static pos_t make_pos_t(id_t id, bool is_reverse, off_t offset) {
        return std::make_tuple(id, is_reverse, offset);
    }

    /// Decode code_type as pos_t.
    static pos_t decode(code_type pos) {
        return make_pos_t(pos >> ID_OFFSET, pos & REV_MASK, pos & OFF_MASK);
    }

//------------------------------------------------------------------------------

private:
    void copy(const MinimizerIndex& source);
    void clear(size_t i);   // Deletes the pointer at hash_table[i].
    void clear();           // Deletes all pointers in the hash table.

    // Find the hash table offset for the key.
    size_t find_offset(key_type key) const;

    // Insert (key, pos) to hash_table[offset], which is assumed to be empty.
    // Rehashing may be necessary.
    void insert(key_type key, code_type pos, size_t offset);

    // Add pos to the list of occurrences of key at hash_table[offset].
    // The occurrences may be deleted.
    void append(key_type key, code_type pos, size_t offset);

    // Does the list of occurrences at hash_table[offset] contain pos?
    bool contains(size_t offset, code_type pos) const;

    // Double the size of the hash table.
    void rehash();

    // A separate copy of Thomas Wang's hash function for 64-bit integers.
    static size_t hash(key_type key) {
        key = (~key) + (key << 21); // key = (key << 21) - key - 1;
        key = key ^ (key >> 24);
        key = (key + (key << 3)) + (key << 8); // key * 265
        key = key ^ (key >> 14);
        key = (key + (key << 2)) + (key << 4); // key * 21
        key = key ^ (key >> 28);
        key = key + (key << 31);
        return key;
    }
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_MINIMIZER_HPP_INCLUDED
