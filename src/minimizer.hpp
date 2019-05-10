#ifndef VG_MINIMIZER_HPP_INCLUDED
#define VG_MINIMIZER_HPP_INCLUDED

/** \file 
 * A minimizer index prototype with minimal VG dependencies.
 */

#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "types.hpp"

namespace vg {

//------------------------------------------------------------------------------

/**
 * A class that implements the minimizer index as a hash table mapping kmers to sets of pos_t.
 * The hash table uses quadratic probing with power-of-two size.
 * We encode kmers using 2 bits/character and take wang_hash_64() of the encoding. A minimizer
 * is the kmer with the smallest hash in a window of w consecutive kmers and their reverse
 * complements.
 *
 * Index versions:
 *
 *   1  The initial version.
 *
 *   2  Minimizer selection is based on hashes instead of lexicographic order. A sequence and
 *      its reverse complement have the same minimizers, reducing index size by 50%. Not
 *      compatible with version 1.
 *
 *   3  Construction-time hit cap is no longer used. Compatible with version 2.
 */
class MinimizerIndex {
public:
    typedef std::uint64_t key_type;
    typedef std::uint64_t code_type;
    typedef std::uint32_t offset_type;

    // Public constants.
    // Any key is smaller than NO_KEY and NO_VALUE maps to an empty pos_t.
    constexpr static size_t    KMER_LENGTH      = 21;
    constexpr static size_t    WINDOW_LENGTH    = 11;
    constexpr static size_t    KMER_MAX_LENGTH  = 31;
    constexpr static size_t    INITIAL_CAPACITY = 1024;
    constexpr static double    MAX_LOAD_FACTOR  = 0.77;
    constexpr static key_type  NO_KEY           = std::numeric_limits<key_type>::max();
    constexpr static code_type NO_VALUE         = 0;

    union value_type {
        code_type               value;
        std::vector<code_type>* pointer;
    };

    typedef std::pair<key_type, value_type> cell_type;

    static cell_type empty_cell() { return cell_type(NO_KEY, { NO_VALUE }); }

    struct minimizer_type {
        key_type    key;        // Encoded minimizer.
        size_t      hash;       // Hash of the minimizer.
        offset_type offset;     // First/last offset of the kmer for forward/reverse complement.
        bool        is_reverse; // The minimizer is the reverse complement of the kmer.

        /// Is the minimizer empty?
        bool empty() const { return (this->key == NO_KEY); }

        /// Sort by (offset, !is_reverse). When the offsets are equal, a reverse complement
        /// minimizer is earlier in the sequence than a forward minimizer.
        bool operator<(const minimizer_type& another) const {
            return ((this->offset < another.offset) ||
                    (this->offset == another.offset && this->is_reverse > another.is_reverse));
        }

        bool operator==(const minimizer_type& another) const {
            return (this->key == another.key && this->offset == another.offset && this->is_reverse == another.is_reverse);
        }
    };

    struct Header {
        std::uint32_t tag, version;
        std::uint64_t flags;
        std::uint64_t k, w;
        std::uint64_t keys, capacity, max_keys;
        std::uint64_t values;
        std::uint64_t unused1; // This used to be max_occs.
        std::uint64_t unique;
        std::uint64_t unused2; // This used to be frequent.

        constexpr static std::uint32_t TAG = 0x31513151;
        constexpr static std::uint32_t VERSION = 3;
        constexpr static std::uint32_t MIN_VERSION = 2;

        Header();
        Header(size_t kmer_length, size_t window_length);
        void sanitize();
        bool check() const;

        bool operator==(const Header& another) const;
        bool operator!=(const Header& another) const { return !(this->operator==(another)); }
    };

//------------------------------------------------------------------------------

    /// Constructs an index with the default parameters.
    MinimizerIndex();

    /// Constructs an index with the specified parameter values.
    MinimizerIndex(size_t kmer_length, size_t window_length);

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

    /// Returns the minimizer in the window specified by the iterators. If no minimizer
    /// exists (e.g. because all kmers contain invalid characters), the return value is
    /// an empty minimizer.
    minimizer_type minimizer(std::string::const_iterator begin, std::string::const_iterator end) const;

    /// Returns all minimizers in the string specified by the iterators. The return
    /// value is a vector of minimizers sorted by their offsets.
    std::vector<minimizer_type> minimizers(std::string::const_iterator begin, std::string::const_iterator end) const;

    /// Returns all minimizers in the string. The return value is a vector of
    /// minimizers sorted by their offsets.
    std::vector<minimizer_type> minimizers(const std::string& str) const {
        return this->minimizers(str.begin(), str.end());
    }

    /// Inserts the position into the index, using minimizer.key as the key and
    /// minimizer.hash as its hash. Does not insert empty minimizers or positions.
    /// The offset of the position will be truncated to fit in REV_OFFSET bits.
    /// Use minimizer() or minimizers() to get the minimizer and valid_offset() to check
    /// if the offset fits in the available space.
    /// The position should match the orientation of the minimizer: a path label
    /// starting from the position should have the minimizer as its prefix.
    void insert(const minimizer_type& minimizer, const pos_t& pos);

    /// Returns the sorted set of occurrences of the minimizer.
    /// Use minimizer() or minimizers() to get the minimizer.
    /// If the minimizer is in reverse orientation, use reverse_base_pos() to reverse
    /// the reported occurrences.
    std::vector<pos_t> find(const minimizer_type& minimizer) const;

    /// Returns the occurrence count of the minimizer.
    /// Use minimizer() or minimizers() to get the minimizer.
    size_t count(const minimizer_type& minimizer) const;

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

    /// Is the offset small enough to fit in the low-order bits of the encoding?
    static bool valid_offset(const pos_t& pos) {
        return (offset(pos) <= OFF_MASK);
    }

    /// Encode pos_t as code_type.
    static code_type encode(const pos_t& pos) {
        return (static_cast<code_type>(id(pos)) << ID_OFFSET) |
               (static_cast<code_type>(is_rev(pos)) << REV_OFFSET) |
               (static_cast<code_type>(offset(pos)) & OFF_MASK);
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

    // Find the hash table offset for the key with the given hash value.
    size_t find_offset(key_type key, size_t hash) const;

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
};

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_MINIMIZER_HPP_INCLUDED
