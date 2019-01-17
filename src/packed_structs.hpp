//
//  packed_structs.hpp
//  
// Contains implementations of classic data structures converted into
// bit-packed integer vectors
//

#ifndef VG_PACKED_STRUCTS_HPP_INCLUDED
#define VG_PACKED_STRUCTS_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <algorithm>
#include <vector>
#include "sdsl/int_vector.hpp"

namespace vg {
    
using namespace std;
    
/*
 * A dynamic integer vector that maintains integers in bit-compressed form.
 * Automatically adjusts bit-width for entries depending on input data.
 */
class PackedVector {
public:
    /// Constructor (starts empty)
    PackedVector();
    
    /// Move constructors
    PackedVector(PackedVector&& other) = default;
    PackedVector& operator=(PackedVector&& other) = default;
        
    /// Destructor
    ~PackedVector();
        
    /// Set the i-th value
    inline void set(const size_t& i, const uint64_t& value);
        
    /// Returns the i-th value
    inline uint64_t get(const size_t& i) const;
        
    /// Add a value to the end
    inline void append(const uint64_t& value);
        
    /// Remove the last value
    inline void pop();
    
    /// Either shrink the vector or grow the vector to the new size. New
    /// entries created by growing are filled with 0.
    inline void resize(size_t new_size);
        
    /// Returns the number of values
    inline size_t size() const;
        
    /// Returns true if there are no entries and false otherwise
    inline bool empty() const;

    /// Clears the backing vector
    inline void clear();
        
private:
        
    // the underlying vector representation
    sdsl::int_vector<> vec;
    // tracker for number of values
    size_t filled = 0;
    // geometric expansion factor
    static const double factor;
};

/*
 * A dynamic integer vector that provides better compression when values in the
 * integer vector either 1) do not vary much from their neighbors or 2) are 0.
 * Compression is also optimized for vectors that are mostly increasing.
 */
class PagedVector {
public:
    
    /// Construct and set page size (starts empty)
    PagedVector(size_t page_size);
    
    /// Move constructors
    PagedVector(PagedVector&& other) = default;
    PagedVector& operator=(PagedVector&& other) = default;
    
    // Destructor
    ~PagedVector();
    
    /// Set the i-th value
    inline void set(const size_t& i, const uint64_t& value);
    
    /// Returns the i-th value
    inline uint64_t get(const size_t& i) const;
    
    /// Add a value to the end
    inline void append(const uint64_t& value);
    
    /// Remove the last value
    inline void pop();
    
    /// Either shrink the vector or grow the vector to the new size. New
    /// entries created by growing are filled with 0.
    inline void resize(size_t new_size);
    
    /// Returns the number of values
    inline size_t size() const;
    
    /// Returns true if there are no entries and false otherwise
    inline bool empty() const;
    
    /// Clears the backing vector
    inline void clear();
    
private:
    
    PagedVector();
    
    inline uint64_t to_diff(const uint64_t& value, const uint64_t& page) const;
    inline uint64_t from_diff(const uint64_t& diff, const uint64_t& page) const;
    
    // TODO: is there a way to const this and still allow copy/move constructors?
    size_t page_size = 64;
    
    // The number of entries filled so far
    size_t filled = 0;
    
    // Evenly spaced entries from the vector
    PackedVector anchors;
    // All entries in the vector expressed as a difference from the preceding page value
    vector<PackedVector> pages;
};

/*
 * A deque implementation that maintains integers in bit-compressed form, with the bit
 * width automatically adjusted to the entries.
 */
class PackedDeque {
public:
    PackedDeque(void);
    ~PackedDeque(void);
    
    /// Move constructors
    PackedDeque(PackedDeque&& other) = default;
    PackedDeque& operator=(PackedDeque&& other) = default;
    
    /// Set the i-th value
    inline void set(const size_t& i, const uint64_t& value);
    
    /// Returns the i-th value
    inline uint64_t get(const size_t& i) const;
    
    /// Add a value to the front
    inline void append_front(const uint64_t& value);
    
    /// Add a value to the back
    inline void append_back(const uint64_t& value);
    
    /// Remove the front value
    inline void pop_front();
    
    /// Remove the back value
    inline void pop_back();
    
    /// Returns the number of values
    inline size_t size() const;
    
    /// Returns true if there are no entries and false otherwise
    inline bool empty() const;
    
    /// Empty the contents
    inline void clear();
    
private:
    
    inline void contract();
    
    inline size_t internal_index(const size_t& i) const;
    
    PackedVector vec;
    
    size_t begin_idx = 0;
    size_t filled = 0;
    static const double factor;
};
    
/*
 * A splay-tree implementation that stores keys, values, and pointers in
 * bit-compressed form.
 */
class PackedSplayTree {

public:
    PackedSplayTree(void);
    ~PackedSplayTree(void);
    
    /// Move constructors
    PackedSplayTree(PackedSplayTree&& other) = default;
    PackedSplayTree& operator=(PackedSplayTree&& other) = default;
    
    /// Insert a key-value pair. If the key already exists, the current value
    /// will be replaced with the given value.
    void insert(const size_t& key, const size_t& value);
    
    /// Erase the key-value pair associated with the key. If the key does not
    /// exist, do nothing.
    void erase(const size_t& key);
    
    /// Returns true if there are no entries, otherwise false.
    bool empty() const;
    
    /// Returns the number of entries.
    size_t size() const;
    
    /// Returns a handle to the key-value pair associated with a key, or 0 if
    /// there the key does not exist.
    size_t find(const size_t& key) const;
    
    /// Returns a handle to the key-value pair with the largest key that is less-than
    /// or equal to the given key, or 0 if the given key is less than the minimum or
    /// the tree is empty.
    size_t first_lower(const size_t& key) const;
    
    /// Returns the handle to the key-value pair with the next-smallest key to the
    /// given handle.
    size_t next(const size_t& x) const;
    
    /// Returns the key of a handle.
    inline size_t get_key(const size_t& x) const;
    
    /// Returns the value of a handle.
    inline size_t get_value(const size_t& x) const;

        
private:
    const static size_t NODE_SIZE = 5;
    const static int64_t KEY_OFFSET = 0;
    const static int64_t VALUE_OFFSET = 1;
    const static int64_t PARENT_OFFSET = 2;
    const static int64_t LEFT_CHILD_OFFSET = 3;
    const static int64_t RIGHT_CHILD_OFFSET = 4;
        
    PackedVector tree;
    size_t root = 0;
    size_t num_nodes = 0;
        
    inline size_t get_parent(size_t x) const;
    
    inline size_t get_left(size_t x) const;
    
    inline size_t get_right(size_t x) const;
        
    inline void set_key(size_t x, size_t val);
        
    inline void set_value(size_t x, size_t val);
        
    inline void set_left(size_t x, size_t y);
        
    inline void set_right(size_t x, size_t y);
        
    inline void set_parent(size_t x, size_t y);
        
    void left_rotate(size_t x);
        
    void right_rotate(size_t x);
        
    void splay(size_t x);
        
    void replace(size_t u, size_t v );
        
    size_t subtree_minimum(size_t u) const;
    
    size_t subtree_maximum(size_t u) const;
    
    size_t add_node(const size_t& key, const size_t& value);
        
    void delete_node(size_t x);
    
    void print_topology(std::ostream& out) const;
    void print_vector(std::ostream& out) const;
};
    
    
    
    
/// Inline functions
    
inline void PackedVector::set(const size_t& i, const uint64_t& value) {
    assert(i < filled);
        
    uint8_t width = vec.width();
    uint64_t mask = std::numeric_limits<uint64_t>::max() << width;
    while (mask & value) {
        width++;
        mask = std::numeric_limits<uint64_t>::max() << width;
    }
        
    if (width > vec.width()) {
        sdsl::int_vector<> wider_vec;
        wider_vec.width(width);
        wider_vec.resize(vec.size());
        for (size_t i = 0; i < filled; i++) {
            wider_vec[i] = vec[i];
        }
        vec = std::move(wider_vec);
    }
        
    vec[i] = value;
}
    
inline uint64_t PackedVector::get(const size_t& i) const {
    assert(i < filled);
    return vec[i];
}
    
inline void PackedVector::append(const uint64_t& value) {
    resize(filled + 1);
    set(filled - 1, value);
}
    
inline void PackedVector::pop() {
    resize(filled - 1);
}
    
inline void PackedVector::resize(size_t new_size) {
    if (new_size < filled) {
        size_t shrink_capacity = vec.size() / (factor * factor);
        if (new_size < shrink_capacity) {
            sdsl::int_vector<> tmp;
            tmp.width(vec.width());
            tmp.resize(new_size);
            for (size_t i = 0; i < new_size; i++) {
                tmp[i] = vec[i];
            }
            vec = std::move(tmp);
        }
    }
    else if (new_size > vec.size()) {
        size_t new_capacity = std::max<size_t>(size_t(vec.size() * factor) + 1, new_size);
        sdsl::int_vector<> tmp;
        tmp.width(vec.width());
        tmp.resize(new_capacity);
        for (size_t i = 0; i < filled; i++) {
            tmp[i] = vec[i];
        }
        vec = std::move(tmp);
    }
    filled = new_size;
}
    
inline size_t PackedVector::size() const {
    return filled;
}
    
inline bool PackedVector::empty() const {
    return filled == 0;
}

inline void PackedVector::clear() {
    vec.resize(0);
    vec.width(1);
    filled = 0;
}
    
inline bool PackedSplayTree::empty() const {
    return root == 0;
}
    
inline size_t PackedSplayTree::size() const {
    return num_nodes;
}
    
inline size_t PackedSplayTree::get_key(const size_t& x) const {
    return tree.get((x - 1) * NODE_SIZE + KEY_OFFSET);
}

inline size_t PackedSplayTree::get_value(const size_t& x) const {
    return tree.get((x - 1) * NODE_SIZE + VALUE_OFFSET);
}

inline size_t PackedSplayTree::get_parent(size_t x) const {
    return tree.get((x - 1) * NODE_SIZE + PARENT_OFFSET);
}

inline size_t PackedSplayTree::get_left(size_t x) const {
    return tree.get((x - 1) * NODE_SIZE + LEFT_CHILD_OFFSET);
}

inline size_t PackedSplayTree::get_right(size_t x) const {
    return tree.get((x - 1) * NODE_SIZE + RIGHT_CHILD_OFFSET);
}
    
inline void PackedSplayTree::set_key(size_t x, size_t val) {
    tree.set((x - 1) * NODE_SIZE + KEY_OFFSET, val);
}
    
inline void PackedSplayTree::set_value(size_t x, size_t val) {
    tree.set((x - 1) * NODE_SIZE + VALUE_OFFSET, val);
}
    
inline void PackedSplayTree::set_left(size_t x, size_t y) {
    tree.set((x - 1) * NODE_SIZE + LEFT_CHILD_OFFSET, y);
}
    
inline void PackedSplayTree::set_right(size_t x, size_t y) {
    tree.set((x - 1) * NODE_SIZE + RIGHT_CHILD_OFFSET, y);
}
    
inline void PackedSplayTree::set_parent(size_t x, size_t y) {
    tree.set((x - 1) * NODE_SIZE + PARENT_OFFSET, y);
}

inline size_t PackedDeque::internal_index(const size_t& i) const {
    assert(i < filled);
    return i < vec.size() - begin_idx ? begin_idx + i : i - (vec.size() - begin_idx);
}

inline void PackedDeque::set(const size_t& i, const uint64_t& value) {
    return vec.set(internal_index(i), value);
}

inline uint64_t PackedDeque::get(const size_t& i) const {
    return vec.get(internal_index(i));
}

inline void PackedDeque::append_front(const uint64_t& value) {
    if (filled == vec.size()) {
        size_t new_capacity = size_t(factor * vec.size()) + 1;
        PackedVector new_vec;
        new_vec.resize(new_capacity);
        
        new_vec.set(0, value);
        for (size_t i = 0; i < filled; i++) {
            new_vec.set(i + 1, get(i));
        }
        
        vec = std::move(new_vec);
        begin_idx = 0;
    }
    else {
        if (begin_idx == 0) {
            begin_idx = vec.size() - 1;
        }
        else {
            begin_idx--;
        }
        vec.set(begin_idx, value);
    }
    
    filled++;
}

inline void PackedDeque::append_back(const uint64_t& value) {
    if (filled == vec.size()) {
        size_t new_capacity = size_t(factor * vec.size()) + 1;
        PackedVector new_vec;
        new_vec.resize(new_capacity);
        
        for (size_t i = 0; i < filled; i++) {
            new_vec.set(i, get(i));
        }
        new_vec.set(filled, value);
        
        vec = std::move(new_vec);
        begin_idx = 0;
        filled++;
    }
    else {
        filled++;
        vec.set(internal_index(filled - 1), value);
    }
}
    
inline void PackedDeque::contract() {
    size_t shrink_capacity = vec.size() / (factor * factor);
    if (filled <= shrink_capacity) {
        PackedVector new_vec;
        new_vec.resize(filled);
        for (size_t i = 0; i < filled; i++) {
            new_vec.set(i, get(i));
        }
        
        vec = std::move(new_vec);
        begin_idx = 0;
    }
}

inline void PackedDeque::pop_front() {
    begin_idx++;
    if (begin_idx == vec.size()) {
        begin_idx = 0;
    }
    filled--;
    contract();
}

inline void PackedDeque::pop_back() {
    filled--;
    contract();
}

inline size_t PackedDeque::size() const {
    return filled;
}

inline bool PackedDeque::empty() const {
    return filled == 0;
}
    
inline void PackedDeque::clear() {
    vec.clear();
    filled = 0;
    begin_idx = 0;
}
    
inline void PagedVector::set(const size_t& i, const uint64_t& value) {
    assert(i < filled);
    uint64_t anchor = anchors.get(i / page_size);
    if (anchor == 0) {
        // this page does not have a non-zero anchor yet, use this one
        anchors.set(i / page_size, value);
        anchor = value;
    }
    pages[i / page_size].set(i % page_size, to_diff(value, anchor));
}

inline uint64_t PagedVector::get(const size_t& i) const {
    assert(i < filled);
    return from_diff(pages[i / page_size].get(i % page_size),
                     anchors.get(i / page_size));
}

inline void PagedVector::append(const uint64_t& value) {
    if (filled % page_size == 0) {
        // init a new page and a new anchor
        pages.emplace_back();
        pages.back().resize(page_size);
        anchors.append(0);
    }
    
    // use the logic in set to choose anchor and diff
    filled++;
    set(filled - 1, value);
}

inline void PagedVector::pop() {
    filled--;
    if (filled % page_size == 0) {
        // we've emptied a page, remove it
        pages.pop_back();
        anchors.pop();
    }
}

inline void PagedVector::resize(size_t new_size) {
    // how many pages does this require?
    size_t num_pages = new_size > 0 ? (new_size - 1) / page_size + 1 : 0;
    
    anchors.resize(num_pages);
    // add pages if necessary
    while (num_pages > pages.size()) {
        pages.emplace_back();
        pages.back().resize(page_size);
    }
    // remove pages if necessary
    pages.resize(num_pages);
    
    filled = new_size;
}

inline size_t PagedVector::size() const {
    return filled;
}

inline bool PagedVector::empty() const {
    return filled == 0;
}

inline void PagedVector::clear() {
    pages.clear();
    anchors.clear();
    filled = 0;
}
    
inline uint64_t PagedVector::to_diff(const uint64_t& value, const uint64_t& anchor) const {
    // leaves 0 unchanged, encodes other values as a difference from the anchor value
    // with a bijection into the positive integers as follows:
    // difference  0  1  2  3 -1  4  5  6  7 -2  8  9 10 11 -3 ...
    // integer     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 ...
    // the goal here is use smaller integers to maintain low bit-width, allowing 0 as
    // a sentinel. the bijection is biased to encode positive differences as smaller bit-
    // width integers since anchors are taken from the beginning of their page in the
    // vector and we expect most vectors to be mostly increasing
    
    if (value == 0) {
        return 0;
    }
    else if (value >= anchor) {
        uint64_t raw_diff = value - anchor;
        return raw_diff + raw_diff / 4 + 1;
    }
    else {
        return 5 * (anchor - value);
    }
}

inline uint64_t PagedVector::from_diff(const uint64_t& diff, const uint64_t& anchor) const {
    if (diff == 0) {
        return 0;
    }
    else if (diff % 5 == 0) {
        return anchor - diff / 5;
    }
    else {
        return anchor + diff - diff / 5 - 1;
    }
}
}



#endif /* dynamic_structs_hpp */
