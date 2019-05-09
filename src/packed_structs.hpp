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
#include <iostream>
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
    
    /// Construct from contents in a stream
    PackedVector(istream& in);
    
    /// Move constructor
    PackedVector(PackedVector&& other) = default;
    /// Move assignment operator
    PackedVector& operator=(PackedVector&& other) = default;
        
    /// Destructor
    ~PackedVector();
    
    /// Clear current contents and load from contents in a stream
    void deserialize(istream& in);
    
    /// Output contents to a stream
    void serialize(ostream& out) const ;
    
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
    inline void resize(const size_t& new_size);
    
    /// If necessary, expand capacity so that the given number of entries can
    /// be included in the vector without reallocating. Never shrinks capacity.
    inline void reserve(const size_t& future_size);
        
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
 * Compression is also optimized for vectors that are mostly (but not necessarily
 * exclusively) increasing.
 */
class PagedVector {
public:
    
    /// Construct and set page size (starts empty)
    PagedVector(size_t page_size);
    
    /// Construct from contents in a stream
    PagedVector(istream& in);
    
    /// Move constructor
    PagedVector(PagedVector&& other) = default;
    /// Move assignment operator
    PagedVector& operator=(PagedVector&& other) = default;
    
    // Destructor
    ~PagedVector();
    
    /// Clear current contents and load from contents in a stream
    void deserialize(istream& in);
    
    /// Output contents to a stream
    void serialize(ostream& out) const ;
    
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
    inline void resize(const size_t& new_size);
    
    /// If necessary, expand capacity so that the given number of entries can
    /// be included in the vector without reallocating. Never shrinks capacity.
    inline void reserve(const size_t& future_size);
    
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
    /// Construct empty
    PackedDeque(void);
    /// Construct from contents in a stream
    PackedDeque(istream& in);
    
    /// Move constructor
    PackedDeque(PackedDeque&& other) = default;
    /// Move assignment operator
    PackedDeque& operator=(PackedDeque&& other) = default;
    
    /// Destructor
    ~PackedDeque(void);
    
    /// Clear current contents and load from contents in a stream
    void deserialize(istream& in);
    
    /// Output contents to a stream
    void serialize(ostream& out) const ;
    
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
    
    /// If necessary, expand capacity so that the given number of entries can
    /// be included in the deque without reallocating. Never shrinks capacity.
    inline void reserve(const size_t& future_size);
    
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
    
    
    
/// Inline functions
    
/////////////////////
/// PackedVector
/////////////////////
    
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
    
inline void PackedVector::resize(const size_t& new_size) {
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
        reserve(new_capacity);
    }
    filled = new_size;
}
    
inline void PackedVector::reserve(const size_t& future_size) {
    if (future_size > vec.size()) {
        sdsl::int_vector<> tmp;
        tmp.width(vec.width());
        tmp.resize(future_size);
        for (size_t i = 0; i < filled; i++) {
            tmp[i] = vec[i];
        }
        vec = std::move(tmp);
    }
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
    
/////////////////////
/// PackedDeque
/////////////////////
    
    
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
    
inline void PackedDeque::reserve(const size_t& future_size) {
    if (future_size > vec.size()) {
        PackedVector new_vec;
        new_vec.resize(future_size);
        
        for (size_t i = 0; i < filled; i++) {
            new_vec.set(i, get(i));
        }
        vec = std::move(new_vec);
        begin_idx = 0;
    }
}

inline void PackedDeque::append_front(const uint64_t& value) {
    // expand capacity if necessary
    if (filled == vec.size()) {
        size_t new_capacity = size_t(factor * vec.size()) + 1;
        reserve(new_capacity);
    }
    
    // update the pointer to the front
    if (begin_idx == 0) {
        begin_idx = vec.size() - 1;
    }
    else {
        begin_idx--;
    }
    // update the pointer to the back
    filled++;
    
    // set the value
    vec.set(internal_index(0), value);
}

inline void PackedDeque::append_back(const uint64_t& value) {
    // expand capacity if necessary
    if (filled == vec.size()) {
        size_t new_capacity = size_t(factor * vec.size()) + 1;
        reserve(new_capacity);
    }
    
    // update the pointer to the back
    filled++;
    
    // set the value
    vec.set(internal_index(filled - 1), value);
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
    // update the pointer to the beginning
    begin_idx++;
    if (begin_idx == vec.size()) {
        begin_idx = 0;
    }
    // update the pointer to the end
    filled--;
    
    // shrink if necessary
    contract();
}

inline void PackedDeque::pop_back() {
    // update the pointer to the end
    filled--;
    
    // shrink if necessary
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
    
/////////////////////
/// PagedVector
/////////////////////
    
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
    if (filled == pages.size() * page_size) {
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
    while (filled + page_size <= pages.size() * page_size) {
        // the final page is unused now, remove it
        pages.pop_back(); // TODO: this won't resize since it's an STL vector
        anchors.pop();
    }
}

inline void PagedVector::resize(const size_t& new_size) {
    if (new_size < filled) {
        // shrink down to the number of pages we would need
        size_t num_pages = new_size == 0 ? 0 : (new_size - 1) / page_size + 1;
        anchors.resize(num_pages);
        pages.resize(num_pages);
    }
    else if (new_size > filled) {
        // make sure we capacity for this many elements
        reserve(new_size);
    }
    filled = new_size;
}
    
inline void PagedVector::reserve(const size_t& future_size) {
    if (future_size > pages.size() * page_size) {
        // how many pages does this require?
        size_t num_pages = (future_size - 1) / page_size + 1;
        // note: we don't need to worry about underflow b/c previous condition
        // implies future_size > 0
        
        // expand anchor and pages vectors out to the capacity of the number of pages
        anchors.reserve(num_pages);
        pages.reserve(num_pages);
        
        // add the anchors and fixed-width pages in this
        anchors.resize(num_pages);
        while (num_pages > pages.size()) {
            pages.emplace_back();
            pages.back().resize(page_size);
        }
    }
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
    // with a reversible mapping into the positive integers as follows:
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
    // convert backward from the transformation described in to_diff
    
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
