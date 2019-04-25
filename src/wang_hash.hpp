#ifndef VG_WANG_HASH_HPP_INCLUDED
#define VG_WANG_HASH_HPP_INCLUDED

#include <cstddef>

namespace vg {

/// Thomas Wang's integer hash function. In many implementations, std::hash
/// is identity function for integers, which leads to performance issues.
inline size_t wang_hash_64(size_t key) {
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

}   // namespace vg

#endif
