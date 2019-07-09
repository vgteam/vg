#ifndef VG_HASH_MAP_HPP_INCLUDED
#define VG_HASH_MAP_HPP_INCLUDED

#include <cstdint>
#include <tuple>
#include <type_traits>

#include "wang_hash.hpp"

// Comment out to use sparse_hash_map and sparse_hash_set instead of
// dense_hash_map and dense_hash_set.
//#define USE_DENSE_HASH

#ifdef USE_DENSE_HASH
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#else
#include <sparsepp/spp.h>
#endif


// http://stackoverflow.com/questions/4870437/pairint-int-pair-as-key-of-unordered-map-issue#comment5439557_4870467
// https://github.com/Revolutionary-Games/Thrive/blob/fd8ab943dd4ced59a8e7d1e4a7b725468b7c2557/src/util/pair_hash.h
// taken from boost
#ifndef OVERLOAD_PAIR_HASH
#define OVERLOAD_PAIR_HASH
namespace std {
namespace
{
    
    // Code from boost
    // Reciprocal of the golden ratio helps spread entropy
    //     and handles duplicates.
    // See Mike Seymour in magic-numbers-in-boosthash-combine:
    //     http://stackoverflow.com/questions/4948780
    
    template <class T>
    inline void hash_combine(size_t& seed, T const& v)
    {
        seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }
    
    // Recursive template code derived from Matthieu M.
    template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
    struct HashValueImpl
    {
        static void apply(size_t& seed, Tuple const& tuple)
        {
            HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
            hash_combine(seed, std::get<Index>(tuple));
        }
    };
    
    template <class Tuple>
    struct HashValueImpl<Tuple,0>
    {
        static void apply(size_t& seed, Tuple const& tuple)
        {
            hash_combine(seed, std::get<0>(tuple));
        }
    };
}
    
template <typename A, typename B>
struct hash<pair<A,B> > {
    size_t operator()(const pair<A,B>& x) const {
        size_t hash_val = std::hash<A>()(x.first);
        hash_combine(hash_val, x.second);
        return hash_val;
    }
};

// from http://stackoverflow.com/questions/7110301/generic-hash-for-tuples-in-unordered-map-unordered-set
template <typename ... TT>
struct hash<std::tuple<TT...>>
{
    size_t
    operator()(std::tuple<TT...> const& tt) const
    {
        size_t seed = 0;
        HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
        return seed;
    }
    
};
}
#endif  // OVERLOAD_PAIR_HASH


namespace vg {

// We need this second type for enable_if-based specialization
template<typename T, typename ImplementationMatched = void>
struct wang_hash;

// We can hash pointers
template<typename T>
struct wang_hash<T*> {
    size_t operator()(const T* pointer) const {
        return wang_hash_64(reinterpret_cast<size_t>(pointer));
    }
};

// We can hash any integer that can be implicitly widened to size_t.
// This covers 32 bit ints (which we need to be able to hash on Mac) and 64 bit ints
// This also coveres bools.
// See <https://stackoverflow.com/a/42679086>
template<typename T>
struct wang_hash<T, typename std::enable_if<std::is_integral<T>::value>::type> {
    size_t operator()(const T& x) const {
        static_assert(sizeof(T) <= sizeof(size_t), "widest hashable type is size_t");
        return wang_hash_64(static_cast<size_t>(x));
    }
};

// We can hash pairs
template<typename A, typename B>
struct wang_hash<std::pair<A, B>> {
    size_t operator()(const std::pair<A, B>& x) const {
        size_t hash_val = wang_hash<A>()(x.first);
        hash_val ^= wang_hash<B>()(x.second) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        return hash_val;
    }
};


// Replacements for std::unordered_map.

template<typename K, typename V>
#ifdef USE_DENSE_HASH
class hash_map : public google::dense_hash_map<K, V, wang_hash<K>>
#else
class hash_map : public spp::sparse_hash_map<K, V, wang_hash<K>>
#endif
{
#ifdef USE_DENSE_HASH
public:
    hash_map() {
        this->set_empty_key(-1);
    }
#endif
};

template<typename K, typename V>
#ifdef USE_DENSE_HASH
class string_hash_map : public google::dense_hash_map<K, V>
#else
class string_hash_map : public spp::sparse_hash_map<K, V>
#endif
{
#ifdef USE_DENSE_HASH
public:
    string_hash_map() {
        this->set_empty_key(" ");
    }
#endif
};

template<typename K, typename V>
#ifdef USE_DENSE_HASH
class pair_hash_map : public google::dense_hash_map<K, V, wang_hash<K>>
#else
class pair_hash_map : public spp::sparse_hash_map<K, V, wang_hash<K>>
#endif
{
#ifdef USE_DENSE_HASH
public:
    pair_hash_map() {
        this->set_empty_key(K(-1, -1));
    }
#endif
};

template<typename K, typename V>
#ifdef USE_DENSE_HASH
class hash_map<K*, V> : public google::dense_hash_map<K*, V, wang_hash<K*>>
#else
class hash_map<K*, V> : public spp::sparse_hash_map<K*, V, wang_hash<K*>>
#endif
{
#ifdef USE_DENSE_HASH
public:
    hash_map() {
        this->set_empty_key((K*)(~0));
    }
#endif
};


// Replacements for std::unordered_set.

template<typename K>
#ifdef USE_DENSE_HASH
class hash_set : public google::dense_hash_set<K, wang_hash<K>>
#else
class hash_set : public spp::sparse_hash_set<K, wang_hash<K>>
#endif
    {
#ifdef USE_DENSE_HASH
public:
    hash_set() {
        this->set_empty_key(-1);
    }
#endif
};

template<typename K>
#ifdef USE_DENSE_HASH
class string_hash_set : public google::dense_hash_set<K>
#else
class string_hash_set : public spp::sparse_hash_set<K>
#endif
{
#ifdef USE_DENSE_HASH
public:
    string_hash_set() {
        this->set_empty_key(" ");
    }
#endif
};

template<typename K>
#ifdef USE_DENSE_HASH
class pair_hash_set : public google::dense_hash_set<K, wang_hash<K>>
#else
class pair_hash_set : public spp::sparse_hash_set<K, wang_hash<K>>
#endif
{
#ifdef USE_DENSE_HASH
public:
    pair_hash_set() {
        this->set_empty_key(K(-1, -1));
    }
#endif
};

template<typename K>
#ifdef USE_DENSE_HASH
class hash_set<K*> : public google::dense_hash_set<K*, wang_hash<K*>>
#else
class hash_set<K*> : public spp::sparse_hash_set<K*, wang_hash<K*>>
#endif
{
#ifdef USE_DENSE_HASH
public:
    hash_set() {
        this->set_empty_key((K*)(~0));
    }
#endif
};


}   // namespace vg

#endif
