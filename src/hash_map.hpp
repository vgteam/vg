#ifndef VG_HASH_MAP_HPP_INCLUDED
#define VG_HASH_MAP_HPP_INCLUDED

#include "sparsehash/sparse_hash_map"
#include "sparsehash/dense_hash_map"

#include <tuple>

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
    inline void hash_combine(std::size_t& seed, T const& v)
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
#endif

namespace vg {

using google::sparse_hash_map;
using google::dense_hash_map;

// comment out to use sparse_hash_map instead
#define USE_DENSE_HASH

template<typename K, typename V>
#ifdef USE_DENSE_HASH
    class hash_map : public dense_hash_map<K,V>
#else
    class hash_map : public sparse_hash_map<K,V>
#endif
    {
public:
    hash_map() {
#ifdef USE_DENSE_HASH
        this->set_empty_key(-1);
#endif
        this->set_deleted_key(-2);
    }
};

template<typename K, typename V>
#ifdef USE_DENSE_HASH
    class string_hash_map : public dense_hash_map<K,V>
#else
    class string_hash_map : public sparse_hash_map<K,V>
#endif
{
public:
    string_hash_map() {
#ifdef USE_DENSE_HASH
        this->set_empty_key(" ");
#endif
        this->set_deleted_key("");
    }
};

template<typename K, typename V>
#ifdef USE_DENSE_HASH
    class pair_hash_map : public dense_hash_map<K,V,std::hash<K> >
#else
    class pair_hash_map : public sparse_hash_map<K,V,std::hash<K> >
#endif
{
public:
    pair_hash_map() {
#ifdef USE_DENSE_HASH
        this->set_empty_key(K(-1, -1));
#endif
        this->set_deleted_key(K(-2, -2));
    }
};

template<typename K, typename V>
#ifdef USE_DENSE_HASH
    class hash_map<K*,V> : public dense_hash_map<K*,V>
#else
    class hash_map<K*,V> : public sparse_hash_map<K*,V>
#endif
{
public:
    hash_map() {
#ifdef USE_DENSE_HASH
        this->set_empty_key((K*)(~0));
#endif
        this->set_deleted_key((K*)(0));
    }
};

}

#endif
