#ifndef HASH_MAP_H
#define HASH_MAP_H

#include "sparsehash/sparse_hash_map"
#include "sparsehash/dense_hash_map"


// http://stackoverflow.com/questions/4870437/pairint-int-pair-as-key-of-unordered-map-issue#comment5439557_4870467
// https://github.com/Revolutionary-Games/Thrive/blob/fd8ab943dd4ced59a8e7d1e4a7b725468b7c2557/src/util/pair_hash.h
// taken from boost
#ifndef OVERLOAD_PAIR_HASH
#define OVERLOAD_PAIR_HASH
namespace std {
template <typename A, typename B>
struct hash<pair<A,B> > {
    size_t operator()(const pair<A,B>& x) const {
        std::size_t hashA = std::hash<A>()(x.first);
        std::size_t hashB = std::hash<B>()(x.second);
        return hashA ^ (hashB + 0x9e3779b9 + (hashA << 6) + (hashA >> 2));
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
