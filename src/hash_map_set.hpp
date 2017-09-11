#ifndef XG_HASH_MAP_SET_H
#define XG_HASH_MAP_SET_H

#include "sparsehash/sparse_hash_map"
#include "sparsehash/sparse_hash_set"

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

namespace xg {

using google::sparse_hash_map;
using google::sparse_hash_set;

template<typename K, typename V>
    class hash_map : public sparse_hash_map<K,V>
    {
public:
    hash_map() {
        this->set_deleted_key(-2);
    }
};

template<typename K, typename V>
    class string_hash_map : public sparse_hash_map<K,V>
{
public:
    string_hash_map() {
        this->set_deleted_key("");
    }
};

template<typename K, typename V>
    class pair_hash_map : public sparse_hash_map<K,V,std::hash<K> >
{
public:
    pair_hash_map() {
        this->set_deleted_key(K(-2, -2));
    }
};

template<typename K, typename V>
    class hash_map<K*,V> : public sparse_hash_map<K*,V>
{
public:
    hash_map() {
        this->set_deleted_key((K*)(0));
    }
};



template<typename K>
class hash_set : public sparse_hash_set<K>
    {
public:
    hash_set() {
        this->set_deleted_key(-2);
    }
};

template<typename K>
class string_hash_set : public sparse_hash_set<K>
{
public:
    string_hash_set() {
        this->set_deleted_key("");
    }
};

template<typename K>
    class pair_hash_set : public sparse_hash_set<K,std::hash<K> >
{
public:
    pair_hash_set() {
        this->set_deleted_key(K(-2, -2));
    }
};

template<typename K>
    class hash_set<K*> : public sparse_hash_set<K*>
{
public:
    hash_set() {
        this->set_deleted_key((K*)(0));
    }
};

}

#endif
