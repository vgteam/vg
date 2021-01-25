#ifndef STRUCTURES_SPARSE_UNION_FIND_HPP_INCLUDED
#define STRUCTURES_SPARSE_UNION_FIND_HPP_INCLUDED

#include <vector>
#include <unordered_set>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include "structures/union_find.hpp"

namespace structures {

using namespace std;

class SparseUnionFind : public UnionFind{

    public:
    vector<size_t> IDs;

    SparseUnionFind(size_t size, bool include_children, vector<size_t> node_ids);

    find_group(size_t i);

    private:
    struct UFNode;
    vector<UFNode> uf_nodes;
    bool include_children;
    unordered_map<size_t, size_t> sparse_to_dense;//incoming 
    unordered_map<size_t, size_t> dense_to_sparse;//outgoing



};

}
#endif