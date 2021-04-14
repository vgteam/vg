#ifndef STRUCTURES_SPARSE_UNION_FIND_HPP_INCLUDED
#define STRUCTURES_SPARSE_UNION_FIND_HPP_INCLUDED

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <structures/union_find.hpp>

namespace vg {

    using namespace std;
    using namespace structures;

    class SparseUnionFind : public UnionFind{

        public:
        vector<size_t> node_ids;
        unordered_map<size_t, size_t> sparse_to_dense;//incoming 
        unordered_map<size_t, size_t> dense_to_sparse;//outgoing
        

        SparseUnionFind(bool include_children, vector<size_t> node_ids);

        /// Destructor
        ~SparseUnionFind();
        
        /// Returns the number of indices in the UnionFind
        size_t size();
        
        /// Returns the group ID that index i belongs to (can change after calling union)
        size_t find_group(size_t i);
        
        /// Merges the group containing index i with the group containing index j
        void union_groups(size_t i, size_t j);
        
        /// Returns the size of the group containing index i
        size_t group_size(size_t i);
        
        /// Returns a vector of the indices in the same group as index i
        vector<size_t> group(size_t i);
        
        /// Returns all of the groups, each in a separate vector
        vector<vector<size_t>> all_groups();

        private:
        bool include_children;

    };
  

}
#endif