
#include "structures/sparse_union_find.hpp"

namespace structures {

    using namespace std;

    SparseUnionFind::SparseUnionFind(){size_t size, bool include_children, vector<size_t> node_ids) :
        IDs(IDs), include_children(include_children), node_ids(node_ids){
        
        uf_nodes.reserve(size);
        for (size_t i = 0; i < size; i++) {
            uf_nodes.emplace_back(i);
            sparse_to_dense[node_ids[i]] = i; 
            dense_to_sparse[i]= node_ids[i]; 
        }
        
        
    }


    SparseUnionFind::find_group(size_t i){
       size_t dense_group_id = UnionFind::find_group(size_t i);
       
       size_t sparse_group_id = dense_to_sparse[dense_group_id];

       return sparse_group_id;

    }

    vector<size_t> UnionFind::group(size_t i) {
        assert(include_children);
        vector<size_t> to_return;//contains dense node_ids
        // go to head of group
        vector<size_t> stack{find_group(i)};
        // traverse tree downwards to find all indices in group
        while (!stack.empty()) {
            size_t curr = stack.back();
            stack.pop_back();
            to_return.push_back(curr);
            unordered_set<size_t>& children = uf_nodes[curr].children;
            for (size_t child : children) {
                stack.push_back(child);
            }
        }

        vector<size_t> sparse_to_return;
        //traverse to_return and retrieve ids from dense_to_sparse
        for(size_t i =0; i < to_return.size(); i++){
            //iterate through vector 
            size_t node_to_lookup = to_return[i];
            //do a lookup
            size_t translated_node_id= dense_to_sparse[node_to_lookup];
            //push to vector
            sparse_to_return.push_back(translated_node_id);
        }
        return sparse_to_return;
    }

    vector<vector<size_t>> UnionFind::all_groups() {
        vector<vector<size_t>> to_return(uf_nodes.size());
        for (size_t i = 0; i < uf_nodes.size(); i++) {
            to_return[find_group(i)].push_back(i);
        }
        auto new_end = std::remove_if(to_return.begin(), to_return.end(),
                                    [](const vector<size_t>& grp) { return grp.empty(); });
        to_return.resize(new_end - to_return.begin());

        vector<vector<size_t>> sparse_to_return;

        //translate from dense to sparse
        for(size_t i = 0; i<to_return.size(); i++){
            vector<size_t> sparse_groups;
            for(size_t j = 0; j<to_return[i].size(); j++){
                size_t dense_node = to_return[i][j];
                //lookup
                size_t translated_node = dense_to_sparse[dense_node];
                //push back 
                sparse_groups.push_back(translated_node);
            }
            sparse_to_return.push_back(sparse_groups);

        }
        return sparse_to_return;
    }


};