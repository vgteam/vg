
#include "sparse_union_find.hpp"

namespace vg {

    using namespace std;

    SparseUnionFind::SparseUnionFind(bool include_children, vector<size_t> node_ids) :
        include_children(include_children), node_ids(node_ids), UnionFind(node_ids.size(), include_children) {
        for (size_t i = 0; i < node_ids.size(); i++) {
            sparse_to_dense[node_ids[i]] = i; 
            dense_to_sparse[i]= node_ids[i]; 
        }        
    }
    SparseUnionFind::~SparseUnionFind(){
        //nothing to do
    }

    size_t SparseUnionFind::find_group(size_t i){
       size_t dense_group_id = UnionFind::find_group(i);
       
       size_t sparse_group_id = dense_to_sparse[dense_group_id];

       return sparse_group_id;

    }

    vector<size_t> SparseUnionFind::group(size_t i) {
        
        vector<size_t> to_return = UnionFind::group(sparse_to_dense.at(i));
       
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

    vector<vector<size_t>> SparseUnionFind::all_groups() {
        vector<vector<size_t>> to_return = UnionFind::all_groups();
        
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


}