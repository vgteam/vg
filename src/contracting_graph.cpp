
#include "contracting_graph.hpp"
#include <structures/union_find.hpp>


namespace vg{
    using namespace std;


    ContractingGraph::ContractingGraph(Graph graph, int n_nodes, vector<int> adj_list[], structures::UnionFind uf)
        :graph(graph), n_nodes(n_nodes), adj(adj), uf(uf){

    }
    unordered_map<size_t, size_t> ContractingGraph::get_edges(size_t group_num){

        vector<size_t> group_nodes = get_nodes(group_num);

                
    }
    //returns a vector of the indices in the same group as index i
    vector<size_t> ContractingGraph::get_nodes(size_t node_i){

         return uf.group(node_i); 
    }



}