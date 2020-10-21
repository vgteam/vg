#ifndef CONTRACTINGGRAPH_HPP
#define CONTRACTINGGRAPH_HPP

#include "algorithms/min_cut_graph.hpp"


namespace vg{
    using namespace std;
    using vg::algorithms::Graph;

    class ContractingGraph{
            
            Graph graph;
            int n_nodes;
            vector<int> adj[];
            structures::UnionFind uf;

            public:
            ContractingGraph(Graph graph, int n_nodes, vector<int> adj[], structures::UnionFind uf);

            unordered_map<size_t, size_t>  get_edges(size_t group_num);

            vector<size_t> get_nodes(size_t node_at_index_i);
            
        };



}




#endif