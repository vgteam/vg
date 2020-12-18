#ifndef CONTRACTINGGRAPH_HPP
#define CONTRACTINGGRAPH_HPP

#include "algorithms/min_cut_graph.hpp"
#include <structures/union_find.hpp>

namespace vg{
    using namespace std;
    using vg::algorithms::Graph;
    using namespace structures;

class ContractingGraph{   
    Graph graph;
    size_t V;   
    
    
public:
    UnionFind uf = UnionFind(V, true);

    ContractingGraph(Graph graph, size_t V);

    void contract(size_t random_node, size_t other_node);

    unordered_map<size_t, size_t>  get_edges(size_t group_num);


    vector<size_t> get_nodes();
    //looping through original nodes and sending them to find_group() and determining which nodes are heads
    //return all heads 
    //if node = find_group(node) then its a head otherwise not 

    

    
        
};


}




#endif