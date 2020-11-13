#include "contracting_graph.hpp"
#include <structures/union_find.hpp>


namespace vg{
    using namespace std;
    using namespace structures;


    ContractingGraph::ContractingGraph(Graph graph, int n_nodes, vector<vector<int>> adj_list, UnionFind uf)
        :graph(graph), n_nodes(n_nodes), adj_list(adj_list), uf(uf){
            

    }

    unordered_map<size_t, size_t> ContractingGraph::get_edges(size_t group_num){
        
        //create container to keep group and edge_totals
        unordered_map<size_t, size_t> group_edges; 

        //get contents of group 
        vector<size_t> group_nodes = uf.group(group_num);
        

        // the connecting nodes in adjacency list corresponding to node[group_num]
        // note: will remove the group_nodes from connecting nodes (contracted edges)
        // vector<int> connecting_nodes  = adj_list[group_num];

        
        //if an adj_node exists in group_nodes then it is a contracted edge, we treat it as a special case  
        for(size_t i = 0; i<group_nodes.size(); i++){
            int member = group_nodes[i];
            for (size_t j = 0; j < graph.nodes[member].edges.size(); j++){
                int connecting_node = graph.nodes[member].edges[j].other;
                size_t connecting_node_group_id = uf.find_group(connecting_node);

                if(connecting_node_group_id == group_num){
                    continue;
                }else{

                //if it doesn't exist add, otherwise add weight to existing value
                group_edges[connecting_node] += graph.nodes[connecting_node].edges[j].weight;

                }
            }

        }

        return group_edges;
        
    }
    
    vector<size_t> ContractingGraph::get_nodes(){
        
        //holds indices of nodes that are heads in the contracting graph 
        vector<size_t> heads;
        
        //loop through graph nodes and determine which nodes are heads of the group
        for(int i =0; i < graph.nodes.size(); i++){
            size_t group_head = uf.find_group(i);
            if(i == group_head){
                heads.push_back(i);

            }
        }
       
        return heads;
    }



}