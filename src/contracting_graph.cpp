#include "contracting_graph.hpp"
#include <structures/union_find.hpp>
#include <iostream>


#define debug
namespace vg{
    using namespace std;
    using namespace structures;

    ContractingGraph::ContractingGraph(Graph graph, size_t V)
        :graph(graph), V(V){

    }
  
    void ContractingGraph::contract(size_t random_node, size_t other_node){

        uf.union_groups(random_node, other_node);

    }
    unordered_map<size_t, size_t> ContractingGraph::get_edges(size_t group_num){
        
        //create container to keep group and edge_totals
        unordered_map<size_t, size_t> group_edges; 

        //get contents of group 
        vector<size_t> group_nodes = uf.group(group_num);

        //if an adj_node exists in group_nodes then it is a contracted edge, we treat it as a special case  
        for(size_t i = 0; i<group_nodes.size(); i++){
            int member = group_nodes[i];
#ifdef debug
            cout << "============================================================================= " << endl;    
            cout << "group num " << group_num+1 << endl;       
#endif
            
            if(graph.nodes[member].edges.size() == 0){
                break;
            }
            for (size_t j = 0; j < graph.nodes[member].edges.size(); j++){
#ifdef debug
            cout << "============================================================================= " << endl;    
            cout << "node member from super group " << member+1 << endl;       
            cout << "node member other " << graph.nodes[member].edges[j].other+1 << endl;  
#endif
                size_t connecting_node = graph.nodes[member].edges[j].other;

                // check if the connecting node is contracted with other nodes
                size_t connecting_node_group_id = uf.find_group(connecting_node);
                if(connecting_node_group_id == group_num){
                    cout << "continue" <<endl;
                    continue;
                }else{
#ifdef debug   
            cout << "connecting node group id " << connecting_node_group_id << endl;
            cout << "node edges prev " << group_edges[connecting_node] << endl;  
            cout << "node edges new" << graph.nodes[connecting_node].edges[j].weight << endl;  
#endif
                //if it doesn't exist add, otherwise add weight to existing value
                group_edges[connecting_node] += graph.nodes[connecting_node].edges[j].weight;
#ifdef debug     
            
            cout << "node edges total" << group_edges[connecting_node+1] << endl;  
            cout << "============================================================================= " << endl;  
#endif
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