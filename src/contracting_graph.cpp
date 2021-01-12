#include "contracting_graph.hpp"
#include <structures/union_find.hpp>
#include <iostream>


// #define debug
namespace vg{
    using namespace std;
    using namespace structures;

    ContractingGraph::ContractingGraph(Graph graph, size_t V)
        :graph(graph), V(V){
#ifdef debug
            cout << "original graph" <<endl; 
            cout << "============================================================================= " << endl;
                   
            for(size_t i = 0; i < graph.nodes.size(); i++){
                cout << "node: "<< i  << ", weight: "<<graph.nodes[i].weight << endl;
                for (size_t j = 0; j < graph.nodes[i].edges.size(); j++){
                    cout << "edge "<< i  << "->" << graph.nodes[i].edges[j].other <<", weight: " << graph.nodes[i].edges[j].weight << endl;
                }
            }
            cout << "============================================================================= " << endl;
#endif

    }
  
    void ContractingGraph::contract(size_t random_node, size_t other_node){

        uf.union_groups(random_node, other_node);

    }
    unordered_map<size_t, size_t> ContractingGraph::get_edges(size_t group_num){

        //create container to keep group and edge_totals
        unordered_map<size_t, size_t> group_edges; 

        //get contents of group 
        vector<size_t> group_nodes = uf.group(group_num);
#ifdef debug
            cout << "============================================================================= " << endl;    
            cout << "group num " << group_num << endl;       
#endif

        //if an adj_node exists in group_nodes then it is a contracted edge, we treat it as a special case  
        for(size_t i = 0; i<group_nodes.size(); i++){
            int member = group_nodes[i];

            
            if(graph.nodes[member].edges.size() == 0){
                cout<<"disconnected edge" << endl;
                break;
            }
            for (size_t j = 0; j < graph.nodes[member].edges.size(); j++){
#ifdef debug
            cout << "============================================================================= " << endl;    
            cout << "edge " << member  << "->"<< graph.nodes[member].edges[j].other  << endl;        
#endif
                size_t connecting_node = graph.nodes[member].edges[j].other;

                // check if the connecting node is contracted with other nodes
                size_t connecting_node_group_id = uf.find_group(connecting_node);

                // avoid double counting edges btween members of group num aka contracted edge
                if(connecting_node_group_id == group_num){
#ifdef debug     
                    cout << "continue" <<endl;
#endif      
                    continue;
                }else{
#ifdef debug   
            cout << i << j << endl;
            cout << "connecting node" << connecting_node <<endl;
            cout << "connecting node group id " << connecting_node_group_id << endl;
            cout << "node edges prev " << group_edges[connecting_node_group_id] << endl;  
            cout << "node edges new " << graph.nodes[member].edges[j].weight << endl;  
#endif
                //if it doesn't exist add, otherwise add weight to existing value
                // group edges indexed by connecting node group number 
                group_edges[connecting_node_group_id] += graph.nodes[member].edges[j].weight;
#ifdef debug     
            
            cout << "node edges total group id " << group_edges[connecting_node_group_id] << endl;  
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

    vector<vector<size_t>> ContractingGraph::get_disjoint_sets(){
        vector<vector<size_t>> to_return;
        to_return = uf.all_groups();
        return to_return;
    }



}