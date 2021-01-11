/**
 * \file min_cut_graph.cpp
 *
 * Contains implementation of min_cut_graph function
 */
#include <random>
#include "min_cut_graph.hpp"
#include <stdio.h> 
#include <stdlib.h>
#include <vector>
#include <structures/union_find.hpp>
#include "../contracting_graph.hpp"
#include <random>
#include <iostream>
#define debug

namespace vg {
    namespace algorithms {

        using namespace std;

        pair<vector<vector<size_t>>, size_t> kargers_min_cut(Graph graph, const int n_iterations, const int seed, size_t V) {
 
     
            minstd_rand0 random_engine(seed);
            ContractingGraph cg(graph, V); 

            unordered_map<size_t, size_t> cgraph_total_edge_weights;
            pair<vector<vector<size_t>>, size_t> to_return;
            
            //check for empty or unconnected graphs
            if (V == 0 || V==1){
                //there is no min-cut to calculate
                return to_return;
                
            }
            //handle graphs with <= 2 nodes, with and without edges
            if (V==2){
                //check that both nodes have edges
                if (graph.nodes[0].edges.empty() == false && graph.nodes[1].edges.empty() == false){
                    vector<vector<size_t>> disjoint_sets;
                    //disjoint sets will just be two sets, each containing one node
                    //using index starting at 0 for nodes
                    vector<size_t> supernode0 = {0};
                    vector<size_t> supernode1 = {1};
                    disjoint_sets.push_back(supernode0);
                    disjoint_sets.push_back(supernode1);

                    
                    //assumes weights from node 1->node 2 and node 2->node 1 are equal
                    size_t weight_of_cut = graph.nodes[0].edges[0].weight;
                    to_return = make_pair(disjoint_sets, weight_of_cut);
                }else{
                    //not  connected graph
                    //to_return is empty
                    return to_return;
                    
                }

            }

            // get nodes will return the heads of all the nodes
            // at first call all nodes will be heads
            vector<size_t> super_nodes = cg.get_nodes();
            //get the total edge weights for super_nodes in contracted graph
            unordered_map<size_t, unordered_map<size_t, size_t>> contracted_graph;
            for (int i =0; i < super_nodes.size(); i++){
                    //get total edge weights for each super node
                    unordered_map<size_t, size_t> supernode_edge_weights = cg.get_edges(super_nodes[i]);
                    contracted_graph[super_nodes[i]] = supernode_edge_weights; 
#ifdef debug
            cout << "============================================================================= " << endl;    
            cout << "supernode edge weights for node " << super_nodes[i]   << endl;
            
            
                
#endif
            // tally up the total weights of incident edges
            int total_weight = 0;
            for(pair<size_t, size_t > element: supernode_edge_weights){
                    total_weight += element.second;
#ifdef debug
                    cout << element.first   << ":" <<element.second <<endl;
                    cout << "============================================================================= " << endl;
#endif
            }
                
                    //add total edge weight for each super node
                    //ex: 2: {1:1, 4:14}, sum up 1  4 and put that in the total edge weight list for node #2
                    cgraph_total_edge_weights[super_nodes[i]] = total_weight;
 #ifdef debug
            cout << "============================================================================= " << endl;    
            cout << "supernode edge weights total " << total_weight <<  " for node " << super_nodes[i]   << endl;
            
            cout << "============================================================================= " << endl;
                
#endif                   

            }

#ifdef debug
            cout << "============================================================================= " << endl;    
            cout << "number of nodes: " << V << endl;
            for (int i =0; i <V; i++){
                cout << "node "<< super_nodes[i]   << endl;
                
            }
            cout << "============================================================================= " << endl;
                
#endif 


            
            //for graphs with >2 nodes
            //assumes the graph is connected, acyclic, graph 
            while(V > 2){
                
                // create a vector with the node weights 
                vector<int> node_weights;
                size_t node_num;
                for (int i =0; i <V; i++){
                    node_num = super_nodes[i];
                    node_weights.push_back(cgraph_total_edge_weights[node_num]);
                }

                // create a discrete distrbution with node weights 
                discrete_distribution<int> nodes_distribution(begin(node_weights), end(node_weights));
                
                //pick an node proportional to its total weight from incident edges
                //and chooses uniformly from duplicates
                //discrete distribution returns an index
                int random_weight_idx = nodes_distribution(random_engine);
                size_t random_node = super_nodes[random_weight_idx];


                //get the edge weights of random node
                vector <size_t> rand_ew; 
                unordered_map<size_t, size_t> rand_node_edges = contracted_graph[random_node];
                for_each(rand_node_edges.begin(), rand_node_edges.end() , [&](pair<size_t, size_t > element){
                        //push back the weights 
                         rand_ew.push_back(element.second);
                });


                // create a discrete distrbution with edge weights 
                discrete_distribution<int> edges_distribution(begin(rand_ew), end(rand_ew));

                //pick a random edge weight proportional to its value
                int other_node;
                int random_edge_weight_idx = edges_distribution(random_engine);

                int count = 0;
                for (pair<size_t, size_t> element : rand_node_edges){
                    //iterate through the unordered_map rand_node:{other_node: edge_weight} 
                    //use the idx to access other node
                    if (count==random_edge_weight_idx){
                        other_node = element.first;
                        break;
                    }
                    count++;
                }

                //contract edge between random node and other node 
                cg.contract(random_node, other_node);

                //get nodes after contraction 
                super_nodes = cg.get_nodes();

                //calculate new number of supernodes left in contracted graph
                V = super_nodes.size();

                // will hold most up-to-date contracted graphs weights
                //clear it to recalculate and update
                cgraph_total_edge_weights.clear();

                //update contracted graph and total edge weights
                for (int i =0; i < super_nodes.size(); i++){
                    //get total edge weights for each super node
                    unordered_map<size_t, size_t> supernode_edge_weights = cg.get_edges(super_nodes[i]);
                    contracted_graph[super_nodes[i]] = supernode_edge_weights; 

                    // tally up the total weights of incident edges
                    int total_weight = 0;
                    for(pair<size_t, size_t > element: supernode_edge_weights){
                        total_weight += element.second;
                    }
                
                    //add total edge weight for each super node
                    //ex: 2: {1:1, 4:14}, sum up 1+14 and put that in the total edge weight list for node #2
                    cgraph_total_edge_weights[super_nodes[i]] = total_weight;

                }
#ifdef debug
                cout << "============================================================================= " << endl;
                cout << "random node " << random_node   << " and other node " << other_node   << " have been unioned" << endl;
                cout << "number of nodes after union: " << V << endl;
                for (int i =0; i <V; i++){
                    cout << "node "<< super_nodes[i]   << endl;
                    
                }
                cout << "============================================================================= " << endl;

                
#endif  
#ifdef debug
                cout << "============================================================================= " << endl;
                cout << "Printing contracted graph total edge weights " << endl;
                for (pair<size_t, size_t> element: cgraph_total_edge_weights){
                        cout << "node "<<element.first   << " : "<< element.second <<endl;
                }
                cout << "============================================================================= " << endl;
                
#endif              

                // check the new size of super nodes after contraction
                // if we have only two we can break out of the while loop
                if(V == 2){
                    vector<vector<size_t>> disjoint_sets = cg.get_disjoint_sets();
                    //compute the min cut of graph which is equal to the total edge weights of two supernodes 
                    size_t weight_of_cut;
                    for (pair<size_t, size_t> element: cgraph_total_edge_weights){
                        weight_of_cut = element.second;
                    }

                    to_return = make_pair(disjoint_sets, weight_of_cut);
                }
                
            }
            
            
        // or send back a pair containing contracted_graph, cgraph_total_edge_weights
        return to_return;    
        } 

        pair<vector<vector<size_t>>, size_t> compute_min_cut(Graph graph, const int n_iterations, const int seed, size_t V){

            // compute min-cut twice and choose the min-cut with least total graph weights
            //the minimum total edge weight of graph will give us the min-cut
            const int seed2 = 3;

            //TODO: generate seeds in here or send two seeds
            pair<vector<vector<size_t>>, size_t> to_return;
            pair<vector<vector<size_t>>, size_t> min_cut1 = kargers_min_cut(graph, n_iterations, seed, V);            
            pair<vector<vector<size_t>>, size_t> min_cut2 = kargers_min_cut(graph, n_iterations, seed2, V);
            if (min_cut1.second = 0 || min_cut2.second == 0 ){
                // if pair is empty pair.first and pair.second will both be initialized to 0 during contruction
                //return empty container
                return to_return;
            }
            if (min_cut1.second < min_cut2.second){
                to_return= min_cut1;

            }else{
                to_return= min_cut2;
            }
            cout << "min_cut " << to_return.second <<endl;
            return to_return;

        }
         
    }
}