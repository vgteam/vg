/**
 * \file min_cut_graph.cpp
 *
 * Contains implementation of min_cut_graph function
 */

#include "min_cut_graph.hpp"
#include <stdio.h> 
#include <stdlib.h>
#include <vector>
#include <structures/union_find.hpp>
#include "../contracting_graph.hpp"
#include <random>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <functional>

// #define debug_min_decomp
// #define debug_kargers_min
// #define debug_compute_min

namespace vg {
    namespace algorithms {

        using namespace std;

        pair<vector<unordered_set<size_t>>, size_t> kargers_min_cut(Graph graph, const int seed) {
 
            size_t V = graph.get_size();
            minstd_rand0 random_engine(seed);
            ContractingGraph cg(graph); 
            unordered_map<size_t, size_t> cgraph_total_edge_weights;
            pair<vector<unordered_set<size_t>>, size_t> to_return;

            vector<size_t> node_ids = graph.get_node_ids();
            // check for graph containing a node without an edge
            for (auto& id : node_ids){
                if(graph.get_node_by_id(id).edges.size() <=0){
                    //return empty container
#ifdef debug_kargers_min
                    cout << "============================================================================= " << endl;    
                    cout << "Disconnected graph " << endl; 
                    cout << "Node " <<id << "has 0 edges"<<endl;         
#endif
                    return to_return;
                }
            }

            //check for empty graph or 1-node graph
            if (V == 0 || V==1){
#ifdef debug_kargers_min
                cout << "============================================================================= " << endl;    
                cout << "Empty or 1-node graph "  << endl;         
#endif
                //there is no min-cut to calculate
                //return empty container
                return to_return;
                
            }
            //graph with exactly 2-nodes
            if (V==2){
#ifdef debug_kargers_min
                cout << "============================================================================= " << endl;
                cout << "2-node graph "  << endl; 
#endif
                //both nodes will have edges since we tested for unconnected graph above
                vector<unordered_set<size_t>> disjoint_sets;
                //disjoint sets will just be two sets, each containing one node
                //using index starting at 0 for nodes
                
                
                vector<size_t> vnodes = graph.get_node_ids();
                for(auto& id : vnodes){
#ifdef debug_kargers_min            
                   
                    // node id
                    cout << "node "<<id  << endl; 
#endif                 
                }    
                
                
                unordered_set<size_t> supernode0 = {vnodes[1]};
                unordered_set<size_t> supernode1 = {vnodes[0]};
                disjoint_sets.push_back(supernode0);
                disjoint_sets.push_back(supernode1);

                //assumes weights from node 0->node 1 and node 1->node 0 are equal
                size_t node_id = vnodes[0];
                size_t weight_of_cut = graph.get_node_by_id(node_id).edges[0].weight;
                to_return = make_pair(disjoint_sets, weight_of_cut);
                
            }

            // get nodes will return the heads of all the nodes
            // at first call all nodes will be heads
            vector<size_t> super_nodes = cg.get_nodes();
            //get the total edge weights for super_nodes in contracted graph
           
            for (int i =0; i < super_nodes.size(); i++){
                    //get total edge weights for each super node
                    unordered_map<size_t, size_t> supernode_edge_weights = cg.get_edges(super_nodes[i]);
                
#ifdef debug_kargers_min
                    cout << "============================================================================= " << endl;    
                    cout << "supernode edge weights for node " << super_nodes[i]   << endl;         
#endif
            // tally up the total weights of incident edges
            int total_weight = 0;
            for(pair<size_t, size_t > element: supernode_edge_weights){
                    total_weight += element.second;
#ifdef debug_kargers_min
                    cout << element.first   << ":" <<element.second <<endl;
#endif
            }
                
                    //add total edge weight for each super node
                    //ex: 2: {1:1, 4:14}, sum up 1  4 and put that in the total edge weight list for node #2
                    cgraph_total_edge_weights[super_nodes[i]] = total_weight;
 #ifdef debug_kargers_min
            cout << "============================================================================= " << endl;    
            cout << "supernode edge weights total " << total_weight <<  " for node " << super_nodes[i]   << endl;
            
            cout << "============================================================================= " << endl;
                
#endif                   

            }

#ifdef debug_kargers_min
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
                unordered_map<size_t, size_t> rand_node_edges = cg.get_edges(random_node);
                
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

                    // tally up the total weights of incident edges
                    int total_weight = 0;
                    for(pair<size_t, size_t > element: supernode_edge_weights){
                        total_weight += element.second;
                    }
                
                    //add total edge weight for each super node
                    //ex: 2: {1:1, 4:14}, sum up 1+14 and put that in the total edge weight list for node #2
                    cgraph_total_edge_weights[super_nodes[i]] = total_weight;

                }
#ifdef debug_kargers_min
                cout << "============================================================================= " << endl;
                cout << "random node " << random_node   << " and other node " << other_node   << " have been unioned" << endl;
                cout << "number of nodes after union: " << V << endl;
                for (int i =0; i <V; i++){
                    cout << "node "<< super_nodes[i]   << endl;
                    
                }
                cout << "============================================================================= " << endl;

                
#endif  
#ifdef debug_kargers_min
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
                    vector<vector<size_t>> disjoint_vector = cg.get_disjoint_sets();
 #ifdef debug_kargers_min   
                cout << "============================================================================= " << endl;
                cout << "vector"<<endl;
                for(int i = 0; i < 2; i++){

                    for ( int j = 0; j< disjoint_vector[i].size(); j++ ){
                            cout << "vector " <<i <<"element " <<  disjoint_vector[i][j] << endl;

                    }
                }
                cout << "============================================================================= " << endl;
                    
#endif
                    // make vector into unordered_set
                    unordered_set<size_t> disjoint_set1(disjoint_vector[0].begin(), disjoint_vector[0].end());
                    unordered_set<size_t> disjoint_set2(disjoint_vector[1].begin(), disjoint_vector[1].end());                                     
                    vector<unordered_set<size_t>> disjoint_sets;
                    disjoint_sets.push_back(disjoint_set1);
                    disjoint_sets.push_back(disjoint_set2);                   

                    //compute the min cut of graph which is equal to the total edge weights of two supernodes 
                    size_t weight_of_cut;
                    for (pair<size_t, size_t> element: cgraph_total_edge_weights){
                        weight_of_cut = element.second;
                    }
#ifdef debug_kargers_min
                cout << "============================================================================= " << endl;
                cout << "Weight of cut " << weight_of_cut<< endl;
                cout << "============================================================================= " << endl;
                
#endif
                    to_return = make_pair(disjoint_sets, weight_of_cut);

                }
                
            }
            
            
        // or send back a pair containing min_cut, disjoint_sets
        return to_return;    
        } 

        pair<vector<unordered_set<size_t>>, size_t> compute_min_cut(Graph graph, const int seed){

            // compute min-cut twice and choose the min-cut with least total graph weights
            //the minimum total edge weight of graph will give us the min-cut
            const int seed2 = seed+1;

            //TODO: generate seeds in here or send two seeds
            pair<vector<unordered_set<size_t>>, size_t> to_return;
            pair<vector<unordered_set<size_t>>, size_t> min_cut1 = kargers_min_cut(graph, seed);          
            pair<vector<unordered_set<size_t>>, size_t> min_cut2 = kargers_min_cut(graph, seed2);

            if (min_cut1.second == 0 || min_cut2.second == 0 ){
                // if pair is empty pair.first and pair.second will both be initialized to 0 during contruction
                //return empty container
#ifdef debug_compute_min
                cout << "============================================================================= " << endl;
                cout << "RETURNING EMPTY MINCUT"  <<endl;
                cout << "============================================================================= " << endl;
                
#endif
                return to_return;
            }
#ifdef debug_compute_min   
            cout << "============================================================================= " << endl;
            cout << "MIN-CUT-GRAPH: set"<<endl;
            cout << "disjoint sets in min cut 1 " << endl;
            vector<unordered_set<size_t>> disjoint_set = min_cut1.first;
            for (size_t i = 0;  i < disjoint_set.size(); i++){
                for (auto& x:disjoint_set[i] ) {
                cout << "MCG set "<< i << "has" << x <<endl;
                }

            }
            cout << "cut is : " << min_cut1.second <<endl;
            cout << "============================================================================= " << endl;
            cout << "disjoint sets in min cut 2 " << endl;
            vector<unordered_set<size_t>> disjoint_set2 = min_cut2.first;
            for (size_t i = 0;  i < disjoint_set2.size(); i++){
                for (auto& x:disjoint_set2[i] ) {
                cout << "MCG set "<< i << "has" << x <<endl;
                }

            }
            cout << "cut is : " << min_cut2.second <<endl;
            cout << "============================================================================= " << endl;


#endif
            if (min_cut1.second < min_cut2.second){
                to_return= min_cut1;

            }else{
                to_return= min_cut2;
            }
            
            return to_return;

        }
        vector<unordered_set<size_t>> min_cut_decomposition(Graph graph, const int seed){

            vector<unordered_set<size_t>> Gamma;

            const int rand_seed = seed;

            function<void(Graph)> recurse = [&](Graph graph){
#ifdef debug_min_decomp
                cout << "============================================================================= " << endl;
                cout << "MIN-CUT-DECOMPOSITION"  <<endl;
                cout << "============================================================================= " << endl;
#endif               
                size_t V = graph.get_size();

                //base case
                //singleton nodes won't be added to Gamma
                if(V <= 2){
#ifdef debug_min_decomp             
                    cout << "Base Case"  <<endl;
                    cout << "V= " << V <<endl;
           
#endif
                    return;
                }
                pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, rand_seed);
                vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
               
                if(disjoint_sets.empty()==1){
#ifdef debug_min_decomp             
                    cout << "Empty sets"  <<endl;
#endif
                    return;
                }
#ifdef debug_min_decomp             
                cout << "set size "<<disjoint_sets[0].size() << ", " << disjoint_sets[1].size() <<endl;
#endif
                //push sets into Gamma at each iteration, only if size >=2
                if(disjoint_sets[0].size() >=2){
                    Gamma.push_back(disjoint_sets[0]);
                }
                if(disjoint_sets[1].size() >=2){
                    Gamma.push_back(disjoint_sets[1]);
                }
                
                

                //build the subgraphs from disjoint sets
                vector<Graph> subgraph(2);
                vector<size_t> node_ids = graph.get_node_ids();
                for(size_t h =0; h < disjoint_sets.size(); h++){

                    for(auto& id : node_ids){
                    // if node from original graph is in the disjoint set
                        if (disjoint_sets[h].count(id)==1){
                     
                            size_t node_weight =0;
                            Node node;
                            Edge edge;
                            // check if any edges connect to other nodes in the disjoint set
                            for(size_t j =0; j < graph.get_node_by_id(id).edges.size(); j++){
                                
                                if (disjoint_sets[h].count(graph.get_node_by_id(id).edges[j].other)==1){
                                    edge.other = graph.get_node_by_id(id).edges[j].other;
                                    edge.weight = graph.get_node_by_id(id).edges[j].weight;
                                    node_weight += edge.weight;
                                    node.edges.push_back(edge);                                   


                                }
                            }
                            //tally the node weight using the edges
                            node.weight = node_weight; 
                            subgraph[h].add_node(id, node);

                        
                        }
                    
                    }

                }
                
                recurse(subgraph[0]);

                recurse(subgraph[1]);


            };
            
            recurse(graph);
#ifdef debug_min_decomp
            for(size_t i = 0; i < Gamma.size(); i++){
                for (auto& x:Gamma[i] ) {
                    cout << "Gamma " << i <<"has " << x <<endl; 
                }

            }
            
#endif
            return Gamma;
        }

         
    }
}