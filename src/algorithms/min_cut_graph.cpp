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
// #define debug

namespace vg {
    namespace algorithms {

        using namespace std;




        pair<vector<vector<size_t>>, size_t> kargers_min_cut(Graph graph, const int n_iterations, const int seed, int V) {
 
     
            minstd_rand0 random_engine(seed);
            structures::UnionFind uf(V, true);
            
            

            ContractingGraph cg(graph, uf); 
            unordered_map<size_t, size_t> cgraph_total_edge_weights;
            pair<vector<vector<size_t>>, size_t> to_return;
            
            //TODO: add case for V <= 2 from beginning 
            while(V > 2){
                

                // get nodes will return the heads of all the nodes
                // at first call all nodes will be heads
                vector<size_t> super_nodes = cg.get_nodes();
                V = super_nodes.size();


                //get the total edge weights for super_nodes in contracted graph
                unordered_map<size_t, unordered_map<size_t, size_t>> contracted_graph;
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
#ifdef debug
                cout << "random node " << random_node << endl;
                
#endif

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
#ifdef debug
                cout << "random edge idx " << random_edge_weight_idx << endl;
                cout << "other node  " << other_node << endl;
                size_t g_size_random =  uf.group_size(random_node);
                size_t g_size_other =  uf.group_size(other_node);
                cout << "group size random" << g_size_random <<endl;
                cout << "group size other" << g_size_other <<endl;
                
#endif
    
                //contract edge between random node and other node 
                uf.union_groups(random_node, other_node);
                
                
                
                // get the size of the current groups - will tell us how many supernodes
                vector<vector<size_t>> all_gs = uf.all_groups();
                V = all_gs.size();
#ifdef debug

               
                // cout << "g size" << all_gs.size() << endl;
                // for (size_t i = 0; i < all_gs.size(); i++){
                //     cout << "g_i " << i << endl;
                //     vector<size_t> v = all_gs[i];

                //     for(size_t j = 0; j < v.size(); j++)
                //         cout << "m_j " << v[j] <<endl;
                // }

                // g_size_random =  uf.group_size(random_node);
                // g_size_other =  uf.group_size(other_node);
                // cout << "group size random" << g_size_random <<endl;
                // cout << "group size other" << g_size_other <<endl;
                // vector<size_t> group = uf.group(random_node);
                // for (size_t i = 0; i < group.size(); i++){

                //     cout << "group i member "<< group[i] <<endl;
                // }
                //call get_nodes again to see if the contraction worked
                vector<size_t> sp = cg.get_nodes();
                
                
#endif
                // check the new size of super nodes after contraction
                // if we have only two we can break out of the while loop
                if(V == 2){
                    //compute the min cut of graph which is equal to the total edge weights of two supernodes 
                    size_t min_cut;
                    for (pair<size_t, size_t> element: cgraph_total_edge_weights){
                        min_cut += element.second;
                    }

                    vector<vector<size_t>> disjoint_sets;
                    for(size_t i = 0; i<V; i++){
                        vector<size_t> d_set;
                        for(pair<size_t, size_t> element: contracted_graph[i]){
                            d_set.push_back(element.first);
                        }
                        disjoint_sets.push_back(d_set);

                    }
                    to_return = make_pair(disjoint_sets, min_cut);
                    break;
                }
                
                // will hold most up-to-date contracted graphs weights
                cgraph_total_edge_weights.clear();
#ifdef debug
                cout << "breaking after one loop" << endl;
                break;
#endif
            }
            
            
        // or send back a pair containing contracted_graph, cgraph_total_edge_weights
        return to_return;    
        } 

        pair<vector<vector<size_t>>, size_t> compute_min_cut(Graph graph, const int n_iterations, const int seed, int V){

            // compute min-cut twice and choose the min-cut with least total graph weights
            //the minimum total edge weight of graph will give us the min-cut
            const int seed2 = 3;

            //TODO: generate seeds in here or send two seeds
            
            pair<vector<vector<size_t>>, size_t> to_return;
            pair<vector<vector<size_t>>, size_t> min_cut1 = kargers_min_cut(graph, n_iterations, seed, V);
            pair<vector<vector<size_t>>, size_t> min_cut2 = kargers_min_cut(graph, n_iterations, seed2, V);

            if (min_cut1.second < min_cut2.second){
                to_return = min_cut1;

            }else{
                to_return = min_cut2;

            }

            return to_return;

        }
         
    }
}