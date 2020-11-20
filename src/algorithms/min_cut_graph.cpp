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
//#define debug

namespace vg {
    namespace algorithms {

        using namespace std;

        #ifdef debug
                cerr << "some debug statement" << endl; 
        #endif


        structures::UnionFind kargers_min_cut(Graph graph, const int n_iterations, const int seed, int V) {
 
     
            minstd_rand0 random_engine(seed);
            structures::UnionFind uf(V, true);
            
            

            ContractingGraph cg(graph,uf); 
            unordered_map<size_t, size_t> cgraph_total_edge_weights;
            
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
                    size_t total_weight = 0;
                    for_each(supernode_edge_weights.begin(), supernode_edge_weights.end() , [&](pair<size_t, size_t > element){
                        total_weight += element.second;
                    });
                    //add total edge weight for each super node
                    //ex: 2: {1:1, 4:14}, sum up 1+14 and put that in the total edge weight list for node #2
                    cgraph_total_edge_weights[super_nodes[i]] = total_weight;

                }

                // create a vector with the node weights 
                vector<int> node_weights;
                int node_num;
                for (int i =0; i <V; i++){
                    node_num = super_nodes[i];
                    node_weights.push_back(cgraph_total_edge_weights[node_num]);
                }
                // create a discrete distrbution with node weights 
                discrete_distribution<int> nodes_distribution(begin(node_weights), end(node_weights));
                
                //pick an node proportional to its total weight from incident edges
                //will choose the weight first and will get the node that maps to it after
                size_t incident_edge_totweight = nodes_distribution(random_engine);

                //check if there is more than one node with equal chosen weight
                vector <int> dup_nodes = find_dups(cgraph_total_edge_weights, incident_edge_totweight);
     
                int random_node;
                if(dup_nodes.size()>1){
                    //pick nodes with duplicate weights with equal likelihood 
                    random_node = sample_dups_uniformly(dup_nodes,random_engine);
                }else{
                    random_node = dup_nodes[0];
                }
                

                //get the edge weights of random node
                vector <int> rand_ew; 
                unordered_map<size_t, size_t> rand_node_edges = contracted_graph[random_node];
                for_each(rand_node_edges.begin(), rand_node_edges.end() , [&](pair<size_t, size_t > element){
                        //push back the weights 
                         rand_ew.push_back(element.second);
                });

                // create a discrete distrbution with edge weights 
                discrete_distribution<int> edges_distribution(begin(rand_ew), end(rand_ew));

                //pick a random edge weight proportional to its value
                size_t random_edge_weight = edges_distribution(random_engine);

                vector<int> dup_edges = find_dups(rand_node_edges, random_edge_weight);

                //if there are edges with duplicate total weights, sample from them uniformly
                int other_node;
                if(dup_edges.size()>1){
                    other_node = sample_dups_uniformly(dup_edges,random_engine);
                }else{
                    other_node = dup_edges[0];
                }

                //contract edge between random node and other node 
                uf.union_groups(random_node, other_node);

                
                // check the new size of super nodes after contraction
                // if we have only two we can break out of the while loop
                vector<size_t> super_nodes = cg.get_nodes();
                if(super_nodes.size() == 2){
                    break;
                }

                // will hold most up-to-date contracted graphs weights
                cgraph_total_edge_weights.clear();

            }
        // or send back a pair containing min_cut_of_cgraph and uf objects so we can get the disjoint sets
        // can also send back the disjoint set instead of the object
        return uf;    
        } 

        vector<vector<int>> compute_min_cut(Graph graph, const int n_iterations, const int seed, int V){

            // compute min-cut twice and choose the min-cut with least total graph weights
            //the minimum total edge weight of graph will give us the min-cut
            structures::UnionFind uf1 = kargers_min_cut(graph, n_iterations, seed, V);
            structures::UnionFind uf1 = kargers_min_cut(graph, n_iterations, seed, V);

            //two dijoint sets are nodes contracted in each of supernodes
            //we will use these sets to get out of local optima
            //might need to send back a pair <int graph_min_cut, vector<vector<int>> disjoint_sets> for unit testing





        }

        size_t sample_dups_uniformly(vector<int> to_check, minstd_rand0 random_engine){
            int size = to_check.size();
            uniform_int_distribution<int> distribution(0,size-1);
            int node_index = distribution(random_engine);

            size_t random_node = to_check[node_index];

            return random_node;

        }

        vector<int> find_dups(unordered_map<size_t, size_t> to_search, size_t num_to_find){
            vector<int> to_return;
            for_each(to_search.begin(), to_search.end() , [&](pair<size_t, size_t > element){
                    if(element.second == num_to_find){
                        to_return.push_back(element.first);
                    }
                    
            });
            return to_return;


        }
         
    }
}