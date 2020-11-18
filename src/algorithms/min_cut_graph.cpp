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
                // will hold most up-to-date contracted graphs weights
                cgraph_total_edge_weights.clear();

                // get nodes will return the heads of all the nodes
                // at first call all nodes will be heads
                vector<size_t> super_nodes = cg.get_nodes();
                V = super_nodes.size();

                // create a vector with the node weights 
                vector<int> node_weights;
                int node_num;
                for (int i =0; i <V; i++){
                    node_num = super_nodes[i];
                    node_weights.push_back(graph.nodes[node_num].weight);
                }
                // create a uniform discrete distrbution with node weights 
                discrete_distribution<int> nodes_distribution(begin(node_weights), end(node_weights));
                //pick a random node proportional to its weight
                size_t random_node = nodes_distribution(random_engine);

                // get the group_id for the contracted node
                //initially this will be equal to the node 
                size_t group_id = uf.find_group(random_node);

                //create a vector with edges of random node
                vector<int> edge_weights;
                for (int j =0; j <graph.nodes[random_node].edges.size(); j++){
                    edge_weights.push_back(graph.nodes[random_node].edges[j].weight);
                }
                // create a uniform discrete distrbution with edge weights 
                discrete_distribution<int> edges_distribution(begin(edge_weights), end(edge_weights));
                //pick a random edge in that node
                size_t random_edge = edges_distribution(random_engine);

                //get the other node of the edge
                size_t other_node = graph.nodes[random_node].edges[random_edge].other;

                //contract edge between random node and other node 
                uf.union_groups(random_node, other_node);

                //get the total edge weights for super_nodes in contracted graph
                for (int i =0; i < super_nodes.size(); i++){
                    //get total edge weights for each super node
                    unordered_map<size_t, size_t> supernode_edge_weights = cg.get_edges(super_nodes[i]);
                    size_t total_weight = 0;
                    for_each(supernode_edge_weights.begin(), supernode_edge_weights.end() , [&](pair<size_t, size_t > element){
                        total_weight += element.second;
                    });
                    //add total edge weight for each super node
                    //ex: 2: {1:1, 4:14}, sum up 1+14 and put that in the total edge weight list for node #2
                    cgraph_total_edge_weights[super_nodes[i]] = total_weight;

                }
                // check the new size of super nodes after contraction
                // if we have only two we can break out of the while loop
                vector<size_t> super_nodes = cg.get_nodes();
                if(super_nodes.size() == 2){
                    break;
                }

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
         
    }
}