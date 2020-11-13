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


        int64_t kargers_min_cut(Graph graph, const int n_iterations, const int seed, int V) {
 
     
            minstd_rand0 random_engine(seed);
            structures::UnionFind uf(V, true);
            vector<unordered_map<size_t, size_t>> group_total_edge_weights;

            // // make adjacency list
            vector<vector<int> > adj;
            // for(int i = 0; i< V; i++){
            //     for(int j = 0; j< graph.nodes.at(i).edges.size(); j++ ){
            //         adj[i].push_back(graph.nodes.at(i).edges[j].other);

            //     }

            // }
            

            ContractingGraph cg(graph, V, adj, uf); 
            // pick a random edge and contract it
            while(V > 2){

                // create a vector with the node weights 
                vector<int> node_weights;
                for (int i =0; i <V; i++){
                    node_weights.push_back(graph.nodes[i].weight);
                }
                // create a uniform discrete distrbution with node weights 
                discrete_distribution<int> nodes_distribution(begin(node_weights), end(node_weights));
                //pick a random node proportional to its weight
                size_t random_node = nodes_distribution(random_engine);
                
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

                //subtract the contracted edge weight from the incident edges total weight 
                //we don't want to update the original graph
                // graph.nodes[random_node].weight = graph.nodes[random_node].weight - graph.nodes[random_node].edges[random_edge].weight;
                // graph.nodes[other_node].weight = graph.nodes[other_node].weight - graph.nodes[other_node].edges[random_edge].weight;
                
                // have to check how many nodes we have and will break out of while loop when we only have two nodes in contracted graph
                vector<size_t> super_nodes = cg.get_nodes();
                V = super_nodes.size();

                // get the group_id for the contracted node 
                size_t group_id = uf.find_group(random_node);

                //get the total edge weights for contracted graph
                //how do we know the group_id once we finish the loop?
                //maybe we don't want to keep a vector outside and only keep the last result?
                group_total_edge_weights.push_back(cg.get_edges(group_id));


                //the minimum total edge weight of graph will give us the min-cut

                //two dijoint sets are the nodes inside each super-node or head

                

            }

            
            

        return 0;    
        } 
         
    }
}