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
            int V = graph.nodes.size();
            minstd_rand0 random_engine(seed);
            structures::UnionFind uf(V, true);
            

            // make adjacency list
            vector<int> adj[V];
            for(int i = 0; i< V; i++){
                for(int j = 0; j< graph.nodes.at(i).edges.size(); j++ ){
                    adj[i].push_back(graph.nodes.at(i).edges[j].other);

                }

            }

            ContractingGraph cg(graph, V, adj, uf); 
            // pick a random edge and contract it
            while(V > 2){
                //pick a random node 
                uniform_int_distribution<int> distribution(0, V-1);  
                size_t random_node = distribution(random_engine);
                
                //pick a random edge in that node
                int E = graph.nodes[random_node].edges.size();
                uniform_int_distribution<int> distribution(0, E-1); 
                size_t random_edge = distribution(random_engine);

                //get the other node of the edge
                size_t other_node = graph.nodes[random_node].edges[random_edge].other;

                //get the random_edge located between random_node and other_node 
                uf.union_groups(random_node, other_node);
                size_t group_id = uf.find_group(random_node);
                unordered_map<size_t, size_t> group_edges = cg.get_edges(group_id);
                

            }

            
            

        return 0;    
        } 
       
    

        
        
        
    
    }
}