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
//#define debug

namespace vg {
    namespace algorithms {

        using namespace std;

        #ifdef debug
                cerr << "some debug statement" << endl; 
        #endif
        int kargers_min_cut(vector<Node> graph , const int n_iterations, const int seed) {

            
        } 
        // // a structure to represent a weighted edge in graph  
        // class Edge  
        // {  
        //     public: 
        //     int src, dest, weight;  
        // };  
        
        // // a structure to represent a connected, undirected  
        // // and weighted graph  
        // class Graph  
        // {  
        //     public: 
        //     // V-> Number of vertices, E-> Number of edges  
        //     int V, E;  
        
        //     // graph is represented as an array of edges.  
        //     // Since the graph is undirected, the edge  
        //     // from src to dest is also edge from dest  
        //     // to src. Both are counted as 1 edge here.  
        //     Edge* edge;  
        // };  

        // Graph* create_graph(int V, int E) {
        //     Graph* graph = new Graph; 
        //     graph->V = V; 
        //     graph->E = E; 
        //     graph->edge = new Edge[E]; 
        //     return graph;

        // }
      

        // // Compare two edges according to their weights.  
        // // Used in qsort() for sorting an array of edges  
        // int myComp(const void* a, const void* b)  
        // {  
        //     Edge* a1 = (Edge*)a;  
        //     Edge* b1 = (Edge*)b;  
        //     return a1->weight > b1->weight;  
        // }

        // contraction(Graph* graph, int vertices_length, minstd_rand0& random_engine ){
            
        //     //base case: 
        //     //if there are exactly two nodes left, stop
        //     //the edges crossing those nodes form a cut
        //     // return (C = (S,V-S)) equal to 1
            
        //     // if get_length(graph) < 6, use brute force for min-cut 
            
            
        //     while (get_length(graph) > 6 ){
        //         int start = discrete_uniform_sample(random_engine);
                
            
        //     //else pick a random edge, contract it, repeat 


        //     }

        // }
        // kerger_stein_min_cut(){
        //     //call constraction t log n times
        //     // return the cut that is the smallest among those constructed in above step

        // }
        // // choose a start and end edge randomly and uniformly from the graph
        // discrete_uniform_sample(minstd_rand0& random_engine){

        //     number_of_edges = get_num_edges();

        //     uniform_int_distribution<int> distribution(0, number_of_edges-1);  
        //     int random_num = distribution(random_engine);
        // }

        
        // get_num_edges(){}

        // get_length(Graph* graph){}

        // // Given an edfe (u,v) in a multigraph, we can contract u,v as follows:
        // // delete all edges between u,v
        // // replace u and v with a new "supernode" uv
        // // replace a;; edges incident to u or v with edges incident to the supernode 
        // contract_edges(){}
        
        
        
    
    }
}