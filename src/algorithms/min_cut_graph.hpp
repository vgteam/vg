#ifndef VG_ALGORITHMS_FIND_MIN_CUT_IN_GRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_FIND_MIN_CUT_IN_GRAPH_HPP_INCLUDED

/**
 * \file min_cut_graph.hpp
 * A randomized, probabilistic algorithm that 
 * finds the min cut on an undirected, weighted graph
 * that holds the indices of snarls
 */


#include <vector>
#include <utility>
#include <unordered_map>
#include <map>
#include <unordered_set>


namespace vg {
    namespace algorithms {

    using namespace std;
       
        struct Edge{
            int other; //node at other end 
            int weight;
        };

        struct Node{
            int weight;
            vector<Edge> edges;
        };
        struct Graph {
            unordered_map<size_t,Node> nodes;

        }; 

        pair<vector<unordered_set<size_t>>, size_t> kargers_min_cut(Graph graph, const int n_iterations, const int seed);

        pair<vector<unordered_set<size_t>>, size_t> compute_min_cut(Graph graph, const int n_iterations, const int seed);

        // Assumption: handles one connected component at a time
        vector<unordered_set<size_t>> min_cut_decomposition(Graph graph, const int n_iterations, const int seed);
        
        
    }
}

#endif