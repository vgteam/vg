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


namespace vg {
namespace algorithms {

using namespace std;

class find_min_cut{ 
    vector<Node> graph;
    const int n_iterations;
    const int seed;
    mutable minstd_rand0 random_engine;

public:
    
    find_min_cut(vector<Node> graph, const int n_iterations, const int seed); 

    int kargers_min_cut(vector<Node> graph)const;    

    

};

struct Edge{
    int other;
    int weight;
};

struct Node{
    int weight;
    vector<Edge> edges;
};
using Graph = vector<Node>;
}
}

#endif