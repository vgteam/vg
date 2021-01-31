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

            inline vector<size_t> get_node_ids(){
                vector<size_t> node_ids;
                for (auto& id_and_node : nodes){
                    size_t node_id = id_and_node.first; 
                    node_ids.push_back(node_id);
                }
                return node_ids;

            }

            inline size_t get_size(){
                return nodes.size();
            }

            inline Node get_node_by_id(size_t node_id){
                return nodes.at(node_id); 
            }

            inline void add_node(size_t id, Node node){
                nodes.emplace(id, node);
            }

        }; 

        pair<vector<unordered_set<size_t>>, size_t> kargers_min_cut(Graph graph, const int n_iterations, const int seed);

        pair<vector<unordered_set<size_t>>, size_t> compute_min_cut(Graph graph, const int n_iterations, const int seed);

        // Assumption: handles one connected component at a time
        vector<unordered_set<size_t>> min_cut_decomposition(Graph graph, const int n_iterations, const int seed);
        
        
    }
}

#endif