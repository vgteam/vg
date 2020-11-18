/// \file min_cut_graph.cpp
///  
/// unit tests for min_cut_graph construction and utility functions
///
#include <xg.hpp>
#include <stdio.h>
#include <iostream>
#include "../algorithms/min_cut_graph.hpp"
#include <vector>
#include "catch.hpp"

// #define debug

namespace vg {
    namespace unittest {
        
        using namespace std;
        using vg::algorithms::Graph;
        using vg::algorithms::Edge;
        using vg::algorithms::Node;
        using vg::algorithms::compute_min_cut;

        const int seed = 0;
        const int n_iterations = 100;

        TEST_CASE("Min_cut1") {
            
            SECTION("Test1: Can find a min-cut on a 4 node graph") {
                /* Let us create following undirected, weighted graph  
                      e1 10  
                    1--------2  
                    | \      |  
                e2 6|   |5 e5 |15  
                    |      \ |  e4
                    3--------4  
                       e3 4 
                */
                Graph graph; 
                Edge edge12, edge21, edge13, edge31, edge14, edge41, edge34, edge43,edge24, edge42; //naming convention edge:source:destination
                Node node1, node2, node3, node4;
                int V = 4; //size of nodes
                //weights
                edge12.weight = 10;
                edge21.weight = 10;
                edge13.weight = 6;
                edge31.weight = 6;
                edge14.weight = 5;
                edge41.weight = 5;
                edge34.weight = 4;
                edge43.weight = 4;
                edge24.weight = 15;
                edge42.weight = 15;

                //other node
                edge12.other = 2;
                edge21.other = 1;
                edge13.other = 3;
                edge31.other = 1;
                edge14.other = 4;
                edge41.other = 1;
                edge34.other = 4;
                edge43.other = 3;
                edge24.other = 4;
                edge42.other = 2;

                //node 1
                node1.edges.push_back(edge12);
                node1.edges.push_back(edge13);
                node1.edges.push_back(edge14);
                node1.weight = 21;


                //node 2
                node2.edges.push_back(edge21);
                node2.edges.push_back(edge24);
                node2.weight = 25;

                //node 3
                node3.edges.push_back(edge31);
                node3.edges.push_back(edge34);
                node3.weight = 10;

                //node 4
                node4.edges.push_back(edge41);
                node4.edges.push_back(edge42);
                node4.edges.push_back(edge43);
                node4.weight = 24;

                graph.nodes.push_back(node1);
                graph.nodes.push_back(node2);
                graph.nodes.push_back(node3);
                graph.nodes.push_back(node4);
                
            
                //Karger's min-cut
                vector<vector<int>> disjoint_sets= compute_min_cut(graph, n_iterations, seed, V);
                
                //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
                // REQUIRE(mincut == 2);

/*Notes
for unit test we assume weights have been assigned using match(graph) that assigns weights a score using its consistency with the haplotype phasing


*/


                




  
#ifdef debug
                //cerr << "this is the node address"<< &n1 << endl;
#endif
                
                  
            }
 
        }

    }

}


