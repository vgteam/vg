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
        using vg::algorithms::kargers_min_cut;
        const int seed = 0;
        const int n_iterations = 100;

        TEST_CASE("Test1") {
            
            SECTION("Test1: ") {
                /* Let us create following undirected, weighted graph  
                      e1 10  
                    1--------2  
                    | \      |  
                e2 6|   |5 e5 |15  
                    |      \ |  e4
                    3--------4  
                       e3 4 
                */
            //    struct Edge{ 
            //        int other; //destination node 
            //        //source node will be implied because edge will be contained as an edge from source node
            //        int weight;
            //     }; 
                
            //     struct Node { 
            //        int weight; //weight of incident edges
            //        vector<Edge> edges; 
            //     }; 
            //     using Graph = vector<Node>;

                Graph graph; 
                Edge edge12, edge21, edge13, edge31, edge14, edge41, edge34, edge43,edge24, edge42; //naming convention edge:source:destination
                Node node1, node2, node3, node4;

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

                graph.push_back(node1);
                graph.push_back(node2);
                graph.push_back(node3);
                graph.push_back(node4);

                //Karger's min-cut
                int mincut = kargers_min_cut(graph, n_iterations, seed);
                
                //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
                REQUIRE(mincut == 2);

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


