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

#define debug

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
                size_t V = 4; //size of nodes
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
                //edge12.other is edge1->2 but we use the index starting from 0
                //so, edge12 = edge1->edge1
                //TODO: fix the indices
                edge12.other = 1;
                edge21.other = 0;
                edge13.other = 2;
                edge31.other = 0;
                edge14.other = 3;
                edge41.other = 0;
                edge34.other = 3;
                edge43.other = 2;
                edge24.other = 3;
                edge42.other = 1;

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
                pair<vector<vector<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
                vector<vector<size_t>> disjoint_sets = to_recv.first;
                size_t mincut = to_recv.second;
                
                //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
                // REQUIRE(mincut == 2);

  
#ifdef debug
                //cerr << "this is the node address"<< &n1 << endl;
#endif
                
                  
            }
 
        }

            TEST_CASE("Min_cut2") {
            
            SECTION("Test2: Can find a min-cut on a 9 node graph") {
       
                Graph graph; 
                Edge edge12, edge13, edge21, edge23, edge31, edge32, edge34, edge35, edge43, edge53, edge56, edge58, edge65, edge67, edge69, edge76, edge78, edge85, edge87, edge89, edge96, edge98; //naming convention edge:source:destination
                Node node1, node2, node3, node4, node5, node6, node7, node8, node9;
                size_t V = 9; //size of nodes
                
                //weights
                edge12.weight = 1;
                edge13.weight = 2;
                edge21.weight = 1;
                edge23.weight = 4;
                edge31.weight = 2;
                edge32.weight = 4;
                edge34.weight = 3;
                edge35.weight = 5;
                edge43.weight = 3;
                edge53.weight = 5;
                edge56.weight = 4;
                edge58.weight = 1;
                edge65.weight = 4;
                edge67.weight = 7;
                edge69.weight = 5;
                edge76.weight = 7;
                edge78.weight = 2;
                edge85.weight = 1;
                edge87.weight = 2;
                edge89.weight = 1;
                edge96.weight = 5;
                edge98.weight= 1;

                //other 
                //edge12.other is edge1->2 but we use the index starting from 0
                edge12.other = 1;
                edge13.other = 2;
                edge21.other = 0;
                edge23.other = 2;
                edge31.other = 0;
                edge32.other = 1;
                edge34.other = 3;
                edge35.other = 4;
                edge43.other = 2;
                edge53.other = 2;
                edge56.other = 5;
                edge58.other = 7;
                edge65.other = 4;
                edge67.other = 6;
                edge69.other = 8;
                edge76.other = 5;
                edge78.other = 7;
                edge85.other = 4;
                edge87.other = 6;
                edge89.other = 8;
                edge96.other = 5;
                edge98.other = 7;


                //node 1
                node1.edges.push_back(edge12);
                node1.edges.push_back(edge13);
                node1.weight = 3;


                //node 2
                node2.edges.push_back(edge21);
                node2.edges.push_back(edge23);
                node2.weight = 5;

                //node 3
                node3.edges.push_back(edge31);
                node3.edges.push_back(edge32);
                node3.edges.push_back(edge34);
                node3.edges.push_back(edge35);
                node3.weight = 14;

                //node 4
                node4.edges.push_back(edge43);
                node4.weight = 3;

                //node5
                node5.edges.push_back(edge53);
                node5.edges.push_back(edge56);
                node5.edges.push_back(edge58);
                node5.weight = 10;

                //node6
                node6.edges.push_back(edge65);
                node6.edges.push_back(edge67);
                node6.edges.push_back(edge69);
                node6.weight = 16;
                
                //node7
                node7.edges.push_back(edge76);
                node7.edges.push_back(edge78);
                node7.weight = 9;
                
                //node8
                node8.edges.push_back(edge85);
                node8.edges.push_back(edge87);
                node8.edges.push_back(edge89);
                node8.weight = 4;

                //node9
                node9.edges.push_back(edge96);
                node9.edges.push_back(edge98);
                node9.weight = 6;
                
                graph.nodes.push_back(node1);
                graph.nodes.push_back(node2);
                graph.nodes.push_back(node3);
                graph.nodes.push_back(node4);
                graph.nodes.push_back(node5);
                graph.nodes.push_back(node6);
                graph.nodes.push_back(node7);
                graph.nodes.push_back(node8);
                graph.nodes.push_back(node9);
            
                //Karger's min-cut
                pair<vector<vector<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
                vector<vector<size_t>> disjoint_sets = to_recv.first;
                size_t mincut = to_recv.second;
                
                //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
                REQUIRE(mincut == 1);

            }
 
        }
        TEST_CASE("Min_cut3") {
            
            SECTION("Test3: Can find a min-cut on a 2 node graph") {
 
                Graph graph; 
                Edge edge12, edge21; //naming convention edge:source:destination
                Node node1, node2;
                size_t V = 2; //size of nodes

                //weights
                edge12.weight = 10;
                edge21.weight = 10;

                //other node
                //edge12.other is edge1->2 but we use the index starting from 0
                //so, edge12 = edge1->edge1
                //TODO: fix the indices
                edge12.other = 1;
                edge21.other = 0;

                //node 1
                node1.edges.push_back(edge12);
                node1.weight = 10;


                //node 2
                node2.edges.push_back(edge21);
                node2.weight = 10;


                graph.nodes.push_back(node1);
                graph.nodes.push_back(node2);

                
                //Karger's min-cut
                pair<vector<vector<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
                vector<vector<size_t>> disjoint_sets = to_recv.first;
                size_t mincut = to_recv.second;
                
                //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
                REQUIRE(mincut == 10);
            }
 
        }
        TEST_CASE("Min_cut4") {
                
            SECTION("Test4: Can find a min-cut on a 0 node graph, 0 edges") {

                Graph graph; 
                size_t V = 0; //size of nodes

                //Karger's min-cut
                pair<vector<vector<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
                vector<vector<size_t>> disjoint_sets = to_recv.first;
                size_t mincut = to_recv.second;
                
                //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
                REQUIRE(mincut == 0);
            }
    
        }
        TEST_CASE("Min_cut5") {
                
            SECTION("Test5: Can find a min-cut on a 2 node graph, 0 edges") {

                Graph graph; 
                Node node1, node2;
                size_t V = 2; //size of nodes

                //node 1
                node1.weight = 0;


                //node 2
                node2.weight = 0;


                graph.nodes.push_back(node1);
                graph.nodes.push_back(node2);

                //Karger's min-cut
                pair<vector<vector<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
                vector<vector<size_t>> disjoint_sets = to_recv.first;
                size_t mincut = to_recv.second;
                
                //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
                REQUIRE(mincut == 0);
            }
    
        }
        TEST_CASE("Min_cut6") {
                
            SECTION("Test6: Can find a min-cut on a 4 node graph, 0 edges") {

                Graph graph; 
                Node node1, node2, node3, node4;
                size_t V = 4; //size of nodes

                //node 1
                node1.weight = 0;

                //node 2
                node2.weight = 0;

                //node 3
                node1.weight = 0;

                //node 4
                node2.weight = 0;


                graph.nodes.push_back(node1);
                graph.nodes.push_back(node2);
                graph.nodes.push_back(node3);
                graph.nodes.push_back(node4);


                //Karger's min-cut
                pair<vector<vector<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
                vector<vector<size_t>> disjoint_sets = to_recv.first;
                size_t mincut = to_recv.second;
#ifdef debug
                cout << "mincut"<< mincut << endl;
#endif                
                //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
                REQUIRE(mincut == 0);
            }
    
        }

    }

}


