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
#include <unordered_set>

#define debugs

namespace vg {
    namespace unittest {
        
        using namespace std;
        using vg::algorithms::Graph;
        using vg::algorithms::Edge;
        using vg::algorithms::Node;
        using vg::algorithms::compute_min_cut;

        const int seed = 0;
        const int n_iterations = 100;
        // vg test "[Test1]"
        TEST_CASE("Can find a min-cut on a 4 node graph","[Min-cut-graph][Test1]") {
            
            
            /* Let us create following undirected, weighted graph  
                    e0 10  
                0--------1  
                | \      |  
            e1 6|   |5 e4 |15  
                |      \ |  e3
                2--------3  
                    e2 4 
            */
            Graph graph; 
            Edge edge01,edge10,edge02,edge20,edge03,edge30,edge23,edge32,edge13,edge31; //naming convention edge:source:destination
            Node node0, node1, node2, node3;
            size_t V = 4; //size of nodes
            //weights
            edge01.weight = 10;
            edge10.weight = 10;
            edge02.weight = 6;
            edge20.weight = 6;
            edge03.weight = 5;
            edge30.weight = 5;
            edge23.weight = 4;
            edge32.weight = 4;
            edge13.weight = 15;
            edge31.weight = 15;

            

            //other node
            //edge01.other 1, edge0->1 
            edge01.other = 1;
            edge10.other = 0;
            edge02.other = 2;
            edge20.other = 0;
            edge03.other = 3;
            edge30.other = 0;
            edge23.other = 3;
            edge32.other = 2;
            edge13.other = 3;
            edge31.other = 1;

            //node 0
            node0.edges.push_back(edge01);
            node0.edges.push_back(edge02);
            node0.edges.push_back(edge03);
            node0.weight = 21;


            //node 1
            node1.edges.push_back(edge10);
            node1.edges.push_back(edge13);
            node1.weight = 25;

            //node 2
            node2.edges.push_back(edge20);
            node2.edges.push_back(edge23);
            node2.weight = 10;

            //node 3
            node3.edges.push_back(edge30);
            node3.edges.push_back(edge31);
            node3.edges.push_back(edge32);
            node3.weight = 24;

            graph.nodes.push_back(node0);
            graph.nodes.push_back(node1);
            graph.nodes.push_back(node2);
            graph.nodes.push_back(node3);
            
        
            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;

            
          
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 10);
            REQUIRE(disjoint_sets.size() == 2);

            unordered_set<size_t> set1 = {1,2,3};
            unordered_set<size_t> set2 = {0};
            //     REQUIRE(disjoint_sets[1].size() == 3);
            // for (vector<size_t>::iterator it = set1.begin() ; it != set1.end(); ++it){
            //     REQUIRE(disjoint_sets[1].count(*it)==1);
            //     REQUIRE(disjoint_sets[1].count(*it)==1);
            //     REQUIRE(disjoint_sets[1].count(*it)==1);
                
            //     // cout << "dj1 "<< *it << ": " <<disjoint_sets[1].count(*it) <<endl; 
            // }
            // REQUIRE(disjoint_sets[1].size() == 3);

            for (auto& x:disjoint_sets[1] ) {
                
                cout << "set1 has" << x <<endl;
            }

            for (auto& x:disjoint_sets[0] ) {
                cout << "set2 has " << x <<endl;
                
            }
                


            
 
        }

        TEST_CASE("Can find a min-cut on a 9 node graph", "[Min-cut-graph][Test2]") {
            
            
       
            Graph graph; 
            Edge edge01, edge02, edge10, edge12, edge20, edge21, edge23, edge24, edge32, edge42, edge45, edge47, edge54, edge56, edge58, edge65, edge67, edge74, edge76, edge78, edge85, edge87; //naming convention edge:source:destination
            Node node0, node1, node2, node3, node4, node5, node6, node7, node8;
            size_t V = 9; //size of nodes
            
            //weights
            edge01.weight = 100;
            edge02.weight = 200;
            edge10.weight = 100;
            edge12.weight = 400;
            edge20.weight = 200;
            edge21.weight = 400;
            edge23.weight = 300;
            edge24.weight = 5;
            edge32.weight = 300;
            edge42.weight = 5;
            edge45.weight = 400;
            edge47.weight = 100;
            edge54.weight = 400;
            edge56.weight = 700;
            edge58.weight = 500;
            edge65.weight = 700;
            edge67.weight = 200;
            edge74.weight = 100;
            edge76.weight = 200;
            edge78.weight = 100;
            edge85.weight = 500;
            edge87.weight= 100;

            //other 
            edge01.other = 1;
            edge02.other = 2;
            edge10.other = 0;
            edge12.other = 2;
            edge20.other = 0;
            edge21.other = 1;
            edge23.other = 3;
            edge24.other = 4;
            edge32.other = 2;
            edge42.other = 2;
            edge45.other = 5;
            edge47.other = 7;
            edge54.other = 4;
            edge56.other = 6;
            edge58.other = 8;
            edge65.other = 5;
            edge67.other = 7;
            edge74.other = 4;
            edge76.other = 6;
            edge78.other = 8;
            edge85.other = 5;
            edge87.other = 7;


            //node 0
            node0.edges.push_back(edge01);
            node0.edges.push_back(edge02);
            node0.weight = 3;


            //node 1
            node1.edges.push_back(edge10);
            node1.edges.push_back(edge12);
            node1.weight = 5;

            //node 2
            node2.edges.push_back(edge20);
            node2.edges.push_back(edge21);
            node2.edges.push_back(edge23);
            node2.edges.push_back(edge24);
            node2.weight = 14;

            //node 3
            node3.edges.push_back(edge32);
            node3.weight = 3;

            //node4
            node4.edges.push_back(edge42);
            node4.edges.push_back(edge45);
            node4.edges.push_back(edge47);
            node4.weight = 10;

            //node5
            node5.edges.push_back(edge54);
            node5.edges.push_back(edge56);
            node5.edges.push_back(edge58);
            node5.weight = 16;
            
            //node6
            node6.edges.push_back(edge65);
            node6.edges.push_back(edge67);
            node6.weight = 9;
            
            //node7
            node7.edges.push_back(edge74);
            node7.edges.push_back(edge76);
            node7.edges.push_back(edge78);
            node7.weight = 4;

            //node8
            node8.edges.push_back(edge85);
            node8.edges.push_back(edge87);
            node8.weight = 6;
            
            graph.nodes.push_back(node0);
            graph.nodes.push_back(node1);
            graph.nodes.push_back(node2);
            graph.nodes.push_back(node3);
            graph.nodes.push_back(node4);
            graph.nodes.push_back(node5);
            graph.nodes.push_back(node6);
            graph.nodes.push_back(node7);
            graph.nodes.push_back(node8);
        
            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;
            
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 5);
            REQUIRE(disjoint_sets.size() == 2);

            
 
        }
        TEST_CASE("Can find a min-cut on a 2 node graph", "[Min-cut-graph][Test3]") {
            
            
 
            Graph graph; 
            Edge edge01, edge10; //naming convention edge:source:destination
            Node node0, node1;
            size_t V = 2; //size of nodes

            //weights
            edge01.weight = 10;
            edge10.weight = 10;

            //other node
            edge01.other = 1;
            edge10.other = 0;

            //node 0
            node0.edges.push_back(edge01);
            node0.weight = 10;


            //node 1
            node1.edges.push_back(edge10);
            node1.weight = 10;


            graph.nodes.push_back(node0);
            graph.nodes.push_back(node1);

            
            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;
            
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 10);
            REQUIRE(disjoint_sets.size() == 2);
            
 
        }
        TEST_CASE("Can find a min-cut on a 0 node graph, 0 edges","[Min-cut-graph][Test4]") {
                
            
            Graph graph; 
            size_t V = 0; //size of nodes

            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;
            
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 0);
            REQUIRE(disjoint_sets.size() == 0);
            REQUIRE(disjoint_sets.empty() == true);
            
    
        }
        TEST_CASE("Can find a min-cut on a 2 node graph, 0 edges", "[Min-cut-graph][Test5]") {
                
            
            Graph graph; 
            Node node0, node1;
            size_t V = 2; //size of nodes

            //node 0
            node0.weight = 0;


            //node 1
            node1.weight = 0;


            graph.nodes.push_back(node0);
            graph.nodes.push_back(node1);

            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;
            
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 0);
            REQUIRE(disjoint_sets.size() == 0);
            
    
        }
        TEST_CASE("Can find a min-cut on a 4 node graph, 0 edges", "[Min-cut-graph][Test6]") {
                
            Graph graph; 
            Node node0, node1, node2, node3;
            size_t V = 4; //size of nodes

            //node 0
            node0.weight = 0;

            //node 1
            node1.weight = 0;

            //node 2
            node2.weight = 0;

            //node 3
            node3.weight = 0;


            graph.nodes.push_back(node0);
            graph.nodes.push_back(node1);
            graph.nodes.push_back(node2);
            graph.nodes.push_back(node3);


            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;           
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 0);
            REQUIRE(disjoint_sets.size() == 0);
            REQUIRE(disjoint_sets.empty() == true);
    
        }
        TEST_CASE("Min-cut decomposition works on connected component of 4 nodes", "[Min-cut-graph][MCG-Test7]") {
               /* Let us create following undirected, weighted graph  
                    e0 10  
                0--------1  
                | \      |  
            e1 6|   |5 e4 |15  
                |      \ |  e3
                2--------3  
                    e2 4 
            */
            Graph graph; 
            Edge edge01,edge10,edge02,edge20,edge03,edge30,edge23,edge32,edge13,edge31; //naming convention edge:source:destination
            Node node0, node1, node2, node3;
            size_t V = 4; //size of nodes
            //weights
            edge01.weight = 10;
            edge10.weight = 10;
            edge02.weight = 6;
            edge20.weight = 6;
            edge03.weight = 5;
            edge30.weight = 5;
            edge23.weight = 4;
            edge32.weight = 4;
            edge13.weight = 15;
            edge31.weight = 15;

            

            //other node
            //edge01.other 1, edge0->1 
            edge01.other = 1;
            edge10.other = 0;
            edge02.other = 2;
            edge20.other = 0;
            edge03.other = 3;
            edge30.other = 0;
            edge23.other = 3;
            edge32.other = 2;
            edge13.other = 3;
            edge31.other = 1;

            //node 0
            node0.edges.push_back(edge01);
            node0.edges.push_back(edge02);
            node0.edges.push_back(edge03);
            node0.weight = 21;


            //node 1
            node1.edges.push_back(edge10);
            node1.edges.push_back(edge13);
            node1.weight = 25;

            //node 2
            node2.edges.push_back(edge20);
            node2.edges.push_back(edge23);
            node2.weight = 10;

            //node 3
            node3.edges.push_back(edge30);
            node3.edges.push_back(edge31);
            node3.edges.push_back(edge32);
            node3.weight = 24;

            graph.nodes.push_back(node0);
            graph.nodes.push_back(node1);
            graph.nodes.push_back(node2);
            graph.nodes.push_back(node3);
            
        
            //Karger's min-cut
            // min-cut-decom(Graph){

            //     //pair = compute-min_cut(graph)
            //     // base case V == 2 
            //     // get disjoint sets (2)
            //     // conver each disjoint set into a graph using the original graph
            //     // min-cut-decomp (graph1)
            //     // min-cut-decomp (graph2)
            // }
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, n_iterations, seed, V);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;

            
          
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 10);
            REQUIRE(disjoint_sets.size() == 2); 
            
    
        }

    }

}


