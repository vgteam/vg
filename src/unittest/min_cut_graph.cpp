/// \file min_cut_graph.cpp
///  
/// unit tests for min_cut_graph construction and utility functions
///
#include <algorithm>
#include <random>
#include <xg.hpp>
#include <stdio.h>
#include <iostream>
#include "../algorithms/min_cut_graph.hpp"
#include <vector>
#include "catch.hpp"
#include <unordered_set>
#include "sparse_union_find.hpp"
#include <stdlib.h>     
#include <time.h>
#include <chrono>  
#include <iterator>     
// #define debug

namespace vg {
    namespace unittest {
        
        using namespace std;
        using namespace std::chrono; 
        using vg::algorithms::Graph;
        using vg::algorithms::Edge;
        using vg::algorithms::Node;
        using vg::algorithms::compute_min_cut;

        const int seed = 0;
        const int n_iterations = 100;
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

            graph.add_node(0,node0);
            graph.add_node(1,node1);
            graph.add_node(2,node2);
            graph.add_node(3,node3);
            

            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, seed);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;

            //check for reasonable sizes
            size_t max = graph.get_size();
            REQUIRE(disjoint_sets[0].size() < max); // sets are a reasonable size not exceeding number of nodes
            REQUIRE(disjoint_sets[1].size() < max);

            // compute-min-cut will always return 2 sets
            // then min-graph-decomp will remove singleton sets
            REQUIRE(mincut ==10);
            REQUIRE(disjoint_sets.size() ==2); 

           
            //check for disjoint-ness
            //check1: check that first two are disjoint - these are the sets from compute-min-cut
            vector<size_t> v1(disjoint_sets[0].begin(), disjoint_sets[0].end());
            vector<size_t> v2(disjoint_sets[1].begin(), disjoint_sets[1].end());
            sort(v1.begin(),v1.end());
            sort(v2.begin(),v2.end());
            int i =0; 
            int j =0;
            bool is_disjoint = true;
            while (i < v1.size() && j< v2.size() )
            {
                if (v1[i] < v2[j])
                    i++;
                else if (v2[j] < v1[i])
                    j++;
                else /* if *v1[i] == v2[j] */
                    is_disjoint =  false;
            }
            REQUIRE(is_disjoint);
            //check2: check that others are subsets of either the first two
            // this unittest does not recurse (not using min-cut-decomposition yet)


            
        }

        TEST_CASE("Can find a min-cut on a 9 node graph", "[Min-cut-graph][Test2]") {
            
            
       
            Graph graph; 
            Edge edge01, edge02, edge10, edge12, edge20, edge21, edge23, edge24, edge32, edge42, edge45, edge47, edge54, edge56, edge58, edge65, edge67, edge74, edge76, edge78, edge85, edge87; //naming convention edge:source:destination
            Node node0, node1, node2, node3, node4, node5, node6, node7, node8;
            
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
            node0.weight = 300;


            //node 1
            node1.edges.push_back(edge10);
            node1.edges.push_back(edge12);
            node1.weight = 500;

            //node 2
            node2.edges.push_back(edge20);
            node2.edges.push_back(edge21);
            node2.edges.push_back(edge23);
            node2.edges.push_back(edge24);
            node2.weight = 905;

            //node 3
            node3.edges.push_back(edge32);
            node3.weight = 300;

            //node4
            node4.edges.push_back(edge42);
            node4.edges.push_back(edge45);
            node4.edges.push_back(edge47);
            node4.weight = 505;

            //node5
            node5.edges.push_back(edge54);
            node5.edges.push_back(edge56);
            node5.edges.push_back(edge58);
            node5.weight = 1600;
            
            //node6
            node6.edges.push_back(edge65);
            node6.edges.push_back(edge67);
            node6.weight = 900;
            
            //node7
            node7.edges.push_back(edge74);
            node7.edges.push_back(edge76);
            node7.edges.push_back(edge78);
            node7.weight = 400;

            //node8
            node8.edges.push_back(edge85);
            node8.edges.push_back(edge87);
            node8.weight = 600;
            
            graph.add_node(0,node0);
            graph.add_node(1,node1);
            graph.add_node(2,node2);
            graph.add_node(3,node3);
            graph.add_node(4,node4);
            graph.add_node(5,node5);
            graph.add_node(6,node6);
            graph.add_node(7,node7);
            graph.add_node(8,node8);
        
            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, seed);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;
            
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 5);
            REQUIRE(disjoint_sets.size() == 2);
            REQUIRE(disjoint_sets[0].size() < graph.get_size());
            REQUIRE(disjoint_sets[1].size() < graph.get_size());

            //check for disjoint-ness
            //check1: check that first two are disjoint - these are the sets from compute-min-cut
            vector<size_t> v1(disjoint_sets[0].begin(), disjoint_sets[0].end());
            vector<size_t> v2(disjoint_sets[1].begin(), disjoint_sets[1].end());
            sort(v1.begin(),v1.end());
            sort(v2.begin(),v2.end());
            int i =0; 
            int j =0;
            bool is_disjoint = true;
            while (i < v1.size() && j< v2.size() )
            {
                if (v1[i] < v2[j])
                    i++;
                else if (v2[j] < v1[i])
                    j++;
                else /* if *v1[i] == v2[j] */
                    is_disjoint =  false;
            }
            REQUIRE(is_disjoint);
            //check2: check that others are subsets of either the first two
            // this unittest does not recurse (not using min-cut-decomposition yet)



        }
        TEST_CASE("Can find a min-cut on a 2 node graph", "[Min-cut-graph][Test3]") {
            
            
 
            Graph graph; 
            Edge edge01, edge10; //naming convention edge:source:destination
            Node node0, node1;

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


            graph.add_node(0,node0);
            graph.add_node(1,node1);

            
            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, seed);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;

            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 10);
            REQUIRE(disjoint_sets.size() == 2);
            REQUIRE(disjoint_sets[0].size() < graph.get_size());
            REQUIRE(disjoint_sets[1].size() < graph.get_size());

            //check for disjoint-ness
            //check1: check that first two are disjoint - these are the sets from compute-min-cut
            vector<size_t> v1(disjoint_sets[0].begin(), disjoint_sets[0].end());
            vector<size_t> v2(disjoint_sets[1].begin(), disjoint_sets[1].end());
            sort(v1.begin(),v1.end());
            sort(v2.begin(),v2.end());
            int i =0; 
            int j =0;
            bool is_disjoint = true;
            while (i < v1.size() && j< v2.size() )
            {
                if (v1[i] < v2[j])
                    i++;
                else if (v2[j] < v1[i])
                    j++;
                else /* if *v1[i] == v2[j] */
                    is_disjoint =  false;
            }
            REQUIRE(is_disjoint);
            //check2: check that others are subsets of either the first two
            // this unittest does not recurse (not using min-cut-decomposition yet)


            
 
        }
        TEST_CASE("Can find a min-cut on a 0 node graph with 0 edges","[Min-cut-graph][Test4]") {
                
            
            Graph graph; 

            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, seed);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;
            
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 0);
            REQUIRE(disjoint_sets.size() == 0);
            REQUIRE(disjoint_sets.empty() == true);


            
    
        }
        TEST_CASE("Can find a min-cut on a 2 node graph with 0 edges", "[Min-cut-graph][Test5]") {
                
            
            Graph graph; 
            Node node0, node1;

            //node 0
            node0.weight = 0;


            //node 1
            node1.weight = 0;


            graph.add_node(0,node0);
            graph.add_node(1,node1);

            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, seed);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;
            
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 0);
            REQUIRE(disjoint_sets.size() == 0);
            REQUIRE(disjoint_sets.empty() == true);
    
        }
        TEST_CASE("Can find a min-cut on a 4 node graph with 0 edges", "[Min-cut-graph][Test6]") {
                
            Graph graph; 
            Node node0, node1, node2, node3;


            //node 0
            node0.weight = 0;

            //node 1
            node1.weight = 0;

            //node 2
            node2.weight = 0;

            //node 3
            node3.weight = 0;


            graph.add_node(0,node0);
            graph.add_node(1,node1);
            graph.add_node(2,node2);
            graph.add_node(3,node3);


            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, seed);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second;           
            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(mincut == 0);
            REQUIRE(disjoint_sets.size() == 0);
            REQUIRE(disjoint_sets.empty() == true);
    
        }
        TEST_CASE("Min-cut decomposition works on connected component of 4 nodes and ouput correct Gamma (set of subsets of snarls)", "[Min-cut-graph][MCG-Test7]") {
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

            graph.add_node(0,node0);
            graph.add_node(1,node1);
            graph.add_node(2,node2);
            graph.add_node(3,node3);
            
            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, seed);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second; 

            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(disjoint_sets.size() == 2);
            REQUIRE(disjoint_sets[0].size() < graph.get_size());
            REQUIRE(disjoint_sets[1].size() < graph.get_size());

            //check for disjoint-ness
            //check1: check that first two are disjoint - these are the sets from compute-min-cut
            vector<size_t> v1(disjoint_sets[0].begin(), disjoint_sets[0].end());
            vector<size_t> v2(disjoint_sets[1].begin(), disjoint_sets[1].end());
            sort(v1.begin(),v1.end());
            sort(v2.begin(),v2.end());
            int i =0; 
            int j =0;
            bool is_disjoint = true;
            while (i < v1.size() && j< v2.size() )
            {
                if (v1[i] < v2[j])
                    i++;
                else if (v2[j] < v1[i])
                    j++;
                else /* if *v1[i] == v2[j] */
                    is_disjoint =  false;
            }
            REQUIRE(is_disjoint);
            //check2: check that others are subsets of either the first two
            // this unittest does not recurse (not using min-cut-decomposition yet)

            //call min-cut-decomposition
            vector<unordered_set<size_t>> gamma =  min_cut_decomposition(graph, seed);
            
            //convert unordered sets into vectors
            for(unordered_set<size_t> set: gamma){
                
                // turn set into vector 
                vector<size_t> vsub(set.begin(), set.end());
                // sort 
                sort(vsub.begin(), vsub.end());

                //check for equality with first two subsets
                bool equal_to_v1 = vsub.size() == v1.size() && equal(v1.begin(), v1.end(), vsub.begin());
                bool equal_to_v2 = vsub.szie() == v2.size() && equal(v2.begin(), v2.end(), vsub.begin());

                bool pass_equal_to_either = true;
                //if the subset is equal to either of the first two subsets, skip it
                if(equal_to_v1 && equal_to_v2){
                    //equal to both 
                    pass_equal_to_either = false;
                    REQUIRE(pass_equal_to_either);
                }else if(equal_to_v1 || equal_to_v2){
                    //equal to one
                    continue;
                }else{
                    //only check if subset when not equal to first subsets
                    vector<size_t> diff_result_v1(v1.size());
                    vector<size_t> diff_result_v2(v2.size());
                    bool subset_of_v1 = false;
                    bool subset_of_v2 = false;
                    vector<size_t>::iterator it1;
                    vector<size_t>::iterator it2;
                    //equal to neither, check if a subset of either
                    //check subset against v1 & v2
                    //has to be proper subset
                    if(v1.size() > vsub.size()){
                        it1 = set_difference(v1.begin(), v1.end(), vsub.begin(), vsub.end(), diff_result_v1.begin());
                        diff_result_v1.resize(it1-diff_result_v1.begin()); //resize the results vector
                        
                        if(diff_result_v1.size() > 0 && diff_result_v1.size() < v1.size() ){
                            subset_of_v1 = true;
#ifdef debug
                            //check that there are differences between the set and subset
                            //and subset is proper
                            cerr <<"vsub.size() " << vsub.size() << endl;
                            cerr <<"v1.size " << v1.size() << endl;
                            cerr << "The v1 difference has " << (diff_result_v1.size()) << " elements:\n";
                            for(size_t element: diff_result_v1){
                                cerr << element << " ";
                            }
                            cerr << endl;
#endif
                        }
                        
                    }
                    //has to be a proper subset
                    if(v2.size() > vsub.size()){
                        it2 = set_difference(v2.begin(), v2.end(), vsub.begin(), vsub.end(), diff_result_v2.begin());
                        diff_result_v2.resize(it2-diff_result_v2.begin()); //resize the results vector
                        if(diff_result_v2.size() > 0 && diff_result_v2.size() < v2.size()){
                            subset_of_v2 = true;
#ifdef debug
                            cerr <<"vsub.size() " << vsub.size() << endl;
                            cerr <<"v2.size " << v2.size() << endl;
                            cerr << "The v2 difference has " << (diff_result_v2.size()) << " elements:\n";
                            for(size_t element: diff_result_v2){
                                cerr << element <<" ";
                            }
                            cerr << endl;
#endif
                         }
                    }
                    //check that is a subset of one set only
                    bool pass_match_one_only = false;
                    if(subset_of_v1 == false || subset_of_v2 == false){
                        pass_match_one_only = true;
                    }
                    REQUIRE(pass_match_one_only);
                    

                    //require that the set is subset of either of the first two, but not both
                    REQUIRE(subset_of_v1 != subset_of_v2);


                }
                
                
            }

        }

        TEST_CASE("SparseUnionFind can work with non-consecutive node ids", "[SparseUnionFind][MCG-Test8]") {


            /* Let us create following undirected, weighted graph  
                    e0 10  
                5--------6  
                | \      |  
            e1 6|   |5 e4 |15  
                |      \ |  e3
                7--------8  
                    e2 4 
            */
            Graph graph; 
            Edge edge01,edge10,edge02,edge20,edge03,edge30,edge23,edge32,edge13,edge31; //naming convention edge:source:destination
            Node node0, node1, node2, node3;

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

            //id 5
            node0.edges.push_back(edge01);
            node0.edges.push_back(edge02);
            node0.edges.push_back(edge03);
            node0.weight = 21;


            //id 6
            node1.edges.push_back(edge10);
            node1.edges.push_back(edge13);
            node1.weight = 25;

            //id 7
            node2.edges.push_back(edge20);
            node2.edges.push_back(edge23);
            node2.weight = 10;

            //id 8
            node3.edges.push_back(edge30);
            node3.edges.push_back(edge31);
            node3.edges.push_back(edge32);
            node3.weight = 24;

            graph.add_node(5,node0);
            graph.add_node(6,node1);
            graph.add_node(7,node2);
            graph.add_node(8,node3);
            

            SparseUnionFind suf = SparseUnionFind(true, graph.get_node_ids());

#ifdef debug
                cout << "suf has " << suf.size() << " nodes" <<endl;
#endif
                REQUIRE(suf.size() == 4);
                vector<vector<size_t>> all = suf.all_groups(); 
                for(int i = 0 ; i < all.size(); i++){
                    for(int j = 0 ; j < all[i].size(); j++){
#ifdef debug
                        cout <<"group " <<all[i][j] << " has " << suf.group_size(all[i][j]) << " members"<<endl;
#endif
                    }
                        
                }

                suf.union_groups(5, 6);
#ifdef debug
                cout << "union between 5,6" <<endl;
                cout << "suf now has " << suf.all_groups().size() << " groups" <<endl;
#endif
                REQUIRE(suf.all_groups().size() ==3);

                vector<vector<size_t>> new_all = suf.all_groups(); 
                for(int i = 0 ; i < new_all.size(); i++){
                    for(int j = 0 ; j < new_all[i].size(); j++){
#ifdef debug
                        cout <<"group " <<new_all[i][j] << " has " << suf.group_size(new_all[i][j]) << " members"<<endl;
#endif
                    }
                        
                }
                size_t group_5_id = suf.find_group(5);
                size_t group_6_id = suf.find_group(6);
#ifdef debug
                cout << "group 5 is now in " << group_5_id <<endl;
                cout << "group 6 is now in " << group_6_id <<endl;
#endif

                size_t group_5_size = suf.group_size(5);
                size_t group_6_size = suf.group_size(6);
#ifdef debug
                cout << "group 5 size " << group_5_size <<endl;
                cout << "group 6 size " << group_6_size <<endl;
#endif
                
                if(5 != group_5_id){
                    vector<size_t> group_members = suf.group(5);
                    for (size_t i = 0 ; i < group_members.size(); i++){
#ifdef debug
                        cout << "super node is "<< group_5_id <<endl;
                        cout << "member "<< group_members[i] <<endl;
#endif
                    }
                }
                if(6 != group_6_id){
                    vector<size_t> group_members = suf.group(6);
                    for (size_t i = 0 ; i < group_members.size(); i++){
#ifdef debug
                        cout << "super node is "<< group_6_id <<endl;
                        cout << "member "<< group_members[i] <<endl;
#endif
                    }
                }

                suf.union_groups(7, 6);
#ifdef debug
                cout << "union between 7,6" <<endl;
                cout << "suf now has " << suf.all_groups().size() << " groups" <<endl;
#endif
                REQUIRE(suf.all_groups().size() ==2);

                size_t group_7_id = suf.find_group(7);
                size_t group_six_id = suf.find_group(6);
#ifdef debug
                cout << "group 7 is now in " << group_7_id <<endl;
                cout << "group 6 is now in " << group_six_id <<endl;
#endif
                vector<vector<size_t>> g = suf.all_groups(); 
                for(int i = 0 ; i < g.size(); i++){
                    for(int j = 0 ; j < g[i].size(); j++){
#ifdef debug
                        cout <<"group " <<g[i][j] << " has " << suf.group_size(g[i][j]) << " members"<<endl;
#endif
                    }
                        
                }

    
        }
        TEST_CASE("min-cut-decompomposition works on a 9 node graph", "[Min-cut-graph][MCG-Test9]") {
            
            
       
            Graph graph; 
            Edge edge01, edge02, edge10, edge12, edge20, edge21, edge23, edge24, edge32, edge42, edge45, edge47, edge54, edge56, edge58, edge65, edge67, edge74, edge76, edge78, edge85, edge87; //naming convention edge:source:destination
            Node node0, node1, node2, node3, node4, node5, node6, node7, node8;
            
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
            node0.weight = 300;


            //node 1
            node1.edges.push_back(edge10);
            node1.edges.push_back(edge12);
            node1.weight = 500;

            //node 2
            node2.edges.push_back(edge20);
            node2.edges.push_back(edge21);
            node2.edges.push_back(edge23);
            node2.edges.push_back(edge24);
            node2.weight = 905;

            //node 3
            node3.edges.push_back(edge32);
            node3.weight = 300;

            //node4
            node4.edges.push_back(edge42);
            node4.edges.push_back(edge45);
            node4.edges.push_back(edge47);
            node4.weight = 505;

            //node5
            node5.edges.push_back(edge54);
            node5.edges.push_back(edge56);
            node5.edges.push_back(edge58);
            node5.weight = 1600;
            
            //node6
            node6.edges.push_back(edge65);
            node6.edges.push_back(edge67);
            node6.weight = 900;
            
            //node7
            node7.edges.push_back(edge74);
            node7.edges.push_back(edge76);
            node7.edges.push_back(edge78);
            node7.weight = 400;

            //node8
            node8.edges.push_back(edge85);
            node8.edges.push_back(edge87);
            node8.weight = 600;
            
            graph.add_node(0,node0);
            graph.add_node(1,node1);
            graph.add_node(2,node2);
            graph.add_node(3,node3);
            graph.add_node(4,node4);
            graph.add_node(5,node5);
            graph.add_node(6,node6);
            graph.add_node(7,node7);
            graph.add_node(8,node8);
        
            //Karger's min-cut
            pair<vector<unordered_set<size_t>>, size_t> to_recv= compute_min_cut(graph, seed);
            vector<unordered_set<size_t>> disjoint_sets = to_recv.first;
            size_t mincut = to_recv.second; 

            //returns the number of minimum edge cuts required to devide graph into two disjoint connected subgraphs 
            REQUIRE(disjoint_sets.size() == 2);
            REQUIRE(disjoint_sets[0].size() < graph.get_size());
            REQUIRE(disjoint_sets[1].size() < graph.get_size());

            //check for disjoint-ness
            //check1: check that first two are disjoint - these are the sets from compute-min-cut
            vector<size_t> v1(disjoint_sets[0].begin(), disjoint_sets[0].end());
            vector<size_t> v2(disjoint_sets[1].begin(), disjoint_sets[1].end());
            sort(v1.begin(),v1.end());
            sort(v2.begin(),v2.end());
            int i =0; 
            int j =0;
            bool is_disjoint = true;
            while (i < v1.size() && j< v2.size() )
            {
                if (v1[i] < v2[j])
                    i++;
                else if (v2[j] < v1[i])
                    j++;
                else /* if *v1[i] == v2[j] */
                    is_disjoint =  false;
            }
            REQUIRE(is_disjoint);
            //check2: check that others are subsets of either the first two
            // this unittest does not recurse (not using min-cut-decomposition yet)

            //call min-cut-decomposition
            vector<unordered_set<size_t>> gamma =  min_cut_decomposition(graph, seed);
            
            //convert unordered sets into vectors
            for(unordered_set<size_t> set: gamma){
                
                // turn set into vector 
                vector<size_t> vsub(set.begin(), set.end());
                // sort 
                sort(vsub.begin(), vsub.end());

                //check for equality with first two subsets
                bool equal_to_v1 = vsub.size() == v1.size() && equal(v1.begin(), v1.end(), vsub.begin());
                bool equal_to_v2 = vsub.size() == v2.size() && equal(v2.begin(), v2.end(), vsub.begin());

                bool pass_equal_to_either = true;
                //if the subset is equal to either of the first two subsets, skip it
                if(equal_to_v1 && equal_to_v2){
                    //equal to both 
                    pass_equal_to_either = false;
                    REQUIRE(pass_equal_to_either);
                }else if(equal_to_v1 || equal_to_v2){
                    //equal to one
                    continue;
                }else{
                    //only check if subset when not equal to first subsets
                    vector<size_t> diff_result_v1(v1.size());
                    vector<size_t> diff_result_v2(v2.size());
                    bool subset_of_v1 = false;
                    bool subset_of_v2 = false;
                    vector<size_t>::iterator it1;
                    vector<size_t>::iterator it2;
                    //equal to neither, check if a subset of either
                    //check subset against v1 & v2
                    //has to be proper subset
                    if(v1.size() > vsub.size()){
                        it1 = set_difference(v1.begin(), v1.end(), vsub.begin(), vsub.end(), diff_result_v1.begin());
                        diff_result_v1.resize(it1-diff_result_v1.begin()); //resize the results vector
                        
                        if(diff_result_v1.size() > 0 && diff_result_v1.size() < v1.size() ){
                            subset_of_v1 = true;
#ifdef debug
                            //check that there are differences between the set and subset
                            //and subset is proper
                            cerr <<"vsub.size() " << vsub.size() << endl;
                            cerr <<"v1.size " << v1.size() << endl;
                            cerr << "The v1 difference has " << (diff_result_v1.size()) << " elements:\n";
                            for(size_t element: diff_result_v1){
                                cerr << element << " ";
                            }
                            cerr << endl;
#endif
                        }
                        
                    }
                    //has to be a proper subset
                    if(v2.size() > vsub.size()){
                        it2 = set_difference(v2.begin(), v2.end(), vsub.begin(), vsub.end(), diff_result_v2.begin());
                        diff_result_v2.resize(it2-diff_result_v2.begin()); //resize the results vector
                        if(diff_result_v2.size() > 0 && diff_result_v2.size() < v2.size()){
                            subset_of_v2 = true;
#ifdef debug
                            cerr <<"vsub.size() " << vsub.size() << endl;
                            cerr <<"v2.size " << v2.size() << endl;
                            cerr << "The v2 difference has " << (diff_result_v2.size()) << " elements:\n";
                            for(size_t element: diff_result_v2){
                                cerr << element <<" ";
                            }
                            cerr << endl;
#endif
                         }
                    }
                    //check that is a subset of one set only
                    bool pass_match_one_only = false;
                    if(subset_of_v1 == false || subset_of_v2 == false){
                        pass_match_one_only = true;
                    }
                    REQUIRE(pass_match_one_only);
                    

                    //require that the set is subset of either of the first two, but not both
                    REQUIRE(subset_of_v1 != subset_of_v2);


                }
            }
            

        }
        TEST_CASE("min-cut-decompomposition works on a 1000 node graph", "[Min-cut-graph][MCG-Test10]") {
            Graph graph;
            size_t max_nodes = 1000;
            for(size_t i = 0; i<max_nodes; i++){
                
                //assign edge a random edge weight using a rand num generator 
                int random_weight = 0;  
                // 0 weights not allowed  
                while(random_weight == 0){
                    random_weight = rand() % 1000; 
                }
// #ifdef debug
//                 cout << "============================================================================= " << endl;
//                 cout << "random weight" << random_weight << endl;        
//                 cout << "============================================================================= " << endl;
// #endif  
                //assign it an `other` -> next for forward , and prev for backward edge
                //case 1 : first node
                if(i == 0 ){
                    Edge forward_edge;
                    forward_edge.weight = random_weight;
                    size_t next_node = i+1;
                    forward_edge.other = next_node;
                    Node node;
                    node.edges.push_back(forward_edge);
                    node.weight = random_weight; //only has one edge 
                    graph.add_node(i,node);
                }
                //case 2: last node
                else if (i ==max_nodes -1){
                    Edge backward_edge;
                    size_t prev_node_id = i-1;
                    backward_edge.other = prev_node_id;
                    Node prev_node = graph.get_node_by_id(prev_node_id);
                    //get the prev node edge weight going into current weight by checking which of the edges has `other` equal to current 
                    size_t prev_node_weight = graph.get_weight_using_other(prev_node, i);
                    backward_edge.weight = prev_node_weight;
                    Node node;
                    node.edges.push_back(backward_edge);
                    node.weight = prev_node_weight;
                    graph.add_node(i,node);
                }

                //case 3: middle node
                else{
                    Edge forward_edge, backward_edge;//outward edges from node <-()->
                    Node node;
                   
                    forward_edge.weight = random_weight;
                    size_t prev_node_id = i-1;
                    size_t next_node_id = i+1;
                    backward_edge.other = prev_node_id;
                    forward_edge.other = next_node_id;
                    Node prev_node = graph.get_node_by_id(prev_node_id);
                    size_t prev_node_weight = graph.get_weight_using_other(prev_node, i);
//  #ifdef debug
//                 cout << "============================================================================= " << endl;
//                 cout << "prev weight" << prev_node_weight << endl;        
//                 cout << "============================================================================= " << endl;
// #endif 
                    backward_edge.weight = prev_node_weight;
                    node.edges.push_back(backward_edge);
                    node.edges.push_back(forward_edge);
                    node.weight = prev_node_weight + random_weight;
                    graph.add_node(i,node);

                }
                
            }
           
            // Get starting timepoint 
            auto start = high_resolution_clock::now(); 
            //call min-cut-decomposition
            vector<unordered_set<size_t>> to_recv =  min_cut_decomposition(graph, seed);
            // Get ending timepoint 
            auto stop = high_resolution_clock::now(); 
            auto duration = duration_cast<microseconds>(stop - start); 
#ifdef debug
            cout << "Size of Gamma "<< to_recv.size() << endl;
            for(size_t i = 0; i < to_recv.size(); i++){
                for (auto& x:to_recv[i] ) {

                cout << "set "<<i<<" has" << x <<endl;
                }
            }
            cout << "Time taken by function: " << duration.count() / 1000000 << " seconds" <<endl;
#endif 
            
        }

    }

}


