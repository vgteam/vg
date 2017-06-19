//
//  suffix_tree.cpp
//  
//
//  Created by Jordan Eizenga on 1/25/17.
//
//

#include <stdio.h>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include "utility.hpp"

#include "catch.hpp"

using namespace std;

namespace vg {
    namespace unittest {
        
        vector<pair<size_t, size_t>> random_unions(size_t size) {
            
            int num_pairs = size * size;
            
            random_device rd;
            default_random_engine gen(rd());
            uniform_int_distribution<int> distr(0, num_pairs);
            
            vector<pair<size_t, size_t>> pairs;
            vector<size_t> all_options(num_pairs);
            
            for (size_t i = 0; i < num_pairs; i++) {
                all_options[i] = i;
            }
            
            std::shuffle(all_options.begin(), all_options.end(), gen);
            
            int num_pairs_to_select = distr(gen);
            
            for (int i = 0; i < num_pairs_to_select; i++) {
                pairs.emplace_back(all_options[i] / size, all_options[i] % size);
            }
            
            return pairs;
        }
        
        TEST_CASE("Union-Find builds without crashing", "[unionfind]") {
            
            UnionFind union_find(100);
            
        }
        
        TEST_CASE("Union-Find produces correct groups in a hand-selected case", "[unionfind]") {
            
            UnionFind union_find(10);
            
            SECTION("Union-Find maintains invariants after a single union") {
                
                REQUIRE(union_find.find(0) != union_find.find(1));
                
                REQUIRE(union_find.group_size(0) == 1);
                REQUIRE(union_find.group_size(1) == 1);
                
                union_find.union(0, 1);
                
                REQUIRE(union_find.find(0) == union_find.find(1));
                
                REQUIRE(union_find.group_size(0) == 2);
                REQUIRE(union_find.group_size(1) == 2);
            }
            
            SECTION("Union-Find maintains invariants after several more unions") {
                
                union_find.union(2, 3);
                union_find.union(3, 4);
                
                union_find.union(5, 6);
                
                REQUIRE(union_find.group_size(4) == 3);
                REQUIRE(union_find.find(2) == union_find.find(3));
                REQUIRE(union_find.find(5) == union_find.find(6));
            }
            
            SECTION("Union-Find maintains invariants after a null union") {
                
                union_find.union(2, 4);
                
                REQUIRE(union_find.group_size(4) == 3);
                REQUIRE(union_find.find(2) == union_find.find(3));
                REQUIRE(union_find.find(3) == union_find.find(4));
            }
            
            SECTION("Union-Find can correctly extract a group") {
                
                vector<size_t> correct_group{2, 3, 4};
                vector<size_t> group = union_find.group(3);
                std::sort(group.begin(), group.end());
                REQUIRE(group.size() == correct_group.size());
                REQUIRE(std::equal(group.begin(), group.end(), correct_group.begin()));
            }
        }
        
        TEST_CASE("Union-Find produces identical results when operations are performed in different orders", "[unionfind]") {
            
            SECTION("Single group extraction and all group extraction can be performed in either order") {
                
                UnionFind union_find_1(10);
                UnionFind union_find_2(10);
                
                union_find_1.union(0, 1);
                union_find_1.union(2, 1);
                union_find_1.union(3, 2);
                union_find_1.union(4, 5);
                union_find_1.union(7, 6);
                union_find_1.union(7, 8);
                union_find_1.union(7, 9);
                
                union_find_2.union(0, 1);
                union_find_2.union(2, 1);
                union_find_2.union(3, 2);
                union_find_2.union(4, 5);
                union_find_2.union(7, 6);
                union_find_2.union(7, 8);
                union_find_2.union(7, 9);
                
                vector<size_t> group_of_0_1 = union_find_1.group(0);
                vector<size_t> group_of_4_1 = union_find_1.group(4);
                vector<size_t> group_of_9_1 = union_find_1.group(9);
                
                vector<vector<size_t>> all_groups_1 = union_find_1.all_groups();
                
                vector<vector<size_t>> all_groups_2 = union_find_2.all_groups();
                
                vector<size_t> group_of_0_2 = union_find_2.group(0);
                vector<size_t> group_of_4_2 = union_find_2.group(4);
                vector<size_t> group_of_9_2 = union_find_2.group(9);
                
                vector<size_t> group_of_0_1_from_all, group_of_4_1_from_all, group_of_9_1_from_all;
                vector<size_t> group_of_0_2_from_all, group_of_4_2_from_all, group_of_9_2_from_all;
                
                for (vector<size_t>& vec : all_groups_1) {
                    for (size_t i : vec) {
                        if (i == 0) {
                            group_of_0_1_from_all = vec;
                        }
                        else if (i == 4) {
                            group_of_4_1_from_all = vec;
                        }
                        else if (i == 9) {
                            group_of_9_1_from_all = vec;
                        }
                    }
                }
                
                for (vector<size_t>& vec : all_groups_2) {
                    for (size_t i : vec) {
                        if (i == 0) {
                            group_of_0_2_from_all = vec;
                        }
                        else if (i == 4) {
                            group_of_4_2_from_all = vec;
                        }
                        else if (i == 9) {
                            group_of_9_2_from_all = vec;
                        }
                    }
                }
                
                std::sort(group_of_0_1.begin(), group_of_0_1.end());
                std::sort(group_of_4_1.begin(), group_of_4_1.end());
                std::sort(group_of_9_1.begin(), group_of_9_1.end());
                
                std::sort(group_of_0_2.begin(), group_of_0_2.end());
                std::sort(group_of_4_2.begin(), group_of_4_2.end());
                std::sort(group_of_9_2.begin(), group_of_9_2.end());
                
                std::sort(group_of_0_1_from_all.begin(), group_of_0_1_from_all.end());
                std::sort(group_of_4_1_from_all.begin(), group_of_4_1_from_all.end());
                std::sort(group_of_9_1_from_all.begin(), group_of_9_1_from_all.end());
                
                std::sort(group_of_0_2_from_all.begin(), group_of_0_2_from_all.end());
                std::sort(group_of_4_2_from_all.begin(), group_of_4_2_from_all.end());
                std::sort(group_of_9_2_from_all.begin(), group_of_9_2_from_all.end());
                
                REQUIRE(group_of_0_1.size() == group_of_0_2.size());
                REQUIRE(group_of_0_1.size() == group_of_0_1_from_all.size());
                REQUIRE(group_of_0_1.size() == group_of_0_2_from_all.size());
                
                REQUIRE(group_of_4_1.size() == group_of_4_2.size());
                REQUIRE(group_of_4_1.size() == group_of_4_1_from_all.size());
                REQUIRE(group_of_4_1.size() == group_of_4_2_from_all.size());
                
                REQUIRE(group_of_9_1.size() == group_of_9_2.size());
                REQUIRE(group_of_9_1.size() == group_of_9_1_from_all.size());
                REQUIRE(group_of_9_1.size() == group_of_9_2_from_all.size());
                
                REQUIRE(std::equal(group_of_0_1.begin(), group_of_0_1.end(), group_of_0_1_from_all.begin()));
                REQUIRE(std::equal(group_of_0_1.begin(), group_of_0_1.end(), group_of_0_2_from_all.begin()));
                REQUIRE(std::equal(group_of_0_1.begin(), group_of_0_1.end(), group_of_0_2.begin()));
                
                REQUIRE(std::equal(group_of_4_1.begin(), group_of_4_1.end(), group_of_4_1_from_all.begin()));
                REQUIRE(std::equal(group_of_4_1.begin(), group_of_4_1.end(), group_of_4_2_from_all.begin()));
                REQUIRE(std::equal(group_of_4_1.begin(), group_of_4_1.end(), group_of_4_2.begin()));
                
                REQUIRE(std::equal(group_of_9_1.begin(), group_of_9_1.end(), group_of_9_1_from_all.begin()));
                REQUIRE(std::equal(group_of_9_1.begin(), group_of_9_1.end(), group_of_9_2_from_all.begin()));
                REQUIRE(std::equal(group_of_9_1.begin(), group_of_9_1.end(), group_of_9_2.begin()));
            }
            
            SECTION("Group extraction and group size queries can be performed in either order") {
                UnionFind union_find_1(10);
                UnionFind union_find_2(10);
                
                union_find_1.union(0, 1);
                union_find_1.union(2, 1);
                union_find_1.union(3, 2);
                union_find_1.union(4, 5);
                union_find_1.union(7, 6);
                union_find_1.union(7, 8);
                union_find_1.union(7, 9);
                
                union_find_2.union(0, 1);
                union_find_2.union(2, 1);
                union_find_2.union(3, 2);
                union_find_2.union(4, 5);
                union_find_2.union(7, 6);
                union_find_2.union(7, 8);
                union_find_2.union(7, 9);
                
                REQUIRE(union_find_1.group_size(0) == union_find_2.group(0).size());
                REQUIRE(union_find_2.group_size(0) == union_find_1.group(0).size());
                REQUIRE(union_find_1.group_size(0) == union_find_1.group(0).size());
                REQUIRE(union_find_2.group_size(0) == union_find_2.group(0).size());
                
                REQUIRE(union_find_1.group_size(3) == union_find_2.group(3).size());
                REQUIRE(union_find_2.group_size(3) == union_find_1.group(3).size());
                REQUIRE(union_find_1.group_size(3) == union_find_1.group(3).size());
                REQUIRE(union_find_2.group_size(3) == union_find_2.group(3).size());
                
                REQUIRE(union_find_1.group_size(9) == union_find_2.group(9).size());
                REQUIRE(union_find_2.group_size(9) == union_find_1.group(9).size());
                REQUIRE(union_find_1.group_size(9) == union_find_1.group(9).size());
                REQUIRE(union_find_2.group_size(9) == union_find_2.group(9).size());
            }
        }
        
        TEST_CASE("Union-Find maintains invariants in many randomized sequences of unions", "[unionfind]") {
            
            for (size_t repetition = 0; repetition < 3000; repetition++) {
                
                UnionFind union_find(50);
                
                vector<pair<size_t, size_t>> unions = random_unions(union_find.size());
                
                for (pair<size_t, size_t> idxs : unions) {
                    union_find.union(idxs.first, idxs.second);
                }
                
                vector<vector<size_t>> groups_direct(union_find.size());
                vector<vector<size_t>> groups_from_all(union_find.size());
                // alternate the order these operations are done in just in case the
                // internal changes during find functions affects the outcome
                if (repetition % 2 == 0) {
                    for (size_t i = 0; i < union_find.size(); i++) {
                        groups_direct[i] = union_find.group(i);
                    }
                    
                    vector<vector<size_t>> groups = union_find.all_groups();
                    for (vector<size_t>& group : groups) {
                        for (size_t i : group) {
                            groups_from_all[i] = group;
                        }
                    }
                }
                else {
                    vector<vector<size_t>> groups = union_find.all_groups();
                    for (vector<size_t>& group : groups) {
                        for (size_t i : group) {
                            groups_from_all[i] = group;
                        }
                    }
                    
                    for (size_t i = 0; i < union_find.size(); i++) {
                        groups_direct[i] = union_find.group(i);
                    }
                }
                
                for (size_t i = 0; i < union_find.size(); i++) {
                    
                    std::sort(groups_from_all[i].begin(), groups_from_all[i].end());
                    std::sort(groups_direct[i].begin(), groups_direct[i].end());
                    
                    if (!std::equal(groups_from_all[i].begin(), groups_from_all[i].end(),
                                    groups_direct[i].begin())) {
                        cerr << "FAILURE mismatch in repetition " << repetition << " in the following groups containing " << i << endl;
                        cerr << "computed directly: ";
                        for (size_t j : groups_direct[i]) {
                            cerr << j << " ";
                        }
                        cerr << endl;
                        cerr << "computed in batch: ";
                        for (size_t j : groups_from_all[i]) {
                            cerr << j << " ";
                        }
                        cerr << endl;
                        cerr << "groups formed by unions: "
                        for (pair<size_t, size_t> idxs : unions) {
                            cerr << "(" << idxs.first << "," << idxs.second << ") ";
                        }
                        cerr << endl;
                    }
                    
                    REQUIRE(union_find.group_size(i) == groups_from_all[i].size());
                    REQUIRE(union_find.group_size(i) == groups_direct[i].size());
                    REQUIRE(groups_from_all[i].size() == groups_direct[i].size());
                    
                    REQUIRE(std::equal(groups_from_all[i].begin(), groups_from_all[i].end(),
                                       groups_direct[i].begin()));
                }
            }
        }
    }
}
