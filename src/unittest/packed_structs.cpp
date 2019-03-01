/// \file packed_structs.cpp
///  
/// Unit tests for bit-packed data structures
///

#include <iostream>
#include <random>
#include <unordered_set>

#include "../json2pb.h"
#include "../packed_structs.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

    TEST_CASE("PackedVector acts the same as STL vector", "[packed]") {
        
        enum vec_op_t {SET = 0, GET = 1, APPEND = 2, POP = 3};
        
        random_device rd;
        default_random_engine prng(rd());
        uniform_int_distribution<int> op_distr(0, 3);
        
        int num_runs = 1000;
        int num_ops = 200;
        int gets_per_op = 5;
        int sets_per_op = 5;
        int appends_per_op = 3;
        int pops_per_op = 1;
        
        for (size_t i = 0; i < num_runs; i++) {
            
            uint64_t next_val = 0;
            
            vector<uint64_t> std_vec;
            PackedVector dyn_vec;
            
            for (size_t j = 0; j < num_ops; j++) {
                
                vec_op_t op = (vec_op_t) op_distr(prng);
                switch (op) {
                    case SET:
                        if (!std_vec.empty()) {
                            for (size_t k = 0; k < sets_per_op; k++) {
                                size_t idx = prng() % dyn_vec.size();
                                std_vec[idx] = next_val;
                                dyn_vec.set(idx, next_val);
                                next_val++;
                            }
                        }
                        
                        break;
                        
                    case GET:
                        if (!std_vec.empty()) {
                            for (size_t k = 0; k < gets_per_op; k++) {
                                size_t idx = prng() % dyn_vec.size();
                                REQUIRE(std_vec[idx] == dyn_vec.get(idx));
                                next_val++;
                            }
                        }
                        
                        break;
                        
                    case APPEND:
                        for (size_t k = 0; k < appends_per_op; k++) {
                            std_vec.push_back(next_val);
                            dyn_vec.append(next_val);
                            next_val++;
                        }
                        
                        break;
                        
                    case POP:
                        if (!std_vec.empty()) {
                            for (size_t k = 0; k < pops_per_op; k++) {
                                std_vec.pop_back();
                                dyn_vec.pop();
                            }
                        }
                        
                        break;
                        
                    default:
                        break;
                }
                
                REQUIRE(std_vec.empty() == dyn_vec.empty());
                REQUIRE(std_vec.size() == dyn_vec.size());
            }
        }
    }
    
    TEST_CASE("PagedVector acts the same as STL vector", "[packed]") {
        
        enum vec_op_t {SET = 0, GET = 1, APPEND = 2, POP = 3};
        std::random_device rd;
        std::default_random_engine prng(rd());
        std::uniform_int_distribution<int> op_distr(0, 3);
        std::uniform_int_distribution<int> page_distr(1, 5);
        std::uniform_int_distribution<int> val_distr(0, 100);
        
        int num_runs = 1000;
        int num_ops = 200;
        int gets_per_op = 5;
        int sets_per_op = 5;
        int appends_per_op = 3;
        int pops_per_op = 1;
        
        for (size_t i = 0; i < num_runs; i++) {
            
            uint64_t next_val = val_distr(prng);
            
            std::vector<uint64_t> std_vec;
            PagedVector dyn_vec(page_distr(prng));
            
            for (size_t j = 0; j < num_ops; j++) {
                
                vec_op_t op = (vec_op_t) op_distr(prng);
                switch (op) {
                    case SET:
                        if (!std_vec.empty()) {
                            for (size_t k = 0; k < sets_per_op; k++) {
                                size_t idx = prng() % dyn_vec.size();
                                std_vec[idx] = next_val;
                                dyn_vec.set(idx, next_val);
                                next_val = val_distr(prng);
                            }
                        }
                        
                        break;
                        
                    case GET:
                        if (!std_vec.empty()) {
                            for (size_t k = 0; k < gets_per_op; k++) {
                                size_t idx = prng() % dyn_vec.size();
                                REQUIRE(std_vec[idx] == dyn_vec.get(idx));
                                next_val = val_distr(prng);
                            }
                        }
                        
                        break;
                        
                    case APPEND:
                        for (size_t k = 0; k < appends_per_op; k++) {
                            std_vec.push_back(next_val);
                            dyn_vec.append(next_val);
                            next_val = val_distr(prng);
                        }
                        
                        break;
                        
                    case POP:
                        if (!std_vec.empty()) {
                            for (size_t k = 0; k < pops_per_op; k++) {
                                std_vec.pop_back();
                                dyn_vec.pop();
                            }
                        }
                        
                        break;
                        
                    default:
                        break;
                }
                
                REQUIRE(std_vec.empty() == dyn_vec.empty());
                REQUIRE(std_vec.size() == dyn_vec.size());
            }
        }
    }
    
    TEST_CASE("PackedDeque acts the same as STL deque", "[packed]") {
        
        enum deque_op_t {SET = 0, GET = 1, APPEND_LEFT = 2, POP_LEFT = 3, APPEND_RIGHT = 4, POP_RIGHT = 5};
        std::random_device rd;
        std::default_random_engine prng(rd());
        std::uniform_int_distribution<int> op_distr(0, 5);
        
        int num_runs = 1000;
        int num_ops = 200;
        int gets_per_op = 5;
        int sets_per_op = 5;
        int appends_per_op = 3;
        int pops_per_op = 1;
        
        for (size_t i = 0; i < num_runs; i++) {
            
            uint64_t next_val = 0;
            
            std::deque<uint64_t> std_deq;
            PackedDeque suc_deq;
            
            for (size_t j = 0; j < num_ops; j++) {
                
                deque_op_t op = (deque_op_t) op_distr(prng);
                switch (op) {
                    case SET:
                        if (!std_deq.empty()) {
                            for (size_t k = 0; k < sets_per_op; k++) {
                                size_t idx = prng() % std_deq.size();
                                std_deq[idx] = next_val;
                                suc_deq.set(idx, next_val);
                                next_val++;
                            }
                        }
                        
                        break;
                        
                    case GET:
                        if (!std_deq.empty()) {
                            for (size_t k = 0; k < gets_per_op; k++) {
                                size_t idx = prng() % std_deq.size();
                                REQUIRE(std_deq[idx] == suc_deq.get(idx));
                                next_val++;
                            }
                        }
                        
                        break;
                        
                    case APPEND_LEFT:
                        for (size_t k = 0; k < appends_per_op; k++) {
                            std_deq.push_front(next_val);
                            suc_deq.append_front(next_val);
                            next_val++;
                        }
                        
                        break;
                        
                    case POP_LEFT:
                        for (size_t k = 0; k < pops_per_op && !std_deq.empty(); k++) {
                            std_deq.pop_front();
                            suc_deq.pop_front();
                        }
                        
                        break;
                        
                    case APPEND_RIGHT:
                        for (size_t k = 0; k < appends_per_op; k++) {
                            std_deq.push_back(next_val);
                            suc_deq.append_back(next_val);
                            next_val++;
                        }
                        
                        break;
                        
                    case POP_RIGHT:
                        for (size_t k = 0; k < pops_per_op && !std_deq.empty(); k++) {
                            std_deq.pop_back();
                            suc_deq.pop_back();
                        }
                        
                        break;
                        
                    default:
                        break;
                }
                
                REQUIRE(std_deq.empty() == suc_deq.empty());
                REQUIRE(std_deq.size() == suc_deq.size());
            }
        }
    }
    
    TEST_CASE("PackedSplayTree acts the same as STL map", "[packed]") {
        
        enum tree_op_t {INSERT = 0, ERASE = 1, ACCESS = 2};
        std::random_device rd;
        std::default_random_engine prng(rd());
        std::uniform_int_distribution<int> op_distr(0, 2);
        
        int num_runs = 1000;
        int num_ops = 100;
        
        int inserts_per_op = 3;
        int erases_per_op = 1;
        int accesses_per_op = 5;
        
        size_t value_bank_size = num_ops * inserts_per_op + 1;
        size_t max_value = 10 * value_bank_size;
        
        for (size_t i = 0; i < num_runs; i++) {
            
            std::vector<size_t> bank(max_value);
            for (size_t j = 0; j < bank.size(); j++) {
                // only even so we can use lower bound on odds
                bank[j] = j * 2;
            }
            std::shuffle(bank.begin(), bank.end(), prng);
            bank.resize(value_bank_size);
            size_t bank_idx = 0;
            size_t erase_idx = 0;
            
            std::map<size_t, size_t> std_map;
            PackedSplayTree suc_map;
            
            std::unordered_set<size_t> erased;
            
            for (size_t j = 0; j < num_ops; j++) {
                
                tree_op_t op = (tree_op_t) op_distr(prng);
                switch (op) {
                    case INSERT:
                        for (size_t k = 0; k < inserts_per_op; k++) {
                            std_map[bank[bank_idx]] = bank_idx;
                            suc_map.insert(bank[bank_idx], bank_idx);
                            bank_idx++;
                        }
                        
                        break;
                        
                    case ERASE:
                        if (!std_map.empty()) {
                            for (size_t k = 0; k < erases_per_op; k++) {
                                std_map.erase(bank[erase_idx]);
                                suc_map.erase(bank[erase_idx]);
                                erased.insert(bank[erase_idx]);
                                erase_idx++;
                            }
                        }
                        
                        break;
                        
                    case ACCESS:
                        if (!std_map.empty()) {
                            for (size_t k = 0; k < accesses_per_op; k++) {
                                size_t access_val = bank[prng() % bank_idx];
                                size_t x = suc_map.find(access_val);
                                auto xi = std_map.find(access_val);
                                
                                
                                
                                if (erased.count(access_val)) {
                                    REQUIRE(x == 0);
                                }
                                else {
                                    REQUIRE(suc_map.get_key(x) == xi->first);
                                    REQUIRE(suc_map.get_value(x) == xi->second);
                                }
                                
                                
                                size_t y = suc_map.first_lower(access_val);
                                auto yi = std_map.lower_bound(access_val);
                                
                                if (!erased.count(access_val)) {
                                    REQUIRE(suc_map.get_key(y) == yi->first);
                                    REQUIRE(suc_map.get_value(y) == yi->second);
                                }
                                
                                size_t z = suc_map.first_lower(access_val + 1);
                                auto zi = std_map.lower_bound(access_val + 1);
                                
                                if (!erased.count(access_val)) {
                                    zi--;
                                    REQUIRE(suc_map.get_key(z) == zi->first);
                                    REQUIRE(suc_map.get_value(z) == zi->second);
                                }
                                
                                if (x != 0) {
                                    size_t w = suc_map.next(x);
                                    auto wi = xi;
                                    wi++;
                                    if (wi == std_map.end()) {
                                        REQUIRE(w == 0);
                                    }
                                    else {
                                        REQUIRE(suc_map.get_key(w) == wi->first);
                                        REQUIRE(suc_map.get_value(w) == wi->second);
                                    }
                                }
                            }
                        }
                        
                        break;
                        
                    default:
                        break;
                }
                
                REQUIRE(std_map.empty() == suc_map.empty());
                REQUIRE(std_map.size() == suc_map.size());
            }
        }
    }
}
}
        
