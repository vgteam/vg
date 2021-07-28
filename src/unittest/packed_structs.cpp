/// \file packed_structs.cpp
///  
/// Unit tests for bit-packed data structures
///

#include <iostream>
#include <random>
#include <unordered_set>
#include <sstream>

#include "vg/io/json2pb.h"
#include "catch.hpp"
#include "randomness.hpp"

#include <bdsg/internal/packed_structs.hpp>

namespace vg {
namespace unittest {
using namespace std;

    TEST_CASE("PackedVector acts the same as STL vector", "[packed]") {
        
        enum vec_op_t {SET = 0, GET = 1, APPEND = 2, POP = 3, SERIALIZE = 4};
        
        default_random_engine prng(test_seed_source());
        uniform_int_distribution<int> op_distr(0, 4);
        
        int num_runs = 1000;
        int num_ops = 200;
        int gets_per_op = 5;
        int sets_per_op = 5;
        int appends_per_op = 3;
        int pops_per_op = 1;
        
        for (size_t i = 0; i < num_runs; i++) {
            
            uint64_t next_val = 0;
            
            vector<uint64_t> std_vec;
            bdsg::PackedVector<> dyn_vec;
            
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
                        
                    case SERIALIZE:
                    {
                        stringstream strm;
                        
                        dyn_vec.serialize(strm);
                        strm.seekg(0);
                        bdsg::PackedVector<> copy_vec(strm);
                        
                        REQUIRE(copy_vec.size() == dyn_vec.size());
                        for (size_t i = 0; i < copy_vec.size(); i++) {
                            REQUIRE(copy_vec.get(i) == dyn_vec.get(i));
                        }
                        break;
                    }
                        
                    default:
                        break;
                }
                
                REQUIRE(std_vec.empty() == dyn_vec.empty());
                REQUIRE(std_vec.size() == dyn_vec.size());
            }
        }
    }
    
    TEST_CASE("PagedVector acts the same as STL vector", "[packed]") {
        
        enum vec_op_t {SET = 0, GET = 1, APPEND = 2, POP = 3, SERIALIZE = 4};
        std::default_random_engine prng(test_seed_source());
        std::uniform_int_distribution<int> op_distr(0, 4);
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
            bdsg::PagedVector<> dyn_vec;
            
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
                        
                    case SERIALIZE:
                    {
                        stringstream strm;
                        
                        dyn_vec.serialize(strm);
                        strm.seekg(0);
                        bdsg::PagedVector<> copy_vec(strm);
                        
                        REQUIRE(copy_vec.size() == dyn_vec.size());
                        for (size_t i = 0; i < copy_vec.size(); i++) {
                            REQUIRE(copy_vec.get(i) == dyn_vec.get(i));
                        }
                        break;
                    }
                        
                    default:
                        break;
                }
                
                REQUIRE(std_vec.empty() == dyn_vec.empty());
                REQUIRE(std_vec.size() == dyn_vec.size());
            }
        }
    }
    
    TEST_CASE("PackedDeque acts the same as STL deque", "[packed]") {
        
        enum deque_op_t {SET = 0, GET = 1, APPEND_LEFT = 2, POP_LEFT = 3, APPEND_RIGHT = 4, POP_RIGHT = 5, SERIALIZE = 6};
        std::default_random_engine prng(test_seed_source());
        std::uniform_int_distribution<int> op_distr(0, 6);
        
        int num_runs = 1000;
        int num_ops = 200;
        int gets_per_op = 5;
        int sets_per_op = 5;
        int appends_per_op = 3;
        int pops_per_op = 1;
        
        for (size_t i = 0; i < num_runs; i++) {
            
            uint64_t next_val = 0;
            
            std::deque<uint64_t> std_deq;
            bdsg::PackedDeque<> suc_deq;
            
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
                        
                    case SERIALIZE:
                    {
                        stringstream strm;
                        
                        suc_deq.serialize(strm);
                        strm.seekg(0);
                        bdsg::PackedDeque<> copy_deq(strm);
                        
                        REQUIRE(copy_deq.size() == suc_deq.size());
                        for (size_t i = 0; i < copy_deq.size(); i++) {
                            REQUIRE(copy_deq.get(i) == suc_deq.get(i));
                        }
                        break;
                    }
                        
                    default:
                        break;
                }
                
                REQUIRE(std_deq.empty() == suc_deq.empty());
                REQUIRE(std_deq.size() == suc_deq.size());
            }
        }
    }
}
}
        
