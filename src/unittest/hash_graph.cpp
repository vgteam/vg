/// \file hash_graph.cpp
///  
/// Unit tests for the HashGraph class.
///

#include <iostream>
#include <sstream>
#include <random>

#include "vg/io/json2pb.h"
#include "random_graph.hpp"
#include "randomness.hpp"

#include "bdsg/hash_graph.hpp"

#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;
    
    TEST_CASE("HashGraph serialization works on randomized graphs", "[hashgraph]") {
        
        int num_graphs = 100;
        int seq_length = 200;
        int num_variants = 30;
        int long_var_length = 10;
        
        default_random_engine gen(test_seed_source());
        uniform_int_distribution<int> circ_distr(0, 1);
        
        for (int i = 0; i < num_graphs; i++) {
            bdsg::HashGraph graph;
            random_graph(seq_length, long_var_length, num_variants, &graph);
            
            graph.for_each_path_handle([&](const path_handle_t& path) {
                graph.set_circularity(path, circ_distr(gen));
            });
            
            stringstream strm;
            graph.serialize(strm);
            strm.seekg(0);
            bdsg::HashGraph loaded(strm);
            REQUIRE(handlealgs::are_equivalent_with_paths(&graph, &loaded));
        }
    }
}
}
        
