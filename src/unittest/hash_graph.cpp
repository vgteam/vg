/// \file hash_graph.cpp
///  
/// Unit tests for the HashGraph class.
///

#include <iostream>
#include <sstream>

#include "../json2pb.h"
#include "../hash_graph.hpp"
#include "random_graph.hpp"

#include "algorithms/are_equivalent.hpp"

#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;
    
    TEST_CASE("HashGraph serialization works on randomized graphs", "[hashgraph]") {
        
        int num_graphs = 100;
        int seq_length = 200;
        int num_variants = 30;
        int long_var_length = 10;
        
        for (int i = 0; i < num_graphs; i++) {
            HashGraph graph;
            random_graph(seq_length, long_var_length, num_variants, &graph);
            
            stringstream strm;
            graph.serialize(strm);
            strm.seekg(0);
            HashGraph loaded(strm);
            
            REQUIRE(algorithms::are_equivalent_with_paths(&graph, &loaded));
        }
    }
}
}
        
