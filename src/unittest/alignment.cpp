/// \file alignment.cpp
///  
/// unit tests for Alignments and their utility functions
///

#include <iostream>
#include <string>
#include "../json2pb.h"
#include "../vg.pb.h"
#include "../alignment.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Alignment simplification removes deletions on the edges of Mappings", "[alignment]") {
    
    string alignment_string = R"(
        {
            "sequence": "A",
            "path": {"mapping": [
                {
                    "position": {"node_id": 1},
                    "edit": [
                        {"from_length": 1, "to_length": 1},
                        {"from_length": 1}
                    ]
                },
                {
                    "position": {"node_id": 2},
                    "edit": [
                        {"from_length": 1, "to_length": 1}
                    ]
                }
            ]}
        }
    )";
    
    Alignment a;
    json2pb(a, alignment_string.c_str(), alignment_string.size());
        
    auto simple = simplify(a);
    REQUIRE(simple.path().mapping_size() == 2);
    REQUIRE(simple.path().mapping(0).edit_size() == 1);
    
}
   
}
}
        
