/// \file mapping.cpp
///  
/// unit tests for Mappings and their utility functions
///

#include <iostream>
#include <string>
#include "../json2pb.h"
#include "../vg.pb.h"
#include "../path.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Mapping simplification keeps pure deletion edits", "[mapping]") {
    
    string mapping_string = R"({"edit": [{"from_length": 1, "to_length": 1}, {"from_length": 1}]})";
    
    Mapping m;
    json2pb(m, mapping_string.c_str(), mapping_string.size());
    
    auto simple = simplify(m);
    
    REQUIRE(simple.edit_size() == m.edit_size());
    
}
   
}
}
        
