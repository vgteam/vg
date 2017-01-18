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

TEST_CASE("Mapping simplification fails when mapping position is not set", "[mapping]") {
    
    // We can't remove leading deletions if there's no position to update
    string mapping_string = R"({"edit": [{"from_length": 1}, {"from_length": 1, "to_length": 1}]})";
    
    Mapping m;
    json2pb(m, mapping_string.c_str(), mapping_string.size());
    
    REQUIRE_THROWS(simplify(m));
    
}

TEST_CASE("Mapping adjacent edit merging keeps leading deletions", "[mapping]") {
    
    string mapping_string = R"({"edit": [{"from_length": 1}, {"from_length": 1, "to_length": 1}]})";
    
    Mapping m;
    json2pb(m, mapping_string.c_str(), mapping_string.size());
    
    auto merged = merge_adjacent_edits(m);
    
    REQUIRE(merged.edit_size() == 2);
    
}

TEST_CASE("Mapping adjacent edit merging keeps trailing deletions", "[mapping]") {
    
    string mapping_string = R"({"edit": [{"from_length": 1, "to_length": 1}, {"from_length": 1}]})";
    
    Mapping m;
    json2pb(m, mapping_string.c_str(), mapping_string.size());
    
    auto merged = merge_adjacent_edits(m);
    
    REQUIRE(merged.edit_size() == 2);
    
}
   
}
}
        
