/// \file path.cpp
///  
/// unit tests for Paths and their utility functions
///

#include <iostream>
#include <string>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../path.hpp"
#include "../vg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Path simplification tolerates adjacent insertions and deletions", "[path]") {

    string path_string = R"(
        {
            "mapping": [
                {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "68"}},
                {"edit": [{"sequence": "AAGG", "to_length": 4}, {"from_length": 3}], "position": {"node_id": "67"}},
                {"edit": [{"from_length": 17, "to_length": 17}], "position": {"node_id": "66"}}
            ]
        }
    )";

    Path path;
    json2pb(path, path_string.c_str(), path_string.size());

    auto simple = simplify(path);
    
    std::cerr << pb2json(simple) << std::endl;

    // We need to still touch all the nodes after simplification.
    REQUIRE(simple.mapping_size() == 3);
    REQUIRE(simple.mapping(0).position().node_id() == 68);
    REQUIRE(simple.mapping(1).position().node_id() == 67);
    REQUIRE(simple.mapping(2).position().node_id() == 66);

}

}
}
