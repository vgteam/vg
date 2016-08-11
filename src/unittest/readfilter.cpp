/**
 * unittest/readfilter.cpp: test cases for the read filtering/transforming logic
 */

#include "catch.hpp"
#include "readfilter.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("reads with ambiguous ends can be trimmed", "[filter]") {
    
    // Build a toy graph
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "GATTAC"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "AAAAA"},
            {"id": 4, "sequence": "CATTAG"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 3},
            {"from": 1, "to": 3},
            {"from": 3, "to": 4}            
        ]
    }
    
    )";
    
    // Load it into Protobuf
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    
    // Pass it over to XG
    xg::XG index(chunk);
    
    // Make a ReadFilter;
    ReadFilter filter;
    
    // Configure it for end trimming
    filter.defray_length = 3;
    
    SECTION("A read with an ambiguous end should be trimmed") {
        // Build the read
        const string read_json = R"(
        
        {
            "sequence": "GATTACAA",
            "quality": "ICAgICAgICA=",
            "path": {
                "mapping": [
                    {"position": {"node_id": 1}, "edit": [{"from_length": 6, "to_length": 6}]},
                    {"position": {"node_id": 2}, "edit": [{"from_length": 1, "to_length": 1}]},
                    {"position": {"node_id": 3}, "edit": [{"from_length": 1, "to_length": 1}]}
                ]
            }
        }
        
        )";
        
        // Fluff it up into a real read
        Alignment ambiguous;
        json2pb(ambiguous, read_json.c_str(), read_json.size());
        
        
    }

}


}
}
