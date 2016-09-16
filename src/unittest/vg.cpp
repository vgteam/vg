/**
 * unittest/vg.cpp: test cases for vg::VG methods
 */

#include "catch.hpp"
#include "vg.hpp"

namespace vg {
namespace unittest {

using namespace std;

// Turn a JSON string into a VG graph
VG string_to_graph(const string& json) {
    VG graph;
    Graph chunk;
    json2pb(chunk, json.c_str(), json.size());
    graph.merge(chunk);
    
    return graph;
}

TEST_CASE("is_acyclic should return whether the graph is acyclic", "[vg][cycles]") {
    
    SECTION("a tiny DAG should be acyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == true);
    }
    
    SECTION("a tiny cyclic graph should be cyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 2, "to": 1}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == false);
    }
    
    SECTION("a tiny cyclic graph using from_start and to_end should be cyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 2, "from_start": true, "to_end": true}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == false);
    }
    
    SECTION("a tiny cyclic graph using from_start and to_end the other way should be cyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"}
            ],
            "edge": [
                {"from": 2, "to": 1},
                {"from": 2, "to": 1, "from_start": true, "to_end": true}
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == false);
    }
    
    SECTION("a nontrivial DAG should be acyclic") {
        const string graph_json = R"(
        
        {
            "node": [
                {"id": 1, "sequence": "G"},
                {"id": 2, "sequence": "A"},
                {"id": 3, "sequence": "T"},
                {"id": 4, "sequence": "GGG"},
                {"id": 5, "sequence": "T"},
                {"id": 6, "sequence": "A"},
                {"id": 7, "sequence": "C"},
                {"id": 8, "sequence": "A"},
                {"id": 9, "sequence": "A"}
            ],
            "edge": [
                {"from": 1, "to": 2},
                {"from": 1, "to": 6},
                {"from": 2, "to": 3},
                {"from": 2, "to": 4},
                {"from": 3, "to": 5},
                {"from": 4, "to": 5},
                {"from": 5, "to": 6},
                {"from": 6, "to": 7},
                {"from": 6, "to": 8},
                {"from": 7, "to": 9},
                {"from": 8, "to": 9}
                
            ]
        }
    
        )";
        
        VG graph = string_to_graph(graph_json);
        
        REQUIRE(graph.is_acyclic() == true);
    }
}

}
}
