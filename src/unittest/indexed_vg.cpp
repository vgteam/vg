/// \file indexed_vg.cpp
///  
/// unit tests for the vg-file-backed handle graph implementation

#include <iostream>
#include "json2pb.h"
#include "vg.pb.h"
#include "../indexed_vg.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

/// Save the given Graph into a temporary vg file and return the file name
string save_vg_file(const Graph& graph) {
    string filename = temp_file::create();
    VG vg;
    vg.extend(graph);
    // Serialize with a small chunk size so we can exercise the indexing
    vg.serialize_to_file(filename, 10);
    return filename;
}

using namespace std;
    
TEST_CASE("An IndexedVG can be created for a single node", "[handle][indexed-vg]") {

    // This 1-node graph is sorted
    string graph_json = R"({
        "node": [{"id": 1, "sequence": "GATTACA"}],
        "path": [
            {"name": "ref", "mapping": [
                {"position": {"node_id": 1}, "edit": [{"from_length": 7, "to_length": 7}]}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());

    // Save a vg file
    string vg_name = save_vg_file(proto_graph);
    
    // Load the file up and build and save the index
    IndexedVG indexed(vg_name);
    
    for (auto& node : proto_graph.node()) {
        // For each node in the graph we saved
        
        // Get a handle
        handle_t handle = indexed.get_handle(node.id());
        
        // Make sure it has the right ID
        REQUIRE(indexed.get_id(handle) == node.id());
        
        // Make sure it has the right initial orientation
        REQUIRE(!indexed.get_is_reverse(handle));
        
        // Make sure flipping works
        handle_t flipped = indexed.flip(handle);
        REQUIRE(indexed.get_id(flipped) == node.id());
        REQUIRE(indexed.get_is_reverse(flipped));
        REQUIRE(indexed.flip(flipped) == handle);
        
        // Make sure the length is correct
        REQUIRE(indexed.get_length(handle) == node.sequence().size());
        
        // Make sure the sequence is correct
        REQUIRE(indexed.get_sequence(handle) == node.sequence());
        
        // Since this graph is empty, make sure there are no edges
        REQUIRE(indexed.get_degree(handle, false) == 0);
        REQUIRE(indexed.get_degree(handle, true) == 0);
    }
    
}

}

}
