/// \file indexed_vg.cpp
///  
/// unit tests for the vg-file-backed handle graph implementation

#include <iostream>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../indexed_vg.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../algorithms/id_sort.hpp"
#include "random_graph.hpp"
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
    
    REQUIRE(proto_graph.node_size() == 1);

    // Save a vg file
    string vg_name = save_vg_file(proto_graph);
    
    {
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
    
    // Clean up after the index
    temp_file::remove(vg_name);
    temp_file::remove(vg_name + ".vgi");
}

TEST_CASE("IndexedVG works on random graphs", "[handle][indexed-vg]") {
    for (size_t trial = 0; trial < 5; trial++) {
        // Make a bunch of random graphs
        VG random;
        random_graph(300, 3, 30, &random);
        
        // Sort each by ID
        random.id_sort();
        
        string filename = temp_file::create();
        random.serialize_to_file(filename, 10);
        
        {
        
            // Load the file up and build and save the index
            IndexedVG indexed(filename);
            
            unordered_set<handle_t> indexed_handles;
            
            random.for_each_handle([&](const handle_t& node) {
                // For each node in the graph we saved
                
                // Get a handle in the index
                handle_t handle = indexed.get_handle(random.get_id(node));
                
                // Remember it for later
                indexed_handles.insert(handle);
                
                // Make sure it has the right ID
                REQUIRE(indexed.get_id(handle) == random.get_id(node));
                
                // Make sure it has the right initial orientation
                REQUIRE(!indexed.get_is_reverse(handle));
                
                // Make sure flipping works
                handle_t flipped = indexed.flip(handle);
                REQUIRE(indexed.get_id(flipped) == random.get_id(node));
                REQUIRE(indexed.get_is_reverse(flipped));
                REQUIRE(indexed.flip(flipped) == handle);
                
                // Make sure the length is correct
                REQUIRE(indexed.get_length(handle) == random.get_length(node));
                
                // Make sure the sequence is correct
                REQUIRE(indexed.get_sequence(handle) == random.get_sequence(node));
                REQUIRE(indexed.get_sequence(indexed.flip(handle)) == random.get_sequence(random.flip(node)));
                
                // Make sure degrees agree
                REQUIRE(indexed.get_degree(handle, false) == random.get_degree(node, false));
                REQUIRE(indexed.get_degree(handle, true) == random.get_degree(node, true));
                
                // Make sure the edges match up under all conditions
                for (bool flip_it : {false, true}) {
                    for (bool go_left : {false, true}) {
                        unordered_set<handle_t> random_edge_partners;
                    
                        random.follow_edges(flip_it ? random.flip(node) : node, go_left, [&](const handle_t& other) {
                            random_edge_partners.insert(other);
                        });
                        
                        indexed.follow_edges(flip_it ? indexed.flip(handle) : handle, go_left, [&](const handle_t& other) {
                            // Get the corresponding random graph handle
                            handle_t random_handle = random.get_handle(indexed.get_id(other), indexed.get_is_reverse(other));
                            
                            REQUIRE(random_edge_partners.count(random_handle));
                            random_edge_partners.erase(random_handle);
                        });
                        
                        // Make sure we erased all the edges we added
                        REQUIRE(random_edge_partners.empty());
                    }
                }
                
                
                
            });
            
            // Make sure basic statistics agree
            REQUIRE(indexed.get_node_count() == random.get_node_count());
            REQUIRE(indexed.min_node_id() == random.min_node_id());
            REQUIRE(indexed.max_node_id() == random.max_node_id());
            
            // Make sure we can loop over the handles OK
            indexed.for_each_handle([&](const handle_t& handle) {
                REQUIRE(indexed_handles.count(handle));
                indexed_handles.erase(handle);
            });
            REQUIRE(indexed_handles.size() == 0);
            
        }
        
        // Clean up after the index
        temp_file::remove(filename);
        temp_file::remove(filename + ".vgi");
        
    }
}

}

}
