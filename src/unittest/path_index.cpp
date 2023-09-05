/// \file path_index.cpp
///  
/// Unit tests for the PathIndex class, which indexes paths for random access.
///

#include <iostream>
#include <string>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../path_index.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

// Build a toy graph
const string path_index_graph_1 = R"(
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
            {"id": 9, "sequence": "ACACA"}
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
            
        ],
        "path": [
            {"name": "cool", "mapping": [
                {"position": {"node_id": 1}},
                {"position": {"node_id": 2}},
                {"position": {"node_id": 4}},
                {"position": {"node_id": 5}},
                {"position": {"node_id": 6}},
                {"position": {"node_id": 8}},
                {"position": {"node_id": 9}}
            ]}
        ]
    }
    )";


TEST_CASE("PathIndex can be created", "[pathindex]") {
    
    // Load the graph
    Graph graph;
    json2pb(graph, path_index_graph_1.c_str(), path_index_graph_1.size());
    
    // Make it into a VG
    VG to_index;
    to_index.extend(graph);
    
    // Make a PathIndex
    PathIndex index(to_index, "cool", true);
    
    SECTION("PathIndex has the right string") {
        REQUIRE(index.sequence == "GAGGGTAAACACA");
    }
}

TEST_CASE("PathIndex translation can change a node ID", "[pathindex]") {

    // Load the graph
    Graph graph;
    json2pb(graph, path_index_graph_1.c_str(), path_index_graph_1.size());
    
    // Make it into a VG
    VG to_index;
    to_index.extend(graph);
    
    // Make a PathIndex
    PathIndex index(to_index, "cool", true);
    
    SECTION("Before translation, node 6 is at position 6") {
        REQUIRE(index.at_position(6).node == 6);
    }
    
    // Make a Translation to change node 6 to node 99
    Translation t;
    auto* from_mapping = t.mutable_from()->add_mapping();
    from_mapping->mutable_position()->set_node_id(6);
    auto* from_edit = from_mapping->add_edit();
    from_edit->set_from_length(1);
    from_edit->set_to_length(1);
    auto* to_mapping = t.mutable_to()->add_mapping();
    to_mapping->mutable_position()->set_node_id(99);
    auto* to_edit = to_mapping->add_edit();
    to_edit->set_from_length(1);
    to_edit->set_to_length(1);
    
    // Apply it
    index.apply_translation(t);
    
    SECTION("After translation, node 99 is at position 6") {
        REQUIRE(index.at_position(6).node == 99);
    }
    
}

TEST_CASE("PathIndex translation can divide a node", "[pathindex]") {
    
    // Load the graph
    Graph graph;
    json2pb(graph, path_index_graph_1.c_str(), path_index_graph_1.size());
    
    // Make it into a VG
    VG to_index;
    to_index.extend(graph);
    
    // Make a PathIndex
    PathIndex index(to_index, "cool", true);
    
    SECTION("Before translation, node 4 is at positions 2-4") {
        REQUIRE(index.at_position(2).node == 4);
        REQUIRE(index.at_position(3).node == 4);
        REQUIRE(index.at_position(4).node == 4);
    }
    
    // Make a Translation to change node 4 (GGG) to nodes 1337 (G) and 1338 (GG)
    Translation t;
    auto* from_mapping = t.mutable_from()->add_mapping();
    from_mapping->mutable_position()->set_node_id(4);
    auto* from_edit = from_mapping->add_edit();
    from_edit->set_from_length(3);
    from_edit->set_to_length(3);
    
    auto* to_mapping_1 = t.mutable_to()->add_mapping();
    to_mapping_1->mutable_position()->set_node_id(1337);
    auto* to_edit_1 = to_mapping_1->add_edit();
    to_edit_1->set_from_length(1);
    to_edit_1->set_to_length(1);
    
    auto* to_mapping_2 = t.mutable_to()->add_mapping();
    to_mapping_2->mutable_position()->set_node_id(1338);
    auto* to_edit_2 = to_mapping_2->add_edit();
    to_edit_2->set_from_length(2);
    to_edit_2->set_to_length(2);
    
    // Apply it
    index.apply_translation(t);
    
    SECTION("After translation, node 1337 is at position 2 forward") {
        REQUIRE(index.at_position(2).node == 1337);
        REQUIRE(index.at_position(2).is_end == false);
        REQUIRE(index.node_length(index.find_position(2)) == 1);
    }
    
    SECTION("After translation, node 1338 is at positions 3-4 forward") {
        REQUIRE(index.at_position(3).node == 1338);
        REQUIRE(index.at_position(3).is_end == false);
        REQUIRE(index.node_length(index.find_position(3)) == 2);
        REQUIRE(index.at_position(4).node == 1338);
        REQUIRE(index.at_position(4).is_end == false);
        REQUIRE(index.node_length(index.find_position(4)) == 2);
    }
    
}

TEST_CASE("PathIndex translation can create reverse strand mappings", "[pathindex]") {
    
    // Load the graph
    Graph graph;
    json2pb(graph, path_index_graph_1.c_str(), path_index_graph_1.size());
    
    // Make it into a VG
    VG to_index;
    to_index.extend(graph);
    
    // Make a PathIndex
    PathIndex index(to_index, "cool", true);
    
    SECTION("Before translation, node 4 is at positions 2-4") {
        REQUIRE(index.at_position(2).node == 4);
        REQUIRE(index.at_position(3).node == 4);
        REQUIRE(index.at_position(4).node == 4);
    }
    
    // Make a Translation to change node 4 (GGG) to nodes 1337 (G) and 1338 (GG)
    Translation t;
    auto* from_mapping = t.mutable_from()->add_mapping();
    from_mapping->mutable_position()->set_node_id(4);
    auto* from_edit = from_mapping->add_edit();
    from_edit->set_from_length(3);
    from_edit->set_to_length(3);
    
    auto* to_mapping_1 = t.mutable_to()->add_mapping();
    to_mapping_1->mutable_position()->set_node_id(1337);
    to_mapping_1->mutable_position()->set_is_reverse(true);
    auto* to_edit_1 = to_mapping_1->add_edit();
    to_edit_1->set_from_length(1);
    to_edit_1->set_to_length(1);
    
    auto* to_mapping_2 = t.mutable_to()->add_mapping();
    to_mapping_2->mutable_position()->set_node_id(1338);
    to_mapping_2->mutable_position()->set_is_reverse(true);
    auto* to_edit_2 = to_mapping_2->add_edit();
    to_edit_2->set_from_length(2);
    to_edit_2->set_to_length(2);
    
    // Apply it
    index.apply_translation(t);
    
    SECTION("After translation, node 1337 is at position 2 in reverse") {
        REQUIRE(index.at_position(2).node == 1337);
        REQUIRE(index.at_position(2).is_end == true);
        REQUIRE(index.node_length(index.find_position(2)) == 1);
    }
    
    SECTION("After translation, node 1338 is at positions 3-4 in reverse") {
        REQUIRE(index.at_position(3).node == 1338);
        REQUIRE(index.at_position(3).is_end == true);
        REQUIRE(index.node_length(index.find_position(3)) == 2);
        REQUIRE(index.at_position(4).node == 1338);
        REQUIRE(index.at_position(4).is_end == true);
        REQUIRE(index.node_length(index.find_position(4)) == 2);
    }
    
}

TEST_CASE("PathIndex translation can handle translations articulated for the reverse strand", "[pathindex]") {
    
    // Load the graph
    Graph graph;
    json2pb(graph, path_index_graph_1.c_str(), path_index_graph_1.size());
    
    // Make it into a VG
    VG to_index;
    to_index.extend(graph);
    
    // Make a PathIndex
    PathIndex index(to_index, "cool", true);
    
    SECTION("Before translation, node 4 is at positions 2-4") {
        REQUIRE(index.at_position(2).node == 4);
        REQUIRE(index.at_position(3).node == 4);
        REQUIRE(index.at_position(4).node == 4);
    }
    
    // Make a Translation to change node 4 (GGG) to nodes 1337 (G) and 1338 (GG)
    // But do it in reverse on both paths
    Translation t;
    auto* from_mapping = t.mutable_from()->add_mapping();
    from_mapping->mutable_position()->set_node_id(4);
    from_mapping->mutable_position()->set_is_reverse(true);
    auto* from_edit = from_mapping->add_edit();
    from_edit->set_from_length(3);
    from_edit->set_to_length(3);
    
    auto* to_mapping_1 = t.mutable_to()->add_mapping();
    to_mapping_1->mutable_position()->set_node_id(1338);
    to_mapping_1->mutable_position()->set_is_reverse(true);
    auto* to_edit_1 = to_mapping_1->add_edit();
    to_edit_1->set_from_length(2);
    to_edit_1->set_to_length(2);
    
    auto* to_mapping_2 = t.mutable_to()->add_mapping();
    to_mapping_2->mutable_position()->set_node_id(1337);
    to_mapping_2->mutable_position()->set_is_reverse(true);
    auto* to_edit_2 = to_mapping_2->add_edit();
    to_edit_2->set_from_length(1);
    to_edit_2->set_to_length(1);
    
    // Apply it
    index.apply_translation(t);
    
    // Everything should come out forward
    
    SECTION("After translation, node 1337 is at position 2 forward") {
        REQUIRE(index.at_position(2).node == 1337);
        REQUIRE(index.at_position(2).is_end == false);
        REQUIRE(index.node_length(index.find_position(2)) == 1);
    }
    
    SECTION("After translation, node 1338 is at positions 3-4 forward") {
        REQUIRE(index.at_position(3).node == 1338);
        REQUIRE(index.at_position(3).is_end == false);
        REQUIRE(index.node_length(index.find_position(3)) == 2);
        REQUIRE(index.at_position(4).node == 1338);
        REQUIRE(index.at_position(4).is_end == false);
        REQUIRE(index.node_length(index.find_position(4)) == 2);
    }
    
}

TEST_CASE("PathIndex translation can divide the last node", "[pathindex]") {
    
    // Load the graph
    Graph graph;
    json2pb(graph, path_index_graph_1.c_str(), path_index_graph_1.size());
    
    // Make it into a VG
    VG to_index;
    to_index.extend(graph);
    
    // Make a PathIndex
    PathIndex index(to_index, "cool", true);
    
    SECTION("Before translation, node 9 is at positions 8-12") {
        for (size_t i = 8; i < 13; i++) {
            REQUIRE(index.at_position(i).node == 9);
        }
    }
    
    // Make a Translation to change node 9 (ACACA) to nodes 20 (ACA) and 21 (CA)
    Translation t;
    auto* from_mapping = t.mutable_from()->add_mapping();
    from_mapping->mutable_position()->set_node_id(9);
    auto* from_edit = from_mapping->add_edit();
    from_edit->set_from_length(5);
    from_edit->set_to_length(5);
    
    auto* to_mapping_1 = t.mutable_to()->add_mapping();
    to_mapping_1->mutable_position()->set_node_id(20);
    auto* to_edit_1 = to_mapping_1->add_edit();
    to_edit_1->set_from_length(3);
    to_edit_1->set_to_length(3);
    
    auto* to_mapping_2 = t.mutable_to()->add_mapping();
    to_mapping_2->mutable_position()->set_node_id(21);
    auto* to_edit_2 = to_mapping_2->add_edit();
    to_edit_2->set_from_length(2);
    to_edit_2->set_to_length(2);
    
    // Apply it
    index.apply_translation(t);
    
    // Everything should come out forward
    
    SECTION("After translation, node 20 is at positions 8-10 forward") {
        for (size_t i = 8; i < 11; i++) {
            REQUIRE(index.at_position(i).node == 20);
            REQUIRE(index.at_position(i).is_end == false);
            REQUIRE(index.node_length(index.find_position(i)) == 3);
        }
    }
    
    SECTION("After translation, node 21 is at positions 11-12 forward") {
        for (size_t i = 11; i < 13; i++) {
            REQUIRE(index.at_position(i).node == 21);
            REQUIRE(index.at_position(i).is_end == false);
            REQUIRE(index.node_length(index.find_position(i)) == 2);
        }
    }
    
}
   
}
}
        
