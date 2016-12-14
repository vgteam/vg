/** \file
 *
 * Unit tests for the PhaseDuplicator, which replaces complex subgraphs with
 * observed traversals.
 */

#include <iostream>
#include "../phase_duplicator.hpp"
#include "../json2pb.h"

#include "catch.hpp"

namespace vg {
namespace unittest {

// Build a toy graph
const string duplicator_graph_1 = R"(
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
            
        ],
        "path": [
            {"name": "hint", "mapping": [
                {"position": {"node_id": 1}},
                {"position": {"node_id": 6}},
                {"position": {"node_id": 8}},
                {"position": {"node_id": 9}}
            ]}
        ]
    }
    )";

TEST_CASE("PhaseDuplicator can emit a single traversal", "[phaseduplicator][indexing]") {
        
    // Make a graph
    
    // Make an actual graph
    Graph graph;
    json2pb(graph, duplicator_graph_1.c_str(), duplicator_graph_1.size());
            
    // Make an XG index of the graph
    xg::XG index(graph);
    
    // Add in a thread
    vector<xg::XG::thread_t> threads {{
        {1, false},
        {2, false},
        {3, false},
        {5, false},
        {6, false},
        {7, false},
        {9, false}
    }};
    index.insert_threads_into_dag(threads);
    
    // Make a duplicator
    PhaseDuplicator duplicator(index);
    
    // Define a subgraph to duplicate out
    set<id_t> to_duplicate {2, 3, 4, 5};
    // And the next ID to use
    id_t next_id = 10;
    
    // Do the duplication
    pair<Graph, vector<Translation>> result(duplicator.duplicate(to_duplicate, next_id));
    Graph& duped = result.first;
    vector<Translation>& translations = result.second;
    
#ifdef debug
    cerr << pb2json(duped) << endl;
#endif
    
    // Check the result
    SECTION("the duplicated graph should have 3 nodes") {
        REQUIRE(duped.node_size() == 3);
        
        SECTION("the nodes should be A, T, and T") {
            bool a_found = false;
            for (size_t i = 0; i < duped.node_size(); i++) {
                auto& node = duped.node(i);
                if (node.sequence() == "A") {
                    REQUIRE(!a_found);
                    a_found = true;
                } else {
                    REQUIRE(node.sequence() == "T");
                }
            }
            REQUIRE(a_found);
        }
    }
    
    SECTION("the duplicated graph should have 4 edges") {
        REQUIRE(duped.edge_size() == 4);
    }
    
    SECTION("there should be 1 translation") {
        REQUIRE(translations.size() == 1);
        
        auto& translation = translations.front();
        
        SECTION("the translation should have 3 mappings in each path") {
            REQUIRE(translation.from().mapping_size() == 3);
            REQUIRE(translation.to().mapping_size() == 3);
            
            SECTION("the from mappings should trace the original path") {
                REQUIRE(translation.from().mapping(0).position().node_id() == 2);
                REQUIRE(translation.from().mapping(0).position().is_reverse() == false);
                REQUIRE(translation.from().mapping(1).position().node_id() == 3);
                REQUIRE(translation.from().mapping(1).position().is_reverse() == false);
                REQUIRE(translation.from().mapping(2).position().node_id() == 5);
                REQUIRE(translation.from().mapping(2).position().is_reverse() == false);
            }
            
            SECTION("the to mappings should also spell out ATT") {
                // We need to index the duplicated nodes
                map<id_t, Node> id_to_node;
                for (size_t i = 0; i < duped.node_size(); i++) {
                    id_to_node[duped.node(i).id()] = duped.node(i);
                }
                
                for (size_t i = 0; i < translation.to().mapping_size(); i++) {
                    // All the mappings have to be to forward copies of these nodes.
                    auto& mapping = translation.to().mapping(i);
                    REQUIRE(id_to_node.count(mapping.position().node_id()));
                    REQUIRE(mapping.position().is_reverse() == false);
                }
                
                // And when we look at the nodes they need the right sequences
                REQUIRE(id_to_node[translation.to().mapping(0).position().node_id()].sequence() == "A");
                REQUIRE(id_to_node[translation.to().mapping(1).position().node_id()].sequence() == "T");
                REQUIRE(id_to_node[translation.to().mapping(2).position().node_id()].sequence() == "T");
            }
        }
        
    }
    
            
}

TEST_CASE("PhaseDuplicator can deduplicate identical traversals", "[phaseduplicator][indexing]") {
        
    // Make a graph
    
    // Make an actual graph
    Graph graph;
    json2pb(graph, duplicator_graph_1.c_str(), duplicator_graph_1.size());
            
    // Make an XG index of the graph
    xg::XG index(graph);
    
    // Add in a thread
    vector<xg::XG::thread_t> threads {{
        {1, false},
        {2, false},
        {3, false},
        {5, false},
        {6, false},
        {7, false},
        {9, false}
    }, {
        {1, false},
        {2, false},
        {3, false},
        {5, false},
        {6, false},
        {7, false},
        {9, false}
    }, {
        {9, true},
        {7, true},
        {6, true},
        {5, true},
        {3, true},
        {2, true},
        {1, true}
    }};
    index.insert_threads_into_dag(threads);
    
    // Make a duplicator
    PhaseDuplicator duplicator(index);
    
    // Define a subgraph to duplicate out
    set<id_t> to_duplicate {2, 3, 4, 5};
    // And the next ID to use
    id_t next_id = 10;
    
    // Do the duplication
    pair<Graph, vector<Translation>> result(duplicator.duplicate(to_duplicate, next_id));
    Graph& duped = result.first;
    vector<Translation>& translations = result.second;
    
#ifdef debug
    cerr << pb2json(duped) << endl;
#endif
    
    // Check the result
    SECTION("the duplicated graph should have 3 nodes") {
        REQUIRE(duped.node_size() == 3);
        
        SECTION("the nodes should be A, T, and T") {
            bool a_found = false;
            for (size_t i = 0; i < duped.node_size(); i++) {
                auto& node = duped.node(i);
                if (node.sequence() == "A") {
                    REQUIRE(!a_found);
                    a_found = true;
                } else {
                    REQUIRE(node.sequence() == "T");
                }
            }
            REQUIRE(a_found);
        }
    }
    
    SECTION("the duplicated graph should have 4 edges") {
        REQUIRE(duped.edge_size() == 4);
    }
}

TEST_CASE("PhaseDuplicator can produce internal traversals", "[phaseduplicator][indexing]") {
        
    // Make a graph
    
    // Make an actual graph
    Graph graph;
    json2pb(graph, duplicator_graph_1.c_str(), duplicator_graph_1.size());
            
    // Make an XG index of the graph
    xg::XG index(graph);
    
    // Add in a thread
    vector<xg::XG::thread_t> threads {{
        {1, false},
        {2, false},
        {3, false},
        {5, false},
        {6, false},
        {7, false},
        {9, false}
    }, {
        {1, false},
        {2, false},
        {3, false},
        {5, false},
        {6, false},
        {7, false},
        {9, false}
    }, {
        {4, true}
    }};
    index.insert_threads_into_dag(threads);
    
    // Make a duplicator
    PhaseDuplicator duplicator(index);
    
    // Define a subgraph to duplicate out
    set<id_t> to_duplicate {2, 3, 4, 5};
    // And the next ID to use
    id_t next_id = 10;
    
    // Do the duplication
    pair<Graph, vector<Translation>> result(duplicator.duplicate(to_duplicate, next_id));
    Graph& duped = result.first;
    vector<Translation>& translations = result.second;
    
#ifdef debug
    cerr << pb2json(duped) << endl;
#endif
    
    // Check the result
    SECTION("the duplicated graph should have 4 nodes") {
        REQUIRE(duped.node_size() == 4);
        
        SECTION("the nodes should be A, T, T, and GGG") {
            bool a_found = false;
            bool ggg_found = false;
            for (size_t i = 0; i < duped.node_size(); i++) {
                auto& node = duped.node(i);
                if (node.sequence() == "A") {
                    REQUIRE(!a_found);
                    a_found = true;
                } else if (node.sequence() == "GGG") {
                    REQUIRE(!ggg_found);
                    ggg_found = true;
                } else {
                    REQUIRE(node.sequence() == "T");
                }
            }
            REQUIRE(a_found);
            REQUIRE(ggg_found);
        }
    }
    
    SECTION("the duplicated graph should have 4 edges") {
        REQUIRE(duped.edge_size() == 4);
    }
}
    
}
}
