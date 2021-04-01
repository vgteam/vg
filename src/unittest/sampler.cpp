/// \file sampler.cpp
///  
/// unit tests for the Sampler

#include <iostream>
#include <unordered_set>
#include <utility>

#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../sampler.hpp"
#include "../xg.hpp"
#include "../vg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
    
TEST_CASE( "Sampler can sample from a 1-node graph", "[sampler]" ) {
    
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
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    // Define a sampler    
    Sampler sampler(&xg_index, 1337);
    
    SECTION( "Can sample all bases in both directions" ) {
        
        unordered_set<pair<size_t, bool>> seen;
        
        for (size_t i = 0; i < 100; i++) {
            // Sample a bunch of alignments of 1 base
            Alignment aln = sampler.alignment(1);
            
            // Make sure we got 1 base
            REQUIRE(aln.sequence().size() == 1);
            
            // And that it's on the right node
            REQUIRE(aln.path().mapping(0).position().node_id() == 1);
            
            // Remember what offset and orientation we got.
            seen.emplace(aln.path().mapping(0).position().offset(), aln.path().mapping(0).position().is_reverse());
        }
        
        // We need to see all 7 bases in both orientations.
        REQUIRE(seen.size() == 7 * 2);
        
    }
    
    SECTION( "Can sample all bases in both directions from a path" ) {
        
        // Same as above except we do this
        sampler.set_source_paths({"ref"}, {}, {}, {});
        
        unordered_set<pair<size_t, bool>> seen;
        
        for (size_t i = 0; i < 100; i++) {
            // Sample a bunch of alignments of 1 base
            Alignment aln = sampler.alignment(1);
            
            // Make sure we got 1 base
            REQUIRE(aln.sequence().size() == 1);
            
            // And that it's on the right node
            REQUIRE(aln.path().mapping(0).position().node_id() == 1);
            
            // Remember what offset and orientation we got.
            seen.emplace(aln.path().mapping(0).position().offset(), aln.path().mapping(0).position().is_reverse());
        }
        
        // We need to see all 7 bases in both orientations.
        REQUIRE(seen.size() == 7 * 2);
        
        
    }
    
}

TEST_CASE( "position_at works", "[sampler]" ) {

    string graph_json = R"({
        "node": [
            {"id": 1, "sequence": "GA"},
            {"id": 2, "sequence": "ACC"},
            {"id": 3, "sequence": "ACA"}
        ],
        "edge": [
            {"from": 1, "to": 2, "to_end": true},
            {"from": 2, "to": 2},
            {"from": 2, "to": 3, "from_start": true}
        ],
        "path": [
            {"name": "ref", "mapping": [
                {"rank": 1, "position": {"node_id": 1}, "edit": [{"from_length": 2, "to_length": 2}]},
                {"rank": 2, "position": {"node_id": 2, "is_reverse": true}, "edit": [{"from_length": 3, "to_length": 3}]},
                {"rank": 3, "position": {"node_id": 2, "is_reverse": true}, "edit": [{"from_length": 3, "to_length": 3}]},
                {"rank": 4, "position": {"node_id": 3}, "edit": [{"from_length": 3, "to_length": 3}]}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    SECTION( "position_at works on the forward strand of a forward node" ) {
        REQUIRE(position_at(&xg_index, "ref", 0, false) == make_pos_t(1, false, 0));
        REQUIRE(position_at(&xg_index, "ref", 1, false) == make_pos_t(1, false, 1));
        
        REQUIRE(position_at(&xg_index, "ref", 8, false) == make_pos_t(3, false, 0));
        REQUIRE(position_at(&xg_index, "ref", 9, false) == make_pos_t(3, false, 1));
        REQUIRE(position_at(&xg_index, "ref", 10, false) == make_pos_t(3, false, 2));
    }
    
    SECTION( "position_at works on the reverse strand of a forward node" ) {
        REQUIRE(position_at(&xg_index, "ref", 0, true) == make_pos_t(1, true, 1));
        REQUIRE(position_at(&xg_index, "ref", 1, true) == make_pos_t(1, true, 0));
        
        REQUIRE(position_at(&xg_index, "ref", 8, true) == make_pos_t(3, true, 2));
        REQUIRE(position_at(&xg_index, "ref", 9, true) == make_pos_t(3, true, 1));
        REQUIRE(position_at(&xg_index, "ref", 10, true) == make_pos_t(3, true, 0));
    }
    
    SECTION( "position_at works on the reverse strand of a reverse node" ) {
        REQUIRE(position_at(&xg_index, "ref", 2, false) == make_pos_t(2, true, 0));
        REQUIRE(position_at(&xg_index, "ref", 3, false) == make_pos_t(2, true, 1));
        REQUIRE(position_at(&xg_index, "ref", 4, false) == make_pos_t(2, true, 2));
        
        REQUIRE(position_at(&xg_index, "ref", 5, false) == make_pos_t(2, true, 0));
        REQUIRE(position_at(&xg_index, "ref", 6, false) == make_pos_t(2, true, 1));
        REQUIRE(position_at(&xg_index, "ref", 7, false) == make_pos_t(2, true, 2));
    }
    
    SECTION( "position_at works on the forward strand of a reverse node" ) {
        REQUIRE(position_at(&xg_index, "ref", 2, true) == make_pos_t(2, false, 2));
        REQUIRE(position_at(&xg_index, "ref", 3, true) == make_pos_t(2, false, 1));
        REQUIRE(position_at(&xg_index, "ref", 4, true) == make_pos_t(2, false, 0));
        
        REQUIRE(position_at(&xg_index, "ref", 5, true) == make_pos_t(2, false, 2));
        REQUIRE(position_at(&xg_index, "ref", 6, true) == make_pos_t(2, false, 1));
        REQUIRE(position_at(&xg_index, "ref", 7, true) == make_pos_t(2, false, 0));
    }



}

TEST_CASE( "Sampler can sample from a loop-containing path", "[sampler]" ) {
    
    string graph_json = R"({
        "node": [
            {"id": 1, "sequence": "GA"},
            {"id": 2, "sequence": "T"},
            {"id": 3, "sequence": "ACA"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 2},
            {"from": 2, "to": 3}
        ],
        "path": [
            {"name": "ref", "mapping": [
                {"rank": 1, "position": {"node_id": 1}, "edit": [{"from_length": 2, "to_length": 2}]},
                {"rank": 2, "position": {"node_id": 2}, "edit": [{"from_length": 1, "to_length": 1}]},
                {"rank": 3, "position": {"node_id": 2}, "edit": [{"from_length": 1, "to_length": 1}]},
                {"rank": 4, "position": {"node_id": 3}, "edit": [{"from_length": 3, "to_length": 3}]}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    // Define a sampler    
    Sampler sampler(&xg_index, 1337);
    
    SECTION( "Can sample entire path" ) {
        
        // Same as above except we do this
        sampler.set_source_paths({"ref"}, {}, {}, {});
        
        unordered_set<string> found;
        
        for (size_t i = 0; i < 100; i++) {
            // Sample a bunch of alignments of up to 7 bases
            Alignment aln = sampler.alignment(7);
            found.insert(aln.sequence());
        }
        
        // We need to see the proper ref path all the way through.
        REQUIRE(found.count("GATTACA") == 1);
        // In both orientations
        REQUIRE(found.count("TGTAATC") == 1);
        
        // And that we don't have the version that skips the loop
        REQUIRE(found.count("GATACA") == 0);
        REQUIRE(found.count("TGTATC") == 0);
        
        
    }
}

TEST_CASE( "Sampler can across reversing edges", "[sampler]" ) {
    
    string graph_json = R"({
        "node": [
            {"id": 1, "sequence": "GA"},
            {"id": 2, "sequence": "ACC"},
            {"id": 3, "sequence": "ACA"}
        ],
        "edge": [
            {"from": 1, "to": 2, "to_end": true},
            {"from": 2, "to": 2},
            {"from": 2, "to": 3, "from_start": true}
        ],
        "path": [
            {"name": "ref", "mapping": [
                {"rank": 1, "position": {"node_id": 1}, "edit": [{"from_length": 2, "to_length": 2}]},
                {"rank": 2, "position": {"node_id": 2, "is_reverse": true}, "edit": [{"from_length": 3, "to_length": 3}]},
                {"rank": 3, "position": {"node_id": 2, "is_reverse": true}, "edit": [{"from_length": 3, "to_length": 3}]},
                {"rank": 4, "position": {"node_id": 3}, "edit": [{"from_length": 3, "to_length": 3}]}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    // Define a sampler    
    Sampler sampler(&xg_index, 1337);
    
    SECTION( "Can sample entire path" ) {
        
        // Same as above except we do this
        sampler.set_source_paths({"ref"}, {}, {}, {});
        
        unordered_set<string> found;
        
        for (size_t i = 0; i < 100; i++) {
            // Sample a bunch of alignments
            Alignment aln = sampler.alignment(11);
            found.insert(aln.sequence());
        }
        
        // We need to see the proper ref path all the way through.
        REQUIRE(found.count("GAGGTGGTACA") == 1);
        // In both orientations
        REQUIRE(found.count("TGTACCACCTC") == 1);
        
        // And that we don't have the version that skips the loop
        REQUIRE(found.count("GATGGACA") == 0);
        REQUIRE(found.count("TGTACCTC") == 0);
        
        
    }
}

}

}

    
