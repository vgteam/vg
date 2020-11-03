/**
 * unittest/genotypekit.cpp: test cases for genotypekit modular genotyper pieces
 */

#include "catch.hpp"
#include "../genotyper.hpp"
#include "../snarls.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../traversal_finder.hpp"

namespace vg {
namespace unittest {

TEST_CASE("traversals can be found from reads", "[genotyper]") {
    
    // Build a toy graph
    // We will have snarls 1 to 6 and 6 to 9
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
            
        ],
        "path": [
            {"name": "hint", "mapping": [
                {"position": {"node_id": 1}, "rank" : 1 },
                {"position": {"node_id": 6}, "rank" : 2 },
                {"position": {"node_id": 8}, "rank" : 3 },
                {"position": {"node_id": 9}, "rank" : 4 }
            ]}
        ]
    }
    
    )";
    
    // Make an actual graph
    VG graph;
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    graph.merge(chunk);
    
    // Find the snarls
    SnarlManager manager = CactusSnarlFinder(graph).find_snarls();
    
    // Pull out the 1 to 6 one and the 6 to 9 one.
    auto snarls = manager.top_level_snarls();
    REQUIRE(snarls.size() == 2);
    REQUIRE(snarls[0]->start().node_id() == 1);
    REQUIRE(snarls[0]->end().node_id() == 6);
    REQUIRE(snarls[1]->start().node_id() == 6);
    REQUIRE(snarls[1]->end().node_id() == 9);
    
    // Make a Genotyper
    Genotyper genotyper;
    // Turn off recurrence checking
    genotyper.min_recurrence = 0;
    
    SECTION("Forward reads produce paths") {
    
        // Make an aligned read
        Alignment fwd = graph.align("GATTACA");
        
        // Make sure it is in the graph for lookup
        fwd.mutable_path()->set_name("read1");
        graph.paths.extend(fwd.path());

        // Make a map to hold it
        map<string, const Alignment*> reads_by_name{{"read1", &fwd}};

        // Get paths through the 1 to 6 snarl supported by the read's alignment
        PathRestrictedTraversalFinder path_trav_finder(graph, manager, reads_by_name, 0);
        vector<SnarlTraversal> paths = path_trav_finder.find_traversals(*snarls[0]);
        
        // We should cross the snarl the easy forward way.
        REQUIRE(paths.size() == 1);
        REQUIRE(paths[0].visit_size() == 5);
        size_t i = 0;
        REQUIRE(paths[0].visit(i++).node_id() == 1);
        REQUIRE(paths[0].visit(i++).node_id() == 2);
        REQUIRE(paths[0].visit(i++).node_id() == 3);
        REQUIRE(paths[0].visit(i++).node_id() == 5);
        REQUIRE(paths[0].visit(i++).node_id() == 6);
        
    }
    
    SECTION("Backward reads produce paths") {
    
        // Make an aligned read
        Alignment rev = graph.align("GATTACA");
        
        // Make it be on the reverse strand of the graph
        rev = reverse_complement_alignment(rev, [&](id_t id) {
            return graph.get_node(id)->sequence().size();
        });
        
        // Make sure it is in the graph for lookup
        rev.mutable_path()->set_name("read1");
        graph.paths.extend(rev.path());

        // Make an augmented graph
        AugmentedGraph aug;
        aug.graph = graph;
        vector<Alignment> alignments = {rev};
        aug.augment_from_alignment_edits(alignments);

        // Get paths through the 1 to 6 snarl supported by the read's alignment
        ReadRestrictedTraversalFinder read_trav_finder(aug, manager, 0);
        vector<SnarlTraversal> paths = read_trav_finder.find_traversals(*snarls[0]);
        
        // We should cross the snarl the easy forward way.
        REQUIRE(paths.size() == 1);
        REQUIRE(paths[0].visit_size() == 5);
        size_t i = 0;
        REQUIRE(paths[0].visit(i++).node_id() == 1);
        REQUIRE(paths[0].visit(i++).node_id() == 2);
        REQUIRE(paths[0].visit(i++).node_id() == 3);
        REQUIRE(paths[0].visit(i++).node_id() == 5);
        REQUIRE(paths[0].visit(i++).node_id() == 6);
        
    }
    
    
    
}


}
}
