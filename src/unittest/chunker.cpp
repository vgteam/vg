/**
 * unittest/chunker.cpp: test cases for the read filtering/transforming logic
 */

#include "catch.hpp"
#include "chunker.hpp"
#include "readfilter.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("basic graph chunking", "[chunk]") {
    
    // Build a toy graph (copied from readfilter test)
    const string graph_json = R"(
    
    {
        "node": [
            {"id": 1, "sequence": "GATTAC"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "AAAAA"},
            {"id": 4, "sequence": "CATTAG"},
            {"id": 5, "sequence": "TAGTAG"},
            {"id": 6, "sequence": "TAG"},
            {"id": 7, "sequence": "AGATA"},
            {"id": 8, "sequence": "TTT"},
            {"id": 9, "sequence": "AAA"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 2, "to": 3},
            {"from": 1, "to": 3},
            {"from": 3, "to": 4},
            {"from": 4, "to": 5},
            {"from": 6, "to": 4, "from_start": true, "to_end": true},
            {"from": 4, "to": 6, "from_start": true},
            {"from": 5, "to": 7},
            {"from": 6, "to": 7},
            {"from": 7, "to": 8, "to_end": true},
            {"from": 7, "to": 9},
            {"from": 7, "to": 5, "to_end": true}
        ],
        "path": [
            {"name": "x", "mapping": [
                 {"position": {"node_id": 1}, "rank": 1}, 
                 {"position": {"node_id": 2}, "rank": 2}, 
                 {"position": {"node_id": 3}, "rank": 3}, 
                 {"position": {"node_id": 4}, "rank": 4}, 
                 {"position": {"node_id": 5}, "rank": 5}, 
                 {"position": {"node_id": 7}, "rank": 6}, 
                 {"position": {"node_id": 9}, "rank": 7}]},
            {"name": "y", "mapping": [
                 {"position": {"node_id": 1}, "rank": 1}, 
                 {"position": {"node_id": 3}, "rank": 2}, 
                 {"position": {"node_id": 4}, "rank": 3}, 
                 {"position": {"node_id": 6}, "rank": 4}, 
                 {"position": {"node_id": 7}, "rank": 5}, 
                 {"position": {"node_id": 5, "is_reverse": true}, "rank": 6}, 
                 {"position": {"node_id": 4, "is_reverse": true}, "rank": 7}]},
            {"name": "z", "mapping": [
                 {"position": {"node_id": 1}, "rank": 1}, 
                 {"position": {"node_id": 3}, "rank": 2}, 
                 {"position": {"node_id": 4}, "rank": 3}, 
                 {"position": {"node_id": 6}, "rank": 4}, 
                 {"position": {"node_id": 7}, "rank": 5}, 
                 {"position": {"node_id": 5, "is_reverse": true}, "rank": 6}, 
                 {"position": {"node_id": 4, "is_reverse": true}, "rank": 7},
                 {"position": {"node_id": 6}, "rank": 8},
                 {"position": {"node_id": 7}, "rank": 9}, 
                 {"position": {"node_id": 5, "is_reverse": true}, "rank": 10}, 
                 {"position": {"node_id": 4, "is_reverse": true}, "rank": 11},
                 {"position": {"node_id": 6}, "rank": 12},
                 {"position": {"node_id": 7}, "rank": 13}, 
                 {"position": {"node_id": 5, "is_reverse": true}, "rank": 14}, 
                 {"position": {"node_id": 4, "is_reverse": true}, "rank": 15}]}

        ] 
    }
    
    )";
    
    // Load it into Protobuf
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    
    // Pass it over to XG
    xg::XG index(chunk);

    PathChunker chunker(&index);

    SECTION("Extract whole graph as chunk") {

        // Note: regions are 1-based inclusive
        Region region = {"x", 1, 32};
        VG subgraph;
        Region out_region;
        chunker.extract_subgraph(region, 1, false, subgraph, out_region);

        REQUIRE(subgraph.node_count() == 9);
        REQUIRE(subgraph.edge_count() == 12);
        REQUIRE(out_region.start == 1);
    }

    SECTION("Extract partial graph as chunk") {
        
        Region region = {"x", 9, 16};
        VG subgraph;
        Region out_region;        
        chunker.extract_subgraph(region, 1, false, subgraph, out_region);

        REQUIRE(subgraph.node_count() == 6);
        REQUIRE(subgraph.edge_count() == 7);
        REQUIRE(out_region.start == 1);
    }

    SECTION("Extract partial graph as chunk via id range") {
        
        Region region = {"x", 3, 6};
        VG subgraph;
        Region out_region;        
        chunker.extract_id_range(region.start, region.end, 1, false, subgraph, out_region);
        
        REQUIRE(subgraph.node_count() == 7);
        REQUIRE(subgraph.edge_count() == 10);
        REQUIRE(out_region.start == 1);
    }

    SECTION("Extract partial graph as chunk with exact node boundary") {

        Region region = {"x", 8, 16};
        VG subgraph;
        Region out_region;        
        chunker.extract_subgraph(region, 1, false, subgraph, out_region);

        REQUIRE(subgraph.node_count() == 6);
        REQUIRE(subgraph.edge_count() == 7);
        REQUIRE(out_region.start == 1);
    }

    SECTION("Extract whole graph via harder path") {

        Region region = {"y", 1, 37};
        VG subgraph;
        Region out_region;
        chunker.extract_subgraph(region, 1, false, subgraph, out_region);

        REQUIRE(subgraph.node_count() == 9);
        REQUIRE(subgraph.edge_count() == 12);
        REQUIRE(out_region.start == 1);

    }

    SECTION("Extract partial graph via harder path") {

        Region region = {"y", 15, 35};
        VG subgraph;
        Region out_region;
        chunker.extract_subgraph(region, 1, false, subgraph, out_region);

        REQUIRE(subgraph.node_count() == 7);
        REQUIRE(subgraph.edge_count() == 9);
        REQUIRE(out_region.start == 7);
    }
    
    SECTION("Extract whole graph via cyclic path") {
        
        Region region = {"z", 1, 60};
        VG subgraph;
        Region out_region;        
        chunker.extract_subgraph(region, 1, false, subgraph, out_region);
        
        REQUIRE(subgraph.node_count() == 9);
        REQUIRE(subgraph.edge_count() == 12);
        REQUIRE(out_region.start == 1);
        
    }
    
    SECTION("Partial graph via cyclic path") {
        
        Region region = {"z", 36, 59};
        VG subgraph;
        Region out_region;        
        chunker.extract_subgraph(region, 1, false, subgraph, out_region);
        
        REQUIRE(subgraph.node_count() == 7);
        REQUIRE(subgraph.edge_count() == 9);
        REQUIRE(out_region.start == 7);
        
    }

}


}
}
