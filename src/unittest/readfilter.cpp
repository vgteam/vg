/**
 * unittest/readfilter.cpp: test cases for the read filtering/transforming logic
 */

#include "catch.hpp"
#include "readfilter.hpp"
#include "xg.hpp"

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
            {"from": 5, "to": 7},
            {"from": 6, "to": 7},
            {"from": 7, "to": 8, "to_end": true},
            {"from": 7, "to": 9}       
        ]
    }
    
    )";
    
    // Load it into Protobuf
    Graph chunk;
    json2pb(chunk, graph_json.c_str(), graph_json.size());
    
    // Pass it over to XG
    xg::XG index;
    index.from_path_handle_graph(VG(chunk));
    
    // Make a ReadFilter;
    ReadFilter<Alignment> filter;
    filter.graph = &index;
    
    SECTION("A read that is not ambiguous is not trimmed") {
        // Build the read
        const string read_json = R"(
        
        {
            "sequence": "GATTACAAAAAA",
            "quality": "MDAwMDAwMDAwMDAw",
            "path": {
                "mapping": [
                    {"position": {"node_id": 1}, "edit": [{"from_length": 6, "to_length": 6}]},
                    {"position": {"node_id": 2}, "edit": [{"from_length": 1, "to_length": 1}]},
                    {"position": {"node_id": 3}, "edit": [{"from_length": 5, "to_length": 5}]}
                ]
            }
        }
        
        )";
        
        // Fluff it up into a real read
        Alignment ambiguous;
        json2pb(ambiguous, read_json.c_str(), read_json.size());
        
        // It has to say it didn't trim it
        REQUIRE(filter.trim_ambiguous_ends(ambiguous, 10) == false);
        // And it has to not mess it up
        REQUIRE(ambiguous.sequence() == "GATTACAAAAAA");
        REQUIRE(ambiguous.quality() == "000000000000");
        REQUIRE(ambiguous.path().mapping_size() == 3);
        REQUIRE(ambiguous.path().mapping(2).position().node_id() == 3);
        REQUIRE(ambiguous.path().mapping(2).edit_size() == 1);
        REQUIRE(ambiguous.path().mapping(2).edit(0).from_length() == 5);
        REQUIRE(ambiguous.path().mapping(2).edit(0).to_length() == 5);
        
        
    }
    
    SECTION("A read with an ambiguous end should be trimmed to the unambiguous part") {
        // Build the read
        const string read_json = R"(
        
        {
            "sequence": "GATTACAA",
            "quality": "MTIzNDU2Nzg=",
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
        
        // It has to say it trimmed it
        REQUIRE(filter.trim_ambiguous_ends(ambiguous, 10) == true);
        // And it has to trim it to the right thing
        REQUIRE(ambiguous.sequence() == "GATTAC");
        REQUIRE(ambiguous.quality() == "123456");
        REQUIRE(ambiguous.path().mapping_size() == 1);
        REQUIRE(ambiguous.path().mapping(0).position().node_id() == 1);
        REQUIRE(ambiguous.path().mapping(0).edit_size() == 1);
        REQUIRE(ambiguous.path().mapping(0).edit(0).from_length() == 6);
        REQUIRE(ambiguous.path().mapping(0).edit(0).to_length() == 6);
        
        
    }
    
    SECTION("A read with an ambiguous start should be trimmed to the unambiguous part") {
        // Build the read
        const string read_json = R"(
        
        {
            "sequence": "TTGTAATC",
            "quality": "MTIzNDU2Nzg=",
            "path": {
                "mapping": [
                    {"position": {"node_id": 3, "is_reverse": true, "offset": 4}, "edit": [{"from_length": 1, "to_length": 1}]},
                    {"position": {"node_id": 2, "is_reverse": true}, "edit": [{"from_length": 1, "to_length": 1}]},
                    {"position": {"node_id": 1, "is_reverse": true}, "edit": [{"from_length": 6, "to_length": 6}]}
                ]
            }
        }
        
        )";
        
        // Fluff it up into a real read
        Alignment ambiguous;
        json2pb(ambiguous, read_json.c_str(), read_json.size());
        
        // It has to say it trimmed it
        REQUIRE(filter.trim_ambiguous_ends(ambiguous, 10) == true);
        // And it has to trim it to the right thing
        REQUIRE(ambiguous.sequence() == "GTAATC");
        REQUIRE(ambiguous.quality() == "345678");
        REQUIRE(ambiguous.path().mapping_size() == 1);
        REQUIRE(ambiguous.path().mapping(0).position().node_id() == 1);
        REQUIRE(ambiguous.path().mapping(0).position().is_reverse() == true);
        REQUIRE(ambiguous.path().mapping(0).edit_size() == 1);
        REQUIRE(ambiguous.path().mapping(0).edit(0).from_length() == 6);
        REQUIRE(ambiguous.path().mapping(0).edit(0).to_length() == 6);
        
        
    }
    
    SECTION("A read with both an ambiguous start and an ambiguous end should be trimmed to the unambiguous part") {
        // Build the read
        const string read_json = R"(
        
        {
            "sequence": "TAGTAGAGATAAAA",
            "quality": "QkxBQkxBMDAwMDBMT0w=",
            "path": {
                "mapping": [
                    {"position": {"node_id": 4, "offset": 3}, "edit": [{"from_length": 3, "to_length": 3}]},
                    {"position": {"node_id": 6}, "edit": [{"from_length": 3, "to_length": 3}]},
                    {"position": {"node_id": 7}, "edit": [{"from_length": 5, "to_length": 5}]},
                    {"position": {"node_id": 8, "is_reverse": true}, "edit": [{"from_length": 3, "to_length": 3}]}
                ]
            }
        }
        
        )";
        
        // Fluff it up into a real read
        Alignment ambiguous;
        json2pb(ambiguous, read_json.c_str(), read_json.size());
        
        // It has to say it trimmed it
        REQUIRE(filter.trim_ambiguous_ends(ambiguous, 10) == true);
        // And it has to trim it to the right thing
        REQUIRE(ambiguous.sequence() == "AGATA");
        REQUIRE(ambiguous.quality() == "00000");
        REQUIRE(ambiguous.path().mapping_size() == 1);
        REQUIRE(ambiguous.path().mapping(0).position().node_id() == 7);
        REQUIRE(ambiguous.path().mapping(0).position().is_reverse() == false);
        REQUIRE(ambiguous.path().mapping(0).edit_size() == 1);
        REQUIRE(ambiguous.path().mapping(0).edit(0).from_length() == 5);
        REQUIRE(ambiguous.path().mapping(0).edit(0).to_length() == 5);
        
        
    }
    
    SECTION("A read with no edits in its mappings is rejected") {
        // Build the read
        const string read_json = R"(
        
        {
            "sequence": "TAGTAGAGATAAAA",
            "quality": "QkxBQkxBMDAwMDBMT0w=",
            "path": {
                "mapping": [
                    {"position": {"node_id": 4, "offset": 3}},
                    {"position": {"node_id": 6}},
                    {"position": {"node_id": 7}},
                    {"position": {"node_id": 8, "is_reverse": true}}
                ]
            }
        }
        
        )";
        
        // Fluff it up into a real read
        Alignment ambiguous;
        json2pb(ambiguous, read_json.c_str(), read_json.size());
        
        // It has to explode instead of trimming it. Mappings are always
        // required to be specified now.
        REQUIRE_THROWS(filter.trim_ambiguous_ends(ambiguous, 10));
        
    }
    
    SECTION("A read with differences from the reference should be trimmed to the unambiguous part") {
        // Build the read
        const string read_json = R"(
        
        {
            "sequence": "GATTACCAG",
            "quality": "MTIzNDU2Nzg5",
            "path": {
                "mapping": [
                    {"position": {"node_id": 1}, "edit": [{"from_length": 6, "to_length": 6}]},
                    {"position": {"node_id": 2}, "edit": [{"from_length": 1, "to_length": 1, "sequence": "C"}]},
                    {"position": {"node_id": 3}, "edit": [
                        {"from_length": 1, "to_length": 1}, 
                        {"from_length": 1, "to_length": 1, "sequence": "G"}
                    ]}
                ]
            }
        }
        
        )";
        
        // Fluff it up into a real read
        Alignment ambiguous;
        json2pb(ambiguous, read_json.c_str(), read_json.size());
        
        // It has to say it trimmed it
        REQUIRE(filter.trim_ambiguous_ends(ambiguous, 10) == true);
        // And it has to trim it to the right thing
        REQUIRE(ambiguous.sequence() == "GATTAC");
        REQUIRE(ambiguous.quality() == "123456");
        REQUIRE(ambiguous.path().mapping_size() == 1);
        REQUIRE(ambiguous.path().mapping(0).position().node_id() == 1);
        REQUIRE(ambiguous.path().mapping(0).edit_size() == 1);
        REQUIRE(ambiguous.path().mapping(0).edit(0).from_length() == 6);
        REQUIRE(ambiguous.path().mapping(0).edit(0).to_length() == 6);
        
        
    }

}


}
}
