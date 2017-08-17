/// \file multipath_mapper.cpp
///  
/// unit tests for the multipath mapper

#include <iostream>
#include "json2pb.h"
#include "vg.pb.h"
#include "../multipath_mapper.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
    
TEST_CASE( "MultipathMapper can map to a one-node graph", "[multipath][mapping][multipathmapper]" ) {
    
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
    
    // Configure GCSA temp directory to the system temp directory
    gcsa::TempFile::setDirectory(find_temp_dir());
    // And make it quiet
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Make pointers to fill in
    gcsa::GCSA* gcsaidx = nullptr;
    gcsa::LCPArray* lcpidx = nullptr;
    
    // Build the GCSA index
    graph.build_gcsa_lcp(gcsaidx, lcpidx, 16, false, false, 3);
    
    // Build the xg index
    xg::XG xg_index(proto_graph);
    
    // Make a multipath mapper to map against the graph.
    MultipathMapper mapper(&xg_index, gcsaidx, lcpidx);
    
    SECTION( "MultipathMapper can map a short fake read" ) {

        // Make the alignment        
        string read = "GAT";
        Alignment aln;
        aln.set_sequence(read);
        
        // Have a list to fill with results
        list<MultipathAlignment> results;
        
        // Align for just one alignment
        mapper.multipath_map(aln, results, 1);
        
        SECTION("there should be one alignment") {
            REQUIRE(results.size() == 1);
            auto& aligned = results.front();
            
            SECTION("which should have one subpath") {
                REQUIRE(aligned.subpath_size() == 1);
                auto& subpath = aligned.subpath(0);
                
                SECTION("which should have one mapping") {
                    REQUIRE(subpath.path().mapping_size() == 1);
                    auto& mapping = subpath.path().mapping(0);
                    
                    SECTION("which should be at the start of node 1") {
                        REQUIRE(mapping.position().node_id() == 1);
                        REQUIRE(mapping.position().offset() == 0);
                        REQUIRE(mapping.position().is_reverse() == false);
                    }
                    
                    SECTION("which should have one edit") {
                        REQUIRE(mapping.edit_size() == 1);
                        auto& edit = mapping.edit(0);
                        
                        SECTION("which should be a length 3 perfect match") {
                            REQUIRE(edit.from_length() == 3);
                            REQUIRE(edit.to_length() == 3);
                            REQUIRE(edit.sequence() == "");
                        }
                    }
                }
            }
        }
    }
    
    SECTION( "MultipathMapper can map two tiny paired reads" ) {
    
        // Here are two reads in opposing, inward-facing directions
        Alignment read1, read2;
        read1.set_sequence("GAT");
        read2.set_sequence("TGT");
        
        // Have a list to fill with results
        list<pair<MultipathAlignment, MultipathAlignment>> results;
        
        // Align for just one pair of alignments
        mapper.multipath_map_paired(read1, read2, 10, results, 1);
        
        SECTION("there should be one pair of alignments") {
            REQUIRE(results.size() == 1);
            auto& aligned1 = results.front().first;
            auto& aligned2 = results.front().second;
            
            SECTION("each read should have one subpath") {
                REQUIRE(aligned1.subpath_size() == 1);
                REQUIRE(aligned2.subpath_size() == 1);
                auto& subpath1 = aligned1.subpath(0);
                auto& subpath2 = aligned2.subpath(0);
                
                SECTION("each should have one mapping") {
                    REQUIRE(subpath1.path().mapping_size() == 1);
                    REQUIRE(subpath2.path().mapping_size() == 1);
                    auto& mapping1 = subpath1.path().mapping(0);
                    auto& mapping2 = subpath2.path().mapping(0);
                    
                    SECTION("the first should be at the start of node 1") {
                        REQUIRE(mapping1.position().node_id() == 1);
                        REQUIRE(mapping1.position().offset() == 0);
                        REQUIRE(mapping1.position().is_reverse() == false);
                    }
                    
                    SECTION("the second should be at the end of node 1, reversed") {
                        REQUIRE(mapping2.position().node_id() == 1);
                        REQUIRE(mapping2.position().offset() == 0);
                        REQUIRE(mapping2.position().is_reverse() == true);
                    }
                    
                    SECTION("each should have one edit") {
                        REQUIRE(mapping1.edit_size() == 1);
                        REQUIRE(mapping2.edit_size() == 1);
                        auto& edit1 = mapping1.edit(0);
                        auto& edit2 = mapping2.edit(0);
                        
                        SECTION("which should be a length 3 perfect match") {
                            REQUIRE(edit1.from_length() == 3);
                            REQUIRE(edit1.to_length() == 3);
                            REQUIRE(edit1.sequence() == "");
                            REQUIRE(edit2.from_length() == 3);
                            REQUIRE(edit2.to_length() == 3);
                            REQUIRE(edit2.sequence() == "");
                        }
                    }
                }
            }
        }
    
    }
    
    // Clean up the GCSA/LCP index
    delete gcsaidx;
    delete lcpidx;
}

}

}
