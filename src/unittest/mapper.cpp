/// \file mapper.cpp
///  
/// unit tests for the mapper

#include <iostream>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include <bdsg/hash_graph.hpp>
#include "../mapper.hpp"
#include "xg.hpp"
#include "../build_index.hpp"
#include "catch.hpp"
#include "../algorithms/alignment_path_offsets.hpp"

namespace vg {
namespace unittest {
    
TEST_CASE( "Mapper can map to a one-node graph", "[mapping][mapper]" ) {
    
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
    
    // Make GCSA quiet
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Make pointers to fill in
    gcsa::GCSA* gcsaidx = nullptr;
    gcsa::LCPArray* lcpidx = nullptr;
    
    // Build the GCSA index
    build_gcsa_lcp(graph, gcsaidx, lcpidx, 16, 3);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    // Make a multipath mapper to map against the graph.
    Mapper mapper(&xg_index, gcsaidx, lcpidx);
    
    SECTION( "Mapper can map a short fake read" ) {

        // Make the alignment        
        string read = "GAT";
        Alignment aln;
        aln.set_sequence(read);
        
        // Align for just one alignment
        auto results = mapper.align_multi(aln);
        
        SECTION("there should be one alignment") {
            REQUIRE(results.size() == 1);
            auto& aligned = results.front();
            
                
            SECTION("which should have one mapping") {
                REQUIRE(aligned.path().mapping_size() == 1);
                auto& mapping = aligned.path().mapping(0);
                
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
    
    SECTION( "Mapper can map two tiny paired reads" ) {
    
        // Here are two reads in opposing, inward-facing directions
        Alignment read1, read2;
        read1.set_sequence("GAT");
        read2.set_sequence("TGT");
        
        // Align for just one pair of alignments
        bool resolve_later = false;
        auto results = mapper.align_paired_multi(read1, read2, resolve_later);
        
        SECTION("there should be one pair of alignments") {
            REQUIRE(results.first.size() == 1);
            REQUIRE(results.second.size() == 1);
            auto& aligned1 = results.first.front();
            auto& aligned2 = results.second.front();
            
                
            SECTION("each should have one mapping") {
                REQUIRE(aligned1.path().mapping_size() == 1);
                REQUIRE(aligned2.path().mapping_size() == 1);
                auto& mapping1 = aligned1.path().mapping(0);
                auto& mapping2 = aligned2.path().mapping(0);
                
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
    
    // Clean up the GCSA/LCP index
    delete gcsaidx;
    delete lcpidx;
}

TEST_CASE( "Mapper finds optimal mapping for read starting with node-border MEM", "[mapping][mapper]" ) {
    
    // We have a node 9999 in here to bust some MEM we don't want, to trigger the condition we are trying to test
    string graph_json = R"(
    {"node":[{"id":1430,"sequence":"GAGATCGTGCTACCGCACTCCATGCACTCTAG"},
    {"id":1428,"sequence":"C"},
    {"id":1429,"sequence":"T"},
    {"id":1431,"sequence":"CCTGGGCAACAGAACGAGATGCTGTC"},
    {"id":1427,"sequence":"AGGTTGGAGTGAGC"},
    {"id":1432,"sequence":"ACA"},
    {"id":1433,"sequence":"ACAACAACAACAACAA"},
    {"id":1426,"sequence":"GAGGCAGGAAAATCACTTGAACCGGGAGGCGG"},
    {"id":1434,"sequence":"T"},
    {"id":1435,"sequence":"C"},
    {"id":1424,"sequence":"C"},
    {"id":1425,"sequence":"T"},
    {"id":1436,"sequence":"AACAACAACAA"},
    {"id":1423,"sequence":"TCGGGAGGC"},
    {"id":1437,"sequence":"T"},
    {"id":1438,"sequence":"C"},
    {"id":1421,"sequence":"T"},
    {"id":1422,"sequence":"C"},
    {"id":1439,"sequence":"AACAACAACAACAA"},
    {"id":1420,"sequence":"TA"},
    {"id":1440,"sequence":"A"},
    {"id":1441,"sequence":"C"},
    {"id":1418,"sequence":"A"},
    {"id":1419,"sequence":"C"},
    {"id":1442,"sequence":"AA"},
    {"id":1417,"sequence":"TGCCTGTAATCCCAG"},
    {"id":1443,"sequence":"C"},
    {"id":1444,"sequence":"A"},
    {"id":1416,"sequence":"AAATACAAGTATTAGCCAGGCATTGTGGCAGG"},
    {"id":1445,"sequence":"TTCTCACATCTAAAACAGAGTTCCTGGTTCCA"},
    {"id":9999,"sequence":"TGTTGTTGTTGTTGTGTTNCTCATTTCGTTGCCAAT"}
    ],"edge":[{"to":1431,"from":1430},
    {"to":1430,"from":1428},
    {"to":1430,"from":1429},
    {"to":1432,"from":1431},
    {"to":1433,"from":1431},
    {"to":1428,"from":1427},
    {"to":1429,"from":1427},
    {"to":1433,"from":1432},
    {"to":1434,"from":1433},
    {"to":1435,"from":1433},
    {"to":1427,"from":1426},
    {"to":1436,"from":1434},
    {"to":1436,"from":1435},
    {"to":1426,"from":1424},
    {"to":1426,"from":1425},
    {"to":1437,"from":1436},
    {"to":1438,"from":1436},
    {"to":1424,"from":1423},
    {"to":1425,"from":1423},
    {"to":1439,"from":1437},
    {"to":1439,"from":1438},
    {"to":1423,"from":1421},
    {"to":1423,"from":1422},
    {"to":1440,"from":1439},
    {"to":1441,"from":1439},
    {"to":1421,"from":1420},
    {"to":1422,"from":1420},
    {"to":1442,"from":1440},
    {"to":1442,"from":1441},
    {"to":1420,"from":1418},
    {"to":1420,"from":1419},
    {"to":1443,"from":1442},
    {"to":1444,"from":1442},
    {"to":1418,"from":1417},
    {"to":1419,"from":1417},
    {"to":1445,"from":1443},
    {"to":1445,"from":1444},
    {"to":1417,"from":1416}],"path":[{"name":"17","mapping":[{"position":{"node_id":1416},"rank":1040},
    {"position":{"node_id":1417},"rank":1041},
    {"position":{"node_id":1419},"rank":1042},
    {"position":{"node_id":1420},"rank":1043},
    {"position":{"node_id":1422},"rank":1044},
    {"position":{"node_id":1423},"rank":1045},
    {"position":{"node_id":1425},"rank":1046},
    {"position":{"node_id":1426},"rank":1047},
    {"position":{"node_id":1427},"rank":1048},
    {"position":{"node_id":1429},"rank":1049},
    {"position":{"node_id":1430},"rank":1050},
    {"position":{"node_id":1431},"rank":1051},
    {"position":{"node_id":1433},"rank":1052},
    {"position":{"node_id":1435},"rank":1053},
    {"position":{"node_id":1436},"rank":1054},
    {"position":{"node_id":1438},"rank":1055},
    {"position":{"node_id":1439},"rank":1056},
    {"position":{"node_id":1441},"rank":1057},
    {"position":{"node_id":1442},"rank":1058},
    {"position":{"node_id":1444},"rank":1059},
    {"position":{"node_id":1445},"rank":1060}]}]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make GCSA quiet
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Make pointers to fill in
    gcsa::GCSA* gcsaidx = nullptr;
    gcsa::LCPArray* lcpidx = nullptr;
    
    // Build the GCSA index
    build_gcsa_lcp(graph, gcsaidx, lcpidx, 16, 3);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    // Make a multipath mapper to map against the graph.
    Mapper mapper(&xg_index, gcsaidx, lcpidx);
    mapper.set_alignment_scores(1, 4, 6, 1, 5);
    mapper.mem_reseed_length = 8;
    
    SECTION( "Mapper can map a read starting with a node-border MEM" ) {
    
        Alignment aln;
        aln.set_sequence("TTGTTGTTGTTGTTGTGTTGTTGGGACAGCCATCTCATTTCGTTGCCAATGCTAGAGTGCATGGAG"
                         "TGCGGTAGCACGATCTCAGCTCACTCCAACATCCGCCTCTCGGTTCAAGAGATTTTCCTGGTTCAG"
                         "CCTCCCGAGTAGCTGGAT");
                         
        // Map to the graph
        auto results = mapper.align_multi(aln);
        
        // We want one alignment.
        REQUIRE(results.size() == 1);

        // We want its first mapping to be to 1433 like the optimal alignment.
        REQUIRE(results.front().path().mapping_size() >= 1);
        REQUIRE(results.front().path().mapping(0).position().node_id() == 1436);
    
    }
    
    // Clean up the GCSA/LCP index
    delete gcsaidx;
    delete lcpidx;
    
}

TEST_CASE( "Mapper can annotate positions correctly on both strands", "[mapper][annotation]" ) {
    
    // This node is 73 bp long
    string graph_json = R"(
    {"node":[
        {"id": 1, "sequence": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"}
    ],
    "path":[
        {"name":"x","mapping":[
            {"position":{"node_id":1},"rank":1}
        ]}
    ]}
    )";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make GCSA quiet
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Make pointers to fill in
    gcsa::GCSA* gcsaidx = nullptr;
    gcsa::LCPArray* lcpidx = nullptr;
    
    // Build the GCSA index
    build_gcsa_lcp(graph, gcsaidx, lcpidx, 16, 3);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    // Make a multipath mapper to map against the graph.
    Mapper mapper(&xg_index, gcsaidx, lcpidx);
    
    SECTION( "Mapper can annotate on the forward strand" ) {
        
        // Load up a forward-strand alignment
        string aln_json = R"(
            {"sequence": "A", "path": {"mapping": [
                {"position": {"node_id": 1, "offset": 5}, "edit": [
                    {"from_length": 1, "to_length": 1}
                ]}
            ]}}
        )";
        Alignment aln;
        json2pb(aln, aln_json.c_str(), aln_json.size());
        
        // Annotate it
        algorithms::annotate_with_initial_path_positions(*mapper.xindex, aln);
                        
        // It should have one refpos
        REQUIRE(aln.refpos_size() == 1);
        
        // It should be on the correct path
        REQUIRE(aln.refpos(0).name() == "x");
        
        // It should be at the correct position
        REQUIRE(aln.refpos(0).offset() == 5);
        
        // It should be on the correct strand
        REQUIRE(aln.refpos(0).is_reverse() == false);
    }
    
    SECTION( "Mapper can annotate on the reverse strand" ) {
        
        // Load up a reverse-strand alignment which is the reverse complement of the forward strand one
        string aln_json = R"(
            {"sequence": "T", "path": {"mapping": [
                {"position": {"node_id": 1, "is_reverse": true, "offset": 67}, "edit": [
                    {"from_length": 1, "to_length": 1}
                ]}
            ]}}
        )";
        Alignment aln;
        json2pb(aln, aln_json.c_str(), aln_json.size());
        
        // Annotate it
        algorithms::annotate_with_initial_path_positions(*mapper.xindex, aln);
                        
        // It should have one refpos
        REQUIRE(aln.refpos_size() == 1);
        
        // It should be on the correct path
        REQUIRE(aln.refpos(0).name() == "x");
        
        // It should be at the correct position
        REQUIRE(aln.refpos(0).offset() == 5);
        
        // It should be on the correct strand
        REQUIRE(aln.refpos(0).is_reverse() == true);
    }
    
    
    SECTION( "Mapper can annotate multi-base paths on the forward strand" ) {
        
        // Load up a forward-strand alignment
        string aln_json = R"(
            {"sequence": "AAAAA", "path": {"mapping": [
                {"position": {"node_id": 1, "offset": 5}, "edit": [
                    {"from_length": 5, "to_length": 5}
                ]}
            ]}}
        )";
        Alignment aln;
        json2pb(aln, aln_json.c_str(), aln_json.size());
        
        // Annotate it
        algorithms::annotate_with_initial_path_positions(*mapper.xindex, aln);
                        
        // It should have one refpos
        REQUIRE(aln.refpos_size() == 1);
        
        // It should be on the correct path
        REQUIRE(aln.refpos(0).name() == "x");
        
        // It should be at the correct position
        REQUIRE(aln.refpos(0).offset() == 5);
        
        // It should be on the correct strand
        REQUIRE(aln.refpos(0).is_reverse() == false);
    }
    
    SECTION( "Mapper can annotate multi-base paths on the reverse strand" ) {
        
        // Load up a reverse-strand alignment which is the reverse complement of the forward strand one
        string aln_json = R"(
            {"sequence": "TTTTT", "path": {"mapping": [
                {"position": {"node_id": 1, "is_reverse": true, "offset": 63}, "edit": [
                    {"from_length": 5, "to_length": 5}
                ]}
            ]}}
        )";
        Alignment aln;
        json2pb(aln, aln_json.c_str(), aln_json.size());
        
        // Annotate it
        algorithms::annotate_with_initial_path_positions(*mapper.xindex, aln);
                        
        // It should have one refpos
        REQUIRE(aln.refpos_size() == 1);
        
        // It should be on the correct path
        REQUIRE(aln.refpos(0).name() == "x");
        
        // It should be at the correct position
        REQUIRE(aln.refpos(0).offset() == 5);
        
        // It should be on the correct strand
        REQUIRE(aln.refpos(0).is_reverse() == true);
    }
    
    SECTION( "Mapper can annotate multi-mapping paths on the forward strand" ) {
        
        // Load up a forward-strand alignment
        string aln_json = R"(
            {"sequence": "AAAA", "path": {"mapping": [
                {"position": {"node_id": 1, "offset": 5}, "edit": [
                    {"from_length": 2, "to_length": 2}
                ]},
                {"position": {"node_id": 1, "offset": 7}, "edit": [
                    {"from_length": 1},
                    {"from_length": 2, "to_length": 2}
                ]}
            ]}}
        )";
        Alignment aln;
        json2pb(aln, aln_json.c_str(), aln_json.size());
        
        // Annotate it
        algorithms::annotate_with_initial_path_positions(*mapper.xindex, aln);
                        
        // It should have one refpos
        REQUIRE(aln.refpos_size() == 1);
        
        // It should be on the correct path
        REQUIRE(aln.refpos(0).name() == "x");
        
        // It should be at the correct position
        REQUIRE(aln.refpos(0).offset() == 5);
        
        // It should be on the correct strand
        REQUIRE(aln.refpos(0).is_reverse() == false);
    }
    
    SECTION( "Mapper can annotate multi-mapping paths on the reverse strand" ) {
        
        // Load up a reverse-strand alignment which is not quite the reverse
        // complement of the forward strand one. The deletion is different.
        string aln_json = R"(
            {"sequence": "TTTT", "path": {"mapping": [
                {"position": {"node_id": 1, "is_reverse": true, "offset": 63}, "edit": [
                    {"from_length": 2, "to_length": 2}
                ]},
                {"position": {"node_id": 1, "is_reverse": true, "offset": 65}, "edit": [
                    {"from_length": 1},
                    {"from_length": 2, "to_length": 2}
                ]}
            ]}}
        )";
        Alignment aln;
        json2pb(aln, aln_json.c_str(), aln_json.size());
        
        // Annotate it
        algorithms::annotate_with_initial_path_positions(*mapper.xindex, aln);
                        
        // It should have one refpos
        REQUIRE(aln.refpos_size() == 1);
        
        // It should be on the correct path
        REQUIRE(aln.refpos(0).name() == "x");
        
        // It should be at the correct position
        REQUIRE(aln.refpos(0).offset() == 5);
        
        // It should be on the correct strand
        REQUIRE(aln.refpos(0).is_reverse() == true);
    }
    
    // Clean up the GCSA/LCP index
    delete gcsaidx;
    delete lcpidx;
    
}

TEST_CASE( "Mapper can walk paths from fan-out MEM algorithm", "[mapping][mapper][mem]" ) {
    
    
    bdsg::HashGraph graph;
    auto h1 = graph.create_handle("GATTGGACACCCATAGC");
    auto h2 = graph.create_handle("TGGCCAC");
    auto h3 = graph.create_handle("AGT");
    auto h4 = graph.create_handle("AGT");
    auto h5 = graph.create_handle("AGT");
    auto h6 = graph.create_handle("AGT");
    auto h7 = graph.create_handle("AGT");
    auto h8 = graph.create_handle("AGTCA");
    auto h9 = graph.create_handle("C");
    auto h10 = graph.create_handle("C");
    auto h11 = graph.create_handle("C");
    
    graph.create_edge(h1, h2);
    graph.create_edge(h2, h3);
    graph.create_edge(h2, h4);
    graph.create_edge(h2, h5);
    graph.create_edge(h2, h6);
    graph.create_edge(h2, h7);
    graph.create_edge(h2, h8);
    graph.create_edge(h3, h9);
    graph.create_edge(h4, h10);
    graph.create_edge(h5, h11);
    
    // Make GCSA quiet
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Make pointers to fill in
    gcsa::GCSA* gcsaidx = nullptr;
    gcsa::LCPArray* lcpidx = nullptr;
    
    // Build the GCSA index
    build_gcsa_lcp(graph, gcsaidx, lcpidx, 16, 3);
    
    // Build the xg index
    xg::XG xg_index;
    xg_index.from_path_handle_graph(graph);
    
    // Make a multipath mapper to map against the graph.
    Mapper mapper(&xg_index, gcsaidx, lcpidx);
    
    SECTION("Path within a node") {
        
        //         GATTGGACACCCATAGC
        //            ...X....
        string seq = "TGGGCACC";
        
        deque<pair<string::const_iterator, char>> fanout_breaks;
        fanout_breaks.emplace_back(seq.begin() + 3, 'A');
        
        gcsa::node_type pos = gcsa::Node::encode(graph.get_id(h1), 3, false);
        
        auto path = mapper.walk_fanout_path(seq.begin(), seq.end(),
                                            fanout_breaks, pos);
        
        REQUIRE(path.size() == 1);
        REQUIRE(id(path[0]) == graph.get_id(h1));
        REQUIRE(offset(path[0]) == 3);
        REQUIRE(is_rev(path[0]) == false);
    }
    
    SECTION("Path across two nodes") {
        
        //                 |
        //GATTGGACACCCATAGCTGGCCAC
        //              ..X.X....
        string seq =   "AGTTTGCCA";
        
        deque<pair<string::const_iterator, char>> fanout_breaks;
        fanout_breaks.emplace_back(seq.begin() + 2, 'C');
        fanout_breaks.emplace_back(seq.begin() + 4, 'G');
        
        gcsa::node_type pos = gcsa::Node::encode(graph.get_id(h1), 14, false);
        
        auto path = mapper.walk_fanout_path(seq.begin(), seq.end(),
                                            fanout_breaks, pos);
        
        REQUIRE(path.size() == 2);
        REQUIRE(id(path[0]) == graph.get_id(h1));
        REQUIRE(offset(path[0]) == 14);
        REQUIRE(is_rev(path[0]) == false);
        REQUIRE(id(path[1]) == graph.get_id(h2));
        REQUIRE(offset(path[1]) == 0);
        REQUIRE(is_rev(path[1]) == false);
    }
    
    SECTION("Path across two nodes that requires backtracking") {
        
        //                |
        //         TGGCCACAGTCA
        //            ..X.X..X.
        string seq = "CCGCGGTTA";
        
        deque<pair<string::const_iterator, char>> fanout_breaks;
        fanout_breaks.emplace_back(seq.begin() + 2, 'A');
        fanout_breaks.emplace_back(seq.begin() + 4, 'A');
        fanout_breaks.emplace_back(seq.begin() + 7, 'C');
        
        gcsa::node_type pos = gcsa::Node::encode(graph.get_id(h2), 3, false);
        
        auto path = mapper.walk_fanout_path(seq.begin(), seq.end(),
                                            fanout_breaks, pos);
        
        REQUIRE(path.size() == 2);
        REQUIRE(id(path[0]) == graph.get_id(h2));
        REQUIRE(offset(path[0]) == 3);
        REQUIRE(is_rev(path[0]) == false);
        REQUIRE(id(path[1]) == graph.get_id(h8));
        REQUIRE(offset(path[1]) == 0);
        REQUIRE(is_rev(path[1]) == false);
    }
    
    delete gcsaidx;
    delete lcpidx;
}

}
}
