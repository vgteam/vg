/// \file multipath_mapper.cpp
///  
/// unit tests for the multipath mapper

#include <iostream>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../multipath_mapper.hpp"
#include "../build_index.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

// We define a child class to expose all the protected stuff for testing
class TestMultipathMapper : public MultipathMapper {
public:
    using MultipathMapper::MultipathMapper;
    using MultipathMapper::multipath_map_internal;
    using MultipathMapper::attempt_unpaired_multipath_map_of_pair;
    using MultipathMapper::align_to_cluster_graphs;
    using MultipathMapper::align_to_cluster_graph_pairs;
    using MultipathMapper::query_cluster_graphs;
    using MultipathMapper::multipath_align;
    using MultipathMapper::strip_full_length_bonuses;
    using MultipathMapper::sort_and_compute_mapping_quality;
    using MultipathMapper::set_read_coverage;
};

TEST_CASE( "MultipathMapper::read_coverage works", "[multipath][mapping][multipathmapper]" ) {

    // Make up some MEMs
    // GCSA range_type is just a pair of [start, end], so we can fake them.
    
    // This will actually own the MEMs
    vector<MaximalExactMatch> mems;
    
    // This will hold our MEMs and their start positions in the imaginary graph.
    MultipathMapper::clustergraph_t cluster_graph;
    pair<vector<pair<const MaximalExactMatch*, pos_t>>, double>& mem_hits = get<1>(cluster_graph);
    mem_hits.second = 1.0;
    
    // We need a fake read
    string read("GATTACA");
    
    SECTION("No hits cover no bases") {
        TestMultipathMapper::set_read_coverage(cluster_graph);
        REQUIRE(get<2>(cluster_graph) == 0);
    }
    
    SECTION("A hit of zero length covers 0 bases") {
        // Make a MEM hit
        mems.emplace_back(read.begin(), read.begin(), make_pair(5, 5), 1);
        // Drop it on some arbitrary node
        mem_hits.first.emplace_back(&mems.back(), make_pos_t(999, false, 3));
        
        TestMultipathMapper::set_read_coverage(cluster_graph);
        REQUIRE(get<2>(cluster_graph) == 0);
    }
    
    SECTION("A hit that one-past-the-ends just after its start position covers 1 base") {
        // Make a MEM hit
        mems.emplace_back(read.begin(), read.begin() + 1, make_pair(5, 5), 1);
        // Drop it on some arbitrary node
        mem_hits.first.emplace_back(&mems.back(), make_pos_t(999, false, 3));
        
        TestMultipathMapper::set_read_coverage(cluster_graph);
        REQUIRE(get<2>(cluster_graph) == 1);
    }
    
    SECTION("A hit that ends at the end of the string covers the string") {
        // Make a MEM hit
        mems.emplace_back(read.begin(), read.end(), make_pair(5, 5), 1);
        // Drop it on some arbitrary node
        mem_hits.first.emplace_back(&mems.back(), make_pos_t(999, false, 3));
        
        TestMultipathMapper::set_read_coverage(cluster_graph);
        REQUIRE(get<2>(cluster_graph) == read.size());
    }
    
    SECTION("Two overlapping hits still only cover the string once") {
        // Make a MEM hit
        mems.emplace_back(read.begin(), read.begin() + 5, make_pair(5, 5), 1);
        // Drop it on some arbitrary node
        // Do it again
        mems.emplace_back(read.begin() + 2, read.end(), make_pair(6, 6), 1);
        
        // Point into vector *after* it's done.
        mem_hits.first.emplace_back(&mems[0], make_pos_t(999, false, 3));
        mem_hits.first.emplace_back(&mems[1], make_pos_t(888, false, 3));
        
        TestMultipathMapper::set_read_coverage(cluster_graph);
        REQUIRE(get<2>(cluster_graph) == read.size());
    }
    
    SECTION("Two abutting hits cover the string") {
        // Make a MEM hit
        mems.emplace_back(read.begin(), read.begin() + 3, make_pair(5, 5), 1);
        // Do it again
        mems.emplace_back(read.begin() + 3, read.end(), make_pair(6, 6), 1);

        // Point into vector *after* it's done.
        mem_hits.first.emplace_back(&mems[0], make_pos_t(999, false, 3));
        mem_hits.first.emplace_back(&mems[1], make_pos_t(888, false, 3));
        
        TestMultipathMapper::set_read_coverage(cluster_graph);
        REQUIRE(get<2>(cluster_graph) == read.size());
    }

}

TEST_CASE( "MultipathMapper::query_cluster_graphs works", "[multipath][mapping][multipathmapper]" ) {
    
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
    xg_index.from_path_handle_graph(VG(proto_graph));
    
    // Make a multipath mapper to map against the graph.
    TestMultipathMapper mapper(&xg_index, gcsaidx, lcpidx);
    
    // Make an Alignment that we're pretending we're doing
    Alignment aln;
    aln.set_sequence("GATTACA");
    const string& read = aln.sequence();
    
    // Make some MEMs and their clusters
    vector<MaximalExactMatch> mems;
    vector<MultipathMapper::memcluster_t> clusters;

    // Remember the important types:
    // vector<clustergraph_t>
    // using clustergraph_t = tuple<VG*, memcluster_t, size_t>;
    // using memcluster_t = vector<pair<const MaximalExactMatch*, pos_t>>;
    
    SECTION("no MEMs produce no graphs") {
    
        auto results = mapper.query_cluster_graphs(aln, mems, clusters);
        
        REQUIRE(results.size() == 0);
    
    }
    
    SECTION("a single MEM produces one graph with the whole node") {
        
        // Make MEMs
        mems.emplace_back(read.begin(), read.begin() + 3, make_pair(5, 5), 1);
        // Fill in MEM hits with node_types
        mems.back().nodes.push_back(gcsa::Node::encode(1, 0));
        
        // Make a cluster
        clusters.resize(1);
        clusters.back().first.emplace_back(&mems.at(0), make_pos_t(1, false, 0));
        
        REQUIRE(mems.size() == 1);
        
        // Now get the graph
        auto results = mapper.query_cluster_graphs(aln, mems, clusters);
        
        // We have one graph
        REQUIRE(results.size() == 1);
        // It has one node
        REQUIRE(get<0>(results[0])->get_node_count() == 1);
        // It contains the one MEM we fed in
        REQUIRE(get<1>(results[0]).first.size() == 1);
        MultipathMapper::memcluster_t& assigned_mems = get<1>(results[0]);
        const MaximalExactMatch* mem = assigned_mems.first[0].first;
        pos_t where = assigned_mems.first[0].second;
        REQUIRE(mem == &mems.back());
        REQUIRE(where == make_pos_t(1, false, 0));
        // It covers 3 bases from that MEM
        REQUIRE(get<2>(results[0]) == 3);
    
    }
    
    SECTION("two MEMs close together in one cluster make one graph") {
    
        // We will use two fake overlapping MEMs
        mems.emplace_back(read.begin(), read.begin() + 3, make_pair(5, 5), 1);
        mems.back().nodes.push_back(gcsa::Node::encode(1, 0));
        mems.emplace_back(read.begin() + 3, read.end(), make_pair(6, 6), 1);
        mems.back().nodes.push_back(gcsa::Node::encode(1, 3));
        clusters.resize(1);
        clusters.back().first.emplace_back(&mems[0], make_pos_t(1, false, 0));
        clusters.back().first.emplace_back(&mems[1], make_pos_t(1, false, 3));
        
        REQUIRE(mems.size() == 2);

        // Now get the graph
        auto results = mapper.query_cluster_graphs(aln, mems, clusters);
        
        // We have one graph
        REQUIRE(results.size() == 1);
        // It has one node
        REQUIRE(get<0>(results[0])->get_node_count() == 1);
        // It came from two MEM hits
        REQUIRE(get<1>(results[0]).first.size() == 2);
        // They are hits of the two MEMs we fed in at the right places
        set<pair<const MaximalExactMatch*, pos_t>> found{get<1>(results[0]).first.begin(), get<1>(results[0]).first.end()};
        set<pair<const MaximalExactMatch*, pos_t>> wanted{clusters.back().first.begin(), clusters.back().first.end()};
        REQUIRE(found == wanted);
        // It covers all 7 bases, like the MEMs do together
        REQUIRE(get<2>(results[0]) == 7);
    }
    
    
    SECTION("two MEMs close together in two clusters make one graph") {
    
        // We will use two fake overlapping MEMs
        mems.emplace_back(read.begin(), read.begin() + 3, make_pair(5, 5), 1);
        mems.back().nodes.push_back(gcsa::Node::encode(1, 0));
        mems.emplace_back(read.begin() + 3, read.end(), make_pair(6, 6), 1);
        mems.back().nodes.push_back(gcsa::Node::encode(1, 3));
        clusters.resize(2);
        clusters.front().first.emplace_back(&mems[0], make_pos_t(1, false, 0));
        clusters.back().first.emplace_back(&mems[1], make_pos_t(1, false, 3));
        
        REQUIRE(mems.size() == 2);
        
        // Now get the graph
        auto results = mapper.query_cluster_graphs(aln, mems, clusters);
        
        // We have one graph
        REQUIRE(results.size() == 1);
        // It has one node
        REQUIRE(get<0>(results[0])->get_node_count() == 1);
        // It came from two MEM hits
        REQUIRE(get<1>(results[0]).first.size() == 2);
        // They are hits of the two MEMs we fed in at the right places
        set<pair<const MaximalExactMatch*, pos_t>> found{get<1>(results[0]).first.begin(), get<1>(results[0]).first.end()};
        set<pair<const MaximalExactMatch*, pos_t>> wanted{clusters.front().first[0], clusters.back().first[0]};
        REQUIRE(found == wanted);
        // It covers all 7 bases, like the MEMs do together
        REQUIRE(get<2>(results[0]) == 7);
    
    }
}
    
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
    MultipathMapper mapper(&xg_index, gcsaidx, lcpidx);
    // Lower the max mapping quality so that it thinks it can find unambiguous mappings of
    // short sequences
    mapper.max_mapping_quality = 10;
    
    SECTION( "MultipathMapper can map a short fake read" ) {

        // Make the alignment        
        string read = "GAT";
        Alignment aln;
        aln.set_sequence(read);
        
        // Have a list to fill with results
        vector<multipath_alignment_t> results;
        
        // Align for just one alignment
        mapper.multipath_map(aln, results);
        
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
    
    // Give it a fragment length distribution so it doesn't try to learn one
    mapper.force_fragment_length_distr(4, 2);
    
    SECTION( "MultipathMapper can map two tiny paired reads" ) {
    
        // Here are two reads on the same strand, as the mapper requires.
        Alignment read1, read2;
        read1.set_sequence("GAT");
        read2.set_sequence("ACA");
        
        // Have a list to fill with results
        vector<pair<multipath_alignment_t, multipath_alignment_t>> results;
        vector<pair<Alignment, Alignment>> buffer;
        
        // Align for just one pair of alignments
        mapper.multipath_map_paired(read1, read2, results, buffer);
        
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
                    
                    SECTION("the second should be further along the same node") {
                        REQUIRE(mapping2.position().node_id() == 1);
                        REQUIRE(mapping2.position().offset() == 4);
                        REQUIRE(mapping2.position().is_reverse() == false);
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

TEST_CASE( "MultipathMapper can work on a bigger graph", "[multipath][mapping][multipathmapper]" ) {
    
    string graph_json = R"({"node":[{"sequence":"CTTCTCATCCCTCCTCAAGGGCCTTTAACTACTCCACATCCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAGATG","id":12},{"sequence":"A","id":2},{"sequence":"CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTT","id":3},{"sequence":"TTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAATGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTT","id":4},{"sequence":"CAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGAACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTA","id":5},{"sequence":"CTGACAGCAGAGGTTTGCTGTGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTT","id":6},{"sequence":"ACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCTTCTGGCTTTTCATTGTCAAC","id":7},{"sequence":"ACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCTCACCGTGATCTTCAAGTTTGAAAATTGCATCTCAAATCTA","id":8},{"sequence":"AGACCCAGAGGGCTCACCCAGAGTCGAGGCTCAAGGACAGCTCTCCTTTGTGTCCAGAGTGTATACGATGTAACTCTGTTCGGGCACTGGTGAAAGATAA","id":9},{"sequence":"CAGAGGAAATGCCTGGCTTTTTATCAGAACATGTTTCCAAGCTTATCCCTTTTCCCAGCTCTCCTTGTCCCTCCCAAGATCTCTTCACTGGCCTCTTATC","id":10},{"sequence":"TTTACTGTTACCAAATCTTTCCAGAAGCTGCTCTTTCCCTCAATTGTTCATTTGTCTTCTTGTCCAGGAATGAACCACTGCTCTCTTCTTGTCAGATCAG","id":11}],"path":[{"name":"x","mapping":[{"position":{"node_id":3},"edit":[{"from_length":100,"to_length":100}],"rank":1},{"position":{"node_id":4},"edit":[{"from_length":100,"to_length":100}],"rank":2},{"position":{"node_id":5},"edit":[{"from_length":100,"to_length":100}],"rank":3},{"position":{"node_id":6},"edit":[{"from_length":100,"to_length":100}],"rank":4},{"position":{"node_id":7},"edit":[{"from_length":100,"to_length":100}],"rank":5},{"position":{"node_id":8},"edit":[{"from_length":100,"to_length":100}],"rank":6},{"position":{"node_id":9},"edit":[{"from_length":100,"to_length":100}],"rank":7},{"position":{"node_id":10},"edit":[{"from_length":100,"to_length":100}],"rank":8},{"position":{"node_id":11},"edit":[{"from_length":100,"to_length":100}],"rank":9},{"position":{"node_id":12},"edit":[{"from_length":100,"to_length":100}],"rank":10},{"position":{"node_id":2},"edit":[{"from_length":1,"to_length":1}],"rank":11}]}],"edge":[{"from":12,"to":2},{"from":3,"to":4},{"from":4,"to":5},{"from":5,"to":6},{"from":6,"to":7},{"from":7,"to":8},{"from":8,"to":9},{"from":9,"to":10},{"from":10,"to":11},{"from":11,"to":12}]})";
    
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
    xg_index.from_path_handle_graph(VG(proto_graph));
    
    // Make a multipath mapper to map against the graph.
    TestMultipathMapper mapper(&xg_index, gcsaidx, lcpidx);
    // Lower the max mapping quality so that it thinks it can find unambiguous mappings of
    // short sequences
    mapper.max_mapping_quality = 10;
    
    
    SECTION( "topologically_order_subpaths works within a node" ) {    
        string aln_json = R"(
        {
            "subpath": [
                {"path": {"mapping": [
                    {"position": {"node_id": 3, "offset": 11}, "edit": [{"from_length": 10, "to_length": 10}]}
                ]}},
                {"path": {"mapping": [
                    {"position": {"node_id": 3, "offset": 0}, "edit": [{"from_length": 10, "to_length": 10}]}
                ]}, "next": [0]}
            ]
        }
        )";
    
        MultipathAlignment disordered_pb;
        json2pb(disordered_pb, aln_json.c_str(), aln_json.size());
        multipath_alignment_t disordered;
        from_proto_multipath_alignment(disordered_pb, disordered);
        
        topologically_order_subpaths(disordered);
        
        REQUIRE(disordered.subpath_size() == 2);
        // First subpath (used to be last) first
        REQUIRE(disordered.subpath(0).path().mapping(0).position().node_id() == 3);
        REQUIRE(disordered.subpath(0).path().mapping(0).position().offset() == 0);
        // Then second
        REQUIRE(disordered.subpath(1).path().mapping(0).position().node_id() == 3);
        REQUIRE(disordered.subpath(1).path().mapping(0).position().offset() == 11);
        
    }
    
    SECTION( "topologically_order_subpaths works between nodes" ) {
        string aln_json = R"(
        {
            "subpath": [
                {"path": {"mapping": [
                    {"position": {"node_id": 4}, "edit": [{"from_length": 10, "to_length": 10}]}
                ]}, "next": [1]},
                {"path": {"mapping": [
                    {"position": {"node_id": 4, "offset": 11}, "edit": [{"from_length": 10, "to_length": 10}]}
                ]}},
                {"path": {"mapping": [
                    {"position": {"node_id": 3, "offset": 5}, "edit": [{"from_length": 10, "to_length": 10}]}
                ]}, "next": [0]}
            ]
        }
        )";
    
        MultipathAlignment disordered_pb;
        json2pb(disordered_pb, aln_json.c_str(), aln_json.size());
        multipath_alignment_t disordered;
        from_proto_multipath_alignment(disordered_pb, disordered);
        
        topologically_order_subpaths(disordered);
        
        REQUIRE(disordered.subpath_size() == 3);
        // First subpath (used to be last) first
        REQUIRE(disordered.subpath(0).path().mapping(0).position().node_id() == 3);
        REQUIRE(disordered.subpath(0).path().mapping(0).position().offset() == 5);
        // Then second
        REQUIRE(disordered.subpath(1).path().mapping(0).position().node_id() == 4);
        REQUIRE(disordered.subpath(1).path().mapping(0).position().offset() == 0);
        // Then third
        REQUIRE(disordered.subpath(2).path().mapping(0).position().node_id() == 4);
        REQUIRE(disordered.subpath(2).path().mapping(0).position().offset() == 11);
    
    }
    
    SECTION( "MultipathMapper buffers pairs that don't map the first time" ) {
        // Here are two reads on the same strand
        Alignment read1, read2;
        read1.set_sequence("CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAAC");
        read2.set_sequence("TCCTT");
        
        // Have a list to fill with results
        vector<pair<multipath_alignment_t, multipath_alignment_t>> results;
        vector<pair<Alignment, Alignment>> buffer;
        
        // Align for just one pair of alignments
        mapper.multipath_map_paired(read1, read2, results, buffer);
        
        // The second read was ambiguous so we should have buffered this read.
        REQUIRE(results.empty());
        REQUIRE(buffer.size() == 1);
        
        // Did we buffer the right thing?
        REQUIRE(buffer[0].first.sequence() == read1.sequence());
        REQUIRE(buffer[0].second.sequence() == read2.sequence());
    }
    
    SECTION( "MultipathMapper does not buffer unambiguous pairs" ) {
        // Here are two reads on the same strand
        Alignment read1, read2;
        read1.set_sequence("CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAAC");
        read2.set_sequence("TCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCT");
        
        // Have a list to fill with results
        vector<pair<multipath_alignment_t, multipath_alignment_t>> results;
        vector<pair<Alignment, Alignment>> buffer;
        
        // Align for just one pair of alignments
        mapper.multipath_map_paired(read1, read2, results, buffer);
        
        // The pair was not ambiguous so we should have produced a result.
        REQUIRE(results.size() == 1);
        REQUIRE(buffer.empty());
        
        // Did we map the right thing?
        REQUIRE(results[0].first.sequence() == read1.sequence());
        REQUIRE(results[0].second.sequence() == read2.sequence());
        
        // Did we get alignments?
        REQUIRE(results[0].first.subpath_size() > 0);
        REQUIRE(results[0].second.subpath_size() > 0);
    }
    
    // Fix the distribution
    mapper.force_fragment_length_distr(100, 200);
    
    SECTION( "MultipathMapper aligns ambiguous pairs after the fragment length distribution is set" ) {
        // Here are two reads on the same strand
        Alignment read1, read2;
        read1.set_sequence("CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAAC");
        read2.set_sequence("TCCTT");
        
        // Have a list to fill with results
        vector<pair<multipath_alignment_t, multipath_alignment_t>> results;
        vector<pair<Alignment, Alignment>> buffer;
        
        // Align for just one pair of alignments
        mapper.multipath_map_paired(read1, read2, results, buffer);
        
        // The distribution has been estimated so we should have produced a result.
        REQUIRE(results.size() == 1);
        REQUIRE(buffer.empty());
        
        // Did we map the right thing?
        REQUIRE(results[0].first.sequence() == read1.sequence());
        REQUIRE(results[0].second.sequence() == read2.sequence());
        
        // Did we get alignments?
        REQUIRE(results[0].first.subpath_size() > 0);
        REQUIRE(results[0].second.subpath_size() > 0);
        
        // TODO: this error check is no longer accurate now that I've implemented the fragment length
        // distribution into the mapping quality
        // But it should have MAPQ 0 for the second, ambiguously-placed read.
        //REQUIRE(results[0].second.mapping_quality() == 0);
        // TODO: It also zeros MAPQ for the first read; is that smart?
    }
    
    // Clean up the GCSA/LCP index
    delete gcsaidx;
    delete lcpidx;
}

}

}
