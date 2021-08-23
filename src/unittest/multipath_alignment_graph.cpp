/// \file multipath_alignment_graph.cpp
///  
/// unit tests for the multipath mapper's MultipathAlignmentGraph

#include <iostream>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../cactus_snarl_finder.hpp"
#include "../multipath_alignment_graph.hpp"
#include "../min_distance.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

TEST_CASE( "MultipathAlignmentGraph::align handles tails correctly", "[multipath][mapping][multipathalignmentgraph]" ) {

    string graph_json = R"({
        "node": [
            {"id": 1, "sequence": "GATT"},
            {"id": 2, "sequence": "A"},
            {"id": 3, "sequence": "G"},
            {"id": 4, "sequence": "C"},
            {"id": 5, "sequence": "A"},
            {"id": 6, "sequence": "G"},
            {"id": 7, "sequence": "A"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 3},
            {"from": 2, "to": 4},
            {"from": 3, "to": 4},
            {"from": 4, "to": 5},
            {"from": 4, "to": 6},
            {"from": 5, "to": 7},
            {"from": 6, "to": 7}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG vg;
    vg.extend(proto_graph);
    
    // Make snarls on it
    CactusSnarlFinder bubble_finder(vg);
    SnarlManager snarl_manager = bubble_finder.find_snarls();
    MinimumDistanceIndex dist_index(&vg, &snarl_manager);
    
    // We need a fake read
    string read("GATTACAA");
    
    // Pack it into an Alignment.
    // Note that we need to use the Alignment's copy for getting iterators for the MEMs.
    Alignment query;
    query.set_sequence(read);
    
    // Make an Aligner to use for the actual aligning and the scores
    Aligner aligner;
        
    // Make an identity projection translation
    auto identity = MultipathAlignmentGraph::create_identity_projection_trans(vg);
    
    // Make up a fake MEM
    // GCSA range_type is just a pair of [start, end], so we can fake them.
    
    // This will actually own the MEMs
    vector<MaximalExactMatch> mems;
    
    // This will hold our MEMs and their start positions in the imaginary graph.
    // Note that this is also a memcluster_t
    pair<vector<pair<const MaximalExactMatch*, pos_t>>, double> mem_hits;
    mem_hits.second = 1.0;
    
    
    /*
    SECTION ("Works with right tail only") {
    
        // Make a MEM hit over all of node 1's sequence
        mems.emplace_back(query.sequence().begin(), query.sequence().begin() + 4, make_pair(5, 5), 1);
        // Drop it on node 1 where it should sit
        mem_hits.first.emplace_back(&mems.back(), make_pos_t(1, false, 0));
        
        // Make the MultipathAlignmentGraph to test
        MultipathAlignmentGraph mpg(vg, mem_hits, identity);
        
        // Make the output multipath_alignment_t
        multipath_alignment_t out;
        
        SECTION("Tries multiple traversals of snarls in tails") {
        
            // Generate 2 fake tail anchors
            mpg.synthesize_tail_anchors(query, vg, &aligner, 1, 2, false, 100, 0.0);
            
            // Cut new anchors on snarls
            mpg.resect_snarls_from_paths(&snarl_manager, identity, 5);
            
            // Make it align, with alignments per gap/tail
            mpg.align(query, vg, &aligner, true, 2, false, 100, 0.0, 5, out);
            
            // Make sure to topologically sort the resulting alignment. TODO: Should
            // the MultipathAlignmentGraph guarantee this for us by construction?
            topologically_order_subpaths(out);
            
            // Make sure it worked at all
            REQUIRE(out.sequence() == read);
            REQUIRE(out.subpath_size() > 0);
            
            // Get the top 10 optimal alignments
            vector<Alignment> opt = optimal_alignments(out, 100);
            
            // Convert to a set of vectors of visited node IDs
            set<vector<id_t>> got;
            for(size_t i = 0; i < opt.size(); i++) {
                auto& aln = opt[i];
                
                vector<id_t> ids;
                for (auto& mapping : aln.path().mapping()) {
                    ids.push_back(mapping.position().node_id());
                }
                
                // Save each list of visited IDs for checking.
                got.insert(ids);
            }
            
            // Make sure all combinations are there
            REQUIRE(got.count({1, 2, 4, 6, 7}));
            REQUIRE(got.count({1, 2, 4, 5, 7}));
            REQUIRE(got.count({1, 3, 4, 6, 7}));
            REQUIRE(got.count({1, 3, 4, 5, 7}));
            
        }
        
        SECTION("Handles tails when anchors for them are not generated") {
        
            // Make it align, with alignments per gap/tail
            mpg.align(query, vg, &aligner, true, 2, false, 100, 0.0, 5, out);
            
            // Make sure to topologically sort the resulting alignment. TODO: Should
            // the MultipathAlignmentGraph guarantee this for us by construction?
            topologically_order_subpaths(out);
            
            // Make sure it worked at all
            REQUIRE(out.sequence() == read);
            REQUIRE(out.subpath_size() > 0);
            
            // Get the top optimal alignments
            vector<Alignment> opt = optimal_alignments(out, 100);
            
            // Convert to a set of vectors of visited node IDs
            set<vector<id_t>> got;
            for(size_t i = 0; i < opt.size(); i++) {
                auto& aln = opt[i];
                
                vector<id_t> ids;
                for (auto& mapping : aln.path().mapping()) {
                    ids.push_back(mapping.position().node_id());
                }
                
                // Save each list of visited IDs for checking.
                got.insert(ids);
            }
            
            // Make sure the best alignment is there.
            REQUIRE(got.count({1, 2, 4, 5, 7}));
        
        }
        
    }
    */
    
    SECTION ("Works with both tails at once") {
    
        // Make a MEM hit over all of node 4's sequence
        mems.emplace_back(query.sequence().begin() + 5, query.sequence().begin() + 6, make_pair(5, 5), 1);
        // Drop it on node 4 where it should sit
        mem_hits.first.emplace_back(&mems.back(), make_pos_t(4, false, 0));
        
        // Make the MultipathAlignmentGraph to test
        vector<size_t> provenance;
        MultipathAlignmentGraph mpg(vg, mem_hits, identity, provenance);
        
        // Make the output multipath_alignment_t
        multipath_alignment_t out;
        
        
        SECTION("Tries multiple traversals of snarls in tails") {
        
            // Generate 2 fake tail anchors
            mpg.synthesize_tail_anchors(query, vg, &aligner, 1, 2, false, 100, 0.0);
            
            // Cut new anchors on snarls
            mpg.resect_snarls_from_paths(&snarl_manager, nullptr,
                                         MultipathAlignmentGraph::create_projector(identity), 5);
            
            // Make it align, with alignments per gap/tail
            mpg.align(query, vg, &aligner, true, 2, false, 100, 0.0, false, 0, 5, out);
            
            // Make sure to topologically sort the resulting alignment. TODO: Should
            // the MultipathAlignmentGraph guarantee this for us by construction?
            topologically_order_subpaths(out);
            // Make sure it worked at all
            REQUIRE(out.sequence() == read);
            REQUIRE(out.subpath_size() > 0);
            
            // Get the top optimal alignments
            vector<Alignment> opt = optimal_alignments(out, 100);
            
            // Convert to a set of vectors of visited node IDs
            set<vector<id_t>> got;
            for(size_t i = 0; i < opt.size(); i++) {
                auto& aln = opt[i];
                
                vector<id_t> ids;
                for (auto& mapping : aln.path().mapping()) {
                    ids.push_back(mapping.position().node_id());
                }
                
                // Save each list of visited IDs for checking.
                got.insert(ids);
            }
            
            // Make sure all combinations are there
            REQUIRE(got.count({1, 2, 4, 6, 7}));
            REQUIRE(got.count({1, 2, 4, 5, 7}));
            REQUIRE(got.count({1, 3, 4, 6, 7}));
            REQUIRE(got.count({1, 3, 4, 5, 7}));
        }
        
        
        SECTION("Handles tails when anchors for them are not generated") {
            // Make it align, with alignments per gap/tail
            mpg.align(query, vg, &aligner, true, 2, false, 100, 0.0, false, 0, 5, out);
            
            // Make sure to topologically sort the resulting alignment. TODO: Should
            // the MultipathAlignmentGraph guarantee this for us by construction?
            topologically_order_subpaths(out);
            
            // Make sure it worked at all
            REQUIRE(out.sequence() == read);
            REQUIRE(out.subpath_size() > 0);
            
            // Get the top optimal alignments
            vector<Alignment> opt = optimal_alignments(out, 100);
            
            // Convert to a set of vectors of visited node IDs
            set<vector<id_t>> got;
            for(size_t i = 0; i < opt.size(); i++) {
                auto& aln = opt[i];
                
                vector<id_t> ids;
                for (auto& mapping : aln.path().mapping()) {
                    ids.push_back(mapping.position().node_id());
                }
                
                // Save each list of visited IDs for checking.
                got.insert(ids);
            }
            
            // Make sure the best alignment is there.
            REQUIRE(got.count({1, 2, 4, 5, 7}));
        }
        
    }
    
    
    SECTION ("Works with left tail only") {
    
        // Make a MEM hit over all of node 7's sequence
        mems.emplace_back(query.sequence().begin() + 7, query.sequence().begin() + 8, make_pair(5, 5), 1);
        // Drop it on node 7 where it should sit
        mem_hits.first.emplace_back(&mems.back(), make_pos_t(7, false, 0));
        
        // Make the MultipathAlignmentGraph to test
        vector<size_t> provenance;
        MultipathAlignmentGraph mpg(vg, mem_hits, identity, provenance);
        
        // Make the output multipath_alignment_t
        multipath_alignment_t out;
        
        SECTION("Tries multiple traversals of snarls in tails") {
            // Generate 2 fake tail anchors
            mpg.synthesize_tail_anchors(query, vg, &aligner, 1, 2, false, 100, 0.0);
            
            // Cut new anchors on snarls
            mpg.resect_snarls_from_paths(nullptr, &dist_index,
                                         MultipathAlignmentGraph::create_projector(identity), 5);
            
            // Make it align, with alignments per gap/tail
            mpg.align(query, vg, &aligner, true, 2, false, 100, 0.0, false, 0, 5, out);
            
            // Make sure to topologically sort the resulting alignment. TODO: Should
            // the MultipathAlignmentGraph guarantee this for us by construction?
            topologically_order_subpaths(out);
            
            // Make sure it worked at all
            REQUIRE(out.sequence() == read);
            REQUIRE(out.subpath_size() > 0);
            
            // Get the top optimal alignments
            vector<Alignment> opt = optimal_alignments(out, 100);
            
            // Convert to a set of vectors of visited node IDs
            set<vector<id_t>> got;
            for(size_t i = 0; i < opt.size(); i++) {
                auto& aln = opt[i];
                
                vector<id_t> ids;
                for (auto& mapping : aln.path().mapping()) {
                    ids.push_back(mapping.position().node_id());
                }
                
                // Save each list of visited IDs for checking.
                got.insert(ids);
            }
            
            // Make sure all combinations are there
            REQUIRE(got.count({1, 2, 4, 6, 7}));
            REQUIRE(got.count({1, 2, 4, 5, 7}));
            REQUIRE(got.count({1, 3, 4, 6, 7}));
            REQUIRE(got.count({1, 3, 4, 5, 7}));
            
        }
        
        SECTION("Handles tails when anchors for them are not generated") {
            // Make it align, with alignments per gap/tail
            mpg.align(query, vg, &aligner, true, 2, false, 100, 0.0, false, 0, 5, out);
            
            // Make sure to topologically sort the resulting alignment. TODO: Should
            // the MultipathAlignmentGraph guarantee this for us by construction?
            topologically_order_subpaths(out);
            
            // Make sure it worked at all
            REQUIRE(out.sequence() == read);
            REQUIRE(out.subpath_size() > 0);
            
            // Get the top optimal alignments
            vector<Alignment> opt = optimal_alignments(out, 100);
            
            // Convert to a set of vectors of visited node IDs
            set<vector<id_t>> got;
            for(size_t i = 0; i < opt.size(); i++) {
                auto& aln = opt[i];
                
                vector<id_t> ids;
                for (auto& mapping : aln.path().mapping()) {
                    ids.push_back(mapping.position().node_id());
                }
                
                // Save each list of visited IDs for checking.
                got.insert(ids);
            }
            
            // Make sure the best alignment is there.
            REQUIRE(got.count({1, 2, 4, 5, 7})); 
        }
        
    }
    
        
}

TEST_CASE("Corresponding path lengths can be computed", "[multipathalignmentgraph]") {
    
    path_t path;
    auto m0 = path.add_mapping();
    auto e0 = m0->add_edit();
    e0->set_from_length(5);
    e0->set_to_length(5);
    auto e1 = m0->add_edit();
    e1->set_from_length(0);
    e1->set_to_length(5);
    auto m1 = path.add_mapping();
    auto e2 = m1->add_edit();
    e2->set_from_length(5);
    e2->set_to_length(0);
    auto e3 = m1->add_edit();
    e3->set_from_length(5);
    e3->set_to_length(5);
    e3->set_sequence("AAAAA");
    
    
    REQUIRE(corresponding_to_length(path, 0, false) == 0);
    REQUIRE(corresponding_from_length(path, 0, false) == 0);
    REQUIRE(corresponding_to_length(path, 1, false) == 1);
    REQUIRE(corresponding_from_length(path, 1, false) == 1);
    REQUIRE(corresponding_to_length(path, 5, false) == 5);
    REQUIRE(corresponding_from_length(path, 5, false) == 5);
    REQUIRE(corresponding_to_length(path, 7, false) == 10);
    REQUIRE(corresponding_from_length(path, 7, false) == 5);
    REQUIRE(corresponding_to_length(path, 11, false) == 11);
    REQUIRE(corresponding_from_length(path, 11, false) == 11);
    REQUIRE(corresponding_to_length(path, 15, false) == 15);
    REQUIRE(corresponding_from_length(path, 15, false) == 15);
    
    REQUIRE(corresponding_to_length(path, 0, true) == 0);
    REQUIRE(corresponding_from_length(path, 0, true) == 0);
    REQUIRE(corresponding_to_length(path, 1, true) == 1);
    REQUIRE(corresponding_from_length(path, 1, true) == 1);
    REQUIRE(corresponding_to_length(path, 5, true) == 5);
    REQUIRE(corresponding_from_length(path, 5, true) == 5);
    REQUIRE(corresponding_from_length(path, 7, true) == 10);
    REQUIRE(corresponding_to_length(path, 7, true) == 5);
    REQUIRE(corresponding_to_length(path, 11, true) == 11);
    REQUIRE(corresponding_from_length(path, 11, true) == 11);
    REQUIRE(corresponding_to_length(path, 15, true) == 15);
    REQUIRE(corresponding_from_length(path, 15, true) == 15);
    
}

}
}
