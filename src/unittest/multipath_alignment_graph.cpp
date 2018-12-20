/// \file multipath_alignment_graph.cpp
///  
/// unit tests for the multipath mapper's MultipathAlignmentGraph

#include <iostream>
#include "json2pb.h"
#include "vg.pb.h"
#include "../multipath_alignment_graph.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

TEST_CASE( "MultipathAlignmentGraph::align tries multiple traversals of snarls in tails", "[multipath][mapping][multipathalignmentgraph]" ) {

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
    
    // We need a fake read
    string read("GATTACAA");
    
    // Pack it into an Alignment.
    // Note that we need to use the Alignment's copy for getting iterators for the MEMs.
    Alignment query;
    query.set_sequence(read);
    
    // Make up a fake MEM
    // GCSA range_type is just a pair of [start, end], so we can fake them.
    
    // This will actually own the MEMs
    vector<MaximalExactMatch> mems;
    
    // This will hold our MEMs and their start positions in the imaginary graph.
    // Note that this is also a memcluster_t
    vector<pair<const MaximalExactMatch*, pos_t>> mem_hits;
    
    // Make a MEM hit over all of node 1's sequence
    mems.emplace_back(query.sequence().begin(), query.sequence().begin() + 4, make_pair(5, 5), 1);
    // Drop it on node 1 where it should sit
    mem_hits.emplace_back(&mems.back(), make_pos_t(1, false, 0));
    
    // Make an Aligner to use for the actual aligning and the scores
    Aligner aligner;
    
    // Make an identity projection translation
    auto identity = MultipathAlignmentGraph::create_identity_projection_trans(vg);
    
    // Make the MultipathAlignmentGraph to test
    MultipathAlignmentGraph mpg(vg, mem_hits, identity);
    
    // Make the output MultipathAlignment
    MultipathAlignment out;
    
    // Generate 2 fake tail anchors
    mpg.synthesize_tail_anchors(query, vg, &aligner, 2, false);
    
    mpg.to_dot(cerr, &query);
    
    // Cut new anchors on snarls
    mpg.resect_snarls_from_paths(&snarl_manager, identity, 5);
    
    mpg.to_dot(cerr, &query);
    
    // Make it align, with alignments per gap/tail
    mpg.align(query, vg, &aligner, true, 2, false, 5, out);
    
    cerr << pb2json(out) << endl;
    
    // Make sure it worked at all
    REQUIRE(out.sequence() == read);
    REQUIRE(out.subpath_size() > 0);
    
    // Get the top 10 optimal alignments
    vector<Alignment> opt = optimal_alignments(out, 10);
    
    // Convert to a set of vectors of visited node IDs
    set<vector<id_t>> got;
    for(auto& aln : opt) {
        cerr << "Linearization: " << pb2json(aln.path()) << endl;
        vector<id_t> ids;
        for (auto& mapping : aln.path().mapping()) {
            ids.push_back(mapping.position().node_id());
        }
    }
    
    // Make sure all combinations are there
    REQUIRE(got.count({1, 2, 4, 6, 7}));
    REQUIRE(got.count({1, 2, 4, 5, 7}));
    REQUIRE(got.count({1, 3, 4, 6, 7}));
    REQUIRE(got.count({1, 3, 4, 5, 7}));
    
}

}

}
