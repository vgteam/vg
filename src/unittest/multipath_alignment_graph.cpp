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
            {"id": 4, "sequence": "CA"}
        ],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 3},
            {"from": 2, "to": 4},
            {"from": 3, "to": 4}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG vg;
    vg.extend(proto_graph);
    
    // We need a fake read
    string read("GATTACA");
    
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
    
    // Make a MEM hit
    mems.emplace_back(query.sequence().begin(), ++query.sequence().begin(), make_pair(5, 5), 1);
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
    
    // Make it align
    mpg.align(query, vg, &aligner, true, 4, false, 5, out);
    
    cerr << pb2json(out) << endl;
    
    // Make sure it worked at all
    REQUIRE(out.sequence() == read);
    REQUIRE(out.subpath_size() > 0);
}

}

}
