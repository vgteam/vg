//
//  primers.cpp
//
//  Unit tests for primer filter 
//

#include <stdio.h>
#include <iostream>
#include <regex>
#include <fstream>
#include <vector>
#include <sstream>
#include <set>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "catch.hpp"
#include "random_graph.hpp"
#include "randomness.hpp"
#include "../snarl_distance_index.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../genotypekit.hpp"
#include "../traversal_finder.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>
#include "xg.hpp"
#include "../primer_filter.hpp"

namespace vg {
namespace unittest {

using namespace std;
using namespace Primer_finder;

    TEST_CASE( "filter simple primers",
                "[primers]" ) {
        
        SnarlDistanceIndex distance_index;
        unique_ptr<handlegraph::PathPositionHandleGraph> graph;
        string snarl_index_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.dist";
        string xg_graph_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.xg";
        string primers_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.primer3.out";
        distance_index.deserialize(snarl_index_path);
        graph = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_graph_path);
        Primer_finder primer_finder(graph, "y", &distance_index);
        primer_finder.load_primers(primers_path);
        
        SECTION("Loads the correct number of primer pairs") {
            REQUIRE(primer_finder.get_primer_pairs().size() == 5);
        }

        SECTION("Loads in the sequences correctly") {
            vector<Primer_pair> primer_pairs = primer_finder.get_primer_pairs();
            REQUIRED(primer_pairs[0].left_primer.sequence == 'TGCCTGGCATAGAGGAAAGC');
            REQUIRED(primer_pairs[0].left_primer.position == 362);
            REQUIRED(primer_pairs[0].left_primer.length == 20);
            REQUIRED(primer_pairs[0].right_primer.sequence == 'GCCAGAAGAGCCTCAAGGAG');
            REQUIRED(primer_pairs[0].right_primer.position == 466);
            REQUIRED(primer_pairs[0].right_primer.length == 20);
        }

    }
}
}                   