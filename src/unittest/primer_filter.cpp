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

    TEST_CASE( "filter simple primers",
                "[primer_filter]" ) {
        
        SnarlDistanceIndex distance_index;
        unique_ptr<handlegraph::PathPositionHandleGraph> graph;
        string snarl_index_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.dist";
        string xg_graph_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.xg";
        string primers_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.primer3.out";
        distance_index.deserialize(snarl_index_path);
        graph = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_graph_path);
        PrimerFinder primer_finder(graph, "y", &distance_index);
        primer_finder.load_primers(primers_path);
        
        SECTION("Loads the correct number of primer pairs") {
            REQUIRE(primer_finder.get_primer_pairs().size() == 5);
        }

        SECTION("Loads and processes the primers correctly") {
            primer_finder.add_primer_pair(9, 14, 20, 22, 0, 20);
            const vector<PrimerPair>& primer_pairs = primer_finder.get_primer_pairs();
            const vector<PrimerPair>& selected_primer_pairs = primer_finder.get_selected_primer_pairs();
            const PrimerPair& pair_0 = primer_pairs[0]; // 1st set of primers read from primer3 output. No variation in either primers.
            const PrimerPair& pair_5 = primer_pairs[5]; // made up set of primers. Variation in both priemrs.
                
            SECTION("Check for basic primer attributes") {
                REQUIRE(pair_0.left_primer.sequence == "TGCCTGGCATAGAGGAAAGC");
                REQUIRE(pair_0.left_primer.position == 362);
                REQUIRE(pair_0.left_primer.length == 20);
                REQUIRE(pair_0.right_primer.sequence == "GCCAGAAGAGCCTCAAGGAG");
                REQUIRE(pair_0.right_primer.position == 466);
                REQUIRE(pair_0.right_primer.length == 20);
                REQUIRE(pair_5.left_primer.sequence == "AGCCAGACAAATCTGGGTTC");
                REQUIRE(pair_5.left_primer.position == 181);
                REQUIRE(pair_5.left_primer.length == 20);
                REQUIRE(pair_5.right_primer.sequence == "AGATAATTAAACTGAAGTTC");
                REQUIRE(pair_5.right_primer.position == 260);
                REQUIRE(pair_5.right_primer.length == 20);
            }

            SECTION("Check for minimum and maximum distance") {
                REQUIRE(pair_0.linear_product_size == 124);
                REQUIRE(pair_0.min_product_size == 124);
                REQUIRE(pair_0.max_product_size == 124);
                REQUIRE(pair_5.linear_product_size == 99);
                REQUIRE(pair_5.min_product_size == 97);
                REQUIRE(pair_5.max_product_size == 100);
            }

            SECTION("Check that primers are mapped to correct nodes") {
                vector<size_t> pair_0_left_primer_nodes {27, 8};
                for (size_t i = 0; i < pair_0.left_primer.mapped_nodes_ids.size()-1; i++) {
                    REQUIRE(pair_0.left_primer.mapped_nodes_ids[i] == pair_0_left_primer_nodes[i]);
                }

                vector<size_t> pair_0_right_primer_nodes {33, 34};
                for (size_t i = 0; i < pair_0.right_primer.mapped_nodes_ids.size()-1; i++) {
                    REQUIRE(pair_0.right_primer.mapped_nodes_ids[i] == pair_0_right_primer_nodes[i]);
                }

                vector<size_t> pair_5_left_primer_nodes {9, 11, 12, 14, 15, 17};
                for (size_t i = 0; i < pair_5.left_primer.mapped_nodes_ids.size()-1; i++) {
                    REQUIRE(pair_5.left_primer.mapped_nodes_ids[i] == pair_5_left_primer_nodes[i]);
                }

                vector<size_t> pair_5_right_primer_nodes {22, 24, 25};
                for (size_t i = 0; i < pair_5.right_primer.mapped_nodes_ids.size()-1; i++) {
                    REQUIRE(pair_5.right_primer.mapped_nodes_ids[i] == pair_5_right_primer_nodes[i]);
                }
            }

            SECTION("Check for variation at primer sites") {
                REQUIRE(primer_pairs.size() == 6);
                REQUIRE(selected_primer_pairs.size() == 5);
            }

        }

    }
}
}                   