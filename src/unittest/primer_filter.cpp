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
#include "../recombinator.hpp"

namespace vg {
namespace unittest {

using namespace std;


    TEST_CASE( "filter simple primers",
                "[primer_filter]" ) {
        
        SnarlDistanceIndex distance_index;
        unique_ptr<handlegraph::PathPositionHandleGraph> graph;
        gbwtgraph::GBWTGraph gbwt_graph;
        gbwt::GBWT gbwt_index;
        gbwt::FastLocate r_index;
        string snarl_index_path = "test/primers/y.dist";
        string xg_graph_path = "test/primers/y.xg";
        distance_index.deserialize(snarl_index_path);
        graph = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_graph_path);
        load_r_index(r_index, "test/primers/y.ri");
        load_gbz(gbwt_index, gbwt_graph, "test/primers/y.giraffe.gbz");
        gbwt_graph.set_gbwt(gbwt_index);
        r_index.setGBWT(gbwt_index);
        
        SECTION("template_position=0") {
            string primers_path = "test/primers/y.primer3_with_ref_pos.out";
            ifstream file_handle(primers_path);
            PrimerFinder primer_finder(graph, &distance_index, file_handle, gbwt_graph, gbwt_index, r_index);

            SECTION("Loads the correct number of chromosomes") {
                REQUIRE(primer_finder.total_reference_paths() == 1);
            }

            SECTION("Loads the correct number of primer pairs") {
                REQUIRE(primer_finder.get_primer_pairs_of_chrom("y").size() == 5);
            }

            SECTION("Loads and processes the primers correctly") {
                primer_finder.add_primer_pair("y", 9, 14, 20, 22, 0, 20); // made up data, variation both at primers and in product
                primer_finder.add_primer_pair("y", 31, 0, 15, 34, 1, 15); // made up data, no variation at primers or in product

                // Correct primer attributes
                const vector<string> left_primers_sequences {
                    "TGCCTGGCATAGAGGAAAGC", "GAGTCGAGGCTCAAGGACAG", "CAGAGTCGAGGCTCAAGGAC",
                    "GAGGCTCAAGGACAGCTCTC", "TCCAGAAGCTGCTCTTTCCC", "AGCCAGACAAATCTGGGTTC",
                    "CAACTGGTAGTTACT"
                };

                const vector<size_t> left_primers_positions {
                    362, 620, 618, 625, 819, 181, 388
                };

                const vector<size_t> left_primers_lengths {
                    20, 20, 20, 20, 20, 20, 15
                };

                const vector<size_t> left_primers_nodes_count {
                    2, 1, 1, 2, 2, 6, 1
                };

                const vector<string> right_primers_sequences {
                    "GCCAGAAGAGCCTCAAGGAG", "AGGAGAGCTGGGAAAAGGGA", "AGGAGAGCTGGGAAAAGGGA",
                    "AGGAGAGCTGGGAAAAGGGA", "GCCTGGGTAGCTTTGGATGT", "AGATAATTAAACTGAAGTTC",
                    "GTTGACAATGAAAAG"
                };

                const vector<size_t> right_primers_positions {
                    466, 745, 745, 745, 935, 260, 485
                };

                const vector<size_t> right_primers_lengths {
                    20, 20, 20, 20, 20, 20, 15
                };

                const vector<size_t> right_primers_nodes_count {
                    2, 1, 1, 1, 2, 3, 1
                };

                const vector<size_t> min_product_sizes {
                    124, 142, 144, 137, 136, 99, 112
                };

                const vector<size_t> max_product_sizes {
                    124, 145, 147, 140, 137, 99, 112
                };

                const vector<size_t> linear_product_sizes {
                    124, 145, 147, 140, 136, 99, 112
                };
                
                const vector<double> variation_level {
                    1.0, 1.0, 1.0, 1.0, 1.0, 0.33333, 1.0
                };


                const vector<PrimerPair>& primer_pairs = primer_finder.get_primer_pairs_of_chrom("y");
                
                REQUIRE(primer_pairs.size() == left_primers_sequences.size());
                for (size_t i = 0; i < primer_pairs.size(); ++i) {
                    REQUIRE(left_primers_nodes_count[i]  == primer_pairs[i].left_primer.mapped_nodes_ids.size());
                    REQUIRE(left_primers_sequences[i]    == primer_pairs[i].left_primer.sequence);
                    REQUIRE(left_primers_positions[i]    == primer_pairs[i].left_primer.position_chromosome);
                    REQUIRE(left_primers_lengths[i]      == primer_pairs[i].left_primer.length);
                    REQUIRE(right_primers_nodes_count[i] == primer_pairs[i].right_primer.mapped_nodes_ids.size());
                    REQUIRE(right_primers_sequences[i]   == primer_pairs[i].right_primer.sequence);
                    REQUIRE(right_primers_positions[i]   == primer_pairs[i].right_primer.position_chromosome);
                    REQUIRE(right_primers_lengths[i]     == primer_pairs[i].right_primer.length);
                    REQUIRE(linear_product_sizes[i]      == primer_pairs[i].linear_product_size);
                    REQUIRE(min_product_sizes[i]         == primer_pairs[i].min_product_size);
                    REQUIRE(max_product_sizes[i]         == primer_pairs[i].max_product_size);
                    REQUIRE(abs(variation_level[i] - primer_pairs[i].variation_level) <= 0.0001);
                }

                SECTION("Check that primers are assigned with correct nodes") {
                    vector<size_t> pair_0_left_primer_nodes {27, 28};
                    for (size_t i = 0; i < primer_pairs[0].left_primer.mapped_nodes_ids.size(); i++) {
                        REQUIRE(primer_pairs[0].left_primer.mapped_nodes_ids[i] == pair_0_left_primer_nodes[i]);
                    }

                    vector<size_t> pair_0_right_primer_nodes {33, 34};
                    for (size_t i = 0; i < primer_pairs[0].right_primer.mapped_nodes_ids.size(); i++) {
                        REQUIRE(primer_pairs[0].right_primer.mapped_nodes_ids[i] == pair_0_right_primer_nodes[i]);
                    }

                    vector<size_t> pair_5_left_primer_nodes {9, 11, 12, 14, 15, 17};
                    for (size_t i = 0; i < primer_pairs[5].left_primer.mapped_nodes_ids.size(); i++) {
                        REQUIRE(primer_pairs[5].left_primer.mapped_nodes_ids[i] == pair_5_left_primer_nodes[i]);
                    }

                    vector<size_t> pair_5_right_primer_nodes {22, 24, 25};
                    for (size_t i = 0; i < primer_pairs[5].right_primer.mapped_nodes_ids.size(); i++) {
                        REQUIRE(primer_pairs[5].right_primer.mapped_nodes_ids[i] == pair_5_right_primer_nodes[i]);
                    }
                }

            }
        }

        SECTION("template_position=11") {
            string primers_path = "test/primers/y.primer3_with_ref_pos_11.out";
            ifstream file_handle(primers_path);
            PrimerFinder primer_finder(graph, &distance_index, file_handle, gbwt_graph, gbwt_index, r_index);

            SECTION("Loads the correct number of chromosomes") {
                REQUIRE(primer_finder.total_reference_paths() == 1);
            }

            SECTION("Loads the correct number of primer pairs") {
                REQUIRE(primer_finder.get_primer_pairs_of_chrom("y").size() == 5);
            }

            SECTION("Loads and processes the primers correctly") {
                primer_finder.add_primer_pair("y", 9, 14, 20, 22, 0, 20); // made up data, variation both at primers and in product
                primer_finder.add_primer_pair("y", 31, 0, 15, 34, 1, 15); // made up data, no variation at primers or in product

                // Correct primer attributes
                const vector<string> left_primers_sequences {
                    "TGCCTGGCATAGAGGAAAGC", "GAGTCGAGGCTCAAGGACAG", "CAGAGTCGAGGCTCAAGGAC",
                    "GAGGCTCAAGGACAGCTCTC", "TCCAGAAGCTGCTCTTTCCC", "AGCCAGACAAATCTGGGTTC",
                    "CAACTGGTAGTTACT"
                };

                const vector<size_t> left_primers_positions {
                    362, 620, 618, 625, 819, 181, 388
                };

                const vector<size_t> left_primers_lengths {
                    20, 20, 20, 20, 20, 20, 15
                };

                const vector<size_t> left_primers_nodes_count {
                    2, 1, 1, 2, 2, 6, 1
                };

                const vector<string> right_primers_sequences {
                    "GCCAGAAGAGCCTCAAGGAG", "AGGAGAGCTGGGAAAAGGGA", "AGGAGAGCTGGGAAAAGGGA",
                    "AGGAGAGCTGGGAAAAGGGA", "GCCTGGGTAGCTTTGGATGT", "AGATAATTAAACTGAAGTTC",
                    "GTTGACAATGAAAAG"
                };

                const vector<size_t> right_primers_positions {
                    466, 745, 745, 745, 935, 260, 485
                };

                const vector<size_t> right_primers_lengths {
                    20, 20, 20, 20, 20, 20, 15
                };

                const vector<size_t> right_primers_nodes_count {
                    2, 1, 1, 1, 2, 3, 1
                };

                const vector<size_t> min_product_sizes {
                    124, 142, 144, 137, 136, 99, 112
                };

                const vector<size_t> max_product_sizes {
                    124, 145, 147, 140, 137, 99, 112
                };

                const vector<size_t> linear_product_sizes {
                    124, 145, 147, 140, 136, 99, 112
                };
                
                const vector<double> variation_level {
                    1.0, 1.0, 1.0, 1.0, 1.0, 0.33333, 1.0
                };


                const vector<PrimerPair>& primer_pairs = primer_finder.get_primer_pairs_of_chrom("y");
                
                REQUIRE(primer_pairs.size() == left_primers_sequences.size());
                for (size_t i = 0; i < primer_pairs.size(); ++i) {
                    REQUIRE(left_primers_nodes_count[i]  == primer_pairs[i].left_primer.mapped_nodes_ids.size());
                    REQUIRE(left_primers_sequences[i]    == primer_pairs[i].left_primer.sequence);
                    REQUIRE(left_primers_positions[i]    == primer_pairs[i].left_primer.position_chromosome);
                    REQUIRE(left_primers_lengths[i]      == primer_pairs[i].left_primer.length);
                    REQUIRE(right_primers_nodes_count[i] == primer_pairs[i].right_primer.mapped_nodes_ids.size());
                    REQUIRE(right_primers_sequences[i]   == primer_pairs[i].right_primer.sequence);
                    REQUIRE(right_primers_positions[i]   == primer_pairs[i].right_primer.position_chromosome);
                    REQUIRE(right_primers_lengths[i]     == primer_pairs[i].right_primer.length);
                    REQUIRE(linear_product_sizes[i]      == primer_pairs[i].linear_product_size);
                    REQUIRE(min_product_sizes[i]         == primer_pairs[i].min_product_size);
                    REQUIRE(max_product_sizes[i]         == primer_pairs[i].max_product_size);
                    REQUIRE(abs(variation_level[i] - primer_pairs[i].variation_level) <= 0.0001);
                }

                SECTION("Check that primers are assigned with correct nodes") {
                    vector<size_t> pair_0_left_primer_nodes {27, 28};
                    for (size_t i = 0; i < primer_pairs[0].left_primer.mapped_nodes_ids.size(); i++) {
                        REQUIRE(primer_pairs[0].left_primer.mapped_nodes_ids[i] == pair_0_left_primer_nodes[i]);
                    }

                    vector<size_t> pair_0_right_primer_nodes {33, 34};
                    for (size_t i = 0; i < primer_pairs[0].right_primer.mapped_nodes_ids.size(); i++) {
                        REQUIRE(primer_pairs[0].right_primer.mapped_nodes_ids[i] == pair_0_right_primer_nodes[i]);
                    }

                    vector<size_t> pair_5_left_primer_nodes {9, 11, 12, 14, 15, 17};
                    for (size_t i = 0; i < primer_pairs[5].left_primer.mapped_nodes_ids.size(); i++) {
                        REQUIRE(primer_pairs[5].left_primer.mapped_nodes_ids[i] == pair_5_left_primer_nodes[i]);
                    }

                    vector<size_t> pair_5_right_primer_nodes {22, 24, 25};
                    for (size_t i = 0; i < primer_pairs[5].right_primer.mapped_nodes_ids.size(); i++) {
                        REQUIRE(primer_pairs[5].right_primer.mapped_nodes_ids[i] == pair_5_right_primer_nodes[i]);
                    }
                }

            }
        }
    }
}
}