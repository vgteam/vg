//
//  primers.cpp
//
//  Unit tests for primer filter 
//

#include <stdio.h>
#include <iostream>
#include <regex>
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
#include <bdsg/overlays/overlay_helper.hpp>
#include <gbwtgraph/index.h>
#include "xg.hpp"
#include "../primer_filter.hpp"
#include "../recombinator.hpp"
#include "../haplotype_indexer.hpp"

namespace vg {
namespace unittest {

using namespace std;


    TEST_CASE( "filter simple primers",
                "[primer_filter]" ) {
        
        // This is the same as test/primers/y.giraffe.gbz
        std::string graph_json = R"(
            {"edge": [{"from": "1", "to": "2"}, {"from": "2", "to": "3"}, {"from": "3", "to": "4"}, {"from": "3", "to": "5"}, {"from": "4", "to": "6"}, {"from": "5", "to": "6"}, {"from": "6", "to": "7"}, {"from": "7", "to": "8"}, {"from": "8", "to": "9"}, {"from": "9", "to": "10"}, {"from": "9", "to": "11"}, {"from": "10", "to": "12"}, {"from": "11", "to": "12"}, {"from": "12", "to": "13"}, {"from": "12", "to": "14"}, {"from": "13", "to": "15"}, {"from": "14", "to": "15"}, {"from": "15", "to": "16"}, {"from": "15", "to": "17"}, {"from": "16", "to": "17"}, {"from": "17", "to": "18"}, {"from": "17", "to": "19"}, {"from": "18", "to": "19"}, {"from": "19", "to": "20"}, {"from": "19", "to": "21"}, {"from": "20", "to": "21"}, {"from": "21", "to": "22"}, {"from": "22", "to": "23"}, {"from": "22", "to": "24"}, {"from": "23", "to": "25"}, {"from": "24", "to": "25"}, {"from": "25", "to": "26"}, {"from": "26", "to": "27"}, {"from": "27", "to": "28"}, {"from": "28", "to": "29"}, {"from": "28", "to": "30"}, {"from": "29", "to": "31"}, {"from": "30", "to": "31"}, {"from": "31", "to": "32"}, {"from": "32", "to": "33"}, {"from": "33", "to": "34"}, {"from": "34", "to": "35"}, {"from": "34", "to": "36"}, {"from": "35", "to": "36"}, {"from": "36", "to": "37"}, {"from": "37", "to": "38"}, {"from": "38", "to": "39"}, {"from": "39", "to": "40"}, {"from": "39", "to": "41"}, {"from": "40", "to": "42"}, {"from": "41", "to": "42"}, {"from": "42", "to": "43"}, {"from": "43", "to": "44"}, {"from": "44", "to": "45"}, {"from": "44", "to": "46"}, {"from": "45", "to": "46"}, {"from": "46", "to": "47"}, {"from": "47", "to": "48"}, {"from": "47", "to": "49"}, {"from": "48", "to": "49"}, {"from": "49", "to": "50"}, {"from": "49", "to": "51"}, {"from": "50", "to": "52"}, {"from": "51", "to": "53"}, {"from": "52", "to": "54"}, {"from": "53", "to": "54"}, {"from": "54", "to": "55"}, {"from": "55", "to": "56"}, {"from": "56", "to": "57"}, {"from": "57", "to": "58"}, {"from": "58", "to": "59"}, {"from": "58", "to": "60"}, {"from": "59", "to": "60"}, {"from": "60", "to": "61"}, {"from": "60", "to": "62"}, {"from": "61", "to": "62"}, {"from": "62", "to": "63"}, {"from": "62", "to": "64"}, {"from": "63", "to": "64"}, {"from": "64", "to": "65"}, {"from": "65", "to": "66"}], "node": [{"id": "1", "sequence": "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTA"}, {"id": "2", "sequence": "TTATATTCCAACTCTCTGGTTCCTGGTGCTAT"}, {"id": "3", "sequence": "GTGTAA"}, {"id": "4", "sequence": "T"}, {"id": "5", "sequence": "C"}, {"id": "6", "sequence": "TAGTAATGGTAATGGATATGTTGGGCTTTTTT"}, {"id": "7", "sequence": "CTTTGATTTATTTGAAGTGACGTTTGACAATC"}, {"id": "8", "sequence": "TATCACTAGGGGTAATGTGGGGAAATGGAAAG"}, {"id": "9", "sequence": "AATACAAGATTTGGA"}, {"id": "10", "sequence": "C"}, {"id": "11", "sequence": "G"}, {"id": "12", "sequence": "CCAGAC"}, {"id": "13", "sequence": "T"}, {"id": "14", "sequence": "A"}, {"id": "15", "sequence": "AA"}, {"id": "16", "sequence": "A"}, {"id": "17", "sequence": "TCTGGGTTCAAAT"}, {"id": "18", "sequence": "C"}, {"id": "19", "sequence": "CTCACTTTGCCACATATTAGC"}, {"id": "20", "sequence": "C"}, {"id": "21", "sequence": "ATGTGACTTTGAACAAGTTAGTTAATCTCTCT"}, {"id": "22", "sequence": "GAACTTCAGTTTAATTA"}, {"id": "23", "sequence": "C"}, {"id": "24", "sequence": "T"}, {"id": "25", "sequence": "CTCTAATATGGAGATGATACTACTGACAGCAG"}, {"id": "26", "sequence": "AGGTTTGCTGTGAAGATTAAATTAGGTGATGC"}, {"id": "27", "sequence": "TTGTAAAGCTCAGGGAATAGTGCCTGGCATAG"}, {"id": "28", "sequence": "AGGAAAGCCTCTG"}, {"id": "29", "sequence": "G"}, {"id": "30", "sequence": "A"}, {"id": "31", "sequence": "CAACTGGTAGTTACTGTTATTTACTATGAATC"}, {"id": "32", "sequence": "CTCACCTTCCTTGACTTCTTGAAACATTTGGC"}, {"id": "33", "sequence": "TATTGACCTCTTTCCTCCTTGAGGCTCTTCTG"}, {"id": "34", "sequence": "GCTTTTCATTGTCAACACAGTCAAC"}, {"id": "35", "sequence": "G"}, {"id": "36", "sequence": "CTCAATACAAGGGACATTAGGATTGGCAGTAG"}, {"id": "37", "sequence": "CTCAGAGATCTCTCTGCTCACCGTGATCTTCA"}, {"id": "38", "sequence": "AGTTTGAAAATTGCATCTCAAATCTAAGACCC"}, {"id": "39", "sequence": "AGA"}, {"id": "40", "sequence": "A"}, {"id": "41", "sequence": "G"}, {"id": "42", "sequence": "GGCTCACCCAGAGTCGAGGCTCAAGGACAGCT"}, {"id": "43", "sequence": "CTCCTTTGTGTCCAGAGTGTATACGATGTAAC"}, {"id": "44", "sequence": "TCTG"}, {"id": "45", "sequence": "TT"}, {"id": "46", "sequence": "CGGGCACTGGTGAAAGATAACAGAGGAAATGC"}, {"id": "47", "sequence": "CTGGCTTTTTATCAGA"}, {"id": "48", "sequence": "A"}, {"id": "49", "sequence": "CAT"}, {"id": "50", "sequence": "A"}, {"id": "51", "sequence": "G"}, {"id": "52", "sequence": "C"}, {"id": "53", "sequence": "T"}, {"id": "54", "sequence": "TTCCAAGCTTATCCCTTTTCCCAGCTCTCCTT"}, {"id": "55", "sequence": "GTCCCTCCCAAGATCTCTTCACTGGCCTCTTA"}, {"id": "56", "sequence": "TCTTTACTGTTACCAAATCTTTCCAGAAGCTG"}, {"id": "57", "sequence": "CTCTTTCCCTCAATTGTTCATTTGTCTTCTTG"}, {"id": "58", "sequence": "TCCAGGAATGAACCACT"}, {"id": "59", "sequence": "G"}, {"id": "60", "sequence": "GCTCTCTTCTTGTCAGATCAGCTT"}, {"id": "61", "sequence": "A"}, {"id": "62", "sequence": "CTCATCCC"}, {"id": "63", "sequence": "T"}, {"id": "64", "sequence": "CCTCAAGGGCCTTTAACTACTCCACATCCAAA"}, {"id": "65", "sequence": "GCTACCCAGGCCATTTTAAGTTTCCTGTGGAC"}, {"id": "66", "sequence": "TAAGGACAAAGGTGCGGGGAGATGA"}], "path": [{"mapping": [{"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "1"}, "rank": "1"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "2"}, "rank": "2"}, {"edit": [{"from_length": 6, "to_length": 6}], "position": {"node_id": "3"}, "rank": "3"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "4"}, "rank": "4"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "6"}, "rank": "5"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "7"}, "rank": "6"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "8"}, "rank": "7"}, {"edit": [{"from_length": 15, "to_length": 15}], "position": {"node_id": "9"}, "rank": "8"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "10"}, "rank": "9"}, {"edit": [{"from_length": 6, "to_length": 6}], "position": {"node_id": "12"}, "rank": "10"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "13"}, "rank": "11"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "15"}, "rank": "12"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "16"}, "rank": "13"}, {"edit": [{"from_length": 13, "to_length": 13}], "position": {"node_id": "17"}, "rank": "14"}, {"edit": [{"from_length": 21, "to_length": 21}], "position": {"node_id": "19"}, "rank": "15"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "21"}, "rank": "16"}, {"edit": [{"from_length": 17, "to_length": 17}], "position": {"node_id": "22"}, "rank": "17"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "23"}, "rank": "18"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "25"}, "rank": "19"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "26"}, "rank": "20"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "27"}, "rank": "21"}, {"edit": [{"from_length": 13, "to_length": 13}], "position": {"node_id": "28"}, "rank": "22"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "29"}, "rank": "23"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "31"}, "rank": "24"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "32"}, "rank": "25"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "33"}, "rank": "26"}, {"edit": [{"from_length": 25, "to_length": 25}], "position": {"node_id": "34"}, "rank": "27"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "36"}, "rank": "28"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "37"}, "rank": "29"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "38"}, "rank": "30"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "39"}, "rank": "31"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "40"}, "rank": "32"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "42"}, "rank": "33"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "43"}, "rank": "34"}, {"edit": [{"from_length": 4, "to_length": 4}], "position": {"node_id": "44"}, "rank": "35"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "46"}, "rank": "36"}, {"edit": [{"from_length": 16, "to_length": 16}], "position": {"node_id": "47"}, "rank": "37"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "49"}, "rank": "38"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "50"}, "rank": "39"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "52"}, "rank": "40"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "54"}, "rank": "41"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "55"}, "rank": "42"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "56"}, "rank": "43"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "57"}, "rank": "44"}, {"edit": [{"from_length": 17, "to_length": 17}], "position": {"node_id": "58"}, "rank": "45"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "59"}, "rank": "46"}, {"edit": [{"from_length": 24, "to_length": 24}], "position": {"node_id": "60"}, "rank": "47"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "61"}, "rank": "48"}, {"edit": [{"from_length": 8, "to_length": 8}], "position": {"node_id": "62"}, "rank": "49"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "64"}, "rank": "50"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "65"}, "rank": "51"}, {"edit": [{"from_length": 25, "to_length": 25}], "position": {"node_id": "66"}, "rank": "52"}], "name": "1#0#y#0"}, {"mapping": [{"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "1"}, "rank": "1"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "2"}, "rank": "2"}, {"edit": [{"from_length": 6, "to_length": 6}], "position": {"node_id": "3"}, "rank": "3"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "4"}, "rank": "4"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "6"}, "rank": "5"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "7"}, "rank": "6"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "8"}, "rank": "7"}, {"edit": [{"from_length": 15, "to_length": 15}], "position": {"node_id": "9"}, "rank": "8"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "11"}, "rank": "9"}, {"edit": [{"from_length": 6, "to_length": 6}], "position": {"node_id": "12"}, "rank": "10"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "13"}, "rank": "11"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "15"}, "rank": "12"}, {"edit": [{"from_length": 13, "to_length": 13}], "position": {"node_id": "17"}, "rank": "13"}, {"edit": [{"from_length": 21, "to_length": 21}], "position": {"node_id": "19"}, "rank": "14"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "21"}, "rank": "15"}, {"edit": [{"from_length": 17, "to_length": 17}], "position": {"node_id": "22"}, "rank": "16"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "23"}, "rank": "17"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "25"}, "rank": "18"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "26"}, "rank": "19"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "27"}, "rank": "20"}, {"edit": [{"from_length": 13, "to_length": 13}], "position": {"node_id": "28"}, "rank": "21"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "29"}, "rank": "22"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "31"}, "rank": "23"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "32"}, "rank": "24"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "33"}, "rank": "25"}, {"edit": [{"from_length": 25, "to_length": 25}], "position": {"node_id": "34"}, "rank": "26"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "36"}, "rank": "27"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "37"}, "rank": "28"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "38"}, "rank": "29"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "39"}, "rank": "30"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "41"}, "rank": "31"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "42"}, "rank": "32"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "43"}, "rank": "33"}, {"edit": [{"from_length": 4, "to_length": 4}], "position": {"node_id": "44"}, "rank": "34"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "45"}, "rank": "35"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "46"}, "rank": "36"}, {"edit": [{"from_length": 16, "to_length": 16}], "position": {"node_id": "47"}, "rank": "37"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "49"}, "rank": "38"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "50"}, "rank": "39"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "52"}, "rank": "40"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "54"}, "rank": "41"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "55"}, "rank": "42"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "56"}, "rank": "43"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "57"}, "rank": "44"}, {"edit": [{"from_length": 17, "to_length": 17}], "position": {"node_id": "58"}, "rank": "45"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "59"}, "rank": "46"}, {"edit": [{"from_length": 24, "to_length": 24}], "position": {"node_id": "60"}, "rank": "47"}, {"edit": [{"from_length": 8, "to_length": 8}], "position": {"node_id": "62"}, "rank": "48"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "64"}, "rank": "49"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "65"}, "rank": "50"}, {"edit": [{"from_length": 25, "to_length": 25}], "position": {"node_id": "66"}, "rank": "51"}], "name": "1#1#y#0"}, {"mapping": [{"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "1"}, "rank": "1"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "2"}, "rank": "2"}, {"edit": [{"from_length": 6, "to_length": 6}], "position": {"node_id": "3"}, "rank": "3"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "5"}, "rank": "4"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "6"}, "rank": "5"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "7"}, "rank": "6"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "8"}, "rank": "7"}, {"edit": [{"from_length": 15, "to_length": 15}], "position": {"node_id": "9"}, "rank": "8"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "11"}, "rank": "9"}, {"edit": [{"from_length": 6, "to_length": 6}], "position": {"node_id": "12"}, "rank": "10"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "14"}, "rank": "11"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "15"}, "rank": "12"}, {"edit": [{"from_length": 13, "to_length": 13}], "position": {"node_id": "17"}, "rank": "13"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "18"}, "rank": "14"}, {"edit": [{"from_length": 21, "to_length": 21}], "position": {"node_id": "19"}, "rank": "15"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "20"}, "rank": "16"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "21"}, "rank": "17"}, {"edit": [{"from_length": 17, "to_length": 17}], "position": {"node_id": "22"}, "rank": "18"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "24"}, "rank": "19"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "25"}, "rank": "20"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "26"}, "rank": "21"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "27"}, "rank": "22"}, {"edit": [{"from_length": 13, "to_length": 13}], "position": {"node_id": "28"}, "rank": "23"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "30"}, "rank": "24"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "31"}, "rank": "25"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "32"}, "rank": "26"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "33"}, "rank": "27"}, {"edit": [{"from_length": 25, "to_length": 25}], "position": {"node_id": "34"}, "rank": "28"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "35"}, "rank": "29"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "36"}, "rank": "30"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "37"}, "rank": "31"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "38"}, "rank": "32"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "39"}, "rank": "33"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "41"}, "rank": "34"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "42"}, "rank": "35"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "43"}, "rank": "36"}, {"edit": [{"from_length": 4, "to_length": 4}], "position": {"node_id": "44"}, "rank": "37"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "45"}, "rank": "38"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "46"}, "rank": "39"}, {"edit": [{"from_length": 16, "to_length": 16}], "position": {"node_id": "47"}, "rank": "40"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "48"}, "rank": "41"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "49"}, "rank": "42"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "51"}, "rank": "43"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "53"}, "rank": "44"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "54"}, "rank": "45"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "55"}, "rank": "46"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "56"}, "rank": "47"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "57"}, "rank": "48"}, {"edit": [{"from_length": 17, "to_length": 17}], "position": {"node_id": "58"}, "rank": "49"}, {"edit": [{"from_length": 24, "to_length": 24}], "position": {"node_id": "60"}, "rank": "50"}, {"edit": [{"from_length": 8, "to_length": 8}], "position": {"node_id": "62"}, "rank": "51"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "63"}, "rank": "52"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "64"}, "rank": "53"}, {"edit": [{"from_length": 32, "to_length": 32}], "position": {"node_id": "65"}, "rank": "54"}, {"edit": [{"from_length": 25, "to_length": 25}], "position": {"node_id": "66"}, "rank": "55"}], "name": "y"}]}
        )";

        // Load graph with haplotypes from JSON
        VG vg_graph;
        Graph chunk;
        json2pb(chunk, graph_json.c_str(), graph_json.size());
        // Need to use extend() to bring along the paths.
        vg_graph.extend(chunk);

        REQUIRE(vg_graph.has_path("y"));

        // Make PathPositionHandleGraph
        bdsg::PathPositionOverlayHelper overlay_helper;
        PathPositionHandleGraph* graph = overlay_helper.apply(&vg_graph);

        REQUIRE(graph->has_path("y"));

        // Make GBWT
        // First build it uncompressed
        HaplotypeIndexer haplotype_indexer;
        std::unique_ptr<gbwt::DynamicGBWT> dynamic_gbwt = haplotype_indexer.build_gbwt(*graph);
        // Then compress it
        gbwt::GBWT gbwt_index(*dynamic_gbwt);
        dynamic_gbwt.reset();
        
        // Make GBWTGraph
        gbwtgraph::GBWTGraph gbwt_graph(gbwt_index, *graph);
        
        // Make R index
        gbwt::FastLocate r_index(gbwt_index);
        
        // Make distance index
        IntegratedSnarlFinder snarl_finder(*graph);
        SnarlDistanceIndex distance_index;
        fill_in_distance_index(&distance_index, graph, &snarl_finder);

        SECTION("template_position=0") {
            std::string primer3_with_ref_pos = R"(SEQUENCE_ID=y|gene|feature|0
SEQUENCE_TEMPLATE=CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAATGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGAACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTACTGACAGCAGAGGTTTGCTGTGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTTACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCTTCTGGCTTTTCATTGTCAACACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCTCACCGTGATCTTCAAGTTTGAAAATTGCATCTCAAATCTAAGACCCAGAGGGCTCACCCAGAGTCGAGGCTCAAGGACAGCTCTCCTTTGTGTCCAGAGTGTATACGATGTAACTCTGTTCGGGCACTGGTGAAAGATAACAGAGGAAATGCCTGGCTTTTTATCAGAACATGTTTCCAAGCTTATCCCTTTTCCCAGCTCTCCTTGTCCCTCCCAAGATCTCTTCACTGGCCTCTTATCTTTACTGTTACCAAATCTTTCCAGAAGCTGCTCTTTCCCTCAATTGTTCATTTGTCTTCTTGTCCAGGAATGAACCACTGCTCTCTTCTTGTCAGATCAGCTTCTCATCCCTCCTCAAGGGCCTTTAACTACTCCACATCCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAGATGA
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=75-150
PRIMER_EXPLAIN_FLAG=1
PRIMER_LEFT_EXPLAIN=considered 4586, GC content failed 41, low tm 3329, high tm 93, high hairpin stability 6, long poly-x seq 15, ok 1102
PRIMER_RIGHT_EXPLAIN=considered 4585, GC content failed 40, low tm 3257, high tm 106, long poly-x seq 15, ok 1167
PRIMER_PAIR_EXPLAIN=considered 106, unacceptable product size 101, ok 5
PRIMER_LEFT_NUM_RETURNED=5
PRIMER_RIGHT_NUM_RETURNED=5
PRIMER_INTERNAL_NUM_RETURNED=0
PRIMER_PAIR_NUM_RETURNED=5
PRIMER_PAIR_0_PENALTY=0.214768
PRIMER_LEFT_0_PENALTY=0.107017
PRIMER_RIGHT_0_PENALTY=0.107752
PRIMER_LEFT_0_SEQUENCE=TGCCTGGCATAGAGGAAAGC
PRIMER_RIGHT_0_SEQUENCE=GCCAGAAGAGCCTCAAGGAG
PRIMER_LEFT_0=362,20
PRIMER_RIGHT_0=485,20
PRIMER_LEFT_0_TM=60.107
PRIMER_RIGHT_0_TM=60.108
PRIMER_LEFT_0_GC_PERCENT=55.000
PRIMER_RIGHT_0_GC_PERCENT=60.000
PRIMER_LEFT_0_SELF_ANY_TH=19.30
PRIMER_RIGHT_0_SELF_ANY_TH=0.00
PRIMER_LEFT_0_SELF_END_TH=0.00
PRIMER_RIGHT_0_SELF_END_TH=0.00
PRIMER_LEFT_0_HAIRPIN_TH=31.71
PRIMER_RIGHT_0_HAIRPIN_TH=37.39
PRIMER_LEFT_0_END_STABILITY=3.5100
PRIMER_RIGHT_0_END_STABILITY=3.6900
PRIMER_PAIR_0_COMPL_ANY_TH=6.57
PRIMER_PAIR_0_COMPL_END_TH=4.13
PRIMER_PAIR_0_PRODUCT_SIZE=124
PRIMER_PAIR_0_PRODUCT_TM=81.8
PRIMER_PAIR_1_PENALTY=0.351214
PRIMER_LEFT_1_PENALTY=0.172352
PRIMER_RIGHT_1_PENALTY=0.178861
PRIMER_LEFT_1_SEQUENCE=GAGTCGAGGCTCAAGGACAG
PRIMER_RIGHT_1_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_1=620,20
PRIMER_RIGHT_1=764,20
PRIMER_LEFT_1_TM=59.828
PRIMER_RIGHT_1_TM=60.179
PRIMER_LEFT_1_GC_PERCENT=60.000
PRIMER_RIGHT_1_GC_PERCENT=55.000
PRIMER_LEFT_1_SELF_ANY_TH=15.95
PRIMER_RIGHT_1_SELF_ANY_TH=0.00
PRIMER_LEFT_1_SELF_END_TH=0.00
PRIMER_RIGHT_1_SELF_END_TH=0.00
PRIMER_LEFT_1_HAIRPIN_TH=35.67
PRIMER_RIGHT_1_HAIRPIN_TH=0.00
PRIMER_LEFT_1_END_STABILITY=3.5100
PRIMER_RIGHT_1_END_STABILITY=4.2000
PRIMER_PAIR_1_COMPL_ANY_TH=0.00
PRIMER_PAIR_1_COMPL_END_TH=0.00
PRIMER_PAIR_1_PRODUCT_SIZE=145
PRIMER_PAIR_1_PRODUCT_TM=83.5
PRIMER_PAIR_2_PENALTY=0.351214
PRIMER_LEFT_2_PENALTY=0.172352
PRIMER_RIGHT_2_PENALTY=0.178861
PRIMER_LEFT_2_SEQUENCE=CAGAGTCGAGGCTCAAGGAC
PRIMER_RIGHT_2_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_2=618,20
PRIMER_RIGHT_2=764,20
PRIMER_LEFT_2_TM=59.828
PRIMER_RIGHT_2_TM=60.179
PRIMER_LEFT_2_GC_PERCENT=60.000
PRIMER_RIGHT_2_GC_PERCENT=55.000
PRIMER_LEFT_2_SELF_ANY_TH=15.95
PRIMER_RIGHT_2_SELF_ANY_TH=0.00
PRIMER_LEFT_2_SELF_END_TH=13.47
PRIMER_RIGHT_2_SELF_END_TH=0.00
PRIMER_LEFT_2_HAIRPIN_TH=37.94
PRIMER_RIGHT_2_HAIRPIN_TH=0.00
PRIMER_LEFT_2_END_STABILITY=3.8500
PRIMER_RIGHT_2_END_STABILITY=4.2000
PRIMER_PAIR_2_COMPL_ANY_TH=0.00
PRIMER_PAIR_2_COMPL_END_TH=0.00
PRIMER_PAIR_2_PRODUCT_SIZE=147
PRIMER_PAIR_2_PRODUCT_TM=83.6
PRIMER_PAIR_3_PENALTY=0.354392
PRIMER_LEFT_3_PENALTY=0.175531
PRIMER_RIGHT_3_PENALTY=0.178861
PRIMER_LEFT_3_SEQUENCE=GAGGCTCAAGGACAGCTCTC
PRIMER_RIGHT_3_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_3=625,20
PRIMER_RIGHT_3=764,20
PRIMER_LEFT_3_TM=59.824
PRIMER_RIGHT_3_TM=60.179
PRIMER_LEFT_3_GC_PERCENT=60.000
PRIMER_RIGHT_3_GC_PERCENT=55.000
PRIMER_LEFT_3_SELF_ANY_TH=10.25
PRIMER_RIGHT_3_SELF_ANY_TH=0.00
PRIMER_LEFT_3_SELF_END_TH=0.00
PRIMER_RIGHT_3_SELF_END_TH=0.00
PRIMER_LEFT_3_HAIRPIN_TH=37.05
PRIMER_RIGHT_3_HAIRPIN_TH=0.00
PRIMER_LEFT_3_END_STABILITY=3.2000
PRIMER_RIGHT_3_END_STABILITY=4.2000
PRIMER_PAIR_3_COMPL_ANY_TH=26.57
PRIMER_PAIR_3_COMPL_END_TH=26.57
PRIMER_PAIR_3_PRODUCT_SIZE=140
PRIMER_PAIR_3_PRODUCT_TM=83.2
PRIMER_PAIR_4_PENALTY=0.360353
PRIMER_LEFT_4_PENALTY=0.326264
PRIMER_RIGHT_4_PENALTY=0.034089
PRIMER_LEFT_4_SEQUENCE=TCCAGAAGCTGCTCTTTCCC
PRIMER_RIGHT_4_SEQUENCE=GCCTGGGTAGCTTTGGATGT
PRIMER_LEFT_4=819,20
PRIMER_RIGHT_4=954,20
PRIMER_LEFT_4_TM=59.674
PRIMER_RIGHT_4_TM=60.034
PRIMER_LEFT_4_GC_PERCENT=55.000
PRIMER_RIGHT_4_GC_PERCENT=55.000
PRIMER_LEFT_4_SELF_ANY_TH=1.00
PRIMER_RIGHT_4_SELF_ANY_TH=0.00
PRIMER_LEFT_4_SELF_END_TH=0.00
PRIMER_RIGHT_4_SELF_END_TH=0.00
PRIMER_LEFT_4_HAIRPIN_TH=34.56
PRIMER_RIGHT_4_HAIRPIN_TH=0.00
PRIMER_LEFT_4_END_STABILITY=3.9700
PRIMER_RIGHT_4_END_STABILITY=3.0600
PRIMER_PAIR_4_COMPL_ANY_TH=13.72
PRIMER_PAIR_4_COMPL_END_TH=10.53
PRIMER_PAIR_4_PRODUCT_SIZE=136
PRIMER_PAIR_4_PRODUCT_TM=83.6
=
)";
            std::stringstream file_handle(primer3_with_ref_pos);
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
            std::string primer3_with_ref_pos_11 = R"(SEQUENCE_ID=y|gene|feature|11
SEQUENCE_TEMPLATE=TGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAATGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGAACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTACTGACAGCAGAGGTTTGCTGTGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTTACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCTTCTGGCTTTTCATTGTCAACACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCTCACCGTGATCTTCAAGTTTGAAAATTGCATCTCAAATCTAAGACCCAGAGGGCTCACCCAGAGTCGAGGCTCAAGGACAGCTCTCCTTTGTGTCCAGAGTGTATACGATGTAACTCTGTTCGGGCACTGGTGAAAGATAACAGAGGAAATGCCTGGCTTTTTATCAGAACATGTTTCCAAGCTTATCCCTTTTCCCAGCTCTCCTTGTCCCTCCCAAGATCTCTTCACTGGCCTCTTATCTTTACTGTTACCAAATCTTTCCAGAAGCTGCTCTTTCCCTCAATTGTTCATTTGTCTTCTTGTCCAGGAATGAACCACTGCTCTCTTCTTGTCAGATCAGCTTCTCATCCCTCCTCAAGGGCCTTTAACTACTCCACATCCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAGATGA
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=75-150
PRIMER_EXPLAIN_FLAG=1
PRIMER_LEFT_EXPLAIN=considered 4531, GC content failed 41, low tm 3277, high tm 93, high hairpin stability 6, long poly-x seq 15, ok 1099
PRIMER_RIGHT_EXPLAIN=considered 4530, GC content failed 40, low tm 3203, high tm 106, long poly-x seq 15, ok 1166
PRIMER_PAIR_EXPLAIN=considered 106, unacceptable product size 101, ok 5
PRIMER_LEFT_NUM_RETURNED=5
PRIMER_RIGHT_NUM_RETURNED=5
PRIMER_INTERNAL_NUM_RETURNED=0
PRIMER_PAIR_NUM_RETURNED=5
PRIMER_PAIR_0_PENALTY=0.214768
PRIMER_LEFT_0_PENALTY=0.107017
PRIMER_RIGHT_0_PENALTY=0.107752
PRIMER_LEFT_0_SEQUENCE=TGCCTGGCATAGAGGAAAGC
PRIMER_RIGHT_0_SEQUENCE=GCCAGAAGAGCCTCAAGGAG
PRIMER_LEFT_0=351,20
PRIMER_RIGHT_0=474,20
PRIMER_LEFT_0_TM=60.107
PRIMER_RIGHT_0_TM=60.108
PRIMER_LEFT_0_GC_PERCENT=55.000
PRIMER_RIGHT_0_GC_PERCENT=60.000
PRIMER_LEFT_0_SELF_ANY_TH=19.30
PRIMER_RIGHT_0_SELF_ANY_TH=0.00
PRIMER_LEFT_0_SELF_END_TH=0.00
PRIMER_RIGHT_0_SELF_END_TH=0.00
PRIMER_LEFT_0_HAIRPIN_TH=31.71
PRIMER_RIGHT_0_HAIRPIN_TH=37.39
PRIMER_LEFT_0_END_STABILITY=3.5100
PRIMER_RIGHT_0_END_STABILITY=3.6900
PRIMER_PAIR_0_COMPL_ANY_TH=6.57
PRIMER_PAIR_0_COMPL_END_TH=4.13
PRIMER_PAIR_0_PRODUCT_SIZE=124
PRIMER_PAIR_0_PRODUCT_TM=81.8
PRIMER_PAIR_1_PENALTY=0.351214
PRIMER_LEFT_1_PENALTY=0.172352
PRIMER_RIGHT_1_PENALTY=0.178861
PRIMER_LEFT_1_SEQUENCE=GAGTCGAGGCTCAAGGACAG
PRIMER_RIGHT_1_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_1=609,20
PRIMER_RIGHT_1=753,20
PRIMER_LEFT_1_TM=59.828
PRIMER_RIGHT_1_TM=60.179
PRIMER_LEFT_1_GC_PERCENT=60.000
PRIMER_RIGHT_1_GC_PERCENT=55.000
PRIMER_LEFT_1_SELF_ANY_TH=15.95
PRIMER_RIGHT_1_SELF_ANY_TH=0.00
PRIMER_LEFT_1_SELF_END_TH=0.00
PRIMER_RIGHT_1_SELF_END_TH=0.00
PRIMER_LEFT_1_HAIRPIN_TH=35.67
PRIMER_RIGHT_1_HAIRPIN_TH=0.00
PRIMER_LEFT_1_END_STABILITY=3.5100
PRIMER_RIGHT_1_END_STABILITY=4.2000
PRIMER_PAIR_1_COMPL_ANY_TH=0.00
PRIMER_PAIR_1_COMPL_END_TH=0.00
PRIMER_PAIR_1_PRODUCT_SIZE=145
PRIMER_PAIR_1_PRODUCT_TM=83.5
PRIMER_PAIR_2_PENALTY=0.351214
PRIMER_LEFT_2_PENALTY=0.172352
PRIMER_RIGHT_2_PENALTY=0.178861
PRIMER_LEFT_2_SEQUENCE=CAGAGTCGAGGCTCAAGGAC
PRIMER_RIGHT_2_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_2=607,20
PRIMER_RIGHT_2=753,20
PRIMER_LEFT_2_TM=59.828
PRIMER_RIGHT_2_TM=60.179
PRIMER_LEFT_2_GC_PERCENT=60.000
PRIMER_RIGHT_2_GC_PERCENT=55.000
PRIMER_LEFT_2_SELF_ANY_TH=15.95
PRIMER_RIGHT_2_SELF_ANY_TH=0.00
PRIMER_LEFT_2_SELF_END_TH=13.47
PRIMER_RIGHT_2_SELF_END_TH=0.00
PRIMER_LEFT_2_HAIRPIN_TH=37.94
PRIMER_RIGHT_2_HAIRPIN_TH=0.00
PRIMER_LEFT_2_END_STABILITY=3.8500
PRIMER_RIGHT_2_END_STABILITY=4.2000
PRIMER_PAIR_2_COMPL_ANY_TH=0.00
PRIMER_PAIR_2_COMPL_END_TH=0.00
PRIMER_PAIR_2_PRODUCT_SIZE=147
PRIMER_PAIR_2_PRODUCT_TM=83.6
PRIMER_PAIR_3_PENALTY=0.354392
PRIMER_LEFT_3_PENALTY=0.175531
PRIMER_RIGHT_3_PENALTY=0.178861
PRIMER_LEFT_3_SEQUENCE=GAGGCTCAAGGACAGCTCTC
PRIMER_RIGHT_3_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_3=614,20
PRIMER_RIGHT_3=753,20
PRIMER_LEFT_3_TM=59.824
PRIMER_RIGHT_3_TM=60.179
PRIMER_LEFT_3_GC_PERCENT=60.000
PRIMER_RIGHT_3_GC_PERCENT=55.000
PRIMER_LEFT_3_SELF_ANY_TH=10.25
PRIMER_RIGHT_3_SELF_ANY_TH=0.00
PRIMER_LEFT_3_SELF_END_TH=0.00
PRIMER_RIGHT_3_SELF_END_TH=0.00
PRIMER_LEFT_3_HAIRPIN_TH=37.05
PRIMER_RIGHT_3_HAIRPIN_TH=0.00
PRIMER_LEFT_3_END_STABILITY=3.2000
PRIMER_RIGHT_3_END_STABILITY=4.2000
PRIMER_PAIR_3_COMPL_ANY_TH=26.57
PRIMER_PAIR_3_COMPL_END_TH=26.57
PRIMER_PAIR_3_PRODUCT_SIZE=140
PRIMER_PAIR_3_PRODUCT_TM=83.2
PRIMER_PAIR_4_PENALTY=0.360353
PRIMER_LEFT_4_PENALTY=0.326264
PRIMER_RIGHT_4_PENALTY=0.034089
PRIMER_LEFT_4_SEQUENCE=TCCAGAAGCTGCTCTTTCCC
PRIMER_RIGHT_4_SEQUENCE=GCCTGGGTAGCTTTGGATGT
PRIMER_LEFT_4=808,20
PRIMER_RIGHT_4=943,20
PRIMER_LEFT_4_TM=59.674
PRIMER_RIGHT_4_TM=60.034
PRIMER_LEFT_4_GC_PERCENT=55.000
PRIMER_RIGHT_4_GC_PERCENT=55.000
PRIMER_LEFT_4_SELF_ANY_TH=1.00
PRIMER_RIGHT_4_SELF_ANY_TH=0.00
PRIMER_LEFT_4_SELF_END_TH=0.00
PRIMER_RIGHT_4_SELF_END_TH=0.00
PRIMER_LEFT_4_HAIRPIN_TH=34.56
PRIMER_RIGHT_4_HAIRPIN_TH=0.00
PRIMER_LEFT_4_END_STABILITY=3.9700
PRIMER_RIGHT_4_END_STABILITY=3.0600
PRIMER_PAIR_4_COMPL_ANY_TH=13.72
PRIMER_PAIR_4_COMPL_END_TH=10.53
PRIMER_PAIR_4_PRODUCT_SIZE=136
PRIMER_PAIR_4_PRODUCT_TM=83.6
=
)";
            
            std::stringstream file_handle(primer3_with_ref_pos_11);
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
        SECTION("template_position=11, no path name") {
            string primers_path = "primers/y.primer3_with_ref_pos_11.nopath.out";
            std::string primer3_with_ref_pos_11_nopath = R"(SEQUENCE_ID=x|gene|feature|11
SEQUENCE_TEMPLATE=TGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAATGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGAACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTACTGACAGCAGAGGTTTGCTGTGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTTACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCTTCTGGCTTTTCATTGTCAACACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCTCACCGTGATCTTCAAGTTTGAAAATTGCATCTCAAATCTAAGACCCAGAGGGCTCACCCAGAGTCGAGGCTCAAGGACAGCTCTCCTTTGTGTCCAGAGTGTATACGATGTAACTCTGTTCGGGCACTGGTGAAAGATAACAGAGGAAATGCCTGGCTTTTTATCAGAACATGTTTCCAAGCTTATCCCTTTTCCCAGCTCTCCTTGTCCCTCCCAAGATCTCTTCACTGGCCTCTTATCTTTACTGTTACCAAATCTTTCCAGAAGCTGCTCTTTCCCTCAATTGTTCATTTGTCTTCTTGTCCAGGAATGAACCACTGCTCTCTTCTTGTCAGATCAGCTTCTCATCCCTCCTCAAGGGCCTTTAACTACTCCACATCCAAAGCTACCCAGGCCATTTTAAGTTTCCTGTGGACTAAGGACAAAGGTGCGGGGAGATGA
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=75-150
PRIMER_EXPLAIN_FLAG=1
PRIMER_LEFT_EXPLAIN=considered 4531, GC content failed 41, low tm 3277, high tm 93, high hairpin stability 6, long poly-x seq 15, ok 1099
PRIMER_RIGHT_EXPLAIN=considered 4530, GC content failed 40, low tm 3203, high tm 106, long poly-x seq 15, ok 1166
PRIMER_PAIR_EXPLAIN=considered 106, unacceptable product size 101, ok 5
PRIMER_LEFT_NUM_RETURNED=5
PRIMER_RIGHT_NUM_RETURNED=5
PRIMER_INTERNAL_NUM_RETURNED=0
PRIMER_PAIR_NUM_RETURNED=5
PRIMER_PAIR_0_PENALTY=0.214768
PRIMER_LEFT_0_PENALTY=0.107017
PRIMER_RIGHT_0_PENALTY=0.107752
PRIMER_LEFT_0_SEQUENCE=TGCCTGGCATAGAGGAAAGC
PRIMER_RIGHT_0_SEQUENCE=GCCAGAAGAGCCTCAAGGAG
PRIMER_LEFT_0=351,20
PRIMER_RIGHT_0=474,20
PRIMER_LEFT_0_TM=60.107
PRIMER_RIGHT_0_TM=60.108
PRIMER_LEFT_0_GC_PERCENT=55.000
PRIMER_RIGHT_0_GC_PERCENT=60.000
PRIMER_LEFT_0_SELF_ANY_TH=19.30
PRIMER_RIGHT_0_SELF_ANY_TH=0.00
PRIMER_LEFT_0_SELF_END_TH=0.00
PRIMER_RIGHT_0_SELF_END_TH=0.00
PRIMER_LEFT_0_HAIRPIN_TH=31.71
PRIMER_RIGHT_0_HAIRPIN_TH=37.39
PRIMER_LEFT_0_END_STABILITY=3.5100
PRIMER_RIGHT_0_END_STABILITY=3.6900
PRIMER_PAIR_0_COMPL_ANY_TH=6.57
PRIMER_PAIR_0_COMPL_END_TH=4.13
PRIMER_PAIR_0_PRODUCT_SIZE=124
PRIMER_PAIR_0_PRODUCT_TM=81.8
PRIMER_PAIR_1_PENALTY=0.351214
PRIMER_LEFT_1_PENALTY=0.172352
PRIMER_RIGHT_1_PENALTY=0.178861
PRIMER_LEFT_1_SEQUENCE=GAGTCGAGGCTCAAGGACAG
PRIMER_RIGHT_1_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_1=609,20
PRIMER_RIGHT_1=753,20
PRIMER_LEFT_1_TM=59.828
PRIMER_RIGHT_1_TM=60.179
PRIMER_LEFT_1_GC_PERCENT=60.000
PRIMER_RIGHT_1_GC_PERCENT=55.000
PRIMER_LEFT_1_SELF_ANY_TH=15.95
PRIMER_RIGHT_1_SELF_ANY_TH=0.00
PRIMER_LEFT_1_SELF_END_TH=0.00
PRIMER_RIGHT_1_SELF_END_TH=0.00
PRIMER_LEFT_1_HAIRPIN_TH=35.67
PRIMER_RIGHT_1_HAIRPIN_TH=0.00
PRIMER_LEFT_1_END_STABILITY=3.5100
PRIMER_RIGHT_1_END_STABILITY=4.2000
PRIMER_PAIR_1_COMPL_ANY_TH=0.00
PRIMER_PAIR_1_COMPL_END_TH=0.00
PRIMER_PAIR_1_PRODUCT_SIZE=145
PRIMER_PAIR_1_PRODUCT_TM=83.5
PRIMER_PAIR_2_PENALTY=0.351214
PRIMER_LEFT_2_PENALTY=0.172352
PRIMER_RIGHT_2_PENALTY=0.178861
PRIMER_LEFT_2_SEQUENCE=CAGAGTCGAGGCTCAAGGAC
PRIMER_RIGHT_2_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_2=607,20
PRIMER_RIGHT_2=753,20
PRIMER_LEFT_2_TM=59.828
PRIMER_RIGHT_2_TM=60.179
PRIMER_LEFT_2_GC_PERCENT=60.000
PRIMER_RIGHT_2_GC_PERCENT=55.000
PRIMER_LEFT_2_SELF_ANY_TH=15.95
PRIMER_RIGHT_2_SELF_ANY_TH=0.00
PRIMER_LEFT_2_SELF_END_TH=13.47
PRIMER_RIGHT_2_SELF_END_TH=0.00
PRIMER_LEFT_2_HAIRPIN_TH=37.94
PRIMER_RIGHT_2_HAIRPIN_TH=0.00
PRIMER_LEFT_2_END_STABILITY=3.8500
PRIMER_RIGHT_2_END_STABILITY=4.2000
PRIMER_PAIR_2_COMPL_ANY_TH=0.00
PRIMER_PAIR_2_COMPL_END_TH=0.00
PRIMER_PAIR_2_PRODUCT_SIZE=147
PRIMER_PAIR_2_PRODUCT_TM=83.6
PRIMER_PAIR_3_PENALTY=0.354392
PRIMER_LEFT_3_PENALTY=0.175531
PRIMER_RIGHT_3_PENALTY=0.178861
PRIMER_LEFT_3_SEQUENCE=GAGGCTCAAGGACAGCTCTC
PRIMER_RIGHT_3_SEQUENCE=AGGAGAGCTGGGAAAAGGGA
PRIMER_LEFT_3=614,20
PRIMER_RIGHT_3=753,20
PRIMER_LEFT_3_TM=59.824
PRIMER_RIGHT_3_TM=60.179
PRIMER_LEFT_3_GC_PERCENT=60.000
PRIMER_RIGHT_3_GC_PERCENT=55.000
PRIMER_LEFT_3_SELF_ANY_TH=10.25
PRIMER_RIGHT_3_SELF_ANY_TH=0.00
PRIMER_LEFT_3_SELF_END_TH=0.00
PRIMER_RIGHT_3_SELF_END_TH=0.00
PRIMER_LEFT_3_HAIRPIN_TH=37.05
PRIMER_RIGHT_3_HAIRPIN_TH=0.00
PRIMER_LEFT_3_END_STABILITY=3.2000
PRIMER_RIGHT_3_END_STABILITY=4.2000
PRIMER_PAIR_3_COMPL_ANY_TH=26.57
PRIMER_PAIR_3_COMPL_END_TH=26.57
PRIMER_PAIR_3_PRODUCT_SIZE=140
PRIMER_PAIR_3_PRODUCT_TM=83.2
PRIMER_PAIR_4_PENALTY=0.360353
PRIMER_LEFT_4_PENALTY=0.326264
PRIMER_RIGHT_4_PENALTY=0.034089
PRIMER_LEFT_4_SEQUENCE=TCCAGAAGCTGCTCTTTCCC
PRIMER_RIGHT_4_SEQUENCE=GCCTGGGTAGCTTTGGATGT
PRIMER_LEFT_4=808,20
PRIMER_RIGHT_4=943,20
PRIMER_LEFT_4_TM=59.674
PRIMER_RIGHT_4_TM=60.034
PRIMER_LEFT_4_GC_PERCENT=55.000
PRIMER_RIGHT_4_GC_PERCENT=55.000
PRIMER_LEFT_4_SELF_ANY_TH=1.00
PRIMER_RIGHT_4_SELF_ANY_TH=0.00
PRIMER_LEFT_4_SELF_END_TH=0.00
PRIMER_RIGHT_4_SELF_END_TH=0.00
PRIMER_LEFT_4_HAIRPIN_TH=34.56
PRIMER_RIGHT_4_HAIRPIN_TH=0.00
PRIMER_LEFT_4_END_STABILITY=3.9700
PRIMER_RIGHT_4_END_STABILITY=3.0600
PRIMER_PAIR_4_COMPL_ANY_TH=13.72
PRIMER_PAIR_4_COMPL_END_TH=10.53
PRIMER_PAIR_4_PRODUCT_SIZE=136
PRIMER_PAIR_4_PRODUCT_TM=83.6
=
)";

            std::stringstream file_handle(primer3_with_ref_pos_11_nopath);
            
            // Make a minimizer index
            gbwtgraph::DefaultMinimizerIndex minimizer_index(29, 11, false);
            // And zipcode collection
            ZipCodeCollection oversized_zipcodes;
           
            // Map node id to what gets stored in the payload - either the zipcode or index into oversized_zipcodes
            hash_map<vg::id_t, gbwtgraph::Payload> node_id_to_payload;
            node_id_to_payload.reserve(gbwt_graph.max_node_id() - gbwt_graph.min_node_id());
            
            // TODO: Use a function to deduplicate with minimizer_main
            gbwtgraph::index_haplotypes(gbwt_graph, minimizer_index, [&](const pos_t& pos) -> gbwtgraph::Payload {
                gbwtgraph::Payload payload = MIPayload::NO_CODE;

                #pragma omp critical 
                {
                    //If we've already seen this node before, then return the saved payload
                    if (node_id_to_payload.count(id(pos))) {
                        payload =  node_id_to_payload[id(pos)];
                    }
                }
                if (payload != MIPayload::NO_CODE) {
                    return payload;
                }
               

                ZipCode zipcode;
                zipcode.fill_in_zipcode(distance_index, pos);

                payload = zipcode.get_payload_from_zip();
                if (payload != MIPayload::NO_CODE) {
                    //If the zipcode is small enough to store in the payload
                    #pragma omp critical 
                    {
                        node_id_to_payload.emplace(id(pos), payload);
                    }
                    return payload;
                } else {
                    //Otherwise, if they are being saved, add the zipcode to the oversized zipcode list
                    //And remember the zipcode

                    //Fill in the decoder to be saved too
                    zipcode.fill_in_full_decoder();
                    
                    #pragma omp critical 
                    {
                        oversized_zipcodes.emplace_back(zipcode);
                        size_t zip_index = oversized_zipcodes.size() - 1;
                        payload= {0, zip_index};
                        node_id_to_payload.emplace(id(pos), payload);
                    }
                    return payload;
                }
            });

            MinimizerMapper giraffe_mapper(gbwt_graph, minimizer_index, &distance_index, &oversized_zipcodes);
            PrimerFinder primer_finder(graph, &distance_index, file_handle, gbwt_graph, gbwt_index, r_index, &giraffe_mapper);

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
