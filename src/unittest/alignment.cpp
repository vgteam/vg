/// \file alignment.cpp
///  
/// unit tests for Alignments and their utility functions
///

#include <iostream>
#include <string>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../alignment.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Default alignment simplification removes deletions on the edges of Mappings", "[alignment]") {

    string alignment_string = R"(
        {
            "sequence": "A",
            "path": {"mapping": [
                {
                    "position": {"node_id": 1},
                    "edit": [
                        {"from_length": 1, "to_length": 1},
                        {"from_length": 1}
                    ]
                },
                {
                    "position": {"node_id": 2},
                    "edit": [
                        {"from_length": 1, "to_length": 1}
                    ]
                }
            ]}
        }
    )";

    Alignment a;
    json2pb(a, alignment_string.c_str(), alignment_string.size());

    auto simple = simplify(a);
    REQUIRE(simple.path().mapping_size() == 2);
    REQUIRE(simple.path().mapping(0).edit_size() == 1);

}

TEST_CASE("Non-trim alignment simplification does not remove deletions on the edges of Mappings", "[alignment]") {

    string alignment_string = R"(
        {
            "sequence": "A",
            "path": {"mapping": [
                {
                    "position": {"node_id": 1},
                    "edit": [
                        {"from_length": 1, "to_length": 1},
                        {"from_length": 1}
                    ]
                },
                {
                    "position": {"node_id": 2},
                    "edit": [
                        {"from_length": 1, "to_length": 1}
                    ]
                }
            ]}
        }
    )";

    Alignment a;
    json2pb(a, alignment_string.c_str(), alignment_string.size());

    auto simple = simplify(a, false);
    REQUIRE(simple.path().mapping_size() == 2);
    REQUIRE(simple.path().mapping(0).edit_size() == 2);

}

TEST_CASE("Alignment simplification handles unaligned alignments", "[alignment]") {
    
    string alignment_string = R"(
        {"sequence": "TTGCTTTTTCCATCGCCAATTTCTGCAGCATTAATTTATTCTGGTTGTCATTTACAGTATATCAAGAGTTTAACCAGGAAATATTCCAAAGAAGGTTATTTTACCCATATACAGTACTTCATATGTCCAAAATTGGCCATTTGATGCACCTCTACCATCTATTTTTCCTTTCATGTGGCTTATATAACCATTGCTTGGATATTGATTTTACCAAAAGAGGGCGAATCCATATAATTATCTGGTATCGATGGAATGGTTATGGTAATTCAGAATATTTTTCTTATTATTGGTCGCTAATTGGCTCTTTTTATTTGGAGTTTCAAACAACAAGAAGAAGTCCGAGAATTCGTGCAATCTACTTTTGGAGTATTAAGCCATTTATTAATTTATCTTATAACATGTTGATGGACCATTCATATATGACCGACCCTAGCCCTAAACGAGTTGGTTTGCCTGCCGCCTACCAAATCATGTCATAATGGTAATAGTTGCATGCTTCTTTTTTCTGGTTAGTTTAAATTCGTTTTACTACCTAAATTCTTCATTAATAAAACTTTTACCCAAAATCACTTATTTGCATGATGGGCTCATTGTCCACAGAATGGGTGAAAAAGCCATGTCTGTGTTTTAGTTTTTCCTTTTCCTACATAAATATTTTAATAGACTACTATTGCATTTTTAATAAATAATATTCATGCCTGCATCTTTGAATAAAAGAAAAGAAACATTCGTTATTTCATCCTCTTACTTGATTGCTATAGAGCTATTTGCTTTCAAAGCTCTAGGTCCCACTGTCACGACACATACTATACTCTCACCGGTTCATTTTTCAACAGTTCATAAGGGCGATCATATTCTTTAACCTCACCGGATCCATCACAAATTATCTGTCGTAATCGATAACAGATCTCAACGATGTGCAATAGTAAGAATTGTGTTTTTATTAAACTCACTTCTTAATAATACCCTGAT", "path": {"mapping": [{"edit": [{"to_length": 793, "sequence": "TTGCTTTTTCCATCGCCAATTTCTGCAGCATTAATTTATTCTGGTTGTCATTTACAGTATATCAAGAGTTTAACCAGGAAATATTCCAAAGAAGGTTATTTTACCCATATACAGTACTTCATATGTCCAAAATTGGCCATTTGATGCACCTCTACCATCTATTTTTCCTTTCATGTGGCTTATATAACCATTGCTTGGATATTGATTTTACCAAAAGAGGGCGAATCCATATAATTATCTGGTATCGATGGAATGGTTATGGTAATTCAGAATATTTTTCTTATTATTGGTCGCTAATTGGCTCTTTTTATTTGGAGTTTCAAACAACAAGAAGAAGTCCGAGAATTCGTGCAATCTACTTTTGGAGTATTAAGCCATTTATTAATTTATCTTATAACATGTTGATGGACCATTCATATATGACCGACCCTAGCCCTAAACGAGTTGGTTTGCCTGCCGCCTACCAAATCATGTCATAATGGTAATAGTTGCATGCTTCTTTTTTCTGGTTAGTTTAAATTCGTTTTACTACCTAAATTCTTCATTAATAAAACTTTTACCCAAAATCACTTATTTGCATGATGGGCTCATTGTCCACAGAATGGGTGAAAAAGCCATGTCTGTGTTTTAGTTTTTCCTTTTCCTACATAAATATTTTAATAGACTACTATTGCATTTTTAATAAATAATATTCATGCCTGCATCTTTGAATAAAAGAAAAGAAACATTCGTTATTTCATCCTCTTACTTGATTGCTATAGAGCTATTTGCTTTCAAAGCTCTAGGTCCCACT"}, {"to_length": 18, "sequence": "GTCACGACACATACTATA"}, {"to_length": 161, "sequence": "CTCTCACCGGTTCATTTTTCAACAGTTCATAAGGGCGATCATATTCTTTAACCTCACCGGATCCATCACAAATTATCTGTCGTAATCGATAACAGATCTCAACGATGTGCAATAGTAAGAATTGTGTTTTTATTAAACTCACTTCTTAATAATACCCTGAT"}]}]}}
    )";
    
    Alignment a;
    json2pb(a, alignment_string.c_str(), alignment_string.size());
        
    auto simple = simplify(a);
    REQUIRE(simple.path().mapping_size() == 1);
    REQUIRE(simple.path().mapping(0).edit_size() == 1);
    
}

TEST_CASE("Alignment trimming works even on unaligned reads", "[alignment]") {
    
    string alignment_string = R"(
        {"sequence": "TTGCTTTTTCCATCGCCAATTTCTGCAGCATTAATTTATTCTGGTTGTCATTTACAGTATATCAAGAGTTTAACCAGGAAATATTCCAAAGAAGGTTATTTTACCCATATACAGTACTTCATATGTCCAAAATTGGCCATTTGATGCACCTCTACCATCTATTTTTCCTTTCATGTGGCTTATATAACCATTGCTTGGATATTGATTTTACCAAAAGAGGGCGAATCCATATAATTATCTGGTATCGATGGAATGGTTATGGTAATTCAGAATATTTTTCTTATTATTGGTCGCTAATTGGCTCTTTTTATTTGGAGTTTCAAACAACAAGAAGAAGTCCGAGAATTCGTGCAATCTACTTTTGGAGTATTAAGCCATTTATTAATTTATCTTATAACATGTTGATGGACCATTCATATATGACCGACCCTAGCCCTAAACGAGTTGGTTTGCCTGCCGCCTACCAAATCATGTCATAATGGTAATAGTTGCATGCTTCTTTTTTCTGGTTAGTTTAAATTCGTTTTACTACCTAAATTCTTCATTAATAAAACTTTTACCCAAAATCACTTATTTGCATGATGGGCTCATTGTCCACAGAATGGGTGAAAAAGCCATGTCTGTGTTTTAGTTTTTCCTTTTCCTACATAAATATTTTAATAGACTACTATTGCATTTTTAATAAATAATATTCATGCCTGCATCTTTGAATAAAAGAAAAGAAACATTCGTTATTTCATCCTCTTACTTGATTGCTATAGAGCTATTTGCTTTCAAAGCTCTAGGTCCCACTGTCACGACACATACTATACTCTCACCGGTTCATTTTTCAACAGTTCATAAGGGCGATCATATTCTTTAACCTCACCGGATCCATCACAAATTATCTGTCGTAATCGATAACAGATCTCAACGATGTGCAATAGTAAGAATTGTGTTTTTATTAAACTCACTTCTTAATAATACCCTGAT", "path": {"mapping": [{"edit": [{"to_length": 793, "sequence": "TTGCTTTTTCCATCGCCAATTTCTGCAGCATTAATTTATTCTGGTTGTCATTTACAGTATATCAAGAGTTTAACCAGGAAATATTCCAAAGAAGGTTATTTTACCCATATACAGTACTTCATATGTCCAAAATTGGCCATTTGATGCACCTCTACCATCTATTTTTCCTTTCATGTGGCTTATATAACCATTGCTTGGATATTGATTTTACCAAAAGAGGGCGAATCCATATAATTATCTGGTATCGATGGAATGGTTATGGTAATTCAGAATATTTTTCTTATTATTGGTCGCTAATTGGCTCTTTTTATTTGGAGTTTCAAACAACAAGAAGAAGTCCGAGAATTCGTGCAATCTACTTTTGGAGTATTAAGCCATTTATTAATTTATCTTATAACATGTTGATGGACCATTCATATATGACCGACCCTAGCCCTAAACGAGTTGGTTTGCCTGCCGCCTACCAAATCATGTCATAATGGTAATAGTTGCATGCTTCTTTTTTCTGGTTAGTTTAAATTCGTTTTACTACCTAAATTCTTCATTAATAAAACTTTTACCCAAAATCACTTATTTGCATGATGGGCTCATTGTCCACAGAATGGGTGAAAAAGCCATGTCTGTGTTTTAGTTTTTCCTTTTCCTACATAAATATTTTAATAGACTACTATTGCATTTTTAATAAATAATATTCATGCCTGCATCTTTGAATAAAAGAAAAGAAACATTCGTTATTTCATCCTCTTACTTGATTGCTATAGAGCTATTTGCTTTCAAAGCTCTAGGTCCCACT"}, {"to_length": 18, "sequence": "GTCACGACACATACTATA"}, {"to_length": 161, "sequence": "CTCTCACCGGTTCATTTTTCAACAGTTCATAAGGGCGATCATATTCTTTAACCTCACCGGATCCATCACAAATTATCTGTCGTAATCGATAACAGATCTCAACGATGTGCAATAGTAAGAATTGTGTTTTTATTAAACTCACTTCTTAATAATACCCTGAT"}]}]}}
    )";
    
    Alignment a;
    json2pb(a, alignment_string.c_str(), alignment_string.size());
        
    auto aln = simplify(a);
    aln = strip_from_start(aln, 50);
    aln = strip_from_end(aln, 50);

    REQUIRE(a.sequence().size() - 100 == aln.sequence().size());
    
}
    
TEST_CASE("Alignment normalization behaves as expected","[alignment]") {
    string alignment_string = R"(
    {"sequence":"ATNNNNANCT","path":{"mapping":[{"position":{"node_id":"1"},"edit":[{"from_length":"1","to_length":"1"},{"from_length":"1","to_length":"1"},{"from_length":"1","to_length":"1","sequence":"N"},{"from_length":"1","to_length":"1"},{"from_length":"2","to_length":"2","sequence":"NN"},{"from_length":"2","to_length":"2","sequence":"AN"},{"to_length":"1","sequence":"C"},{"to_length":"1","sequence":"T"},{"from_length":"1"},{"from_length":"1"}]}]}}
    )";
    
    string normalized_string = R"(
    {"sequence":"ATNNNNANCT","path":{"mapping":[{"position":{"node_id":"1"},"edit":[{"from_length":"2","to_length":"2"},{"from_length":"4","to_length":"4","sequence":"NNNN"},{"from_length":"1","to_length":"1","sequence":"A"},{"from_length":"1","to_length":"1","sequence":"N"},{"to_length":"2","sequence":"CT"},{"from_length":"2"}]}]}}
    )";
    
    Alignment aln;
    json2pb(aln, alignment_string.c_str(), alignment_string.size());
    Alignment target;
    json2pb(target, normalized_string.c_str(), normalized_string.size());
    
    normalize_alignment(aln);
    
    REQUIRE(aln.path().mapping_size() == target.path().mapping_size());
    for (size_t i = 0; i < aln.path().mapping_size(); i++) {
        REQUIRE(aln.path().mapping(i).edit_size() == target.path().mapping(i).edit_size());
        for (size_t j = 0; j < target.path().mapping(i).edit_size(); j++) {
            REQUIRE(aln.path().mapping(i).edit(j).from_length() == target.path().mapping(i).edit(j).from_length());
            REQUIRE(aln.path().mapping(i).edit(j).to_length() == target.path().mapping(i).edit(j).to_length());
            REQUIRE(aln.path().mapping(i).edit(j).sequence() == target.path().mapping(i).edit(j).sequence());
        }
    }
}
    
    
TEST_CASE("Target to alignment extraction", "[target-to-aln]") {
    
    VG vg;
    
    Node* n0 = vg.create_node("CGA");
    Node* n1 = vg.create_node("TTGG");
    Node* n2 = vg.create_node("CCGT");
    Node* n3 = vg.create_node("C");
    Node* n4 = vg.create_node("GT");
    Node* n5 = vg.create_node("GATAA");
    Node* n6 = vg.create_node("CGG");
    Node* n7 = vg.create_node("ACA");
    Node* n8 = vg.create_node("GCCG");
    Node* n9 = vg.create_node("A");
    Node* n10 = vg.create_node("C");
    Node* n11 = vg.create_node("G");
    Node* n12 = vg.create_node("T");
    Node* n13 = vg.create_node("A");
    Node* n14 = vg.create_node("C");
    Node* n15 = vg.create_node("C");
    
    vg.create_edge(n0, n1);
    vg.create_edge(n2, n0, true, true);
    vg.create_edge(n1, n3);
    vg.create_edge(n2, n3);
    vg.create_edge(n3, n4, false, true);
    vg.create_edge(n4, n5, true, false);
    vg.create_edge(n5, n6);
    vg.create_edge(n8, n6, false, true);
    vg.create_edge(n6, n7, false, true);
    vg.create_edge(n7, n9, true, true);
    vg.create_edge(n9, n10, true, false);
    vg.create_edge(n10, n11, false, false);
    vg.create_edge(n12, n11, false, true);
    vg.create_edge(n13, n12, false, false);
    vg.create_edge(n14, n13, true, false);
    vg.create_edge(n15, n14, true, true);
    
    Graph graph = vg.graph;
    
    Path* path = graph.add_path();
    path->set_name("path");
    Mapping* mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(n0->id());
    mapping->set_rank(1);
    mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(n2->id());
    mapping->set_rank(2);
    mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(n3->id());
    mapping->set_rank(3);
    mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(n4->id());
    mapping->mutable_position()->set_is_reverse(true);
    mapping->set_rank(4);
    mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(n5->id());
    mapping->set_rank(5);
    mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(n6->id());
    mapping->set_rank(6);
    mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(n8->id());
    mapping->mutable_position()->set_is_reverse(true);
    mapping->set_rank(7);
    
    xg::XG xg_index;
    xg_index.from_path_handle_graph(VG(graph));
    
    SECTION("Subpath getting gives us the expected 1bp alignment") {
        Alignment target = target_alignment(&xg_index, "path", 1, 2, "feature", false);
        REQUIRE(alignment_from_length(target) == 2 - 1);
    }
    
    SECTION("Subpath getting gives us the expected 10bp alignment") {
        Alignment target = target_alignment(&xg_index, "path", 10, 20, "feature", false);
        REQUIRE(alignment_from_length(target) == 20 - 10);
    }
    
    SECTION("Subpath getting gives us the expected 14bp alignment") {
        Alignment target = target_alignment(&xg_index, "path", 0, 14, "feature", false);
        REQUIRE(alignment_from_length(target) == 14);
    }
    
    SECTION("Subpath getting gives us the expected 21bp alignment") {
        Alignment target = target_alignment(&xg_index, "path", 0, 21, "feature", false);
        REQUIRE(alignment_from_length(target) == 21);
    }
    
    SECTION("Subpath getting gives us the expected inverted 7bp alignment") {
        Alignment target = target_alignment(&xg_index, "path", 0, 7, "feature", true);
        REQUIRE(alignment_from_length(target) == 7);
        REQUIRE(target.path().mapping(0).position().node_id() == n2->id());
        REQUIRE(target.path().mapping(1).position().node_id() == n0->id());
        REQUIRE(target.path().mapping(0).position().is_reverse() == true);
        REQUIRE(target.path().mapping(1).position().is_reverse() == true);
    }
    
}

}
}
