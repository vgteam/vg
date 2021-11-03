/// \file alignment.cpp
///  
/// unit tests for Alignments and their utility functions
///

#include <iostream>
#include <string>
#include "vg/io/json2pb.h"
#include <bdsg/hash_graph.hpp>
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

TEST_CASE("Alignments can be left-shifted", "[left-shift]") {
    
    HashGraph g;
    
    //                                                              v We have a substitution here
    handle_t n1 = g.create_handle("AAAACATTAGCATTAGCATTAGCATTAGCATTATCATTAGAAAA");
    // This will be node 1
    
    SECTION("when a deletion is in a homopolymer") {
        Alignment aln;
        json2pb(aln, R"(
            {"path": {"mapping": [
                {"position": {"node_id": "1"}, "edit": {[
                    {"from": 3, "to": 3},
                    {"from": 1, "to": 0},
                    {"from": 1, "to": 1}
                ]}}
            ]}}
        )");
        
        left_shift_alignment_in_place(&aln, g);
        
        REQUIRE(aln.path().mapping(0).edit_size() == 2);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 1);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(0).sequence() == ""); 
        REQUIRE(aln.path().mapping(0).edit(1).from_length() == 4);
        REQUIRE(aln.path().mapping(0).edit(1).to_length() == 4);
        REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
    }
    
    SECTION("when deletions in a homopolymer must merge") {
        Alignment aln;
        json2pb(aln, R"(
            {"path": {"mapping": [
                {"position": {"node_id": "1"}, "edit": {[
                    {"from": 1, "to": 1},
                    {"from": 1, "to": 0},
                    {"from": 1, "to": 1},
                    {"from": 1, "to": 0},
                    {"from": 1, "to": 1}
                ]}}
            ]}}
        )");
        
        left_shift_alignment_in_place(&aln, g);
        
        REQUIRE(aln.path().mapping(0).edit_size() == 2);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 2);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
        REQUIRE(aln.path().mapping(0).edit(1).from_length() == 3);
        REQUIRE(aln.path().mapping(0).edit(1).to_length() == 3);
        REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
    }
    
    SECTION("when an insertion is in a homopolymer") {
        Alignment aln;
        json2pb(aln, R"(
            {"path": {"mapping": [
                {"position": {"node_id": "1"}, "edit": {[
                    {"from": 4, "to": 4},
                    {"from": 0, "to": 1, "sequence": "A"},
                    {"from": 1, "to": 1}
                ]}}
            ]}}
        )");
        
        left_shift_alignment_in_place(&aln, g);
        
        REQUIRE(aln.path().mapping(0).edit_size() == 2);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 1);
        REQUIRE(aln.path().mapping(0).edit(0).sequence() == "A");
        REQUIRE(aln.path().mapping(0).edit(1).from_length() == 5);
        REQUIRE(aln.path().mapping(0).edit(1).to_length() == 5);
        REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
    }
    
    SECTION("when insertions in a homopolymer must merge") {
        Alignment aln;
        json2pb(aln, R"(
            {"path": {"mapping": [
                {"position": {"node_id": "1"}, "edit": {[
                    {"from": 2, "to": 2},
                    {"from": 0, "to": 1, "sequence": "A"},
                    {"from": 2, "to": 2},
                    {"from": 0, "to": 1, "sequence": "A"},
                    {"from": 1, "to": 1}
                ]}}
            ]}}
        )");
        
        left_shift_alignment_in_place(&aln, g);
        
        REQUIRE(aln.path().mapping(0).edit_size() == 2);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 2);
        REQUIRE(aln.path().mapping(0).edit(0).sequence() == "AA");
        REQUIRE(aln.path().mapping(0).edit(1).from_length() == 5);
        REQUIRE(aln.path().mapping(0).edit(1).to_length() == 5);
        REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
    }
    
    SECTION("when an insertion in a homopolymer must merge with an immovable insertion") {
        Alignment aln;
        json2pb(aln, R"(
            {"path": {"mapping": [
                {"position": {"node_id": "1"}, "edit": {[
                    {"from": 2, "to": 2},
                    {"from": 0, "to": 1, "sequence": "G"},
                    {"from": 2, "to": 2},
                    {"from": 0, "to": 1, "sequence": "A"},
                    {"from": 1, "to": 1}
                ]}}
            ]}}
        )");
        
        left_shift_alignment_in_place(&aln, g);
        
        REQUIRE(aln.path().mapping(0).edit_size() == 3);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 2);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 2);
        REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
        REQUIRE(aln.path().mapping(0).edit(1).from_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(1).to_length() == 2);
        REQUIRE(aln.path().mapping(0).edit(1).sequence() == "GA");
        REQUIRE(aln.path().mapping(0).edit(2).from_length() == 3);
        REQUIRE(aln.path().mapping(0).edit(2).to_length() == 3);
        REQUIRE(aln.path().mapping(0).edit(2).sequence() == "");
    }

    SECTION("when the last matching repeat is deleted") {
        Alignment aln;
        json2pb(aln, R"(
            {"path": {"mapping": [
                {"position": {"node_id": "1"}, "edit": {[
                    {"from": 22, "to": 22},
                    {"from": 6, "to": 0},
                    {"from": 16, "to": 16}
                ]}}
            ]}}
        )");
        
        left_shift_alignment_in_place(&aln, g);
        
        REQUIRE(aln.path().mapping(0).edit_size() == 3);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 4);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 4);
        REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
        REQUIRE(aln.path().mapping(0).edit(1).from_length() == 6);
        REQUIRE(aln.path().mapping(0).edit(1).to_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
        REQUIRE(aln.path().mapping(0).edit(2).from_length() == 34);
        REQUIRE(aln.path().mapping(0).edit(2).to_length() == 34);
        REQUIRE(aln.path().mapping(0).edit(2).sequence() == "");
    }
    
    SECTION("when the next matching repeat is deleted and the previous repeat made to match") {
        // TODO: should this one even be caught??? Does Freebayes catch it?
        Alignment aln;
        json2pb(aln, R"(
            {"path": {"mapping": [
                {"position": {"node_id": "1"}, "edit": {[
                    {"from": 28, "to": 28},
                    {"from": 5, "to": 5},
                    {"from": 1, "to": 1, "sequence": "G"},
                    {"from": 6, "to": 0},
                    {"from": 4, "to": 4}
                ]}}
            ]}}
        )");
        
        left_shift_alignment_in_place(&aln, g);
        
        // Comes out as a deletion of the repeat that was actually distinctive
        REQUIRE(aln.path().mapping(0).edit_size() == 3);
        REQUIRE(aln.path().mapping(0).edit(0).from_length() == 28);
        REQUIRE(aln.path().mapping(0).edit(0).to_length() == 28);
        REQUIRE(aln.path().mapping(0).edit(0).sequence() == "");
        REQUIRE(aln.path().mapping(0).edit(1).from_length() == 6);
        REQUIRE(aln.path().mapping(0).edit(1).to_length() == 0);
        REQUIRE(aln.path().mapping(0).edit(1).sequence() == "");
        REQUIRE(aln.path().mapping(0).edit(2).from_length() == 10);
        REQUIRE(aln.path().mapping(0).edit(2).to_length() == 10);
        REQUIRE(aln.path().mapping(0).edit(2).sequence() == "");
    }
    
}
    

TEST_CASE("consolidate_ID_runs merges runs of adjacent I's and D's in cigars", "[alignment][surject]") {
    
    vector<pair<int, char>> cigar{
        make_pair(2, 'D'),
        make_pair(1, 'I'),
        make_pair(4, 'D'),
        make_pair(1, 'M'),
        make_pair(3, 'I'),
        make_pair(5, 'D'),
        make_pair(1, 'I')
    };

    consolidate_ID_runs(cigar);
    REQUIRE(cigar.size() == 5);
    bool consolidated_1 = ((cigar[0] == make_pair(6, 'D') && cigar[1] == make_pair(1, 'I'))
                           || (cigar[0] == make_pair(1, 'I') && cigar[1] == make_pair(6, 'D')));
    bool consolidated_2 = ((cigar[3] == make_pair(5, 'D') && cigar[4] == make_pair(4, 'I'))
                           || (cigar[3] == make_pair(4, 'I') && cigar[4] == make_pair(5, 'D')));
    
    REQUIRE(consolidated_1);
    REQUIRE(consolidated_2);
}

TEST_CASE("Inter-alignment distance computation for HTS output formats matches BWA", "[alignment]") {
    // See https://github.com/vgteam/vg/issues/3078. We want to match BWA on
    // these straightforward, fully-matching reads.
    auto lengths = compute_template_lengths(10206220, {{151, 'M'}}, 10206662, {{151, 'M'}});
    REQUIRE(lengths.first == 593);
    REQUIRE(lengths.second == -593);
}

TEST_CASE("CIGAR generation forces adjacent insertions and deletions to obey GATK's constraints", "[alignment]") {
    // See https://github.com/vgteam/vg/issues/3080
    vector<pair<int, char>> cigar;
    
    SECTION("DID becomes DI") {
        append_cigar_operation(1, 'D', cigar);
        append_cigar_operation(5, 'I', cigar);
        append_cigar_operation(2, 'D', cigar);
        
        REQUIRE(cigar.size() == 2);
        REQUIRE(cigar[0].first == 3);
        REQUIRE(cigar[0].second == 'D');
        REQUIRE(cigar[1].first == 5);
        REQUIRE(cigar[1].second == 'I');
    }
    
    SECTION("MMIDIIDDIMM becomes MDIM") {
        append_cigar_operation(1, 'M', cigar);
        append_cigar_operation(1, 'M', cigar);
        append_cigar_operation(1, 'I', cigar);
        append_cigar_operation(1, 'D', cigar);
        append_cigar_operation(1, 'I', cigar);
        append_cigar_operation(1, 'I', cigar);
        append_cigar_operation(1, 'D', cigar);
        append_cigar_operation(1, 'D', cigar);
        append_cigar_operation(1, 'I', cigar);
        append_cigar_operation(1, 'M', cigar);
        append_cigar_operation(1, 'M', cigar);
        
        REQUIRE(cigar.size() == 4);
        REQUIRE(cigar[0].first == 2);
        REQUIRE(cigar[0].second == 'M');
        REQUIRE(cigar[1].first == 3);
        REQUIRE(cigar[1].second == 'D');
        REQUIRE(cigar[2].first == 4);
        REQUIRE(cigar[2].second == 'I');
        REQUIRE(cigar[3].first == 2);
        REQUIRE(cigar[3].second == 'M');
    }
    
}

}
}
