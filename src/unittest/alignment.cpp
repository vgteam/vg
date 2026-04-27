/// \file alignment.cpp
///  
/// unit tests for Alignments and their utility functions
///

#include <iostream>
#include <string>
#include <google/protobuf/util/message_differencer.h>

#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../alignment.hpp"
#include "../vg.hpp"
#include "xg.hpp"
#include "../handle.hpp"
#include "../annotation.hpp"
#include "catch.hpp"

#include "bdsg/hash_graph.hpp"

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
    path_handle_t path_handle = xg_index.get_path_handle("path");
    
    SECTION("Subpath getting gives us the expected 1bp alignment") {
        Alignment target = target_alignment(&xg_index, path_handle, 1, 2, "feature", false);
        REQUIRE(alignment_from_length(target) == 2 - 1);
    }
    
    SECTION("Subpath getting gives us the expected 10bp alignment") {
        Alignment target = target_alignment(&xg_index, path_handle, 10, 20, "feature", false);
        REQUIRE(alignment_from_length(target) == 20 - 10);
    }
    
    SECTION("Subpath getting gives us the expected 14bp alignment") {
        Alignment target = target_alignment(&xg_index, path_handle, 0, 14, "feature", false);
        REQUIRE(alignment_from_length(target) == 14);
    }
    
    SECTION("Subpath getting gives us the expected 21bp alignment") {
        Alignment target = target_alignment(&xg_index, path_handle, 0, 21, "feature", false);
        REQUIRE(alignment_from_length(target) == 21);
    }
    
    SECTION("Subpath getting gives us the expected inverted 7bp alignment") {
        Alignment target = target_alignment(&xg_index, path_handle, 0, 7, "feature", true);
        REQUIRE(alignment_from_length(target) == 7);
        REQUIRE(target.path().mapping(0).position().node_id() == n2->id());
        REQUIRE(target.path().mapping(1).position().node_id() == n0->id());
        REQUIRE(target.path().mapping(0).position().is_reverse() == true);
        REQUIRE(target.path().mapping(1).position().is_reverse() == true);
    }
    
}

TEST_CASE("simplify_cigar merges runs of adjacent I's and D's in cigars", "[alignment][surject]") {
    
    vector<pair<int, char>> cigar{
        make_pair(2, 'D'),
        make_pair(1, 'I'),
        make_pair(4, 'D'),
        make_pair(1, 'M'),
        make_pair(3, 'I'),
        make_pair(5, 'D'),
        make_pair(1, 'I')
    };

    simplify_cigar(cigar);
    REQUIRE(cigar.size() == 5);
    bool consolidated_1 = ((cigar[0] == make_pair(6, 'D') && cigar[1] == make_pair(1, 'I'))
                           || (cigar[0] == make_pair(1, 'I') && cigar[1] == make_pair(6, 'D')));
    bool consolidated_2 = ((cigar[3] == make_pair(5, 'D') && cigar[4] == make_pair(4, 'I'))
                           || (cigar[3] == make_pair(4, 'I') && cigar[4] == make_pair(5, 'D')));
    
    REQUIRE(consolidated_1);
    REQUIRE(consolidated_2);
}

TEST_CASE("simplify_cigar merges runs of adjacent operations and removes empty operations", "[alignment][surject]") {
    
    vector<pair<int, char>> cigar{
        make_pair(2, 'S'),
        make_pair(1, 'M'),
        make_pair(1, 'M'),
        make_pair(0, 'D'),
        make_pair(1, 'I'),
        make_pair(1, 'M')
    };
    
    simplify_cigar(cigar);
    REQUIRE(cigar.size() == 4);
    REQUIRE(cigar[0] == make_pair(2, 'S'));
    REQUIRE(cigar[1] == make_pair(2, 'M'));
    REQUIRE(cigar[2] == make_pair(1, 'I'));
    REQUIRE(cigar[3] == make_pair(1, 'M'));
}

TEST_CASE("Template length for HTS output formats matches BWA", "[alignment][tlen]") {
    // See https://github.com/vgteam/vg/issues/3078. We want to match BWA on
    // these straightforward, fully-matching reads.
    auto lengths = compute_template_lengths(10206220, {{151, 'M'}}, 10206662, {{151, 'M'}});
    REQUIRE(lengths.first == 593);
    REQUIRE(lengths.second == -593);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 1 starts inside read 2", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(1000, {{151, 'M'}}, 900, {{151, 'M'}});
    REQUIRE(lengths.first == -251);
    REQUIRE(lengths.second == 251);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 2 starts inside read 1", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(900, {{151, 'M'}}, 1000, {{151, 'M'}});
    REQUIRE(lengths.first == 251);
    REQUIRE(lengths.second == -251);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 1 is contained in read 2", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(1000, {{10, 'M'}}, 900, {{151, 'M'}});
    REQUIRE(lengths.first == -151); // Read 1 is the rightmost since it starts latest
    REQUIRE(lengths.second == 151);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 2 is contained in read 1", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(900, {{151, 'M'}}, 1000, {{10, 'M'}});
    REQUIRE(lengths.first == 151); // Read 1 is the leftmost since it starts earliest
    REQUIRE(lengths.second == -151);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 1 is a prefix of read 2", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(1000, {{10, 'M'}}, 1000, {{151, 'M'}});
    REQUIRE(lengths.first == 151); // Read 1 is the leftmost since it ends earliest
    REQUIRE(lengths.second == -151);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 2 is a prefix of read 1", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(1000, {{151, 'M'}}, 1000, {{10, 'M'}});
    REQUIRE(lengths.first == -151); // Read 1 is the rightmost since it ends latest
    REQUIRE(lengths.second == 151);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 1 is a suffix of read 2", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(1141, {{10, 'M'}}, 1000, {{151, 'M'}});
    REQUIRE(lengths.first == -151); // Read 1 is the rightmost since it starts latest
    REQUIRE(lengths.second == 151);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 2 is a suffix of read 1", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(1000, {{151, 'M'}}, 1141, {{10, 'M'}});
    REQUIRE(lengths.first == 151); // Read 1 is the leftmost since it starts earliest
    REQUIRE(lengths.second == -151);
}

TEST_CASE("Template length for HTS output formats is outermost to outermost when read 1 and read 2 have identical coordinates", "[alignment][tlen]") {
    auto lengths = compute_template_lengths(1000, {{151, 'M'}}, 1000, {{151, 'M'}});
    // Signs are arbitrary but must be opposite
    REQUIRE(std::abs(lengths.first) == 151);
    REQUIRE(std::abs(lengths.second) == 151);
    REQUIRE(lengths.first + lengths.second == 0);
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

TEST_CASE("Conversion to GAF removes an unused final node", "[alignment][gaf]") {
    VG graph;
    graph.create_node("GATTACA", 1);
    graph.create_node("CAT", 2);
    graph.create_node("GATTA", 3);
    graph.create_edge(graph.get_handle(1, false), graph.get_handle(2, false));
    graph.create_edge(graph.get_handle(2, false), graph.get_handle(3, false));

    Alignment aln;
    {
        aln.set_sequence("TACACTTAC");
        aln.set_name("softclip-at-end");
        Path* path = aln.mutable_path();

        // TACA to 1:3
        Mapping* mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(1);
        mapping->mutable_position()->set_is_reverse(false);
        mapping->mutable_position()->set_offset(3);
        mapping->set_rank(1);
        Edit* edit = mapping->add_edit();
        edit->set_from_length(4);
        edit->set_to_length(4);

        // CTT to 2:0
        mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(2);
        mapping->mutable_position()->set_is_reverse(false);
        mapping->mutable_position()->set_offset(0);
        mapping->set_rank(2);
        edit = mapping->add_edit();
        edit->set_from_length(1);
        edit->set_to_length(1);
        edit = mapping->add_edit();
        edit->set_from_length(1);
        edit->set_to_length(1);
        edit->set_sequence("T");
        edit = mapping->add_edit();
        edit->set_from_length(1);
        edit->set_to_length(1);

        // AC to 3:0 as a softclip
        mapping = path->add_mapping();
        mapping->mutable_position()->set_node_id(3);
        mapping->mutable_position()->set_is_reverse(false);
        mapping->mutable_position()->set_offset(0);
        mapping->set_rank(3);
        edit = mapping->add_edit();
        edit->set_from_length(0);
        edit->set_to_length(2);
        edit->set_sequence("AC");
    }

    // The unused final node should be removed and path_end should be set correctly.
    gafkluge::GafRecord gaf = alignment_to_gaf(graph, aln, nullptr, true, false, false);
    SECTION("GAF record is correct") {
        REQUIRE(gaf.query_name == aln.name());
        REQUIRE(gaf.query_length == aln.sequence().size());
        REQUIRE(gaf.query_start == 0);
        REQUIRE(gaf.query_end == aln.sequence().size());
        REQUIRE(gaf.path.size() == 2);
        REQUIRE(gaf.path[0].name == "1");
        REQUIRE(gaf.path[0].is_reverse == false);
        REQUIRE(gaf.path[1].name == "2");
        REQUIRE(gaf.path[1].is_reverse == false);
        size_t path_length = graph.get_length(graph.get_handle(1)) + graph.get_length(graph.get_handle(2));
        REQUIRE(gaf.path_length == path_length);
        REQUIRE(gaf.path_start == aln.path().mapping(0).position().offset());
        REQUIRE(gaf.path_end == path_length);
        std::string difference_string = gaf.opt_fields["cs"].second;
        std::string true_difference_string = ":5*AT:1+AC";
        REQUIRE(difference_string == true_difference_string);
    }
}

TEST_CASE("Supplementary alignments can identified and processed", "[alignment][supplementary]") {

    string seq = "CAAGTTCTGCTTTCTGCGAT";
    string qual(seq.size(), 40);
    string name = "read";
    string group = "group";
    string samp = "samp";

    Alignment primary, supp1, supp2, irrelevant;
    for (auto aln :  {&primary, &supp1, &supp2, &irrelevant}) {
        aln->set_sequence(seq);
        aln->set_quality(qual);
        aln->set_name(name);
        aln->set_read_group(group);
        aln->set_sample_name(samp);
    }

    {
        // CAAGTTCTGCTTTCTGCGAT
        // ......CTGCTTTC......
        auto path = primary.mutable_path();
        auto m = path->add_mapping();
        m->set_rank(1);
        auto p = m->mutable_position();
        p->set_node_id(1);
        p->set_offset(4);
        p->set_is_reverse(false);
        auto e1 = m->add_edit();
        e1->set_from_length(0);
        e1->set_to_length(6);
        e1->set_sequence("CAAGTT");
        auto e2 = m->add_edit();
        e2->set_from_length(8);
        e2->set_to_length(8);
        auto e3 = m->add_edit();
        e3->set_from_length(0);
        e3->set_to_length(6);
        e3->set_sequence("TGCGAT");
        primary.set_score(8);
        primary.set_mapping_quality(40);
    }

    {
        // C-AAGTTCTGCTTTCTGCGAT
        // CAAAGGT.............
        auto path = supp1.mutable_path();
        auto m = path->add_mapping();
        m->set_rank(1);
        auto p = m->mutable_position();
        p->set_node_id(4);
        p->set_offset(2);
        p->set_is_reverse(true);
        auto e1 = m->add_edit();
        e1->set_from_length(1);
        e1->set_to_length(1);
        auto e2 = m->add_edit();
        e2->set_from_length(1);
        e2->set_to_length(0);
        auto e3 = m->add_edit();
        e3->set_from_length(3);
        e3->set_to_length(3);
        auto e4 = m->add_edit();
        e4->set_from_length(1);
        e4->set_to_length(1);
        e4->set_sequence("T");
        auto e5 = m->add_edit();
        e5->set_from_length(1);
        e5->set_to_length(1);
        auto e6 = m->add_edit();
        e6->set_from_length(0);
        e6->set_to_length(14);
        e6->set_sequence("CTGCTTTCTGCGAT");
        supp1.set_score(5);
        supp1.set_mapping_quality(15);
    }

    {
        // CAAGTTCTGCTTTCTGCGAT
        // ..............TGC-AT
        auto path = supp2.mutable_path();
        auto m1 = path->add_mapping();
        m1->set_rank(1);
        auto p1 = m1->mutable_position();
        p1->set_node_id(10);
        p1->set_offset(4);
        p1->set_is_reverse(false);
        auto e0 = m1->add_edit();
        e0->set_from_length(0);
        e0->set_to_length(14);
        e0->set_sequence("CAAGTTCTGCTTTC");
        auto e1 = m1->add_edit();
        e1->set_from_length(3);
        e1->set_to_length(3);
        auto m2 = path->add_mapping();
        m2->set_rank(2);
        auto p2 = m2->mutable_position();
        p2->set_node_id(11);
        p2->set_offset(0);
        p2->set_is_reverse(true);
        auto e2 = m2->add_edit();
        e2->set_from_length(0);
        e2->set_to_length(1);
        e2->set_sequence("G");
        auto e3 = m2->add_edit();
        e3->set_from_length(2);
        e3->set_to_length(2);
        supp2.set_score(5);
        supp2.set_mapping_quality(10);
    }

    {
        // CAAGTTCTGCTTTCTGCGAT
        // .......TGCTTT.......
        auto path = irrelevant.mutable_path();
        auto m = path->add_mapping();
        m->set_rank(1);
        auto p = m->mutable_position();
        p->set_node_id(20);
        p->set_offset(0);
        p->set_is_reverse(false);
        auto e1 = m->add_edit();
        e1->set_from_length(0);
        e1->set_to_length(7);
        e1->set_sequence("CAAGTTC");
        auto e2 = m->add_edit();
        e2->set_from_length(6);
        e2->set_to_length(6);
        auto e3 = m->add_edit();
        e3->set_from_length(0);
        e3->set_to_length(7);
        e3->set_sequence("CTGCGAT");
        irrelevant.set_score(5);
        irrelevant.set_mapping_quality(1);
    }

    bdsg::HashGraph graph;
    auto h1 = graph.create_handle("ATATCTGCTTTC", 1);
    auto h4 = graph.create_handle(reverse_complement("GGCAAAGGT"), 4);
    auto h10 = graph.create_handle("CCGCTGC", 10);
    auto h11 = graph.create_handle(reverse_complement("ATTC"), 11);
    auto h20 = graph.create_handle("TGCTTT", 20);
    graph.create_edge(h1, graph.flip(h11));

    SECTION("Supplementaries can be correctly identified among a list of Alignments") {

        vector<Alignment> alns{primary, supp1, supp2, irrelevant};

        size_t min_size = 5;
        size_t separation = 1;
        double read_coverage = 0.95;
        double score_fraction = 0.5;

        SECTION("Can identify a complete partition of the read") {

            auto supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            sort(supps.begin(), supps.end());
            REQUIRE(supps == vector<size_t>{1, 2});
        }

        SECTION("Respects the mininmum size constraint") {
            min_size = 6;
            auto supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 2);
            min_size = 7;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.empty());
        }

        SECTION("Respects the mininmum score constraint") {
            score_fraction = 5.0 / 8.0;
            auto supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 2);
            score_fraction = 0.7;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.empty());
        }

        SECTION("Respects read coverage constraint") {

            read_coverage = 1.0;
            auto supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 2);

            auto m = alns[2].mutable_path()->mutable_mapping(0);

            // move an aligned base into softclip
            auto e1 = m->mutable_edit(0);
            auto e2 = m->mutable_edit(1);
            e1->set_sequence("CAAGTTCTGCTTTCT");
            e1->set_to_length(15);
            e2->set_from_length(2);
            e2->set_to_length(2);

            read_coverage = 0.95;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 2);

            read_coverage = 0.96;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 0);
        }

        SECTION("Respects separation constraint") {

            separation = 0;
            auto supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 2);

            auto m = alns[2].mutable_path()->mutable_mapping(0);

            // move an aligned base into softclip
            auto e1 = m->mutable_edit(0);
            auto e2 = m->mutable_edit(1);
            e1->set_sequence("CAAGTTCTGCTTTCT");
            e1->set_to_length(15);
            e2->set_from_length(2);
            e2->set_to_length(2);

            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 0);

            separation = 1;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 2);

            // add a base of overlap
            e1->set_sequence("CAAGTTCTGCTTT");
            e1->set_to_length(13);
            e2->set_from_length(4);
            e2->set_to_length(4);
            m->mutable_position()->set_offset(3);

            separation = 0;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 0);

            separation = 1;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 2);
            
            // return it to the original state
            e1->set_sequence("CAAGTTCTGCTTTC");
            e1->set_to_length(14);
            e2->set_from_length(3);
            e2->set_to_length(3);

            // split it into two alignments

            alns.emplace_back(alns[2]);
            auto& split1 = alns[2];
            auto& split2 = alns.back();
            split1.mutable_path()->mutable_mapping()->RemoveLast();
            auto m2 = split2.mutable_path()->mutable_mapping(0);
            *m2 = split2.path().mapping(1);
            split2.mutable_path()->mutable_mapping()->RemoveLast();

            // fix up the softclips
            auto e3 = m->add_edit();
            e3->set_to_length(0);
            e3->set_to_length(3);
            e3->set_sequence("GAT");
            auto e4 = m2->mutable_edit(0);
            e4->set_from_length(0);
            e4->set_to_length(18);
            e4->set_sequence("CAAGTTCTGCTTTCTGCG");

            min_size = 1;
            score_fraction = 0.1;

            separation = 1;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 3);

            separation = 0;
            supps = identify_supplementaries(alns, read_coverage, separation, score_fraction, min_size);
            REQUIRE(supps.size() == 0);
        }
    }
}

TEST_CASE("Unaligned sequences survive round-trip to GAF", "[alignment][gaf]") {
    VG graph;
    Alignment aln;
    {
        aln.set_sequence("TACACTTAC");
        aln.set_name("unaligned");
    }

    gafkluge::GafRecord gaf = alignment_to_gaf(graph, aln, nullptr, true, false, false);
    SECTION("GAF record is correct") {
        REQUIRE(gaf.query_name == aln.name());
        REQUIRE(gaf.query_length == aln.sequence().size());
        REQUIRE(gaf.query_start == 0);
        REQUIRE(gaf.query_end == aln.sequence().size());
        REQUIRE(gaf.path.size() == 0);
        REQUIRE(gaf.opt_fields.find("cs") != gaf.opt_fields.end());
        std::string expected_cs = "+" + aln.sequence();
        REQUIRE(gaf.opt_fields["cs"].second == expected_cs);
    }

    Alignment roundtrip;
    gaf_to_alignment(graph, gaf, roundtrip);
    SECTION("Round-tripped alignment is correct") {
        REQUIRE(roundtrip.name() == aln.name());
        REQUIRE(roundtrip.sequence() == aln.sequence());
        REQUIRE(roundtrip.path().mapping_size() == 0);
    }
}

TEST_CASE("Alignments can be left/right-aligned", "[alignment]") {

    auto paths_equivalent = [](const Path& path1, const Path& path2) {
        REQUIRE(path1.mapping_size() == path2.mapping_size());
        for (size_t i = 0; i < path1.mapping_size(); ++i) {
            const auto& mapping1 = path1.mapping(i);
            const auto& mapping2 = path2.mapping(i);
            REQUIRE(mapping1.position().node_id() == mapping2.position().node_id());
            REQUIRE(mapping1.position().is_reverse() == mapping2.position().is_reverse());
            REQUIRE(mapping1.position().offset() == mapping2.position().offset());
            REQUIRE(mapping1.edit_size() == mapping2.edit_size());
            for (size_t j = 0; j < mapping1.edit_size(); ++j) {
                const auto& edit1 = mapping1.edit(j);
                const auto& edit2 = mapping2.edit(j);
                REQUIRE(edit1.from_length() == edit2.from_length());
                REQUIRE(edit1.to_length() == edit2.to_length());
                REQUIRE(edit1.sequence() == edit2.sequence());
            }
        }
    };

    auto check_left_right_adj_consistency = [&](const Alignment& aln, const HandleGraph& graph) {
        
        auto get_length = [&](nid_t node_id) -> int64_t { return graph.get_length(graph.get_handle(node_id)); };

        for (bool left_adj : {true, false}) {

            auto fwd = aln;
            auto rev = reverse_complement_alignment(fwd, get_length);
            
            normalize_indel_adjustment(fwd, left_adj, graph);
            normalize_indel_adjustment(rev, !left_adj, graph);

            reverse_complement_alignment_in_place(&rev, get_length);

            paths_equivalent(fwd.path(), rev.path());
        }
    };

    SECTION("Within a single node") {

        bdsg::HashGraph graph;

        auto h = graph.create_handle("TACACACT");

        Alignment aln;
        aln.set_sequence("TACACT");
        auto path = aln.mutable_path();

        auto m = path->add_mapping();
        auto p = m->mutable_position();
        p->set_node_id(graph.get_id(h));
        p->set_is_reverse(false);
        p->set_offset(0);

        auto e1 = m->add_edit();
        e1->set_from_length(3);
        e1->set_to_length(3);

        auto e2 = m->add_edit();
        e2->set_from_length(2);
        e2->set_to_length(0);

        auto e3 = m->add_edit();
        e3->set_from_length(3);
        e3->set_to_length(3);

        SECTION("Consistency") {
            check_left_right_adj_consistency(aln, graph);
        }

        SECTION("Deletion, left"){

            normalize_indel_adjustment(aln, true, graph);

            const auto& path = aln.path();
            REQUIRE(path.mapping_size() == 1);
            const auto& mapping = path.mapping(0);
            REQUIRE(mapping.edit_size() == 3);
            REQUIRE(make_pos_t(mapping.position()) == make_pos_t(graph.get_id(h), false, 0));
            REQUIRE(mapping.edit(0).from_length() == 1);
            REQUIRE(mapping.edit(0).to_length() == 1);
            REQUIRE(mapping.edit(0).sequence().empty());
            REQUIRE(mapping.edit(1).from_length() == 2);
            REQUIRE(mapping.edit(1).to_length() == 0);
            REQUIRE(mapping.edit(1).sequence().empty());
            REQUIRE(mapping.edit(2).from_length() == 5);
            REQUIRE(mapping.edit(2).to_length() == 5);
            REQUIRE(mapping.edit(2).sequence().empty());
        }

        SECTION("Deletion, right") {

            normalize_indel_adjustment(aln, false, graph);

            const auto& path = aln.path();
            REQUIRE(path.mapping_size() == 1);
            const auto& mapping = path.mapping(0);
            REQUIRE(mapping.edit_size() == 3);
            REQUIRE(make_pos_t(mapping.position()) == make_pos_t(graph.get_id(h), false, 0));
            REQUIRE(mapping.edit(0).from_length() == 5);
            REQUIRE(mapping.edit(0).to_length() == 5);
            REQUIRE(mapping.edit(0).sequence().empty());
            REQUIRE(mapping.edit(1).from_length() == 2);
            REQUIRE(mapping.edit(1).to_length() == 0);
            REQUIRE(mapping.edit(1).sequence().empty());
            REQUIRE(mapping.edit(2).from_length() == 1);
            REQUIRE(mapping.edit(2).to_length() == 1);
            REQUIRE(mapping.edit(2).sequence().empty());
        }

        // switch to insertion
        aln.set_sequence("TACACACACT");

        e1->set_from_length(5);
        e1->set_to_length(5);

        e2->set_from_length(0);
        e2->set_to_length(2);
        e2->set_sequence("AC");

        SECTION("Consistency") {
            check_left_right_adj_consistency(aln, graph);
        }

        SECTION("Insertion, left"){

            

            normalize_indel_adjustment(aln, true, graph);

            const auto& path = aln.path();
            REQUIRE(path.mapping_size() == 1);
            const auto& mapping = path.mapping(0);
            REQUIRE(mapping.edit_size() == 3);
            REQUIRE(make_pos_t(mapping.position()) == make_pos_t(graph.get_id(h), false, 0));
            REQUIRE(mapping.edit(0).from_length() == 1);
            REQUIRE(mapping.edit(0).to_length() == 1);
            REQUIRE(mapping.edit(0).sequence().empty());
            REQUIRE(mapping.edit(1).from_length() == 0);
            REQUIRE(mapping.edit(1).to_length() == 2);
            REQUIRE(mapping.edit(1).sequence() == "AC");
            REQUIRE(mapping.edit(2).from_length() == 7);
            REQUIRE(mapping.edit(2).to_length() == 7);
            REQUIRE(mapping.edit(2).sequence().empty());
        }

        SECTION("Insertion, right"){

            normalize_indel_adjustment(aln, false, graph);

            const auto& path = aln.path();
            REQUIRE(path.mapping_size() == 1);
            const auto& mapping = path.mapping(0);
            REQUIRE(mapping.edit_size() == 3);
            REQUIRE(make_pos_t(mapping.position()) == make_pos_t(graph.get_id(h), false, 0));
            REQUIRE(mapping.edit(0).from_length() == 7);
            REQUIRE(mapping.edit(0).to_length() == 7);
            REQUIRE(mapping.edit(0).sequence().empty());
            REQUIRE(mapping.edit(1).from_length() == 0);
            REQUIRE(mapping.edit(1).to_length() == 2);
            REQUIRE(mapping.edit(1).sequence() == "AC");
            REQUIRE(mapping.edit(2).from_length() == 1);
            REQUIRE(mapping.edit(2).to_length() == 1);
            REQUIRE(mapping.edit(2).sequence().empty());
        }
    }

    SECTION("Across multiple nodes") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("TTTTACGACGACG");
        auto h2 = graph.create_handle("A");
        auto h3 = graph.create_handle("CG");
        auto h4 = graph.create_handle("AC");
        auto h5 = graph.create_handle("GAC");
        auto h6 = graph.create_handle("G");
        auto h7 = graph.create_handle("ACGTTTT");

        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        graph.create_edge(h5, h6);
        graph.create_edge(h6, h7);

        SECTION("Insertion"){

            Alignment aln;
            aln.set_sequence("TTACGACGACGACGACGACGACGACGACGTT");

            auto path = aln.mutable_path();

            auto m1 = path->add_mapping();
            auto p1 = m1->mutable_position();
            p1->set_node_id(graph.get_id(h1));
            p1->set_is_reverse(false);
            p1->set_offset(2);
            auto e11 = m1->add_edit();
            e11->set_from_length(5);
            e11->set_to_length(5);
            auto e12 = m1->add_edit();
            e12->set_from_length(0);
            e12->set_to_length(6);
            auto e13 = m1->add_edit();
            e13->set_from_length(6);
            e13->set_to_length(6);

            auto m2 = path->add_mapping();
            auto p2 = m2->mutable_position();
            p2->set_node_id(graph.get_id(h2));
            p2->set_is_reverse(false);
            p2->set_offset(0);
            auto e2 = m2->add_edit();
            e2->set_from_length(1);
            e2->set_to_length(1);

            auto m3 = path->add_mapping();
            auto p3 = m3->mutable_position();
            p3->set_node_id(graph.get_id(h3));
            p3->set_is_reverse(false);
            p3->set_offset(0);
            auto e3 = m3->add_edit();
            e3->set_from_length(2);
            e3->set_to_length(2);

            auto m4 = path->add_mapping();
            auto p4 = m4->mutable_position();
            p4->set_node_id(graph.get_id(h4));
            p4->set_is_reverse(false);
            p4->set_offset(0);
            auto e4 = m4->add_edit();
            e4->set_from_length(2);
            e4->set_to_length(2);

            auto m5 = path->add_mapping();
            auto p5 = m5->mutable_position();
            p5->set_node_id(graph.get_id(h5));
            p5->set_is_reverse(false);
            p5->set_offset(0);
            auto e5 = m5->add_edit();
            e5->set_from_length(3);
            e5->set_to_length(3);

            auto m6 = path->add_mapping();
            auto p6 = m6->mutable_position();
            p6->set_node_id(graph.get_id(h6));
            p6->set_is_reverse(false);
            p6->set_offset(0);
            auto e6 = m6->add_edit();
            e6->set_from_length(1);
            e6->set_to_length(1);

            auto m7 = path->add_mapping();
            auto p7 = m7->mutable_position();
            p7->set_node_id(graph.get_id(h7));
            p7->set_is_reverse(false);
            p7->set_offset(0);
            auto e7 = m7->add_edit();
            e7->set_from_length(5);
            e7->set_to_length(5);

            SECTION("Consistency") {
                check_left_right_adj_consistency(aln, graph);
            }

            normalize_indel_adjustment(aln, false, graph);

            REQUIRE(path->mapping_size() == 7);

            REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 2));
            REQUIRE(path->mapping(0).edit_size() == 1);
            REQUIRE(path->mapping(0).edit(0).from_length() == 11);
            REQUIRE(path->mapping(0).edit(0).to_length() == 11);
            REQUIRE(path->mapping(0).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(1).position()) == pos_t(graph.get_id(h2), false, 0));
            REQUIRE(path->mapping(1).edit_size() == 1);
            REQUIRE(path->mapping(1).edit(0).from_length() == 1);
            REQUIRE(path->mapping(1).edit(0).to_length() == 1);
            REQUIRE(path->mapping(1).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(2).position()) == pos_t(graph.get_id(h3), false, 0));
            REQUIRE(path->mapping(2).edit_size() == 1);
            REQUIRE(path->mapping(2).edit(0).from_length() == 2);
            REQUIRE(path->mapping(2).edit(0).to_length() == 2);
            REQUIRE(path->mapping(2).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(3).position()) == pos_t(graph.get_id(h4), false, 0));
            REQUIRE(path->mapping(3).edit_size() == 1);
            REQUIRE(path->mapping(3).edit(0).from_length() == 2);
            REQUIRE(path->mapping(3).edit(0).to_length() == 2);
            REQUIRE(path->mapping(3).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(4).position()) == pos_t(graph.get_id(h5), false, 0));
            REQUIRE(path->mapping(4).edit_size() == 1);
            REQUIRE(path->mapping(4).edit(0).from_length() == 3);
            REQUIRE(path->mapping(4).edit(0).to_length() == 3);
            REQUIRE(path->mapping(4).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(5).position()) == pos_t(graph.get_id(h6), false, 0));
            REQUIRE(path->mapping(5).edit_size() == 1);
            REQUIRE(path->mapping(5).edit(0).from_length() == 1);
            REQUIRE(path->mapping(5).edit(0).to_length() == 1);
            REQUIRE(path->mapping(5).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(6).position()) == pos_t(graph.get_id(h7), false, 0));
            REQUIRE(path->mapping(6).edit_size() == 3);
            REQUIRE(path->mapping(6).edit(0).from_length() == 3);
            REQUIRE(path->mapping(6).edit(0).to_length() == 3);
            REQUIRE(path->mapping(6).edit(0).sequence().empty());
            REQUIRE(path->mapping(6).edit(1).from_length() == 0);
            REQUIRE(path->mapping(6).edit(1).to_length() == 6);
            REQUIRE(path->mapping(6).edit(1).sequence() == "ACGACG");
            REQUIRE(path->mapping(6).edit(2).from_length() == 2);
            REQUIRE(path->mapping(6).edit(2).to_length() == 2);
            REQUIRE(path->mapping(6).edit(2).sequence().empty());

            
        }

        SECTION("Deletion") {
            
            Alignment aln;
            aln.set_sequence("TTACGACGACGACGACGTT");

            auto path = aln.mutable_path();

            auto m1 = path->add_mapping();
            auto p1 = m1->mutable_position();
            p1->set_node_id(graph.get_id(h1));
            p1->set_is_reverse(false);
            p1->set_offset(2);
            auto e11 = m1->add_edit();
            e11->set_from_length(3);
            e11->set_to_length(3);
            auto e12 = m1->add_edit();
            e12->set_from_length(6);
            e12->set_to_length(0);
            auto e13 = m1->add_edit();
            e13->set_from_length(2);
            e13->set_to_length(2);

            auto m2 = path->add_mapping();
            auto p2 = m2->mutable_position();
            p2->set_node_id(graph.get_id(h2));
            p2->set_is_reverse(false);
            p2->set_offset(0);
            auto e2 = m2->add_edit();
            e2->set_from_length(1);
            e2->set_to_length(1);

            auto m3 = path->add_mapping();
            auto p3 = m3->mutable_position();
            p3->set_node_id(graph.get_id(h3));
            p3->set_is_reverse(false);
            p3->set_offset(0);
            auto e3 = m3->add_edit();
            e3->set_from_length(2);
            e3->set_to_length(2);

            auto m4 = path->add_mapping();
            auto p4 = m4->mutable_position();
            p4->set_node_id(graph.get_id(h4));
            p4->set_is_reverse(false);
            p4->set_offset(0);
            auto e4 = m4->add_edit();
            e4->set_from_length(2);
            e4->set_to_length(2);

            auto m5 = path->add_mapping();
            auto p5 = m5->mutable_position();
            p5->set_node_id(graph.get_id(h5));
            p5->set_is_reverse(false);
            p5->set_offset(0);
            auto e5 = m5->add_edit();
            e5->set_from_length(3);
            e5->set_to_length(3);

            auto m6 = path->add_mapping();
            auto p6 = m6->mutable_position();
            p6->set_node_id(graph.get_id(h6));
            p6->set_is_reverse(false);
            p6->set_offset(0);
            auto e6 = m6->add_edit();
            e6->set_from_length(1);
            e6->set_to_length(1);

            auto m7 = path->add_mapping();
            auto p7 = m7->mutable_position();
            p7->set_node_id(graph.get_id(h7));
            p7->set_is_reverse(false);
            p7->set_offset(0);
            auto e7 = m7->add_edit();
            e7->set_from_length(5);
            e7->set_to_length(5);

            SECTION("Consistency") {
                check_left_right_adj_consistency(aln, graph);
            }

            normalize_indel_adjustment(aln, false, graph);

            REQUIRE(path->mapping_size() == 7);

            REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 2));
            REQUIRE(path->mapping(0).edit_size() == 1);
            REQUIRE(path->mapping(0).edit(0).from_length() == 11);
            REQUIRE(path->mapping(0).edit(0).to_length() == 11);
            REQUIRE(path->mapping(0).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(1).position()) == pos_t(graph.get_id(h2), false, 0));
            REQUIRE(path->mapping(1).edit_size() == 1);
            REQUIRE(path->mapping(1).edit(0).from_length() == 1);
            REQUIRE(path->mapping(1).edit(0).to_length() == 1);
            REQUIRE(path->mapping(1).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(2).position()) == pos_t(graph.get_id(h3), false, 0));
            REQUIRE(path->mapping(2).edit_size() == 1);
            REQUIRE(path->mapping(2).edit(0).from_length() == 2);
            REQUIRE(path->mapping(2).edit(0).to_length() == 2);
            REQUIRE(path->mapping(2).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(3).position()) == pos_t(graph.get_id(h4), false, 0));
            REQUIRE(path->mapping(3).edit_size() == 1);
            REQUIRE(path->mapping(3).edit(0).from_length() == 2);
            REQUIRE(path->mapping(3).edit(0).to_length() == 2);
            REQUIRE(path->mapping(3).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(4).position()) == pos_t(graph.get_id(h5), false, 0));
            REQUIRE(path->mapping(4).edit_size() == 2);
            REQUIRE(path->mapping(4).edit(0).from_length() == 1);
            REQUIRE(path->mapping(4).edit(0).to_length() == 1);
            REQUIRE(path->mapping(4).edit(1).from_length() == 2);
            REQUIRE(path->mapping(4).edit(1).to_length() == 0);
            REQUIRE(path->mapping(4).edit(1).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(5).position()) == pos_t(graph.get_id(h6), false, 0));
            REQUIRE(path->mapping(5).edit_size() == 1);
            REQUIRE(path->mapping(5).edit(0).from_length() == 1);
            REQUIRE(path->mapping(5).edit(0).to_length() == 0);
            REQUIRE(path->mapping(5).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(6).position()) == pos_t(graph.get_id(h7), false, 0));
            REQUIRE(path->mapping(6).edit_size() == 2);
            REQUIRE(path->mapping(6).edit(0).from_length() == 3);
            REQUIRE(path->mapping(6).edit(0).to_length() == 0);
            REQUIRE(path->mapping(6).edit(0).sequence().empty());
            REQUIRE(path->mapping(6).edit(1).from_length() == 2);
            REQUIRE(path->mapping(6).edit(1).to_length() == 2);
            REQUIRE(path->mapping(6).edit(1).sequence().empty());
        }
    }

    SECTION("Does not perform incorrect shifts") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("ACGTACGTACGT");

        SECTION("Unmovable deletion"){

            Alignment aln;
            aln.set_sequence("ACGTCGTACGT");

            auto path = aln.mutable_path();
            auto m = path->add_mapping();
            auto p = m->mutable_position();
            p->set_node_id(graph.get_id(h1));
            p->set_is_reverse(false);
            p->set_offset(0);
            auto e1 = m->add_edit();
            e1->set_from_length(4);
            e1->set_to_length(4);
            auto e2 = m->add_edit();
            e2->set_from_length(1);
            e2->set_to_length(0);
            auto e3 = m->add_edit();
            e3->set_from_length(7);
            e3->set_to_length(7);

            auto aln_copy = aln;

            for (bool left : {true, false}) {
                normalize_indel_adjustment(aln, left, graph);
                paths_equivalent(aln.path(), aln_copy.path());
            }
        }

        SECTION("Unmovable insertion"){

            Alignment aln;
            aln.set_sequence("ACGTACGATACGT");

            auto path = aln.mutable_path();
            auto m = path->add_mapping();
            auto p = m->mutable_position();
            p->set_node_id(graph.get_id(h1));
            p->set_is_reverse(false);
            p->set_offset(0);
            auto e1 = m->add_edit();
            e1->set_from_length(7);
            e1->set_to_length(7);
            auto e2 = m->add_edit();
            e2->set_from_length(0);
            e2->set_to_length(1);
            e2->set_sequence("A");
            auto e3 = m->add_edit();
            e3->set_from_length(5);
            e3->set_to_length(5);

            auto aln_copy = aln;

            for (bool left : {true, false}) {
                normalize_indel_adjustment(aln, left, graph);
                paths_equivalent(aln.path(), aln_copy.path());
            }
        }        
    }

    SECTION("Merge insertions/deletions as adjusting") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("TTTTACGACGACG");
        auto h2 = graph.create_handle("A");
        auto h3 = graph.create_handle("CG");
        auto h4 = graph.create_handle("AC");
        auto h5 = graph.create_handle("GAC");
        auto h6 = graph.create_handle("G");
        auto h7 = graph.create_handle("ACGTTTT");

        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        graph.create_edge(h5, h6);
        graph.create_edge(h6, h7);

        SECTION("Insertion"){

            Alignment aln;
            aln.set_sequence("TTACGACGACGACGACGACGACGACGACGTT");

            auto path = aln.mutable_path();

            auto m1 = path->add_mapping();
            auto p1 = m1->mutable_position();
            p1->set_node_id(graph.get_id(h1));
            p1->set_is_reverse(false);
            p1->set_offset(2);
            auto e11 = m1->add_edit();
            e11->set_from_length(5);
            e11->set_to_length(5);
            auto e12 = m1->add_edit();
            e12->set_from_length(0);
            e12->set_to_length(3);
            e12->set_sequence("ACG");
            auto e13 = m1->add_edit();
            e13->set_from_length(6);
            e13->set_to_length(6);

            auto m2 = path->add_mapping();
            auto p2 = m2->mutable_position();
            p2->set_node_id(graph.get_id(h2));
            p2->set_is_reverse(false);
            p2->set_offset(0);
            auto e2 = m2->add_edit();
            e2->set_from_length(1);
            e2->set_to_length(1);

            auto m3 = path->add_mapping();
            auto p3 = m3->mutable_position();
            p3->set_node_id(graph.get_id(h3));
            p3->set_is_reverse(false);
            p3->set_offset(0);
            auto e31 = m3->add_edit();
            e31->set_from_length(2);
            e31->set_to_length(2);
            auto e32 = m3->add_edit();
            e32->set_from_length(0);
            e32->set_to_length(3);
            e32->set_sequence("ACG");

            auto m4 = path->add_mapping();
            auto p4 = m4->mutable_position();
            p4->set_node_id(graph.get_id(h4));
            p4->set_is_reverse(false);
            p4->set_offset(0);
            auto e4 = m4->add_edit();
            e4->set_from_length(2);
            e4->set_to_length(2);

            auto m5 = path->add_mapping();
            auto p5 = m5->mutable_position();
            p5->set_node_id(graph.get_id(h5));
            p5->set_is_reverse(false);
            p5->set_offset(0);
            auto e5 = m5->add_edit();
            e5->set_from_length(3);
            e5->set_to_length(3);

            auto m6 = path->add_mapping();
            auto p6 = m6->mutable_position();
            p6->set_node_id(graph.get_id(h6));
            p6->set_is_reverse(false);
            p6->set_offset(0);
            auto e6 = m6->add_edit();
            e6->set_from_length(1);
            e6->set_to_length(1);

            auto m7 = path->add_mapping();
            auto p7 = m7->mutable_position();
            p7->set_node_id(graph.get_id(h7));
            p7->set_is_reverse(false);
            p7->set_offset(0);
            auto e7 = m7->add_edit();
            e7->set_from_length(5);
            e7->set_to_length(5);

            auto aln_copy = aln;

            normalize_indel_adjustment(aln, false, graph);

            REQUIRE(path->mapping_size() == 7);

            REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 2));
            REQUIRE(path->mapping(0).edit_size() == 1);
            REQUIRE(path->mapping(0).edit(0).from_length() == 11);
            REQUIRE(path->mapping(0).edit(0).to_length() == 11);
            REQUIRE(path->mapping(0).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(1).position()) == pos_t(graph.get_id(h2), false, 0));
            REQUIRE(path->mapping(1).edit_size() == 1);
            REQUIRE(path->mapping(1).edit(0).from_length() == 1);
            REQUIRE(path->mapping(1).edit(0).to_length() == 1);
            REQUIRE(path->mapping(1).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(2).position()) == pos_t(graph.get_id(h3), false, 0));
            REQUIRE(path->mapping(2).edit_size() == 1);
            REQUIRE(path->mapping(2).edit(0).from_length() == 2);
            REQUIRE(path->mapping(2).edit(0).to_length() == 2);
            REQUIRE(path->mapping(2).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(3).position()) == pos_t(graph.get_id(h4), false, 0));
            REQUIRE(path->mapping(3).edit_size() == 1);
            REQUIRE(path->mapping(3).edit(0).from_length() == 2);
            REQUIRE(path->mapping(3).edit(0).to_length() == 2);
            REQUIRE(path->mapping(3).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(4).position()) == pos_t(graph.get_id(h5), false, 0));
            REQUIRE(path->mapping(4).edit_size() == 1);
            REQUIRE(path->mapping(4).edit(0).from_length() == 3);
            REQUIRE(path->mapping(4).edit(0).to_length() == 3);
            REQUIRE(path->mapping(4).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(5).position()) == pos_t(graph.get_id(h6), false, 0));
            REQUIRE(path->mapping(5).edit_size() == 1);
            REQUIRE(path->mapping(5).edit(0).from_length() == 1);
            REQUIRE(path->mapping(5).edit(0).to_length() == 1);
            REQUIRE(path->mapping(5).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(6).position()) == pos_t(graph.get_id(h7), false, 0));
            REQUIRE(path->mapping(6).edit_size() == 3);
            REQUIRE(path->mapping(6).edit(0).from_length() == 3);
            REQUIRE(path->mapping(6).edit(0).to_length() == 3);
            REQUIRE(path->mapping(6).edit(0).sequence().empty());
            REQUIRE(path->mapping(6).edit(1).from_length() == 0);
            REQUIRE(path->mapping(6).edit(1).to_length() == 6);
            REQUIRE(path->mapping(6).edit(1).sequence() == "ACGACG");
            REQUIRE(path->mapping(6).edit(2).from_length() == 2);
            REQUIRE(path->mapping(6).edit(2).to_length() == 2);
            REQUIRE(path->mapping(6).edit(2).sequence().empty());

            {
                // move the insert to the start of a node
                aln_copy.mutable_path()->mutable_mapping(2)->mutable_edit()->RemoveLast();
                auto e41 = aln_copy.mutable_path()->mutable_mapping(3)->mutable_edit(0);
                auto e42 = aln_copy.mutable_path()->mutable_mapping(3)->add_edit();
                *e42 = *e41;
                e41->set_to_length(3);
                e41->set_from_length(0);
                e41->set_sequence("ACG");

                normalize_indel_adjustment(aln_copy, false, graph);

                paths_equivalent(aln.path(), aln_copy.path());
            }
        }

        SECTION("Deletion") {
            
            Alignment aln;
            aln.set_sequence("TTACGACGACGACGACGTT");

            auto path = aln.mutable_path();

            auto m1 = path->add_mapping();
            auto p1 = m1->mutable_position();
            p1->set_node_id(graph.get_id(h1));
            p1->set_is_reverse(false);
            p1->set_offset(2);
            auto e11 = m1->add_edit();
            e11->set_from_length(9);
            e11->set_to_length(9);
            auto e12 = m1->add_edit();
            e12->set_from_length(2);
            e12->set_to_length(0);

            auto m2 = path->add_mapping();
            auto p2 = m2->mutable_position();
            p2->set_node_id(graph.get_id(h2));
            p2->set_is_reverse(false);
            p2->set_offset(0);
            auto e2 = m2->add_edit();
            e2->set_from_length(1);
            e2->set_to_length(0);

            auto m3 = path->add_mapping();
            auto p3 = m3->mutable_position();
            p3->set_node_id(graph.get_id(h3));
            p3->set_is_reverse(false);
            p3->set_offset(0);
            auto e3 = m3->add_edit();
            e3->set_from_length(2);
            e3->set_to_length(2);

            auto m4 = path->add_mapping();
            auto p4 = m4->mutable_position();
            p4->set_node_id(graph.get_id(h4));
            p4->set_is_reverse(false);
            p4->set_offset(0);
            auto e4 = m4->add_edit();
            e4->set_from_length(2);
            e4->set_to_length(2);

            auto m5 = path->add_mapping();
            auto p5 = m5->mutable_position();
            p5->set_node_id(graph.get_id(h5));
            p5->set_is_reverse(false);
            p5->set_offset(0);
            auto e51 = m5->add_edit();
            e51->set_from_length(2);
            e51->set_to_length(2);
            auto e52 = m5->add_edit();
            e52->set_from_length(1);
            e52->set_to_length(0);

            auto m6 = path->add_mapping();
            auto p6 = m6->mutable_position();
            p6->set_node_id(graph.get_id(h6));
            p6->set_is_reverse(false);
            p6->set_offset(0);
            auto e6 = m6->add_edit();
            e6->set_from_length(1);
            e6->set_to_length(0);

            auto m7 = path->add_mapping();
            auto p7 = m7->mutable_position();
            p7->set_node_id(graph.get_id(h7));
            p7->set_is_reverse(false);
            p7->set_offset(0);
            auto e71 = m7->add_edit();
            e71->set_from_length(1);
            e71->set_to_length(0);
            auto e72 = m7->add_edit();
            e72->set_from_length(4);
            e72->set_to_length(4);

            normalize_indel_adjustment(aln, false, graph);

            REQUIRE(path->mapping_size() == 7);

            REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 2));
            REQUIRE(path->mapping(0).edit_size() == 1);
            REQUIRE(path->mapping(0).edit(0).from_length() == 11);
            REQUIRE(path->mapping(0).edit(0).to_length() == 11);
            REQUIRE(path->mapping(0).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(1).position()) == pos_t(graph.get_id(h2), false, 0));
            REQUIRE(path->mapping(1).edit_size() == 1);
            REQUIRE(path->mapping(1).edit(0).from_length() == 1);
            REQUIRE(path->mapping(1).edit(0).to_length() == 1);
            REQUIRE(path->mapping(1).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(2).position()) == pos_t(graph.get_id(h3), false, 0));
            REQUIRE(path->mapping(2).edit_size() == 1);
            REQUIRE(path->mapping(2).edit(0).from_length() == 2);
            REQUIRE(path->mapping(2).edit(0).to_length() == 2);
            REQUIRE(path->mapping(2).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(3).position()) == pos_t(graph.get_id(h4), false, 0));
            REQUIRE(path->mapping(3).edit_size() == 1);
            REQUIRE(path->mapping(3).edit(0).from_length() == 2);
            REQUIRE(path->mapping(3).edit(0).to_length() == 2);
            REQUIRE(path->mapping(3).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(4).position()) == pos_t(graph.get_id(h5), false, 0));
            REQUIRE(path->mapping(4).edit_size() == 2);
            REQUIRE(path->mapping(4).edit(0).from_length() == 1);
            REQUIRE(path->mapping(4).edit(0).to_length() == 1);
            REQUIRE(path->mapping(4).edit(1).from_length() == 2);
            REQUIRE(path->mapping(4).edit(1).to_length() == 0);
            REQUIRE(path->mapping(4).edit(1).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(5).position()) == pos_t(graph.get_id(h6), false, 0));
            REQUIRE(path->mapping(5).edit_size() == 1);
            REQUIRE(path->mapping(5).edit(0).from_length() == 1);
            REQUIRE(path->mapping(5).edit(0).to_length() == 0);
            REQUIRE(path->mapping(5).edit(0).sequence().empty());

            REQUIRE(make_pos_t(path->mapping(6).position()) == pos_t(graph.get_id(h7), false, 0));
            REQUIRE(path->mapping(6).edit_size() == 2);
            REQUIRE(path->mapping(6).edit(0).from_length() == 3);
            REQUIRE(path->mapping(6).edit(0).to_length() == 0);
            REQUIRE(path->mapping(6).edit(0).sequence().empty());
            REQUIRE(path->mapping(6).edit(1).from_length() == 2);
            REQUIRE(path->mapping(6).edit(1).to_length() == 2);
            REQUIRE(path->mapping(6).edit(1).sequence().empty());
        }
    }

    SECTION("Insertions and deletions can cancel each other out") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("TTTTACGACGACG");
        auto h2 = graph.create_handle("A");
        auto h3 = graph.create_handle("CG");
        auto h4 = graph.create_handle("AC");
        auto h5 = graph.create_handle("GAC");
        auto h6 = graph.create_handle("G");
        auto h7 = graph.create_handle("ACGTTTT");

        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        graph.create_edge(h5, h6);
        graph.create_edge(h6, h7);

        Alignment aln;
        aln.set_sequence("TTACGACGACGACGACGACGACGTT");

        auto path = aln.mutable_path();

        auto m1 = path->add_mapping();
        auto p1 = m1->mutable_position();
        p1->set_node_id(graph.get_id(h1));
        p1->set_is_reverse(false);
        p1->set_offset(2);
        auto e11 = m1->add_edit();
        e11->set_from_length(5);
        e11->set_to_length(5);
        auto e12 = m1->add_edit();
        e12->set_from_length(6);
        e12->set_to_length(0);

        auto m2 = path->add_mapping();
        auto p2 = m2->mutable_position();
        p2->set_node_id(graph.get_id(h2));
        p2->set_is_reverse(false);
        p2->set_offset(0);
        auto e2 = m2->add_edit();
        e2->set_from_length(1);
        e2->set_to_length(1);

        auto m3 = path->add_mapping();
        auto p3 = m3->mutable_position();
        p3->set_node_id(graph.get_id(h3));
        p3->set_is_reverse(false);
        p3->set_offset(0);
        auto e3 = m3->add_edit();
        e3->set_from_length(2);
        e3->set_to_length(2);

        auto m4 = path->add_mapping();
        auto p4 = m4->mutable_position();
        p4->set_node_id(graph.get_id(h4));
        p4->set_is_reverse(false);
        p4->set_offset(0);
        auto e4 = m4->add_edit();
        e4->set_from_length(2);
        e4->set_to_length(2);

        auto m5 = path->add_mapping();
        auto p5 = m5->mutable_position();
        p5->set_node_id(graph.get_id(h5));
        p5->set_is_reverse(false);
        p5->set_offset(0);
        auto e51 = m5->add_edit();
        e51->set_from_length(0);
        e51->set_to_length(6);
        e51->set_sequence("GACGAC");
        auto e52 = m5->add_edit();
        e52->set_from_length(3);
        e52->set_to_length(3);

        auto m6 = path->add_mapping();
        auto p6 = m6->mutable_position();
        p6->set_node_id(graph.get_id(h6));
        p6->set_is_reverse(false);
        p6->set_offset(0);
        auto e6 = m6->add_edit();
        e6->set_from_length(1);
        e6->set_to_length(1);

        auto m7 = path->add_mapping();
        auto p7 = m7->mutable_position();
        p7->set_node_id(graph.get_id(h7));
        p7->set_is_reverse(false);
        p7->set_offset(0);
        auto e7 = m7->add_edit();
        e7->set_from_length(5);
        e7->set_to_length(5);

        check_left_right_adj_consistency(aln, graph);

        normalize_indel_adjustment(aln, false, graph);

        REQUIRE(path->mapping_size() == 7);

        REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 2));
        REQUIRE(path->mapping(0).edit_size() == 1);
        REQUIRE(path->mapping(0).edit(0).from_length() == 11);
        REQUIRE(path->mapping(0).edit(0).to_length() == 11);
        REQUIRE(path->mapping(0).edit(0).sequence().empty());

        REQUIRE(make_pos_t(path->mapping(1).position()) == pos_t(graph.get_id(h2), false, 0));
        REQUIRE(path->mapping(1).edit_size() == 1);
        REQUIRE(path->mapping(1).edit(0).from_length() == 1);
        REQUIRE(path->mapping(1).edit(0).to_length() == 1);
        REQUIRE(path->mapping(1).edit(0).sequence().empty());

        REQUIRE(make_pos_t(path->mapping(2).position()) == pos_t(graph.get_id(h3), false, 0));
        REQUIRE(path->mapping(2).edit_size() == 1);
        REQUIRE(path->mapping(2).edit(0).from_length() == 2);
        REQUIRE(path->mapping(2).edit(0).to_length() == 2);
        REQUIRE(path->mapping(2).edit(0).sequence().empty());

        REQUIRE(make_pos_t(path->mapping(3).position()) == pos_t(graph.get_id(h4), false, 0));
        REQUIRE(path->mapping(3).edit_size() == 1);
        REQUIRE(path->mapping(3).edit(0).from_length() == 2);
        REQUIRE(path->mapping(3).edit(0).to_length() == 2);
        REQUIRE(path->mapping(3).edit(0).sequence().empty());

        REQUIRE(make_pos_t(path->mapping(4).position()) == pos_t(graph.get_id(h5), false, 0));
        REQUIRE(path->mapping(4).edit_size() == 1);
        REQUIRE(path->mapping(4).edit(0).from_length() == 3);
        REQUIRE(path->mapping(4).edit(0).to_length() == 3);
        REQUIRE(path->mapping(4).edit(0).sequence().empty());

        REQUIRE(make_pos_t(path->mapping(5).position()) == pos_t(graph.get_id(h6), false, 0));
        REQUIRE(path->mapping(5).edit_size() == 1);
        REQUIRE(path->mapping(5).edit(0).from_length() == 1);
        REQUIRE(path->mapping(5).edit(0).to_length() == 1);
        REQUIRE(path->mapping(5).edit(0).sequence().empty());

        REQUIRE(make_pos_t(path->mapping(6).position()) == pos_t(graph.get_id(h7), false, 0));
        REQUIRE(path->mapping(6).edit_size() == 1);
        REQUIRE(path->mapping(6).edit(0).from_length() == 5);
        REQUIRE(path->mapping(6).edit(0).to_length() == 5);
        REQUIRE(path->mapping(6).edit(0).sequence().empty());
    }

    SECTION("Left/right-aligning indels shifting appropriately handles softclips") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("TTAAAAAATT");
        auto h2 = graph.create_handle("TTAC");
        auto h3 = graph.create_handle("ACAC");
        auto h4 = graph.create_handle("ACTT");

        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);

        Alignment aln;

        aln.set_sequence("GGGAAAAAAAAGGG");

        auto path = aln.mutable_path();

        auto m = path->add_mapping();
        auto p = m->mutable_position();
        p->set_node_id(graph.get_id(h1));
        p->set_is_reverse(false);
        p->set_offset(2);

        auto e1 = m->add_edit();
        e1->set_from_length(0);
        e1->set_to_length(3);
        e1->set_sequence("GGG");

        auto e2 = m->add_edit();
        e2->set_from_length(3);
        e2->set_to_length(3);

        auto e3 = m->add_edit();
        e3->set_from_length(0);
        e3->set_to_length(2);
        e3->set_sequence("AA");

        auto e4 = m->add_edit();
        e4->set_from_length(3);
        e4->set_to_length(3);

        auto e5 = m->add_edit();
        e5->set_from_length(0);
        e5->set_to_length(3);
        e5->set_sequence("GGG");

        SECTION("Shift insertion right") {

            normalize_indel_adjustment(aln, false, graph);

            REQUIRE(path->mapping_size() == 1);

            REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 2));
            REQUIRE(path->mapping(0).edit_size() == 3);

            REQUIRE(path->mapping(0).edit(0).from_length() == 0);
            REQUIRE(path->mapping(0).edit(0).to_length() == 3);
            REQUIRE(path->mapping(0).edit(0).sequence() == "GGG");

            REQUIRE(path->mapping(0).edit(1).from_length() == 6);
            REQUIRE(path->mapping(0).edit(1).to_length() == 6);
            REQUIRE(path->mapping(0).edit(1).sequence().empty());

            REQUIRE(path->mapping(0).edit(2).from_length() == 0);
            REQUIRE(path->mapping(0).edit(2).to_length() == 5);
            REQUIRE(path->mapping(0).edit(2).sequence() == "AAGGG");
        }

        SECTION("Shift insertion left") {

            normalize_indel_adjustment(aln, true, graph);

            REQUIRE(path->mapping_size() == 1);

            REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 2));
            REQUIRE(path->mapping(0).edit_size() == 3);

            REQUIRE(path->mapping(0).edit(0).from_length() == 0);
            REQUIRE(path->mapping(0).edit(0).to_length() == 5);
            REQUIRE(path->mapping(0).edit(0).sequence() == "GGGAA");

            REQUIRE(path->mapping(0).edit(1).from_length() == 6);
            REQUIRE(path->mapping(0).edit(1).to_length() == 6);
            REQUIRE(path->mapping(0).edit(1).sequence().empty());

            REQUIRE(path->mapping(0).edit(2).from_length() == 0);
            REQUIRE(path->mapping(0).edit(2).to_length() == 3);
            REQUIRE(path->mapping(0).edit(2).sequence() == "GGG");
        }

        // convert to a deletino
        aln.set_sequence("GGGAAAAGGG");

        e2->set_from_length(2);
        e2->set_to_length(2);

        e3->set_from_length(2);
        e3->set_to_length(0);
        e3->set_sequence("");

        e4->set_from_length(2);
        e4->set_to_length(2);

        SECTION("Shift deletion to the right") {

            normalize_indel_adjustment(aln, false, graph);

            REQUIRE(path->mapping_size() == 1);

            REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 2));
            REQUIRE(path->mapping(0).edit_size() == 3);

            REQUIRE(path->mapping(0).edit(0).from_length() == 0);
            REQUIRE(path->mapping(0).edit(0).to_length() == 3);
            REQUIRE(path->mapping(0).edit(0).sequence() == "GGG");

            REQUIRE(path->mapping(0).edit(1).from_length() == 4);
            REQUIRE(path->mapping(0).edit(1).to_length() == 4);
            REQUIRE(path->mapping(0).edit(1).sequence().empty());

            REQUIRE(path->mapping(0).edit(2).from_length() == 0);
            REQUIRE(path->mapping(0).edit(2).to_length() == 3);
            REQUIRE(path->mapping(0).edit(2).sequence() == "GGG");
        } 

        SECTION("Shift deletion to the left") {

            normalize_indel_adjustment(aln, true, graph);

            REQUIRE(path->mapping_size() == 1);

            REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h1), false, 4));
            REQUIRE(path->mapping(0).edit_size() == 3);

            REQUIRE(path->mapping(0).edit(0).from_length() == 0);
            REQUIRE(path->mapping(0).edit(0).to_length() == 3);
            REQUIRE(path->mapping(0).edit(0).sequence() == "GGG");

            REQUIRE(path->mapping(0).edit(1).from_length() == 4);
            REQUIRE(path->mapping(0).edit(1).to_length() == 4);
            REQUIRE(path->mapping(0).edit(1).sequence().empty());

            REQUIRE(path->mapping(0).edit(2).from_length() == 0);
            REQUIRE(path->mapping(0).edit(2).to_length() == 3);
            REQUIRE(path->mapping(0).edit(2).sequence() == "GGG");
        }

        SECTION("Shifting deletion over multiple nodes") {

            Alignment aln;

            aln.set_sequence("GGGACACGGG");

            auto path = aln.mutable_path();

            auto m1 = path->add_mapping();
            auto p1 = m1->mutable_position();
            p1->set_node_id(graph.get_id(h2));
            p1->set_is_reverse(false);
            p1->set_offset(2);

            auto e11 = m1->add_edit();
            e11->set_from_length(0);
            e11->set_to_length(3);
            e11->set_sequence("GGG");

            auto e12 = m1->add_edit();
            e12->set_from_length(2);
            e12->set_to_length(2);

            auto m2 = path->add_mapping();
            auto p2 = m2->mutable_position();
            p2->set_node_id(graph.get_id(h3));
            p2->set_is_reverse(false);
            p2->set_offset(0);

            auto e2 = m2->add_edit();
            e2->set_from_length(4);
            e2->set_to_length(0);

            auto m3 = path->add_mapping();
            auto p3 = m3->mutable_position();
            p3->set_node_id(graph.get_id(h4));
            p3->set_is_reverse(false);
            p3->set_offset(0);

            auto e31 = m3->add_edit();
            e31->set_from_length(2);
            e31->set_to_length(2);

            auto e32 = m3->add_edit();
            e32->set_from_length(0);
            e32->set_to_length(3);
            e32->set_sequence("GGG");

            SECTION("To the right") {

                normalize_indel_adjustment(aln, false, graph);

                REQUIRE(path->mapping_size() == 2);

                REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h2), false, 2));
                REQUIRE(path->mapping(0).edit_size() == 2);

                REQUIRE(path->mapping(0).edit(0).from_length() == 0);
                REQUIRE(path->mapping(0).edit(0).to_length() == 3);
                REQUIRE(path->mapping(0).edit(0).sequence() == "GGG");

                REQUIRE(path->mapping(0).edit(1).from_length() == 2);
                REQUIRE(path->mapping(0).edit(1).to_length() == 2);
                REQUIRE(path->mapping(0).edit(1).sequence().empty());

                REQUIRE(make_pos_t(path->mapping(1).position()) == pos_t(graph.get_id(h3), false, 0));
                REQUIRE(path->mapping(1).edit_size() == 2);

                REQUIRE(path->mapping(1).edit(0).from_length() == 2);
                REQUIRE(path->mapping(1).edit(0).to_length() == 2);

                REQUIRE(path->mapping(1).edit(1).from_length() == 0);
                REQUIRE(path->mapping(1).edit(1).to_length() == 3);
                REQUIRE(path->mapping(1).edit(1).sequence() == "GGG");
            }

            SECTION("To the left") {

                normalize_indel_adjustment(aln, true, graph);

                REQUIRE(path->mapping_size() == 2);

                REQUIRE(make_pos_t(path->mapping(0).position()) == pos_t(graph.get_id(h3), false, 2));
                REQUIRE(path->mapping(0).edit_size() == 2);

                REQUIRE(path->mapping(0).edit(0).from_length() == 0);
                REQUIRE(path->mapping(0).edit(0).to_length() == 3);
                REQUIRE(path->mapping(0).edit(0).sequence() == "GGG");

                REQUIRE(path->mapping(0).edit(1).from_length() == 2);
                REQUIRE(path->mapping(0).edit(1).to_length() == 2);
                REQUIRE(path->mapping(0).edit(1).sequence().empty());

                REQUIRE(make_pos_t(path->mapping(1).position()) == pos_t(graph.get_id(h4), false, 0));
                REQUIRE(path->mapping(1).edit_size() == 2);

                REQUIRE(path->mapping(1).edit(0).from_length() == 2);
                REQUIRE(path->mapping(1).edit(0).to_length() == 2);

                REQUIRE(path->mapping(1).edit(1).from_length() == 0);
                REQUIRE(path->mapping(1).edit(1).to_length() == 3);
                REQUIRE(path->mapping(1).edit(1).sequence() == "GGG");
            }
        }
    }
}

TEST_CASE("Left alignment properly handles cases involving adjacent insertions and deletions") {

    SECTION("Runs are consolidated before shifting") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("GGAG");

        Alignment aln;

        aln.set_sequence("GGTTAG");

        auto path = aln.mutable_path();

        auto m = path->add_mapping();
        auto p = m->mutable_position();
        p->set_node_id(graph.get_id(h1));
        p->set_is_reverse(false);
        p->set_offset(0);

        auto e1 = m->add_edit();
        e1->set_from_length(2);
        e1->set_to_length(2);

        auto e2 = m->add_edit();
        e2->set_from_length(0);
        e2->set_to_length(1);
        e2->set_sequence("T");

        auto e3 = m->add_edit();
        e3->set_from_length(1);
        e3->set_to_length(0);

        auto e4 = m->add_edit();
        e4->set_from_length(0);
        e4->set_to_length(2);
        e4->set_sequence("TA");

        auto e5 = m->add_edit();
        e5->set_from_length(1);
        e5->set_to_length(1);

        normalize_indel_adjustment(aln, false, graph);

        REQUIRE(m->edit_size() == 3);
        REQUIRE(m->edit(0).from_length() == 2);
        REQUIRE(m->edit(0).to_length() == 2);
        REQUIRE(m->edit(0).sequence() == "");
        REQUIRE(m->edit(1).from_length() == 0);
        REQUIRE(m->edit(1).to_length() == 2);
        REQUIRE(m->edit(1).sequence() == "TT");
        REQUIRE(m->edit(2).from_length() == 2);
        REQUIRE(m->edit(2).to_length() == 2);
        REQUIRE(m->edit(2).sequence() == "");
    }

    SECTION("Runs are consolidated before shifting with opposite configuration") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("GGTTAG");

        Alignment aln;

        aln.set_sequence("GGAG");

        auto path = aln.mutable_path();

        auto m = path->add_mapping();
        auto p = m->mutable_position();
        p->set_node_id(graph.get_id(h1));
        p->set_is_reverse(false);
        p->set_offset(0);

        auto e1 = m->add_edit();
        e1->set_from_length(2);
        e1->set_to_length(2);

        auto e2 = m->add_edit();
        e2->set_from_length(1);
        e2->set_to_length(0);

        auto e3 = m->add_edit();
        e3->set_from_length(0);
        e3->set_to_length(1);
        e3->set_sequence("A");

        auto e4 = m->add_edit();
        e4->set_from_length(2);
        e4->set_to_length(0);

        auto e5 = m->add_edit();
        e5->set_from_length(1);
        e5->set_to_length(1);

        normalize_indel_adjustment(aln, false, graph);

        REQUIRE(m->edit_size() == 3);
        REQUIRE(m->edit(0).from_length() == 2);
        REQUIRE(m->edit(0).to_length() == 2);
        REQUIRE(m->edit(0).sequence() == "");
        REQUIRE(m->edit(1).from_length() == 2);
        REQUIRE(m->edit(1).to_length() == 0);
        REQUIRE(m->edit(1).sequence() == "");
        REQUIRE(m->edit(2).from_length() == 2);
        REQUIRE(m->edit(2).to_length() == 2);
        REQUIRE(m->edit(2).sequence() == "");
    }


    SECTION("Left alignment can continue a released shift after cancellation") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("ATGTGGG");

        Alignment aln;

        aln.set_sequence("ATGTGG");

        auto path = aln.mutable_path();

        auto m = path->add_mapping();
        auto p = m->mutable_position();
        p->set_node_id(graph.get_id(h1));
        p->set_is_reverse(false);
        p->set_offset(0);

        auto e1 = m->add_edit();
        e1->set_from_length(1);
        e1->set_to_length(1);

        auto e2 = m->add_edit();
        e2->set_from_length(2);
        e2->set_to_length(0);

        auto e3 = m->add_edit();
        e3->set_from_length(2);
        e3->set_to_length(2);

        auto e4 = m->add_edit();
        e4->set_from_length(0);
        e4->set_to_length(1);
        e4->set_sequence("T");

        auto e5 = m->add_edit();
        e5->set_from_length(2);
        e5->set_to_length(2);

        normalize_indel_adjustment(aln, false, graph);

        REQUIRE(m->edit_size() == 1);
        REQUIRE(m->edit(0).from_length() == 6);
        REQUIRE(m->edit(0).to_length() == 6);
        REQUIRE(m->edit(0).sequence() == "");
    }
}

TEST_CASE("Left alignment produces valid results on difficult cases") {

    SECTION("Empirical case 1") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("T", 38183332);
        auto h2 = graph.create_handle("TC", 38183333);
        auto h3 = graph.create_handle("C", 38183334);
        auto h4 = graph.create_handle("A", 38183335);
        auto h5 = graph.create_handle("T", 38183337);
        auto h6 = graph.create_handle("T", 38183338);
        auto h7 = graph.create_handle("C", 38183339);
        auto h8 = graph.create_handle("C", 38183340);
        auto h9 = graph.create_handle("A", 38183341);
        auto h10 = graph.create_handle("T", 38183343);
        auto h11 = graph.create_handle("TC", 38183344);
        auto h12 = graph.create_handle("C", 38183345);
        auto h13 = graph.create_handle("A", 38183347);
        auto h14 = graph.create_handle("T", 38183348);
        auto h15 = graph.create_handle("TCC", 38183349);
        auto h16 = graph.create_handle("A", 38183350);
        auto h17 = graph.create_handle("T", 38183353);
        auto h18 = graph.create_handle("T", 38183355);
        auto h19 = graph.create_handle("C", 38183356);
        auto h20 = graph.create_handle("C", 38183357);
        auto h21 = graph.create_handle("ATT", 38183359);
        auto h22 = graph.create_handle("C", 38183360);
        auto h23 = graph.create_handle("C", 38183361);
        auto h24 = graph.create_handle("A", 38183362);
        auto h25 = graph.create_handle("TTC", 38183364);
        auto h26 = graph.create_handle("C", 38183365);
        auto h27 = graph.create_handle("A", 38183366);
        auto h28 = graph.create_handle("TT", 38183368);
        auto h29 = graph.create_handle("C", 38183369);
        auto h30 = graph.create_handle("CA", 38183371);
        auto h31 = graph.create_handle("C", 38183372);
        auto h32 = graph.create_handle("T", 38183374);
        auto h33 = graph.create_handle("C", 38183375);
        auto h34 = graph.create_handle("G", 38183379);
        auto h35 = graph.create_handle("GG", 38183380);
        auto h36 = graph.create_handle("T", 38183381);
        auto h37 = graph.create_handle("T", 38183382);
        auto h38 = graph.create_handle("G", 38183385);
        auto h39 = graph.create_handle("A", 38183387);
        auto h40 = graph.create_handle("T", 38183388);
        auto h41 = graph.create_handle("T", 38183389);
        auto h42 = graph.create_handle("C", 38183390);
        auto h43 = graph.create_handle("C", 38183391);
        auto h44 = graph.create_handle("A", 38183392);
        auto h45 = graph.create_handle("T", 38183393);
        auto h46 = graph.create_handle("T", 38183395);
        auto h47 = graph.create_handle("C", 38183396);
        auto h48 = graph.create_handle("C", 38183398);
        auto h49 = graph.create_handle("A", 38183400);
        auto h50 = graph.create_handle("T", 38183401);
        auto h51 = graph.create_handle("T", 38183403);
        auto h52 = graph.create_handle("C", 38183404);
        auto h53 = graph.create_handle("A", 38183405);
        auto h54 = graph.create_handle("A", 38183409);
        auto h55 = graph.create_handle("TTC", 38183410);
        auto h56 = graph.create_handle("C", 38183411);
        auto h57 = graph.create_handle("G", 38183412);
        auto h58 = graph.create_handle("T", 38183414);
        auto h59 = graph.create_handle("T", 38183415);
        auto h60 = graph.create_handle("C", 38183416);
        auto h61 = graph.create_handle("C", 38183417);
        auto h62 = graph.create_handle("G", 38183421);
        auto h63 = graph.create_handle("T", 38183423);
        auto h64 = graph.create_handle("T", 38183424);
        auto h65 = graph.create_handle("CC", 38183425);
        auto h66 = graph.create_handle("G", 38183426);
        auto h67 = graph.create_handle("TT", 38183428);
        auto h68 = graph.create_handle("CC", 38183429);
        auto h69 = graph.create_handle("G", 38183430);
        auto h70 = graph.create_handle("T", 38183432);
        auto h71 = graph.create_handle("T", 38183433);
        auto h72 = graph.create_handle("C", 38183434);
        auto h73 = graph.create_handle("C", 38183435);
        auto h74 = graph.create_handle("A", 38183437);
        auto h75 = graph.create_handle("T", 38183438);
        auto h76 = graph.create_handle("T", 38183439);
        auto h77 = graph.create_handle("C", 38183440);
        auto h78 = graph.create_handle("C", 38183442);
        auto h79 = graph.create_handle("A", 38183443);
        auto h80 = graph.create_handle("T", 38183444);
        auto h81 = graph.create_handle("T", 38183445);
        auto h82 = graph.create_handle("C", 38183446);
        auto h83 = graph.create_handle("C", 38183449);
        auto h84 = graph.create_handle("A", 38183451);
        auto h85 = graph.create_handle("TT", 38183452);
        auto h86 = graph.create_handle("C", 38183453);
        auto h87 = graph.create_handle("C", 38183454);
        auto h88 = graph.create_handle("A", 38183455);
        auto h89 = graph.create_handle("T", 38183457);
        auto h90 = graph.create_handle("T", 38183458);
        auto h91 = graph.create_handle("T", 38183459);
        auto h92 = graph.create_handle("C", 38183462);
        auto h93 = graph.create_handle("AT", 38183463);
        auto h94 = graph.create_handle("T", 38183464);
        auto h95 = graph.create_handle("C", 38183466);
        auto h96 = graph.create_handle("C", 38183467);
        auto h97 = graph.create_handle("A", 38183469);
        auto h98 = graph.create_handle("T", 38183470);
        auto h99 = graph.create_handle("T", 38183472);
        auto h100 = graph.create_handle("G", 38183473);
        auto h101 = graph.create_handle("CA", 38183475);
        auto h102 = graph.create_handle("TTC", 38183476);
        auto h103 = graph.create_handle("C", 38183477);
        auto h104 = graph.create_handle("A", 38183479);
        auto h105 = graph.create_handle("C", 38183480);
        auto h106 = graph.create_handle("T", 38183482);
        auto h107 = graph.create_handle("C", 38183483);
        auto h108 = graph.create_handle("G", 38183485);
        auto h109 = graph.create_handle("G", 38183486);
        auto h110 = graph.create_handle("G", 38183487);
        auto h111 = graph.create_handle("T", 38183488);
        auto h112 = graph.create_handle("TG", 38183489);
        auto h113 = graph.create_handle("A", 38183490);
        auto h114 = graph.create_handle("T", 38183491);
        auto h115 = graph.create_handle("T", 38183492);
        auto h116 = graph.create_handle("CC", 38183493);
        auto h117 = graph.create_handle("A", 38183494);
        auto h118 = graph.create_handle("A", 38183496);
        auto h119 = graph.create_handle("A", 38183497);
        auto h120 = graph.create_handle("G", 38183498);
        auto h121 = graph.create_handle("C", 38183499);
        auto h122 = graph.create_handle("A", 38183500);
        auto h123 = graph.create_handle("T", 38183503);
        auto h124 = graph.create_handle("T", 38183504);
        auto h125 = graph.create_handle("C", 38183506);
        auto h126 = graph.create_handle("C", 38183507);
        auto h127 = graph.create_handle("A", 38183509);
        auto h128 = graph.create_handle("T", 38183510);
        auto h129 = graph.create_handle("TC", 38183518);
        auto h130 = graph.create_handle("C", 38183519);
        auto h131 = graph.create_handle("A", 38183521);
        auto h132 = graph.create_handle("T", 38183522);
        auto h133 = graph.create_handle("T", 38183523);
        auto h134 = graph.create_handle("CC", 38183524);
        auto h135 = graph.create_handle("A", 38183525);
        auto h136 = graph.create_handle("T", 38183526);
        auto h137 = graph.create_handle("T", 38183527);
        auto h138 = graph.create_handle("CC", 38183528);
        auto h139 = graph.create_handle("A", 38183529);
        auto h140 = graph.create_handle("T", 38183530);
        auto h141 = graph.create_handle("T", 38183532);
        auto h142 = graph.create_handle("C", 38183533);
        auto h143 = graph.create_handle("C", 38183541);
        auto h144 = graph.create_handle("A", 38183543);
        auto h145 = graph.create_handle("T", 38183544);
        auto h146 = graph.create_handle("T", 38183546);
        auto h147 = graph.create_handle("C", 38183547);
        auto h148 = graph.create_handle("C", 38183548);
        auto h149 = graph.create_handle("A", 38183549);
        auto h150 = graph.create_handle("TT", 38183550);
        auto h151 = graph.create_handle("C", 38183551);
        auto h152 = graph.create_handle("C", 38183552);
        auto h153 = graph.create_handle("A", 38183553);
        auto h154 = graph.create_handle("T", 38183555);
        auto h155 = graph.create_handle("T", 38183558);
        auto h156 = graph.create_handle("C", 38183559);
        auto h157 = graph.create_handle("CATTC", 38183560);
        auto h158 = graph.create_handle("C", 38183561);
        auto h159 = graph.create_handle("A", 38183563);
        auto h160 = graph.create_handle("C", 38183564);
        auto h161 = graph.create_handle("T", 38183567);
        auto h162 = graph.create_handle("C", 38183568);
        auto h163 = graph.create_handle("G", 38183569);
        auto h164 = graph.create_handle("G", 38183572);
        auto h165 = graph.create_handle("G", 38183574);
        auto h166 = graph.create_handle("T", 38183575);
        auto h167 = graph.create_handle("T", 38183576);
        auto h168 = graph.create_handle("G", 38183578);
        auto h169 = graph.create_handle("ATT", 38183579);
        auto h170 = graph.create_handle("C", 38183588);
        auto h171 = graph.create_handle("C", 38183593);
        auto h172 = graph.create_handle("A", 38183594);
        auto h173 = graph.create_handle("T", 38183595);
        auto h174 = graph.create_handle("T", 38183596);
        auto h175 = graph.create_handle("C", 38183598);
        auto h176 = graph.create_handle("C", 38183599);
        auto h177 = graph.create_handle("AT", 38183601);
        auto h178 = graph.create_handle("T", 38183602);
        auto h179 = graph.create_handle("C", 38183603);
        auto h180 = graph.create_handle("C", 38183604);
        auto h181 = graph.create_handle("A", 38183605);
        auto h182 = graph.create_handle("T", 38183606);
        auto h183 = graph.create_handle("TCC", 38183607);
        auto h184 = graph.create_handle("A", 38183608);
        auto h185 = graph.create_handle("T", 38183609);
        auto h186 = graph.create_handle("T", 38183610);
        auto h187 = graph.create_handle("C", 38183612);
        auto h188 = graph.create_handle("TG", 38183613);
        auto h189 = graph.create_handle("GAA", 38183614);
        auto h190 = graph.create_handle("C", 38183617);
        auto h191 = graph.create_handle("AAT", 38183618);
        auto h192 = graph.create_handle("C", 38183619);
        auto h193 = graph.create_handle("C", 38183620);
        auto h194 = graph.create_handle("A", 38183621);
        auto h195 = graph.create_handle("T", 38183622);
        auto h196 = graph.create_handle("TCC", 38183623);
        auto h197 = graph.create_handle("A", 38183624);
        auto h198 = graph.create_handle("T", 38183625);
        auto h199 = graph.create_handle("TC", 38183626);
        auto h200 = graph.create_handle("C", 38183627);
        auto h201 = graph.create_handle("AT", 38183628);
        auto h202 = graph.create_handle("T", 38183629);
        auto h203 = graph.create_handle("C", 38183631);
        auto h204 = graph.create_handle("C", 38183634);
        auto h205 = graph.create_handle("A", 38183635);
        auto h206 = graph.create_handle("T", 38183637);
        auto h207 = graph.create_handle("T", 38183638);
        auto h208 = graph.create_handle("C", 38183639);
        auto h209 = graph.create_handle("C", 38183640);
        auto h210 = graph.create_handle("A", 38183642);
        auto h211 = graph.create_handle("C", 38183643);
        auto h212 = graph.create_handle("T", 38183646);
        auto h213 = graph.create_handle("C", 38183647);
        auto h214 = graph.create_handle("CATT", 38183673);
        auto h215 = graph.create_handle("C", 38183674);
        auto h216 = graph.create_handle("C", 38183679);
        auto h217 = graph.create_handle("A", 38183680);
        auto h218 = graph.create_handle("T", 38183682);
        auto h219 = graph.create_handle("T", 38183683);
        auto h220 = graph.create_handle("C", 38183684);
        auto h221 = graph.create_handle("C", 38183686);
        auto h222 = graph.create_handle("A", 38183687);
        auto h223 = graph.create_handle("A", 38183689);
        auto h224 = graph.create_handle("T", 38183691);
        auto h225 = graph.create_handle("C", 38183692);
        auto h226 = graph.create_handle("C", 38183693);

        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        graph.create_edge(h5, h6);
        graph.create_edge(h6, h7);
        graph.create_edge(h7, h8);
        graph.create_edge(h8, h9);
        graph.create_edge(h9, h10);
        graph.create_edge(h10, h11);
        graph.create_edge(h11, h12);
        graph.create_edge(h12, h13);
        graph.create_edge(h13, h14);
        graph.create_edge(h14, h15);
        graph.create_edge(h15, h16);
        graph.create_edge(h16, h17);
        graph.create_edge(h17, h18);
        graph.create_edge(h18, h19);
        graph.create_edge(h19, h20);
        graph.create_edge(h20, h21);
        graph.create_edge(h21, h22);
        graph.create_edge(h22, h23);
        graph.create_edge(h23, h24);
        graph.create_edge(h24, h25);
        graph.create_edge(h25, h26);
        graph.create_edge(h26, h27);
        graph.create_edge(h27, h28);
        graph.create_edge(h28, h29);
        graph.create_edge(h29, h30);
        graph.create_edge(h30, h31);
        graph.create_edge(h31, h32);
        graph.create_edge(h32, h33);
        graph.create_edge(h33, h34);
        graph.create_edge(h34, h35);
        graph.create_edge(h35, h36);
        graph.create_edge(h36, h37);
        graph.create_edge(h37, h38);
        graph.create_edge(h38, h39);
        graph.create_edge(h39, h40);
        graph.create_edge(h40, h41);
        graph.create_edge(h41, h42);
        graph.create_edge(h42, h43);
        graph.create_edge(h43, h44);
        graph.create_edge(h44, h45);
        graph.create_edge(h45, h46);
        graph.create_edge(h46, h47);
        graph.create_edge(h47, h48);
        graph.create_edge(h48, h49);
        graph.create_edge(h49, h50);
        graph.create_edge(h50, h51);
        graph.create_edge(h51, h52);
        graph.create_edge(h52, h53);
        graph.create_edge(h53, h54);
        graph.create_edge(h54, h55);
        graph.create_edge(h55, h56);
        graph.create_edge(h56, h57);
        graph.create_edge(h57, h58);
        graph.create_edge(h58, h59);
        graph.create_edge(h59, h60);
        graph.create_edge(h60, h61);
        graph.create_edge(h61, h62);
        graph.create_edge(h62, h63);
        graph.create_edge(h63, h64);
        graph.create_edge(h64, h65);
        graph.create_edge(h65, h66);
        graph.create_edge(h66, h67);
        graph.create_edge(h67, h68);
        graph.create_edge(h68, h69);
        graph.create_edge(h69, h70);
        graph.create_edge(h70, h71);
        graph.create_edge(h71, h72);
        graph.create_edge(h72, h73);
        graph.create_edge(h73, h74);
        graph.create_edge(h74, h75);
        graph.create_edge(h75, h76);
        graph.create_edge(h76, h77);
        graph.create_edge(h77, h78);
        graph.create_edge(h78, h79);
        graph.create_edge(h79, h80);
        graph.create_edge(h80, h81);
        graph.create_edge(h81, h82);
        graph.create_edge(h82, h83);
        graph.create_edge(h83, h84);
        graph.create_edge(h84, h85);
        graph.create_edge(h85, h86);
        graph.create_edge(h86, h87);
        graph.create_edge(h87, h88);
        graph.create_edge(h88, h89);
        graph.create_edge(h89, h90);
        graph.create_edge(h90, h91);
        graph.create_edge(h91, h92);
        graph.create_edge(h92, h93);
        graph.create_edge(h93, h94);
        graph.create_edge(h94, h95);
        graph.create_edge(h95, h96);
        graph.create_edge(h96, h97);
        graph.create_edge(h97, h98);
        graph.create_edge(h98, h99);
        graph.create_edge(h99, h100);
        graph.create_edge(h100, h101);
        graph.create_edge(h101, h102);
        graph.create_edge(h102, h103);
        graph.create_edge(h103, h104);
        graph.create_edge(h104, h105);
        graph.create_edge(h105, h106);
        graph.create_edge(h106, h107);
        graph.create_edge(h107, h108);
        graph.create_edge(h108, h109);
        graph.create_edge(h109, h110);
        graph.create_edge(h110, h111);
        graph.create_edge(h111, h112);
        graph.create_edge(h112, h113);
        graph.create_edge(h113, h114);
        graph.create_edge(h114, h115);
        graph.create_edge(h115, h116);
        graph.create_edge(h116, h117);
        graph.create_edge(h117, h118);
        graph.create_edge(h118, h119);
        graph.create_edge(h119, h120);
        graph.create_edge(h120, h121);
        graph.create_edge(h121, h122);
        graph.create_edge(h122, h123);
        graph.create_edge(h123, h124);
        graph.create_edge(h124, h125);
        graph.create_edge(h125, h126);
        graph.create_edge(h126, h127);
        graph.create_edge(h127, h128);
        graph.create_edge(h128, h129);
        graph.create_edge(h129, h130);
        graph.create_edge(h130, h131);
        graph.create_edge(h131, h132);
        graph.create_edge(h132, h133);
        graph.create_edge(h133, h134);
        graph.create_edge(h134, h135);
        graph.create_edge(h135, h136);
        graph.create_edge(h136, h137);
        graph.create_edge(h137, h138);
        graph.create_edge(h138, h139);
        graph.create_edge(h139, h140);
        graph.create_edge(h140, h141);
        graph.create_edge(h141, h142);
        graph.create_edge(h142, h143);
        graph.create_edge(h143, h144);
        graph.create_edge(h144, h145);
        graph.create_edge(h145, h146);
        graph.create_edge(h146, h147);
        graph.create_edge(h147, h148);
        graph.create_edge(h148, h149);
        graph.create_edge(h149, h150);
        graph.create_edge(h150, h151);
        graph.create_edge(h151, h152);
        graph.create_edge(h152, h153);
        graph.create_edge(h153, h154);
        graph.create_edge(h154, h155);
        graph.create_edge(h155, h156);
        graph.create_edge(h156, h157);
        graph.create_edge(h157, h158);
        graph.create_edge(h158, h159);
        graph.create_edge(h159, h160);
        graph.create_edge(h160, h161);
        graph.create_edge(h161, h162);
        graph.create_edge(h162, h163);
        graph.create_edge(h163, h164);
        graph.create_edge(h164, h165);
        graph.create_edge(h165, h166);
        graph.create_edge(h166, h167);
        graph.create_edge(h167, h168);
        graph.create_edge(h168, h169);
        graph.create_edge(h169, h170);
        graph.create_edge(h170, h171);
        graph.create_edge(h171, h172);
        graph.create_edge(h172, h173);
        graph.create_edge(h173, h174);
        graph.create_edge(h174, h175);
        graph.create_edge(h175, h176);
        graph.create_edge(h176, h177);
        graph.create_edge(h177, h178);
        graph.create_edge(h178, h179);
        graph.create_edge(h179, h180);
        graph.create_edge(h180, h181);
        graph.create_edge(h181, h182);
        graph.create_edge(h182, h183);
        graph.create_edge(h183, h184);
        graph.create_edge(h184, h185);
        graph.create_edge(h185, h186);
        graph.create_edge(h186, h187);
        graph.create_edge(h187, h188);
        graph.create_edge(h188, h189);
        graph.create_edge(h189, h190);
        graph.create_edge(h190, h191);
        graph.create_edge(h191, h192);
        graph.create_edge(h192, h193);
        graph.create_edge(h193, h194);
        graph.create_edge(h194, h195);
        graph.create_edge(h195, h196);
        graph.create_edge(h196, h197);
        graph.create_edge(h197, h198);
        graph.create_edge(h198, h199);
        graph.create_edge(h199, h200);
        graph.create_edge(h200, h201);
        graph.create_edge(h201, h202);
        graph.create_edge(h202, h203);
        graph.create_edge(h203, h204);
        graph.create_edge(h204, h205);
        graph.create_edge(h205, h206);
        graph.create_edge(h206, h207);
        graph.create_edge(h207, h208);
        graph.create_edge(h208, h209);
        graph.create_edge(h209, h210);
        graph.create_edge(h210, h211);
        graph.create_edge(h211, h212);
        graph.create_edge(h212, h213);
        graph.create_edge(h213, h214);
        graph.create_edge(h214, h215);
        graph.create_edge(h215, h216);
        graph.create_edge(h216, h217);
        graph.create_edge(h217, h218);
        graph.create_edge(h218, h219);
        graph.create_edge(h219, h220);
        graph.create_edge(h220, h221);
        graph.create_edge(h221, h222);
        graph.create_edge(h222, h223);
        graph.create_edge(h223, h224);
        graph.create_edge(h224, h225);
        graph.create_edge(h225, h226);


        string alignment_string = R"({"path": {"mapping": [{"edit": [{"sequence": "CGAGTGCTGGGGAATGTAAT", "to_length": 20}, {"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183693"}, "rank": "1"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183692"}, "rank": "2"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183691"}, "rank": "3"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "38183689"}, "rank": "4"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183687"}, "rank": "5"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183686"}, "rank": "6"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183684"}, "rank": "7"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183683"}, "rank": "8"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183682"}, "rank": "9"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183680"}, "rank": "10"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183679"}, "rank": "11"}, {"edit": [{"from_length": 1, "sequence": "C", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183674"}, "rank": "12"}, {"edit": [{"from_length": 4, "to_length": 4}], "position": {"is_reverse": true, "node_id": "38183673"}, "rank": "13"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183647"}, "rank": "14"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183646"}, "rank": "15"}, {"edit": [{"from_length": 1, "sequence": "A", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183643"}, "rank": "16"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183642"}, "rank": "17"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183640"}, "rank": "18"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183639"}, "rank": "19"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183638"}, "rank": "20"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183637"}, "rank": "21"}, {"edit": [{"from_length": 1, "to_length": 1}, {"sequence": "CATCC", "to_length": 5}], "position": {"is_reverse": true, "node_id": "38183635"}, "rank": "22"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183634"}, "rank": "1"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183631"}, "rank": "2"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183629"}, "rank": "3"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183628"}, "rank": "4"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183627"}, "rank": "5"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183626"}, "rank": "6"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183625"}, "rank": "7"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183624"}, "rank": "8"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"is_reverse": true, "node_id": "38183623"}, "rank": "9"}, {"edit": [{"from_length": 1, "sequence": "G", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183622"}, "rank": "10"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183621"}, "rank": "11"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183620"}, "rank": "12"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183619"}, "rank": "13"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "38183618"}, "rank": "14"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183617"}, "rank": "15"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "38183614"}, "rank": "16"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "38183613"}, "rank": "17"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183612"}, "rank": "18"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183610"}, "rank": "19"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183609"}, "rank": "20"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183608"}, "rank": "21"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"is_reverse": true, "node_id": "38183607"}, "rank": "22"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183606"}, "rank": "23"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183605"}, "rank": "24"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183604"}, "rank": "25"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183603"}, "rank": "26"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183602"}, "rank": "27"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183601"}, "rank": "28"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183599"}, "rank": "29"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183598"}, "rank": "30"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183596"}, "rank": "31"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183595"}, "rank": "32"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183594"}, "rank": "33"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183593"}, "rank": "34"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183588"}, "rank": "35"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"is_reverse": true, "node_id": "38183579"}, "rank": "36"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183578"}, "rank": "37"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183576"}, "rank": "38"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183575"}, "rank": "39"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183574"}, "rank": "40"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183572"}, "rank": "41"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183569"}, "rank": "42"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183568"}, "rank": "43"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183567"}, "rank": "44"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183564"}, "rank": "45"}, {"edit": [{"from_length": 1, "sequence": "C", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183563"}, "rank": "46"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183561"}, "rank": "47"}, {"edit": [{"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 4, "to_length": 4}], "position": {"is_reverse": true, "node_id": "38183560"}, "rank": "48"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183559"}, "rank": "49"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183558"}, "rank": "50"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183555"}, "rank": "51"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183553"}, "rank": "52"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183552"}, "rank": "53"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183551"}, "rank": "54"}, {"edit": [{"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "G", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183550"}, "rank": "55"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183549"}, "rank": "56"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183548"}, "rank": "57"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183547"}, "rank": "58"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183546"}, "rank": "59"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183544"}, "rank": "60"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183543"}, "rank": "61"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183541"}, "rank": "62"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183533"}, "rank": "63"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183532"}, "rank": "64"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183530"}, "rank": "65"}, {"edit": [{"from_length": 1, "sequence": "A", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183529"}, "rank": "66"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183528"}, "rank": "67"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183527"}, "rank": "68"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183526"}, "rank": "69"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183525"}, "rank": "70"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183524"}, "rank": "71"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183523"}, "rank": "72"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183522"}, "rank": "73"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183521"}, "rank": "74"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183519"}, "rank": "75"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183518"}, "rank": "76"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183510"}, "rank": "77"}, {"edit": [{"sequence": "CAACCCGAA", "to_length": 9}, {"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183509"}, "rank": "1"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183507"}, "rank": "2"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183506"}, "rank": "3"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183504"}, "rank": "4"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183503"}, "rank": "5"}, {"edit": [{"from_length": 1, "sequence": "T", "to_length": 1}], "position": {"node_id": "38183500"}, "rank": "6"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183499"}, "rank": "7"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "38183498"}, "rank": "8"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "38183497"}, "rank": "9"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "38183496"}, "rank": "10"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183494"}, "rank": "11"}, {"edit": [{"from_length": 1, "to_length": 1}, {"sequence": "TAATGGAGAGTAAGGGAGTTGAATA", "to_length": 25}, {"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183493"}, "rank": "12"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183492"}, "rank": "13"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183491"}, "rank": "14"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183490"}, "rank": "15"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183489"}, "rank": "16"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183488"}, "rank": "17"}, {"edit": [{"from_length": 1, "sequence": "T", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183487"}, "rank": "18"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183486"}, "rank": "19"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183485"}, "rank": "20"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183483"}, "rank": "21"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183482"}, "rank": "22"}, {"edit": [{"from_length": 1, "sequence": "A", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183480"}, "rank": "23"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183479"}, "rank": "24"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183477"}, "rank": "25"}, {"edit": [{"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183476"}, "rank": "26"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183475"}, "rank": "27"}, {"edit": [{"from_length": 1, "sequence": "G", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183473"}, "rank": "28"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183472"}, "rank": "29"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183470"}, "rank": "30"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183469"}, "rank": "31"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183467"}, "rank": "32"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183466"}, "rank": "33"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183464"}, "rank": "34"}, {"edit": [{"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "C", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183463"}, "rank": "35"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183462"}, "rank": "36"}, {"edit": [{"from_length": 1, "sequence": "G", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183459"}, "rank": "37"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183458"}, "rank": "38"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183457"}, "rank": "39"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183455"}, "rank": "40"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183454"}, "rank": "41"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183453"}, "rank": "42"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183452"}, "rank": "43"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183451"}, "rank": "44"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183449"}, "rank": "45"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183446"}, "rank": "46"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183445"}, "rank": "47"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183444"}, "rank": "48"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183443"}, "rank": "49"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183442"}, "rank": "50"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183440"}, "rank": "51"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183439"}, "rank": "52"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183438"}, "rank": "53"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183437"}, "rank": "54"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183435"}, "rank": "55"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183434"}, "rank": "56"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183433"}, "rank": "57"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183432"}, "rank": "58"}, {"edit": [{"from_length": 1, "sequence": "T", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183430"}, "rank": "59"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183429"}, "rank": "60"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183428"}, "rank": "61"}, {"edit": [{"from_length": 1, "sequence": "T", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183426"}, "rank": "62"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183425"}, "rank": "63"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183424"}, "rank": "64"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183423"}, "rank": "65"}, {"edit": [{"from_length": 1, "sequence": "T", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183421"}, "rank": "66"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183417"}, "rank": "67"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183416"}, "rank": "68"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183415"}, "rank": "69"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183414"}, "rank": "70"}, {"edit": [{"from_length": 1, "sequence": "T", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183412"}, "rank": "71"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183411"}, "rank": "72"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"is_reverse": true, "node_id": "38183410"}, "rank": "73"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183409"}, "rank": "74"}, {"edit": [{"from_length": 1, "sequence": "G", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183405"}, "rank": "75"}, {"edit": [{"from_length": 1, "sequence": "C", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183404"}, "rank": "76"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183403"}, "rank": "77"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183401"}, "rank": "78"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183400"}, "rank": "79"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183398"}, "rank": "80"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183396"}, "rank": "81"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183395"}, "rank": "82"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183393"}, "rank": "83"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183392"}, "rank": "84"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183391"}, "rank": "85"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183390"}, "rank": "86"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183389"}, "rank": "87"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183388"}, "rank": "88"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183387"}, "rank": "89"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183385"}, "rank": "90"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183382"}, "rank": "91"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183381"}, "rank": "92"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183380"}, "rank": "93"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183379"}, "rank": "94"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183375"}, "rank": "95"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183374"}, "rank": "96"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183372"}, "rank": "97"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183371"}, "rank": "98"}, {"edit": [{"from_length": 1, "sequence": "C", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183369"}, "rank": "99"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183368"}, "rank": "100"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183366"}, "rank": "101"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183365"}, "rank": "102"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"is_reverse": true, "node_id": "38183364"}, "rank": "103"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183362"}, "rank": "104"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183361"}, "rank": "105"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183360"}, "rank": "106"}, {"edit": [{"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 1, "sequence": "A", "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183359"}, "rank": "107"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183357"}, "rank": "108"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183356"}, "rank": "109"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183355"}, "rank": "110"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183353"}, "rank": "111"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183350"}, "rank": "112"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"is_reverse": true, "node_id": "38183349"}, "rank": "113"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183348"}, "rank": "114"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183347"}, "rank": "115"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183345"}, "rank": "116"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183344"}, "rank": "117"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183343"}, "rank": "118"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183341"}, "rank": "119"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183340"}, "rank": "120"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183339"}, "rank": "121"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183338"}, "rank": "122"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183337"}, "rank": "123"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183335"}, "rank": "124"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "38183334"}, "rank": "125"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"is_reverse": true, "node_id": "38183333"}, "rank": "126"}, {"edit": [{"from_length": 1, "to_length": 1}, {"sequence": "ACTACCCAAATGGAATGGAATGTAATGGACTGTAAAGGAATTGAATAGAATCAATCCGAATGTAATGGAATGGAATGGAACGGAATG", "to_length": 87}], "position": {"is_reverse": true, "node_id": "38183332"}, "rank": "127"}]}, "sequence": "CGAGTGCTGGGGAATGTAATGGAATGGAATGCAATGGAATGGAATCATCCGGAATGGAATGGAGTGGAATGGAATGGAATGGAATGGAATGGAATGGAATCAACCCGAGCGCAATGGAATGGAGTGGAATGGAAAGGAATGGAATGGAACAACCCGAATGGAATGGAATGTAATGGAGAGTAAGGGAGTTGAATAGAATCAATCCGAATGTAATGGAATGGAACGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGCAATGGAATGGAATCAACCCGAGTGCAATGGAATGGAGAGGAATGGAATGGAATGGAATGGAAACTACCCAAATGGAATGGAATGTAATGGACTGTAAAGGAATTGAATAGAATCAATCCGAATGTAATGGAATGGAATGGAACGGAATG"})";

        Alignment aln;
        json2pb(aln, alignment_string.c_str(), alignment_string.size());

        for (bool left : {false, true}) {
            normalize_indel_adjustment(aln, left, graph);
            auto validity = alignment_is_valid(aln, &graph, true);
            cerr << validity.bad_mapping_index << ", " << validity.bad_edit_index << ", " << validity.bad_read_position << ", " << validity.message << endl;
            REQUIRE((bool) validity);
        }
    }

    SECTION("Empirical case 2") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("TTCCATTCCATTCCA");
        auto h2 = graph.create_handle("TTCCATTCCATTCCATTCC");
        auto h3 = graph.create_handle("ATTC");
        auto h4 = graph.create_handle("CA");
        auto h5 = graph.create_handle("G");
        auto h6 = graph.create_handle("T");
        auto h7 = graph.create_handle("TG");
        auto h8 = graph.create_handle("ATTCCATT");
        auto h9 = graph.create_handle("G");
        auto h10 = graph.create_handle("CATTCC");
        auto h11 = graph.create_handle("A");
        auto h12 = graph.create_handle("TTC");
        auto h13 = graph.create_handle("C");
        auto h14 = graph.create_handle("A");
        auto h15 = graph.create_handle("TTCCAT");
        auto h16 = graph.create_handle("TCCATTCC");
        auto h17 = graph.create_handle("G");
        auto h18 = graph.create_handle("TTCC");
        auto h19 = graph.create_handle("A");
        auto h20 = graph.create_handle("TTC");
        auto h21 = graph.create_handle("G");
        auto h22 = graph.create_handle("T");
        auto h23 = graph.create_handle("GAA");
        auto h24 = graph.create_handle("CATTCCATTC");
        auto h25 = graph.create_handle("C");
        auto h26 = graph.create_handle("ATTCCAT");
        auto h27 = graph.create_handle("TCCAT");
        auto h28 = graph.create_handle("TCCATTCCTTTCCACTCGGGTTGATTCCATT");
        auto h29 = graph.create_handle("CCAT");
        auto h30 = graph.create_handle("TCCATTCCTTTCCATTCCATTCCATTCCGTTCCACTCGGCTTGATTCCATTCCATTCCATTCCATTTTTTCCAATCCACTCGGGTTGATTCCATTCTATTCCATTCCATTCCAGTTG");

        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h3, h4);
        graph.create_edge(h4, h5);
        graph.create_edge(h5, h6);
        graph.create_edge(h6, h7);
        graph.create_edge(h7, h8);
        graph.create_edge(h8, h9);
        graph.create_edge(h9, h10);
        graph.create_edge(h10, h11);
        graph.create_edge(h11, h12);
        graph.create_edge(h12, h13);
        graph.create_edge(h13, h14);
        graph.create_edge(h14, h15);
        graph.create_edge(h15, h16);
        graph.create_edge(h16, h17);
        graph.create_edge(h17, h18);
        graph.create_edge(h18, h19);
        graph.create_edge(h19, h20);
        graph.create_edge(h20, h21);
        graph.create_edge(h21, h22);
        graph.create_edge(h22, h23);
        graph.create_edge(h23, h24);
        graph.create_edge(h24, h25);
        graph.create_edge(h25, h26);
        graph.create_edge(h26, h27);
        graph.create_edge(h27, h28);
        graph.create_edge(h28, h29);
        graph.create_edge(h29, h30);


        string alignment_string = R"({"path": {"mapping": [{"edit": [{"sequence": "ATTCCATTCCACTTGGCGTGATTCA", "to_length": 25}, {"from_length": 9, "to_length": 9}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 1, "to_length": 1}], "position": {"node_id": "1", "offset": "4"}, "rank": "1"}, {"edit": [{"from_length": 6, "to_length": 6}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 11, "to_length": 11}], "position": {"node_id": "2"}, "rank": "2"}, {"edit": [{"from_length": 4, "to_length": 4}], "position": {"node_id": "3"}, "rank": "3"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "4"}, "rank": "4"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "5"}, "rank": "5"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "6"}, "rank": "6"}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "7"}, "rank": "7"}, {"edit": [{"from_length": 4, "to_length": 4}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 3, "to_length": 3}], "position": {"node_id": "8"}, "rank": "8"}, {"edit": [{"from_length": 1, "sequence": "C", "to_length": 1}], "position": {"node_id": "9"}, "rank": "9"}, {"edit": [{"from_length": 5, "to_length": 5}, {"from_length": 1, "sequence": "T", "to_length": 1}], "position": {"node_id": "10"}, "rank": "10"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "11"}, "rank": "11"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "12"}, "rank": "12"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "13"}, "rank": "13"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "14"}, "rank": "14"}, {"edit": [{"from_length": 5, "to_length": 5}, {"from_length": 1}], "position": {"node_id": "15"}, "rank": "15"}, {"edit": [{"from_length": 8}], "position": {"node_id": "16"}, "rank": "16"}, {"edit": [{"from_length": 1}], "position": {"node_id": "17"}, "rank": "17"}, {"edit": [{"from_length": 4, "to_length": 4}], "position": {"node_id": "18"}, "rank": "1"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "19"}, "rank": "2"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"node_id": "20"}, "rank": "3"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "21"}, "rank": "4"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"is_reverse": true, "node_id": "22"}, "rank": "5"}, {"edit": [{"from_length": 3, "to_length": 3}], "position": {"is_reverse": true, "node_id": "23"}, "rank": "6"}, {"edit": [{"from_length": 10, "to_length": 10}], "position": {"node_id": "24"}, "rank": "7"}, {"edit": [{"from_length": 1, "to_length": 1}], "position": {"node_id": "25"}, "rank": "8"}, {"edit": [{"from_length": 3, "to_length": 3}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 3, "to_length": 3}], "position": {"node_id": "26"}, "rank": "9"}, {"edit": [{"from_length": 3, "to_length": 3}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 1, "to_length": 1}], "position": {"node_id": "27"}, "rank": "10"}, {"edit": [{"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 5, "to_length": 5}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 5, "to_length": 5}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 13, "to_length": 13}], "position": {"node_id": "28"}, "rank": "11"}, {"edit": [{"from_length": 4, "to_length": 4}], "position": {"node_id": "29"}, "rank": "12"}, {"edit": [{"from_length": 19, "to_length": 19}, {"sequence": "GTTCATTCAATTCCATTCTATTCCATTCCAATGGATTCCATTCCACTTCAT", "to_length": 51}], "position": {"node_id": "30"}, "rank": "13"}]},  "sequence": "ATTCCATTCCACTTGGCGTGATTCAATTCCATTCAATTCCATATCATTCCATTCCATTCCAGTTGATTCTATTCCATTCTATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTGCATTCCTTTCAATTCCATTTCACTCAGGTTGATTCCATTCCATTCCATTCCTTTCCATTCCAGTTCATTCAATTCCATTCTATTCCATTCCAATGGATTCCATTCCACTTCAT"})";

        Alignment aln;
        json2pb(aln, alignment_string.c_str(), alignment_string.size());

        for (bool left : {false, true}) {
            normalize_indel_adjustment(aln, left, graph);
            auto validity = alignment_is_valid(aln, &graph, true);
            REQUIRE((bool) validity);
        }
    }

    SECTION("Empirical case 3") {

        bdsg::HashGraph graph;

        auto h1 = graph.create_handle("GCATTCTCAGTAAGAGGTTGTTATGTTTGCAATCATCTCACAGACTTGAACCTTTCTTTTGATAGAGCAGTGTTGAAATGCACTTTTTGTAGAATCTGCATGTGTTCATTTGGAGCGCTTTGTTGCCAGTGGTGGAAAAAGAAATATCTTCATATAAAAACTAGATAGAAGCATTTGCAGAAACCCCTTTGTGATGTGTGTGTTCAATTCACAGAGTTTAACCTTTCTTTAGATAGAGCAGTTTTGAAACACTGCCTTTGTAGGATCTGCTTGTGGATATTTGGAGCTCTTTAAGGAATTCGTTGTAAACGGGATAACTTCACATACAAACTAGACAGAAGCATTCTCAGAAACTGCTTTGTGACGTGTGCATTCAACTCACAGAGTTGAACCTTCTTTTTGAGAGATCAATTTTGAAACAGTCTGATTGTATTATCTGCAAGTGGATAATTGGAGCGACTTGAGGCCTAAGATGGACAAGGAGATAGCTTCACATACAAACTAGACAGAAGAATTCTCAGAAACTTCTTTGTGATGTGTGCATTCAACTCATTGTCTAGAACAATTCTTTTGATAGAGCAGTGTTGAAACAGCCTTCTTGTGGAATCTGCAAGTGTTCATTTGGAGAGCTTTGTTGCCTATTCTGGAAAAAGAAACATCTTCACATAAAAACTAGACAGAAGCATTCTCAGAGACTCCTTTTTGATATGTGTGTTCCATTCACAGAGTTGAACCATTCTTTTGATAGAGCAGTTTTGAAACACTGTTTTTGTAGCATCTGCTTGTGGATATT");
        auto h2 = graph.create_handle("GA");
        auto h3 = graph.create_handle("AGCTCTTTCTGGAATTCGTTGTAAACGGGATATCTTCACATACAAACTAGACAGAAGCATTCTCAGAAACTACTTTGTGATGTGTGCATTCAACTCACAGAGTTGAGCCTTCCTTTTGAGAGAGCAGTTTTGAAACAGTCTTTTTGTAGTATCTGCAAGTGGATATTAGGAGCGATTTGAGGCCTATGATGAAAAAGGAAATATCTTCACATACAACGAGACAGAAGCATTCTCAGACACTGCTTTGTGATGTGTGCATTCAACTCACAGTGTTGAAGCTTCCTTTTGAGAGAGCAGTTTTGAAACAGTCTTTTTGTAGTGGCTGCAAGTGGATATTTGGAGCGATTTCAGGCCTGTGATGGAAAAGGAAATATCTTCACATAAAAACTAGACAGAATCATTCTCAGAAACTGCTTTGTGATGTGTGCATTCACCTCACAGAGTGGAACCGTTCTTTTGATAGATCAGTTTTGAAACAGTCTTTTTGTAGGATCTGCAAGTGTTCATGTGGAGTGCTTTGTAGCCAAAGATGAGAACGGAAATATCTTCACGTAAAAACTAGACAGAAGCGTTCTCAGGAACTTCATTGAGATGTGTGCATTCAACTGACAGAGTTGAAACTGTCTTTTGATAGAGCAGTATTGAAACACTCCGTTTGTAGTATCTGGTTGTGGATATTTGGAACTCTTTGAGGAATTCATTGGAAACCGGTATCTTCACATAAAATATAGACCCAAGCATTCTCTGAAATTTATTTCTGATGTGTGCATTCAACTAACAGACGTGAACCTTTTCTTTGACAGAGCAGTGATGAAACACACTTTTTGTGGAATCTGCAAATTTTCATTTCGTGTGCTTTGTTGCCTACGGTTGAAAATAATATCTTCACATAAAAACTAGACAGAAGCATTCTCAGAAACTCTTTGGATGTGTCTGTTCAATTCACAGAGTTAAACCTTTCTTTTTATAGAGCAGTTTTGAAACACTACTTTTGTGGAATGTGCTTGTGGAGATTTGGAACTGT");

        graph.create_edge(h1, h2);
        graph.create_edge(h2, h3);

        string alignment_string = R"({"path": {"mapping": [{"edit": [{"sequence": "AAGCAGATTCTACAAAAGCAG", "to_length": 21}, {"from_length": 17, "to_length": 17}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 7, "to_length": 7}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 9, "to_length": 9}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 25, "to_length": 25}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 3, "to_length": 3}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 6, "to_length": 6}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 4, "to_length": 4}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 15, "to_length": 15}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 15, "to_length": 15}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 14, "to_length": 14}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 1}, {"from_length": 4, "to_length": 4}], "position": {"is_reverse": true, "node_id": "3", "offset": "887"}, "rank": "1"}, {"edit": [{"from_length": 2, "to_length": 2}, {"sequence": "A", "to_length": 1}], "position": {"is_reverse": true, "node_id": "2"}, "rank": "2"}, {"edit": [{"from_length": 11, "to_length": 11}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 5, "to_length": 5}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 8, "to_length": 8}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 20, "to_length": 20}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 9, "to_length": 9}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 9, "to_length": 9}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 11, "to_length": 11}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 5, "to_length": 5}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 4, "to_length": 4}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 18, "to_length": 18}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 18, "to_length": 18}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 66, "to_length": 66}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 3, "to_length": 3}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 4, "to_length": 4}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 19, "to_length": 19}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 7, "to_length": 7}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 15, "to_length": 15}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 8, "to_length": 8}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 5, "to_length": 5}, {"from_length": 1}, {"from_length": 7, "to_length": 7}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 14, "to_length": 14}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 6, "to_length": 6}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 3, "to_length": 3}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 3, "to_length": 3}, {"from_length": 1}, {"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 4, "to_length": 4}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 3, "to_length": 3}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 1, "to_length": 1}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 11, "to_length": 11}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 4, "to_length": 4}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 5, "to_length": 5}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 1, "sequence": "T", "to_length": 1}, {"from_length": 13, "to_length": 13}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 2, "to_length": 2}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 10, "to_length": 10}, {"from_length": 1, "sequence": "G", "to_length": 1}, {"from_length": 3, "to_length": 3}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 15, "to_length": 15}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 6, "to_length": 6}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 4, "to_length": 4}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 16, "to_length": 16}, {"from_length": 1, "sequence": "A", "to_length": 1}, {"from_length": 3, "to_length": 3}, {"from_length": 1, "sequence": "C", "to_length": 1}, {"from_length": 11, "to_length": 11}], "position": {"is_reverse": true, "node_id": "1"}, "rank": "1"}]}, "sequence": "AAGCAGATTCTACAAAAGCAGTGTTTCAAAACTGCTCTATCAAAAGAAAGGCTCAAGTCTGTGAGTTGAATGCACACATCACCAAGCAGTTTCCAAGAAAGCTTCTGTCTAGTTTTTATGTGAAGATATCCAGTTTACAACGAATTTCACAACAGCTTCAAATATCCACAAACAGATTCTACAAAAGCAGTGTTTCAAAACTGCTCTTTCAAAAGAAAGGTTCAACTTTGTGAATGGAAAACACACATCACAAAGGAGTCTCTGAGAATACTTCTGTCTAGTTTTTATTTGAAGATGTTTCTTTTTCCAGAATAGGCAACAAAGCTCTCCAAATGAACACTTGCAGATTCCACAAAAAGACCGTTTAAGCACTGCTCTATCAAAAGAAATGTTCTACACTGAGAGTTGAATGCACACCTCACAAAGCAGTTTTGAGAATGCTTCTGTCTAGTTTTTACGTGAAGATATTTTCTTTTCACCATATGCCTGAAATCGTTTGAAATATCCACTTGCCGATACTACAAAATGACTGTTTCAAAACTGCTCTCTCAAAAGGAAAGTTCAACTCTGTGAGATGAATGAACACATCACAAAGCAGTTTCTAAGACTGCTTCTGTCT"})";

        Alignment aln;
        json2pb(aln, alignment_string.c_str(), alignment_string.size());

        normalize_indel_adjustment(aln, false, graph);

        auto validity = alignment_is_valid(aln, &graph, true);
        REQUIRE((bool) validity);

    }
}

}
}
