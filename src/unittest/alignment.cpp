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

TEST_CASE("Conversion to GAF removes an unused final node", "[alignment]") {
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

}
}
