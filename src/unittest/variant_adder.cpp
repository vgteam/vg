/**
 * \file 
 * unittest/variant_adder.cpp: test cases for the add-variants-to-existing-graph tool
 */

#include "catch.hpp"
#include "../variant_adder.hpp"
#include "../handle.hpp"

#include "../utility.hpp"
#include "../path.hpp"
#include "vg/io/json2pb.h"

#include <vector>
#include <sstream>
#include <iostream>
#include <regex>

using namespace vg::io;

namespace vg {
namespace unittest {

TEST_CASE( "Files with variants but not samples are rejected", "[variantadder]" ) {

    // We'll work on this tiny VCF
    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	5	rs1337	A	G	29	PASS	.	GT
)";

    // Make a stream out of the data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);

    string graph_json = R"({
        "node": [{"id": 1, "sequence": "GATTACA"}],
        "path": [
            {"name": "ref", "mapping": [
                {"position": {"node_id": 1}, "edit": [{"from_length": 7, "to_length": 7}]}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    
    // Make a VariantAdder
    VariantAdder adder(graph);
    // Fail to add the variants to the graph
    REQUIRE_THROWS(adder.add_variants(&vcf));

}

TEST_CASE( "A SNP can be added", "[variantadder]" ) {

    // We'll work on this tiny VCF
    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
ref	5	rs1337	A	G	29	PASS	.	GT	0/1
)";

    // Make a stream out of the data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);

    string graph_json = R"({
        "node": [{"id": 1, "sequence": "GATTACA"}],
        "path": [
            {"name": "ref", "mapping": [
                {"position": {"node_id": 1}, "edit": [{"from_length": 7, "to_length": 7}]}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make a VariantAdder
    VariantAdder adder(graph);
    // Add the variants to the graph
    adder.add_variants(&vcf);

    SECTION("the graph should have 4 nodes") {
        REQUIRE(graph.size() == 4);
    }
    
    SECTION("the graph should have 4 edges") {
        REQUIRE(graph.edge_count() == 4);
    }

}

TEST_CASE( "A relatively long deletion can be added", "[variantadder]" ) {

    // We'll work on this tiny VCF
    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
ref	5	rs1337	AAAAAAAAAAAAAAAAAAAAA	A	29	PASS	.	GT	0/1
)";

    // Make a stream out of the data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);

    string graph_json = R"({
        "node": [{"id": 1, "sequence": "GCGCAAAAAAAAAAAAAAAAAAAAAGCGC"}],
        "path": [
            {"name": "ref", "mapping": [
                {"position": {"node_id": 1}, "edit": [{"from_length": 29, "to_length": 29}]}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make a VariantAdder
    VariantAdder adder(graph);
    // Add the variants to the graph
    adder.add_variants(&vcf);

    SECTION("the graph should have 3 nodes") {
        REQUIRE(graph.size() == 3);
    }
    
    SECTION("the graph should have 3 edges") {
        REQUIRE(graph.edge_count() == 3);
    }

}

TEST_CASE( "A relatively long deletion can be added across a reversing edge", "[variantadder]" ) {

    // We'll work on this tiny VCF
    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
ref	5	rs1337	AAAAAAAAAAAAAAAAAAAAA	A	29	PASS	.	GT	0/1
)";

    // Make a stream out of the data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);

    string graph_json = R"({
        "node": [
            {"id": 1, "sequence": "GCGCAAAAAAAAAAA"},
            {"id": 2, "sequence": "GCGCTTTTTTTTTT"}
        ],
        "edge": [
            {"from": 2, "to": 1, "to_end": true}
        ],
        "path": [
            {"name": "ref", "mapping": [
                {"position": {"node_id": 1}, "edit": [{"from_length": 15, "to_length": 15}], "rank": 1},
                {"position": {"node_id": 2, "is_reverse": true}, "edit": [{"from_length": 14, "to_length": 14}], "rank": 2}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    SECTION ("should work when the graph is as given") {
    
        // Make a VariantAdder
        VariantAdder adder(graph);
        // Add the variants to the graph
        adder.add_variants(&vcf);

        SECTION("the graph should have 4 nodes") {
            REQUIRE(graph.size() == 4);
        }
        
        SECTION("the graph should have 4 edges") {
            REQUIRE(graph.edge_count() == 4);
        }
        
    }
    
    SECTION ("should work when the graph is atomized") {
    
        handlealgs::chop(graph, 1);
        graph.paths.compact_ranks();
    
        // Make a VariantAdder
        VariantAdder adder(graph);
        // Add the variants to the graph
        adder.add_variants(&vcf);

        SECTION("the graph should have 29 nodes") {
            REQUIRE(graph.size() == 29);
        }
        
        SECTION("the graph should have 29 edges (nodes - 1 + deletion)") {
            REQUIRE(graph.edge_count() == 29);
        }
        
    }

}

TEST_CASE( "Structural duplications and duplicative CNVs can be skipped", "[variantadder]" ) {

    // We'll work on this tiny VCF
    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
ref	4	rs1337	CAAAAAAAAAAAAAAAAAAAAA	CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA,CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	29	PASS	SVTYPE=CNV	GT	2/1	0/0
ref	5	rs1337	AAAAAAAAAAAAAAAAAAAAA	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	29	PASS	SVTYPE=DUP	GT	0/0	0/1
)";

    // Make a stream out of the data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);

    string graph_json = R"({
        "node": [{"id": 1, "sequence": "GCGCAAAAAAAAAAAAAAAAAAAAAGCGC"}],
        "path": [
            {"name": "ref", "mapping": [
                {"position": {"node_id": 1}, "edit": [{"from_length": 29, "to_length": 29}]}
            ]}
        ]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make a VariantAdder
    VariantAdder adder(graph);
    adder.skip_structural_duplications = true;
    // Add the variants to the graph. None should actually make it.
    adder.add_variants(&vcf);

    SECTION("the graph should have 1 node") {
        REQUIRE(graph.size() == 1);
    }
    
    SECTION("the graph should have no edges") {
        REQUIRE(graph.edge_count() == 0);
    }

}

TEST_CASE( "The smart aligner works on very large inserts", "[variantadder]" ) {

    string graph_json = R"({
        "node": [{"id": 1, "sequence": "GCGCAAAAAAAAAAAAAAAAAAAAAGCGC"}]
    })";
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make a VariantAdder
    VariantAdder adder(graph);
    vector<handle_t> order = handlealgs::topological_order(&adder.get_graph());
    NodeSide front(adder.get_graph().get_id(order.front()), false);
    NodeSide back(adder.get_graph().get_id(order.back()), true);
    
    // Make a really long insert
    stringstream s;
    s << "GCGCAAAAAAAAAAA";
    for (size_t i = 0; i < 10000; i++) {
        s << "C";
    }
    s << "AAAAAAAAAAGCGC";
    
    Alignment aligned = adder.smart_align(graph, make_pair(front, back), s.str(), 10000 + graph.length());
    
    SECTION("the resulting alignment should have the input string") {
        REQUIRE(aligned.sequence() == s.str());
    }
    
    SECTION("the resulting alignment should have one mapping") {
        REQUIRE(aligned.path().mapping_size() == 1);
        
        auto& m = aligned.path().mapping(0);
        
        SECTION("that mapping should have 3 edits") {
            REQUIRE(m.edit_size() == 3);
            
            auto& match1 = m.edit(0);
            auto& insert = m.edit(1);
            auto& match2 = m.edit(2);
            
            SECTION("the first edit should be a match of the leading ref part") {
                REQUIRE(edit_is_match(match1));
                REQUIRE(match1.from_length() == strlen("GCGCAAAAAAAAAAA"));
            }
            
            SECTION("the second edit should be an insert of the inserted sequence") {
                REQUIRE(edit_is_insertion(insert));
                REQUIRE(insert.to_length() == 10000);
            }
            
            SECTION("the third edit should be a match of the trailing ref part") {
                REQUIRE(edit_is_match(match2));
                REQUIRE(match2.from_length() == strlen("AAAAAAAAAAGCGC"));
            }
        }
    }

}

TEST_CASE( "The smart aligner should use mapping offsets on huge deletions", "[variantadder]" ) {

    string graph_json = R"({
        "node": [
            {"id": 1, "sequence": "GCGCAAAAAAAAAAAAAAAAAAAA"},
            {"id": 2, "sequence": "<10kAs>"},
            {"id": 3, "sequence": "AAAAAAAAAAAAAAAAAAAAGCGC"}],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 3},
            {"from": 2, "to": 3}
        ]
    })";
    
    // Make the graph have lots of As
    stringstream a_stream;
    for(size_t i = 0; i < 10000; i++) {
        a_stream << "A";
    }
    graph_json = regex_replace(graph_json, std::regex("<10kAs>"), a_stream.str());
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make a VariantAdder
    VariantAdder adder(graph);
    vector<handle_t> order = handlealgs::topological_order(&adder.get_graph());
    NodeSide front(adder.get_graph().get_id(order.front()), false);
    NodeSide back(adder.get_graph().get_id(order.back()), true);
    
    // Make a deleted version (only 21 As)
    string deleted = "GCGCAAAAAAAAAAAAAAAAAAAAAGCGC";
    
    // Align between 1 and 3
    auto endpoints = make_pair(front, back);
    Alignment aligned = adder.smart_align(graph, endpoints, deleted, graph.length());
    
    SECTION("the resulting alignment should have the input string") {
        REQUIRE(aligned.sequence() == deleted);
    }
    
    SECTION("the resulting alignment should have two mappings") {
        REQUIRE(aligned.path().mapping_size() == 2);
        
        auto& m1 = aligned.path().mapping(0);
        auto& m2 = aligned.path().mapping(1);
        
        SECTION("each mapping should be a single edit") {
            REQUIRE(m1.edit_size() == 1);
            REQUIRE(m2.edit_size() == 1);
            
            auto& match1 = m1.edit(0);
            auto& match2 = m2.edit(0);
        
                
            SECTION("the first edit should be a match") {
                REQUIRE(edit_is_match(match1));
            }
            
            SECTION("the second edit should be a match") {
                REQUIRE(edit_is_match(match2));
            }
            
            SECTION("the match lengths should sum to the length of the aligned string") {
                REQUIRE(match1.from_length() + match2.from_length() == deleted.size());
            }
            
            SECTION("the first mapping should be at the start of the first node") {
                REQUIRE(m1.position().node_id() == front.node);
                REQUIRE(m1.position().offset() == 0);
                REQUIRE(m1.position().is_reverse() == false);
            }
            
            SECTION("the second mapping should be at the end of the last node") {
                REQUIRE(m2.position().node_id() == back.node);
                REQUIRE(m2.position().offset() == adder.get_graph().get_length(order.back()) - mapping_from_length(m2));
                REQUIRE(m2.position().is_reverse() == false);
            }
        }
    }

}

TEST_CASE( "The smart aligner should find existing huge deletions", "[variantadder]" ) {

    string graph_json = R"({
        "node": [
            {"id": 1, "sequence": "GCGCAAAAAAAAAAAAAAAAAAAA"},
            {"id": 2, "sequence": "<10kAs>"},
            {"id": 3, "sequence": "AAAAAAAAAAAAAAAAAAAAGCGC"}],
        "edge": [
            {"from": 1, "to": 2},
            {"from": 1, "to": 3},
            {"from": 2, "to": 3}
        ]
    })";
    
    // Make the graph have lots of As
    stringstream a_stream;
    for(size_t i = 0; i < 10000; i++) {
        a_stream << "A";
    }
    graph_json = regex_replace(graph_json, std::regex("<10kAs>"), a_stream.str());
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make a VariantAdder
    VariantAdder adder(graph);
    vector<handle_t> order = handlealgs::topological_order(&adder.get_graph());
    NodeSide front(adder.get_graph().get_id(order.front()), false);
    NodeSide back(adder.get_graph().get_id(order.back()), true);
    
    // Make a deleted version (only 21 As)
    string deleted = "GCGCAAAAAAAAAAAAAAAAAAAAAGCGC";
    
    // Align it between 1 and 3
    auto endpoints = make_pair(front, back);
    Alignment aligned = adder.smart_align(graph, endpoints, deleted, graph.length());
    
    SECTION("the resulting alignment should have the input string") {
        REQUIRE(aligned.sequence() == deleted);
    }
    
    SECTION("the resulting alignment should have two mappings") {
        REQUIRE(aligned.path().mapping_size() == 2);
        
        auto& m1 = aligned.path().mapping(0);
        auto& m2 = aligned.path().mapping(1);
        
        SECTION("each mapping should be a single edit") {
            REQUIRE(m1.edit_size() == 1);
            REQUIRE(m2.edit_size() == 1);
            
            auto& match1 = m1.edit(0);
            auto& match2 = m2.edit(0);
        
                
            SECTION("the first edit should be a match") {
                REQUIRE(edit_is_match(match1));
            }
            
            SECTION("the second edit should be a match") {
                REQUIRE(edit_is_match(match2));
            }
            
            SECTION("the match lengths should sum to the length of the aligned string") {
                REQUIRE(match1.from_length() + match2.from_length() == deleted.size());
            }
            
            SECTION("the first mapping should be at the start of the first node") {
                REQUIRE(m1.position().node_id() == front.node);
                REQUIRE(m1.position().offset() == 0);
                REQUIRE(m1.position().is_reverse() == false);
            }
            
            SECTION("the second mapping should be at the end of the last node") {
                REQUIRE(m2.position().node_id() == back.node);
                REQUIRE(m2.position().offset() == adder.get_graph().get_length(order.back()) - mapping_from_length(m2));
                REQUIRE(m2.position().is_reverse() == false);
            }
        }
    }

}

TEST_CASE( "The smart aligner should use deletion edits on medium deletions", "[variantadder]" ) {

    string graph_json = R"({
        "node": [{"id": 1, "sequence": "GCGC<100As>GCGC"}]
    })";
    
    // Make the graph have lots of As
    stringstream a_stream;
    for(size_t i = 0; i < 100; i++) {
        a_stream << "A";
    }
    graph_json = regex_replace(graph_json, std::regex("<100As>"), a_stream.str());
    
    // Load the JSON
    Graph proto_graph;
    json2pb(proto_graph, graph_json.c_str(), graph_json.size());
    
    // Make it into a VG
    VG graph;
    graph.extend(proto_graph);
    
    // Make a VariantAdder
    VariantAdder adder(graph);
    
    // Make a deleted version (only 21 As)
    string deleted = "GCGCAAAAAAAAAAAAAAAAAAAAAGCGC";
    
    auto endpoints = make_pair(NodeSide(1, false), NodeSide(1, true));
    Alignment aligned = adder.smart_align(graph, endpoints, deleted, graph.length());
    
    SECTION("the resulting alignment should have the input string") {
        REQUIRE(aligned.sequence() == deleted);
    }
    
    SECTION("the resulting alignment should have one or more mappings") {
        
        REQUIRE(aligned.path().mapping_size() >= 1);
        
        auto& mFirst = aligned.path().mapping(0);
        auto& mLast = aligned.path().mapping(aligned.path().mapping_size() - 1);
        
        SECTION("the first and last edit should be a match") {
            REQUIRE(edit_is_match(mFirst.edit(0)));
            REQUIRE(edit_is_match(mLast.edit(mLast.edit_size() - 1)));
        }
        
        int64_t total_matches = 0, total_deleted = 0;
        
        SECTION("the mappings should have leading and trailing match edits and internal deletions") {
            bool in_lead_match = true;
            bool in_middle_deletion = false;
            for (size_t i = 0; i < aligned.path().mapping_size(); i++) {
                const auto& mapping = aligned.path().mapping(i);
                for(size_t j = 0; j < aligned.path().mapping(i).edit_size(); j++) {
                    const auto& edit = mapping.edit(j);
                    if (in_lead_match) {
                        if (edit_is_match(edit)) {
                            total_matches += edit.from_length();
                        }
                        else {
                            REQUIRE(edit_is_deletion(edit));
                            total_deleted += edit.from_length();
                            in_lead_match = false;
                            in_middle_deletion = true;
                        }
                    }
                    else if (in_middle_deletion) {
                        if (edit_is_deletion(edit)) {
                            total_deleted += edit.from_length();
                        }
                        else {
                            REQUIRE(edit_is_match(edit));
                            total_matches += edit.from_length();
                            in_middle_deletion = false;
                        }
                    }
                    else {
                        REQUIRE(edit_is_match(edit));
                        total_matches += edit.from_length();
                    }
                }
            }
            REQUIRE(total_deleted == 100 - 21);
            REQUIRE(total_matches == deleted.size());
        }
    }

}


}
}
