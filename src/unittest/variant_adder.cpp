/**
 * \file 
 * unittest/variant_adder.cpp: test cases for the add-variants-to-existing-graph tool
 */

#include "catch.hpp"
#include "../variant_adder.hpp"

#include "../utility.hpp"
#include "../path.hpp"
#include "../json2pb.h"

#include <vector>
#include <sstream>
#include <iostream>

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
    
        graph.dice_nodes(1);
    
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

TEST_CASE( "The smart aligner works on very large inserts", "[variantadder]" ) {

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
    
    // Make a really long insert
    stringstream s;
    s << "GCGCAAAAAAAAAAA";
    for (size_t i = 0; i < 10000; i++) {
        s << "C";
    }
    s << "AAAAAAAAAAGCGC";
    
    Alignment aligned = adder.smart_align(graph, s.str(), 10000 + graph.length());
    
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

}
}
