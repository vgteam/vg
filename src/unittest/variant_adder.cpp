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

}
}
