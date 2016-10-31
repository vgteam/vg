/**
 * unittest/constructor.cpp: test cases for the vg graph constructor tool
 */

#include "catch.hpp"
#include "../constructor.hpp"

#include "../path.hpp"
#include "../json2pb.h"

#include <vector>
#include <sstream>
#include <iostream>

namespace vg {
namespace unittest {

TEST_CASE( "An empty chunk with no variants can be constructed", "[constructor]" ) {
    Constructor constructor;
    
    auto result = constructor.construct_chunk("", "empty", std::vector<vcflib::Variant>());
    
    SECTION("the graph should have no elements") {
        REQUIRE(result.graph.node_size() == 0);
        REQUIRE(result.graph.edge_size() == 0);
        REQUIRE(result.left_ends.empty());
        REQUIRE(result.right_ends.empty());
    }

}

TEST_CASE( "A small linear chunk with no variants can be constructed", "[constructor]" ) {
    Constructor constructor;
    
    auto result = constructor.construct_chunk("GATTACA", "movie", std::vector<vcflib::Variant>());
    
    SECTION("the graph should have one node") {
        REQUIRE(result.graph.node_size() == 1);
        auto& node = result.graph.node(0);
        
        SECTION("the node should have the full sequence") {
            REQUIRE(node.sequence() == "GATTACA");
        }
        
        SECTION("the node should have the correct ID") {
            REQUIRE(node.id() == 1);
        }
        
        SECTION("the node should be the only exposed node on the left") {
            REQUIRE(result.left_ends.count(node.id()));
            REQUIRE(result.left_ends.size() == 1);
        }
        
        SECTION("the node should be the only exposed node on the right") {
            REQUIRE(result.right_ends.count(node.id()));
            REQUIRE(result.right_ends.size() == 1);
        }
    }
    
    SECTION("the graph should have no edges") {
        REQUIRE(result.graph.edge_size() == 0);
    }
    
    SECTION("the graph should have one path") {
        REQUIRE(result.graph.path_size() == 1);
        
        auto& path = result.graph.path(0);
        
        SECTION("the path should have the name we passed in") {
            REQUIRE(path.name() == "movie");
        }
        
        SECTION("the path should have one mapping") {
            REQUIRE(path.mapping_size() == 1);
            
            auto& mapping = path.mapping(0);
            
            SECTION("the mapping should be a full length perfect match") {
                REQUIRE(mapping_is_match(mapping));
                REQUIRE(from_length(mapping) == 7);
            }
            
            SECTION("the mapping should be on the node") {
                REQUIRE(mapping.position().node_id() == 1);
                REQUIRE(mapping.position().offset() == 0);
                REQUIRE(mapping.position().is_reverse() == false);
            }
        }
    }
}

TEST_CASE( "Max node length is respected", "[constructor]" ) {
    Constructor constructor;
    constructor.max_node_size = 4;
    
    auto result = constructor.construct_chunk("GATTACA", "movie", std::vector<vcflib::Variant>());
    
    SECTION("the graph should have two nodes") {
        REQUIRE(result.graph.node_size() == 2);
        auto& node1 = result.graph.node(0);
        auto& node2 = result.graph.node(1);
        
        SECTION("node 1 should have the first part of the sequence") {
            REQUIRE(node1.sequence() == "GATT");
            REQUIRE(node1.id() == 1);
        }
        
        SECTION("node 2 should have the second part of the sequence") {
            REQUIRE(node2.sequence() == "ACA");
            REQUIRE(node2.id() == 2);
        }
        
        SECTION("node 1 should be exposed on the left") {
            REQUIRE(result.left_ends.count(node1.id()));
            REQUIRE(result.left_ends.size() == 1);
        }
        
        SECTION("node 2 should be exposed on the right") {
            REQUIRE(result.right_ends.count(node2.id()));
            REQUIRE(result.right_ends.size() == 1);
        }
    }
    
    SECTION("the graph should have one edge") {
        REQUIRE(result.graph.edge_size() == 1);
        
        auto& edge = result.graph.edge(0);
        
        SECTION("the edge should connect node 1 to node 2") {
            REQUIRE(edge.from() == 1);
            REQUIRE(edge.to() == 2);
            REQUIRE(edge.from_start() == false);
            REQUIRE(edge.to_end() == false);
        }
    }
    
    SECTION("the graph should have one path") {
        REQUIRE(result.graph.path_size() == 1);
        
        auto& path = result.graph.path(0);
        
        SECTION("the path should have the name we passed in") {
            REQUIRE(path.name() == "movie");
        }
        
        SECTION("the path should have two mappings") {
            REQUIRE(path.mapping_size() == 2);
            
            auto& mapping1 = path.mapping(0);
            auto& mapping2 = path.mapping(1);
            
            SECTION("mapping 1 should be a full length perfect match") {
                REQUIRE(mapping_is_match(mapping1));
                REQUIRE(from_length(mapping1) == 4);
            }
            
            SECTION("mapping 2 should be a full length perfect match") {
                REQUIRE(mapping_is_match(mapping2));
                REQUIRE(from_length(mapping2) == 3);
            }
            
            SECTION("mapping 1 should be on node 1") {
                REQUIRE(mapping1.position().node_id() == 1);
                REQUIRE(mapping1.position().offset() == 0);
                REQUIRE(mapping1.position().is_reverse() == false);
            }
            
            SECTION("mapping 2 should be on node 2") {
                REQUIRE(mapping2.position().node_id() == 2);
                REQUIRE(mapping2.position().offset() == 0);
                REQUIRE(mapping2.position().is_reverse() == false);
            }
        }
    }
}

TEST_CASE( "A SNP can be constructed", "[constructor]" ) {

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

    std::stringstream vcf_stream(vcf_data);

    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);
    
    std::vector<vcflib::Variant> variants;
    vcflib::Variant var;
    while (vcf.getNextVariant(var)) {
        // Load up the VCF
        // Make sure to correct it to 0-based
        var.position -= 1;
        variants.push_back(var);
    }

    auto ref = "GATTACA";
    
    Constructor constructor;

    // Construct the graph    
    auto result = constructor.construct_chunk(ref, "ref", variants);


    std::cerr << pb2json(result.graph) << std::endl;

    SECTION("the graph should have 4 nodes") {
        REQUIRE(result.graph.node_size() == 4);
        
        // Find all the nodes
        Node before;
        Node after;
        Node snp_ref;
        Node snp_alt;
        
        for (size_t i = 0; i < result.graph.node_size(); i++) {
            auto& node = result.graph.node(i);
            
            if (node.sequence() == "GATT") {
                before = node;
            } else if (node.sequence() == "CA") {
                after = node;
            } else if (node.sequence() == "A") {
                snp_ref = node;
            } else if (node.sequence() == "G") {
                snp_alt = node;
            }
        }
        
        SECTION("before, after, ref, and alt nodes should be present") {
            REQUIRE(before.id() != 0);
            REQUIRE(after.id() != 0);
            REQUIRE(snp_ref.id() != 0);
            REQUIRE(snp_alt.id() != 0);
        }
        
    }
	

}

TEST_CASE( "A small VCF can be constructed", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
20	14370	rs6054257	G	A	29	PASS	.	GT
20	17330	.	T	A	3	q10	.	GT
20	1110696	rs6040355	A	G,T	67	PASS	.	GT
20	1230237	.	T	.	47	PASS	.	GT
20	1234567	microsat1	GTCT	G,GTACT	50	PASS	.	GT
)";

	

}

// TODO: graphs with variants

}
}
