/**
 * \file 
 * unittest/constructor.cpp: test cases for the vg graph constructor tool
 */

#include "catch.hpp"
#include "../constructor.hpp"

#include "../utility.hpp"
#include "../path.hpp"
#include "vg/io/json2pb.h"

#include <vector>
#include <sstream>
#include <iostream>
#include <unordered_map>

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {

TEST_CASE( "An empty chunk with no variants can be constructed", "[constructor]" ) {
    Constructor constructor;
    
    auto result = constructor.construct_chunk("", "empty", std::vector<vcflib::Variant>(), 0);
    
    SECTION("the graph should have no elements") {
        REQUIRE(result.graph.node_size() == 0);
        REQUIRE(result.graph.edge_size() == 0);
        REQUIRE(result.left_ends.empty());
        REQUIRE(result.right_ends.empty());
    }

}

TEST_CASE( "A small linear chunk with no variants can be constructed", "[constructor]" ) {
    Constructor constructor;
    
    auto result = constructor.construct_chunk("GATTACA", "movie", std::vector<vcflib::Variant>(), 0);
    
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

TEST_CASE( "A chunk with lowercase letters has them uppercased", "[constructor]" ) {
    Constructor constructor;
    // Don't warn during testing
    constructor.warn_on_lowercase = false;
    
    auto result = constructor.construct_chunk("cangantan", "lower", std::vector<vcflib::Variant>(), 0);
    
    SECTION("the graph should have one node") {
        REQUIRE(result.graph.node_size() == 1);
        auto& node = result.graph.node(0);
        
        SECTION("the node should have the uppercase sequence") {
            REQUIRE(node.sequence() == "CANGANTAN");
        }
    }
    
}

TEST_CASE( "Max node length is respected", "[constructor]" ) {
    Constructor constructor;
    constructor.max_node_size = 4;
    
    // The reasoning here applies to greedy mode, not old-vg-construct-mimicing mode
    constructor.greedy_pieces = true;
    
    auto result = constructor.construct_chunk("GATTACA", "movie", std::vector<vcflib::Variant>(), 0);
    
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

/**
 * Testing wrapper to build a graph chunk from a VCF string. Adds alt paths by default.
 */
static ConstructedChunk construct_test_chunk(string ref_sequence, string ref_name, string vcf_data) {
    
    // Make a stream out of the data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);
    
    // Fill in this vector of variants
    std::vector<vcflib::Variant> variants;
    vcflib::Variant var;
    while (vcf.getNextVariant(var)) {
        // Make sure to correct each variant's position to 0-based
        //var.position -= 1;
        variants.push_back(var);
    }

    Constructor constructor;
    constructor.alt_paths = true;
    // Make sure we can test the node splitting behavior at reasonable sizes
    constructor.max_node_size = 50;

    // Construct the graph    
    return constructor.construct_chunk(ref_sequence, ref_name, variants, 0);
}

/**
 * Testing wrapper to build a whole graph from a VCF string. Adds alt paths by default.
 */
static Graph construct_test_graph(string fasta_data, string vcf_data, size_t max_node_size,
    bool do_svs, bool use_flat_alts = false) {
    
    // Merge all the graphs we get into this graph
    Graph built;
    
    // Make a stream out of the VCF data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);
    
    // Put it in a vector
    vector<vcflib::VariantCallFile*> vcf_pointers {&vcf};
    
    // We have to write the FASTA to a file
    string fasta_filename = temp_file::create();
    ofstream fasta_stream(fasta_filename);
    fasta_stream << fasta_data;
    fasta_stream.close(); 
    
    // Make a FastaReference out of it
    FastaReference reference;
    reference.open(fasta_filename);
    
    // Put it in a vector
    vector<FastaReference*> fasta_pointers {&reference};
    
    // Make an empty vector of insertion files
    vector<FastaReference*> ins_pointers;
    
    // Make a callback to handle the output
    auto callback = [&](Graph& constructed) {
        // Merge everything that comes out into one graph in memory.
        #pragma omp critical
        built.MergeFrom(constructed);
    };
    
    Constructor constructor;
    constructor.alt_paths = true;
    constructor.do_svs = do_svs;
    constructor.flat = use_flat_alts;
    // Make sure we can test the node splitting behavior at reasonable sizes
    constructor.max_node_size = max_node_size;

    // Construct the graph    
    constructor.construct_graph(fasta_pointers, vcf_pointers, ins_pointers, callback);
    
    // Delete our temporary file
    temp_file::remove(fasta_filename);
    
    // Return the aggregated result
    return built;
}

/**
 * Testing wrapper to build a whole graph from a VCF string into a HandleGraph
 */
static unique_ptr<PathHandleGraph> construct_test_handle_graph(string fasta_data, string vcf_data, 
    size_t max_node_size, bool do_svs, bool use_flat_alts = false) {
  
    // Make an empty graph
    auto build_to = new bdsg::HashGraph();
    // Make a unique_ptr own our graph
    unique_ptr<PathHandleGraph> graph(build_to);
    
    // Make a stream out of the VCF data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);
    
    // Put it in a vector
    vector<vcflib::VariantCallFile*> vcf_pointers {&vcf};
    
    // We have to write the FASTA to a file
    string fasta_filename = temp_file::create();
    ofstream fasta_stream(fasta_filename);
    fasta_stream << fasta_data;
    fasta_stream.close(); 
    
    // Make a FastaReference out of it
    FastaReference reference;
    reference.open(fasta_filename);
    
    // Put it in a vector
    vector<FastaReference*> fasta_pointers {&reference};
    
    // Make an empty vector of insertion files
    vector<FastaReference*> ins_pointers;
    
    Constructor constructor;
    constructor.alt_paths = true;
    constructor.do_svs = do_svs;
    constructor.flat = use_flat_alts;
    // Make sure we can test the node splitting behavior at reasonable sizes
    constructor.max_node_size = max_node_size;

    // Construct the graph    
    constructor.construct_graph(fasta_pointers, vcf_pointers, ins_pointers, build_to);
    
    // Delete our temporary file
    temp_file::remove(fasta_filename);
    
    // Return the result
    return std::move(graph);
}

// Find and copy path in the graph that ends with the given name suffix, or
// produce an empty Path.
static Path find_path_by_suffix(const Graph& g, const std::string& suffix) {
    for (size_t i = 0; i < g.path_size(); i++) {
        auto& path = g.path(i);
        if (path.name().size() >= suffix.size()) {
            // This could be it
            if (std::equal(suffix.rbegin(), suffix.rend(), path.name().rbegin())) {
                // Path does end with suffix
                return path;
            }
        }
    }
    // If we get here no path matches
    return Path();
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

    auto ref = "GATTACA";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);

#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif

    // The graph will have these 4 nodes: before and after the SNP, and the two
    // alts of the SNP.
    Node before;
    Node after;
    Node snp_ref;
    Node snp_alt;
    
    // Find all the nodes. All the other test cases depend on knowing them.
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

    SECTION("the graph should have 4 nodes") {
        REQUIRE(result.graph.node_size() == 4);
        
        SECTION("before, after, ref, and alt nodes should be present") {
            REQUIRE(before.id() != 0);
            REQUIRE(after.id() != 0);
            REQUIRE(snp_ref.id() != 0);
            REQUIRE(snp_alt.id() != 0);
        }
        
        SECTION("the single source should be the very first node, with ID 1") {
            REQUIRE(before.id() == 1);
            REQUIRE(result.left_ends.size() == 1);
            REQUIRE(result.left_ends.count(before.id()) == 1);
            REQUIRE(result.graph.node(0).id() == before.id());
        }
        
        SECTION("the single sink should be the very last node, with ID max_id") {
            REQUIRE(after.id() == result.max_id);
            REQUIRE(result.right_ends.size() == 1);
            REQUIRE(result.right_ends.count(after.id()) == 1);
            REQUIRE(result.graph.node(result.graph.node_size() - 1).id() == after.id());
        }
        
    }
    
    SECTION("the graph should have 4 edges") {
        REQUIRE(result.graph.edge_size() == 4);
        
        // We want to recognize all the edges based on the nodes.
        Edge to_ref;
        Edge from_ref;
        Edge to_alt;
        Edge from_alt;
        
        for (size_t i = 0; i < result.graph.edge_size(); i++) {
            auto& edge = result.graph.edge(i);
            
            // Match each edge we expect
            if (edge.from() == before.id() && edge.to() == snp_ref.id()) {
                to_ref = edge;
            } else if (edge.from() == before.id() && edge.to() == snp_alt.id()) {
                to_alt = edge;
            } else if (edge.from() == snp_ref.id() && edge.to() == after.id()) {
                from_ref = edge;
            } else if (edge.from() == snp_alt.id() && edge.to() == after.id()) {
                from_alt = edge;
            }
        }
        
        SECTION("edges should connect into and out of both ref and alt alleles") {
            // Now check them
            REQUIRE(to_ref.from() == before.id());
            REQUIRE(to_ref.to() == snp_ref.id());
            REQUIRE(!to_ref.from_start());
            REQUIRE(!to_ref.to_end());
            
            REQUIRE(to_alt.from() == before.id());
            REQUIRE(to_alt.to() == snp_alt.id());
            REQUIRE(!to_alt.from_start());
            REQUIRE(!to_alt.to_end());
            
            REQUIRE(from_ref.from() == snp_ref.id());
            REQUIRE(from_ref.to() == after.id());
            REQUIRE(!from_ref.from_start());
            REQUIRE(!from_ref.to_end());
            
            REQUIRE(from_alt.from() == snp_alt.id());
            REQUIRE(from_alt.to() == after.id());
            REQUIRE(!from_alt.from_start());
            REQUIRE(!from_alt.to_end());
        }
        
        
    }
    
    // Find the primary path, and the paths for the two alleles, filling in
    // empty paths if absent.
    Path primary = find_path_by_suffix(result.graph, "ref");
    Path allele0 = find_path_by_suffix(result.graph, "0");
    Path allele1 = find_path_by_suffix(result.graph, "1");
    
    SECTION("only one allele should be allowed to be empty") {
        REQUIRE(primary.mapping_size() > 0);
        REQUIRE(allele0.mapping_size() + allele1.mapping_size() > 0);
    }
        
    SECTION("primary, ref allele, and alt allele paths should be named correctly") {
        REQUIRE(primary.name() == "ref");
        REQUIRE(allele0.name().substr(0, 5) == "_alt_");
        REQUIRE(allele0.name().substr(allele0.name().size() - 2, 2) == "_0");
        REQUIRE(allele1.name().substr(0, 5) == "_alt_");
        REQUIRE(allele1.name().substr(allele1.name().size() - 2, 2) == "_1");
        
        // And the two alleles have to be of the same variant
        REQUIRE(allele0.name().substr(5, allele0.name().size() - (5 + 2)) == 
            allele1.name().substr(5, allele1.name().size() - (5 + 2)));
            
        SECTION("the primary path should trace the reference") {
            REQUIRE(primary.mapping_size() == 3);
            
            REQUIRE(primary.mapping(0).position().node_id() == before.id());
            REQUIRE(primary.mapping(0).position().offset() == 0);
            REQUIRE(primary.mapping(0).position().is_reverse() == false);
            REQUIRE(mapping_is_match(primary.mapping(0)));
            REQUIRE(from_length(primary.mapping(0)) == before.sequence().size());
            
            REQUIRE(primary.mapping(1).position().node_id() == snp_ref.id());
            REQUIRE(primary.mapping(1).position().offset() == 0);
            REQUIRE(primary.mapping(1).position().is_reverse() == false);
            REQUIRE(mapping_is_match(primary.mapping(1)));
            REQUIRE(from_length(primary.mapping(1)) == snp_ref.sequence().size());
            
            REQUIRE(primary.mapping(2).position().node_id() == after.id());
            REQUIRE(primary.mapping(2).position().offset() == 0);
            REQUIRE(primary.mapping(2).position().is_reverse() == false);
            REQUIRE(mapping_is_match(primary.mapping(2)));
            REQUIRE(from_length(primary.mapping(2)) == after.sequence().size());
            
        }
        
        SECTION("the ref allele path should visit the ref allele") {
            REQUIRE(allele0.mapping_size() == 1);
            
            REQUIRE(allele0.mapping(0).position().node_id() == snp_ref.id());
            REQUIRE(allele0.mapping(0).position().offset() == 0);
            REQUIRE(allele0.mapping(0).position().is_reverse() == false);
            REQUIRE(mapping_is_match(allele0.mapping(0)));
            REQUIRE(from_length(allele0.mapping(0)) == snp_ref.sequence().size());
        }
        
        SECTION("the alt allele path should visit the alt allele") {
            REQUIRE(allele1.mapping_size() == 1);
            
            REQUIRE(allele1.mapping(0).position().node_id() == snp_alt.id());
            REQUIRE(allele1.mapping(0).position().offset() == 0);
            REQUIRE(allele1.mapping(0).position().is_reverse() == false);
            REQUIRE(mapping_is_match(allele0.mapping(0)));
            REQUIRE(from_length(allele0.mapping(0)) == snp_alt.sequence().size());
        }
                
        SECTION("the reference path should be path 0") {
            REQUIRE(result.graph.path(0).name() == primary.name());
        }
    }

}

TEST_CASE( "A deletion can be constructed", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	5	rs1337	AC	A	29	PASS	.	GT
)";

    auto ref = "GATTACA";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    // We could build this either GATT,A,C,A or GATTA,C,A.
    // In either case we have as many nodes as edges.
    
    SECTION("the graph should have either 3 or 4 nodes depending on structure") {
        REQUIRE(result.graph.node_size() >= 3);
        REQUIRE(result.graph.node_size() <= 4);
    }
    
    SECTION("the graph should have as many edges as nodes") {
        REQUIRE(result.graph.edge_size() == result.graph.node_size());
    }

    // Find the primary path, and the paths for the two alleles, filling in
    // empty paths if absent.
    Path primary = find_path_by_suffix(result.graph, "ref");
    Path allele0 = find_path_by_suffix(result.graph, "0");
    Path allele1 = find_path_by_suffix(result.graph, "1");
    
    SECTION("only one allele should be allowed to be empty") {
        REQUIRE(primary.mapping_size() > 0);
        REQUIRE(allele0.mapping_size() + allele1.mapping_size() > 0);
    }
    
    SECTION("the path for the alt should not include the deleted sequence") {
        if (allele1.mapping_size() == 0) {
            // This definitely lacks the C
            REQUIRE(true);
        } else {
            for (size_t i = 0; i < allele1.mapping_size(); i++) {
                // Look at all the nodes along the path
                id_t node_id = allele1.mapping(i).position().node_id();
                
                for (size_t j = 0; j < result.graph.node_size(); j++) {
                    // Brute force the whole graph to find the node
                    if(node_id == result.graph.node(j).id()) {
                        // In the node we actually visit, there can't be a "C", since we deleted them all
                        REQUIRE(result.graph.node(j).sequence().find("C") == string::npos);
                    }
                }
            }
        }
    }

}

TEST_CASE( "An insertion can be constructed", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	5	rs1337	A	AC	29	PASS	.	GT
)";

    auto ref = "GATTAA";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    // We could build this either GATT,A,C,A or GATTA,C,A.
    // In either case we have as many nodes as edges.
    
    SECTION("the graph should have either 3 or 4 nodes depending on structure") {
        REQUIRE(result.graph.node_size() >= 3);
        REQUIRE(result.graph.node_size() <= 4);
    }
    
    SECTION("the graph should have as many edges as nodes") {
        REQUIRE(result.graph.edge_size() == result.graph.node_size());
    }

    // Find the primary path, and the paths for the two alleles, filling in
    // empty paths if absent.
    Path primary = find_path_by_suffix(result.graph, "ref");
    Path allele0 = find_path_by_suffix(result.graph, "0");
    Path allele1 = find_path_by_suffix(result.graph, "1");
    
    SECTION("only one allele should be allowed to be empty") {
        REQUIRE(primary.mapping_size() > 0);
        REQUIRE(allele0.mapping_size() + allele1.mapping_size() > 0);
    }
        
    SECTION("the path for the ref should not include the inserted sequence") {
        if (allele0.mapping_size() == 0) {
            // This definitely lacks the C
            REQUIRE(true);
        } else {
            for (size_t i = 0; i < allele0.mapping_size(); i++) {
                // Look at all the nodes along the path
                id_t node_id = allele0.mapping(i).position().node_id();
                
                for (size_t j = 0; j < result.graph.node_size(); j++) {
                    // Brute force the whole graph to find the node
                    if(node_id == result.graph.node(j).id()) {
                        // In the node we actually visit, there can't be a "C", since we inserted the only one
                        REQUIRE(result.graph.node(j).sequence().find("C") == string::npos);
                    }
                }
            }
        }
    }

}

TEST_CASE( "A SNP nested inside a deletion can be constructed", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	2	.	ATGTTCTTCC	A	100	PASS	.	GT
ref	6	.	T	C	100	PASS	.	GT
)";

    auto ref = "GATGTTCTTCCG";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
    // It should be like
    //
    //       /C\
    // GA TGT T CTTCC G
    //   \-----------/
    
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    SECTION("the graph should have 6 nodes") {
        REQUIRE(result.graph.node_size() == 6);
    }
    
    SECTION("the graph should have 7 edges") {
        REQUIRE(result.graph.edge_size() == 7);
    }
}

TEST_CASE( "A variable count repeat can be constructed", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	2	.	GAATT	GAATTAATT,GAATTAATTAATT,G	100	PASS	.	GT
)";

    auto ref = "CGAATTC";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
    // It should be like
    //    /AATTAATT\
    //   /-AATT-----\
    // CG------------AATT-C 
    //   \---------------/
    
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    SECTION("the graph should have 5 nodes") {
        REQUIRE(result.graph.node_size() == 5);
    }
    
    SECTION("the graph should have 7 edges") {
        REQUIRE(result.graph.edge_size() == 7);
    }
}

TEST_CASE( "A merged SNP and indel can be constructed", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	2	rs5839893	TG	GG,T	100	PASS	.	GT
)";

    auto ref = "ATGA";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    
    
    SECTION("the graph should have 5 nodes") {
        REQUIRE(result.graph.node_size() == 5);
    }
    
    SECTION("the graph should have 7 edges") {
        REQUIRE(result.graph.edge_size() == 7);
    }
}


TEST_CASE( "Path names do not depend on chunking", "[constructor]" ) {

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

    string ref = "GATTACA";

    // Make a stream out of the data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);
    
    // Fill in this vector of variants
    std::vector<vcflib::Variant> variants;
    vcflib::Variant var;
    while (vcf.getNextVariant(var)) {
        // Make sure to correct each variant's position to 0-based
        //var.position -= 1;
        variants.push_back(var);
    }

    // Make a constructor
    Constructor constructor;
    constructor.alt_paths = true;
    
    // Construct the graph    
    auto result1 = constructor.construct_chunk(ref, "ref", variants, 0);
    
    // Construct the graph with a slight offset
    auto result2 = constructor.construct_chunk(ref.substr(1), "ref", variants, 1);

    // Get the two sets of path names
    set<string> paths1;
    for(size_t i = 0; i < result1.graph.path_size(); i++) {
        paths1.insert(result1.graph.path(i).name());
    }
    set<string> paths2;
    for(size_t i = 0; i < result2.graph.path_size(); i++) {
        paths2.insert(result2.graph.path(i).name());
    }
    
    REQUIRE(paths1 == paths2);

}


TEST_CASE( "Outer matching sequence is trimmed on inserts", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	3	rs1337	TTC	TTAC	29	PASS	.	GT
)";

    auto ref = "GATTCA";
    
    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    // We insist on building this GATT,A,CA with the minimum insert
    
    SECTION("the graph should have the minimum number of nodes") {
        REQUIRE(result.graph.node_size() == 3);
        
        SECTION("the nodes should be pre-insert, inserted sequence, and post-insert") {
            CHECK(result.graph.node(0).sequence() == "GATT");
            CHECK(result.graph.node(1).sequence() == "A");
            CHECK(result.graph.node(2).sequence() == "CA");
        }
        
        SECTION("the nodes should be numbered 1, 2, 3 in order") {
            CHECK(result.graph.node(0).id() == 1);
            CHECK(result.graph.node(1).id() == 2);
            CHECK(result.graph.node(2).id() == 3);
        }
    }
    
    SECTION("the graph should have the minimum number of edges") {
        REQUIRE(result.graph.edge_size() == 3);
    }

    // Find the primary path, and the paths for the two alleles, filling in
    // empty paths if absent.
    Path primary = find_path_by_suffix(result.graph, "ref");
    Path allele0 = find_path_by_suffix(result.graph, "0");
    Path allele1 = find_path_by_suffix(result.graph, "1");
    
    SECTION("only one allele should be allowed to be empty") {
        REQUIRE(primary.mapping_size() > 0);
        REQUIRE(allele0.mapping_size() + allele1.mapping_size() > 0);
    }
    
    SECTION("the path for the ref allele should be completely empty") {
        CHECK(allele0.mapping_size() == 0);
    }
    
    SECTION("the path for the alt allele should have just one node") {
        CHECK(allele1.mapping_size() == 1);
    }

}

TEST_CASE( "Large deletions are broken appropriately", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	2	rs1337	CATTTAATTATTAATTAAATAATTTAATTTATTTATTATTATAAATTTATTAATATAAATTAAATA	C	29	PASS	.	GT
)";

    auto ref = "GCATTTAATTATTAATTAAATAATTTAATTTATTTATTATTATAAATTTATTAATATAAATTAAATAG";
    
    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    // We insist on building this GATT,A,CA with the minimum insert
    
    SECTION("the graph should have 4 nodes") {
        REQUIRE(result.graph.node_size() == 4);
        
        SECTION("the nodes should be pre-deletion, deleted sequence 1, deleted sequence 2, and post-deletion") {
            CHECK(result.graph.node(0).sequence() == "GC");
            CHECK(result.graph.node(1).sequence() == "ATTTAATTATTAATTAAATAATTTAATTTATTTATTATTATAAATTTATT");
            CHECK(result.graph.node(2).sequence() == "AATATAAATTAAATA");
            CHECK(result.graph.node(3).sequence() == "G");
        }
        
        SECTION("the nodes should be numbered 1, 2, 3, and 4, in order") {
            CHECK(result.graph.node(0).id() == 1);
            CHECK(result.graph.node(1).id() == 2);
            CHECK(result.graph.node(2).id() == 3);
            CHECK(result.graph.node(3).id() == 4);
        }
    }
    
    SECTION("the graph should have 4 edges") {
        REQUIRE(result.graph.edge_size() == 4);
    }

    // Find the primary path, and the paths for the two alleles, filling in
    // empty paths if absent.
    Path primary = find_path_by_suffix(result.graph, "ref");
    Path allele0 = find_path_by_suffix(result.graph, "0");
    Path allele1 = find_path_by_suffix(result.graph, "1");
    
    SECTION("only one allele should be allowed to be empty") {
        REQUIRE(primary.mapping_size() > 0);
        REQUIRE(allele0.mapping_size() + allele1.mapping_size() > 0);
    }
        
    SECTION("the path for the ref allele should have 2 nodes") {
        CHECK(allele0.mapping_size() == 2);
    }
    
    SECTION("the path for the alt allele should be completely empty") {
        CHECK(allele1.mapping_size() == 0);
    }

}


TEST_CASE( "Multiple inserts don't cross-link", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	6	rs10666768	C	CG,CGG,G	100	PASS	.	GT
)";

    auto ref = "GATTACA";
    
    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    // We insist on building this as:
    //       +-C-+-G--+ 
    // GATTA-+   +----+-A
    //       +-G-+-GG-+
    
    SECTION("the graph should have the minimum number of nodes") {
        REQUIRE(result.graph.node_size() == 6);
    }
    
    SECTION("the graph should contain no self loops") {
        for(size_t i = 0; i < result.graph.edge_size(); i++) {
            auto& edge = result.graph.edge(i);
            REQUIRE(edge.from() != edge.to());
        }
    }
    
    SECTION("the graph should have all and only the edges between the 2-way SNP and the 3-way indel") {
        REQUIRE(result.graph.edge_size() == 10);
    }
    
    

}

TEST_CASE( "A combination insertion and deletion gets appropriate alt paths", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	2	rs554978012;rs201121256	GTA	GTATA,G	.	GT
)";

    auto ref = "CGTATACC";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    // We build this as CG,TA,TA,TAC, with both TAs being taken in the longest
    // allele, neither in the shortestallele, and the rightmost one in the
    // reference allele.
    
    SECTION("the graph should have 4 nodes") {
        REQUIRE(result.graph.node_size() == 4);
        
        SECTION("the nodes should be pre-variant, inserted TA, deleted TA, and post-variant") {
            CHECK(result.graph.node(0).sequence() == "CG");
            CHECK(result.graph.node(1).sequence() == "TA");
            CHECK(result.graph.node(2).sequence() == "TA");
            CHECK(result.graph.node(3).sequence() == "TACC");
        }
        
        SECTION("the nodes should be numbered 1, 2, 3, and 4, in order") {
            CHECK(result.graph.node(0).id() == 1);
            CHECK(result.graph.node(1).id() == 2);
            CHECK(result.graph.node(2).id() == 3);
            CHECK(result.graph.node(3).id() == 4);
        }
    }
    
    SECTION("the graph should have 5 edges") {
        REQUIRE(result.graph.edge_size() == 5);
    }

    // Find the primary path, and the paths for the two alleles, filling in
    // empty paths if absent.
    Path primary = find_path_by_suffix(result.graph, "ref");
    Path allele0 = find_path_by_suffix(result.graph, "0");
    Path allele1 = find_path_by_suffix(result.graph, "1");
    Path allele2 = find_path_by_suffix(result.graph, "2");
    
    SECTION("only one allele should be allowed to be empty") {
        REQUIRE(primary.mapping_size() > 0);
        REQUIRE(allele0.mapping_size() + allele1.mapping_size() > 0);
        REQUIRE(allele1.mapping_size() + allele2.mapping_size() > 0);
        REQUIRE(allele0.mapping_size() + allele2.mapping_size() > 0);
    }
        
    SECTION("the path for the reference alt should visit the second TA node") {
        REQUIRE(allele0.mapping_size() == 1);
        REQUIRE(allele0.mapping(0).position().node_id() == 3);
    }
    
    SECTION("the path for the insert alt should visit the second first and second TA nodes") {
        REQUIRE(allele1.mapping_size() == 2);
        REQUIRE(allele1.mapping(0).position().node_id() == 2);
        REQUIRE(allele1.mapping(1).position().node_id() == 3);
    }
    
    SECTION("the path for the delete alt should be empty") {
        REQUIRE(allele2.mapping_size() == 0);
    }

}

TEST_CASE( "An insert with adjacent SNP can be constructed", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	2	rs572383716	T	TA	100	PASS	.	GT    
ref	3	rs76837267	A	T	100	PASS	.	GT
)";

    auto ref = "CTATAC";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
    // It should be like
    // CT,A/-,A/T,AC
    
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif
    
    SECTION("the graph should have 5 nodes") {
        REQUIRE(result.graph.node_size() == 5);
    }
    
    SECTION("the graph should have 7 edges") {
        REQUIRE(result.graph.edge_size() == 7);
    }
}


TEST_CASE( "A VCF with multiple clumps can be constructed", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	1	.	GA	A	29	PASS	.	GT
ref	5	rs1337	AC	A	29	PASS	.	GT
ref	5	.	A	T	29	PASS	.	GT
ref	6	rs1338	C	G	29	PASS	.	GT
ref	11	.	TAG	T	29	PASS	.	GT
)";

    auto ref = "GATTACACATTAG";

    // Build the graph
    auto result = construct_test_chunk(ref, "ref", vcf_data);
    
#ifdef debug
    std::cerr << pb2json(result.graph) << std::endl;
#endif

    SECTION("a leading deletion is recognized") {
        REQUIRE(result.left_ends.size() == 2);
    }
    
    SECTION("a trailing deletion is recognized") {
        REQUIRE(result.right_ends.size() == 2);
    }
    
    // TODO: the center ought to look like this:
    // (A)TT A-+->C-+->ACAT(T)
    //       T-+----/

}

TEST_CASE( "A graph can be constructed to a HandleGraph", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	5	.	A	T	29	PASS	.	GT
)";

    auto fasta_data = R"(>ref
GATTACACATTAG
)";

    // Build the graph
    unique_ptr<PathHandleGraph> result = construct_test_handle_graph(fasta_data, vcf_data, 50, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    SECTION("graph is the correct size") {
        REQUIRE(result->get_node_count() == 4);
    }

}


TEST_CASE( "A VCF and FASTA on two contigs make a graph with a consistent ID space", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref1	1	.	GA	A	29	PASS	.	GT
ref1	5	rs1337	AC	A	29	PASS	.	GT
ref2	5	.	A	T	29	PASS	.	GT
ref2	6	rs1338	C	G	29	PASS	.	GT
ref2	11	.	TAG	T	29	PASS	.	GT
)";

    auto fasta_data = R"(>ref1
GATTACACATTAG
>ref2
GATTACACATTAG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 50, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    SECTION("node IDs are not repeated") {
        set<id_t> seen_ids;
        
        for (size_t i = 0; i < result.node_size(); i++) {
            // Look at each node
            auto& node = result.node(i);
            
            // Make sure its ID hasn't been seen before
            REQUIRE(!seen_ids.count(node.id()));
            seen_ids.insert(node.id());
        }
    }

}

TEST_CASE( "Non-left-shifted variants can be used to construct valid graphs", "[constructor]" ) {

    auto vcf_data = R"(##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
ref	5	.	AAA	AAG,A	50	PASS	.	GT
)";

    auto fasta_data = R"(>ref
AAAAAAAAAAAAA
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 50, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    // Find all the edges
    unordered_set<pair<id_t, id_t>> edges;
    
    for (auto& edge : result.edge()) {
        // All the edges should be end to start
        REQUIRE(!edge.from_start());
        REQUIRE(!edge.to_end());
        
        pair<id_t, id_t> key = make_pair(edge.from(), edge.to());
    
        // Each edge must be unique
        REQUIRE(!edges.count(key));
        
        edges.insert(key);
    }

    SECTION("Paths follow edges") {
        for (auto& path : result.path()) {
            // For each path
            for (size_t i = 1; i < path.mapping_size(); i++) {
                // Scan alogn adjacent pairs of nodes
                id_t prev = path.mapping(i - 1).position().node_id();
                id_t here = path.mapping(i).position().node_id();
                
                // The edge must have been created.
                REQUIRE(edges.count(make_pair(prev, here)));
            }
        }
    }

}

TEST_CASE( "VG handles structural variants as expected"){
    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<DEL>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=DEL;END=29;CIPOS=0,3;CIEND=-3,0	GT)";

    auto vcf_with_alt_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	GCTTGGAAATTTTCTGGAGTT	G	99	PASS	AC=1;NA=1;NS=1;SVTYPE=DEL;END=29;CIPOS=0,3;CIEND=-3,0	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    SECTION("vg ignores symbolic SVs by default"){
        auto result = construct_test_graph(fasta_data, vcf_data, 10, false, false);

        unordered_map<size_t, string> expected_nodes;
        expected_nodes.insert({1, "CAAATAAGGC"});
        expected_nodes.insert({2, "TTGGAAATTT"});
        expected_nodes.insert({3, "TCTGGAGTTC"});
        expected_nodes.insert({4, "TATTATATTC"});
        expected_nodes.insert({5, "CAACTCTCTG"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected_nodes[node.id()]);
        }
    }

    SECTION("vg handles SVs with alt data by default, even if they have SV tags"){
        auto result = construct_test_graph(fasta_data, vcf_with_alt_data, 10, false, false);

        unordered_map<size_t, string> expected_nodes;
        expected_nodes.insert({1, "CAAATAAGG"});
        expected_nodes.insert({2, "CTTGGAAATT"});
        expected_nodes.insert({3, "TTCTGGAGTT"});
        expected_nodes.insert({4, "CTATTATATT"});
        expected_nodes.insert({5, "CCAACTCTCT"});
        expected_nodes.insert({6, "G"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected_nodes[node.id()]);
        }
    }



}



TEST_CASE( "An SV deletion is represented properly" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<DEL>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=DEL;END=29;CIPOS=0,3;CIEND=-3,0	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "CTTGGAAATT"});
        expected.insert({3, "TTCTGGAGTT"});
        expected.insert({4, "CTATTATATT"});
        expected.insert({5, "CCAACTCTCT"});
        expected.insert({6, "G"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 4);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE( "An SV deletion with matching SVLEN is represented properly" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<DEL>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=DEL;END=29;SVLEN=-20;CIPOS=0,3;CIEND=-3,0	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "CTTGGAAATT"});
        expected.insert({3, "TTCTGGAGTT"});
        expected.insert({4, "CTATTATATT"});
        expected.insert({5, "CCAACTCTCT"});
        expected.insert({6, "G"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 4);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE( "An SV deletion with only SVLEN is represented properly" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<DEL>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=DEL;SVLEN=-20;CIPOS=0,3;CIEND=-3,0	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "CTTGGAAATT"});
        expected.insert({3, "TTCTGGAGTT"});
        expected.insert({4, "CTATTATATT"});
        expected.insert({5, "CCAACTCTCT"});
        expected.insert({6, "G"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 4);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE( "An SV deletion with only SPAN is represented properly" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<DEL>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=DEL;SPAN=20;CIPOS=0,3;CIEND=-3,0	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "CTTGGAAATT"});
        expected.insert({3, "TTCTGGAGTT"});
        expected.insert({4, "CTATTATATT"});
        expected.insert({5, "CCAACTCTCT"});
        expected.insert({6, "G"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 4);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE( "An SV deletion with mismatching SVLEN is rejected" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<DEL>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=DEL;END=29;SVLEN=-19;CIPOS=0,3;CIEND=-3,0	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGGC"});
        expected.insert({2, "TTGGAAATTT"});
        expected.insert({3, "TCTGGAGTTC"});
        expected.insert({4, "TATTATATTC"});
        expected.insert({5, "CAACTCTCTG"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE( "A non-SV deletion is represented properly" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	GCTTGGAAATTTTCTGGAGTT	G	99	PASS	AC=1;NA=1;NS=1;	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "CTTGGAAATT"});
        expected.insert({3, "TTCTGGAGTT"});
        expected.insert({4, "CTATTATATT"});
        expected.insert({5, "CCAACTCTCT"});
        expected.insert({6, "G"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 4);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}


TEST_CASE("VG handles SV insertions", "[constructor]"){
    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<INS>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INS;SEQ=ACTG;SVLEN=4;CIPOS=0,3	GT)";

    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);

    SECTION("Nodes are as expected") {

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "ACTG"});
        expected.insert({3, "CTTGGAAATT"});
        expected.insert({4, "TTCTGGAGTT"});
        expected.insert({5, "CTATTATATT"});
        expected.insert({6, "CCAACTCTCT"});
        expected.insert({7, "G"});


        for (size_t i = 0; i < result.node_size(); i++){
            auto& node = result.node(i);
            REQUIRE(node.sequence() == expected[node.id()]);
        }
    }
    
    SECTION("Edges are as expected") {
    
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 3);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        edges_wanted.emplace(6, 7);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE("VG handles SV insertions with both SVLEN and END", "[constructor]"){
    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<INS>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INS;SEQ=ACTG;SVLEN=4;END=9;CIPOS=0,3	GT)";

    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);

    SECTION("Nodes are as expected") {

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "ACTG"});
        expected.insert({3, "CTTGGAAATT"});
        expected.insert({4, "TTCTGGAGTT"});
        expected.insert({5, "CTATTATATT"});
        expected.insert({6, "CCAACTCTCT"});
        expected.insert({7, "G"});


        for (size_t i = 0; i < result.node_size(); i++){
            auto& node = result.node(i);
            REQUIRE(node.sequence() == expected[node.id()]);
        }
    }
    
    SECTION("Edges are as expected") {
    
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 3);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        edges_wanted.emplace(6, 7);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE("VG handles SV insertions with no SVLEN or END", "[constructor]"){
    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<INS>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INS;SEQ=ACTG;CIPOS=0,3	GT)";

    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);

    SECTION("Nodes are as expected") {

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "ACTG"});
        expected.insert({3, "CTTGGAAATT"});
        expected.insert({4, "TTCTGGAGTT"});
        expected.insert({5, "CTATTATATT"});
        expected.insert({6, "CCAACTCTCT"});
        expected.insert({7, "G"});


        for (size_t i = 0; i < result.node_size(); i++){
            auto& node = result.node(i);
            REQUIRE(node.sequence() == expected[node.id()]);
        }
    }
    
    SECTION("Edges are as expected") {
    
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 3);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        edges_wanted.emplace(6, 7);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE("VG rejects SV insertions with bad ENDs", "[constructor]"){
    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<INS>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INS;SEQ=ACTG;SVLEN=4;END=18;CIPOS=0,3	GT)";

    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);

   SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGGC"});
        expected.insert({2, "TTGGAAATTT"});
        expected.insert({3, "TCTGGAGTTC"});
        expected.insert({4, "TATTATATTC"});
        expected.insert({5, "CAACTCTCTG"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE("VG handles SV insertions with misunderstood ENDs set to POS + SVLEN", "[constructor]"){
    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<INS>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INS;SEQ=ACTG;SVLEN=4;END=13;CIPOS=0,3	GT)";

    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);

    SECTION("Nodes are as expected") {

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "ACTG"});
        expected.insert({3, "CTTGGAAATT"});
        expected.insert({4, "TTCTGGAGTT"});
        expected.insert({5, "CTATTATATT"});
        expected.insert({6, "CCAACTCTCT"});
        expected.insert({7, "G"});


        for (size_t i = 0; i < result.node_size(); i++){
            auto& node = result.node(i);
            REQUIRE(node.sequence() == expected[node.id()]);
        }
    }
    
    SECTION("Edges are as expected") {
    
        unordered_set<pair<id_t, id_t>> edges_wanted;
        edges_wanted.emplace(1, 2);
        edges_wanted.emplace(1, 3);
        edges_wanted.emplace(2, 3);
        edges_wanted.emplace(3, 4);
        edges_wanted.emplace(4, 5);
        edges_wanted.emplace(5, 6);
        edges_wanted.emplace(6, 7);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // All the edges should be forward
            REQUIRE(!edge.from_start());
            REQUIRE(!edge.to_end());
            
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_pair(edge.from(), edge.to())));
        }
    }

}

TEST_CASE( "An SV inversion is represented properly" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<INV>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INV;END=29;CIPOS=0,3;CIEND=-3,0	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    // Inversions are like substitutions, so the POS base is included and inverted.

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGG"});
        expected.insert({2, "CTTGGAAATT"});
        expected.insert({3, "TTCTGGAGTT"});
        expected.insert({4, "CTATTATATT"});
        expected.insert({5, "CCAACTCTCT"});
        expected.insert({6, "G"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<tuple<id_t, bool, id_t, bool>> edges_wanted;
        edges_wanted.emplace(1, false, 2, false);
        edges_wanted.emplace(2, false, 3, false);
        edges_wanted.emplace(3, false, 4, false);
        edges_wanted.emplace(4, false, 5, false);
        edges_wanted.emplace(5, false, 6, false);
        edges_wanted.emplace(1, false, 3, true);
        edges_wanted.emplace(2, true, 4, false);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_tuple(edge.from(), edge.from_start(), edge.to(), edge.to_end())));
        }
    }

}

TEST_CASE( "An SV inversion that ends 1 base after something else is constructable" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	9	sv1	N	<INV>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INV;END=29;CIPOS=0,3;CIEND=-3,0	GT
x	11	sv2	N	<DEL>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=DEL;END=40;CIPOS=0,3;CIEND=-3,0	GT
x	28	sv2	T	GATTACA	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INS;SVLEN=7;CIPOS=0,3;CIEND=-3,0	GT
)";


    //          v    inverted      v
    // CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
    //            ^      deleted              ^
    //                            X Replaced and inserted after


    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 100, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    // Inversions are like substitutions, so the POS base is included and inverted.

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        // Order is a bit weird. First part before insertion.
        expected.insert({1, "CAAATAAGG"});
        // Replacing inserted base. TODO: Should we change this to join with the rest of the insert?
        expected.insert({2, "G"});
        // Rest of inserted bases
        expected.insert({3, "ATTACA"});
        // Inversion before deletion
        expected.insert({4, "CT"});
        // Deleted sequence before insertion
        expected.insert({5, "TGGAAATTTTCTGGAG"});
        // Insertion base replaced
        expected.insert({6, "T"});
        // Base in inversion after insertion
        expected.insert({7, "T"});
        // Non-inverted deleted part
        expected.insert({8, "CTATTATATTC"});
        // Part after inversion
        expected.insert({9, "CAACTCTCTG"});
        
        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
}

TEST_CASE( "SV inversions with smart quotes and a lower-case reference parse correctly" , "[constructor]") {

    // Note the smart quotes

    auto vcf_data = R"(##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description=“End position of the variant described in this record”>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=“Type of structural variant”>
##ALT=<ID=IV,Description=“Inversion”>
#CHROM  POS     ID      REF ALT QUAL    FILTER  INFO
x	10	SRR026655.22810753-B	N	<INV>	.	.	END=41;SVTYPE=INV)";

    auto fasta_data = R"(>x
caaataaggcttggaaattttctggagttctattatattccaactctctg
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    // Inversions are like substitutions, so the POS base is included and inverted.

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGGC"});
        expected.insert({2, "TTGGAAATTT"});
        expected.insert({3, "TCTGGAGTTC"});
        expected.insert({4, "TATTATATTC"});
        expected.insert({5, "C"});
        expected.insert({6, "AACTCTCTG"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<tuple<id_t, bool, id_t, bool>> edges_wanted;
        edges_wanted.emplace(1, false, 2, false);
        edges_wanted.emplace(2, false, 3, false);
        edges_wanted.emplace(3, false, 4, false);
        edges_wanted.emplace(4, false, 5, false);
        edges_wanted.emplace(5, false, 6, false);
        edges_wanted.emplace(1, false, 5, true);
        edges_wanted.emplace(2, true, 6, false);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_tuple(edge.from(), edge.from_start(), edge.to(), edge.to_end())));
        }
    }

}


TEST_CASE( "A shorter SV inversion is represented properly" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	10	sv1	N	<INV>	99	PASS	AC=1;NA=1;NS=1;SVTYPE=INV;END=29;CIPOS=0,3;CIEND=-3,0	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

    // Inversions are like substitutions, so the POS base is included and inverted.

    SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGGC"});
        expected.insert({2, "TTGGAAATTT"});
        expected.insert({3, "TCTGGAGTT"});
        expected.insert({4, "CTATTATATT"});
        expected.insert({5, "CCAACTCTCT"});
        expected.insert({6, "G"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<tuple<id_t, bool, id_t, bool>> edges_wanted;
        edges_wanted.emplace(1, false, 2, false);
        edges_wanted.emplace(2, false, 3, false);
        edges_wanted.emplace(3, false, 4, false);
        edges_wanted.emplace(4, false, 5, false);
        edges_wanted.emplace(5, false, 6, false);
        edges_wanted.emplace(1, false, 3, true);
        edges_wanted.emplace(2, true, 4, false);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_tuple(edge.from(), edge.from_start(), edge.to(), edge.to_end())));
        }
    }

}

TEST_CASE( "SVs that are fully base specified are constructed correctly" , "[constructor]") {

    auto vcf_data = R"(##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
x	11	sv1	T	TTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTC	11	PASS	END=10;SVLEN=35;SVTYPE=INS	GT)";

    auto fasta_data = R"(>x
CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG
)";

    // Build the graph
    auto result = construct_test_graph(fasta_data, vcf_data, 10, true, false);
    
#ifdef debug
    std::cerr << pb2json(result) << std::endl;
#endif

	SECTION("nodes are as expected") {
        // Look at each node

        unordered_map<size_t, string> expected;
        expected.insert({1, "CAAATAAGGC"});
		expected.insert({2, "T"});
        expected.insert({3, "TTCTTTCTTT"});
        expected.insert({4, "CTTTCTTTCT"});
        expected.insert({5, "TTCTTTCTTT"});
        expected.insert({6, "CTTTC"});
        expected.insert({7, "TGGAAATTTT"});
        expected.insert({8, "CTGGAGTTCT"});
        expected.insert({9, "ATTATATTCC"});
        expected.insert({10, "AACTCTCTG"});

        for (size_t i = 0; i < result.node_size(); i++) {
            auto& node = result.node(i);
            REQUIRE(node.sequence()==expected[node.id()]);
        }
    }
    
    SECTION("edges are as expected") {
        unordered_set<tuple<id_t, bool, id_t, bool>> edges_wanted;
        edges_wanted.emplace(1, false, 2, false);
        edges_wanted.emplace(2, false, 3, false);
        edges_wanted.emplace(3, false, 4, false);
        edges_wanted.emplace(4, false, 5, false);
        edges_wanted.emplace(5, false, 6, false);
        edges_wanted.emplace(6, false, 7, false);
        edges_wanted.emplace(7, false, 8, false);
        edges_wanted.emplace(8, false, 9, false);
        edges_wanted.emplace(9, false, 10, false);
        edges_wanted.emplace(2, false, 7, false);
        
        // We should have the right number of edges
        REQUIRE(result.edge_size() == edges_wanted.size());
        
        for (auto& edge : result.edge()) {
            // The edge should be expected
            REQUIRE(edges_wanted.count(make_tuple(edge.from(), edge.from_start(), edge.to(), edge.to_end())));
        }
    }

}




}
}
