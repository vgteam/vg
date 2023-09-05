/**
 * \file 
 * unittest/vcf_buffer.cpp: test cases for the VcfBuffer and WindowedVcfBuffer
 * which let us read through VCFs with context.
 */

#include "catch.hpp"
#include "../vcf_buffer.hpp"

#include "../utility.hpp"
#include "../path.hpp"
#include "vg/io/json2pb.h"

#include <vector>
#include <sstream>
#include <iostream>

namespace vg {
namespace unittest {

TEST_CASE( "WindowedVcfBuffer windowing works", "[windowedvcfbuffer][vcf]" ) {

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
ref	7	rs1338	A	G	29	PASS	.	GT
ref	7	rs1339	A	T	29	PASS	.	GT
ref	8	rs1340	A	G	29	PASS	.	GT
ref	17	rs1341	A	G	29	PASS	.	GT
ref	18	rs1342	A	G	29	PASS	.	GT
)";

    // Make a stream out of the data
    std::stringstream vcf_stream(vcf_data);
    
    // Load it up in vcflib
    vcflib::VariantCallFile vcf;
    vcf.open(vcf_stream);
    
    // Make a windowed buffer
    WindowedVcfBuffer buffer(&vcf, 10);
    
    // Definbe some variables for holding variants and their contexts
    vector<vcflib::Variant*> before;
    vector<vcflib::Variant*> after;
    vcflib::Variant* current;
    
    SECTION("no tabix index is found") {
        REQUIRE(!buffer.has_tabix());
    }
    
    SECTION("some variants exist") {
        REQUIRE(buffer.next());
        
        SECTION("the first variant can be retrieved in context") {
            tie(before, current, after) = buffer.get();
            
            REQUIRE(current != nullptr);
            REQUIRE(current->id == "rs1337");
            
            SECTION("there are no variants before the first") {
                REQUIRE(before.empty());
            }
            
            SECTION("there are 3 variants within 10 bp after the first") {
                REQUIRE(after.size() == 3);
            }
        }
        
        SECTION("the second variant can be retrieved in context") {
            REQUIRE(buffer.next());
            tie(before, current, after) = buffer.get();
            
            REQUIRE(current != nullptr);
            REQUIRE(current->id == "rs1338");
            
            SECTION("there is one variant before the second") {
                REQUIRE(before.size() == 1);
                
                SECTION("it is the first variant") {
                    REQUIRE(before.front()->id == "rs1337");
                }
            }
            
            SECTION("there are 3 variants within 10 bp of the second") {
                REQUIRE(after.size() == 3);
            }
            
            SECTION("the second variant can be retrieved without overlaps") {
                tie(before, current, after) = buffer.get_nonoverlapping();
                
                REQUIRE(current != nullptr);
                REQUIRE(current->id == "rs1338");
                
                SECTION("there are only 2 nonoverlapping subsequent variants in range") {
                    REQUIRE(after.size() == 2);
                    REQUIRE(after.front()->id == "rs1340");
                }
            }
            
            SECTION("the fifth variant can be retrieved in context without overlaps") {
                REQUIRE(buffer.next()); // 3rd
                REQUIRE(buffer.next()); // 4th
                REQUIRE(buffer.next()); // 5th
                
                tie(before, current, after) = buffer.get_nonoverlapping();
                
                REQUIRE(current != nullptr);
                REQUIRE(current->id == "rs1341");
                
                SECTION("there are 2 variants before it in range that don't overlap it or each other") {
                    REQUIRE(before.size() == 2);
                }
                
                SECTION("there is 1 variant after it") {
                    REQUIRE(after.size() == 1);
                }
                
                SECTION("iteration can be exhausted") {
                    REQUIRE(buffer.next()); // 6th
                    REQUIRE(!buffer.next()); // No 7th
                    REQUIRE_THROWS(buffer.get());
                }
            }
            
        }
    }
    
}

}
}
