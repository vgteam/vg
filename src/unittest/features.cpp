/** \file
 *
 * Unit tests for the Features API, which lets us attach metadata to Alignments and MultipathAlignments.
 */

#include <iostream>
#include <sstream>
#include "../features.hpp"

#include "catch.hpp"

namespace vg {
namespace unittest {

using namespace std;

TEST_CASE("Tag features work", "[features][alignment]") {

    Alignment aln;
    
    SECTION("Feature starts not present") {
        REQUIRE(!has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
    }
    
    SECTION("Feature is present after addition") {
        add_feature(&aln, FeatureType::HAPLOTYPE_SCORED_TAG);
        
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
        
        SECTION("Feature can be removed again") {
            remove_feature(&aln, FeatureType::HAPLOTYPE_SCORED_TAG);
            
            REQUIRE(!has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
        }
    }
    
    SECTION("Feature can be added with a bool value") {
        SECTION("False does not set it") {
            add_feature(&aln, FeatureType::HAPLOTYPE_SCORED_TAG, false);
            REQUIRE(!has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
        }
        
        SECTION("True sets it") {
            add_feature(&aln, FeatureType::HAPLOTYPE_SCORED_TAG, true);
            REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
        }
    }
    
    SECTION("Feature can be set") {
        add_feature(&aln, FeatureType::HAPLOTYPE_SCORED_TAG);
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
        
        set_feature(&aln, FeatureType::HAPLOTYPE_SCORED_TAG, true);
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
        
        set_feature(&aln, FeatureType::HAPLOTYPE_SCORED_TAG, false);
        REQUIRE(!has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
        
        set_feature(&aln, FeatureType::HAPLOTYPE_SCORED_TAG, true);
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_SCORED_TAG));
    }
}

TEST_CASE("Single-value features work", "[features][alignment]") {
    
    Alignment aln;
    
    SECTION("Features start not present") {
        REQUIRE(!has_feature(aln, FeatureType::HAPLOTYPE_LOGPROB));
    }
    
    SECTION("Feature is present after addition") {
        add_feature(aln, FeatureType::HAPLOTYPE_LOGPROB, 1.0);
        
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_LOGPROB));
        REQUIRE(get_feature(aln, FeatureType::HAPLOTYPE_LOGPROB) == 1.0);
    }
}


}
}

