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
        add_feature(&aln, FeatureType::HAPLOTYPE_LOGPROB, 1.0);
        
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_LOGPROB));
        REQUIRE(get_feature(aln, FeatureType::HAPLOTYPE_LOGPROB) == 1.0);
        
        SECTION("Feature can be removed") {
            remove_feature(&aln, FeatureType::HAPLOTYPE_LOGPROB);
            REQUIRE(!has_feature(aln, FeatureType::HAPLOTYPE_LOGPROB));
        }
    }
    
    SECTION("Feature can be set") {
        set_feature(&aln, FeatureType::HAPLOTYPE_LOGPROB, 1.0);
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_LOGPROB));
        REQUIRE(get_feature(aln, FeatureType::HAPLOTYPE_LOGPROB) == 1.0);
        
        set_feature(&aln, FeatureType::HAPLOTYPE_LOGPROB, 15);
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_LOGPROB));
        REQUIRE(get_feature(aln, FeatureType::HAPLOTYPE_LOGPROB) == (double)15);
        
        set_feature(&aln, FeatureType::HAPLOTYPE_LOGPROB, 2.0);
        REQUIRE(has_feature(aln, FeatureType::HAPLOTYPE_LOGPROB));
        REQUIRE(get_feature(aln, FeatureType::HAPLOTYPE_LOGPROB) == 2.0);
        
        remove_feature(&aln, FeatureType::HAPLOTYPE_LOGPROB);
        REQUIRE(!has_feature(aln, FeatureType::HAPLOTYPE_LOGPROB));
    }
}

TEST_CASE("Multi-value features work", "[features][alignment]") {
    Alignment aln;
    
    SECTION("Features start not present") {
        REQUIRE(!has_feature(aln, FeatureType::SECONDARY_SCORES));
    }
    
    SECTION("Feature is present after addition") {
        add_feature(&aln, FeatureType::SECONDARY_SCORES, 1.0);
        
        REQUIRE(has_feature(aln, FeatureType::SECONDARY_SCORES));
        REQUIRE(get_feature(aln, FeatureType::SECONDARY_SCORES) == 1.0);
        auto vals = get_features(aln, FeatureType::SECONDARY_SCORES);
        REQUIRE(vals.size() == 1);
        REQUIRE(vals[0] == 1.0);
        
        SECTION("Feature can be removed") {
            remove_feature(&aln, FeatureType::SECONDARY_SCORES);
            REQUIRE(!has_feature(aln, FeatureType::SECONDARY_SCORES));
            
            auto vals = get_features(aln, FeatureType::SECONDARY_SCORES);
            REQUIRE(vals.size() == 0);
        }
        
        SECTION("More values can be added") {
            add_feature(&aln, FeatureType::SECONDARY_SCORES, 2.0);
            add_feature(&aln, FeatureType::SECONDARY_SCORES, 3.0);
            
            REQUIRE(get_feature(aln, FeatureType::SECONDARY_SCORES) == 1.0);
            auto vals = get_features(aln, FeatureType::SECONDARY_SCORES);
            REQUIRE(vals.size() == 3);
            REQUIRE(vals[0] == 1.0);
            REQUIRE(vals[1] == 2.0);
            REQUIRE(vals[2] == 3.0);
            
            SECTION("All values can be removed") {
                remove_feature(&aln, FeatureType::SECONDARY_SCORES);
                REQUIRE(!has_feature(aln, FeatureType::SECONDARY_SCORES));
                
                auto vals = get_features(aln, FeatureType::SECONDARY_SCORES);
                REQUIRE(vals.size() == 0);
            }
        }
    }
}


}
}

