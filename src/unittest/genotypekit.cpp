/**
 * unittest/genotypekit.cpp: test cases for genotypekit modular genotyper pieces
 */

#include "catch.hpp"
#include "genotypekit.hpp"

namespace vg {
namespace unittest {

TEST_CASE( "fixed priors can be assigned to genotypes", "[genotype]" ) {
    
    GenotypePriorCalculator* calculator = new FixedGenotypePriorCalculator();
    
    Genotype het;
    het.add_allele(0);
    het.add_allele(1);
    
    Genotype hom_alt;
    hom_alt.add_allele(1);
    hom_alt.add_allele(1);
    
    Genotype hom_ref;
    hom_ref.add_allele(0);
    hom_ref.add_allele(0);
    
    SECTION("homozygote priors should be equal") {
        REQUIRE(calculator->calculate_log_prior(hom_alt) == calculator->calculate_log_prior(hom_ref));
    }
    
    SECTION("homozygotes should be more likely than heterozygotes") {
        REQUIRE(calculator->calculate_log_prior(het) < calculator->calculate_log_prior(hom_ref));
        REQUIRE(calculator->calculate_log_prior(het) < calculator->calculate_log_prior(hom_alt));
    }
    
    delete calculator;
}

}
}
