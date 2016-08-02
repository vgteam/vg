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
    
    SECTION("haploid genotypes should have nonzero prior") {
        Genotype haploid;
        haploid.add_allele(5);
        REQUIRE(calculator->calculate_log_prior(haploid) > prob_to_logprob(0));
    }
    
    SECTION("zero-ploid genotypes should have nonzero prior") {
        Genotype empty;
        REQUIRE(calculator->calculate_log_prior(empty) > prob_to_logprob(0));
    }
    
    SECTION("polyploid genotypes should have nonzero prior") {
        Genotype polyploid;
        for(int i = 0; i < 100; i++) {
            polyploid.add_allele(i);
        }
        REQUIRE(calculator->calculate_log_prior(polyploid) > prob_to_logprob(0));
    }
    
    delete calculator;
}

}
}
