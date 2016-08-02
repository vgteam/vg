#include "genotypekit.hpp"

namespace vg {

using namespace std;


CactusSiteFinder::CactusSiteFinder(VG& graph, const string& hint_path_name): graph(graph), hint_path_name(hint_path_name) {
    // Nothing to do!
}

void CactusSiteFinder::for_each_site_parallel(const function<void(const NestedSite&)>& lambda) {
    
}

double FixedGenotypePriorCalculator::calculate_log_prior(const Genotype& genotype) {
    // Are all the alleles the same?
    bool all_same = true;
    // What is the common allele number (-1 for unset)
    int allele_value = -1;
    for(size_t i = 0; i < genotype.allele_size(); i++) {
        // For each allele in the genotype
        if(allele_value == -1) {
            // On the first one, grab it. Everyone else will have to match.
            allele_value = genotype.allele(i);
        }
        
        if(allele_value != genotype.allele(i)) {
            // There are two distinct allele values in here
            all_same = false;
            break;
        }
    }
    
    // Return the appropriate prior depending on whether the alleles are all the
    // same (homozygous) or not (heterozygous).
    return all_same ? homozygous_prior_ln : heterozygous_prior_ln;
}

}
