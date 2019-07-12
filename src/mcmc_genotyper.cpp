#include "mcmc_genotyper.hpp"


namespace vg {

using namespace std;

MCMCGenotyper::MCMCGenotyper(const SnarlManager& snarls):snarls(snarls){

}

PhasedGenome MCMCGenotyper::run_genotype(const vector<MultipathAlignment>& reads) const{
    // TODO: return a phase genome 
    return PhasedGenome(snarls);


}

}
