#include "mcmc_genotyper.hpp"


namespace vg {

using namespace std;

PhasedGenome MCMCGenotyper::run_genotype(const vector<MultipathAlignment>& reads) const{
    // TODO: return a phase genome 
    return PhasedGenome(snarls);


}

}
