#ifndef VG_MCMC_GENOTYPER_HPP_INCLUDED
#define VG_MCMC_GENOTYPER_HPP_INCLUDED

/** \file
 * mcmc_genotyper.hpp: defines a class that implements mcmc probabilistic model for variant calling.
 */

#include <map>
#include <string>
#include <vg/vg.pb.h>
#include <vector>
#include "phased_genome.hpp"


namespace vg {

using namespace std;
    
/** 
 * This class is a genotyper that uses MCMC to find two optimal paths through the graph given a set of aligned reads. 
 *  
 */ 
class MCMCGenotyper{
    SnarlManager snarls;    

public:

    /** 
     * This method takes as input a collection of mapped reads stored as a vector of multipath alignments and uses 
     * MCMC to find two optimal paths through the graph.
     * Output: phased genome 
     */    
    PhasedGenome run_genotype(const vector<MultipathAlignment>& reads) const;
    
     

};

}



#endif
