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
    SnarlManager& snarls; 
    VG& graph;
    const int n_iterations = 1000;
    mutable minstd_rand0 random_engine;

public:
    
    MCMCGenotyper(SnarlManager& snarls, VG& graph, const int n_iterations); 

    /** 
     * Takes as input a collection of mapped reads stored as a vector of multipath alignments and uses 
     * MCMC to find two optimal paths through the graph.
     * Output: phased genome 
     */    
    PhasedGenome run_genotype(const vector<MultipathAlignment>& reads, const double log_base) const;
    
    /**
     * Represents the poseterior distribution function 
     * returns the posterir probability
     */ 
     double log_target(PhasedGenome& phased_genome, const vector<MultipathAlignment>& reads) const;  

    /**
     * Generates a proposal sample over the desired distrubution
     * returns a sample from the proposal distribution
     */
     tuple<id_t, Snarl*, vector<NodeTraversal> > proposal_sample(PhasedGenome& current) const;
    /**
     * Samples haplotypes randomly using the discrete uniform distribution
     */
     int sample_uniform_haplotypes(minstd_rand0& random_engine, vector<id_t> matched_haplotypes) const;
    
    /**
     * Given a range [a,b] will return a random number uniformly distributed within that range 
     */
     double runif(const double a, const double b) const;

     /**
      * Generate a PhasedGenome to use as an initial value in M-H
      * Uses the two non-alt paths from the linear reference as haplotypes
      */
     PhasedGenome generate_initial_guess(void) const;
};

}



#endif
