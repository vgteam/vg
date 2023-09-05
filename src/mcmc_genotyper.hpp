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
#include "multipath_alignment.hpp"
#include "algorithms/min_cut_graph.hpp"
#include "sparse_union_find.hpp"

namespace vg {

using namespace std;
   
/** 
 * This class is a genotyper that uses MCMC to find two optimal paths through the graph given a set of aligned reads. 
 *  
 */ 
class MCMCGenotyper{
    SnarlManager& snarls; 
    VG& graph;
    const int n_iterations;
    const int seed;
    const int burn_in; //how many iterations we run only original proposal sample dist. 
    const int frequency; //frequency we remake gamma set
    mutable minstd_rand0 random_engine;

public:
    
    MCMCGenotyper(SnarlManager& snarls, VG& graph, const int n_iterations, const int seed, const int burn_in, const int frequency); 

    /** 
     * Takes as input a collection of mapped reads stored as a vector of multipath alignments and uses 
     * MCMC to find two optimal paths through the graph.
     * Output: phased genome 
     */    
    unique_ptr<PhasedGenome> run_genotype(const vector<multipath_alignment_t>& reads, const double log_base) const;
    
    /**
     * Represents the poseterior distribution function 
     * returns the posterir probability
     */ 
     double log_target(unique_ptr<PhasedGenome>& phased_genome, const vector<multipath_alignment_t>& reads) const;  

    /**
     * Generates a proposal sample over the desired distrubution
     * returns a sample from the proposal distribution
     */
     tuple<int, const Snarl*, vector<NodeTraversal> > proposal_sample(unique_ptr<PhasedGenome>& current) const;
    /**
     * Generates a number randomly using the discrete uniform distribution
     */
     int generate_discrete_uniform(minstd_rand0& random_engine, int lower_bound , int upper_bound) const;
    
    /**
     * Given a range [a,b] will return a random number uniformly distributed within that range 
     */
     double generate_continuous_uniform(const double a, const double b) const;

     /**
      * Generate a PhasedGenome to use as an initial value in M-H
      * Uses the two non-alt paths from the linear reference as haplotypes
      */
     unique_ptr<PhasedGenome> generate_initial_guess()const;
     /**
      * Generate a map from a pair of snarls to an edge weight 
      * Uses snarls read API, reads and the optimal_score_on pahsed genome as a 
      * scoring scheme for the edge weight overlapping the snarl pair.
      */
     unordered_map<pair<const Snarl*, const Snarl*>, int32_t> make_snarl_map(const vector<multipath_alignment_t>& reads, PhasedGenome& phased_genome) const;
     /**
      * Generate a graph using the snarl map
      */
     algorithms::Graph make_snarl_graph(unordered_map<pair<const Snarl*, const Snarl*>, int32_t> map) const;

    /**
      * Make a snarl graph with edge weights scored by how well mapped reads support phasing of snarl
      * Use an alternative proposal distribution using sets generated from karger-stein min cut algorithm
      * to escape bottlenecks leading to rapid convergence
      */
     vector<unordered_set<size_t>> karger_stein(const vector<multipath_alignment_t>& reads, PhasedGenome& genome) const;

     unordered_set<size_t> alt_proposal_sample(vector<unordered_set<size_t>>& gamma, PhasedGenome& genome) const;


};

}



#endif
