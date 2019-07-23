#include "mcmc_genotyper.hpp"
#include <algorithm>    // std::min
#include <random> 
#include <chrono>

namespace vg {

    using namespace std;

    MCMCGenotyper::MCMCGenotyper(SnarlManager& snarls, VG& graph, const int n_iterations):snarls(snarls), graph(graph), n_iterations(n_iterations), 
        random_engine(chrono::system_clock::now().time_since_epoch().count()){
            // will want to implement a way for callee to pass in the seed value for predictability in unit testing 
    
    }

    PhasedGenome MCMCGenotyper::run_genotype(const vector<MultipathAlignment>& reads, const double log_base) const{

        // generate initial value 
        PhasedGenome initial_guess = generate_initial_guess();
        PhasedGenome* current = &initial_guess;
    
        for(int i = 1; i<= n_iterations; i++){
        
            // set x_star 
            PhasedGenome x_star = proposal_sample(*current);
            
            // calculate likelihood ratio of posterior distribution 
            double likelihood_ratio = exp(log_base*(log_target(x_star, reads)-log_target(*current, reads)));

            // calculate acceptance probability 
            double acceptance_probability = min(1.0, likelihood_ratio);
        
            if(runif(0.0,1.0) < acceptance_probability){ 
                current = &x_star;
            }else{
                // current value is kept
            }

        }

        return PhasedGenome(snarls);


    }   
    double MCMCGenotyper::log_target(PhasedGenome& phased_genome, const vector<MultipathAlignment>& reads)const{
        
        // Sum of scores given the reads aligned on the haplotype 
        double sum_scores; 
        
        for(MultipathAlignment mp : reads){
            sum_scores += phased_genome.optimal_score_on_genome(mp, graph);
        } 
        return sum_scores;
    }

    PhasedGenome MCMCGenotyper::proposal_sample(PhasedGenome& current)const{
        
        // set a seed for predictability in testing 
        unsigned seed = 1; 
        minstd_rand0 random_engine(seed);
        
        // Sample uniformly within snarls 
        const Snarl* random_snarl = snarls.discrete_uniform_sample(random_engine);

        // Get list of haplotypes that contain the snarl
        vector<int> matched_haplotypes = current.get_haplotypes_with_snarl(random_snarl);

        // pick a matched haplotype uniformly 
        // write algorithm that gets a different traversal through the snarl by uniformly choosing from all possible ways to traverse the snarl
        // look at the video
        // set allele - use random alleles - this methods replaces with new 
        // swap traversals in snarls not between snarls 
        // get allele - save and return 
        // book keeping: tuple (haplotype ID, snarl* (site that we replaced at), get_allele(old allele, which is the old traversal of the snarl))
        // will not return entire genome because it will be too costly
        return PhasedGenome(snarls);

    }
    double MCMCGenotyper::runif(const double a, const double b)const{
        
        uniform_real_distribution<double> distribution(a,b);
        double random_num = distribution(random_engine);
        
        return random_num;

    }
    PhasedGenome MCMCGenotyper::generate_initial_guess()const{
        
        PhasedGenome genome = PhasedGenome(snarls);
        vector<NodeTraversal> haplotype; //will add twice  

        graph.for_each_path_handle([&](const path_handle_t& path){
        // For each path

            if(!Paths::is_alt(graph.get_path_name(path))) {
            // If it isn't an alt path, we want to trace it
  
                for (handle_t handle : graph.scan_path(path)) {
                // For each occurrence from start to end
                    
                    // get the node and the node postion and add to the vector 
                    Node* node = graph.get_node(graph.get_id(handle));
                    bool position = graph.get_is_reverse(handle);
                    haplotype.push_back(NodeTraversal(node,position));
  
                }
            }
        });
        // construct haplotypes
        // haplotype1 = haplotype2
        genome.add_haplotype(haplotype.begin(), haplotype.end());
        genome.add_haplotype(haplotype.begin(), haplotype.end());

        return genome;
    }
}



