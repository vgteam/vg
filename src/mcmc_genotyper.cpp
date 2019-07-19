#include "mcmc_genotyper.hpp"
#include <algorithm>    // std::min
#include <random> 

namespace vg {

    using namespace std;

    MCMCGenotyper::MCMCGenotyper(const SnarlManager& snarls, VG& graph):snarls(snarls), graph(graph){

    }

    PhasedGenome MCMCGenotyper::run_genotype(const vector<MultipathAlignment>& reads) const{

        
        PhasedGenome current = generate_initial_guess();
    
        for(int i = 1; i<= 10000; i++){
        
        // set x_star 
            PhasedGenome x_star = proposal_sample(current);

            double ratio = std::min(1.0, target(x_star, reads)/target(current, reads));
        
            if(runif(0.0,1.0) < ratio){ 
            //accept x_star - what does it mean to accept this PhasedGenome ?
            }else{
            //accept cuurent 
            }

        }

        return PhasedGenome(snarls);


    }   
    double MCMCGenotyper::target(const PhasedGenome& phased_genome, const vector<MultipathAlignment>& reads)const{

    
    }

    PhasedGenome MCMCGenotyper::proposal_sample(const PhasedGenome& current)const{

    }
    double MCMCGenotyper::runif(const double a, const double b)const{
    
        default_random_engine generator;
        uniform_real_distribution<double> distribution(a,b);
        double random_num = distribution(generator);
        
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



