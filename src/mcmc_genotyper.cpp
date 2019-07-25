#include "mcmc_genotyper.hpp"
#include "algorithms/count_walks.hpp"
#include "subgraph.hpp"
#include <algorithm> 
#include <random> 
#include <chrono>
#include <utility>

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
        
        // sum of scores given the reads aligned on the haplotype 
        double sum_scores; 
        
        for(MultipathAlignment mp : reads){
            sum_scores += phased_genome.optimal_score_on_genome(mp, graph);
        } 
        return sum_scores;
    }

    tuple<id_t, Snarl*, vector<NodeTraversal> > MCMCGenotyper::proposal_sample(PhasedGenome& current)const{
        
        // set a seed for predictability in testing 
        unsigned seed = 1; 
        minstd_rand0 random_engine(seed);
        
        // sample uniformly between snarls 
        const Snarl* random_snarl = snarls.discrete_uniform_sample(random_engine);

        // get list of haplotypes that contain the snarl
        vector<id_t> matched_haplotypes = current.get_haplotypes_with_snarl(random_snarl);

        // choose a haplotype uiformly 
        id_t random_haplotype = sample_uniform_haplotypes(random_engine, matched_haplotypes);

        // returns a pair with Nodes and Edges contained in Snarl
        pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarls.deep_contents(random_snarl, graph, false);

        // unpack the pair, we only care about the nodes
        unordered_set<id_t>& ids = contents.first;
            
        // HandleGraph of vg graph
        
        // build a subhandle graph from super graph
        // we only want to enumerate counts through nodes in a snarl not the entire graph
        SubHandleGraph subgraph(&graph);
        
        for (id_t id : ids){
            
            // for each node_id in the snarl, get handle
            // add the handle to the subgraph
            subgraph.add_handle(subgraph.get_handle(id, false));
        }
        
        unordered_map<handle_t, size_t> count_map = algorithms::count_walks_through_nodes(&graph);
        size_t sink_totals = algorithms::get_total();
        
        // bookkeeping: haplotype ID, snarl* (site that we replaced at), get_allele())
        tuple<id_t, Snarl*, vector<NodeTraversal> > unpacked_pg;
        
        //vector<NodeTraversal> allele;
        //push back the nodes in the allele

        // set new allele, replace with previous allele 
        //current.set_allele(random_snarl , allele.begin(), allele.end(), random_haplotype);

        // old traversal of the snarl 
        // get allele at site where replacement was made
        //vector<NodeTraversal> old_allele = current.get_allele(random_snarl, random_haplotype);

        
        
        // get a different traversal through the snarl by uniformly choosing from all possible ways to traverse the snarl
        // set allele - use random alleles - this methods replaces with new 
        // swap traversals in snarls (not between snarls)
 
        return unpacked_pg;

    }
    int MCMCGenotyper::sample_uniform_haplotypes(minstd_rand0& random_engine, vector<id_t> matched_haplotypes) const{
        int number_of_haplotypes = matched_haplotypes.size;
    
        // choose a haplotype randomly using discrete uniform distribution
        uniform_int_distribution<int> distribution(0, number_of_haplotypes);  
        int random_num = distribution(random_engine);

        int random_haplotype = matched_haplotypes[random_num];

        return random_haplotype;
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
        // capture all variables (paths) in scope by reference 

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



