#include "mcmc_genotyper.hpp"
#include "algorithms/count_walks.hpp"
#include "subgraph.hpp"
#include <algorithm> 
#include <random> 
#include <chrono>
#include <utility>
#include "multipath_alignment.hpp"

namespace vg {

    using namespace std;

    MCMCGenotyper::MCMCGenotyper(SnarlManager& snarls, VG& graph, const int n_iterations, const int seed):snarls(snarls), graph(graph), n_iterations(n_iterations), 
        seed(seed), random_engine(seed){
            
    
    }

    unique_ptr<PhasedGenome> MCMCGenotyper::run_genotype(const vector<MultipathAlignment>& reads, const double log_base) const{

        // set a flag for invalid contents so a message it observed 
        bool invalid_contents = false;
        bool return_optimal = false;
        bool swapped = false;

        // generate initial value
        unique_ptr<PhasedGenome> genome = generate_initial_guess();
        
        double max_likelihood = 0.0; //max log likelihood
        double current_likelihood = 0.0; // current log likelihood
        double previous_likelihood = 0.0;

        unique_ptr<PhasedGenome> optimal;
        
        // build markov chain using Metropolis-Hastings
        for(int i = 0; i< n_iterations; i++){
            
            // holds the previous sample allele
            double x_prev = log_target(*genome, reads);
        
            // get contents from proposal_sample
            tuple<int, const Snarl*, vector<NodeTraversal> > to_receive = proposal_sample(*genome);
            int& modified_haplo = get<0>(to_receive);         
            
            // if the to_receive contents are invalid keep the new allele
            // for graphs that do not contain snarls
            if (modified_haplo ==-1){
                invalid_contents = true;
                continue;
            }else{
                const Snarl& modified_site = *get<1>(to_receive); 
                vector<NodeTraversal>& old_allele = get<2>(to_receive); 
 
                // holds new sample allele
                double x_new = log_target(*genome, reads);

                //genome->print_phased_genome();

                // calculate likelihood ratio of posterior distribution 
                double likelihood_ratio = exp(log_base*(x_new - x_prev));
                
                // cerr << "********************************" <<endl;
                // cerr << "total score of new genome "<< x_new << endl;
                // cerr << "total score of prev genome "<< x_prev << endl;
                // cerr << "likelihood ratio            " << likelihood_ratio << endl;
                // cerr << "prev log likelihood         " << current_likelihood << endl;                
                
                if(!swapped){
                    // if we did not swap back then we use the prev likelihood 
                    previous_likelihood = current_likelihood;
                }// otherwise don't update it 

                current_likelihood = previous_likelihood + log_base*(x_new-x_prev);

                // cerr << "current log likelihood      " << current_likelihood << endl;
                // cerr << "max likelihood              " << max_likelihood << endl;
                // cerr<< endl;
                
                
                // genome->print_phased_genome();
                // cerr << "********************************" <<endl;
                if (current_likelihood > max_likelihood){
                    max_likelihood = current_likelihood;
                    optimal = unique_ptr<PhasedGenome>(new PhasedGenome(*genome));
                    return_optimal=true;
                }

                // calculate acceptance probability 
                double acceptance_probability = min(1.0, likelihood_ratio);

                // if u~U(0,1) > alpha, discard new allele and keep previous 
                if(generate_continuous_uniform(0.0,1.0) > acceptance_probability){ 
                    swapped = true;
                    genome->set_allele(modified_site, old_allele.begin(), old_allele.end(), modified_haplo); 
                    // cerr << "************UPDATED HAPLOTYPE********************" <<endl;
                    // genome->print_phased_genome();
                    // cerr << "********************************" <<endl;
                }         
            }
        } 
        if(invalid_contents || !return_optimal){
            // for graphs without snarls 
            //cerr << "return genome"<<endl;
            //cerr << "************FINAL******************************************" <<endl;
            // genome->print_phased_genome();
            return std::move(genome); 
        }else{
            // cerr << "return optimal"<<endl;
            // cerr << "************FINAL******************************************" <<endl;
            // optimal->print_phased_genome();
            return std::move(optimal); 
        }

    }   
    double MCMCGenotyper::log_target(PhasedGenome& phased_genome, const vector<MultipathAlignment>& reads)const{
        
        // sum of scores given the reads aligned on the haplotype 
        int32_t sum_scores = 0; 
        
        // get scores for mp alignments 
        for(MultipathAlignment mp : reads){
            identify_start_subpaths(mp);
            sum_scores += phased_genome.optimal_score_on_genome(mp, graph);
            
        } 
        
        return sum_scores;
    }

    tuple<int, const Snarl*, vector<NodeTraversal> > MCMCGenotyper::proposal_sample(PhasedGenome& current)const{
        // get a different traversal through the snarl by uniformly choosing from all possible ways to traverse the snarl
        
        // bookkeeping: haplotype ID, snarl* (site that we replaced at), get_allele())
        tuple<int, const Snarl*, vector<NodeTraversal> > to_return;

        int& random_haplotype = get<0>(to_return);
        const Snarl*& random_snarl = get<1>(to_return);
        // the random path through the snarl 
        vector<NodeTraversal>& old_allele = get<2>(to_return);
        
        // sample uniformly between snarls 
        random_snarl = snarls.discrete_uniform_sample(random_engine);
        
        // comment this out when testing Test1 - will cause segfault
        // cerr << "this is the boundary nodes of the snarl " << random_snarl->start() << random_snarl->end() <<endl;
        

        if(random_snarl == nullptr){
            cerr << "random_snarl is null " <<endl;
            random_haplotype = -1;
            return to_return;
        }


        // get list of haplotypes that contain the snarl
        vector<id_t> matched_haplotypes = current.get_haplotypes_with_snarl(random_snarl);


        if(matched_haplotypes.empty()){
            cerr<< "matched haplotype is empty"<< endl;
            // cerr << "looking for snarl starting with " << random_snarl->start() << " and snarl ending with " << random_snarl->end() <<endl;
            random_haplotype = -1;
            return to_return;
            
        }

        // choose a haplotype uiformly 
        id_t lower_bound = 0;
        id_t upper_bound = matched_haplotypes.size()-1;

        int random_num = generate_discrete_uniform(random_engine, lower_bound, upper_bound);
        random_haplotype = matched_haplotypes[random_num];


        pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarls.deep_contents(random_snarl, graph, true);

        // unpack the pair, we only care about the node_ids
        unordered_set<id_t>& ids = contents.first;
            
        // enumerate counts through nodes in snarl not the entire graph
        SubHandleGraph subgraph(&graph);
        
        for (id_t id : ids){
            
            // add each node from snarl in super graph to sub graph
            subgraph.add_handle(graph.get_handle(id, false));
        }
        
        
        // create a count_map of the subgraph
        auto count_contents = algorithms::count_walks_through_nodes(&subgraph);
        
        // unpack the count map from the count_contents 
        unordered_map<handle_t, size_t>& count_map = get<1>(count_contents);

        
        
        // create a topological order of sub graph count map
        vector<handle_t> topological_order = algorithms::lazier_topological_order(&subgraph);

        //  we want to get just the sink handle handle
        handle_t start = topological_order.back();  
        handle_t source = topological_order.front();
        
        // start at sink in topological
        bool start_from_sink =true;
        bool not_source = true;

        vector<NodeTraversal> allele;
 
        while(not_source){

            size_t  cum_sum = 0;
            vector<size_t> cumulative_sum;
            vector<size_t> paths_to_take;
            size_t  count = 0;
            vector<handle_t> handles;
            
            subgraph.follow_edges(start, start_from_sink, [&](const handle_t& next) { 
                unordered_map<handle_t, size_t>::iterator it;
                it = count_map.find(next); // find the handle
                count = it->second; // get the count 
                cum_sum += count;
                cumulative_sum.push_back(cum_sum); 
                handles.push_back(next); 
        
            });

            // choose a random path uniformly
            int l_bound = 0;
            int u_bound = cumulative_sum.back()-1;
            int random_num = generate_discrete_uniform(random_engine,l_bound, u_bound);

            // use the random_num to select a random_handle
            int found = 0, prev = 0;
            for (int i = 0; i< cumulative_sum.size() ; i++) {
                // check what range the random_num falls in    
                if (prev <= random_num && random_num < cumulative_sum[i] ){
                    found = i; // will correspond to the index of handles
                    break; 
                }
                prev = cumulative_sum[i];
            } 

            assert(found != -1);

            // start_ptr will point to random handle 
            start = handles[found]; 
            
            // save the random path 
            bool position = subgraph.get_is_reverse(start);
            Node* n = graph.get_node(subgraph.get_id(start));


            // allele should not include boundary nodes of random_snarl
            if(n->id() != random_snarl->start().node_id() && n->id() != random_snarl->end().node_id() ){
                allele.push_back(NodeTraversal(n,position));
            }
            
            

            // check if we are at the source, if so we terminate loop
            if(start == source){
                not_source = false;
            }    
            
        }
        // save old allele so we can swap back to it if we need to 
        old_allele = current.get_allele(*random_snarl, random_haplotype);
        
        // set new allele with random allele, replace with previous allele 
        current.set_allele(*random_snarl , allele.rbegin(), allele.rend(), random_haplotype);
        
 
        return to_return;

    }
    int MCMCGenotyper::generate_discrete_uniform(minstd_rand0& random_engine, int lower_bound , int upper_bound) const{
        
        // choose a number randomly using discrete uniform distribution
        uniform_int_distribution<int> distribution(lower_bound, upper_bound);  
        int random_num = distribution(random_engine);

        return random_num;
    }
    double MCMCGenotyper::generate_continuous_uniform(const double a, const double b)const{
        
        uniform_real_distribution<double> distribution(a,b);
        double random_num = distribution(random_engine);
        
        return random_num;

    }
    unique_ptr<PhasedGenome> MCMCGenotyper::generate_initial_guess()const{
        
        unique_ptr<PhasedGenome> genome(new PhasedGenome(snarls));
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
        genome->add_haplotype(haplotype.begin(), haplotype.end());
        genome->add_haplotype(haplotype.begin(), haplotype.end());

        // index sites
        genome->build_indices();

        return std::move(genome);
    }

}



