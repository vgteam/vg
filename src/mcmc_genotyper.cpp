#include "mcmc_genotyper.hpp"
#include "algorithms/count_walks.hpp"
#include "subgraph.hpp"
#include <algorithm> 
#include <random> 
#include <chrono>
#include <utility>
#include "multipath_alignment.hpp"

// #define debug_mcmc

namespace vg {

    using namespace std;

    MCMCGenotyper::MCMCGenotyper(SnarlManager& snarls, VG& graph, const int n_iterations, const int seed):snarls(snarls), graph(graph), n_iterations(n_iterations), 
        seed(seed), random_engine(seed){
            
    
    }

    unique_ptr<PhasedGenome> MCMCGenotyper::run_genotype(const vector<multipath_alignment_t>& reads, const double log_base) const{

        // set a flag for invalid contents so a message is observed 
        bool invalid_contents = false;
        bool return_optimal = false;

        // generate initial value
        unique_ptr<PhasedGenome> genome = generate_initial_guess();
        
        double max_likelihood = 0.0; 
        double current_likelihood = 0.0; 
        double previous_likelihood = 0.0;
 
        unique_ptr<PhasedGenome> optimal;
        
        // build markov chain using Metropolis-Hastings
        for(int i = 0; i< n_iterations; i++){

            
            // holds the previous sample allele
            double x_prev = log_target(genome, reads);

            // get contents from proposal_sample
            tuple<int, const Snarl*, vector<NodeTraversal> > to_receive = proposal_sample(genome);
            int& modified_haplo = get<0>(to_receive);         
            
            // if the to_receive contents are invalid keep the new allele
            // for graphs that do not contain snarls
            if (modified_haplo ==-1){
                invalid_contents = true;
                continue;
            }else{
                const Snarl& modified_site = *get<1>(to_receive); 
                vector<NodeTraversal>& old_allele = get<2>(to_receive); 
 
                // holds new sample allele score 
                double x_new = log_target(genome, reads);

                double likelihood_ratio = exp(log_base*(x_new - x_prev));
                
                current_likelihood = previous_likelihood + log_base*(x_new-x_prev);
                
                if (current_likelihood > max_likelihood){
                    max_likelihood = current_likelihood;
                    optimal = unique_ptr<PhasedGenome>(new PhasedGenome(*genome));
                    return_optimal=true;
                }
                
                // calculate acceptance probability 
                double acceptance_probability = min(1.0, likelihood_ratio);

                // if u~U(0,1) > alpha, discard new allele and keep previous 
                auto uniform_smpl = generate_continuous_uniform(0.0,1.0);
                if(uniform_smpl > acceptance_probability){ 
                    genome->set_allele(modified_site, old_allele.begin(), old_allele.end(), modified_haplo); 
#ifdef debug_mcmc
                    cerr << "Rejected new allele" <<endl;
                    cerr << "clikelihood " << previous_likelihood <<endl;
                    genome->print_phased_genome();
#endif                    
                }else{     
#ifdef debug_mcmc 
                    cerr << "Accepted new allele" <<endl;
                    cerr << "clikelihood " << current_likelihood <<endl;
                    genome->print_phased_genome();
#endif
                    previous_likelihood = current_likelihood;
                    
                }         
            }
        } 
        if(invalid_contents || !return_optimal){
            // for graphs without snarls 
            return genome; 
        }else{
#ifdef debug_mcmc 
            cerr <<"klikelihood " << max_likelihood <<endl;
            optimal->print_phased_genome();
#endif
            return optimal; 
        }

    }   
    double MCMCGenotyper::log_target(unique_ptr<PhasedGenome>& phased_genome, const vector<multipath_alignment_t>& reads)const{
        
        // sum of scores given the reads aligned on the haplotype 
        int32_t sum_scores = 0; 
        
        // get scores for mp alignments 
        for(const multipath_alignment_t& mp : reads){
            sum_scores += phased_genome->optimal_score_on_genome(mp, graph);
            
        } 
        
        return sum_scores;
    }

    tuple<int, const Snarl*, vector<NodeTraversal> > MCMCGenotyper::proposal_sample(unique_ptr<PhasedGenome>& current)const{
        // get a different traversal through the snarl by uniformly choosing from all possible ways to traverse the snarl
        
        // bookkeeping: haplotype ID, snarl* (site that we replaced at), get_allele())
        tuple<int, const Snarl*, vector<NodeTraversal> > to_return;

        int& random_haplotype = get<0>(to_return);
        const Snarl*& random_snarl = get<1>(to_return);
        // the random path through the snarl 
        vector<NodeTraversal>& old_allele = get<2>(to_return);
        
        // sample uniformly between snarls 
        random_snarl = snarls.discrete_uniform_sample(random_engine);

        if(random_snarl == nullptr){
            random_haplotype = -1;
            return to_return;
        }


        // get list of haplotypes that contain the snarl
        vector<id_t> matched_haplotypes = current->get_haplotypes_with_snarl(random_snarl);


        if(matched_haplotypes.empty()){
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
        old_allele = current->get_allele(*random_snarl, random_haplotype);

#ifdef debug_mcmc 
        cerr << "modifying haplotype " << random_num << endl; 
        for(auto iter = allele.begin(); iter != allele.end(); iter++ ){
            cerr << "new allele: " <<"node " << iter->node->id() << " " << iter->node->sequence() <<endl;

        }
        for(auto iter = old_allele.begin(); iter != old_allele.end(); iter++ ){
            cerr << "old allele: " <<"node " << iter->node->id() <<  " " << iter->node->sequence() <<endl;

        }
#endif
        // set new allele with random allele, replace with previous allele 
        current->set_allele(*random_snarl , allele.rbegin(), allele.rend(), random_haplotype);
        
        
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

        return genome;
    }
    unordered_map<pair<const Snarl*, const Snarl*>, int32_t> MCMCGenotyper::make_snarl_map(const vector<multipath_alignment_t>& reads, unique_ptr<PhasedGenome>& phased_genome) const{
        
        unordered_set<const Snarl*> snarl_set;
        unordered_map<pair<const Snarl*, const Snarl*>, int32_t> map;
        vector<pair<const Snarl*, const Snarl*>> pairs;
        int32_t score_after_swap,score_before_swap,diff_score;

        //loop over reads
        for(const multipath_alignment_t& multipath_aln : reads){
            //for each pair of snarls that touches that read
            for (const auto& subpath : multipath_aln.subpath()) {
                
                if(subpath.has_path()){
                    auto& path = subpath.path();
                    //for every mapping in the path
                    for(size_t i = 0; i < path.mapping_size(); i++){
                        auto& mapping = path.mapping(i);
                        int64_t node_id = mapping.position().node_id();
                        const Snarl* back_snarl = snarls.into_which_snarl(node_id, true); 
                        const Snarl* fwd_snarl = snarls.into_which_snarl(node_id, false);
                        
                        // insert snarls into unordered set with only unique entries 
                        snarl_set.insert(back_snarl);
                        snarl_set.insert(fwd_snarl);
                    }
                    
                }
                
            }

            //get_optimal_score_on_genome(genome_before_swap, read)
            score_before_swap = phased_genome->optimal_score_on_genome(multipath_aln, graph);
            
            //for each pair of snarls - check if both haplotypes visited these snarls, 
            vector<const Snarl*> v(snarl_set.begin(), snarl_set.end());
            for(int i =0; i < v.size(); i++){
                for(int j =i+1; j < v.size(); j++){
                    pairs.push_back(make_pair(v[i], v[j]));

                }
            }

            for(auto snarl_ptr:pairs){
                // TODO: update to consider cyclic paths
                // check that both haplotypes visit each snarl in the snarl pair
                vector<id_t> haplo_ids1 = phased_genome->get_haplotypes_with_snarl(snarl_ptr.first);
                vector<id_t> haplo_ids2 = phased_genome->get_haplotypes_with_snarl(snarl_ptr.second);
                //if so:
                if(haplo_ids1.size()==2 && haplo_ids2.size()==2){
                    int haplotype_1 =0;
                    int haplotype_2 =1;
                    //generate a random uniform number between [0,1]
                    int random_num = generate_discrete_uniform(random_engine, haplotype_1, haplotype_2);
                    //exchange their alleles with each other at one of the snarls (chosen randomly)
                    //TODO: for < 2 or > 2 haplotypes that overlap snarl pair , skip snarl pair for that read 
                    if(random_num == 1){
                        //dereference the ptr
                        const Snarl& snarl_to_swap = *snarl_ptr.first;
                        //exhange alleles at first snarl in pair
                        phased_genome->swap_alleles(snarl_to_swap, haplotype_1, haplotype_2);
                        // get score after swap
                        score_after_swap = phased_genome->optimal_score_on_genome(multipath_aln, graph);
                        //swap back 
                        phased_genome->swap_alleles(snarl_to_swap, haplotype_1, haplotype_2);
                        
                    }else{
                        //dereference the ptr
                        const Snarl& snarl_to_swap = *snarl_ptr.second;
                        //exchange alleles at second snarl in pair
                        phased_genome->swap_alleles(snarl_to_swap, haplotype_1, haplotype_2);
                        // get score after swap
                        score_after_swap = phased_genome->optimal_score_on_genome(multipath_aln, graph);
                        //swap back 
                        phased_genome->swap_alleles(snarl_to_swap, haplotype_1, haplotype_2);
                    }
                    
                    //getcalculate difference of scores between swaps
                    diff_score = score_before_swap - score_after_swap;

                    if(score_before_swap > score_after_swap){
                        map[make_pair(snarl_ptr.first, snarl_ptr.second)] += diff_score;
                    }else{
                        map[make_pair(snarl_ptr.first, snarl_ptr.second)] -= diff_score;
                    }
                    
                    
                }
                    
               
            }
            
        } 
        return  map;
    }

    algorithms::Graph MCMCGenotyper::make_snarl_graph(unordered_map<pair<const Snarl*, const Snarl*>, int32_t> map) const{
        //TODO: find where the SnarlRecord* are being added to deque and store the index in snarls.cpp
        
        algorithms::Graph snarl_graph;

        for(auto snarl_pair_to_weight: map){
            size_t edge_weight = (size_t)snarl_pair_to_weight.second;
            pair<const Snarl*, const Snarl*> snarl_pair = snarl_pair_to_weight.first;
            const Snarl* snarl_1 = snarl_pair.first;
            const Snarl* snarl_2 = snarl_pair.second;
            // skip edge weights that are <1
            if(edge_weight < 1){
                continue;
            }else{
                algorithms::Node snarl_node_1, snarl_node_2;
                algorithms::Edge edge_fwd, edge_back;
                edge_fwd.weight = edge_weight;
                edge_back.weight = edge_weight;
                snarl_node_1.edges.push_back(edge_fwd);
                snarl_node_2.edges.push_back(edge_back);
                snarl_node_1.weight += edge_weight;
                snarl_node_2.weight += edge_weight;

                //node ids, get the ids from SnarlRecord index member 
                size_t snarl_id_1 = snarls.snarl_number(snarl_1);
                size_t snarl_id_2 = snarls.snarl_number(snarl_2);

                edge_fwd.other = snarl_id_2;
                edge_back.other = snarl_id_1;
                snarl_graph.add_node(snarl_id_1, snarl_node_1);
                snarl_graph.add_node(snarl_id_2, snarl_node_2);

            }

        }
        
        return snarl_graph; 

    }


}



