#include "snarls.hpp"
#include "snarl_graph.hpp"
#include <structures/union_find.hpp>
#include "sparse_union_find.hpp"

namespace vg{
    using namespace std;
    using vg::algorithms::Graph;

class SnarlMap{
       
public:
    SnarlMap(SnarlManager& snarls, const vector<multipath_alignment_t>& reads, unique_ptr<PhasedGenome>& phased_genome): 
    snarls(snarls), reads(reads), phased_genome(phased_genome){

    }

    unordered_map<size_t, size_t> id_weight> make_snarl_map(){
        
    //loop over reads
        //for each pair of snarls that touches that read
            //get_optimal_score_on_genome(genome, read)
            //choose random haplotype?
            //swap_allele_on_haplotype_at_given_snarl()
            //get_optimal_score_on_genome(genome_with_swap, read)
            //get_phased_genome_before_swap()
            //get_difference_of_scores()
            //map.emplace(snarl_id, weight);
            /*Notes: 
            if score decreases after swap, contribute (+) score to edge weight 
            For negative or 0 weights, filter them out of the graph
            **/

    }

    Graph make_snarl_graph(unordred_map<pair<size_t, size_t>, size_t> map){

        
    }

    

    
        
};


}