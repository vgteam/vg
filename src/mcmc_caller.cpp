#include "mcmc_caller.hpp"  
#include "graph_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "memoizing_graph.hpp"


namespace vg {

    /**
     * MCMCCaller : Inherits from VCFOutputCaller    
     */         
    MCMCCaller::MCMCCaller(const PathHandleGraph& graph,
                           SnarlManager& snarl_manager,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<size_t>& ref_path_offsets,
                           const vector<size_t>& ref_path_lengths,
                           ostream& out_stream) :
        graph(graph), snarl_manager(snarl_manager), sample_name(sample_name), ref_paths(ref_paths), ref_path_offsets(ref_path_offsets),
        ref_path_lengths(ref_path_lengths), out_stream(out_stream) {
        

    }
    
    MCMCCaller:: ~MCMCCaller(){

    }


    void process_variant(const Snarl& snarl, TraversalFinder& trav_finder, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, const string& ref_path_name) {
            
        // convert traversal to string
        // function<string(const SnarlTraversal&)> trav_string = [&](const SnarlTraversal& trav) {
        // string seq;
        // for (int i = 0; i < trav.visit_size(); ++i) {
        //     seq += graph->get_sequence(graph->get_handle(trav.visit(i).node_id(), trav.visit(i).backward()));
        // }
        // return seq;
        // };
        
        vcflib::Variant out_variant;

        // fill out the rest of the variant
        out_variant.sequenceName = ref_path_name;
        // +1 to convert to 1-based VCF
        // out_variant.position = get_ref_position(snarl, ref_path_name).first + ref_offsets.find(ref_path_name)->second + 1; 
        out_variant.id = std::to_string(snarl.start().node_id()) + "_" + std::to_string(snarl.end().node_id());
        out_variant.filter = "PASS";
        out_variant.updateAlleleIndexes();


        
        // convert to alleles
        // ref path is always allele 0  
        
        // clean up the alleles to not have so man common prefixes
        // flatten_common_allele_ends(out_variant, true);
        // flatten_common_allele_ends(out_variant, false);

        // if (!out_variant.alt.empty()) {
        //     add_variant(out_variant);
        // }

    }



    

    


}  





