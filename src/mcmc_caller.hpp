#ifndef VG_MCMC_CALLER_HPP_INCLUDED
#define VG_MCMC_CALLER_HPP_INCLUDED


#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
#include "snarls.hpp"
#include "graph_caller.hpp"
#include "traversal_finder.hpp"
#include "phased_genome.hpp"
#include "region.hpp"

namespace vg{

    using namespace std;

    class MCMCCaller : public VCFOutputCaller {
    public:    
        PhasedGenome& genome;
        SnarlManager& snarl_manager;
        const string sample_name = "SAMPLE";
        const vector<size_t> ref_path_offsets = {};
        const vector<size_t> ref_path_lengths = {};
        ostream& out_stream; 
        const SnarlTraversal trav;
        

        MCMCCaller(const PathPositionHandleGraph* path_position_handle_graph,
                    PhasedGenome& genome,
                    SnarlManager& snarl_manager,
                    const string& sample_name,
                    const vector<string>& ref_paths,
                    const vector<size_t>& ref_path_offsets,
                    const vector<size_t>& ref_path_lengths,
                    ostream& out_stream = cout );

        virtual ~MCMCCaller(); 
        
        /// Run call_snarl() on every top-level snarl in the manager.
        /// For any that return false, try the children, etc. (when recurse_on_fail true)
        /// Snarls are processed in parallel
        void call_top_level_snarls(bool recurse_on_fail = true) ;

        /// print vcf header
        virtual string vcf_header(const PathPositionHandleGraph& graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const ;
    
    protected:  
        /// path position handle graph
        const PathPositionHandleGraph* path_position_handle_graph;

        /// keep track of the reference paths
        vector<string> ref_paths;

        /// keep track of offsets in the reference paths
        map<string, size_t> ref_offsets; 

        /// Update INFO and FORMAT fields of the called variant
        void update_vcf_info(const Snarl& snarl,
                                 const vector<SnarlTraversal>& traversals,
                                 const vector<int>& genotype,
                                 const string& sample_name,
                                 vcflib::Variant& variant) const;

        /// print a vcf variant 
        void emit_variant(const Snarl& snarl, const vector<int>& genotype, SnarlTraversal ref_trav, const string& ref_path_name, const vector<SnarlTraversal>& haplo_travs) const;

        /// Call a given snarl, and print the output to out_stream
        bool call_snarl(const Snarl& snarl);

        /// check if a site can be handled 
        bool is_traversable(const Snarl& snarl);

        /// get position of reference path 
        pair<size_t, bool> get_ref_position(const Snarl& snarl, const string& ref_path_name) const;

        /// clean up the alleles to not share common prefixes / suffixes
        void flatten_common_allele_ends(vcflib::Variant& variant, bool backward) const;


    };

    


}
#endif 
