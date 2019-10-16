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
#include "region.hpp"

namespace vg{

    using namespace std;

    class MCMCCaller : public VCFOutputCaller {
    public:    
        const PathHandleGraph& graph;
        SnarlManager& snarl_manager;
        const string& sample_name;
        const vector<string>& ref_paths = {};
        const vector<size_t>& ref_path_offsets = {};
        const vector<size_t>& ref_path_lengths = {};
        ostream& out_stream; 
    
    
    
    
        MCMCCaller(const PathHandleGraph& graph,
                    SnarlManager& snarl_manager,
                    const string& sample_name,
                    const vector<string>& ref_paths = {},
                    const vector<size_t>& ref_path_offsets = {},
                    const vector<size_t>& ref_path_lengths = {},
                    ostream& out_stream = cout );

        virtual ~MCMCCaller(); 

        protected:

        /// print a vcf variant 
        void process_variant(const Snarl& snarl, TraversalFinder& trav_finder, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, const string& ref_path_name) ;


        /// Our Genotyper
        //TODO: find out how this is being used to possibly substitute or remove completely (SnarlCaller)


    };

    


}
#endif 