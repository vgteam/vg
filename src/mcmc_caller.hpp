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
#include "traversal_finder.hpp"
#include "region.hpp"

namespace vg{

    using namespace std;

    class MCMCCaller : public VCFOutputCaller{
        
        const PathPositionHandleGraph& graph;
        SnarlManager& snarl_manager;
        const string& sample_name;
        const vector<string>& ref_paths = {};
        const vector<size_t>& ref_path_offsets = {};
        ostream& out_stream; 
    
    
    
    public:
        MCMCCaller(const PathPositionHandleGraph& graph,
                    SnarlManager& snarl_manager,
                    const string& sample_name,
                    const vector<string>& ref_paths = {},
                    const vector<size_t>& ref_path_offsets = {},
                    ostream& out_stream = cout );

        virtual ~MCMCCaller(); 

        protected:

        /// print a vcf variant 
        void process_variant(const Snarl& snarl, TraversalFinder& trav_finder, const vector<SnarlTraversal>& called_traversals,
                      const vector<int>& genotype, const string& ref_path_name) const;


        /// Our Genotyper
        //TODO: find out how this is being used to possibly substitute or remove completely (SnarlCaller)

        /// Our snarls
        SnarlManager& snarl_manager;

        /// Where all output written
        ostream& out_stream;

        /// the graph
        const PathPositionHandleGraph& graph;
        /// non-vg inputs are converted into vg as-needed, at least until we get the
        /// traversal finding ported
        bool is_vg;

        /// The old vg call traversal finder.  It is fairly efficient but daunting to maintain.
        /// We keep it around until a better replacement is implemented.  It is *not* compatible
        /// with the Handle Graph API because it relise on PathIndex.  We convert to VG as
        /// needed in order to use it. 
        RepresentativeTraversalFinder* traversal_finder;
        /// Needed by above (only used when working on vg inputs -- generated on the fly otherwise)
        vector<PathIndex*> path_indexes;

        /// keep track of the reference paths
        vector<string> ref_paths;

        /// keep track of offsets in the reference paths
        map<string, size_t> ref_offsets;

    };

    


}
#endif VG_MCMC_CALLER_HPP_INCLUDED