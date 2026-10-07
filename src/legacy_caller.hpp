#ifndef VG_LEGACY_CALLER_HPP_INCLUDED
#define VG_LEGACY_CALLER_HPP_INCLUDED

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
#include "snarl_caller.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"
#include "graph_caller.hpp"
#include "vcf_output_caller.hpp"

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/**
 * LegacyCaller : Preserves (most of) the old vg call logic by using 
 * the RepresentativeTraversalFinder to recursively find traversals
 * through arbitrary sites.   
 */
class LegacyCaller : public GraphCaller, public VCFOutputCaller {
public:
    LegacyCaller(const PathPositionHandleGraph& graph,
                 SupportBasedSnarlCaller& snarl_caller,
                 SnarlManager& snarl_manager,
                 const string& sample_name,
                 const vector<string>& ref_paths = {},
                 const vector<size_t>& ref_path_offsets = {},
                 const vector<int>& ref_path_ploidies = {});

    virtual ~LegacyCaller();

    virtual bool call_snarl(const Snarl& snarl);

    virtual string vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                              const vector<size_t>& contig_length_overrides = {}) const;

protected:

    /// recursively genotype a snarl
    /// todo: can this be pushed to a more generic class? 
    pair<vector<SnarlTraversal>, vector<int>> top_down_genotype(const Snarl& snarl, TraversalFinder& trav_finder, int ploidy,
                                                                const string& ref_path_name, pair<size_t, size_t> ref_interval) const;
    
    /// we need the reference traversal for VCF, but if the ref is not called, the above method won't find it. 
    SnarlTraversal get_reference_traversal(const Snarl& snarl, TraversalFinder& trav_finder) const;

    /// re-genotype output of top_down_genotype.  it may give slightly different results as
    /// it's working with fully-defined traversals and can exactly determine lengths and supports
    /// it will also make sure the reference traversal is in the beginning of the output
    tuple<vector<SnarlTraversal>, vector<int>, unique_ptr<SnarlCaller::CallInfo>> re_genotype(const Snarl& snarl,
                                                                                              TraversalFinder& trav_finder,
                                                                                              const vector<SnarlTraversal>& in_traversals,
                                                                                              const vector<int>& in_genotype,
                                                                                              int ploidy,
                                                                                              const string& ref_path_name,
                                                                                              pair<size_t, size_t> ref_interval) const;

    /// check if a site can be handled by the RepresentativeTraversalFinder
    bool is_traversable(const Snarl& snarl);

    /// look up a path index for a site and return its name too
    pair<string, PathIndex*> find_index(const Snarl& snarl, const vector<PathIndex*> path_indexes) const;

protected:

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

    /// keep track of ploidies in the reference paths
    map<string, int> ref_ploidies;

    /// Tuning

    /// How many nodes should we be willing to look at on our path back to the
    /// primary path? Keep in mind we need to look at all valid paths (and all
    /// combinations thereof) until we find a valid pair.
    int max_search_depth = 1000;
    /// How many search states should we allow on the DFS stack when searching
    /// for traversals?
    int max_search_width = 1000;
    /// What's the maximum number of bubble path combinations we can explore
    /// while finding one with maximum support?
    size_t max_bubble_paths = 100;

};

}

#endif
