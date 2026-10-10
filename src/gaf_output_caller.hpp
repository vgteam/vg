#ifndef VG_GAF_OUTPUT_CALLER_HPP_INCLUDED
#define VG_GAF_OUTPUT_CALLER_HPP_INCLUDED

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

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/**
 * Helper class for outputing snarl traversals as GAF
 */
class GAFOutputCaller {
public:
    /// The emitter object is created and owned by external forces
    GAFOutputCaller(AlignmentEmitter* emitter, const string& sample_name, const vector<string>& ref_paths,
                    size_t trav_padding);
    virtual ~GAFOutputCaller();

    /// print the GAF traversals
    void emit_gaf_traversals(const PathHandleGraph& graph, const string& snarl_name,
                             const vector<SnarlTraversal>& travs,
                             int64_t ref_trav_idx,
                             const string& ref_path_name, int64_t ref_path_position,
                             const TraversalSupportFinder* support_finder = nullptr);

    /// print the GAF genotype
    void emit_gaf_variant(const PathHandleGraph& graph, const string& snarl_name,
                          const vector<SnarlTraversal>& travs,
                          const vector<int>& genotype,
                          int64_t ref_trav_idx,
                          const string& ref_path_name, int64_t ref_path_position,
                          const TraversalSupportFinder* support_finder = nullptr);
    
    /// pad a traversal with (first found) reference path, adding up to trav_padding to each side
    SnarlTraversal pad_traversal(const PathHandleGraph& graph, const SnarlTraversal& trav) const;
    
protected:
    
    AlignmentEmitter* emitter;

    /// Sample name
    string gaf_sample_name;

    /// Add padding from reference paths to traversals to make them at least this long
    /// (only in emit_gaf_traversals(), not emit_gaf_variant)
    size_t trav_padding = 0;

    /// Reference paths are used to pad out traversals.  If there are none, then first path found is used
    unordered_set<string> ref_paths;

};

}

#endif
