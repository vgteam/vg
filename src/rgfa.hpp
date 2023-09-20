#ifndef VG_RGFA_HPP_INCLUDED
#define VG_RGFA_HPP_INCLUDED

/**
 * \file rgfa.hpp
 *
 * Interface for computing and querying rGFA path covers.
 *
 * An rGFA cover is a set of path fragments (stored as separate paths) in the graph.
 * They are always relative to an existing reference sample (ie GRCh38).
 * rGFA path fragments have a special name format, where sample=rGFA
 *
 * The data structures used in this class are always relative to the original paths
 * in the graph.  The REfERENCE-SENSE fragments that are used to serialize the
 * cover (and for surjection) can be created and loaded, but they are not used
 * beyond that.  
 */

#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"

namespace vg {

using namespace std;

class RGFACover {
public:
    // I'm not sure how much we want to allow this to change.  For now just define it here
    static const string rgfa_sample_name;

    // make an rGFA path name from a subrange of a normal path
    static string make_rgfa_path_name(const string& path_name, int64_t start, int64_t length,
                                      bool specify_subrange_end = true);

    // test if path name is rGFA
    static bool is_rgfa_path_name(const string& path_name);

    // break up the rGFA locus name back into original sample and locus
    static pair<string, string> parse_rgfa_locus_name(const string& locus_name);

public:
    // clear out the existing rGFA cover from the graph.  recommended to run this
    // before compute()
    void clear(MutablePathMutableHandleGraph* graph);

    // compute the rgfa cover from the graph, starting with a given set of reference paths
    void compute(const PathHandleGraph* graph,
                 SnarlManager* snarl_manager,
                 const unordered_set<path_handle_t>& reference_paths,
                 int64_t minimum_length);

    // load the rgfa cover from the graph, assuming it's been computed already and
    // saved in special rgfa paths
    // this function assumes that the rGFA cover paths are exactly consistent
    // with the original paths.
    void load(const PathHandleGraph* graph,
              const unordered_set<path_handle_t>& reference_paths);

    // apply the rgfa cover to a graph (must have been computed first), adding it
    // as a bunch of (reference sense) paths with a special sample name
    void apply(MutablePathMutableHandleGraph* mutable_graph);

    // get the rgfa rank (level) of a given node (0 if on a reference path)
    int64_t get_rank(nid_t node_id) const;

    // get the rgfa step of a given node
    step_handle_t get_step(nid_t node_id) const;

    // get the parent intervals
    pair<const pair<step_handle_t, step_handle_t>*,
         const pair<step_handle_t, step_handle_t>*> get_parent_intervals(const pair<step_handle_t, step_handle_t>& interval) const;

    // get the intervals
    const vector<pair<step_handle_t, step_handle_t>>& get_intervals() const;

    // get an interval from a node
    // return nullptr if node not in an interval
    const pair<step_handle_t, step_handle_t>* get_interval(nid_t node_id) const;

    
protected:

    // compute the cover for the given snarl, by greedily finding the convered rgfa paths through it
    // the cover is added to the two "thread_" structures.
    void compute_snarl(const Snarl& snarl, PathTraversalFinder& path_trav_finder, int64_t minimum_length,
                       vector<pair<step_handle_t, step_handle_t>>& thread_rgfa_intervals,
                       unordered_map<nid_t, int64_t>& thread_node_to_interval);

    // get intervals in travs that are not covered according to this->node_to_interval or
    // the  thread_node_to_interval parameter
    vector<pair<int64_t, int64_t>> get_uncovered_intervals(const vector<step_handle_t>& trav,
                                                           const unordered_map<nid_t, int64_t>& thread_node_to_interval);

    // get the total coverage of a traversal (sum of step lengths)
    int64_t get_coverage(const vector<step_handle_t>& trav, const pair<int64_t, int64_t>& uncovered_interval);

    // make sure all nodes in all rGFA paths are in forward orientation
    // this is always possible because they are, by definition, disjoint
    // this should only be run from inside apply()
    void forwardize_rgfa_paths(MutablePathMutableHandleGraph* mutable_graph);
    
    
protected:

    const PathHandleGraph* graph;

    // intervals are end-exclusive (like bed)
    vector<pair<step_handle_t, step_handle_t>> rgfa_intervals;

    // so rgfa_intervals[0,num_ref_intervals-1] will be all rank-0 reference intervals
    int64_t num_ref_intervals;

    unordered_map<nid_t, int64_t> node_to_interval;
};

/// Export the given VG graph to the given GFA file.
}

#endif
