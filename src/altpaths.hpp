#ifndef VG_ALTPATHS_HPP_INCLUDED
#define VG_ALTPATHS_HPP_INCLUDED

/**
 * \file altpaths.hpp
 *
 * Interface for computing and querying alternative path covers.
 *
 * An altpath cover is a set of path fragments (stored as separate paths) in the graph.
 * They are always relative to an existing reference sample (ie GRCh38 or CHM13).
 * Unlike rGFA paths which use complex metadata embedding, altpaths use a simple naming
 * scheme: {base_path_name}_alt{N}
 *
 * For example, if the reference path is "CHM13#0#chr1", altpaths would be named:
 *   - CHM13#0#chr1_alt1
 *   - CHM13#0#chr1_alt2
 *   - etc.
 *
 * The data structures used in this class are always relative to the original paths
 * in the graph. The REFERENCE-sense fragments that are used to serialize the
 * cover can be created and loaded, but they are not used beyond that.
 */

#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"

namespace vg {

using namespace std;

class AltPathsCover {
public:
    // The suffix used to identify altpaths
    static const string altpath_suffix;  // "_alt"

    // Create an altpath name from a base reference path name and an index.
    // Example: make_altpath_name("CHM13#0#chr1", 1) -> "CHM13#0#chr1_alt1"
    static string make_altpath_name(const string& base_path_name, int64_t alt_index);

    // Test if a path name is an altpath (contains "_alt" suffix followed by digits).
    static bool is_altpath_name(const string& path_name);

    // Parse an altpath name to extract the base reference path name.
    // Returns the original base path name, or the input if not an altpath.
    // Example: parse_base_path("CHM13#0#chr1_alt3") -> "CHM13#0#chr1"
    static string parse_base_path(const string& altpath_name);

    // Parse an altpath name to extract the alt index.
    // Returns -1 if the path is not an altpath.
    // Example: parse_alt_index("CHM13#0#chr1_alt3") -> 3
    static int64_t parse_alt_index(const string& altpath_name);

public:
    // Clear out any existing altpaths from the graph. Recommended to run this
    // before compute().
    void clear(MutablePathMutableHandleGraph* graph);

    // Compute the altpath cover from the graph, starting with a given set of reference paths.
    void compute(const PathHandleGraph* graph,
                 SnarlManager* snarl_manager,
                 const unordered_set<path_handle_t>& reference_paths,
                 int64_t minimum_length);

    // Load existing altpaths from the graph, assuming they've been computed already.
    // The reference_paths should be the rank-0 paths the altpaths extend from.
    void load(const PathHandleGraph* graph,
              const unordered_set<path_handle_t>& reference_paths);

    // Apply the altpath cover to a graph (must have been computed first), adding it
    // as a bunch of REFERENCE-sense paths with the simplified naming scheme.
    void apply(MutablePathMutableHandleGraph* mutable_graph);

    // Get the rank (level) of a given node (0 if on a reference path).
    int64_t get_rank(nid_t node_id) const;

    // Get the step of a given node in its covering interval.
    step_handle_t get_step(nid_t node_id) const;

    // Get the parent intervals (left and right) of a given interval.
    pair<const pair<step_handle_t, step_handle_t>*,
         const pair<step_handle_t, step_handle_t>*> get_parent_intervals(const pair<step_handle_t, step_handle_t>& interval) const;

    // Get all computed intervals.
    const vector<pair<step_handle_t, step_handle_t>>& get_intervals() const;

    // Get an interval from a node. Returns nullptr if node not in an interval.
    const pair<step_handle_t, step_handle_t>* get_interval(nid_t node_id) const;

    // Get the number of reference intervals (rank-0).
    int64_t get_num_ref_intervals() const;

    // Print out a table of statistics.
    void print_stats(ostream& os);

protected:

    // Compute the cover for the given snarl, by greedily finding the covered paths through it.
    // The cover is added to the two "thread_" structures.
    void compute_snarl(const Snarl& snarl, PathTraversalFinder& path_trav_finder, int64_t minimum_length,
                       vector<pair<step_handle_t, step_handle_t>>& thread_altpath_intervals,
                       unordered_map<nid_t, int64_t>& thread_node_to_interval);

    // Get intervals in traversal that are not covered according to this->node_to_interval or
    // the thread_node_to_interval parameter.
    vector<pair<int64_t, int64_t>> get_uncovered_intervals(const vector<step_handle_t>& trav,
                                                           const unordered_map<nid_t, int64_t>& thread_node_to_interval);

    // Add a new interval into the altpath_intervals vector and update the node_to_interval map.
    // If the interval can be merged into an existing, contiguous interval, do that instead.
    // Returns true if a new interval was added, false if an existing interval was updated.
    bool add_interval(vector<pair<step_handle_t, step_handle_t>>& thread_altpath_intervals,
                      unordered_map<nid_t, int64_t>& thread_node_to_interval,
                      const pair<step_handle_t, step_handle_t>& new_interval,
                      bool global = false);

    // add_interval() can delete an existing interval. This requires a full update at the end.
    void defragment_intervals();

    // Get the total coverage of a traversal (sum of step lengths * path count).
    int64_t get_coverage(const vector<step_handle_t>& trav, const pair<int64_t, int64_t>& uncovered_interval);

    // Make sure all nodes in all altpaths are in forward orientation.
    // This is always possible because they are, by definition, disjoint.
    // This should only be run from inside apply().
    void forwardize_altpaths(MutablePathMutableHandleGraph* mutable_graph);

    // Search back to the reference and return <distance, node_id> when found.
    // (here distance is the number of intervals crossed, aka rank)
    // "first" toggles returning the first interval found vs all of them.
    vector<pair<int64_t, nid_t>> get_reference_nodes(nid_t node_id, bool first) const;

protected:

    const PathHandleGraph* graph = nullptr;

    // Intervals are end-exclusive (like BED).
    vector<pair<step_handle_t, step_handle_t>> altpath_intervals;

    // altpath_intervals[0, num_ref_intervals-1] are all rank-0 reference intervals.
    int64_t num_ref_intervals = 0;

    // Map from node ID to interval index.
    unordered_map<nid_t, int64_t> node_to_interval;

    // Counter for generating unique altpath indices per base path.
    // Using mutable so it can be updated in apply() which is logically const for the cover.
    mutable unordered_map<string, int64_t> base_path_alt_counter;

    // Used when selecting traversals to make the greedy cover.
    struct RankedFragment {
        int64_t coverage;
        const string* name;
        int64_t trav_idx;
        pair<int64_t, int64_t> fragment;
        bool operator<(const RankedFragment& f2) const {
            // note: name comparison is flipped because we want to select high coverage / low name
            return this->coverage < f2.coverage || (this->coverage == f2.coverage && *this->name > *f2.name);
        }
    };
};

}

#endif
