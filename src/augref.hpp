#ifndef VG_AUGREF_HPP_INCLUDED
#define VG_AUGREF_HPP_INCLUDED

/**
 * \file augref.hpp
 *
 * Interface for computing and querying augmented reference path covers.
 *
 * An augref cover is a set of path fragments (stored as separate paths) in the graph.
 * They are always relative to an existing reference sample (ie GRCh38 or CHM13).
 * Unlike rGFA paths which use complex metadata embedding, augref paths use a simple naming
 * scheme: {base_path_name}_{N}_alt
 *
 * For example, if the reference path is "CHM13#0#chr1", augref paths would be named:
 *   - CHM13#0#chr1_1_alt
 *   - CHM13#0#chr1_2_alt
 *   - etc.
 *
 * The data structures used in this class are always relative to the original paths
 * in the graph. The REFERENCE-sense fragments that are used to serialize the
 * cover can be created and loaded, but they are not used beyond that.
 */

#include <optional>

#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"

namespace vg {

using namespace std;

class AugRefCover {
public:
    // The suffix used to identify augref paths
    static const string augref_suffix;  // "_alt"

    // Create an augref path name from a base reference path name and an index.
    // Example: make_augref_name("CHM13#0#chr1", 1) -> "CHM13#0#chr1_1_alt"
    static string make_augref_name(const string& base_path_name, int64_t augref_index);

    // Test if a path name is an augref path (contains "_{N}_alt" suffix).
    static bool is_augref_name(const string& path_name);

    // Parse an augref path name to extract the base reference path name.
    // Returns the original base path name, or the input if not an augref path.
    // Example: parse_base_path("CHM13#0#chr1_3_alt") -> "CHM13#0#chr1"
    static string parse_base_path(const string& augref_name);

    // Parse an augref path name to extract the augref index.
    // Returns -1 if the path is not an augref path.
    // Example: parse_augref_index("CHM13#0#chr1_3_alt") -> 3
    static int64_t parse_augref_index(const string& augref_name);

public:
    // Clear out any existing augref paths from the graph. Recommended to run this
    // before compute().
    void clear(MutablePathMutableHandleGraph* graph);

    // Compute the augref cover from the graph, starting with a given set of reference paths.
    void compute(const PathHandleGraph* graph,
                 SnarlManager* snarl_manager,
                 const unordered_set<path_handle_t>& reference_paths,
                 int64_t minimum_length);

    // Load existing augref paths from the graph, assuming they've been computed already.
    // The reference_paths should be the rank-0 paths the augref paths extend from.
    void load(const PathHandleGraph* graph,
              const unordered_set<path_handle_t>& reference_paths);

    // Apply the augref cover to a graph (must have been computed first), adding it
    // as a bunch of REFERENCE-sense paths with the simplified naming scheme.
    // If augref_sample_name is set, base paths are first copied to the new sample,
    // and augref paths are created under the new sample name.
    void apply(MutablePathMutableHandleGraph* mutable_graph);

    // Set the sample name for augref paths. When set, apply() will:
    // 1. Copy base reference paths to this new sample (CHM13#0#chr1 -> new_sample#0#chr1)
    // 2. Create augref paths under the new sample (new_sample#0#chr1_1_alt, etc.)
    void set_augref_sample(const string& sample_name);

    // Get the current augref sample name (empty string if not set).
    const string& get_augref_sample() const;

    // Enable verbose output (coverage summary, etc.)
    void set_verbose(bool verbose);

    // Check if verbose output is enabled.
    bool get_verbose() const;

    // Get the rank (level) of a given node (0 if on a reference path).
    int64_t get_rank(nid_t node_id) const;

    // Get all computed intervals.
    const vector<pair<step_handle_t, step_handle_t>>& get_intervals() const;

    // Get an interval from a node. Returns nullptr if node not in an interval.
    const pair<step_handle_t, step_handle_t>* get_interval(nid_t node_id) const;

    // Get the number of reference intervals (rank-0).
    int64_t get_num_ref_intervals() const;

    // Write a tab-separated table describing augref segments.
    // Each line contains: source_path, source_start, source_end, augref_path_name,
    //                     ref_path, ref_start, ref_end
    // Must be called after compute() and knows what augref path names will be used.
    // If augref_sample is set, uses that for the augref path names.
    void write_augref_segments(ostream& os);

protected:

    // Compute the cover for the given snarl, by greedily finding the covered paths through it.
    // The cover is added to the two "thread_" structures.
    // top_snarl_start/end are the boundary node IDs of the top-level snarl containing this snarl.
    void compute_snarl(const Snarl& snarl, PathTraversalFinder& path_trav_finder, int64_t minimum_length,
                       vector<pair<step_handle_t, step_handle_t>>& thread_augref_intervals,
                       unordered_map<nid_t, int64_t>& thread_node_to_interval,
                       nid_t top_snarl_start, nid_t top_snarl_end,
                       vector<pair<nid_t, nid_t>>& thread_snarl_bounds);

    // Get intervals in traversal that are not covered according to this->node_to_interval or
    // the thread_node_to_interval parameter.
    vector<pair<int64_t, int64_t>> get_uncovered_intervals(const vector<step_handle_t>& trav,
                                                           const unordered_map<nid_t, int64_t>& thread_node_to_interval);

    // Add a new interval into the augref_intervals vector and update the node_to_interval map.
    // If the interval can be merged into an existing, contiguous interval, do that instead.
    // Returns true if a new interval was added, false if an existing interval was updated.
    bool add_interval(vector<pair<step_handle_t, step_handle_t>>& thread_augref_intervals,
                      unordered_map<nid_t, int64_t>& thread_node_to_interval,
                      const pair<step_handle_t, step_handle_t>& new_interval,
                      bool global = false,
                      vector<pair<nid_t, nid_t>>* snarl_bounds_vec = nullptr,
                      pair<nid_t, nid_t> snarl_bounds = {0, 0});

    // add_interval() can delete an existing interval. This requires a full update at the end.
    void defragment_intervals();

    // Remove non-reference intervals shorter than minimum_length, then defragment.
    // Called after all merging is complete so short intervals have had a chance to merge.
    void filter_short_intervals(int64_t minimum_length);

    // Walk forward from start_step on path, comparing each node ID + orientation
    // against other_interval's steps. Returns the new end step if all match, nullopt otherwise.
    optional<step_handle_t> try_extend_forward(step_handle_t start_step, path_handle_t path,
                                                const pair<step_handle_t, step_handle_t>& other_interval);

    // Collect other_interval's node IDs + orientations into a vector (forward walk).
    // Then walk backward from start_step on path, comparing in reverse.
    // Returns the first matching step if all match, nullopt otherwise.
    optional<step_handle_t> try_extend_backward(step_handle_t start_step, path_handle_t path,
                                                 const pair<step_handle_t, step_handle_t>& other_interval);

    // Get the total coverage of a traversal (sum of step lengths * path count).
    int64_t get_coverage(const vector<step_handle_t>& trav, const pair<int64_t, int64_t>& uncovered_interval);

    // Make sure all nodes in all augref paths are in forward orientation.
    // This is always possible because they are, by definition, disjoint.
    // This should only be run from inside apply().
    void forwardize_augref_paths(MutablePathMutableHandleGraph* mutable_graph);

    // Second pass: greedily cover any nodes not covered by snarl traversals.
    // This handles nodes that are outside of snarls or in complex regions
    // where the traversal finder couldn't find good coverage.
    void fill_uncovered_nodes(int64_t minimum_length);

    // Search back to the reference and return <distance, node_id> when found.
    // (here distance is the number of intervals crossed, aka rank)
    // "first" toggles returning the first interval found vs all of them.
    vector<pair<int64_t, nid_t>> get_reference_nodes(nid_t node_id, bool first) const;

    // Debug function: verify that every node in the graph is covered by the augref cover.
    // Prints a summary of coverage statistics to stderr.
    void verify_cover() const;

    const PathHandleGraph* graph = nullptr;

    // Intervals are end-exclusive (like BED).
    vector<pair<step_handle_t, step_handle_t>> augref_intervals;

    // Top-level snarl boundary nodes for each interval, parallel to augref_intervals.
    // (0, 0) sentinel for reference intervals and fill_uncovered_nodes intervals.
    vector<pair<nid_t, nid_t>> interval_snarl_bounds;

    // augref_intervals[0, num_ref_intervals-1] are all rank-0 reference intervals.
    int64_t num_ref_intervals = 0;

    // Map from node ID to interval index.
    unordered_map<nid_t, int64_t> node_to_interval;

    // Counter for generating unique augref indices per base path.
    // Using mutable so it can be updated in apply() which is logically const for the cover.
    mutable unordered_map<string, int64_t> base_path_augref_counter;

    // Optional sample name for augref paths. When set, base paths are copied to this
    // sample and augref paths are created under it.
    string augref_sample_name;

    // Whether to print verbose output (coverage summary, etc.)
    bool verbose = false;

    // When true, rank traversal fragments by name only (ignore coverage).
    // This ensures deterministic output regardless of thread count.
    bool rank_by_name = false;

    // Copy base reference paths to the augref sample.
    // Creates new paths like "new_sample#0#chr1" from "CHM13#0#chr1".
    void copy_base_paths_to_sample(MutablePathMutableHandleGraph* mutable_graph,
                                   const unordered_set<path_handle_t>& reference_paths);

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
