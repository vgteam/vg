#ifndef VG_GREF_HPP_INCLUDED
#define VG_GREF_HPP_INCLUDED

/**
 * \file gref.hpp
 *
 * Interface for computing and querying graph reference path covers.
 *
 * An gref cover is a set of path fragments (stored as separate paths) in the graph.
 * They are always relative to an existing reference sample (ie GRCh38 or CHM13).
 * Unlike rGFA paths which use complex metadata embedding, gref paths use a simple naming
 * scheme: the base reference path name moved into the "gref_" namespace, with the
 * fragments numbered off it as {gref_base_path_name}_{N}_alt
 *
 * For example, if the reference path is "CHM13#0#chr1", the cover is written as:
 *   - gref_CHM13#0#chr1        (a copy of the reference path itself)
 *   - gref_CHM13#0#chr1_1_alt  (fragments hanging off it)
 *   - gref_CHM13#0#chr1_2_alt
 *   - etc.
 *
 * The naming is a convention, not an option: because the prefix lands on the sample,
 * a graph carries both references at once, "is this a gref path/sample" is a prefix
 * test, and the base path a gref path came from is recoverable by dropping the prefix.
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

class GrefCover {
public:
    // The prefix every gref path name carries
    static const string gref_prefix;  // "gref_"

    // The suffix used to identify gref fragments
    static const string gref_suffix;  // "_alt"

    // Name of the gref copy of a path: the same name moved into the gref namespace,
    // keeping its subrange so that subpaths of one contig stay distinct and keep their
    // coordinates.  Example: make_gref_copy_name("CHM13#0#chr1[100-200]")
    //                        -> "gref_CHM13#0#chr1[100-200]"
    // For PanSN names the prefix renames the sample, which is what makes "is this path
    // a gref path" and "which base path is it derived from" answerable from the name
    // alone.  Anything matching a base path against its gref copy must use this.
    static string make_gref_copy_name(const string& path_name);

    // Name a fragment's gref base is built from: make_gref_copy_name() with the subrange
    // dropped as well, since make_gref_name() appends "_{N}_alt" and that has to land on
    // the locus.  Fragments off different subpaths of one contig therefore share a base
    // name; the per-base counter keeps their full names distinct.
    // Example: make_gref_base_name("CHM13#0#chr1[100-200]") -> "gref_CHM13#0#chr1"
    static string make_gref_base_name(const string& base_path_name);

    // Test if a name is in the gref namespace. Works on path names and, because
    // the prefix lands on the sample for PanSN names, on sample names too.
    // Example: is_gref_derived("gref_CHM13#0#chr1_3_alt") -> true
    //          is_gref_derived("gref_CHM13") -> true
    static bool is_gref_derived(const string& name);

    // Take a name back out of the gref namespace. Returns the input unchanged
    // if it isn't in it.
    // Example: strip_gref_prefix("gref_CHM13#0#chr1") -> "CHM13#0#chr1"
    static string strip_gref_prefix(const string& name);

    // Create an gref path name from a base reference path name and an index.
    // Example: make_gref_name("gref_CHM13#0#chr1", 1) -> "gref_CHM13#0#chr1_1_alt"
    static string make_gref_name(const string& base_path_name, int64_t gref_index);

    // Test if a path name is an gref fragment (contains "_{N}_alt" suffix).
    // Copied base contigs are gref paths but not fragments, so they fail this.
    static bool is_gref_name(const string& path_name);

    // Parse an gref path name to extract the base reference path name.
    // Returns the original base path name, or the input if not an gref path.
    // Example: parse_base_path("gref_CHM13#0#chr1_3_alt") -> "gref_CHM13#0#chr1"
    static string parse_base_path(const string& gref_name);

    // Parse an gref path name to extract the gref index.
    // Returns -1 if the path is not an gref path.
    // Example: parse_gref_index("gref_CHM13#0#chr1_3_alt") -> 3
    static int64_t parse_gref_index(const string& gref_name);

public:
    // Clear out any existing gref paths from the graph. Recommended to run this
    // before compute().
    void clear(MutablePathMutableHandleGraph* graph);

    // Compute the gref cover from the graph, starting with a given set of reference paths.
    void compute(const PathHandleGraph* graph,
                 SnarlManager* snarl_manager,
                 const unordered_set<path_handle_t>& reference_paths,
                 int64_t minimum_length);

    // Load existing gref paths from the graph, assuming they've been computed already.
    // The reference_paths should be the rank-0 paths the gref paths extend from.
    void load(const PathHandleGraph* graph,
              const unordered_set<path_handle_t>& reference_paths);

    // Apply the gref cover to a graph (must have been computed first).  Everything
    // lands in the gref namespace:
    // 1. The base reference paths are copied over (CHM13#0#chr1 -> gref_CHM13#0#chr1)
    // 2. The fragments are added alongside them (gref_CHM13#0#chr1_1_alt, etc.)
    // The original paths are left untouched, so the graph carries both the base
    // reference and the gref reference and they can be told apart by name.
    void apply(MutablePathMutableHandleGraph* mutable_graph);

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

    // Write a tab-separated table describing gref segments.
    // Each line contains: source_path, source_start, source_end, gref_path_name,
    //                     ref_path, ref_start, ref_end, level, parent
    //
    // "level" is how many gref hops separate this fragment from a base reference path:
    // 1 means it hangs directly off the reference, 2 that it is only reachable across
    // another fragment, and so on.  It is the same number the VCF reports as INFO/CH.
    // "parent" names the path it hangs off -- a gref copy of a reference path, or another
    // fragment -- so the two columns form a forest and level is derivable by chasing parent.
    // Both are "." for a fragment in a component with no reference path.
    //
    // Note "level" is a hop count over interval adjacency, not snarl containment: it says
    // how far the fragment is from the reference, not how deeply nested it is.  The two
    // agree at level 1, which is the case worth filtering on.
    //
    // Must be called after compute() and knows what gref path names will be used.
    void write_gref_segments(ostream& os);

protected:

    // Compute the cover for the given snarl, by greedily finding the covered paths through it.
    // The cover is added to the two "thread_" structures.
    // top_snarl_start/end are the boundary node IDs of the top-level snarl containing this snarl.
    void compute_snarl(const Snarl& snarl, PathTraversalFinder& path_trav_finder, int64_t minimum_length,
                       vector<pair<step_handle_t, step_handle_t>>& thread_gref_intervals,
                       unordered_map<nid_t, int64_t>& thread_node_to_interval,
                       nid_t top_snarl_start, nid_t top_snarl_end,
                       vector<pair<nid_t, nid_t>>& thread_snarl_bounds);

    // Get intervals in traversal that are not covered according to this->node_to_interval or
    // the thread_node_to_interval parameter.
    vector<pair<int64_t, int64_t>> get_uncovered_intervals(const vector<step_handle_t>& trav,
                                                           const unordered_map<nid_t, int64_t>& thread_node_to_interval);

    // Add a new interval into the gref_intervals vector and update the node_to_interval map.
    // If the interval can be merged into an existing, contiguous interval, do that instead.
    // Returns true if a new interval was added, false if an existing interval was updated.
    bool add_interval(vector<pair<step_handle_t, step_handle_t>>& thread_gref_intervals,
                      unordered_map<nid_t, int64_t>& thread_node_to_interval,
                      const pair<step_handle_t, step_handle_t>& new_interval,
                      bool global = false,
                      vector<pair<nid_t, nid_t>>* snarl_bounds_vec = nullptr,
                      pair<nid_t, nid_t> snarl_bounds = {0, 0});

    // add_interval() can delete an existing interval. This requires a full update at the end.
    void defragment_intervals();

    // Post-insertion cross-path merge: tries to consolidate the interval containing
    // ref_step with any cross-path neighbor at its left or right boundary.
    // Operates on this->gref_intervals and this->node_to_interval.
    void try_cross_path_merge(step_handle_t ref_step);

    // Remove non-reference intervals shorter than minimum_length, then defragment.
    // Called after all merging is complete so short intervals have had a chance to merge.
    void filter_short_intervals(int64_t minimum_length);

    // Walk forward from start_step on path, comparing each node ID + orientation
    // against other_interval's steps. Returns the new end step if all match, nullopt otherwise.
    optional<step_handle_t> try_extend_forward(step_handle_t start_step, path_handle_t path,
                                                const pair<step_handle_t, step_handle_t>& other_interval);

    // Walk backward from start_step on path, comparing each node ID + orientation
    // against other_interval's steps in reverse order.  Both sides are walked
    // lazily, so a mismatch costs only the steps actually compared.
    // Returns the first matching step if all match, nullopt otherwise.
    optional<step_handle_t> try_extend_backward(step_handle_t start_step, path_handle_t path,
                                                 const pair<step_handle_t, step_handle_t>& other_interval);

    // Check if merging two adjacent/overlapping step ranges on the same path
    // would produce a duplicate node ID.  Walks [interval_a.first, interval_b.second).
    bool merge_would_duplicate_node(const pair<step_handle_t, step_handle_t>& interval_a,
                                    const pair<step_handle_t, step_handle_t>& interval_b) const;

    // Fast duplicate check for global fold: walks [ext_start, ext_end) and
    // checks whether any node ID in that range already belongs to
    // target_interval_idx in nti.  The caller is responsible for pre-trimming
    // any shared boundary step (overlap-by-one) from the walk range so that
    // boundary nodes owned by the target are not flagged as duplicates.
    bool extension_would_duplicate_node(const unordered_map<nid_t, int64_t>& nti,
                                        int64_t target_interval_idx,
                                        step_handle_t ext_start, step_handle_t ext_end) const;

    // Unified duplicate-node check: dispatches to extension_would_duplicate_node
    // (O(extension_length)) when global=true, or merge_would_duplicate_node
    // (O(combined_length)) otherwise.  Callers must pre-trim shared boundary
    // steps from [ext_start, ext_end); merge_would_duplicate_node handles the
    // shared boundary naturally via its combined walk.
    bool would_duplicate_node(bool global,
                              const unordered_map<nid_t, int64_t>& nti,
                              int64_t target_idx,
                              step_handle_t ext_start, step_handle_t ext_end,
                              const pair<step_handle_t, step_handle_t>& interval_a,
                              const pair<step_handle_t, step_handle_t>& interval_b) const;

    // Get the total coverage of a traversal (sum of step lengths * path count).
    int64_t get_coverage(const vector<step_handle_t>& trav, const pair<int64_t, int64_t>& uncovered_interval);

    // Second pass: greedily cover any nodes not covered by snarl traversals.
    // This handles nodes that are outside of snarls or in complex regions
    // where the traversal finder couldn't find good coverage.
    void fill_uncovered_nodes(int64_t minimum_length);

    // For every interval, how many intervals must be crossed to reach a rank-0 reference
    // interval, and which interval it was reached from.  Reference intervals are (0, -1);
    // an interval in a component with no reference path is (-1, -1).
    //
    // One breadth-first sweep from all the reference intervals at once, rather than a
    // search per fragment: uncovered nodes are shared conduits, so a per-fragment search
    // re-walks them once per fragment and cannot terminate at all in a component with no
    // reference.  Uses each interval's own two ends, which also sidesteps node_to_interval
    // being last-write-wins where intervals overlap.
    vector<pair<int64_t, int64_t>> compute_interval_levels() const;

    // Search back to the reference and return <distance, node_id> when found.
    // (here distance is the number of intervals crossed, aka rank)
    // "first" toggles returning the first interval found vs all of them.
    vector<pair<int64_t, nid_t>> get_reference_nodes(nid_t node_id, bool first) const;

    // Debug function: verify that every node in the graph is covered by the gref cover.
    // Prints a summary of coverage statistics to stderr.
    void verify_cover(int64_t minimum_length) const;

    const PathHandleGraph* graph = nullptr;

    // Intervals are end-exclusive (like BED).
    vector<pair<step_handle_t, step_handle_t>> gref_intervals;

    // Top-level snarl boundary nodes for each interval, parallel to gref_intervals.
    // (0, 0) sentinel for reference intervals and fill_uncovered_nodes intervals.
    vector<pair<nid_t, nid_t>> interval_snarl_bounds;

    // gref_intervals[0, num_ref_intervals-1] are all rank-0 reference intervals.
    int64_t num_ref_intervals = 0;

    // Map from node ID to interval index.
    unordered_map<nid_t, int64_t> node_to_interval;

    // Counter for generating unique gref indices per base path.
    // Using mutable so it can be updated in apply() which is logically const for the cover.
    unordered_map<string, int64_t> base_path_gref_counter;

    // Whether to print verbose output (coverage summary, etc.)
    bool verbose = false;

    // When true, rank traversal fragments by name only (ignore coverage).
    // This ensures deterministic output regardless of thread count.
    bool rank_by_name = false;

    // Copy the base reference paths into the gref namespace.
    // Creates new paths like "gref_CHM13#0#chr1" from "CHM13#0#chr1".
    void copy_base_paths_to_gref(MutablePathMutableHandleGraph* mutable_graph,
                                 const unordered_set<path_handle_t>& reference_paths);

    // Used when selecting traversals to make the greedy cover.
    struct RankedFragment {
        int64_t coverage;
        int64_t length;
        bool reverse;
        const string* name;
        int64_t trav_idx;
        pair<int64_t, int64_t> fragment;
        bool operator<(const RankedFragment& f2) const {
            // Max-heap, so this is "worse than".  Longest first, because a short fragment
            // that displaces part of a longer one leaves both halves to be filtered away
            // by --min-gref-len -- which loses the sequence altogether.  Orientation only
            // breaks ties between runs of equal length: a fragment taken from a path that
            // walks it backwards has to be flipped to be written, so prefer one that does
            // not, but never at the cost of covering less in one piece.
            // note: name comparison is flipped because we want to select high coverage / low name
            if (this->coverage != f2.coverage) {
                return this->coverage < f2.coverage;
            }
            if (this->length != f2.length) {
                return this->length < f2.length;
            }
            if (this->reverse != f2.reverse) {
                return this->reverse;
            }
            return *this->name > *f2.name;
        }
    };
};

}

#endif
