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

#include <bdsg/snarl_distance_index.hpp>

#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"

namespace vg {

using namespace std;

class GrefCover {
public:
    // Whether the cover draws candidates from snarl traversals as well as from whole-path
    // runs.  A snarl traversal is bounded by its bubble, so it cannot offer a candidate
    // longer than one snarl; fill_uncovered_nodes() enumerates maximal runs along whole
    // source paths and covers everything on its own.  See the comment in compute().
    //
    // Public because compute() then has no use for its SnarlManager, and building one is
    // the single largest cost of `vg paths -u` in both time and peak memory.  Callers test
    // this to decide whether to build one at all: when it is false, compute() accepts
    // nullptr.  Keeping the two decisions on one constant is what stops them drifting.
    static constexpr bool use_snarl_candidates = false;

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
    //
    // distance_index, when non-null, is a decomposition-only snarl distance index over the
    // same graph.  It is not used to build the cover; it is used afterwards to map each
    // fragment to the top-level snarl containing it (see assign_top_level_snarls()).
    void compute(const PathHandleGraph* graph,
                 SnarlManager* snarl_manager,
                 const bdsg::SnarlDistanceIndex* distance_index,
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

    // Summary of how the final cover sits inside the top-level snarl decomposition.
    // Filled in by assign_top_level_snarls(); all zero when no distance index was supplied.
    //
    // The unit of the decomposition is a child of a top-level chain: either a top-level
    // snarl (a bubble hanging off the spine) or a bare spine node.  A fragment is
    // independently processable exactly when all of its nodes sit in one unit.
    struct TopLevelSnarlStats {
        // Shape of the decomposition
        int64_t top_level_chains = 0;
        int64_t top_level_snarls = 0;
        int64_t spine_nodes = 0;             // node children of top-level chains
        int64_t unexpected_chain_children = 0;  // neither snarl nor node: should be 0

        // Reference-anchoring conditions
        int64_t unanchored_snarls = 0;       // top-level snarls with a non-reference bound
        int64_t non_reference_spine_nodes = 0;
        int64_t non_reference_spine_bp = 0;

        // Where the fragments landed
        int64_t fragments = 0;
        int64_t fragments_in_one_snarl = 0;     // wholly inside a single top-level snarl
        int64_t fragments_on_one_spine_node = 0;
        int64_t fragments_spanning_snarls = 0;  // touch >= 2 top-level snarls: THE CLAIM
        int64_t fragments_touching_spine = 0;   // contain at least one spine node
        int64_t fragments_spanning_units = 0;   // touch >= 2 units of any kind
        int64_t fragments_with_unmapped_nodes = 0;  // decomposition did not place a node
        int64_t max_snarls_in_one_fragment = 0;
        int64_t max_units_in_one_fragment = 0;
        // How many distinct top-level snarls the fragments landed in.  Equal to
        // fragments_in_one_snarl when no two fragments share a snarl; far below it would
        // mean the decomposition is coarse, and 1 would mean the mapping is degenerate.
        int64_t distinct_snarls_with_fragments = 0;
        // Nodes the descent placed in two different units.  The decomposition is a
        // partition, so this must be 0; a non-zero value would mean the unit assignment is
        // ambiguous and every count above it is meaningless.
        int64_t nodes_claimed_twice = 0;
    };

    // Map every fragment onto the top-level snarl decomposition carried by distance_index.
    // Fills interval_snarl_bounds with the boundary nodes of the containing top-level snarl
    // for every fragment that lies wholly inside one, leaves {0,0} for the rest, and records
    // the summary returned by get_top_level_snarl_stats().
    void assign_top_level_snarls(const bdsg::SnarlDistanceIndex& distance_index);

    // Summary from the last assign_top_level_snarls() call.
    const TopLevelSnarlStats& get_top_level_snarl_stats() const;

    // The boundary nodes of the top-level snarl containing the given interval, or {0, 0}
    // when the interval is a reference path, when no distance index was supplied, or when
    // the interval is not wholly inside one top-level snarl.  Index is into get_intervals().
    pair<nid_t, nid_t> get_top_level_snarl(int64_t interval_index) const;

    // Write a tab-separated table describing gref segments.
    // Each line contains: source_path, source_start, source_end, gref_path_name,
    //                     ref_path, ref_start, ref_end
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

    // Search back to the reference and return <distance, node_id> when found.
    // (here distance is the number of intervals crossed, aka rank)
    // "first" toggles returning the first interval found vs all of them.
    vector<pair<int64_t, nid_t>> get_reference_nodes(nid_t node_id, bool first) const;

    // Resolve the gref base path name a fragment interval gets named after: the reference
    // path it traces back to, or its own source path when there is no reference in its
    // component.
    //
    // apply() and write_gref_segments() must agree here.  Both number fragments by
    // incrementing a per-base-name counter over the same interval sequence, and
    // write_gref_segments() runs first, so it cannot read the names apply() picks -- it has
    // to predict them.  A single disagreement does not just mislabel one row: it renumbers
    // every later fragment on both of the names involved.  So the prediction is not a
    // reimplementation, it is this same call.
    //
    // When the interval traces back to a reference, *out_ref_path and *out_ref_node are set
    // to the reference path and the reference node the trace landed on.  When it does not,
    // both are left untouched, so initialise out_ref_node to 0 to detect that case.
    string resolve_base_path_name(int64_t interval_index,
                                  path_handle_t* out_ref_path = nullptr,
                                  nid_t* out_ref_node = nullptr) const;

    // Debug function: verify that every node in the graph is covered by the gref cover.
    // Prints a summary of coverage statistics to stderr.
    void verify_cover(int64_t minimum_length) const;

    // Enforce node-disjointness: no two gref intervals may share a node, and no fragment
    // may claim a node belonging to a reference interval.
    //
    // This is a hard invariant, not a quality measure.  Every gref path becomes its own
    // reference contig with its own coordinate space, so a node owned by two of them is the
    // same sequence at two different coordinates -- downstream that is not a worse cover,
    // it is an incoherent one.  node_to_interval cannot enforce it on its own: it maps each
    // node to ONE interval index, so an interval that overwrites another's entry silently
    // wins the index while the loser keeps the step in its own range and still gets written
    // by apply().  Only walking the intervals themselves can detect that, which is what
    // this does.
    //
    // Called at the end of compute() over the final cover, so it holds whatever route the
    // intervals took to get there.  Exits non-zero on violation rather than warning: a
    // cover that violates this should never reach a downstream tool.
    void verify_disjoint() const;

    const PathHandleGraph* graph = nullptr;

    // Intervals are end-exclusive (like BED).
    vector<pair<step_handle_t, step_handle_t>> gref_intervals;

    // Top-level snarl boundary nodes for each interval, parallel to gref_intervals.
    // (0, 0) sentinel for reference intervals, and for any fragment that is not wholly
    // contained in one top-level snarl (or that was computed without a distance index).
    vector<pair<nid_t, nid_t>> interval_snarl_bounds;

    // Summary from assign_top_level_snarls().
    TopLevelSnarlStats top_level_snarl_stats;

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
            // Max-heap, so this is "worse than".  Rank by source path name.
            //
            // This looks arbitrary and is not: name order is *stable across snarls*.  The
            // same path wins in adjacent snarls, so add_interval()'s same-path merge joins
            // their intervals into one long fragment.  Ranking by run length instead is
            // locally optimal and globally destructive -- each snarl independently picks
            // whichever path happens to have the longest run just there, neighbours end up
            // on different paths, and nothing can merge them: try_cross_path_merge() only
            // consolidates paths that walk the stretch identically, which two paths taking
            // different traversals of a bubble by definition do not.
            //
            // note: name comparison is flipped because we want to select high coverage / low name
            return this->coverage < f2.coverage || (this->coverage == f2.coverage && *this->name > *f2.name);
        }
    };
};

}

#endif
